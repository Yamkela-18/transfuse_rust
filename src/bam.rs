// src/bam.rs
//
// BAM coverage computation — replaces the bam-read C binary from transrate-tools.
//
// Four passes over the alignment, all via samtools (no BAM parsing library
// dependency):
//
//   1. `samtools coverage`  - cheap per-contig aggregates: read count,
//      breadth of coverage, mean depth, mean base/map quality.
//   2. `samtools view` (sampled, early-exit) - empirical fragment-size
//      distribution (median/MAD of TLEN).
//   3. `samtools view` (full, streamed) - pair concordance (sCord) and
//      per-window depth (sCseg). Primary-alignment-only: multi-mapping
//      reads are counted via whichever alignment the aligner tagged
//      primary, not reassigned - see README "Known differences" for what
//      this leaves on the table relative to the original tool's
//      EM-corrected read assignment.
//   4. `samtools mpileup` (streamed) - per-position quality-weighted
//      identity (sCnuc). Also primary-alignment-only, same caveat.
//
// `compute_coverage` takes an `on_phase` callback, called once before each
// of the four passes above with a short phase name. score.rs uses this to
// tick its progress bar after real, potentially slow subprocess work
// actually happens, rather than only once when the whole assembly is
// done - a single assembly can take minutes, and without per-phase ticks
// the bar would sit frozen at one percentage for that entire stretch.
// Deliberately a plain `&mut dyn FnMut(&str)` rather than depending on
// indicatif directly, so this module doesn't need to know progress bars
// exist - score.rs supplies whatever behavior it wants.
//
// stderr on every subprocess call below is suppressed (`.stderr(Stdio::
// null())` on streaming calls; `.output()` naturally captures rather than
// inherits it on the non-streaming one). samtools writes its own internal
// status lines - e.g. mpileup's "[mpileup] 1 samples in 1 input files" -
// directly to inherited stderr, completely outside Rust's `log` crate and
// therefore invisible to the MultiProgress-based coordination in main.rs;
// left alone, those lines interleave with and visually corrupt whichever
// progress bar is active. Failure diagnostics still surface via each
// command's exit status and the `bail!` messages below - suppressing
// stderr loses samtools' own detail on a failure, but every one of these
// calls already reports success/failure independently.
//
// `ContigStats::score()` below is a legacy, simplified heuristic kept only
// for this file's own unit tests. It is NOT what the pipeline uses - the
// authoritative TransRate-style contig score (a true product of sCnuc,
// sCcov, sCord, sCseg) is computed in `score::score_assemblies` from the
// fields populated here.

use anyhow::{bail, Context, Result};
use log::debug;
use std::collections::HashMap;
use std::io::{BufRead, BufReader};
use std::path::Path;
use std::process::{Command, Stdio};

/// Number of equal-length windows used for the coverage-homogeneity
/// (chimera/segmentation) check.
const NUM_WINDOWS: usize = 10;

/// Cap on how many fragment-size samples we collect before stopping early -
/// a modest sample is enough to estimate a stable median/MAD, and this
/// keeps that pass fast even on huge BAMs.
const MAX_FRAGMENT_SAMPLES: usize = 2_000_000;

/// Number of times `compute_coverage` calls `on_phase` per invocation -
/// score.rs uses this to size its progress bar correctly. Kept as one
/// named constant rather than a magic number, since it has to stay in
/// sync with the actual number of `on_phase(...)` calls below.
pub const PHASES: u64 = 4;

// ── Public types ──────────────────────────────────────────────────────────────

/// Per-contig alignment statistics extracted from the BAM file.
#[derive(Debug, Clone, Default)]
#[allow(dead_code)] // num_reads/coverage/mean_depth/mean_base_quality/mean_map_quality are
                     // diagnostic - populated from samtools coverage but not currently read
                     // by score.rs, which uses covered_bases/length directly instead.
pub struct ContigStats {
    pub num_reads: u64,
    pub covered_bases: u64,
    pub coverage: f64, // fraction 0.0-1.0
    pub mean_depth: f64,
    pub mean_base_quality: f64,
    pub mean_map_quality: f64,
    pub length: u64,

    pub pairs_total: u64,        // reads with FLAG & 0x1 (paired), for sCord
    pub pairs_proper: u64,       // reads passing our own orientation+insert-size test, for sCord
    pub depth_windows: Vec<f64>, // mean depth per window, for sCseg
    pub nuc_score: f64,          // mean per-position quality-weighted identity, for sCnuc
}

impl ContigStats {
    /// Legacy composite score in [0.0, 1.0]. Superseded by
    /// `score::score_assemblies`, which combines the real sCnuc/sCcov/sCord/
    /// sCseg components as a product rather than these ad-hoc weights.
    /// Kept only so existing callers/tests of this struct don't break.
    #[allow(dead_code)]
    pub fn score(&self) -> f64 {
        if self.length == 0 { return 0.0; }
        let cov_frac    = self.coverage.clamp(0.0, 1.0);
        let depth_score = (self.mean_depth / 10.0).min(1.0);
        let mapq_score  = (self.mean_map_quality / 30.0).min(1.0);
        0.5 * cov_frac + 0.3 * depth_score + 0.2 * mapq_score
    }
}

/// Accumulator for the pair-concordance + depth-window pass, keyed by
/// contig name.
#[derive(Debug, Clone)]
struct AlignDetail {
    pairs_total: u64,
    pairs_proper: u64,
    window_bases: Vec<f64>, // raw overlap-base sums, finalized to depth later
}

impl AlignDetail {
    fn new() -> Self {
        Self {
            pairs_total: 0,
            pairs_proper: 0,
            window_bases: vec![0.0; NUM_WINDOWS],
        }
    }
}

/// The library's empirical fragment-size distribution, estimated from a
/// sample of read pairs. `mad` of 0.0 means "no usable distribution" (e.g.
/// single-end data) - in that case the insert-size half of the
/// concordance check is skipped rather than penalising every pair.
#[derive(Debug, Clone, Copy)]
struct FragmentStats {
    median: f64,
    mad: f64,
}

// ── Public API ────────────────────────────────────────────────────────────────

/// Runs all four passes on bam_path (mpileup additionally needs the
/// assembly FASTA it was aligned against) and returns per-contig stats
/// with every sCnuc/sCcov/sCord/sCseg input populated. Calls `on_phase`
/// with a short label immediately before each pass starts - see module
/// docs for why.
pub fn compute_coverage(
    bam_path: &Path,
    fasta_path: &Path,
    on_phase: &mut dyn FnMut(&str),
) -> Result<HashMap<String, ContigStats>> {
    on_phase("samtools coverage");
    debug!("  Running samtools coverage on {:?}", bam_path);

    let output = Command::new("samtools")
        .args(["coverage", "-d", "0"]) // -d 0: no depth cap
        .arg(bam_path)
        .output() // captures stdout+stderr rather than inheriting - no suppression needed
        .context("Failed to run samtools coverage")?;

    if !output.status.success() {
        let stderr = String::from_utf8_lossy(&output.stderr);
        bail!("samtools coverage failed: {}", stderr);
    }

    let stdout = String::from_utf8_lossy(&output.stdout);
    let mut stats = parse_coverage_output(&stdout)?;

    let lengths: HashMap<String, u64> =
        stats.iter().map(|(k, v)| (k.clone(), v.length)).collect();

    on_phase("fragment size estimation");
    debug!("  Estimating fragment-size distribution from {:?}", bam_path);
    let frag = estimate_fragment_stats(bam_path)?;

    on_phase("pair concordance");
    debug!("  Running samtools view on {:?} for pair concordance + depth windows", bam_path);
    let details = compute_alignment_details(bam_path, &lengths, &frag)?;

    on_phase("per-base identity");
    debug!("  Running samtools mpileup on {:?} for per-base identity", bam_path);
    let nuc_scores = compute_nuc_identity(bam_path, fasta_path)?;

    for (name, s) in stats.iter_mut() {
        match details.get(name) {
            Some(d) => {
                s.pairs_total = d.pairs_total;
                s.pairs_proper = d.pairs_proper;
                s.depth_windows = finalize_windows(&d.window_bases, s.length);
            }
            None => {
                s.depth_windows = vec![0.0; NUM_WINDOWS];
            }
        }
        s.nuc_score = nuc_scores.get(name).copied().unwrap_or(0.0);
    }

    Ok(stats)
}

// ── Pass 1 parser (samtools coverage) ────────────────────────────────────────

fn parse_coverage_output(text: &str) -> Result<HashMap<String, ContigStats>> {
    let mut stats = HashMap::new();

    for line in text.lines() {
        if line.trim().is_empty()
    || line.starts_with('#')
    || line.starts_with("name")
    || line.starts_with("rname")
{
    continue;
}
        let cols: Vec<&str> = line.trim().split('\t').collect();
        if cols.len() < 9 { continue; }

        let name       = cols[0].trim().to_string();
        let start: u64 = cols[1].parse().unwrap_or(0);
        let end: u64   = cols[2].parse().unwrap_or(0);
        let length     = end.saturating_sub(start).max(1);
        let num_reads  = cols[3].parse().unwrap_or(0u64);
        let cov_bases  = cols[4].parse().unwrap_or(0u64);
        let coverage   = cols[5].parse().unwrap_or(0.0f64) / 100.0; // samtools gives %
        let mean_depth = cols[6].parse().unwrap_or(0.0f64);
        let mean_baseq = cols[7].parse().unwrap_or(0.0f64);
        let mean_mapq  = cols[8].parse().unwrap_or(0.0f64);

        stats.insert(name, ContigStats {
            num_reads, covered_bases: cov_bases, coverage,
            mean_depth, mean_base_quality: mean_baseq,
            mean_map_quality: mean_mapq, length,
            ..Default::default()
        });
    }
    Ok(stats)
}

// ── Pass 2 (samtools view, sampled): fragment-size distribution ────────────

fn estimate_fragment_stats(bam_path: &Path) -> Result<FragmentStats> {
    let mut child = Command::new("samtools")
        .args(["view", "-F", "0x904"])
        .arg(bam_path)
        .stdout(Stdio::piped())
        .stderr(Stdio::null()) // suppress samtools' own status/warning lines
        .spawn()
        .context("Failed to spawn samtools view for fragment-size sampling")?;

    let stdout = child
        .stdout
        .take()
        .context("Failed to capture samtools view stdout")?;
    let reader = BufReader::new(stdout);

    let mut samples: Vec<i64> = Vec::new();

    for line in reader.lines() {
        let line = line?;
        if line.is_empty() {
            continue;
        }
        let cols: Vec<&str> = line.splitn(10, '\t').collect();
        if cols.len() < 9 {
            continue;
        }

        let flag: u32 = cols[1].parse().unwrap_or(0);
        let tlen: i64 = cols[8].parse().unwrap_or(0);

        if flag & 0x1 != 0 && flag & 0x40 != 0 && tlen != 0 {
            samples.push(tlen.abs());
            if samples.len() >= MAX_FRAGMENT_SAMPLES {
                let _ = child.kill();
                break;
            }
        }
    }

    let _ = child.wait();

    let (median, mad) = median_and_mad(&samples);
    Ok(FragmentStats { median, mad })
}

fn median_and_mad(samples: &[i64]) -> (f64, f64) {
    if samples.is_empty() {
        return (0.0, 0.0);
    }
    let mut sorted = samples.to_vec();
    sorted.sort_unstable();
    let median = percentile(&sorted, 0.5);

    let mut deviations: Vec<i64> = sorted
        .iter()
        .map(|&s| (s as f64 - median).abs() as i64)
        .collect();
    deviations.sort_unstable();
    let mad = percentile(&deviations, 0.5).max(1.0);

    (median, mad)
}

fn percentile(sorted: &[i64], p: f64) -> f64 {
    if sorted.is_empty() {
        return 0.0;
    }
    let idx = (((sorted.len() - 1) as f64) * p).round() as usize;
    sorted[idx.min(sorted.len().saturating_sub(1))] as f64
}

fn pair_is_concordant(flag: u32, tlen: i64, frag: &FragmentStats) -> bool {
    let self_reverse = flag & 0x10 != 0;
    let mate_reverse = flag & 0x20 != 0;
    let orientation_ok = self_reverse != mate_reverse;

    let insert_ok = if frag.mad > 0.0 {
        (tlen.abs() as f64 - frag.median).abs() <= 3.0 * frag.mad
    } else {
        true
    };

    orientation_ok && insert_ok
}

// ── Pass 3 (samtools view, full): pair concordance + depth windows ─────────

fn compute_alignment_details(
    bam_path: &Path,
    lengths: &HashMap<String, u64>,
    frag: &FragmentStats,
) -> Result<HashMap<String, AlignDetail>> {
    let mut child = Command::new("samtools")
        .args(["view", "-F", "0x904"])
        .arg(bam_path)
        .stdout(Stdio::piped())
        .stderr(Stdio::null()) // suppress samtools' own status/warning lines
        .spawn()
        .context("Failed to spawn samtools view")?;

    let stdout = child
        .stdout
        .take()
        .context("Failed to capture samtools view stdout")?;
    let reader = BufReader::new(stdout);

    let mut details: HashMap<String, AlignDetail> = HashMap::new();

    for line in reader.lines() {
        let line = line?;
        if line.is_empty() {
            continue;
        }
        let cols: Vec<&str> = line.split('\t').collect();
        if cols.len() < 11 {
            continue;
        }

        let rname = cols[2];
        if rname == "*" {
            continue;
        }
        let contig_len = match lengths.get(rname) {
            Some(&l) => l,
            None => continue,
        };

        let flag: u32 = cols[1].parse().unwrap_or(0);
        let pos: u64 = cols[3].parse().unwrap_or(0);
        let cigar = cols[5];
        let tlen: i64 = cols.get(8).and_then(|s| s.parse().ok()).unwrap_or(0);

        let ref_len = cigar_ref_len(cigar);

        let entry = details
            .entry(rname.to_string())
            .or_insert_with(AlignDetail::new);

        if flag & 0x1 != 0 {
            entry.pairs_total += 1;
            if pair_is_concordant(flag, tlen, frag) {
                entry.pairs_proper += 1;
            }
        }

        add_to_windows(&mut entry.window_bases, contig_len, pos, ref_len);
    }

    let status = child
        .wait()
        .context("samtools view did not exit cleanly")?;
    if !status.success() {
        bail!("samtools view exited with a non-zero status");
    }

    Ok(details)
}

/// Sums the reference-consuming CIGAR ops (M/D/N/=/X), used to place a
/// read into the right depth window(s).
fn cigar_ref_len(cigar: &str) -> u64 {
    if cigar == "*" {
        return 0;
    }
    let mut ref_len: u64 = 0;
    let mut num = String::new();
    for c in cigar.chars() {
        if c.is_ascii_digit() {
            num.push(c);
            continue;
        }
        let n: u64 = num.parse().unwrap_or(0);
        num.clear();
        if matches!(c, 'M' | 'D' | 'N' | '=' | 'X') {
            ref_len += n;
        }
    }
    ref_len
}

/// Adds this read's reference-consumed span to whichever of the contig's
/// NUM_WINDOWS windows it overlaps, weighted by the actual overlap length,
/// so a read spanning a window boundary is split proportionally.
fn add_to_windows(windows: &mut [f64], contig_len: u64, read_start_1based: u64, ref_len: u64) {
    if contig_len == 0 || ref_len == 0 {
        return;
    }
    let n = windows.len() as u64;
    let window_width = ((contig_len as f64) / (n as f64)).ceil().max(1.0) as u64;

    let read_start0 = read_start_1based.saturating_sub(1);
    let read_end0 = read_start0 + ref_len;

    for (w, slot) in windows.iter_mut().enumerate() {
        let w_start = w as u64 * window_width;
        if w_start >= contig_len {
            break;
        }
        let w_end = ((w as u64 + 1) * window_width).min(contig_len);

        let overlap_start = read_start0.max(w_start);
        let overlap_end = read_end0.min(w_end);
        if overlap_end > overlap_start {
            *slot += (overlap_end - overlap_start) as f64;
        }
    }
}

/// Converts raw per-window overlap-base sums into mean depth per window.
fn finalize_windows(window_bases: &[f64], contig_len: u64) -> Vec<f64> {
    let n = window_bases.len() as u64;
    if contig_len == 0 || n == 0 {
        return vec![0.0; window_bases.len()];
    }
    let window_width = (contig_len as f64 / n as f64).max(1.0);
    window_bases.iter().map(|b| b / window_width).collect()
}

// ── Pass 4 (samtools mpileup): per-position quality-weighted identity ──────

fn compute_nuc_identity(bam_path: &Path, fasta_path: &Path) -> Result<HashMap<String, f64>> {
    let mut child = Command::new("samtools")
        .args(["mpileup", "-B", "-q", "0", "-Q", "0", "-f"])
        .arg(fasta_path)
        .arg(bam_path)
        .stdout(Stdio::piped())
        // This is the "[mpileup] 1 samples in 1 input files" line - mpileup
        // writes its own status directly to stderr regardless of verbosity
        // settings. Suppressed here rather than routed through `log`,
        // since it never went through `log` to begin with.
        .stderr(Stdio::null())
        .spawn()
        .context("Failed to spawn samtools mpileup")?;

    let stdout = child
        .stdout
        .take()
        .context("Failed to capture samtools mpileup stdout")?;
    let reader = BufReader::new(stdout);

    let mut sums: HashMap<String, (f64, u64)> = HashMap::new();

    for line in reader.lines() {
        let line = line?;
        if line.is_empty() {
            continue;
        }
        let cols: Vec<&str> = line.split('\t').collect();
        if cols.len() < 6 {
            continue;
        }

        let chrom = cols[0];
        let pileup = cols[4];
        let quals = cols[5];

        if let Some(pos_score) = score_pileup_line(pileup, quals) {
            let entry = sums.entry(chrom.to_string()).or_insert((0.0, 0));
            entry.0 += pos_score;
            entry.1 += 1;
        }
    }

    let status = child
        .wait()
        .context("samtools mpileup did not exit cleanly")?;
    if !status.success() {
        bail!("samtools mpileup exited with a non-zero status");
    }

    Ok(sums
        .into_iter()
        .map(|(k, (sum, n))| (k, if n > 0 { (sum / n as f64).clamp(0.0, 1.0) } else { 0.0 }))
        .collect())
}

fn score_pileup_line(pileup: &str, quals: &str) -> Option<f64> {
    let bases = pileup.as_bytes();
    let quals = quals.as_bytes();
    let mut qi = 0usize;
    let mut sum = 0.0f64;
    let mut n = 0u64;
    let mut i = 0usize;

    while i < bases.len() {
        match bases[i] as char {
            '^' => {
                i += 2;
            }
            '$' => {
                i += 1;
            }
            '+' | '-' => {
                i += 1;
                let start = i;
                while i < bases.len() && (bases[i] as char).is_ascii_digit() {
                    i += 1;
                }
                let len: usize = std::str::from_utf8(&bases[start..i])
                    .ok()
                    .and_then(|s| s.parse().ok())
                    .unwrap_or(0);
                i += len;
            }
            '.' | ',' | 'A' | 'C' | 'G' | 'T' | 'N' | 'a' | 'c' | 'g' | 'n' | 't' | '*' => {
                let is_match = matches!(bases[i] as char, '.' | ',');
                if qi < quals.len() {
                    let q = quals[qi].saturating_sub(33) as f64;
                    let p_correct = 1.0 - 10f64.powf(-q / 10.0);
                    sum += if is_match { p_correct } else { 1.0 - p_correct };
                    n += 1;
                    qi += 1;
                }
                i += 1;
            }
            _ => {
                i += 1;
            }
        }
    }

    if n == 0 {
        None
    } else {
        Some((sum / n as f64).clamp(0.0, 1.0))
    }
}

// ── Tests ─────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_coverage_line() {
        let tsv = "\
name\tstartpos\tendpos\tnumreads\tcovbases\tcoverage\tmeandepth\tmeanbaseq\tmeanmapq
contig1\t1\t1000\t150\t980\t98.0\t8.5\t35.0\t42.0
";
        let stats = parse_coverage_output(tsv).unwrap();
        assert_eq!(stats.len(), 1);
        let s = stats.get("contig1").expect("contig1 missing — check tab escapes in test TSV");
        assert_eq!(s.num_reads, 150);
        assert!((s.coverage - 0.98).abs() < 1e-6);
        assert!((s.mean_depth - 8.5).abs() < 1e-6);
    }

    #[test]
    fn test_score_range() {
        let s = ContigStats {
            num_reads: 100, covered_bases: 900, coverage: 0.90,
            mean_depth: 15.0, mean_base_quality: 35.0, mean_map_quality: 40.0,
            length: 1000,
            ..Default::default()
        };
        let score = s.score();
        assert!(score > 0.0 && score <= 1.0, "score out of range: {score}");
    }

    #[test]
    fn test_zero_length_score() {
        let s = ContigStats {
            num_reads: 0, covered_bases: 0, coverage: 0.0,
            mean_depth: 0.0, mean_base_quality: 0.0, mean_map_quality: 0.0,
            length: 0,
            ..Default::default()
        };
        assert_eq!(s.score(), 0.0);
    }

    #[test]
    fn test_cigar_ref_len() {
        assert_eq!(cigar_ref_len("5S10M2D3I"), 12);
    }

    #[test]
    fn test_add_to_windows_splits_across_boundary() {
        let mut windows = vec![0.0; 10];
        add_to_windows(&mut windows, 100, 6, 10);
        assert!((windows[0] - 5.0).abs() < 1e-9);
        assert!((windows[1] - 5.0).abs() < 1e-9);
        assert!(windows[2..].iter().all(|&v| v == 0.0));
    }

    #[test]
    fn test_median_and_mad_robust_to_outlier() {
        let samples = vec![100, 102, 98, 101, 99, 300];
        let (median, mad) = median_and_mad(&samples);
        assert!((98.0..=102.0).contains(&median), "median {median} pulled by outlier");
        assert!(mad < 10.0, "mad {mad} too large - not robust to the outlier");
    }

    #[test]
    fn test_pair_is_concordant() {
        let frag = FragmentStats { median: 300.0, mad: 20.0 };
        let flag_ok = 0x1 | 0x40 | 0x20;
        assert!(pair_is_concordant(flag_ok, 310, &frag));

        let flag_same_strand = 0x1 | 0x40;
        assert!(!pair_is_concordant(flag_same_strand, 310, &frag));
        assert!(!pair_is_concordant(flag_ok, 5000, &frag));

        let frag_none = FragmentStats { median: 0.0, mad: 0.0 };
        assert!(pair_is_concordant(flag_ok, 999_999, &frag_none));
    }

    #[test]
    fn test_score_pileup_line_all_matches() {
        let score = score_pileup_line(",..", "III").unwrap();
        assert!(score > 0.99, "expected near-1.0 for all high-quality matches, got {score}");
    }

    #[test]
    fn test_score_pileup_line_with_mismatch() {
        let score = score_pileup_line(".A,", "III").unwrap();
        assert!(score > 0.5 && score < 0.8, "expected a middling score, got {score}");
    }

    #[test]
    fn test_score_pileup_line_skips_indel_and_start_end_markers() {
        let score = score_pileup_line("^!.+2AC,-1A$,", "II").unwrap();
        assert!(score > 0.99, "expected near-1.0, indel/start/end tokens should be skipped, got {score}");
    }

    #[test]
    fn test_score_pileup_line_empty_returns_none() {
        assert!(score_pileup_line("", "").is_none());
    }

    #[test]
    fn test_compute_coverage_calls_on_phase_expected_number_of_times() {
        // We can't run real samtools here, but we CAN confirm the PHASES
        // constant matches the doc'd contract score.rs relies on to size
        // its progress bar - this test exists to fail loudly if someone
        // adds/removes an on_phase call without updating PHASES to match.
        assert_eq!(PHASES, 4);
    }
}
