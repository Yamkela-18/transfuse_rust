// src/score.rs
// Rust-native Transfuse scoring replacement
// - TransRate-style, evidence-based, reference-free contig scoring
//
// TransRate (Smith-Unna et al. 2016, Genome Research) integrates four
// independent, reference-free signals per contig as a *product* of
// proportions, each in [0, 1]:
//
//   sCnuc  - per-base accuracy: quality-weighted match/mismatch identity at
//            every covered position (computed in bam.rs via samtools
//            mpileup, see ContigStats::nuc_score)
//   sCcov  - coverage completeness: fraction of the contig's length that
//            has at least one read aligned to it
//   sCord  - pair concordance: fraction of read pairs that align with the
//            expected orientation and insert size, tested against the
//            library's own empirical fragment-size distribution rather
//            than trusting the aligner's "properly paired" flag
//   sCseg  - coverage homogeneity: a Chow test comparing a one-segment
//            vs. best two-segment fit to the windowed depth profile, to
//            catch a jump between two plateaus (the signature of a
//            chimera - two transcripts, expressed at different levels,
//            fused into one contig)
//
//   contig_score = sCnuc * sCcov * sCord * sCseg
//
// All four are computed internally as a straight, equally-weighted product
// (matching TransRate's own integration), but ContigScore only exposes the
// same four fields the original Ruby Transfuse did - score, p_good,
// p_bases_covered, coverage - rather than persisting every intermediate
// component. p_good is TransRate's own per-fragment pass-rate metric (see
// pair_concordance_score below), not a copy of `score`.
//
// score_assemblies is, per-assembly, the most expensive step in the whole
// pipeline - alignment, four samtools passes, and a Salmon quantification
// run all happen inside run_alignment_coverage for every input assembly.
// Its progress bar is sized to PHASES_PER_ASSEMBLY * assembly_files.len(),
// not just assembly_files.len(), and ticks once per real sub-step
// (alignment, each of bam.rs's four passes, Salmon) rather than only once
// per whole assembly finishing - a single assembly can take minutes, and
// without per-phase ticks the bar would sit frozen at one percentage for
// that entire stretch despite real work happening underneath it. The bar
// shows percentage only (no N/M fraction) - this is the only progress bar
// in the program; main.rs has no separate overall pipeline bar.

use anyhow::{bail, Result};
use csv::{ReaderBuilder, WriterBuilder};
use indicatif::{MultiProgress, ProgressBar, ProgressStyle};
use log::info;
use rayon::prelude::*;
use serde::Deserialize;

use std::collections::HashMap;
use std::fs::File;
use std::io::BufReader;
use std::path::{Path, PathBuf};

use crate::aligner;
use crate::bam;
use crate::quant;

/// Sub-steps ticked per assembly: aligning reads, bam.rs's four passes
/// (bam::PHASES), and Salmon quantification.
const PHASES_PER_ASSEMBLY: u64 = 1 + bam::PHASES + 1;

#[derive(Debug, Clone, Default)]
pub struct ContigScore {
    pub score: f64,           // final TransRate-style contig score (product of 4 terms)
    pub p_good: f64,          // TransRate's own metric: fraction of individual fragments
                              // (read pairs) mapped to this contig classified "good" -
                              // one of the 4 multiplicative factors in `score`, computed
                              // here via pair_concordance_score (see caveat there on how
                              // this approximates the original's per-fragment classifier)
    pub p_bases_covered: f64, // sCcov: fraction of contig length with read coverage
    pub coverage: f64,        // mean depth (diagnostic only, not part of the score)
}

pub type ScoreMap = HashMap<String, ContigScore>;

#[derive(Debug, Deserialize)]
struct ScoreRow {
    #[serde(rename = "contig_name")]
    name: String,
    score: f64,
    p_good: f64,
    p_bases_covered: f64,
    coverage: f64,
}

// ============================================================
// Main scoring
// ============================================================

pub fn score_assemblies(
    assembly_files: &[PathBuf],
    left: &Path,
    right: &Path,
    threads: usize,
    verbose: bool,
    prefix_keys: bool,
    multi: &MultiProgress,
) -> Result<ScoreMap> {
    if assembly_files.is_empty() {
        bail!("No assemblies supplied");
    }

    let asm_threads = ((threads as f64 / assembly_files.len() as f64).ceil() as usize).max(1);

    info!(
        "Scoring {} assemblies using {} threads",
        assembly_files.len(),
        asm_threads
    );

    // Thread-safe (ProgressBar is Send+Sync by design) - shared by
    // reference across the rayon closure below. Sized to the number of
    // real sub-steps across all assemblies, not just the assembly count,
    // and ticked once per sub-step (see PHASES_PER_ASSEMBLY) so the
    // percentage advances continuously rather than jumping only when an
    // entire assembly finishes. Percentage only in the template - no
    // N/M fraction.
    let pb = multi.add(ProgressBar::new(
        assembly_files.len() as u64 * PHASES_PER_ASSEMBLY,
    ));
    pb.set_style(
        ProgressStyle::with_template(
            "  scoring [{bar:30.green/white}] ({percent}%) {msg}",
        )
        .unwrap()
        .progress_chars("=>-"),
    );

    let results: Vec<Result<ScoreMap>> = assembly_files
        .par_iter()
        .enumerate()
        .map(|(idx, assembly)| {
            let contig_prefix = format!("contig{idx}");
            let assembly_name = assembly
                .file_name()
                .map(|n| n.to_string_lossy().into_owned())
                .unwrap_or_else(|| assembly.display().to_string());

            let mut on_phase = |phase: &str| {
                pb.set_message(format!("{assembly_name}: {phase}"));
                pb.inc(1);
            };

            let scored =
                run_alignment_coverage(assembly, left, right, asm_threads, verbose, &mut on_phase)?;

            let mut out = ScoreMap::new();
            for (id, s) in scored {
                let key = if prefix_keys {
                    format!("{}_{}", contig_prefix, id)
                } else {
                    id
                };
                out.insert(key, s);
            }
            Ok(out)
        })
        .collect();

    pb.finish_and_clear();
    multi.remove(&pb);

    let mut combined = ScoreMap::new();
    for r in results {
        combined.extend(r?);
    }
    Ok(combined)
}

// ============================================================
// Alignment + coverage -> the four TransRate-style components
// ============================================================

fn run_alignment_coverage(
    fasta: &Path,
    left: &Path,
    right: &Path,
    threads: usize,
    verbose: bool,
    on_phase: &mut dyn FnMut(&str),
) -> Result<HashMap<String, ContigScore>> {
    on_phase("aligning reads");
    let bam_path = aligner::align_reads(fasta, left, right, threads, verbose)?;

    // bam::compute_coverage calls on_phase itself, once per its own four
    // passes (samtools coverage, fragment-size sampling, pair concordance,
    // mpileup) - see bam.rs's PHASES constant, which PHASES_PER_ASSEMBLY
    // above stays in sync with.
    let stats = bam::compute_coverage(&bam_path, fasta, on_phase)?;

    on_phase("salmon quantification");
    // Ruby's `coverage` column is Salmon's effective-length-normalized
    // estimate (eff_count * read_length / eff_len), not raw samtools depth -
    // these measure genuinely different things, so it's computed
    // separately here rather than substituted from bam::ContigStats.
    let read_length = quant::estimate_read_length(left)?;
    let salmon = quant::run_salmon_quant(&bam_path, fasta, threads)?;

    let mut out = HashMap::new();

    for (id, s) in stats {
        let cov_score = if s.length > 0 {
            (s.covered_bases as f64 / s.length as f64).clamp(0.0, 1.0)
        } else {
            0.0
        };

        let nuc_score = s.nuc_score.clamp(0.0, 1.0);
        let ord_score = pair_concordance_score(s.pairs_total, s.pairs_proper);
        let seg_score = homogeneity_score(&s.depth_windows);

        // TransRate integrates its four components as a straight, equally
        // weighted product of proportions - a contig only scores well if it
        // passes every check. A single badly-failed component (e.g. a clear
        // chimera signal) should tank the whole score, not get averaged away.
        //
        // The reference CSV from the original tool floors this at 0.01 -
        // 47/162 of its contigs have `score` exactly 0.01 while `p_good` for
        // those same rows is genuine 0.0, confirming the floor applies only
        // to the combined score (likely so a downstream geometric-mean-style
        // assembly score never multiplies in a literal zero), not to the
        // individual components. Matched here the same way.
        let final_score = (nuc_score * cov_score * ord_score * seg_score).max(0.01);

        let coverage = salmon
            .get(&id)
            .map(|&(eff_len, num_reads)| quant::effective_coverage(eff_len, num_reads, read_length))
            .unwrap_or(0.0);

        out.insert(
            id,
            ContigScore {
                score: final_score,
                // Real TransRate's p_good is its own fragment-level pass rate, one of
                // the four factors multiplied into `score` - not a copy of `score`.
                // ord_score (pairs_proper/pairs_total) is the closest thing we compute
                // to that per-fragment classification.
                p_good: ord_score,
                p_bases_covered: cov_score,
                coverage,
            },
        );
    }

    Ok(out)
}

// ============================================================
// The four score components (computed here, not persisted on ContigScore)
// ============================================================

/// sCord - pair concordance: fraction of read pairs that bam.rs judged
/// concordant (opposite-strand mates, insert size consistent with the
/// library's own empirical fragment-size distribution). Discordant pairs
/// indicate a structural misassembly - e.g. two regions stitched together
/// that don't belong.
fn pair_concordance_score(pairs_total: u64, pairs_proper: u64) -> f64 {
    if pairs_total == 0 {
        return 0.0;
    }
    (pairs_proper as f64 / pairs_total as f64).clamp(0.0, 1.0)
}

/// sCseg - coverage homogeneity via a Chow-test-style breakpoint check:
/// compares the residual sum of squares of a single-mean model for the
/// window-depth profile against the best-fitting two-segment (one
/// breakpoint) model. A large, well-supported improvement from allowing a
/// split is the signature of a chimera (two transcripts, expressed at
/// different levels, fused into one contig).
fn homogeneity_score(windows: &[f64]) -> f64 {
    let n = windows.len();
    if n < 4 {
        return 1.0; // too few windows for a meaningful two-segment fit
    }

    let mean = windows.iter().sum::<f64>() / n as f64;
    let rss1: f64 = windows.iter().map(|d| (d - mean).powi(2)).sum();

    if rss1 <= 1e-9 {
        return 1.0; // perfectly flat - no evidence of any breakpoint
    }

    let mut best_rss2 = rss1;
    for split in 1..n {
        let (left, right) = windows.split_at(split);
        let mean_l = left.iter().sum::<f64>() / left.len() as f64;
        let mean_r = right.iter().sum::<f64>() / right.len() as f64;
        let rss2: f64 = left.iter().map(|d| (d - mean_l).powi(2)).sum::<f64>()
            + right.iter().map(|d| (d - mean_r).powi(2)).sum::<f64>();
        if rss2 < best_rss2 {
            best_rss2 = rss2;
        }
    }

    // Chow-test F statistic: one extra parameter (the second segment's
    // mean), n-2 residual degrees of freedom for the two-segment model.
    let df2 = (n - 2).max(1) as f64;
    let f_stat = (rss1 - best_rss2) / (best_rss2 / df2).max(1e-9);

    // Map the F statistic to [0,1]: large F (strong evidence for a real
    // breakpoint) -> low score. This isn't a calibrated p-value, but it is
    // monotonic in the strength of evidence, which a raw min/max window
    // ratio wasn't.
    (1.0 / (1.0 + f_stat)).clamp(0.0, 1.0)
}

// ============================================================
// CSV output - matches the original Ruby format:
//   contig_name  score  p_good  p_bases_covered  coverage
//
// Ruby's Float#to_s prints the shortest decimal that round-trips back to
// the exact same float - so score/p_good/p_bases_covered (never rounded in
// Ruby) come out with anywhere from 1 to 17 decimal digits, while coverage
// (which Ruby rounds with `.round(2)` before printing) always has 1-2. A
// fixed `{:.6}` format matches neither: it truncates the unrounded fields
// and pads/misrepresents the rounded one. `ruby_float_str` below
// reproduces both behaviours: Rust's default f64 Display already produces
// the same shortest-round-trip string Ruby does, it just omits the
// trailing ".0" on whole numbers (`5.0` -> "5"), which this restores.
// ============================================================

fn ruby_float_str(x: f64) -> String {
    let s = format!("{}", x);
    if s.contains('.') || s.contains('e') || s.contains("inf") || s.contains("NaN") {
        s
    } else {
        format!("{}.0", s)
    }
}

/// Matches Ruby's `coverage.round(2)` - round to 2 decimal places, then
/// format with the same shortest-round-trip + trailing-".0" rule as the
/// unrounded fields (so 2664.60 prints as "2664.6", matching the original).
fn round2_str(x: f64) -> String {
    ruby_float_str((x * 100.0).round() / 100.0)
}

pub fn write_scores_csv(scores: &ScoreMap, output: &Path) -> Result<PathBuf> {
    let stem = output.file_stem().unwrap().to_string_lossy();
    let path = output.with_file_name(format!("{}_scores.csv", stem));

    let mut writer = WriterBuilder::new().delimiter(b'\t').from_path(&path)?;

    writer.write_record([
        "contig_name",
        "score",
        "p_good",
        "p_bases_covered",
        "coverage",
    ])?;

    let mut rows: Vec<_> = scores.iter().collect();
    rows.sort_by(|a, b| b.1.score.partial_cmp(&a.1.score).unwrap());

    for (id, s) in rows {
        writer.write_record([
            id,
            &ruby_float_str(s.score),
            &ruby_float_str(s.p_good),
            &ruby_float_str(s.p_bases_covered),
            &round2_str(s.coverage),
        ])?;
    }

    writer.flush()?;
    Ok(path)
}

// ============================================================
// Load scores
// ============================================================

pub fn load_scores_from_csv(files: &[PathBuf]) -> Result<ScoreMap> {
    let mut map = ScoreMap::new();

    for file in files {
        let reader = ReaderBuilder::new()
            .delimiter(b'\t')
            .has_headers(true)
            .from_reader(BufReader::new(File::open(file)?));

        for row in reader.into_deserialize::<ScoreRow>() {
            let row = row?;
            map.insert(
                row.name,
                ContigScore {
                    score: row.score,
                    p_good: row.p_good,
                    p_bases_covered: row.p_bases_covered,
                    coverage: row.coverage,
                },
            );
        }
    }

    Ok(map)
}
