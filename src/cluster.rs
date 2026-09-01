// src/cluster.rs
//
// Sequence clustering — wraps vsearch, equivalent to Transfuse#cluster.
//
// vsearch flags used:
//   --cluster_fast    greedy centroid-based clustering (USEARCH compatible)
//   --id              minimum sequence identity (0.0 – 1.0)
//   --msaout          per-cluster multiple sequence alignment
//   --uc              UC-format cluster membership table
//   --threads         parallelism
//   --notrunclabels   preserve full sequence headers
//
// Returns the path to the .aln MSA file consumed by consensus.rs.
//
// Progress: vsearch reports genuine live percentage progress per internal
// phase (reading the file, sorting by length, k-mer counting, clustering
// itself, sorting clusters, writing output) via carriage-return-updated
// stderr lines - `--quiet` was previously suppressing exactly that output.
// Dropping it and parsing those lines drives a real progress bar rather
// than a synthesized estimate. vsearch resets to 0-100% once per phase,
// not once for the whole run, so the bar visibly restarts a few times
// with a different phase label each time - that's vsearch's actual
// behavior, not a bug in how it's being read here.

use anyhow::{bail, Context, Result};
use indicatif::{MultiProgress, ProgressBar, ProgressStyle};
use log::{debug, info};
use std::io::Read;
use std::path::{Path, PathBuf};
use std::process::{Command, Stdio};

// ── Public API ────────────────────────────────────────────────────────────────

/// Run vsearch clustering on fasta_path at identity threshold.
/// Idempotent: skips if .aln and .clust output files already exist.
pub fn cluster_vsearch(
    fasta_path: &Path,
    identity: f64,
    threads: usize,
    multi: &MultiProgress,
) -> Result<PathBuf> {
    validate_identity(identity)?;

    let stem = fasta_path
        .file_stem()
        .map(|s| s.to_string_lossy().into_owned())
        .unwrap_or_else(|| "combined".into());
    let parent = fasta_path.parent().unwrap_or(Path::new("."));
    let id_str = format!("{:.2}", identity);

    let aln_path   = parent.join(format!("{stem}-{id_str}.aln"));
    let clust_path = parent.join(format!("{stem}-{id_str}.clust"));

    if aln_path.exists() && clust_path.exists() {
        info!("  vsearch output already exists, skipping clustering.");
        return Ok(aln_path);
    }

    info!(
        "  Running vsearch --cluster_fast --id {} --threads {} on {:?}",
        id_str, threads, fasta_path
    );

    let pb = multi.add(ProgressBar::new(100));
    pb.set_style(
        ProgressStyle::with_template("  clustering [{bar:30.blue/white}] ({percent}%) {msg}")
            .unwrap()
            .progress_chars("=>-"),
    );

    // --quiet deliberately dropped (vs. the original): that flag was
    // suppressing exactly the per-phase percentage output needed to drive
    // a real progress bar. Nothing else about the command changes.
    let mut child = Command::new("vsearch")
        .arg("--cluster_fast")
        .arg(fasta_path)
        .arg("--id").arg(&id_str)
        .arg("--msaout").arg(&aln_path)
        .arg("--uc").arg(&clust_path)
        .arg("--threads").arg(threads.to_string())
        .arg("--notrunclabels")
        .stderr(Stdio::piped())
        .spawn()
        .context("Failed to launch vsearch. Is it installed and on PATH?")?;

    let stderr = child
        .stderr
        .take()
        .context("Failed to capture vsearch stderr")?;
    stream_vsearch_progress(stderr, &pb);

    let status = child
        .wait()
        .context("vsearch did not exit cleanly")?;

    pb.finish_and_clear();
    multi.remove(&pb);

    if !status.success() {
        bail!("vsearch clustering failed for {:?}", fasta_path);
    }

    if !aln_path.exists() {
        bail!(
            "vsearch ran but produced no MSA file at {:?}.              The input may be empty or all sequences deduplicated.",
            aln_path
        );
    }

    debug!("  MSA written to {:?}", aln_path);
    debug!("  Cluster table written to {:?}", clust_path);
    Ok(aln_path)
}

// ── vsearch progress parsing ─────────────────────────────────────────────────

/// Reads vsearch's stderr byte-by-byte rather than line-by-line. vsearch
/// updates its progress in place using carriage returns (`\r`), only
/// emitting a real newline once a phase finishes - a standard
/// `BufRead::lines()` call (which splits on `\n`) would buffer an entire
/// phase's worth of `\r`-separated updates into one string and only hand
/// it over once that phase completed, defeating the point of live
/// progress. Splitting on both `\r` and `\n` surfaces each individual
/// update as vsearch actually writes it.
fn stream_vsearch_progress(mut stderr: impl Read, pb: &ProgressBar) {
    let mut buf = Vec::new();
    let mut byte = [0u8; 1];

    loop {
        match stderr.read(&mut byte) {
            Ok(0) => break, // EOF - vsearch closed stderr, it's done
            Ok(_) => {
                if byte[0] == b'\r' || byte[0] == b'\n' {
                    if !buf.is_empty() {
                        handle_vsearch_line(&String::from_utf8_lossy(&buf), pb);
                        buf.clear();
                    }
                } else {
                    buf.push(byte[0]);
                }
            }
            Err(_) => break, // pipe closed or read error - stop, don't fail the whole run over a progress-parsing hiccup
        }
    }

    if !buf.is_empty() {
        handle_vsearch_line(&String::from_utf8_lossy(&buf), pb);
    }
}

fn handle_vsearch_line(line: &str, pb: &ProgressBar) {
    let line = line.trim();
    if line.is_empty() {
        return;
    }

    match extract_trailing_percent(line) {
        Some(pct) => {
            pb.set_position(pct as u64);
            let label = line
                .trim_end_matches(|c: char| c.is_ascii_digit() || c == '%')
                .trim();
            if !label.is_empty() {
                pb.set_message(label.to_string());
            }
        }
        None => {
            // Non-percentage status lines (e.g. summary counts once a
            // phase completes) - still worth showing as the current
            // message, just without moving the bar's position.
            pb.set_message(line.to_string());
        }
    }
}

/// Parses a trailing "NN%" token off the end of a line, e.g. "Clustering
/// 100%" -> Some(100). Returns None if the line doesn't end in a
/// percentage (vsearch prints plenty of non-progress status lines too).
fn extract_trailing_percent(line: &str) -> Option<u32> {
    let line = line.trim_end();
    if !line.ends_with('%') {
        return None;
    }
    let without_pct = &line[..line.len() - 1];
    let digits: String = without_pct
        .chars()
        .rev()
        .take_while(|c| c.is_ascii_digit())
        .collect();
    if digits.is_empty() {
        return None;
    }
    let digits: String = digits.chars().rev().collect();
    digits.parse::<u32>().ok()
}

// ── UC-format cluster table parser ───────────────────────────────────────────
//
// UC format tab-separated columns:
//   0: record type  H=hit, S=centroid, C=cluster summary, N=no-match
//   1: cluster number (0-based)
//   8: query sequence label
//   9: target (centroid) label, or "*" for centroids
//
// Not currently called anywhere in the crate - kept for future consumers
// that need per-cluster membership (e.g. reporting which sequence was the
// centroid vs. a hit). Silenced rather than deleted since it's a complete,
// tested, working parser; revisit if/when something needs this data.

#[derive(Debug, Clone)]
#[allow(dead_code)]
pub struct ClusterMember {
    pub cluster_id: usize,
    pub seq_id: String,
    pub is_centroid: bool,
}

/// Parse a .clust UC-format file into per-cluster membership lists.
#[allow(dead_code)]
pub fn parse_uc_file(path: &Path) -> Result<Vec<Vec<String>>> {
    let content = std::fs::read_to_string(path)
        .with_context(|| format!("Cannot read UC file {:?}", path))?;

    let mut clusters: Vec<Vec<String>> = Vec::new();

    for line in content.lines() {
        if line.is_empty() || line.starts_with('#') { continue; }
        let cols: Vec<&str> = line.trim().split('\t').collect();
        if cols.len() < 10 { continue; }
        let record_type = cols[0].trim();

        if record_type == "C" {
            continue;
        }

        let cluster_idx: usize = cols[1].parse().unwrap_or(0);
        let seq_id = cols[8].trim().to_string();

        while clusters.len() <= cluster_idx {
            clusters.push(Vec::new());
        }
        clusters[cluster_idx].push(seq_id);
    }
    Ok(clusters)
}

// ── Validation ────────────────────────────────────────────────────────────────

fn validate_identity(id: f64) -> Result<()> {
    if id <= 0.0 || id > 1.0 {
        bail!("Sequence identity must be in the range (0.0, 1.0]. Got: {}", id);
    }
    Ok(())
}

// ── Tests ─────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_validate_identity_ok() {
        assert!(validate_identity(1.0).is_ok());
        assert!(validate_identity(0.95).is_ok());
        assert!(validate_identity(0.01).is_ok());
    }

    #[test]
    fn test_validate_identity_bad() {
        assert!(validate_identity(0.0).is_err());
        assert!(validate_identity(1.1).is_err());
        assert!(validate_identity(-0.5).is_err());
    }

    #[test]
    fn test_parse_uc_file() {
        use std::io::Write;
        let uc = "S	0	500	*	*	*	*	*	seq1	*
                  H	0	450	99.0	+	0	0	*	seq2	seq1
                  H	0	480	97.0	+	0	0	*	seq3	seq1
                  S	1	600	*	*	*	*	*	seq4	*
                  C	0	3	*	*	*	*	*	seq1	*
                  C	1	1	*	*	*	*	*	seq4	*
";
        let tmp = tempfile::NamedTempFile::new().unwrap();
        tmp.as_file().write_all(uc.as_bytes()).unwrap();
        let clusters = parse_uc_file(tmp.path()).unwrap();
        assert_eq!(clusters.len(), 2);
        assert_eq!(clusters[0].len(), 3);  // seq1, seq2, seq3
        assert_eq!(clusters[1].len(), 1);  // seq4
    }

    #[test]
    fn test_parse_uc_complex() {
        use std::io::Write;
        let uc = "S	0	500	*	*	*	*	*	seq1	*
                  H	0	480	99.5	+	0	0	=	seq2	seq1
                  H	0	460	98.2	+	0	0	=	seq3	seq1
                  S	1	300	*	*	*	*	*	seq4	*
                  S	2	200	*	*	*	*	*	seq5	*
                  H	2	195	97.0	+	0	0	=	seq6	seq5
                  C	0	3	*	*	*	*	*	seq1	*
                  C	1	1	*	*	*	*	*	seq4	*
                  C	2	2	*	*	*	*	*	seq5	*
";
        let tmp = tempfile::NamedTempFile::new().unwrap();
        tmp.as_file().write_all(uc.as_bytes()).unwrap();
        let clusters = parse_uc_file(tmp.path()).unwrap();
        assert_eq!(clusters.len(), 3);
        assert_eq!(clusters[0].len(), 3);
        assert_eq!(clusters[1].len(), 1);
        assert_eq!(clusters[2].len(), 2);
    }

    #[test]
    fn test_extract_trailing_percent_present() {
        assert_eq!(extract_trailing_percent("Clustering 100%"), Some(100));
        assert_eq!(extract_trailing_percent("Reading file 45%"), Some(45));
        assert_eq!(extract_trailing_percent("Sorting by length 0%"), Some(0));
    }

    #[test]
    fn test_extract_trailing_percent_absent() {
        assert_eq!(extract_trailing_percent("vsearch v2.31.0"), None);
        assert_eq!(extract_trailing_percent("1234567 nt in 5000 seqs"), None);
        assert_eq!(extract_trailing_percent(""), None);
    }

    #[test]
    fn test_stream_vsearch_progress_reads_carriage_return_updates() {
        // Simulates vsearch's own output style: several \r-updated
        // percentages within one phase, then a real \n before the next
        // phase - without byte-level \r splitting, all of "Reading file"
        // would buffer into one string and never surface intermediate ticks.
        let simulated = b"Reading file 0%\rReading file 50%\rReading file 100%\nClustering 0%\rClustering 100%\n";
        let pb = ProgressBar::hidden();
        stream_vsearch_progress(&simulated[..], &pb);
        // final state should reflect the last update seen (Clustering 100%)
        assert_eq!(pb.position(), 100);
    }
}
