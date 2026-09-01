// src/quant.rs
//
// Salmon-based effective-length-normalized coverage - replaces
// ReadMetrics#analyse_expression in the original Ruby Transfuse/TransRate.
//
// Ruby's coverage formula per contig:
//   coverage = eff_count * read_length / eff_len   (0 if eff_len == 0)
//
// `eff_len` (EffectiveLength) and `eff_count` (NumReads) come from Salmon's
// EM-based quantification, run in *alignment-based* mode against the same
// BAM already produced for the other read-mapping metrics - not a second,
// independent quasi-mapping pass. `read_length` is a single representative
// value estimated from the first few thousand reads, matching Ruby's
// ReadMetrics#get_read_length.
//
// This is intentionally a different metric from bam::ContigStats::mean_depth
// (raw samtools depth): the two measure genuinely different things, and
// samtools depth does not reproduce Salmon's effective-length behaviour -
// including its well-known blow-up for contigs shorter than the library's
// fragment-size distribution, which is a property of the *original* tool's
// metric, not a bug to route around.

use anyhow::{bail, Context, Result};
use std::collections::HashMap;
use std::fs;
use std::io::{BufRead, BufReader};
use std::path::Path;
use std::process::{Command, Stdio};

/// Number of FASTQ records sampled to estimate a single representative read
/// length, matching Ruby's ReadMetrics#get_read_length (reads the first
/// 5000 records, keeps the maximum sequence length seen).
const READ_LENGTH_SAMPLE: usize = 5000;

/// Opens a FASTQ file for reading, transparently decompressing it via
/// `gunzip` if it has a .gz extension. We only ever read a small prefix of
/// the file here (see estimate_read_length), so the child process is left
/// to receive SIGPIPE and exit once the reader is dropped - its exit status
/// isn't needed for this estimate.
fn open_fastq(path: &Path) -> Result<Box<dyn BufRead>> {
    if path.extension().and_then(|e| e.to_str()) == Some("gz") {
        let mut child = Command::new("sh")
            .arg("-c")
            .arg(format!("gunzip -c {:?}", path))
            .stdout(Stdio::piped())
            .stderr(Stdio::null())
            .spawn()
            .context("Failed to spawn gunzip")?;
        let stdout = child
            .stdout
            .take()
            .context("Failed to capture gunzip stdout")?;
        Ok(Box::new(BufReader::new(stdout)))
    } else {
        let file =
            fs::File::open(path).with_context(|| format!("Failed to open reads file {:?}", path))?;
        Ok(Box::new(BufReader::new(file)))
    }
}

/// Estimates a representative read length from the first ~5000 records of
/// a FASTQ file (plain or gzipped).
pub fn estimate_read_length(left: &Path) -> Result<u64> {
    let reader = open_fastq(left)?;
    let mut lines = reader.lines();

    let mut max_len: u64 = 0;
    let mut records = 0usize;

    while records < READ_LENGTH_SAMPLE {
        let _name = match lines.next() {
            Some(l) => l?,
            None => break,
        };
        let seq = match lines.next() {
            Some(l) => l?,
            None => break,
        };
        let _plus = match lines.next() {
            Some(l) => l?,
            None => break,
        };
        let _qual = match lines.next() {
            Some(l) => l?,
            None => break,
        };

        max_len = max_len.max(seq.len() as u64);
        records += 1;
    }

    if max_len == 0 {
        bail!(
            "Could not determine read length from {:?} (file empty or malformed)",
            left
        );
    }

    Ok(max_len)
}

/// Runs `samtools sort -n` to produce a name-grouped copy of the BAM for
/// salmon's alignment mode. Salmon's own CLI reference is explicit that
/// `-a/--alignments` needs "a name-grouped BAM of reads aligned to the
/// transcriptome" (mates adjacent), whereas bam.rs's samtools coverage/
/// mpileup passes want coordinate-sorted for positional access - the two
/// requirements conflict, so a separate throwaway copy is sorted here
/// rather than changing the BAM the rest of the pipeline depends on.
fn name_sort_for_salmon(bam_path: &Path, threads: usize) -> Result<std::path::PathBuf> {
    let out_path = bam_path.with_extension("namesorted.bam");

    let status = Command::new("samtools")
        .arg("sort")
        .arg("-n")
        .arg("-@")
        .arg(threads.to_string())
        .arg("-o")
        .arg(&out_path)
        .arg(bam_path)
        .stderr(Stdio::null()) // suppress samtools' own status lines
        .status()
        .context("Failed to run samtools sort -n for salmon input")?;

    if !status.success() {
        bail!("samtools sort -n exited with a non-zero status");
    }

    Ok(out_path)
}

/// Runs `salmon quant` in alignment-based mode against a name-grouped copy
/// of the BAM already produced for the other read-mapping metrics (the
/// content is identical, just re-ordered - not a second independent
/// alignment pass), then parses quant.sf for (EffectiveLength, NumReads)
/// per contig.
///
/// Uses `.output()` rather than `.status()`: Salmon writes its own verbose
/// progress/timing log directly to stderr regardless of verbosity flags
/// (the "INFO salmon::timing: phase complete ..." lines), which - since it
/// never goes through Rust's `log` crate - is invisible to the
/// MultiProgress-based coordination in main.rs and, left inherited, visibly
/// corrupts whichever progress bar is active. `.output()` captures it
/// instead of inheriting it; it's discarded on success and surfaced in the
/// error message if salmon actually fails, so nothing is lost for
/// debugging, just hidden on the happy path.
pub fn run_salmon_quant(
    bam_path: &Path,
    fasta_path: &Path,
    threads: usize,
) -> Result<HashMap<String, (f64, f64)>> {
    let name_sorted = name_sort_for_salmon(bam_path, threads)?;

    let out_dir = bam_path.with_extension("salmon_quant");
    if out_dir.exists() {
        let _ = fs::remove_dir_all(&out_dir);
    }

    let output = Command::new("salmon")
        .args(["quant", "--libType", "A", "--no-version-check"])
        .arg("--alignments")
        .arg(&name_sorted)
        .arg("--targets")
        .arg(fasta_path)
        .arg("--threads")
        .arg(threads.to_string())
        .arg("--output")
        .arg(&out_dir)
        .output()
        .context("Failed to run salmon quant")?;

    // The name-sorted copy is only ever needed for this one call - clean it
    // up regardless of whether salmon succeeded.
    let _ = fs::remove_file(&name_sorted);

    if !output.status.success() {
        let stderr = String::from_utf8_lossy(&output.stderr);
        bail!("salmon quant exited with a non-zero status: {}", stderr);
    }

    let quant_sf = out_dir.join("quant.sf");
    let file = fs::File::open(&quant_sf)
        .with_context(|| format!("Cannot open salmon output {:?}", quant_sf))?;
    let reader = BufReader::new(file);

    let mut out = HashMap::new();

    for (i, line) in reader.lines().enumerate() {
        let line = line?;
        if i == 0 || line.is_empty() {
            continue; // header: Name  Length  EffectiveLength  TPM  NumReads
        }
        let cols: Vec<&str> = line.split('\t').collect();
        if cols.len() < 5 {
            continue;
        }
        let name = cols[0].to_string();
        let eff_len: f64 = cols[2].parse().unwrap_or(0.0);
        let num_reads: f64 = cols[4].parse().unwrap_or(0.0);
        out.insert(name, (eff_len, num_reads));
    }

    Ok(out)
}

/// Ruby's coverage formula: `eff_count * read_length / eff_len`, or 0 if
/// eff_len is 0 (Salmon reports this for transcripts with essentially no
/// assignable read evidence).
pub fn effective_coverage(eff_len: f64, num_reads: f64, read_length: u64) -> f64 {
    if eff_len <= 0.0 {
        0.0
    } else {
        num_reads * read_length as f64 / eff_len
    }
}

// ── Tests ─────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_effective_coverage_normal_case() {
        // 100 reads of length 100 assigned to a transcript with effective
        // length 500 -> (100 * 100) / 500 = 20x.
        let cov = effective_coverage(500.0, 100.0, 100);
        assert!((cov - 20.0).abs() < 1e-9);
    }

    #[test]
    fn test_effective_coverage_zero_eff_len() {
        assert_eq!(effective_coverage(0.0, 50.0, 100), 0.0);
    }

    #[test]
    fn test_effective_coverage_short_contig_inflates() {
        // A very short effective length with even a modest read count
        // produces the same kind of blow-up the original tool exhibits -
        // this is expected behaviour we're matching, not a bug.
        let cov = effective_coverage(2.0, 300.0, 100);
        assert!(cov > 10_000.0);
    }
}
