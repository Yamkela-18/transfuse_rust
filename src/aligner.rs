// src/aligner.rs
//
// Read alignment — replaces SNAP (Scalable Nucleotide Alignment Program).
//
// Uses minimap2 (splice-aware, actively maintained, no glibc version pinning)
// to align paired-end FASTQ reads to a FASTA reference.
//
// minimap2 flags:
//   -ax sr          short-read preset (equivalent to SNAP's paired-end mode)
//   -t <threads>    parallelism
//   --secondary=no  suppress secondary alignments
//
// SAM output piped through:
//   samtools sort -@ <threads> -o <out.bam>
//   samtools index <out.bam>

use anyhow::{bail, Context, Result};
use log::{debug, info};
use std::path::{Path, PathBuf};
use std::process::{Command, Stdio};

// ── Public API ────────────────────────────────────────────────────────────────

/// Align paired-end reads to reference and return path to sorted BAM.
/// Idempotent: skips alignment if BAM already exists.
pub fn align_reads(
    reference: &Path,
    left: &Path,
    right: &Path,
    threads: usize,
    verbose: bool,
) -> Result<PathBuf> {
    let bam_path = bam_path_for(reference);

    if bam_path.exists() {
        info!("  BAM already exists, skipping alignment: {:?}", bam_path);
        return Ok(bam_path);
    }

    info!("  Aligning reads to {:?}", reference);
    info!("    minimap2 -ax sr -t {} --secondary=no | samtools sort", threads);

    // Spawn minimap2, pipe SAM stdout into samtools sort
    let minimap = Command::new("minimap2")
        .args(["-ax", "sr", "-t", &threads.to_string(), "--secondary=no"])
        .arg(reference)
        .arg(left)
        .arg(right)
        .stdout(Stdio::piped())
        .stderr(if verbose { Stdio::inherit() } else { Stdio::null() })
        .spawn()
        .context("Failed to launch minimap2. Is it installed and on PATH?")?;

    let minimap_stdout = minimap.stdout
        .context("Failed to capture minimap2 stdout")?;

    let sort_status = Command::new("samtools")
        .args(["sort", "-@", &threads.to_string(), "-o"])
        .arg(&bam_path)
        .arg("-")
        .stdin(minimap_stdout)
        .stderr(if verbose { Stdio::inherit() } else { Stdio::null() })
        .status()
        .context("Failed to launch samtools sort")?;

    if !sort_status.success() {
        bail!("samtools sort failed for {:?}", reference);
    }

    let index_status = Command::new("samtools")
        .args(["index", "-@", &threads.to_string()])
        .arg(&bam_path)
        .stderr(if verbose { Stdio::inherit() } else { Stdio::null() })
        .status()
        .context("Failed to launch samtools index")?;

    if !index_status.success() {
        bail!("samtools index failed for {:?}", bam_path);
    }

    debug!("  BAM written to {:?}", bam_path);
    Ok(bam_path)
}

// ── Helpers ───────────────────────────────────────────────────────────────────

/// Derive sorted BAM path from reference path.
/// e.g. k31.fa -> k31.sorted.bam
pub fn bam_path_for(reference: &Path) -> PathBuf {
    let stem = reference
        .file_stem()
        .map(|s| s.to_string_lossy().into_owned())
        .unwrap_or_else(|| "ref".into());
    let parent = reference.parent().unwrap_or(Path::new("."));
    parent.join(format!("{stem}.sorted.bam"))
}

// ── Tests ─────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_bam_path_derivation() {
        let ref_path = Path::new("/tmp/k31.fa");
        let bam = bam_path_for(ref_path);
        assert_eq!(bam, PathBuf::from("/tmp/k31.sorted.bam"));
    }
}
