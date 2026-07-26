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

use anyhow::{bail, Context, Result};
use log::{debug, info};
use std::path::{Path, PathBuf};
use std::process::Command;

// ── Public API ────────────────────────────────────────────────────────────────

/// Run vsearch clustering on fasta_path at identity threshold.
/// Idempotent: skips if .aln and .clust output files already exist.
pub fn cluster_vsearch(
    fasta_path: &Path,
    identity: f64,
    threads: usize,
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

    let status = Command::new("vsearch")
        .arg("--cluster_fast")
        .arg(fasta_path)
        .arg("--id").arg(&id_str)
        .arg("--msaout").arg(&aln_path)
        .arg("--uc").arg(&clust_path)
        .arg("--threads").arg(threads.to_string())
        .arg("--notrunclabels")
        .arg("--quiet")
        .status()
        .context("Failed to launch vsearch. Is it installed and on PATH?")?;

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
}