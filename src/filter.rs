// src/filter.rs

use anyhow::{Context, Result};
use log::{debug, info};
use std::path::{Path, PathBuf};

use crate::fasta::{self, FastaRecord};
use crate::score::ScoreMap;

// ── Public API ────────────────────────────────────────────────────────────────

pub fn filter_assemblies(
    assembly_files: &[PathBuf],
    scores: &ScoreMap,
    min_score: f64,
) -> Result<Vec<PathBuf>> {
    if min_score <= 0.0 {
        return Ok(assembly_files.to_vec());
    }

    let mut kept = Vec::new();

    // NOTE: this must be called with the FULL, original assembly list (in
    // the same order used by score::score_assemblies' prefix_keys=true
    // pass), so that index 0 here means the same assembly as "contig0_..."
    // in the ScoreMap keys.
    for (idx, path) in assembly_files.iter().enumerate() {
        let contig_prefix = format!("contig{idx}_");

        let has_passing = scores.iter().any(|(id, score)| {
            id.starts_with(&contig_prefix)
                && score.score >= min_score
                && score.p_bases_covered >= 0.5
        });

        if has_passing {
            kept.push(path.clone());
            debug!("Keeping assembly {:?}", path);
        } else {
            info!(
                "Dropping assembly {:?} \
(no contigs above threshold {})",
                path,
                min_score
            );
        }
    }

    Ok(kept)
}

pub fn write_filtered_fasta(
    input_fasta: &Path,
    scores: &ScoreMap,
    min_score: f64,
    output: &Path,
) -> Result<PathBuf> {
    let records = fasta::load_fasta_ordered(input_fasta)?;
    let total_records = records.len();

    let kept: Vec<FastaRecord> = records
        .into_iter()
        .filter(|rec| {
            match scores.get(rec.id()) {
                Some(score) => {
                    let keep =
                        score.score >= min_score
                            && score.p_bases_covered >= 0.5;

                    if !keep {
                        debug!(
                            "Dropping {} \
(score={:.3}, p_bases_covered={:.3})",
                            rec.id(),
                            score.score,
                            score.p_bases_covered
                        );
                    }

                    keep
                }

                None => {
                    debug!(
                        "Missing score for {}",
                        rec.id()
                    );
                    false
                }
            }
        })
        .collect();

    info!(
        "Final filter: kept {}/{} contigs",
        kept.len(),
        total_records
    );

    fasta::write_fasta_file(output, &kept)
        .with_context(|| {
            format!(
                "Failed to write filtered FASTA to {:?}",
                output
            )
        })?;

    Ok(output.to_path_buf())
}

// ── Tests ─────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use tempfile::tempdir;
    use crate::score::ContigScore;

    #[test]
    fn test_filter_assemblies_keeps_all_when_threshold_zero() {
        let dir = tempdir().unwrap();

        let a = dir.path().join("a.fa");
        let b = dir.path().join("b.fa");

        std::fs::write(&a, ">s1\nACGT\n").unwrap();
        std::fs::write(&b, ">s1\nACGT\n").unwrap();

        let scores: ScoreMap = [
            (
                "contig0_s1".into(),
                ContigScore {
                    score: 0.1,
                    p_good: 0.1,
                    p_bases_covered: 0.1,
                    coverage: 0.5,
                }
            ),
            (
                "contig1_s1".into(),
                ContigScore {
                    score: 0.9,
                    p_good: 0.9,
                    p_bases_covered: 0.9,
                    coverage: 2.0,
                }
            ),
        ]
            .into_iter()
            .collect();

        let result =
            filter_assemblies(&[a, b], &scores, 0.0)
                .unwrap();

        assert_eq!(result.len(), 2);
    }

    #[test]
    fn test_filter_assemblies_drops_below_threshold() {
        let dir = tempdir().unwrap();

        let a = dir.path().join("a.fa");
        let b = dir.path().join("b.fa");

        std::fs::write(&a, ">s1\nACGT\n").unwrap();
        std::fs::write(&b, ">s1\nACGT\n").unwrap();

        let scores: ScoreMap = [
            (
                "contig0_s1".into(),
                ContigScore {
                    score: 0.1,
                    p_good: 0.1,
                    p_bases_covered: 0.1,
                    coverage: 0.5,
                }
            ),
            (
                "contig1_s1".into(),
                ContigScore {
                    score: 0.9,
                    p_good: 0.9,
                    p_bases_covered: 0.9,
                    coverage: 2.0,
                }
            ),
        ]
            .into_iter()
            .collect();

        let result =
            filter_assemblies(
                &[a, b.clone()],
                &scores,
                0.5,
            )
                .unwrap();

        assert_eq!(result.len(), 1);
        assert_eq!(result[0], b);
    }

    #[test]
    fn test_write_filtered_fasta() {
        let dir = tempdir().unwrap();

        let input = dir.path().join("input.fa");
        let output = dir.path().join("output.fa");

        std::fs::write(
            &input,
            ">good\nAAAA\n>bad\nCCCC\n",
        )
            .unwrap();

        let scores: ScoreMap = [
            (
                "good".into(),
                ContigScore {
                    score: 0.9,
                    p_good: 0.9,
                    p_bases_covered: 0.9,
                    coverage: 2.0,
                }
            ),
            (
                "bad".into(),
                ContigScore {
                    score: 0.1,
                    p_good: 0.1,
                    p_bases_covered: 0.1,
                    coverage: 0.5,
                }
            ),
        ]
            .into_iter()
            .collect();

        write_filtered_fasta(
            &input,
            &scores,
            0.5,
            &output,
        )
            .unwrap();

        let records =
            fasta::load_fasta_ordered(&output)
                .unwrap();

        assert_eq!(records.len(), 1);
        assert_eq!(records[0].id(), "good");
    }
}