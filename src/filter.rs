// src/filter.rs
//
// Mirrors Ruby's Transfuse#filter and Transfuse#transrate_consensus's final
// write exactly: every individual contig is screened against the same
// fixed criteria,
//
//     score > 0.01  AND  coverage >= 1
//
// - not a whole-assembly keep/drop decision. Ruby's `filter` runs this once
// per input assembly *before* concatenation, writing each assembly's
// surviving contigs to a `<stem>_filtered.fa` with the contig's name
// rewritten to the same "contig{idx}_<name>" prefix used everywhere else,
// which is what lets a later plain `cat` concatenate them without name
// collisions. `transrate_consensus` runs the identical criteria once more,
// unprefixed, on the re-scored consensus assembly at the very end.
//
// Both entry points render a progress bar registered against the shared
// MultiProgress passed in from main.rs (so they coordinate with the
// overall pipeline bar and with any log output firing while they run):
// filter_assemblies shows assemblies completed, write_filtered_fasta shows
// contigs screened directly, since it already has the full record count
// up front.

use anyhow::{Context, Result};
use indicatif::{MultiProgress, ProgressBar, ProgressStyle};
use log::{debug, info, warn};
use std::path::{Path, PathBuf};

use crate::fasta::{self, FastaRecord};
use crate::score::{ContigScore, ScoreMap};

/// Ruby's filter criteria are fixed constants, not user-configurable:
/// `score > 0.01 and coverage >= 1`. `--min-score` here can raise the score
/// threshold above 0.01 (for stricter filtering) but never lower it or
/// skip the check entirely - matching the original's unconditional filter.
const MIN_SCORE_FLOOR: f64 = 0.01;
const MIN_COVERAGE: f64 = 1.0;

fn passes_filter(score: &ContigScore, min_score: f64) -> bool {
    let threshold = min_score.max(MIN_SCORE_FLOOR);
    score.score > threshold && score.coverage >= MIN_COVERAGE
}

fn filter_progress_bar(multi: &MultiProgress, len: u64, label: &str) -> ProgressBar {
    let pb = multi.add(ProgressBar::new(len));
    pb.set_style(
        ProgressStyle::with_template(&format!(
            "  {label} [{{bar:30.yellow/white}}] {{pos}}/{{len}} {{msg}} ({{percent}}%)"
        ))
        .unwrap()
        .progress_chars("=>-"),
    );
    pb
}

// ── Public API ────────────────────────────────────────────────────────────────

/// Per-contig filtering, run once per input assembly, before pooling and
/// clustering. For each assembly, writes only the contigs that individually
/// pass `score > 0.01 and coverage >= 1` to a new `<stem>_filtered.fa`,
/// with each surviving contig's header rewritten to
/// "contig{idx}_<original_id>" - matching Ruby's `filter` exactly, so a
/// later plain concatenation needs no further renaming.
///
/// Always returns exactly one path per input assembly, in the same order,
/// even if an assembly ends up with zero surviving contigs (matching
/// Ruby's `filtered_files << ...` running unconditionally for every input
/// file) - callers shouldn't assume the returned list means "this many
/// assemblies survived".
pub fn filter_assemblies(
    assembly_files: &[PathBuf],
    scores: &ScoreMap,
    min_score: f64,
    multi: &MultiProgress,
) -> Result<Vec<PathBuf>> {
    let mut filtered_files = Vec::with_capacity(assembly_files.len());

    let pb = filter_progress_bar(multi, assembly_files.len() as u64, "filtering");

    for (idx, file) in assembly_files.iter().enumerate() {
        pb.set_message(
            file.file_name()
                .map(|n| n.to_string_lossy().into_owned())
                .unwrap_or_else(|| file.display().to_string()),
        );

        let contig_prefix = format!("contig{idx}_");
        let stem = file
            .file_stem()
            .map(|s| s.to_string_lossy().into_owned())
            .unwrap_or_else(|| format!("assembly{idx}"));
        let new_path = file.with_file_name(format!("{stem}_filtered.fa"));

        let records = fasta::load_fasta_ordered(file)?;
        let total = records.len();

        let mut kept = Vec::new();
        for rec in records {
            let key = format!("{contig_prefix}{}", rec.id());

            let score = scores
                .get(&key)
                .with_context(|| format!("Can't find '{key}' in scores"))?;

            if passes_filter(score, min_score) {
                kept.push(FastaRecord {
                    header: key,
                    sequence: rec.sequence,
                });
            } else {
                debug!(
                    "  Dropping {key} (score={:.3}, coverage={:.3})",
                    score.score, score.coverage
                );
            }
        }

        info!("  filtering {:?}: kept {}/{} contigs", file, kept.len(), total);
        if kept.is_empty() {
            warn!("  {:?}: no contigs passed the filter - effectively dropped", file);
        }

        fasta::write_fasta_file(&new_path, &kept)
            .with_context(|| format!("Failed to write filtered FASTA to {:?}", new_path))?;

        filtered_files.push(new_path);
        pb.inc(1);
    }

    pb.finish_and_clear();
    multi.remove(&pb);

    Ok(filtered_files)
}

/// Final filter, run once against the re-scored consensus assembly at the
/// very end - mirrors Ruby's `transrate_consensus` final write: the same
/// fixed criteria as above, just against a single already-unprefixed
/// assembly rather than many.
pub fn write_filtered_fasta(
    input_fasta: &Path,
    scores: &ScoreMap,
    min_score: f64,
    output: &Path,
    multi: &MultiProgress,
) -> Result<PathBuf> {
    let records = fasta::load_fasta_ordered(input_fasta)?;
    let total_records = records.len();

    let pb = filter_progress_bar(multi, total_records as u64, "final filter");

    let kept: Vec<FastaRecord> = records
        .into_iter()
        .filter(|rec| {
            pb.inc(1);
            match scores.get(rec.id()) {
                Some(score) => {
                    let keep = passes_filter(score, min_score);
                    if !keep {
                        debug!(
                            "Dropping {} (score={:.3}, coverage={:.3})",
                            rec.id(),
                            score.score,
                            score.coverage
                        );
                    }
                    keep
                }
                None => {
                    debug!("Missing score for {}", rec.id());
                    false
                }
            }
        })
        .collect();

    pb.finish_and_clear();
    multi.remove(&pb);

    info!(
        "Final filter: kept {}/{} contigs",
        kept.len(),
        total_records
    );

    fasta::write_fasta_file(output, &kept)
        .with_context(|| format!("Failed to write filtered FASTA to {:?}", output))?;

    Ok(output.to_path_buf())
}

// ── Tests ─────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use tempfile::tempdir;

    fn score(score: f64, coverage: f64) -> ContigScore {
        ContigScore {
            score,
            p_good: score,
            p_bases_covered: score,
            coverage,
        }
    }

    #[test]
    fn test_filter_assemblies_screens_per_contig_not_per_file() {
        let dir = tempdir().unwrap();

        let a = dir.path().join("a.fa");
        let b = dir.path().join("b.fa");

        // Each assembly has one contig that should survive and one that
        // shouldn't - a whole-assembly filter would keep or drop the whole
        // file; a correct per-contig filter keeps exactly one from each.
        std::fs::write(&a, ">good_a\nACGT\n>bad_a\nTTTT\n").unwrap();
        std::fs::write(&b, ">good_b\nGGGG\n>bad_b\nCCCC\n").unwrap();

        let scores: ScoreMap = [
            ("contig0_good_a".into(), score(0.9, 2.0)),
            ("contig0_bad_a".into(), score(0.005, 2.0)), // fails score
            ("contig1_good_b".into(), score(0.9, 2.0)),
            ("contig1_bad_b".into(), score(0.9, 0.5)),   // fails coverage
        ]
        .into_iter()
        .collect();

        let multi = MultiProgress::new();
        let filtered = filter_assemblies(&[a, b], &scores, 0.0, &multi).unwrap();
        assert_eq!(filtered.len(), 2); // one output per input, always

        let recs_a = fasta::load_fasta_ordered(&filtered[0]).unwrap();
        assert_eq!(recs_a.len(), 1);
        assert_eq!(recs_a[0].id(), "contig0_good_a");

        let recs_b = fasta::load_fasta_ordered(&filtered[1]).unwrap();
        assert_eq!(recs_b.len(), 1);
        assert_eq!(recs_b[0].id(), "contig1_good_b");
    }

    #[test]
    fn test_filter_assemblies_empty_survivors_still_returns_a_path() {
        let dir = tempdir().unwrap();
        let a = dir.path().join("a.fa");
        std::fs::write(&a, ">only\nACGT\n").unwrap();

        let scores: ScoreMap = [("contig0_only".into(), score(0.005, 0.5))]
            .into_iter()
            .collect();

        let multi = MultiProgress::new();
        let filtered = filter_assemblies(&[a], &scores, 0.0, &multi).unwrap();
        assert_eq!(filtered.len(), 1);
        assert!(fasta::load_fasta_ordered(&filtered[0]).unwrap().is_empty());
    }

    #[test]
    fn test_filter_assemblies_errors_on_missing_score() {
        let dir = tempdir().unwrap();
        let a = dir.path().join("a.fa");
        std::fs::write(&a, ">mystery\nACGT\n").unwrap();

        let scores: ScoreMap = ScoreMap::new();
        let multi = MultiProgress::new();
        assert!(filter_assemblies(&[a], &scores, 0.0, &multi).is_err());
    }

    #[test]
    fn test_min_score_can_raise_but_not_lower_the_floor() {
        let dir = tempdir().unwrap();
        let a = dir.path().join("a.fa");
        std::fs::write(&a, ">c\nACGT\n").unwrap();

        // Passes the fixed 0.01 floor but not a stricter user threshold.
        let scores: ScoreMap = [("contig0_c".into(), score(0.05, 2.0))]
            .into_iter()
            .collect();

        let multi = MultiProgress::new();
        let lenient = filter_assemblies(&[a.clone()], &scores, 0.0, &multi).unwrap();
        assert_eq!(fasta::load_fasta_ordered(&lenient[0]).unwrap().len(), 1);

        let strict = filter_assemblies(&[a], &scores, 0.5, &multi).unwrap();
        assert_eq!(fasta::load_fasta_ordered(&strict[0]).unwrap().len(), 0);
    }

    #[test]
    fn test_write_filtered_fasta() {
        let dir = tempdir().unwrap();

        let input = dir.path().join("input.fa");
        let output = dir.path().join("output.fa");

        std::fs::write(&input, ">good\nAAAA\n>bad\nCCCC\n").unwrap();

        let scores: ScoreMap = [
            ("good".into(), score(0.9, 2.0)),
            ("bad".into(), score(0.1, 0.5)),
        ]
        .into_iter()
        .collect();

        let multi = MultiProgress::new();
        write_filtered_fasta(&input, &scores, 0.5, &output, &multi).unwrap();

        let records = fasta::load_fasta_ordered(&output).unwrap();
        assert_eq!(records.len(), 1);
        assert_eq!(records[0].id(), "good");
    }
}
