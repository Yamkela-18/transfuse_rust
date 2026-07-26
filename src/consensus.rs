// src/consensus.rs
//
// Score-guided consensus selection — replaces Transfuse#consensus in transfuse.rb.
//
// Parse vsearch --msaout output, pick the highest-scoring member per cluster,
// write representatives to <output_stem>_cons.fa.
//
// MSA format: clusters separated by "//" lines; sequences in FASTA format
// with gap characters ('-') that we strip before writing.

use anyhow::{Context, Result};
use log::{debug, info, warn};
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::{Path, PathBuf};

use crate::fasta::{write_fasta_record, FastaRecord};
use crate::score::ScoreMap;

// ────────────────────────────────────────────────────────────────

pub fn build_consensus(
    msa_path: &Path,
    scores: &ScoreMap,
    sequences: &HashMap<String, FastaRecord>,
    output_path: &Path,
) -> Result<PathBuf> {
    let cons_path = consensus_path(output_path);

    let file = File::open(msa_path)
        .with_context(|| format!("Cannot open MSA file {:?}", msa_path))?;
    let reader = BufReader::new(file);

    let out_file = File::create(&cons_path)
        .with_context(|| format!("Cannot create consensus FASTA {:?}", cons_path))?;
    let mut writer = BufWriter::new(out_file);

    let mut total_clusters = 0usize;
    let mut singletons     = 0usize;
    let mut multi_member   = 0usize;

    for cluster in parse_msa_clusters(reader)? {
        total_clusters += 1;
        match pick_best(&cluster, scores, sequences) {
            Some(rec) => {
                write_fasta_record(&mut writer, &rec.header, &rec.sequence)?;
                if cluster.len() == 1 { singletons += 1; }
                else {
                    multi_member += 1;
                    debug!("  Cluster of {} -> selected {} (score={:.3})",
                        cluster.len(), rec.id(),
                        scores.get(rec.id()).map(|s| s.score).unwrap_or(0.0));
                }
            }
            None => warn!("  Empty cluster encountered, skipping."),
        }
    }

    writer.flush()?;
    info!("  Consensus: {} clusters ({} singletons, {} merged), written to {:?}",
          total_clusters, singletons, multi_member, cons_path);
    Ok(cons_path)
}

// ── Core selection logic ──────────────────────────────────────────────────────

fn pick_best<'a>(
    cluster: &'a [FastaRecord],
    scores: &ScoreMap,
    sequences: &'a HashMap<String, FastaRecord>,
) -> Option<&'a FastaRecord> {
    if cluster.is_empty() { return None; }

    let best = cluster.iter().max_by(|a, b| {
        let sa = scores
            .get(a.id())
            .map(|s| s.score)
            .unwrap_or(0.0);

        let sb = scores
            .get(b.id())
            .map(|s| s.score)
            .unwrap_or(0.0);

        sa.partial_cmp(&sb)
            .unwrap_or(std::cmp::Ordering::Equal)
    });

    // Prefer original (ungapped) sequence from loaded FASTA HashMap
    best.map(|rec| sequences.get(rec.id()).unwrap_or(rec))
}

// ── MSA parser ────────────────────────────────────────────────────────────────

fn parse_msa_clusters<R: BufRead>(reader: R) -> Result<Vec<Vec<FastaRecord>>> {
    let mut clusters: Vec<Vec<FastaRecord>> = Vec::new();
    let mut current_cluster: Vec<FastaRecord> = Vec::new();
    let mut current_header: Option<String> = None;
    let mut current_seq = String::new();

    for line in reader.lines() {
        let line = line?;
        let line = line.trim_end();

        if line == "//" {
            if let Some(hdr) = current_header.take() {
                current_cluster.push(FastaRecord {
                    header: hdr,
                    sequence: current_seq.replace('-', "").to_uppercase(),
                });
                current_seq = String::new();
            }
            if !current_cluster.is_empty() {
                clusters.push(current_cluster);
                current_cluster = Vec::new();
            }
            continue;
        }

        if let Some(header) = line.strip_prefix('>') {
            if let Some(hdr) = current_header.take() {
                current_cluster.push(FastaRecord {
                    header: hdr,
                    sequence: current_seq.replace('-', "").to_uppercase(),
                });
                current_seq = String::new();
            }
            current_header = Some(header.to_string());
        } else if !line.is_empty() {
            current_seq.push_str(line);
        }
    }

    // Flush final record and cluster (file may not end with //)
    if let Some(hdr) = current_header.take() {
        current_cluster.push(FastaRecord {
            header: hdr,
            sequence: current_seq.replace('-', "").to_uppercase(),
        });
    }
    if !current_cluster.is_empty() {
        clusters.push(current_cluster);
    }
    Ok(clusters)
}

// ── Path helper ───────────────────────────────────────────────────────────────

fn consensus_path(output_path: &Path) -> PathBuf {
    let parent = output_path.parent().unwrap_or(Path::new("."));
    let stem = output_path
        .file_stem()
        .map(|s| s.to_string_lossy().into_owned())
        .unwrap_or_else(|| "merged".into());
    parent.join(format!("{stem}_cons.fa"))
}

// ── Tests ─────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Cursor;
    use crate::score::ContigScore;

    fn make_seqs(ids: &[(&str, &str)]) -> HashMap<String, FastaRecord> {
        ids.iter().map(|(id, seq)| (id.to_string(), FastaRecord {
            header: id.to_string(), sequence: seq.to_string(),
        })).collect()
    }

    fn make_scores(pairs: &[(&str, f64)]) -> ScoreMap {
        pairs
            .iter()
            .map(|(id, s)| {
                (
                    id.to_string(),
                    ContigScore {
                        score: *s,
                        p_good: *s,
                        p_bases_covered: *s,
                        coverage: *s,
                    },
                )
            })
            .collect()
    }

    #[test]
    fn test_parse_single_cluster() {
        let msa = ">seq1
ACGT
>seq2
ACG-
//
";
        let clusters = parse_msa_clusters(BufReader::new(Cursor::new(msa))).unwrap();
        assert_eq!(clusters.len(), 1);
        assert_eq!(clusters[0].len(), 2);
        assert_eq!(clusters[0][1].sequence, "ACG");
    }

    #[test]
    fn test_parse_multiple_clusters() {
        let msa = ">s1
AAAA
>s2
AAA-
//
>s3
GGGG
//
";
        let clusters = parse_msa_clusters(BufReader::new(Cursor::new(msa))).unwrap();
        assert_eq!(clusters.len(), 2);
        assert_eq!(clusters[0].len(), 2);
        assert_eq!(clusters[1].len(), 1);
    }

    #[test]
    fn test_pick_best_by_score() {
        let cluster = vec![
            FastaRecord { header: "low".into(),  sequence: "AAAA".into() },
            FastaRecord { header: "high".into(), sequence: "CCCC".into() },
        ];
        let scores = make_scores(&[("low", 0.3), ("high", 0.9)]);
        let seqs   = make_seqs(&[("low", "AAAA"), ("high", "CCCC")]);
        assert_eq!(pick_best(&cluster, &scores, &seqs).unwrap().id(), "high");
    }

    #[test]
    fn test_pick_best_fallback_no_scores() {
        let cluster = vec![
            FastaRecord { header: "s1".into(), sequence: "AAAA".into() },
            FastaRecord { header: "s2".into(), sequence: "CCCC".into() },
        ];
        assert!(pick_best(&cluster, &HashMap::new(), &HashMap::new()).is_some());
    }

    #[test]
    fn test_build_consensus_full_pipeline() {
        use tempfile::tempdir;
        let dir = tempdir().unwrap();
        let msa_content = ">k31__t1
ACGTACGT
>k41__t1
ACGTACG-
//
>k31__t2
GGGGGGGG
//
";
        let msa = dir.path().join("combined-1.00.aln");
        std::fs::write(&msa, msa_content).unwrap();
        let scores = make_scores(&[("k31__t1",0.9),("k41__t1",0.6),("k31__t2",0.8)]);
        let seqs   = make_seqs(&[("k31__t1","ACGTACGT"),("k41__t1","ACGTACG"),("k31__t2","GGGGGGGG")]);
        let output = dir.path().join("merged.fa");
        let cons = build_consensus(&msa, &scores, &seqs, &output).unwrap();
        assert!(cons.exists());
        let records = crate::fasta::load_fasta_ordered(&cons).unwrap();
        assert_eq!(records.len(), 2);
        assert_eq!(records[0].id(), "k31__t1");
        assert_eq!(records[1].id(), "k31__t2");
    }

    #[test]
    fn test_consensus_path_derivation() {
        let out = Path::new("/tmp/merged.fa");
        assert_eq!(consensus_path(out), PathBuf::from("/tmp/merged_cons.fa"));
    }
}