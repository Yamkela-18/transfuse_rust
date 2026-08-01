# Transfuse (Rust)

A Rust reimplementation of [Transfuse](https://github.com/cboursnell/transfuse) —
a tool that merges multiple *de novo* transcriptome assemblies into one evidence-scored, non-redundant
transcriptome.

If you assembled the same RNA-seq data with several assemblers or parameter
sets and don't want to just pick one and throw the rest away, this is the
"combine them properly" step.

## What it does

1. **Score** every contig in every input assembly against the original reads
   — no reference genome needed.
2. **Filter** out contigs/assemblies that score poorly.
3. **Pool** everything that's left into one candidate FASTA.
4. **Cluster** the pool by sequence identity (VSEARCH) to find contigs
   different assemblers produced for the same underlying transcript.
5. **Build a consensus** per cluster, weighted by which member scored
   highest.
6. **Re-score and filter** the consensus assembly once more.
7. Write the final merged FASTA.

## Requirements

Binary dependencies (checked via `--install`, not bundled):

| Tool | Used for |
|---|---|
| `samtools` (≥ 1.12) | coverage, alignment-detail, and mpileup passes |
| `salmon` | effective-length-normalized coverage quantification |
| `vsearch` | clustering and per-cluster multiple sequence alignment |
| an aligner (see `aligner.rs`) | producing the BAM everything else reads |

Build with:

```bash
cargo build --release
```

## Usage

Score every contig across both assemblies and write `new_scores.csv`
without building a merged output yet:

```bash
cargo run -- --assemblies trinity_output.fasta,rnaspades_output.fasta --left reads.left.fq.gz --right reads.right.fq.gz --output new.csv --score-only
```

Run the full pipeline end-to-end and write the merged, clustered,
consensus-built assembly:

```bash
cargo run -- --assemblies trinity_output.fasta,rnaspades_output.fasta --left reads.left.fq.gz --right reads.right.fq.gz --output new.fasta
```

Both `--left`/`--right` accept gzipped FASTQ (`.fq.gz`) directly — handled
transparently in `quant.rs`'s read-length sampling and by whichever aligner
`aligner.rs` wraps.

| Flag | Description |
|---|---|
| `-a, --assemblies` | Comma-separated assembly FASTA files |
| `-l, --left` / `-r, --right` | Paired-end reads used for the original assemblies |
| `-o, --output` | Path for the merged output FASTA |
| `-t, --threads` | Thread count (default: 1) |
| `-i, --id` | VSEARCH clustering identity threshold (default: 1.0) |
| `-s, --min-score` | Minimum contig score kept after scoring (default: 0.0) |
| `--scores <CSV_FILES>` | Load pre-computed scores instead of re-running the scoring pass |
| `--score-only` | Score everything, write `<output>_scores.csv`, then exit |
| `--install` | Check/install missing binary dependencies |
| `-v, --verbose` | Verbose logging |

A typical two-pass workflow, so re-runs with a different `--min-score` don't
re-score from scratch:

```bash
cargo run -- --assemblies trinity_output.fasta,rnaspades_output.fasta --left reads.left.fq.gz --right reads.right.fq.gz --output new.fasta --score-only

cargo run -- --assemblies trinity_output.fasta,rnaspades_output.fasta --left reads.left.fq.gz --right reads.right.fq.gz --output new.fasta --scores new_scores.csv --min-score 0.05
```

## Scoring methodology

Each contig gets four independent, reference-free proportions in `[0, 1]`,
combined as a straight product — matching TransRate's own integration, so a
single badly-failed check can't get averaged away by three good ones:

```
score = sCnuc × sCcov × sCord × sCseg      (floored at 0.01)
```

| Component | What it measures | How it's computed here |
|---|---|---|
| **sCnuc** | Per-base read-to-contig identity | Quality-weighted match/mismatch voting per position, via `samtools mpileup` |
| **sCcov** | Breadth of coverage | Fraction of contig length with ≥1 read aligned, via `samtools coverage` |
| **sCord** | Read-pair concordance | Fraction of pairs with opposite-strand mates and an insert size consistent with the library's own empirical fragment-size distribution (not the aligner's own "properly paired" flag) |
| **sCseg** | Coverage homogeneity (chimera detection) | A Chow test comparing a one-segment vs. best two-segment fit across 10 coverage windows |

Two more fields round out the output, matching the original tool's CSV shape:

- **`p_good`** — the same per-fragment pass rate as `sCord` above; this is
  TransRate's own distinct metric (one of the four multiplicative factors),
  *not* a copy of `score`.
- **`coverage`** — Salmon's effective-length-normalized estimate
  (`num_reads × read_length / effective_length`), run in alignment-based
  mode against the same BAM used everywhere else. This is deliberately not
  raw sequencing depth — it's what the original tool reports, artifacts
  (including its known blow-up on contigs shorter than the library's
  fragment size) included.

CSV output (`<output>_scores.csv`) mirrors the original's number formatting:
`score`/`p_good`/`p_bases_covered` are printed at full, unrounded precision;
`coverage` is rounded to 2 decimal places — both matching Ruby's
`Float#to_s` behavior rather than a fixed decimal count.

## Known differences from the original tool

This was built by reverse-engineering the original Ruby `transfuse` +
`transrate` gem's behavior and validating against its actual output, not
from its (unavailable to us) internal source in every case. Where it's a
deliberate simplification rather than a bug:

- **sCord** checks orientation + insert size only. The original's `bam-read`
  C tool most likely also folds in per-read mismatch/mapping-quality
  checks, which aren't included here.
- **sCseg** uses a Chow test on 10 coverage windows. The original is
  described as a Dirichlet-based single-vs-multi-distribution test — same
  goal, coarser statistical machinery.
- **No EM-based multi-mapping reassignment.** The original runs Salmon's
  quasi-mapping EM to probabilistically reassign ambiguous reads across
  contigs before scoring. This pipeline scores whatever the aligner already
  picked as each read's primary alignment. Expect the biggest divergence
  from the original on contigs sharing sequence with other isoforms or
  duplicated genes.
- **No assembly-level score.** The original also computes a read-weighted
  geometric-mean assembly score and searches for an "optimal" cutoff
  (`assembly_score` / `assembly_optimal_score` in the Ruby source). Nothing
  here currently aggregates above the per-contig level.

`filter.rs` screens every individual contig against the same fixed
criteria as the original (`score > 0.01 and coverage >= 1`), both before
clustering and again on the final consensus assembly — not a
whole-assembly keep/drop decision.

Despite the above, validation against real reference output showed strong
agreement on the classification that matters most in practice: which
contigs are clearly good vs. clearly bad (159/162 exact floor agreement,
`p_good` correlation ≈ 0.97). Treat absolute score values as a good-faith
approximation, not a guaranteed numerical match.

## Module layout

| File | Responsibility |
|---|---|
| `main.rs` | CLI, dependency checks, pipeline orchestration |
| `aligner.rs` | Runs the read aligner, produces the BAM |
| `bam.rs` | `samtools`-based per-contig alignment stats (coverage, pair concordance, depth windows, per-base identity) |
| `quant.rs` | Salmon-based effective-length coverage |
| `score.rs` | Combines the above into the four TransRate-style components and the final score |
| `filter.rs` | Drops low-scoring assemblies/contigs |
| `fasta.rs` | FASTA I/O, concatenation |
| `cluster.rs` | VSEARCH clustering |
| `consensus.rs` | Score-guided consensus selection per cluster |
| `deps.rs` | Binary dependency checking/installation |

## Testing

```bash
cargo test
```