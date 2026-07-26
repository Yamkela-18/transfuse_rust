// src/main.rs
//
// Transfuse — Rust implementation
// Equivalent to bin/transfuse (Ruby entry point) + lib/transfuse/transfuse.rb
//
// Pipeline:
//   1. Validate CLI args & files
//   2. Check / install binary dependencies (vsearch, minimap2, samtools, salmon)
//   3. Score every assembly with Salmon (quasi-mapping) + samtools coverage
//   4. Filter low-scoring contigs
//   5. Concatenate passing contigs
//   6. Cluster with vsearch --cluster_fast
//   7. Score-guided consensus selection
//   8. Final scoring pass -> write output FASTA

use anyhow::{bail, Context, Result};
use clap::Parser;
use colored::Colorize;
use log::info;
use std::path::PathBuf;

mod aligner;
mod bam;
mod cluster;
mod consensus;
mod deps;
mod fasta;
mod filter;
mod score;
mod quant;


pub const VERSION: &str = env!("CARGO_PKG_VERSION");

// ── CLI definition (replaces Ruby Trollop options) ────────────────────────────

/// Transfuse: intelligently merge multiple de novo transcriptome assemblies.
#[derive(Parser, Debug)]
#[command(
    name    = "transfuse",
    version = VERSION,
    about   = "Merge multiple de novo transcriptome assemblies",
    after_help = "EXAMPLE:
  transfuse -a k31.fa,k41.fa,k51.fa -l reads_1.fq -r reads_2.fq -o merged.fa -t 12"
)]
pub struct Cli {
    /// Assembly FASTA files, comma-separated
    #[arg(short = 'a', long, value_name = "FILES")]
    pub assemblies: Option<String>,

    /// Left (R1) reads in FASTQ format
    #[arg(short = 'l', long, value_name = "FILE")]
    pub left: Option<String>,

    /// Right (R2) reads in FASTQ format
    #[arg(short = 'r', long, value_name = "FILE")]
    pub right: Option<String>,

    /// Path for the merged output FASTA
    #[arg(short = 'o', long, value_name = "FILE")]
    pub output: Option<String>,

    /// Number of threads [default: 1]
    #[arg(short = 't', long, default_value_t = 1, value_name = "INT")]
    pub threads: usize,

    /// Sequence identity threshold for vsearch clustering [default: 1.0]
    #[arg(short = 'i', long, default_value_t = 1.0, value_name = "FLOAT")]
    pub id: f64,

    /// Minimum contig score to keep after scoring [default: 0.0]
    #[arg(short = 's', long, default_value_t = 0.0, value_name = "FLOAT")]
    pub min_score: f64,

    /// Load pre-computed scores from CSV files instead of running Salmon
    #[arg(long, value_name = "CSV_FILES")]
    pub scores: Option<String>,

    /// Score all assemblies, write <output_stem>_scores.csv, then exit.
    /// Use the output CSV with --scores on a subsequent run to skip re-scoring.
    #[arg(long)]
    pub score_only: bool,

    /// Check and install missing binary dependencies, then exit
    #[arg(long)]
    pub install: bool,

    /// Verbose output
    #[arg(short = 'v', long)]
    pub verbose: bool,
}

// ── Entry point ───────────────────────────────────────────────────────────────

fn main() {
    if let Err(e) = run() {
        eprintln!("{} {}", "error:".red().bold(), e);
        for cause in e.chain().skip(1) {
            eprintln!("  {} {}", "caused by:".yellow(), cause);
        }
        std::process::exit(1);
    }
}

fn run() -> Result<()> {
    let cli = Cli::parse();

    let log_level = if cli.verbose { "debug" } else { "info" };
    env_logger::Builder::new()
        .filter_level(log_level.parse().unwrap())
        .format_timestamp(None)
        .format_target(false)
        .init();

    print_banner();

    // 1. Dependency check / install
    info!("{}", "Checking binary dependencies...".cyan());
    let missing = deps::check_dependencies()?;

    if cli.install {
        if missing.is_empty() {
            println!("{}", "All dependencies already installed.".green());
        } else {
            deps::install_dependencies(&missing)?;
            println!("{}", "All dependencies installed successfully.".green());
        }
        return Ok(());
    }

    if !missing.is_empty() {
        let list: Vec<String> = missing.iter()
            .map(|d| format!("  {} >= {}", d.name.yellow(), d.version))
            .collect();
        bail!(
            "Missing binary dependencies:
{}
Run with --install to download them.",
            list.join("
")
        );
    }
    info!("{}", "All dependencies present.".green());

    // 2. File validation
    let assembly_files = require_files(
        cli.assemblies.as_deref(), "assemblies", "--assemblies / -a",
    )?;

    let output_path = cli.output.as_deref()
        .context("Please specify an output file with --output / -o")?;
    let output = PathBuf::from(output_path);
    if output.exists() {
        bail!("Output file {:?} already exists.", output);
    }

    // 3. Score every contig
    let scores = if let Some(csv_list) = cli.scores.as_deref() {
        info!("{}", "Loading scores from CSV files...".cyan());
        score::load_scores_from_csv(&parse_comma_list(csv_list))?
    } else {
        let left  = require_file(cli.left.as_deref(),  "--left / -l")?;
        let right = require_file(cli.right.as_deref(), "--right / -r")?;
        info!("{}", "Scoring assemblies with Salmon + samtools...".cyan());
        // Multiple raw assemblies can share contig names (e.g. "TRINITY_DN5..."
        // appearing in more than one assembler's output), so keys must be
        // prefixed with the assembly file stem to stay unique.
        score::score_assemblies(&assembly_files, &left, &right, cli.threads, cli.verbose, true)?
    };
    info!("Scored {} contigs across all assemblies.", scores.len());
    if cli.score_only {
        let csv_path = score::write_scores_csv(&scores, &output)?;
        println!("\n{} Scores written to {:?}", "✓".green().bold(), csv_path);
        println!("  {} contigs scored", scores.len().to_string().yellow());
        println!("\n  Re-run with: --scores {:?} --min-score <threshold> -o {:?}", csv_path, output);
        return Ok(());
    }

    // 4. Filter low-scoring contigs
    info!("{}", "Filtering contigs by score...".cyan());
    let filtered = filter::filter_assemblies(&assembly_files, &scores, cli.min_score)?;
    info!("Kept {}/{} assemblies after filtering.", filtered.len(), assembly_files.len());
    if filtered.is_empty() {
        bail!("All assemblies removed by score filter. Try lowering --min-score.");
    }

    // 5. Concatenate
    info!("{}", "Concatenating filtered assemblies...".cyan());
    let cat_path = fasta::concatenate_assemblies(&filtered, &assembly_files, &output)?;
    info!("Combined FASTA written to {:?}", cat_path);

    // 6. Load sequences into memory
    info!("{}", "Loading FASTA sequences...".cyan());
    let sequences = fasta::load_fasta(&cat_path)?;
    info!("Loaded {} sequences.", sequences.len());

    // 7. Cluster with vsearch
    info!("{}", "Clustering with vsearch...".cyan());
    let msa_path = cluster::cluster_vsearch(&cat_path, cli.id, cli.threads)?;
    info!("vsearch MSA written to {:?}", msa_path);

    // 8. Score-guided consensus selection
    info!("{}", "Selecting best representative per cluster...".cyan());
    let cons_path = consensus::build_consensus(&msa_path, &scores, &sequences, &output)?;
    info!("Consensus FASTA written to {:?}", cons_path);

    // 9. Final scoring pass
    let final_output = if let Some(left) = cli.left.as_deref() {
        if let Some(right) = cli.right.as_deref() {
            info!("{}", "Running final scoring pass on consensus...".cyan());
            // The consensus FASTA is a single file whose headers were already
            // made globally unique during the earlier prefixed scoring pass
            // (e.g. "trinity_output__TRINITY_DN5_..."), so scoring it again
            // must NOT re-prefix with the file's own stem ("meg_cons"),
            // or every key will fail to match rec.id() during filtering.
            let final_scores = score::score_assemblies(
                &[cons_path.clone()],
                &PathBuf::from(left),
                &PathBuf::from(right),
                cli.threads, cli.verbose,
                false,
            )?;
            filter::write_filtered_fasta(&cons_path, &final_scores, cli.min_score, &output)?
        } else {
            std::fs::copy(&cons_path, &output)
                .with_context(|| format!("Failed to copy {:?} to {:?}", cons_path, output))?;
            output.clone()
        }
    } else {
        std::fs::copy(&cons_path, &output)
            .with_context(|| format!("Failed to copy {:?} to {:?}", cons_path, output))?;
        output.clone()
    };

    // Summary
    let final_seqs = fasta::load_fasta(&final_output)?;
    println!("
{} Merged assembly written to {:?}", "✓".green().bold(), final_output);
    println!("  {} total contigs in merged assembly", final_seqs.len().to_string().yellow());
    Ok(())
}

// ── Helpers ───────────────────────────────────────────────────────────────────

fn print_banner() {
    println!("
{} {}
{}
",
             "transfuse".cyan().bold(),
             format!("v{VERSION}").dimmed(),
             "Merge multiple de novo transcriptome assemblies".dimmed()
    );
}

pub fn parse_comma_list(s: &str) -> Vec<PathBuf> {
    s.split(',').map(|p| PathBuf::from(p.trim())).collect()
}

fn require_files(value: Option<&str>, label: &str, flag: &str) -> Result<Vec<PathBuf>> {
    let s = value.with_context(|| format!("Please provide {label} with {flag}"))?;
    let paths = parse_comma_list(s);
    for p in &paths {
        if !p.exists() { bail!("{} file not found: {:?}", label, p); }
    }
    Ok(paths)
}

fn require_file(value: Option<&str>, flag: &str) -> Result<PathBuf> {
    let s = value.with_context(|| format!("Please provide a file with {flag}"))?;
    let p = PathBuf::from(s);
    if !p.exists() { bail!("File not found: {:?}", p); }
    Ok(p)
}