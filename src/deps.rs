// src/deps.rs
//
// Binary dependency management — replaces the Ruby bindeps gem.
//
// Tools managed:
//   vsearch   — sequence clustering
//   minimap2  — read alignment (replaces SNAP)
//   samtools  — BAM processing + coverage (replaces bam-read)
//   salmon    — RNA-seq quasi-mapping quantification (modern 1.x+)

use anyhow::{bail, Context, Result};
use colored::Colorize;
use log::{debug, info, warn};
use std::path::{Path, PathBuf};
use std::process::Command;

#[derive(Debug, Clone)]
pub struct BinaryDep {
    pub name: &'static str,
    pub version: &'static str,
    pub version_arg: &'static str,
    pub version_pattern: &'static str,
    pub download_url_linux: &'static str,
    pub download_url_macos: &'static str,
}

pub const REQUIRED_DEPS: &[BinaryDep] = &[
    BinaryDep {
        name: "vsearch", version: "2.0",
        version_arg: "--version", version_pattern: "vsearch v",
        download_url_linux: "https://github.com/torognes/vsearch/releases/download/v2.28.1/vsearch-2.28.1-linux-x86_64.tar.gz",
        download_url_macos: "https://github.com/torognes/vsearch/releases/download/v2.28.1/vsearch-2.28.1-macos-aarch64.tar.gz",
    },
    BinaryDep {
        name: "minimap2", version: "2.0",
        version_arg: "--version", version_pattern: "",
        download_url_linux: "https://github.com/lh3/minimap2/releases/download/v2.28/minimap2-2.28_x64-linux.tar.bz2",
        download_url_macos: "",  // brew install minimap2
    },
    BinaryDep {
        name: "samtools", version: "1.0",
        version_arg: "--version", version_pattern: "samtools ",
        download_url_linux: "",  // apt / conda / brew
        download_url_macos: "",
    },
    BinaryDep {
        name: "salmon", version: "1.0",
        version_arg: "--version", version_pattern: "salmon ",
        download_url_linux: "https://github.com/COMBINE-lab/salmon/releases/download/v1.10.3/salmon-1.10.3_linux_x86_64.tar.gz",
        download_url_macos: "https://github.com/COMBINE-lab/salmon/releases/download/v1.10.3/salmon-1.10.3_mac_osx.tar.gz",
    },
];

// ── Public API ────────────────────────────────────────────────────────────────

pub fn check_dependencies() -> Result<Vec<BinaryDep>> {
    let mut missing = Vec::new();
    for dep in REQUIRED_DEPS {
        match probe_binary(dep) {
            Ok(ver)  => info!("  {} {} ✓  (found {})", dep.name.green(), dep.version, ver.dimmed()),
            Err(e)   => { warn!("  {} {} ✗  ({})", dep.name.red(), dep.version, e);
                missing.push(dep.clone()); }
        }
    }
    Ok(missing)
}

pub fn install_dependencies(missing: &[BinaryDep]) -> Result<()> {
    let install_dir = default_install_dir();
    std::fs::create_dir_all(&install_dir)
        .with_context(|| format!("Cannot create install directory {:?}", install_dir))?;
    info!("Installing to {:?}", install_dir);

    for dep in missing {
        let url = if cfg!(target_os = "macos") { dep.download_url_macos }
        else { dep.download_url_linux };
        if url.is_empty() {
            warn!("No automatic download for '{}'. Install via apt, brew, or conda.", dep.name);
            continue;
        }
        info!("Downloading {} from {}", dep.name, url);
        download_and_install(dep.name, url, &install_dir)?;
        info!("{} installed to {:?}", dep.name.green(), install_dir);
    }
    println!("
Add {:?} to your PATH if not already included.", install_dir);
    Ok(())
}

// ── Internals ─────────────────────────────────────────────────────────────────

/// A lenient (major, minor, patch) version, parsed from tool version output
/// (e.g. "2.3.4", "1.19.2", "2.28.1_linux_x86_64" -> (2,28,1)). Missing
/// trailing components default to 0. Ord is derived, so tuple comparison
/// gives correct version ordering (major first, then minor, then patch).
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
struct SemVer(u32, u32, u32);

impl SemVer {
    /// Scans whitespace/punctuation-delimited tokens for the first one that
    /// starts with a digit, then parses as many of major/minor/patch as are
    /// present from it (stopping at the first non-digit, non-'.' char, so
    /// "2.28.1_linux_x86_64" correctly yields "2.28.1").
    fn parse(text: &str) -> Option<Self> {
        for token in text.split(|c: char| c.is_whitespace() || c == ',' || c == '(' || c == ')') {
            let numeric: String = token
                .chars()
                .take_while(|c| c.is_ascii_digit() || *c == '.')
                .collect();
            if !numeric.starts_with(|c: char| c.is_ascii_digit()) {
                continue;
            }
            let mut parts = numeric.split('.');
            let major: u32 = match parts.next().and_then(|s| s.parse().ok()) {
                Some(m) => m,
                None => continue,
            };
            let minor = parts.next().and_then(|s| s.parse().ok()).unwrap_or(0);
            let patch = parts.next().and_then(|s| s.parse().ok()).unwrap_or(0);
            return Some(SemVer(major, minor, patch));
        }
        None
    }
}

impl std::fmt::Display for SemVer {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}.{}.{}", self.0, self.1, self.2)
    }
}

/// Runs `dep`'s version command, confirms the output actually looks like
/// this tool (not just some unrelated binary sharing its name on PATH), and
/// - critically - parses the real version number and checks it against
/// `dep.version` as an actual minimum, not just a display label.
fn probe_binary(dep: &BinaryDep) -> Result<String> {
    let output = Command::new(dep.name)
        .arg(dep.version_arg)
        .output()
        .with_context(|| format!("'{}' not found on PATH", dep.name))?;
    let stdout = String::from_utf8_lossy(&output.stdout);
    let stderr = String::from_utf8_lossy(&output.stderr);
    let combined = format!("{}{}", stdout, stderr);

    if !dep.version_pattern.is_empty() && !combined.contains(dep.version_pattern) {
        bail!(
            "found a '{}' on PATH, but its version output didn't match the expected pattern '{}': {}",
            dep.name,
            dep.version_pattern,
            combined.lines().next().unwrap_or("").trim()
        );
    }

    // Parse the version from right after the matched pattern (e.g. right
    // after "vsearch v" or "salmon ") rather than blindly taking line 0 of
    // the combined output, which may not be the line the version is on.
    let search_from = if dep.version_pattern.is_empty() {
        combined.as_str()
    } else {
        match combined.find(dep.version_pattern) {
            Some(pos) => &combined[pos + dep.version_pattern.len()..],
            None => combined.as_str(),
        }
    };
    let display = search_from.lines().next().unwrap_or("").trim().to_string();
    debug!("  {} version output: {}", dep.name, display);

    let required = SemVer::parse(dep.version).with_context(|| {
        format!(
            "internal error: couldn't parse configured minimum version '{}' for '{}'",
            dep.version, dep.name
        )
    })?;

    match SemVer::parse(search_from) {
        Some(found) if found >= required => Ok(display),
        Some(found) => bail!(
            "found {} {}, but transfuse requires >= {}",
            dep.name, found, dep.version
        ),
        None => bail!(
            "found '{}' on PATH but couldn't parse a version number from its output: {}",
            dep.name, display
        ),
    }
}

fn download_and_install(name: &str, url: &str, dir: &Path) -> Result<()> {
    let tmp = tempfile::tempdir().context("Failed to create temp dir")?;
    let archive_name = url.split('/').last().unwrap_or("archive");
    let archive_path = tmp.path().join(archive_name);

    let status = Command::new("curl").args(["-L", "-o"]).arg(&archive_path).arg(url)
        .status().context("Failed to run curl. Please install curl.")?;
    if !status.success() { bail!("curl failed for {}", url); }

    let ext = archive_name.to_lowercase();
    if ext.ends_with(".tar.gz") || ext.ends_with(".tgz") {
        Command::new("tar").args(["xzf"]).arg(&archive_path).arg("-C").arg(tmp.path())
            .status().context("Failed to run tar")?;
    } else if ext.ends_with(".tar.bz2") {
        Command::new("tar").args(["xjf"]).arg(&archive_path).arg("-C").arg(tmp.path())
            .status().context("Failed to run tar")?;
    } else if ext.ends_with(".zip") {
        Command::new("unzip").arg("-q").arg(&archive_path).arg("-d").arg(tmp.path())
            .status().context("Failed to run unzip")?;
    }

    let binary_path = find_binary_in_dir(tmp.path(), name)?;
    let dest = dir.join(name);
    std::fs::copy(&binary_path, &dest)
        .with_context(|| format!("Failed to copy {:?} to {:?}", binary_path, dest))?;

    #[cfg(unix)]
    {
        use std::os::unix::fs::PermissionsExt;
        let mut perms = std::fs::metadata(&dest)?.permissions();
        perms.set_mode(0o755);
        std::fs::set_permissions(&dest, perms)?;
    }
    Ok(())
}

fn find_binary_in_dir(dir: &Path, name: &str) -> Result<PathBuf> {
    for entry in walkdir(dir) {
        if entry.file_name().to_string_lossy() == name {
            return Ok(entry.path());
        }
    }
    bail!("Could not find '{}' in unpacked archive", name)
}

fn walkdir(dir: &Path) -> Vec<std::fs::DirEntry> {
    let mut results = Vec::new();
    if let Ok(entries) = std::fs::read_dir(dir) {
        for entry in entries.flatten() {
            let path = entry.path();
            if path.is_dir() { results.extend(walkdir(&path)); }
            else { results.push(entry); }
        }
    }
    results
}

fn default_install_dir() -> PathBuf {
    std::env::var("HOME").map(PathBuf::from)
        .unwrap_or_else(|_| PathBuf::from("/tmp"))
        .join(".local").join("bin")
}

// ── Tests ─────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_semver_parse_basic() {
        assert_eq!(SemVer::parse("2.3.4"), Some(SemVer(2, 3, 4)));
        assert_eq!(SemVer::parse("salmon 2.3.4"), Some(SemVer(2, 3, 4)));
        assert_eq!(SemVer::parse("1.19.2"), Some(SemVer(1, 19, 2)));
    }

    #[test]
    fn test_semver_parse_missing_components_default_to_zero() {
        assert_eq!(SemVer::parse("2"), Some(SemVer(2, 0, 0)));
        assert_eq!(SemVer::parse("1.0"), Some(SemVer(1, 0, 0)));
    }

    #[test]
    fn test_semver_parse_stops_at_non_numeric_suffix() {
        // e.g. minimap2's bare "2.26-r1175" and vsearch's "2.28.1_linux_x86_64"
        assert_eq!(SemVer::parse("2.26-r1175"), Some(SemVer(2, 26, 0)));
        assert_eq!(SemVer::parse("2.28.1_linux_x86_64, 15.6GB RAM"), Some(SemVer(2, 28, 1)));
    }

    #[test]
    fn test_semver_parse_finds_version_after_leading_text() {
        // The exact bug scenario: pattern-matched text before the real version.
        assert_eq!(SemVer::parse("Rognes T, Flouri T (2016) vsearch v2.28.1"), Some(SemVer(2016, 0, 0)));
        // Note: this deliberately shows why searching from right after the
        // matched pattern (not the whole raw string) matters - see probe_binary.
    }

    #[test]
    fn test_semver_ordering() {
        assert!(SemVer(2, 3, 4) >= SemVer(1, 0, 0));
        assert!(SemVer(0, 4, 0) < SemVer(1, 0, 0)); // the exact bug this fixes
        assert!(SemVer(1, 19, 2) >= SemVer(1, 0, 0));
        assert!(SemVer(1, 0, 0) >= SemVer(1, 0, 0));
    }

    #[test]
    fn test_semver_parse_no_digits_returns_none() {
        assert_eq!(SemVer::parse("no version here"), None);
    }
}