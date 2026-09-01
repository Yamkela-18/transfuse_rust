#!/bin/bash
#SBATCH --job-name=merge609894
#SBATCH --output=merge609894_%j.out
#SBATCH --error=merge609894_%j.err
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=72:00:00

set -e

cd ~/transfuse_rust

# --------------------------------------------------
# Environment
# --------------------------------------------------

module load anaconda3/2024.10
conda activate transfuse_env

THRESHOLD=0.05

OUTPUT="data/output/SRR609894_new_1_merged_0.05.fasta"

echo "========================================"
echo "Transfuse Merger - SRR609894"
echo "Job ID: $SLURM_JOB_ID"
echo "Node: $(hostname)"
echo "Start: $(date)"
echo "Threshold: $THRESHOLD"
echo "========================================"

mkdir -p data/output

# --------------------------------------------------
# Check dependencies
# --------------------------------------------------

echo "========================================"
echo "Checking dependencies"
echo "========================================"

which vsearch
which minimap2
which samtools
which salmon
which cargo

# --------------------------------------------------
# Automatically find COMPLETE score file
# --------------------------------------------------

echo "========================================"
echo "Searching for score files"
echo "========================================"

SCORE_FILE=""
MAX_RECORDS=0

for file in data/output/*.csv; do

    [ -f "$file" ] || continue

    HEADER=$(head -n 1 "$file")

    # Check that the file is TAB-separated
    # and has the expected Transfuse columns.
    if echo "$HEADER" | grep -q $'\t' &&
       [ "$(echo "$HEADER" | awk -F'\t' '{print $1}')" = "contig_name" ] &&
       [ "$(echo "$HEADER" | awk -F'\t' '{print $2}')" = "score" ] &&
       [ "$(echo "$HEADER" | awk -F'\t' '{print $3}')" = "p_good" ] &&
       [ "$(echo "$HEADER" | awk -F'\t' '{print $4}')" = "p_bases_covered" ] &&
       [ "$(echo "$HEADER" | awk -F'\t' '{print $5}')" = "coverage" ]; then

        RECORDS=$(awk 'END {print NR-1}' "$file")

        echo "Valid score file:"
        echo "  $file"
        echo "  Records: $RECORDS"

        # Select the file with the most records.
        if [ "$RECORDS" -gt "$MAX_RECORDS" ]; then
            MAX_RECORDS=$RECORDS
            SCORE_FILE="$file"
        fi
    fi

done

# --------------------------------------------------
# Make sure a score file was found
# --------------------------------------------------

if [ -z "$SCORE_FILE" ]; then
    echo "ERROR: No valid Transfuse score file found."
    exit 1
fi

echo "========================================"
echo "SELECTED SCORE FILE"
echo "========================================"
echo "$SCORE_FILE"
echo "Records: $MAX_RECORDS"

# --------------------------------------------------
# Verify TAB-separated format
# --------------------------------------------------

echo "========================================"
echo "Checking score format"
echo "========================================"

FIELDS=$(head -n 1 "$SCORE_FILE" | awk -F'\t' '{print NF}')

if [ "$FIELDS" -ne 5 ]; then
    echo "ERROR: Expected 5 TAB-separated columns."
    echo "Found: $FIELDS"
    exit 1
fi

echo "Format: TAB-separated"
echo "Columns: $FIELDS"

head -n 1 "$SCORE_FILE"

# --------------------------------------------------
# Verify column 2 is score
# --------------------------------------------------

SECOND_COLUMN=$(head -n 1 "$SCORE_FILE" | awk -F'\t' '{print $2}')

if [ "$SECOND_COLUMN" != "score" ]; then
    echo "ERROR: Column 2 is not score."
    echo "Column 2 = $SECOND_COLUMN"
    exit 1
fi

echo "Column 2 confirmed as: score"

# --------------------------------------------------
# Score statistics
# --------------------------------------------------

echo "========================================"
echo "Score statistics"
echo "========================================"

TOTAL=$(awk -F'\t' '
NR > 1 {count++}
END {print count+0}
' "$SCORE_FILE")

PASS=$(awk -F'\t' -v threshold="$THRESHOLD" '
NR > 1 && ($2+0) >= threshold {count++}
END {print count+0}
' "$SCORE_FILE")

FAIL=$(awk -F'\t' -v threshold="$THRESHOLD" '
NR > 1 && ($2+0) < threshold {count++}
END {print count+0}
' "$SCORE_FILE")

PERCENT=$(awk -v pass="$PASS" -v total="$TOTAL" '
BEGIN {
    if (total > 0)
        printf "%.2f", (pass/total)*100
    else
        print "0.00"
}')

echo "Total contigs:       $TOTAL"
echo "Score >= $THRESHOLD: $PASS"
echo "Score < $THRESHOLD:  $FAIL"
echo "Passing percentage:  $PERCENT%"

# --------------------------------------------------
# Check assemblies
# --------------------------------------------------

echo "========================================"
echo "Checking assemblies"
echo "========================================"

ASSEMBLIES="data/assemblies/SRR609894_soapdenovo.fasta,data/assemblies/SRR609894_transcripts.fasta,data/assemblies/SRR609894_tri_assembly.okay.fasta,data/assemblies/SRR609894_trinity.Trinity.fasta"

ls -lh \
    data/assemblies/SRR609894_soapdenovo.fasta \
    data/assemblies/SRR609894_transcripts.fasta \
    data/assemblies/SRR609894_tri_assembly.okay.fasta \
    data/assemblies/SRR609894_trinity.Trinity.fasta

# --------------------------------------------------
# Check reads
# --------------------------------------------------

echo "========================================"
echo "Checking paired reads"
echo "========================================"

LEFT="data/trimmed/SRR609894/SRR609894_1_paired.fastq"
RIGHT="data/trimmed/SRR609894/SRR609894_2_paired.fastq"

ls -lh "$LEFT"
ls -lh "$RIGHT"

# --------------------------------------------------
# Run merger
# --------------------------------------------------

echo "========================================"
echo "STARTING TRANSFUSE MERGER"
echo "========================================"

echo "Complete score file:"
echo "$SCORE_FILE"

echo "Threshold:"
echo "$THRESHOLD"

echo "Rust will use the COMPLETE score file."
echo "Rust will apply --min-score $THRESHOLD."
echo "========================================"

START=$(date +%s)

cargo run --release -- \
    --assemblies "$ASSEMBLIES" \
    --left "$LEFT" \
    --right "$RIGHT" \
    --scores "$SCORE_FILE" \
    --min-score "$THRESHOLD" \
    --output "$OUTPUT"

END=$(date +%s)
TIME=$((END - START))

# --------------------------------------------------
# Completion
# --------------------------------------------------

echo "========================================"
echo "MERGER COMPLETED SUCCESSFULLY"
echo "========================================"

echo "End: $(date)"

echo "Runtime:"
echo "$((TIME/3600))h $(((TIME%3600)/60))m $((TIME%60))s"

echo "Final output:"
ls -lh "$OUTPUT"

echo "Merged contigs:"
grep -c '^>' "$OUTPUT"

echo "========================================"
echo "JOB FINISHED"
echo "========================================"
