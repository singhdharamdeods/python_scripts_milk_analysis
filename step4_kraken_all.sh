#!/bin/bash
set -euo pipefail
shopt -s nullglob

# === USER INPUTS ===
DB="/media/hd3/databases/kraken2_db_NEW"   # Kraken2 database path
THREADS=4
OPTS="--use-names --gzip-compressed --confidence 0.01"

# Input: filtered FASTQs
FILTERED_BASE="/home/dharamdeo/other_projects/milk/sept_11_2025/sample/filtered"

# Output: central kraken_out directory
OUTBASE="/home/dharamdeo/other_projects/milk/sept_11_2025/sample/kraken_out"
mkdir -p "$OUTBASE"

# Loop over each *_starttime directory in filtered
for dir in "$FILTERED_BASE"/*_starttime; do
  if [[ ! -d "$dir" ]]; then
    continue
  fi

  subdir=$(basename "$dir")
  OUTDIR="$OUTBASE/$subdir"
  mkdir -p "$OUTDIR"

  echo "Processing filtered directory: $dir"
  echo "Kraken2 outputs will be in: $OUTDIR"

  # Loop through all filtered FASTQs
  for IN in "$dir"/*.fastq.gz; do
    base=$(basename "$IN" .fastq.gz)
    OUT="$OUTDIR/${base}"

    if [[ ! -s "$IN" ]]; then
      echo "WARN: $IN missing/empty, skipping"
      continue
    fi

    echo ">> Kraken2 $base"
    kraken2 --db "$DB" --threads "$THREADS" $OPTS \
      --report "${OUT}.kreport.txt" \
      --output "${OUT}.kraken.txt" \
      --classified-out "${OUT}.classified.fastq" \
      --unclassified-out "${OUT}.unclassified.fastq" \
      "$IN" || { echo "ERROR: kraken2 failed on $base"; continue; }

    # Compress classified/unclassified outputs
    gzip -f "${OUT}.classified.fastq" "${OUT}.unclassified.fastq"
  done
done

echo "All done. Results are in: $OUTBASE"

