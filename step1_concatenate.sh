#!/bin/bash
set -euo pipefail
shopt -s nullglob

# Move into your dataset directory
cd /home/dharamdeo/other_projects/milk/20250911_1300_fastq_pass

# Output directory for concatenated FASTQs
mkdir -p sample

# Loop through all barcode directories
for b in barcode{01..24}; do
  files=( "$b"/*.fastq.gz )
  if (( ${#files[@]} == 0 )); then
    echo "No FASTQs in $b, skipping."
    continue
  fi
  # Concatenate preserving gzip compression
  cat "${files[@]}" > "sample/${b}.fastq.gz"
  echo "Created sample/${b}.fastq.gz"
done

