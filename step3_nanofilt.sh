#!/bin/bash
set -euo pipefail
shopt -s nullglob

cd /home/dharamdeo/other_projects/milk/20250911_1300_fastq_pass/sample

mkdir -p filtered

# Max parallel jobs (4 jobs × 3 threads each = 12 cores)
MAXJOBS=4
THREADS=3

for b in {01..24}; do
  in="barcode${b}.fastq.gz"
  out="filtered/barcode${b}.filtered.fastq.gz"

  if [[ ! -f "$in" ]]; then
    echo "Input $in not found, skipping."
    continue
  fi

  (
    pigz -dc -p $THREADS "$in" \
      | NanoFilt -q 8 --headcrop 15 --length 100 \
      | pigz -c -p $THREADS > "$out"
    echo "Filtered -> $out"
  ) &

  # Limit concurrency
  if [[ $(jobs -r -p | wc -l) -ge $MAXJOBS ]]; then
    wait -n
  fi
done

wait

