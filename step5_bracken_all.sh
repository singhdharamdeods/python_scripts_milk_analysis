#!/bin/bash
set -euo pipefail
shopt -s nullglob

# === USER INPUTS ===
DB="/media/hd3/databases/kraken2_db_NEW"   # Bracken database (same as Kraken2)
BASE="/home/dharamdeo/other_projects/milk/sept_11_2025/sample/kraken_out"

# Bracken parameters
RLEN=100   # read length used to build Bracken DB
THRESH=10  # abundance threshold (%)

# Central output directory
OUTBASE="/home/dharamdeo/other_projects/milk/sept_11_2025/sample/bracken"
mkdir -p "$OUTBASE"

# Loop over all *_starttime directories in kraken_out
for KDIR in "$BASE"/*_starttime; do
  if [[ ! -d "$KDIR" ]]; then
    continue
  fi

  subdir=$(basename "$KDIR")
  OUTDIR="$OUTBASE/$subdir"
  mkdir -p "$OUTDIR"

  echo "Processing Bracken for: $KDIR -> $OUTDIR"

  # Loop over all Kraken report files in this directory
  for KR in "$KDIR"/*.kreport.txt; do
    [[ -s "$KR" ]] || { echo "Skipping empty: $KR"; continue; }
    base=$(basename "$KR" .kreport.txt)

    echo ">> Bracken species for $base"
    bracken -d "$DB" -i "$KR" \
            -o "$OUTDIR/${base}.bracken.S.txt" \
            -w "$OUTDIR/${base}.bracken.S.with_kreport.txt" \
            -l S -r "$RLEN" -t "$THRESH"
  done
done

echo "All Bracken analyses done. Results are in: $OUTBASE"

