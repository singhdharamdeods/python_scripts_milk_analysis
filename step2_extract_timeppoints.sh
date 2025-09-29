#!/bin/bash
set -euo pipefail

# === USER INPUTS ===
SUMMARY=$(ls sequencing_summary_*.txt | head -n 1)  # auto-detect summary file
if [[ ! -f "$SUMMARY" ]]; then
    echo "ERROR: No sequencing_summary_*.txt file found in current directory."
    exit 1
fi

# Generate cutoffs: every 2h (7200s) until 24h (86400s)
CUTOFFS=$(seq 7200 7200 86400)

# Master output dir
TOPDIR="timepoints"
mkdir -p "$TOPDIR"

# CSV file for extracted counts
CSV_FILE="${TOPDIR}/extracted_counts.csv"
echo "Cutoff_s,File,Full_Reads,Subset_Reads,Percent" > "$CSV_FILE"

# === LOOP OVER CUTOFFS ===
for CUTOFF in $CUTOFFS; do
    OUTDIR="${TOPDIR}/${CUTOFF}s_fastq_starttime"
    mkdir -p "$OUTDIR"

    echo "========================================================"
    echo ">> Processing cutoff = ${CUTOFF}s ($(echo "$CUTOFF/3600" | bc)h)"
    echo "========================================================"

    # Extract read IDs that START within cutoff
    awk -F'\t' -v cutoff="$CUTOFF" 'NR>1 && $10 <= cutoff {print $5}' "$SUMMARY" > readIDs_starttime.txt
    echo ">> Found $(wc -l < readIDs_starttime.txt) reads that START within ${CUTOFF}s."

    # Subset FASTQs
    shopt -s nullglob
    FASTQS=( *.fastq.gz )

    if [[ ${#FASTQS[@]} -eq 0 ]]; then
        echo "ERROR: No FASTQ files found in current directory."
        exit 1
    fi

    for f in "${FASTQS[@]}"; do
        if [[ ! -s "$f" ]]; then
            echo "SKIP: $f is missing or empty."
            continue
        fi

        base=$(basename "$f" .fastq.gz)
        out="${OUTDIR}/${base}_${CUTOFF}s_start.fastq.gz"

        echo ">> Subsetting $f -> $out"
        seqtk subseq "$f" readIDs_starttime.txt | gzip > "$out"

        if [[ -s "$out" ]]; then
            subset=$(zgrep -c "^@" "$out")
            full=$(zgrep -c "^@" "$f")
            pct=$(echo "scale=2; ($subset/$full)*100" | bc)

            echo "DONE: $out ($subset reads, ${pct}%)"
            echo "${CUTOFF},${base},${full},${subset},${pct}" >> "$CSV_FILE"
        else
            echo "SKIP: No reads ≤ ${CUTOFF}s found in $f"
            rm -f "$out"
            # still record zero
            full=$(zgrep -c "^@" "$f")
            echo "${CUTOFF},${base},${full},0,0" >> "$CSV_FILE"
        fi
    done

    # === LAST END-TIME REPORT PER BARCODE ===
    echo
    echo "================= LAST END TIME PER BARCODE (${CUTOFF}s) ================="
    awk -F'\t' -v cutoff="$CUTOFF" '
    NR==1 { next }
    $10 <= cutoff {
        end_time = $10 + $11
        bc = $26
        if (bc != "unclassified") {
            if (end_time > maxend[bc]) {
                maxend[bc] = end_time
            }
            count[bc]++
        }
    }
    END {
        printf("%-15s %-15s %-15s %-12s\n", "Barcode", "Reads_≤cutoff", "Last_end_time(s)", "Last_end_time(hh:mm:ss)")
        for (bc in maxend) {
            sec = int(maxend[bc])
            hh = int(sec/3600)
            mm = int((sec%3600)/60)
            ss = int(sec%60)
            printf("%-15s %-15d %-15.2f %02d:%02d:%02d\n", bc, count[bc], maxend[bc], hh, mm, ss)
        }
    }' "$SUMMARY" | sort
    echo "=========================================================================="
done

echo
echo ">> All cutoff summaries saved to: $CSV_FILE"

