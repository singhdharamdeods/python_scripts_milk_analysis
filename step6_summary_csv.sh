#!/usr/bin/env python3
import glob, os, re
import pandas as pd

# === CONFIG ===
BASE = "/home/dharamdeo/other_projects/milk/sept_11_2025/sample/bracken"
OUTDIR = "/home/dharamdeo/other_projects/milk/sept_11_2025/sample/summarycsv"
os.makedirs(OUTDIR, exist_ok=True)

def natural_key(s):
    """Sort like humans (barcode1, barcode2, ... not barcode10 before barcode2)."""
    return [int(t) if t.isdigit() else t for t in re.split(r'(\d+)', s)]

def clean_sample_name(filename):
    """
    Extract sample name, e.g.
    'barcode01_14400s_start.filtered.bracken.S.txt' -> 'barcode01_14400s_start'
    """
    name = os.path.basename(filename)
    name = name.replace(".bracken.S.txt", "")
    name = name.replace(".filtered", "")   # strip .filtered
    return name

def load_bracken_species_counts(path):
    """
    Load a Bracken .bracken.S.txt file and return {species_name: new_est_reads}
    """
    df = pd.read_csv(path, sep="\t", header=0, dtype={"name": str})
    df["name"] = df["name"].astype(str).str.strip()

    # Keep only rows with nonzero estimated reads
    df = df[df["new_est_reads"] > 0]

    return dict(zip(df["name"], df["new_est_reads"]))

def main():
    # Find all Bracken outputs in all timepoint subdirs
    pattern = os.path.join(BASE, "*_starttime", "*.bracken.S.txt")
    files = sorted(glob.glob(pattern))
    if not files:
        print("No .bracken.S.txt files found under", BASE)
        return

    # Build union of species across all samples
    all_species = set()
    per_sample = {}
    for f in files:
        sample = clean_sample_name(f)
        counts = load_bracken_species_counts(f)
        per_sample[sample] = counts
        all_species.update(counts.keys())

    species_sorted = sorted(all_species)
    samples_sorted = sorted(per_sample.keys(), key=natural_key)

    # Counts matrix
    counts_df = pd.DataFrame(index=species_sorted, columns=samples_sorted).fillna(0).astype(int)
    for s in samples_sorted:
        for sp, c in per_sample[s].items():
            counts_df.at[sp, s] = c

    # Percent abundance matrix
    perc_df = counts_df.div(counts_df.sum(axis=0), axis=1) * 100
    perc_df = perc_df.fillna(0)

    # Combine counts and percents into one dataframe
    combined_df = pd.DataFrame(index=species_sorted)
    for s in samples_sorted:
        combined_df[f"{s}_count"] = counts_df[s]
        combined_df[f"{s}_percent"] = perc_df[s]

    combined_df.index.name = "species"
    combined_path = os.path.join(OUTDIR, "merged.bracken_S.counts+percent.csv")
    combined_df.to_csv(combined_path)

    # Presence/absence matrix
    presence_df = (counts_df > 0).astype(int)
    presence_df.index.name = "species"
    presence_path = os.path.join(OUTDIR, "merged.bracken_S.presence.csv")
    presence_df.to_csv(presence_path)

    print(f"Wrote {combined_path}")
    print(f"Wrote {presence_path}")

if __name__ == "__main__":
    main()

