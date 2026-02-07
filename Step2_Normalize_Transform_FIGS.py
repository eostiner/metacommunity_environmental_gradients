#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Step2_Normalize_Transform_FIGS.py
- Reads a taxa-by-samples count matrix (Genus or Taxon as the ID).
- Ignores extra taxonomy columns: Genus, species, Kingdom, Phylum, Class, Order, Family.
- (Optional) basic QC filtering by minimum sample totals and taxon prevalence.
- Writes:
    * raw_counts_passthrough.csv          (filtered raw counts)
    * hellinger_matrix.csv                (Hellinger-transformed)
    * presence_absence.csv                (0/1)
    * qc_summary.txt                      (quick stats + thresholds)
- Plots (PNG):
    * raw_counts_heatmap.png              (log(count+1))
    * hellinger_heatmap.png
    * presence_absence_heatmap.png
    * sample_totals_bar.png               (library sizes)
    * sample_richness_bar.png             (# taxa present)

Usage example:
  python -u Step2_Normalize_Transform_FIGS.py \
    --counts "/path/USE_GENUS_ONLY_matrix_by_trap_MEAN_REL_WeighteAvg_copy.csv" \
    --outdir "/path/Normalized_Matrices" \
    --min-sample-total 0 \
    --min-taxon-prevalence 0
"""

import os
import argparse
import textwrap

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


# --------------------------- Helpers --------------------------- #

def read_counts(path: str) -> pd.DataFrame:
    """
    Reads an abundance matrix (taxa × samples), ignoring extra taxonomy columns.
    Accepts either 'Genus' or 'Taxon' as the identifier column, renames to 'Taxon',
    sets it as the index, and returns numeric sample columns only.
    """
    # Auto-detect delimiter
    with open(path, 'r', encoding='utf-8') as f:
        first_line = f.readline()
    sep = '\t' if first_line.count('\t') > first_line.count(',') else ','

    df = pd.read_csv(path, sep=sep)

    # Identify the main taxon column (Genus or Taxon)
    lower_cols = {c.lower(): c for c in df.columns}
    id_col = lower_cols.get('genus', lower_cols.get('taxon', None))
    if id_col is None:
        raise ValueError(
            "Input matrix must have a column named 'Genus' or 'Taxon'. "
            f"Found columns: {list(df.columns)}"
        )

    # Normalize the ID column name to 'Taxon'
    df = df.rename(columns={id_col: 'Taxon'})

    # Drop taxonomy metadata columns; keep only sample columns
    ignore_cols = {'genus', 'species', 'kingdom', 'phylum', 'class', 'order', 'family', 'taxon'}
    keep_cols = [c for c in df.columns if c.lower() not in ignore_cols]

    # Ensure 'Taxon' + sample columns in order
    df = df[['Taxon'] + [c for c in keep_cols if c != 'Taxon']]

    # Cast sample columns to numeric
    for c in df.columns[1:]:
        df[c] = pd.to_numeric(df[c], errors='coerce').fillna(0)

    # Set taxa as the INDEX so downstream math ignores the ID
    df['Taxon'] = df['Taxon'].astype(str)
    df = df.set_index('Taxon', drop=True)

    # Optionally drop taxa with all-zero rows (safety)
    if df.shape[1] > 0:
        df = df.loc[(df.sum(axis=1) > 0), :]

    return df


def hellinger_transform(counts: pd.DataFrame) -> pd.DataFrame:
    """Hellinger transform: sqrt( counts / sample_total )."""
    colsum = counts.sum(axis=0).astype(float)
    colsum[colsum == 0] = 1.0  # avoid divide-by-zero
    rel = counts.div(colsum, axis=1)
    return np.sqrt(rel)


def presence_absence_matrix(counts: pd.DataFrame) -> pd.DataFrame:
    return (counts > 0).astype(int)


def write_df(df: pd.DataFrame, path: str, index_label: str = "Taxon"):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    df.to_csv(path, index=True, index_label=index_label)


def plot_heatmap(matrix: pd.DataFrame, title: str, cbar_label: str, out_png: str):
    """Simple matplotlib heatmap for large/sparse matrices."""
    plt.figure(figsize=(10, 6))
    img = plt.imshow(matrix.values, aspect="auto", interpolation="nearest")
    plt.title(title)
    plt.xlabel("Samples")
    plt.ylabel("Taxa")
    cbar = plt.colorbar(img)
    cbar.set_label(cbar_label)
    ax = plt.gca()
    n = matrix.shape[1]
    step = max(1, n // 20)
    ax.set_xticks(np.arange(0, n, step))
    ax.set_xticklabels(matrix.columns[::step], rotation=90, fontsize=8)
    ax.set_yticks([])  # too many taxa to label legibly
    plt.tight_layout()
    plt.savefig(out_png, dpi=300)
    plt.close()


def plot_sample_totals(counts: pd.DataFrame, out_png: str):
    totals = counts.sum(axis=0)
    plt.figure(figsize=(10, 4))
    plt.bar(totals.index, totals.values)
    plt.ylabel("Total reads")
    plt.title("Sample library sizes")
    plt.xticks(rotation=90, fontsize=8)
    plt.tight_layout()
    plt.savefig(out_png, dpi=300)
    plt.close()


def plot_sample_richness(counts: pd.DataFrame, out_png: str):
    richness = (counts > 0).sum(axis=0)
    plt.figure(figsize=(10, 4))
    plt.bar(richness.index, richness.values)
    plt.ylabel("# taxa present")
    plt.title("Sample richness (presence/absence)")
    plt.xticks(rotation=90, fontsize=8)
    plt.tight_layout()
    plt.savefig(out_png, dpi=300)
    plt.close()


# --------------------------- CLI --------------------------- #

def get_parser():
    p = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description="Normalize and transform community matrices (Hellinger, presence/absence) and produce QC plots.",
        epilog=textwrap.dedent("""
        Notes:
          * Input must contain a 'Genus' or 'Taxon' column; the script will use it as ID.
          * Extra taxonomy columns (Genus, species, Kingdom, Phylum, Class, Order, Family) are ignored.
          * Use either --counts (original flag) or --community (alias).
        """)
    )
    p.add_argument("--counts", "--community", dest="counts", required=True,
                   help="Path to taxa-by-samples CSV/TSV (with 'Genus' or 'Taxon' ID column).")
    p.add_argument("--outdir", default="./RESULTS/step2",
                   help="Output directory (default: ./RESULTS/step2)")
    p.add_argument("--min-sample-total", type=int, default=0,
                   help="Drop samples with total counts < this value (default: 0 = keep all).")
    p.add_argument("--min-taxon-prevalence", type=int, default=0,
                   help="Drop taxa present in fewer than this many samples (default: 0 = keep all).")
    return p


# --------------------------- Main --------------------------- #

def main():
    args = get_parser().parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    # Read counts (rows = taxa, cols = samples; index = 'Taxon')
    counts = read_counts(args.counts)

    # Basic filtering (optional)
    # 1) sample totals filter
    sample_totals = counts.sum(axis=0)
    keep_samples = sample_totals.index
    if args.min_sample_total > 0:
        keep_samples = sample_totals[sample_totals >= args.min_sample_total].index
        dropped = set(counts.columns) - set(keep_samples)
        if dropped:
            print(f"[INFO] Dropping {len(dropped)} samples with totals < {args.min_sample_total}: "
                  f"{sorted(list(dropped))[:5]}...")
    counts = counts.loc[:, keep_samples]

    # 2) taxon prevalence filter
    if args.min_taxon_prevalence > 0:
        prevalence = (counts > 0).sum(axis=1)
        keep_taxa = prevalence[prevalence >= args.min_taxon_prevalence].index
        dropped_t = set(counts.index) - set(keep_taxa)
        if dropped_t:
            print(f"[INFO] Dropping {len(dropped_t)} taxa with prevalence < {args.min_taxon_prevalence}.")
        counts = counts.loc[keep_taxa, :]

    # Write raw passthrough
    raw_out = os.path.join(args.outdir, "raw_counts_passthrough.csv")
    write_df(counts, raw_out)

    # Hellinger
    hell = hellinger_transform(counts)
    hell_out = os.path.join(args.outdir, "hellinger_matrix.csv")
    write_df(hell, hell_out)

    # Presence/Absence
    pa = presence_absence_matrix(counts)
    pa_out = os.path.join(args.outdir, "presence_absence.csv")
    write_df(pa, pa_out)

    # QC summary
    qc_txt = os.path.join(args.outdir, "qc_summary.txt")
    with open(qc_txt, "w", encoding="utf-8") as fh:
        fh.write("Step2 QC Summary\n")
        fh.write("================\n\n")
        fh.write(f"Input counts: {args.counts}\n")
        fh.write(f"Output dir  : {args.outdir}\n\n")
        fh.write(f"Samples kept: {counts.shape[1]}\n")
        fh.write(f"Taxa kept   : {counts.shape[0]}\n")
        fh.write(f"Min sample total     : {args.min_sample_total}\n")
        fh.write(f"Min taxon prevalence : {args.min_taxon_prevalence}\n\n")
        fh.write("Per-sample totals (first 10):\n")
        fh.write(sample_totals.loc[counts.columns].head(10).to_string() + "\n")

    # --------- Plots --------- #
    # Heatmaps
    plot_heatmap(np.log1p(counts), "Raw Counts (log1p)", "log(count+1)",
                 os.path.join(args.outdir, "raw_counts_heatmap.png"))
    plot_heatmap(hell, "Hellinger Transformed", "Hellinger",
                 os.path.join(args.outdir, "hellinger_heatmap.png"))
    plot_heatmap(pa, "Presence/Absence Matrix", "0/1",
                 os.path.join(args.outdir, "presence_absence_heatmap.png"))

    # Bar charts
    plot_sample_totals(counts, os.path.join(args.outdir, "sample_totals_bar.png"))
    plot_sample_richness(counts, os.path.join(args.outdir, "sample_richness_bar.png"))

    print(f"[OK] Wrote matrices to: {args.outdir}")
    print("[OK] Wrote plots: raw_counts_heatmap.png, hellinger_heatmap.png, presence_absence_heatmap.png, "
          "sample_totals_bar.png, sample_richness_bar.png")


if __name__ == "__main__":
    main()