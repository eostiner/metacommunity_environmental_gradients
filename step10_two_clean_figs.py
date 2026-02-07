#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
step10_two_clean_figs.py
"""

import argparse
import os
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


# ---------- small helpers ----------

def norm_ids(s):
    s = s.astype(str).str.strip()
    s = s.str.replace(r"\s+", "_", regex=True)
    s = s.str.replace(r"[^A-Za-z0-9_\-]", "", regex=True)
    s = s.str.replace(r"_+", "_", regex=True)
    return s


def detect_id_col(df, candidates):
    cols = list(df.columns)
    for c in candidates:
        if c in cols:
            return c
    lower = {c.lower(): c for c in cols}
    for cand in candidates:
        lc = cand.lower()
        if lc in lower:
            return lower[lc]
    raise ValueError(f"Could not find an ID column in metadata. Tried: {candidates}")


def detect_p_col(df):
    candidates = ["p_value", "p.value", "p.value.", "pval", "p"]
    for c in candidates:
        if c in df.columns:
            return c
    return None


def main():
    ap = argparse.ArgumentParser(
        description="Step 10 figures: LCBD per ranch and indicator taxa heatmap."
    )
    ap.add_argument("--step10_dir", required=True)
    ap.add_argument("--meta_file", required=True)
    ap.add_argument("--group_col", default="ranch")
    ap.add_argument("--alpha", type=float, default=0.05)
    ap.add_argument("--top_k", type=int, default=25)
    ap.add_argument("--out_base", default="Step10_Indicators")
    args = ap.parse_args()

    step10_dir = args.step10_dir

    # ---------- paths ----------
    ind_path = os.path.join(step10_dir, "indicators_full.csv")
    lcbd_path = os.path.join(step10_dir, "LCBD.csv")

    if not os.path.exists(ind_path):
        print(f"[ERR] Missing indicators_full.csv in {step10_dir}", file=sys.stderr)
        sys.exit(1)
    if not os.path.exists(lcbd_path):
        print(f"[ERR] Missing LCBD.csv in {step10_dir}", file=sys.stderr)
        sys.exit(1)
    if not os.path.exists(args.meta_file):
        print(f"[ERR] Metadata not found: {args.meta_file}", file=sys.stderr)
        sys.exit(1)

    # ---------- load metadata ----------
    meta = pd.read_csv(args.meta_file)
    id_col = detect_id_col(meta, ["sample_id", "SampleID", "sample", "id", "ID"])
    if args.group_col not in meta.columns:
        raise ValueError(f"[metadata] group_col '{args.group_col}' not found.")

    meta[id_col] = norm_ids(meta[id_col])
    meta[args.group_col] = meta[args.group_col].astype(str)

    # ---- your desired ranch plotting order ----
    desired_ranch_order = ["La_Venada", "El_Mesquite", "Casa_Vieja", "El_Tule"]

    # ---------- Panel A: LCBD ----------
    lcbd = pd.read_csv(lcbd_path)
    if "sample_id" not in lcbd.columns:
        raise ValueError("LCBD.csv must have 'sample_id'.")

    lcbd["sample_id"] = norm_ids(lcbd["sample_id"])

    merged = lcbd.merge(
        meta[[id_col, args.group_col]],
        left_on="sample_id", right_on=id_col,
        how="left"
    )
    merged = merged.dropna(subset=[args.group_col, "LCBD"])

    # Compute mean LCBD + enforce ranch order
    lcbd_stats = merged.groupby(args.group_col)["LCBD"].agg(["mean", "sem"])

    ranches = [r for r in desired_ranch_order if r in lcbd_stats.index]
    if not ranches:
        ranches = list(lcbd_stats.index)

    lmean = lcbd_stats.loc[ranches, "mean"]
    lse = lcbd_stats.loc[ranches, "sem"]

    figA, axA = plt.subplots(figsize=(6, 4), dpi=300)
    x = np.arange(len(ranches))
    axA.bar(x, lmean.values, yerr=lse.values, capsize=4)
    axA.set_xticks(x)
    axA.set_xticklabels(ranches, rotation=20, ha="right")
    axA.set_ylabel("Mean LCBD")
    axA.set_title("Local contributions to beta diversity by ranch")
    figA.tight_layout()

    figA.savefig(f"{args.out_base}_LCBD.png", bbox_inches="tight")
    figA.savefig(f"{args.out_base}_LCBD.pdf", bbox_inches="tight")
    plt.close(figA)

    # ---------- Panel B: indicator taxa ----------
    ind = pd.read_csv(ind_path)
    if "taxon" not in ind.columns:
        raise ValueError("indicators_full.csv must have 'taxon'.")

    p_col = detect_p_col(ind)
    if p_col is None:
        ind["p_tmp"] = 1.0
        p_col = "p_tmp"

    ind_sorted = ind.sort_values(by=p_col)
    sig = ind_sorted[ind_sorted[p_col] <= args.alpha]
    chosen = sig.head(args.top_k) if len(sig) > 0 else ind_sorted.head(args.top_k)

    ranch_cols_raw = [c for c in ind.columns if c.startswith("s.")]
    if not ranch_cols_raw:
        raise ValueError("No 's.<ranch>' columns in indicators_full.csv.")

    col_to_ranch = {c: c[2:] for c in ranch_cols_raw}
    ranch_to_col = {r: c for c, r in col_to_ranch.items()}

    # enforce ranch plotting order
    ranch_names = [r for r in desired_ranch_order if r in ranch_to_col]
    if not ranch_names:
        ranch_names = list(ranch_to_col.keys())

    ranch_cols = [ranch_to_col[r] for r in ranch_names]

    if "stat" not in ind.columns:
        raise ValueError("Missing 'stat' column for IndVal.")

    rows = []
    for _, row in chosen.iterrows():
        entry = {"taxon": row["taxon"], "p": float(row[p_col])}
        for c, rname in zip(ranch_cols, ranch_names):
            val = float(row["stat"]) if row.get(c, 0) > 0 else 0.0
            entry[rname] = val
        rows.append(entry)

    mat = pd.DataFrame(rows).sort_values(by="p")
    taxa = list(mat["taxon"])
    vals = mat[ranch_names].values

    n_taxa = len(taxa)
    fig_h = max(4, 0.35 * n_taxa)
    figB, axB = plt.subplots(figsize=(8, fig_h), dpi=300)

    im = axB.imshow(vals, aspect="auto", interpolation="nearest")

    axB.set_xticks(np.arange(len(ranch_names)))
    axB.set_xticklabels(ranch_names, rotation=25, ha="right")
    axB.set_yticks(np.arange(n_taxa))
    axB.set_yticklabels(taxa)
    axB.set_xlabel("Ranch")
    axB.set_ylabel("Taxon")
    axB.set_title(f"Indicator taxa by ranch (IndVal stat, alpha = {args.alpha})")

    cbar = figB.colorbar(im, ax=axB)
    cbar.set_label("IndVal statistic")

    for i, (_, row) in enumerate(mat.iterrows()):
        for j, rname in enumerate(ranch_names):
            sc = "s." + rname
            src = chosen[chosen["taxon"] == row["taxon"]].iloc[0]
            if sc in chosen.columns and src.get(sc, 0) > 0 and row["p"] <= args.alpha:
                axB.scatter(j, i, s=40, facecolors="none", edgecolors="black", linewidths=0.8)

    figB.tight_layout()
    figB.savefig(f"{args.out_base}_IndValHeatmap.png", bbox_inches="tight")
    figB.savefig(f"{args.out_base}_IndValHeatmap.pdf", bbox_inches="tight")
    plt.close(figB)


if __name__ == "__main__":
    main()