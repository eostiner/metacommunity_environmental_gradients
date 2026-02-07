#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Step10_Genus_Top20_DotTable_pub_fixed.py  (wide-sites compatible, auto mkdir)

Supports two input modes:

A) SAMPLE MODE (original): comm_file rows = samples (first col sample_id)
   + metadata (samples.csv) used to group samples into sites (ranch).
B) WIDE-SITES MODE (new):  comm_file rows = taxa (Genus or Genus/species),
   columns = site replicates (e.g., Casa_Vieja_1, El_Mesquite_2, ...).
   No metadata required.

Key opts:
  --wide_sites yes|no|auto (default auto)
  --comm_file <csv>
  --meta_file <csv>                   # optional when wide-sites is used
  --group_col ranch
  --top_per_site 20
  --metric relative|absolute (default relative = %)
  --ranch_order <sites...>            # x-axis order (low→high elevation)
  --orient_elevation yes|no (default yes)
  --out <basename>
"""

import argparse, os, re, sys
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt

# ---------------- UTILITIES ---------------- #

def norm_ids(s: pd.Series) -> pd.Series:
    s = s.astype(str).str.strip()
    s = s.str.replace(r"\s+", "_", regex=True)
    s = s.str.replace(r"[^A-Za-z0-9_\-]", "", regex=True)
    s = s.str.replace(r"_+", "_", regex=True)
    return s

def detect_col(df, candidates, required=False, ctx=""):
    for c in candidates:
        if c in df.columns:
            return c
    lower = {c.lower(): c for c in df.columns}
    for cand in candidates:
        forms = {cand.lower(), cand.lower().replace(".", "_"), cand.lower().replace("_", ".")}
        for f in forms:
            if f in lower:
                return lower[f]
    if required:
        raise ValueError(f"[{ctx}] Need one of {candidates}; have {list(df.columns)}")
    return None

def scale_sizes(vals, min_size=25, max_size=900):
    v = np.asarray(vals, float)
    v = np.maximum(v, 0)
    if not np.isfinite(v).any() or np.nanmax(v) <= 0:
        return np.full_like(v, (min_size + max_size) / 2.0)
    v = v / np.nanmax(v)
    return min_size + v * (max_size - min_size)

def token_site(name: str) -> str:
    """From column 'Casa_Vieja_1' return 'Casa_Vieja'."""
    tok = str(name)
    if "_" in tok:
        parts = tok.split("_")
        if parts[-1].isdigit():
            return "_".join(parts[:-1])
    return tok

def looks_like_wide_sites(df: pd.DataFrame) -> bool:
    cols = [c for c in df.columns]
    first = cols[0].lower()
    siteish = sum(1 for c in cols[1:] if re.search(r"_\d+$", str(c)) is not None)
    return first in {"genus", "genus_species", "taxon", "species"} or siteish >= 2

# ---------------- MAIN ---------------- #

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--comm_file", required=True)
    ap.add_argument("--meta_file", required=False)
    ap.add_argument("--group_col", default="ranch")
    ap.add_argument("--top_per_site", type=int, default=20)
    ap.add_argument("--metric", choices=["relative", "absolute"], default="relative")
    ap.add_argument("--ranch_order", nargs="*", default=None)
    ap.add_argument("--orient_elevation", choices=["yes", "no"], default="yes")
    ap.add_argument("--wide_sites", choices=["yes", "no", "auto"], default="auto")
    ap.add_argument("--out", default="Step10_Genus_Top20_DotTable")
    ap.add_argument("--dpi", type=int, default=300)
    ap.add_argument("--font_family", default="DejaVu Sans")
    ap.add_argument("--font_scale", type=float, default=1.0)
    args = ap.parse_args()

    # ---- fonts ----
    base = 10.0 * args.font_scale
    mpl.rcParams.update({
        "font.family": args.font_family,
        "font.size": base,
        "axes.titlesize": base * 1.25,
        "axes.labelsize": base * 1.1,
        "xtick.labelsize": base * 0.95,
        "ytick.labelsize": base * 0.95,
        "legend.fontsize": base * 0.95,
        "figure.constrained_layout.use": True,
    })

    # ---- Read community file ----
    comm_raw = pd.read_csv(args.comm_file)
    mode = args.wide_sites
    if mode == "auto":
        mode = "yes" if looks_like_wide_sites(comm_raw) else "no"

    rows = []
    size_label_unit = "%"

    # ==================== WIDE-SITES MODE ==========================
    if mode == "yes":
        taxon_col = detect_col(comm_raw, ["Genus_species","Genus","Taxon","genus_species","genus","taxon"], required=False, ctx="comm")
        if taxon_col is None:
            taxon_col = comm_raw.columns[0]
        taxa_df = comm_raw.copy()
        taxa_df[taxon_col] = taxa_df[taxon_col].astype(str)

        # numeric columns are site data
        for c in taxa_df.columns:
            if c != taxon_col:
                taxa_df[c] = pd.to_numeric(taxa_df[c], errors="coerce")

        site_cols = [c for c in taxa_df.columns if c != taxon_col]
        col_site = {c: token_site(c) for c in site_cols}
        sites_all = sorted(pd.unique(list(col_site.values())), key=lambda x: x)

        # average replicates per site
        site_strengths = {}
        for s in sites_all:
            cols = [c for c in site_cols if col_site[c] == s]
            site_strengths[s] = taxa_df[cols].mean(axis=1)

        site_strengths_df = pd.DataFrame(site_strengths)
        if args.metric == "relative":
            col_sums = site_strengths_df.sum(axis=0).replace(0, np.nan)
            site_strengths_df = site_strengths_df.divide(col_sums, axis=1) * 100.0
            size_label_unit = "%"
        else:
            size_label_unit = ""

        default_baja = ["La_Venada", "El_Mesquite", "Casa_Vieja", "El_Tule"]
        if args.ranch_order and all(s in site_strengths_df.columns for s in args.ranch_order):
            sites = args.ranch_order
        elif set(default_baja).issubset(set(site_strengths_df.columns)):
            sites = default_baja
        else:
            sites = list(site_strengths_df.columns)

        for s in sites:
            ser = site_strengths_df[s].copy()
            ser.index = taxa_df[taxon_col]
            ser = ser.replace([np.inf, -np.inf], np.nan).dropna()
            top = ser.sort_values(ascending=False).head(args.top_per_site)
            for g, v in top.items():
                rows.append({"Site": s, "Genus": str(g), "Strength": float(v)})

    # ==================== SAMPLE MODE ==========================
    else:
        if args.meta_file is None:
            raise ValueError("Sample mode requires --meta_file.")
        comm = comm_raw.copy()
        if "sample_id" not in comm.columns:
            comm = comm.rename(columns={comm.columns[0]: "sample_id"})
        comm["sample_id"] = norm_ids(comm["sample_id"])
        taxa_cols = [c for c in comm.columns if c != "sample_id"]
        comm[taxa_cols] = comm[taxa_cols].apply(pd.to_numeric, errors="coerce").fillna(0.0)

        meta = pd.read_csv(args.meta_file)
        id_col = detect_col(meta, ["sample_id","SampleID","sample","id","ID"], required=True, ctx="metadata")
        meta[id_col] = norm_ids(meta[id_col])
        if args.group_col not in meta.columns:
            raise ValueError(f"[metadata] group_col '{args.group_col}' not found.")

        ov = np.intersect1d(comm["sample_id"].values, meta[id_col].values)
        if len(ov) < 2:
            raise ValueError("Very few overlapping sample IDs. Use --wide_sites yes for site-wise input.")
        comm = comm.set_index("sample_id").loc[ov]
        meta = meta.set_index(id_col).loc[ov]
        groups = meta[args.group_col].astype(str)
        mean_by_site = comm.groupby(groups).mean(numeric_only=True)

        if args.metric == "relative":
            totals = mean_by_site.sum(axis=1).replace(0, np.nan)
            strength = mean_by_site.divide(totals, axis=0) * 100.0
            size_label_unit = "%"
        else:
            strength = mean_by_site.copy()
            size_label_unit = ""

        default_baja = ["La_Venada", "El_Mesquite", "Casa_Vieja", "El_Tule"]
        if args.ranch_order and all(s in strength.index for s in args.ranch_order):
            sites = args.ranch_order
        elif set(default_baja).issubset(set(strength.index)):
            sites = default_baja
        else:
            sites = list(strength.index)

        for s in sites:
            if s not in strength.index:
                continue
            ser = strength.loc[s].copy().replace([np.inf, -np.inf], np.nan).dropna()
            top = ser.sort_values(ascending=False).head(args.top_per_site)
            for g, v in top.items():
                rows.append({"Site": s, "Genus": g, "Strength": float(v)})

    if not rows:
        raise ValueError("No genera selected — check inputs and --top_per_site.")

    df = pd.DataFrame(rows)

    # ======== ORIENTATION & PLOT ==========
    sites = list(dict.fromkeys(df["Site"].tolist()))
    if args.ranch_order and all(s in sites for s in args.ranch_order):
        sites = args.ranch_order
    x_map = {s: i for i, s in enumerate(sites)}

    if args.orient_elevation == "yes":
        site_rank = {s: i for i, s in enumerate(sites)}
        df["site_rank"] = df["Site"].map(site_rank)
        df = df.sort_values(["site_rank", "Strength"], ascending=[False, False]).reset_index(drop=True)
        y_labels = df["Genus"].tolist()
        y_site = df["Site"].tolist()
        sep_positions = [i - 0.5 for i in range(1, len(df)) if df.loc[i, "site_rank"] != df.loc[i - 1, "site_rank"]]
    else:
        df = df.sort_values(["Site", "Strength"], ascending=[True, False]).reset_index(drop=True)
        y_labels, y_site, sep_positions = [], [], []
        running = 0
        for s in sites[:-1]:
            sub = df[df["Site"] == s]
            y_labels.extend(sub["Genus"].tolist())
            y_site.extend(sub["Site"].tolist())
            running += len(sub)
            sep_positions.append(running - 0.5)
        y_labels.extend(df[df["Site"] == sites[-1]]["Genus"].tolist())
        y_site.extend(df[df["Site"] == sites[-1]]["Site"].tolist())

    xvals = [x_map[s] for s in y_site]
    yvals = np.arange(len(y_labels))
    sizes = scale_sizes(df.loc[:len(y_labels)-1, "Strength"].values, min_size=30, max_size=1000)
    palette = ["#3B7EA1", "#F39C12", "#27AE60", "#8E44AD",
               "#E74C3C", "#16A085", "#7F8C8D", "#2C3E50"]
    color_map = {s: palette[i % len(palette)] for i, s in enumerate(sites)}
    colors = [color_map[s] for s in y_site]

    fig_h = max(6, 0.35 * len(yvals))
    fig, ax = plt.subplots(figsize=(10, fig_h))
    ax.scatter(xvals, yvals, s=sizes, c=colors, alpha=0.85, edgecolor="k", linewidth=0.4)
    ax.set_xticks(list(x_map.values()))
    ax.set_xticklabels(list(x_map.keys()), rotation=25, ha="right")
    ax.set_yticks(yvals)
    ax.set_yticklabels(y_labels)
    ax.set_xlabel("Site")
    ax.set_ylabel("Genus")
    mode_label = "elevation-oriented" if args.orient_elevation == "yes" else "site blocks"
    ax.set_title(f"Top {args.top_per_site} genera per site (dot size = {args.metric} strength) — {mode_label}")

    for y in sep_positions:
        ax.axhline(y, color="#bbbbbb", lw=0.6, ls="--", alpha=0.6)

    ref = np.nanpercentile(df["Strength"], [30, 60, 90])
    ref_sizes = scale_sizes(ref, min_size=30, max_size=1000)
    for rs, rv in zip(ref_sizes, ref):
        ax.scatter([], [], s=rs, c="#777777", alpha=0.85, edgecolor="k", linewidth=0.4,
                   label=f"{rv:.1f}{size_label_unit}")
    ax.legend(title="Size ~ strength", bbox_to_anchor=(1.02, 1), loc="upper left", frameon=False)
    fig.tight_layout()

    # --- Create output folder automatically ---
    outdir = os.path.dirname(args.out)
    if outdir:
        os.makedirs(outdir, exist_ok=True)

    png = f"{args.out}.png"
    pdf = f"{args.out}.pdf"
    svg = f"{args.out}.svg"
    fig.savefig(png, dpi=args.dpi, bbox_inches="tight")
    fig.savefig(pdf, bbox_inches="tight")
    fig.savefig(svg, bbox_inches="tight")
    plt.close(fig)
    print(f"[OK] Wrote: {png}, {pdf}, {svg}")

    df.rename(columns={"Site": "site", "Genus": "genus", "Strength": "strength"}).to_csv(f"{args.out}_table.csv", index=False)

if __name__ == "__main__":
    main()