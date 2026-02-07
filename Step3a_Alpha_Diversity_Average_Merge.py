#!/usr/bin/env python3
# Step3a_Alpha_Diversity_Average_Merge.py
# Alpha diversity at trap or collapsed-to-ranch level with clear labels, linear + spline fits,
# and dynamic stats-box placement to avoid covering points.

import os, argparse, textwrap
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy import stats
from scipy.interpolate import UnivariateSpline

# ----------------------------- basic utils ---------------------------------

def ensure_outdir(path: str):
    os.makedirs(path, exist_ok=True)

def detect_col(df: pd.DataFrame, candidates):
    for c in candidates:
        if c in df.columns:
            return c
    lower = {c.lower(): c for c in df.columns}
    for c in candidates:
        if c.lower() in lower:
            return lower[c.lower()]
    raise ValueError(f"Missing any of: {candidates}")

# ----------------------------- I/O helpers ---------------------------------

def _autosep(path: str) -> str:
    with open(path, "r", encoding="utf-8") as f:
        first = f.readline()
    return "\t" if first.count("\t") > first.count(",") else ","

def read_counts_table(path: str) -> pd.DataFrame:
    """Read taxa x samples. Use first column as taxa ID if present."""
    sep = _autosep(path)
    counts = pd.read_csv(path, sep=sep)
    taxon_col = None
    for c in ["Taxon", "taxon", "Genus", "ASV", "OTU", counts.columns[0]]:
        if c in counts.columns:
            taxon_col = c
            break
    if taxon_col is None:
        raise ValueError("No taxon column detected in counts table.")
    counts = counts.set_index(taxon_col)
    counts = counts.apply(pd.to_numeric, errors="coerce").fillna(0)
    return counts

def read_metadata(path: str) -> pd.DataFrame:
    """Standardize to SampleID, ranch, elev_m; keep all cols (for weights)."""
    meta = pd.read_excel(path) if path.lower().endswith((".xlsx", ".xls")) else pd.read_csv(path)
    sample_col = detect_col(meta, ["SampleID", "sample_id", "sample", "Sample"])
    ranch_col  = detect_col(meta, ["ranch", "Ranch"])
    elev_col   = detect_col(meta, ["elev_m", "Elevation", "elevation", "elev"])
    meta = meta.rename(columns={sample_col: "SampleID", ranch_col: "ranch", elev_col: "elev_m"})
    meta["SampleID"] = meta["SampleID"].astype(str)
    meta["ranch"]    = meta["ranch"].astype(str)
    return meta

# ----------------------------- collapsing helpers ---------------------------------

def _geometric_mean(arr, eps=1e-12):
    arr = np.asarray(arr, float) + eps
    return float(np.exp(np.mean(np.log(arr))))

def _harmonic_mean(arr, eps=1e-12):
    arr = np.asarray(arr, float) + eps
    return float(len(arr) / np.sum(1.0 / arr))

def collapse_by_mean(counts: pd.DataFrame,
                     meta: pd.DataFrame,
                     method: str = "arithmetic",
                     weight_col: str | None = None) -> pd.DataFrame:
    """
    Collapse duplicate traps to one column per ranch using:
      arithmetic (default for MEAN_REL), weighted (needs --weight-column),
      median (robust), geometric/harmonic (rare here).
    Returns taxa x ranch.
    """
    sample_ids = [c for c in counts.columns if c in set(meta["SampleID"])]
    if not sample_ids:
        raise ValueError("No overlapping SampleID between counts and metadata.")
    sub = counts[sample_ids].copy()

    # columns: MultiIndex (ranch, sample)
    ranch_for_sample = meta.set_index("SampleID")["ranch"].loc[sample_ids]
    sub.columns = pd.MultiIndex.from_arrays([ranch_for_sample.values, sub.columns],
                                            names=["ranch", "sample"])

    blocks = []
    # future-proof: group on transposed columns
    for ranch, blockT in sub.T.groupby(level="ranch"):
        block = blockT.T  # taxa x samples (replicates of this ranch)
        if method == "arithmetic":
            agg = block.mean(axis=1)
        elif method == "median":
            agg = block.median(axis=1)
        elif method == "geometric":
            agg = block.apply(_geometric_mean, axis=1)
        elif method == "harmonic":
            agg = block.apply(_harmonic_mean, axis=1)
        elif method == "weighted":
            if weight_col is None:
                raise ValueError("weighted mean requested but --weight-column was not provided.")
            samp = block.columns.get_level_values("sample").tolist()
            w = (meta.set_index("SampleID")[weight_col]
                     .astype(float)
                     .reindex(samp))
            w = w.where(np.isfinite(w), 0.0).to_numpy()
            w = w / w.sum() if w.sum() > 0 else np.ones_like(w) / len(w)
            agg = (block * w).sum(axis=1)
        else:
            raise ValueError(f"Unknown collapse-mean: {method}")
        agg.name = ranch
        blocks.append(agg)

    return pd.concat(blocks, axis=1)

# ----------------------------- alpha diversity ---------------------------------

def alpha_metrics_from_counts(counts: pd.DataFrame) -> pd.DataFrame:
    """Compute richness, Shannon, Simpson for each column."""
    totals   = counts.sum(axis=0)
    richness = (counts > 0).sum(axis=0)

    def sh(col):
        s = col.sum()
        if s <= 0: return 0.0
        p = (col / s).to_numpy()
        p = p[p > 0]
        return -float(np.sum(p * np.log(p)))

    def si(col):
        s = col.sum()
        if s <= 0: return 0.0
        p = (col / s).to_numpy()
        return 1.0 - float(np.sum(p * p))

    shannon = counts.apply(sh, axis=0)
    simpson = counts.apply(si, axis=0)

    return pd.DataFrame({
        "SampleID": counts.columns,
        "richness": richness.values,
        "shannon":  shannon.values,
        "simpson":  simpson.values,
        "total_count": totals.values
    })

# ----------------------------- robust fitting & plotting ---------------------------------

def _safe_lm_stats(x, y):
    """Linear fit stats; return NaNs if not enough data or x has no spread."""
    x = np.asarray(x, float); y = np.asarray(y, float)
    m = np.isfinite(x) & np.isfinite(y)
    x, y = x[m], y[m]
    if len(x) < 2 or np.ptp(x) == 0:
        return {"slope": np.nan, "intercept": np.nan, "r2": np.nan, "p": np.nan, "stderr": np.nan}
    res = stats.linregress(x, y)
    return {"slope": res.slope, "intercept": res.intercept, "r2": res.rvalue**2,
            "p": res.pvalue, "stderr": res.stderr}

def _place_stats_box(ax, x, y, text):
    """
    Put the stats box in the upper corner farthest from the highest point.
    If max-y point is on the left half -> place top-right; else top-left.
    """
    x = np.asarray(x, float); y = np.asarray(y, float)
    if len(x) == 0:
        corner = ("left", (0.02, 0.98))
    else:
        i_top = int(np.argmax(y))
        x_med = np.median(x)
        corner = ("right", (0.98, 0.98)) if x[i_top] < x_med else ("left", (0.02, 0.98))
    ha = corner[0]; (ax_x, ax_y) = corner[1]
    ax.text(ax_x, ax_y, text,
            transform=ax.transAxes, va="top", ha=ha,
            bbox=dict(boxstyle="round", facecolor="white", alpha=0.85, edgecolor="0.3"))

def add_linear_and_spline(ax, x, y, label_points=None):
    """
    Big outlined points + bold labels, linear fit, and a visible dashed spline.
    For n<5, spline uses k=2, s=0 so you always see a curve.
    """
    x = np.asarray(x, float)
    y = np.asarray(y, float)
    mask = np.isfinite(x) & np.isfinite(y)
    x, y = x[mask], y[mask]

    labels = None
    if label_points is not None:
        lp = np.asarray(label_points, dtype=object)
        labels = lp[mask]

    # Larger points with outline (avoid hiding under the line)
    ax.scatter(x, y, s=80, zorder=3, color="steelblue", edgecolor="black")

    if len(x) >= 2 and np.ptp(x) > 0:
        # Linear fit (solid blue)
        try:
            coeff = np.polyfit(x, y, 1)
            xx = np.linspace(x.min(), x.max(), 200)
            yy = np.polyval(coeff, xx)
            ax.plot(xx, yy, lw=2, color="#1f77b4", label="Linear fit")
        except Exception:
            pass

        # Spline fit (dashed orange); force visibility for tiny n
        try:
            k = 2 if len(x) < 5 else min(3, len(x)-1)
            s = 0.0 if len(x) < 5 else max(0.0, 0.1 * len(x))
            spl = UnivariateSpline(x, y, k=k, s=s)
            xx = np.linspace(x.min(), x.max(), 200)
            ax.plot(xx, spl(xx), ls="--", lw=2, color="orange", label=f"Spline fit (k={k})")
        except Exception:
            # Fallback: dashed quadratic if spline fails
            if len(x) >= 3:
                try:
                    p2 = np.polyfit(x, y, 2)
                    xx = np.linspace(x.min(), x.max(), 200)
                    ax.plot(xx, np.polyval(p2, xx), ls="--", lw=2, color="orange", label="Quadratic fit")
                except Exception:
                    pass

    # Bold ranch labels
    if labels is not None:
        for xi, yi, lbl in zip(x, y, labels):
            ax.annotate(str(lbl), (xi, yi),
                        textcoords="offset points", xytext=(8, 6),
                        fontsize=11, fontweight="bold", color="darkblue")

# ----------------------------- CLI ---------------------------------

def get_parser():
    p = argparse.ArgumentParser(
        description="Step 3a: Alpha diversity with optional replicate collapse and clear figures",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=textwrap.dedent("""Outputs:
  <outdir>/alpha_diversity.csv
  <outdir>/alpha_summary_by_ranch.csv
  <outdir>/alpha_stats.csv
  <outdir>/figures/*.png""")
    )
    p.add_argument("--counts", required=True, help="Taxa x samples counts CSV/TSV")
    p.add_argument("--meta",   required=True, help="Aligned metadata with SampleID, ranch, elev_m")
    p.add_argument("--outdir", required=True, help="Output directory")
    p.add_argument("--collapse-to-ranch", action="store_true",
                   help="Collapse duplicate traps to a single ranch column before alpha")
    p.add_argument("--collapse-mean",
                   choices=["arithmetic", "weighted", "median", "geometric", "harmonic"],
                   default="arithmetic",
                   help="Averaging method for replicate collapse (default: arithmetic)")
    p.add_argument("--weight-column", default=None,
                   help="Metadata column to weight by when --collapse-mean weighted")
    return p

# ----------------------------- main ---------------------------------

def main():
    args = get_parser().parse_args()
    ensure_outdir(args.outdir)
    figdir = os.path.join(args.outdir, "figures")
    ensure_outdir(figdir)

    counts = read_counts_table(args.counts)
    meta   = read_metadata(args.meta)

    if args.collapse_to_ranch:
        counts = collapse_by_mean(counts, meta, method=args.collapse_mean, weight_col=args.weight_column)

    alpha = alpha_metrics_from_counts(counts)
    df = alpha.merge(meta, on="SampleID", how="left")

    # If collapsed, SampleID == ranch; map elevation by ranch mean
    if args.collapse_to_ranch:
        df["ranch"] = df["SampleID"]
        if "elev_m" in meta.columns:
            elev_by_ranch = meta.groupby("ranch", dropna=False)["elev_m"].mean()
            df["elev_m"] = df["ranch"].map(elev_by_ranch)

    if df.empty:
        raise ValueError("No overlap of SampleID between counts and metadata.")

    out_csv = os.path.join(args.outdir, "alpha_diversity.csv")
    df.sort_values("SampleID").to_csv(out_csv, index=False)

    ranch_summary = (df.groupby("ranch")[["richness","shannon","simpson","total_count"]]
                       .mean().reset_index().sort_values("ranch"))
    ranch_summary.to_csv(os.path.join(args.outdir, "alpha_summary_by_ranch.csv"), index=False)

    # --- Bar figures ---
    plt.figure(figsize=(7,5))
    plt.bar(ranch_summary["ranch"], ranch_summary["richness"], color="#6baed6")
    plt.ylabel("Mean richness")
    plt.title("Mean richness by ranch")
    plt.xticks(rotation=25)
    plt.tight_layout()
    plt.savefig(os.path.join(figdir, "alpha_richness_by_ranch.png"), dpi=200)
    plt.close()

    plt.figure(figsize=(11,5))
    order = df.sort_values("richness", ascending=False)
    plt.bar(order["SampleID"], order["richness"], color="#9ecae1")
    plt.ylabel("Richness (taxa > 0)")
    plt.title("Per-sample richness" + (" (ranch level)" if args.collapse_to_ranch else " (trap level)"))
    plt.xticks(rotation=35, ha="right")
    plt.tight_layout()
    plt.savefig(os.path.join(figdir, "alpha_richness_by_sample.png"), dpi=200)
    plt.close()

    # --- Elevation figures & stats ---
    if "elev_m" in df.columns:
        df_elev = df[np.isfinite(df["elev_m"].to_numpy(dtype=float))].copy()
    else:
        df_elev = pd.DataFrame()

    if not df_elev.empty and df_elev.shape[0] >= 2:
        # 3-panel with big points + labels
        fig, axes = plt.subplots(1,3, figsize=(18,6))
        for ax, metric, label in zip(
            axes,
            ["richness","shannon","simpson"],
            ["Richness","Shannon","Simpson"]
        ):
            ax.scatter(df_elev["elev_m"], df_elev[metric], s=80, color="steelblue", edgecolor="black")
            for _, row in df_elev.iterrows():
                ax.annotate(str(row["ranch"]), (row["elev_m"], row[metric]),
                            textcoords="offset points", xytext=(8,6),
                            fontsize=11, fontweight="bold", color="darkblue")
            ax.set_xlabel("Elevation (m)")
            ax.set_ylabel(label)
            ax.set_title(f"{label} vs Elevation")
        plt.tight_layout()
        plt.savefig(os.path.join(figdir, "alpha_vs_elevation.png"), dpi=200)
        plt.close()

        # Single-metric panels with fits + dynamic stats box
        stats_rows = []
        def one_metric(metric, fname):
            x = df_elev["elev_m"].to_numpy(float)
            y = df_elev[metric].to_numpy(float)
            sdict = _safe_lm_stats(x, y)
            fig, ax = plt.subplots(figsize=(10,7))
            add_linear_and_spline(ax, x, y, label_points=df_elev["ranch"].tolist())
            ax.set_xlabel("elev_m"); ax.set_ylabel(metric); ax.legend(frameon=False)
            txt = (f"slope={sdict['slope'] if pd.notna(sdict['slope']) else 'NA'}\n"
                   f"R²={sdict['r2'] if pd.notna(sdict['r2']) else 'NA'}\n"
                   f"p={sdict['p'] if pd.notna(sdict['p']) else 'NA'}")
            _place_stats_box(ax, x, y, txt)   # <— dynamic placement
            ax.set_title(f"{metric} vs elev_m")
            plt.tight_layout()
            plt.savefig(os.path.join(figdir, fname), dpi=200)
            plt.close()
            stats_rows.append({"metric": metric, "predictor": "elev_m", **sdict})

        one_metric("richness", "richness_vs_elev.png")
        one_metric("shannon",  "shannon_vs_elev.png")
        one_metric("simpson",  "simpson_vs_elev.png")
        pd.DataFrame(stats_rows).to_csv(os.path.join(args.outdir, "alpha_stats.csv"), index=False)

    print(f"[OK] Wrote {out_csv}")
    print(f"[OK] Wrote {os.path.join(args.outdir, 'alpha_summary_by_ranch.csv')}")
    if os.path.exists(os.path.join(args.outdir, 'alpha_stats.csv')):
        print(f"[OK] Wrote {os.path.join(args.outdir, 'alpha_stats.csv')}")
    print(f"[OK] Wrote figures to: {figdir}")
    print("[DONE] Step 3a with optional replicate collapsing.")

if __name__ == "__main__":
    main()