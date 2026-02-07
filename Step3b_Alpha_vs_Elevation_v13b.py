#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# CHANGES vs v13b:
# 1) Normalize SampleID from any of: sampleid, sample_id, sample, id  → SampleID
#    (This fixes KeyError: "['SampleID'] not in index".)
# 2) If no ID column exists, synthesize SampleID = "{ranch}_{rownum}" so plotting still runs.
# 3) Keep the “read Step3a alpha_diversity.csv directly” behavior.

import os, argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy import stats

RANCH_ORDER = ["La Venada", "El Mesquite", "Casa Vieja", "El Tule"]
ELEV_ALIASES = [
    "elev_m","elevation","Elevation","elev","ElevDEM","elev_meters",
    "elev_meters_mean","elev_mean","elev_m_mean","elevation_m","Elevation_m"
]

def ensure_dir(p): os.makedirs(p, exist_ok=True)

def _auto_sep_read(path: str) -> pd.DataFrame:
    with open(path, "r", encoding="utf-8", errors="ignore") as f:
        head = f.read(2048)
    sep = "\t" if head.count("\t") > head.count(",") else ","
    return pd.read_csv(path, sep=sep)

def _norm_col(c: str) -> str:
    c = str(c).replace("\u00a0", " ")
    c = " ".join(c.split()).strip()
    return c.lower().replace(" ", "_").replace("-", "_")

def _norm_df_cols(df: pd.DataFrame) -> pd.DataFrame:
    df = df.copy()
    df.columns = [_norm_col(c) for c in df.columns]
    return df

def normalize_ranch(x) -> str:
    if not isinstance(x, str): return x
    t = x.strip().replace("_"," "); low = t.lower()
    if low in {"el mezquital","mezquital","el mesquite","mesquite"}: return "El Mesquite"
    if low in {"la venada","lavenada"}: return "La Venada"
    if low in {"casa vieja","casavieja","casa viejas","casa_vieja"}: return "Casa Vieja"
    if low in {"el tule","eltule"}: return "El Tule"
    return " ".join(w.capitalize() for w in t.split())

def coerce_elevation(df: pd.DataFrame) -> pd.DataFrame:
    df = df.copy(); cols = set(df.columns); elev_col = None
    for name in ELEV_ALIASES:
        key = _norm_col(name)
        if key in cols:
            elev_col = key; break
    if elev_col is None:
        raise ValueError("No elevation column found. Looked for: " + ", ".join(ELEV_ALIASES))
    df.rename(columns={elev_col:"elev_m"}, inplace=True)
    return df

# ---------- stats/plot (unchanged) ----------
def linear_fit_stats(x, y):
    if np.allclose(np.std(x),0) or len(x) < 2: return (np.nan,np.nan,np.nan,np.nan,np.nan)
    slope, intercept, r, p, stderr = stats.linregress(x, y)
    return slope, intercept, r**2, p, stderr

def ci_band(ax, x, y, slope, intercept):
    n = len(x)
    if n < 3 or np.isnan(slope): return
    x = np.asarray(x); y = np.asarray(y)
    yhat = intercept + slope*x
    resid = y - yhat
    s_err = np.sqrt(np.sum(resid**2)/(n-2))
    x0 = np.linspace(np.min(x), np.max(x), 200)
    xbar = np.mean(x)
    from scipy.stats import t as _t
    se = s_err * np.sqrt(1/n + (x0 - xbar)**2 / np.sum((x - xbar)**2))
    y0 = intercept + slope*x0
    tval = _t.ppf(0.975, df=n-2)
    ax.fill_between(x0, y0 - tval*se, y0 + tval*se, alpha=0.15, linewidth=0, label="95% CI")

def add_quadratic(ax, x, y, xgrid):
    if len(x) < 3: return
    coeff = np.polyfit(x, y, 2)
    ax.plot(xgrid, np.polyval(coeff, xgrid), linestyle="--", linewidth=2.0, color="orange", label="quadratic")

def _scatter_panel(df, xcol, ycol, out_png, title):
    D = df.dropna(subset=[xcol, ycol]).copy()
    if D.empty: return
    fig, ax = plt.subplots(figsize=(6.2, 4.6), dpi=180)
    colors = {"La Venada":"#1f77b4","El Mesquite":"#2ca02c","Casa Vieja":"#9467bd","El Tule":"#d62728"}
    if "ranch" in D.columns:
        for rn, sub in D.groupby("ranch", sort=False):
            ax.scatter(sub[xcol], sub[ycol], s=55, alpha=0.85,
                       edgecolor="white", linewidth=0.6,
                       label=str(rn), color=colors.get(str(rn), "#333333"))
    else:
        ax.scatter(D[xcol], D[ycol], s=60, alpha=0.9, edgecolor="white", linewidth=0.6)
    slope, intercept, r2, p, stderr = linear_fit_stats(D[xcol].values, D[ycol].values)
    xg = np.linspace(D[xcol].min(), D[xcol].max(), 200)
    if not np.isnan(slope):
        ax.plot(xg, intercept + slope*xg, linewidth=2.2, label="linear")
        ci_band(ax, D[xcol].values, D[ycol].values, slope, intercept)
    add_quadratic(ax, D[xcol].values, D[ycol].values, xg)
    ax.set_title(title, pad=10); ax.set_xlabel(xcol.replace("_"," ").title()); ax.set_ylabel(ycol.title())
    ax.legend(frameon=False, fontsize=9, loc="best")
    ax.text(0.98, 0.02, f"N={len(D)}\nslope={slope:.3g}\nR²={r2:.3f}\np={p:.3g}",
            ha="right", va="bottom", transform=ax.transAxes,
            fontsize=9, bbox=dict(boxstyle="round,pad=0.3", fc="#f7f7f7", ec="#ccc"))
    fig.tight_layout(); fig.savefig(out_png); plt.close(fig)

def write_model_row(rows, metric, xcol, D):
    slope, intercept, r2, p, stderr = linear_fit_stats(D[xcol].values, D[metric].values)
    rows.append(dict(metric=metric, predictor=xcol, N=len(D),
                     slope=slope, intercept=intercept, r2=r2, p=p, stderr=stderr))

def make_level_outputs(level_name, df_in, outdir):
    ensure_dir(outdir)
    df_in.to_csv(os.path.join(outdir, f"{level_name}_data_used.csv"), index=False)
    model_rows = []
    predictors = [c for c in ["elev_m","cum_distance_miles","cum_gain_ft"] if c in df_in.columns]
    metrics    = [c for c in ["richness","shannon","simpson"] if c in df_in.columns]
    for m in metrics:
        for x in predictors:
            _scatter_panel(df_in, x, m, os.path.join(outdir, f"{m}_vs_{x}.png"),
                           f"{m.title()} vs {x.replace('_',' ').title()} ({level_name})")
            D = df_in.dropna(subset=[x, m])
            if not D.empty: write_model_row(model_rows, m, x, D)
    pd.DataFrame(model_rows).to_csv(os.path.join(outdir, f"{level_name}_model_summaries.csv"), index=False)

# ---------- NEW: robust ID adoption ----------
def adopt_sample_id(T: pd.DataFrame) -> pd.DataFrame:
    T = T.copy()
    id_candidates = ["sampleid","sample_id","sample","id"]
    found = None
    for c in id_candidates:
        if c in T.columns:
            found = c; break
    if found:
        T.rename(columns={found:"SampleID"}, inplace=True)
    else:
        # synthesize SampleID if none present
        base = T["ranch"].astype(str).fillna("NA") if "ranch" in T.columns else pd.Series(["row"]*len(T))
        T["SampleID"] = [f"{r}_{i+1}" for i, r in enumerate(base)]
    return T

def table_has_ranch_and_elev(T: pd.DataFrame) -> bool:
    cols = set(T.columns)
    has_ranch = ("ranch" in cols) or ("site" in cols) or ("location" in cols)
    has_elev  = any((_norm_col(c) in {_norm_col(e) for e in ELEV_ALIASES}) for c in cols)
    return has_ranch and has_elev

def adopt_ranch_and_elev(T: pd.DataFrame) -> pd.DataFrame:
    T = T.copy()
    if "ranch" not in T.columns:
        for c in ["site","location"]:
            if c in T.columns:
                T.rename(columns={c:"ranch"}, inplace=True); break
    T["ranch"] = T["ranch"].map(normalize_ranch)
    T = coerce_elevation(T)
    T["ranch"] = pd.Categorical(T["ranch"], RANCH_ORDER, ordered=True)
    return T

def build_from_table(table_path):
    T = _norm_df_cols(_auto_sep_read(table_path))
    # metrics required
    if not {"richness","shannon","simpson"} <= set(T.columns):
        raise ValueError("Table missing metrics: richness, shannon, simpson")
    # ranch + elevation (from Step3a table)
    if not table_has_ranch_and_elev(T):
        raise ValueError("Table lacks ranch/elevation. Use --from-step3a-dir or pass a Step3a alpha_diversity.csv.")
    T = adopt_ranch_and_elev(T)
    # normalize SampleID no matter how it’s named
    T = adopt_sample_id(T)
    keep = ["SampleID","ranch","elev_m","richness","shannon","simpson"]
    for c in ["cum_distance_miles","cum_gain_ft"]:
        if c in T.columns: keep.append(c)
    return T[keep].sort_values(["ranch","SampleID"])

def main():
    ap = argparse.ArgumentParser(description="Alpha vs elevation (reads Step3a output directly; robust SampleID).")
    g = ap.add_mutually_exclusive_group(required=True)
    g.add_argument("--table", help="Single table (e.g., Step3a alpha_diversity.csv) that already includes ranch + elevation.")
    g.add_argument("--from-step3a-dir", help="Directory containing Step3a/alpha_diversity.csv")
    g.add_argument("--auto-merge", action="store_true", help="Legacy mode (requires --alpha + --meta)")
    ap.add_argument("--alpha", default=None); ap.add_argument("--meta", default=None)
    ap.add_argument("--outdir", required=True)
    args = ap.parse_args()

    ensure_dir(args.outdir)
    rep_dir = os.path.join(args.outdir, "replicates"); ensure_dir(rep_dir)
    mean_dir = os.path.join(args.outdir, "means");      ensure_dir(mean_dir)

    if args.from_step3a_dir:
        table_path = os.path.join(args.from_step3a_dir, "alpha_diversity.csv")
        if not os.path.isfile(table_path):
            raise SystemExit(f"Could not find Step3a alpha_diversity.csv at: {table_path}")
        df = build_from_table(table_path)
    elif args.table:
        df = build_from_table(args.table)
    else:
        # legacy path (still works, but you likely won’t need it)
        if not args.alpha or not args.meta:
            raise SystemExit("--auto-merge requires --alpha and --meta")
        A = _norm_df_cols(_auto_sep_read(args.alpha))
        M = _norm_df_cols(_auto_sep_read(args.meta))
        if "sampleid" not in A.columns and "sample_id" not in A.columns:
            raise ValueError("alpha file must contain a sample ID column (SampleID / sampleid / sample_id).")
        if "sampleid" in A.columns: A.rename(columns={"sampleid":"SampleID"}, inplace=True)
        if "sample_id" in A.columns: A.rename(columns={"sample_id":"SampleID"}, inplace=True)
        if not {"richness","shannon","simpson"} <= set(A.columns):
            raise ValueError("alpha missing metrics: richness, shannon, simpson")
        ranch_col = None
        for c in ["ranch","site","location"]:
            if c in M.columns: ranch_col = c; break
        if ranch_col is None: raise ValueError("metadata must have ranch/site/location")
        M.rename(columns={ranch_col:"ranch"}, inplace=True)
        M["ranch"] = M["ranch"].map(normalize_ranch)
        M = coerce_elevation(M)
        df = pd.merge(A, M[["SampleID","ranch","elev_m"]], on="SampleID", how="left")
        df["ranch"] = pd.Categorical(df["ranch"], RANCH_ORDER, ordered=True)
        df = df.sort_values(["ranch","SampleID"])

    # replicate-level outputs
    make_level_outputs("replicates", df, rep_dir)
    # ranch means
    means = (df[["ranch","elev_m","richness","shannon","simpson"]].groupby("ranch", as_index=False).mean(numeric_only=True))
    means["ranch"] = pd.Categorical(means["ranch"], RANCH_ORDER, ordered=True)
    means = means.sort_values("ranch")
    make_level_outputs("means", means, mean_dir)
    print("[OK] Wrote replicate & ranch-mean plots and CSVs to:", args.outdir)

if __name__ == "__main__":
    main()