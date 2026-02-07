#!/usr/bin/env python3
"""
Step4_beta_partition_v6.py
Compute beta diversity at sample- and ranch-level with preserved ranch order
AND save ordered heatmaps.

- Presence/absence -> Sørensen (β_sor), Simpson turnover (β_sim), Nestedness (β_sne = β_sor - β_sim)
  using Baselga (2010) definitions from a,b,c counts.
- Hellinger matrix -> Bray–Curtis dissimilarity.
- Outputs both SAMPLE-level and RANCH-level square matrices (CSV), sorted by user ranch order.
- Also writes PNG heatmaps with that order.

Inputs
  --presence : CSV, taxa x samples (first column = taxa)
  --hellinger: CSV, taxa x samples (first column = taxa), Hellinger-transformed abundances
  --metadata : CSV with columns [SampleID, ranch, ...]
  --alias-map: JSON mapping raw->std names (optional)
  --ranch-order: comma-separated ranch order (after aliasing), e.g. "La_Venada,El_Mezquital,Casa_Vieja,El_Tule"
"""

import argparse, json
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.spatial.distance import pdist, squareform
import seaborn as sns
import matplotlib.pyplot as plt

# ---------- IO helpers ----------

def read_matrix_taxa_by_samples(path: str) -> pd.DataFrame:
    """Read a taxa x samples table with first column as taxa index."""
    df = pd.read_csv(path, index_col=0)
    df.columns = [str(c).strip() for c in df.columns]
    return df

def apply_alias_and_standardize_ranch(meta: pd.DataFrame, alias_map: dict) -> pd.DataFrame:
    meta = meta.copy()
    meta.columns = [c.strip() for c in meta.columns]
    # flexible renaming
    rename = {
        "sample_id": "SampleID", "Sample_id": "SampleID",
        "sampleID": "SampleID", "sampleId": "SampleID",
        "Latitude":"lat", "Longitude":"lon", "Longitiude":"lon",
        "Elevation":"elev_m", "elevation":"elev_m", "elevation_m":"elev_m"
    }
    meta = meta.rename(columns={c: rename.get(c, c) for c in meta.columns})
    if "SampleID" not in meta.columns or "ranch" not in meta.columns:
        raise SystemExit(f"[ERROR] metadata must have columns ['SampleID','ranch']; got {list(meta.columns)}")

    # alias mapping (accept keys w/ spaces or underscores)
    if alias_map:
        inv_alias = {k.replace("_"," "): v for k,v in alias_map.items()}
        meta["ranch"] = meta["ranch"].map(lambda x: alias_map.get(x, inv_alias.get(x, x)))
    meta["ranch"] = meta["ranch"].astype(str).str.strip().str.replace(" ", "_", regex=False)
    return meta

# ---------- beta computations ----------

def compute_baselga_components(pres_bool_samples_x_taxa: np.ndarray):
    """
    Baselga (2010) pairwise components from presence/absence:
      a = shared presences
      b = presences unique to i
      c = presences unique to j
      β_sor = (b + c) / (2a + b + c)
      β_sim = min(b, c) / (a + min(b, c))
      β_sne = β_sor - β_sim
    Input shape: (n_samples, n_taxa) as bool/int
    """
    B = pres_bool_samples_x_taxa.astype(bool)
    n = B.shape[0]

    beta_sor = np.zeros((n, n), float)
    beta_sim = np.zeros((n, n), float)

    for i in range(n):
        Bi = B[i]
        for j in range(i+1, n):
            Bj = B[j]
            a = np.logical_and(Bi, Bj).sum()
            b = np.logical_and(Bi, ~Bj).sum()
            c = np.logical_and(~Bi, Bj).sum()

            denom_sor = (2*a + b + c)
            bs = 0.0 if denom_sor == 0 else (b + c) / denom_sor

            m = min(b, c)
            denom_sim = (a + m)
            bt = 0.0 if denom_sim == 0 else m / denom_sim

            beta_sor[i, j] = beta_sor[j, i] = bs
            beta_sim[i, j] = beta_sim[j, i] = bt

    beta_sne = beta_sor - beta_sim
    return beta_sor, beta_sim, beta_sne

def bray_curtis_from_hellinger(hell_taxa_by_samples: pd.DataFrame, order: list[str]) -> pd.DataFrame:
    """Bray–Curtis on Hellinger matrix (taxa x samples). Output square matrix indexed by given order."""
    use = [s for s in order if s in hell_taxa_by_samples.columns]
    X = hell_taxa_by_samples[use].T.values  # samples x taxa
    dm = squareform(pdist(X, metric="braycurtis"))
    return pd.DataFrame(dm, index=use, columns=use)

def order_square(df: pd.DataFrame, order: list[str]) -> pd.DataFrame:
    use = [x for x in order if x in df.index]
    return df.loc[use, use]

# ---------- ordering & plotting ----------

def build_sample_order(sample_ids, meta_df, ranch_order):
    """Return samples ordered by ranch_order, preserving within-ranch original order."""
    s2r = meta_df.set_index("SampleID")["ranch"]
    ordered = []
    # preferred order first
    for r in ranch_order:
        r_samples = [s for s in sample_ids if s in s2r.index and s2r[s] == r]
        ordered.extend(r_samples)
    # any leftovers (ranches not in ranch_order) appended after
    leftovers = [s for s in sample_ids if s not in ordered]
    ordered.extend(leftovers)
    return ordered

def heatmap(df: pd.DataFrame, outpng: Path, title: str, vmin=0, vmax=1, cmap="viridis"):
    plt.figure(figsize=(8, 6))
    sns.heatmap(df, square=True, vmin=vmin, vmax=vmax, cmap=cmap,
                cbar_kws=dict(label="Dissimilarity"))
    plt.title(title)
    plt.tight_layout()
    plt.savefig(outpng, dpi=150)
    plt.close()

# ---------- main ----------

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--presence", required=True, help="presence_absence.csv (taxa x samples)")
    ap.add_argument("--hellinger", required=True, help="hellinger_matrix.csv (taxa x samples)")
    ap.add_argument("--metadata", required=True, help="metadata CSV with SampleID,ranch")
    ap.add_argument("--outdir", required=True, help="output directory")
    ap.add_argument("--alias-map", default="{}", help="JSON mapping of raw->std ranch names")
    ap.add_argument("--ranch-order", default="",
                    help='Comma-separated ranch order after aliasing, e.g. "La_Venada,El_Mezquital,Casa_Vieja,El_Tule"')
    args = ap.parse_args()

    outdir = Path(args.outdir); outdir.mkdir(parents=True, exist_ok=True)
    alias_map = json.loads(args.alias_map) if args.alias_map else {}
    user_ranch_order = [r.strip() for r in args.ranch_order.split(",") if r.strip()]
    # normalize ranch order same way as metadata aliasing would
    user_ranch_order = [alias_map.get(r, r).replace(" ", "_") for r in user_ranch_order]

    # read inputs
    pres = read_matrix_taxa_by_samples(args.presence)
    hell = read_matrix_taxa_by_samples(args.hellinger)
    meta = pd.read_csv(args.metadata)
    meta = apply_alias_and_standardize_ranch(meta, alias_map)

    # align samples across all three inputs
    sample_cols = sorted(set(pres.columns) & set(hell.columns) & set(meta["SampleID"]))
    if not sample_cols:
        raise SystemExit("[ERROR] No overlapping samples across presence/hellinger/metadata.")
    pres = pres[sample_cols].copy()
    hell = hell[sample_cols].copy()
    meta = meta.set_index("SampleID").loc[sample_cols].reset_index()

    # ---------- SAMPLE-LEVEL ----------
    pres_bool = (pres > 0).astype(bool)        # taxa x samples (bools)
    B = pres_bool.T.values                      # samples x taxa

    beta_sor, beta_sim, beta_sne = compute_baselga_components(B)
    idx = sample_cols
    df_sor = pd.DataFrame(beta_sor, index=idx, columns=idx)
    df_sim = pd.DataFrame(beta_sim, index=idx, columns=idx)
    df_sne = pd.DataFrame(beta_sne, index=idx, columns=idx)
    df_bray = bray_curtis_from_hellinger(hell, order=sample_cols)

    # order samples by ranch order
    sample_order = build_sample_order(sample_cols, meta, user_ranch_order if user_ranch_order else list(meta["ranch"].unique()))
    df_sor_o  = order_square(df_sor,  sample_order)
    df_sim_o  = order_square(df_sim,  sample_order)
    df_sne_o  = order_square(df_sne,  sample_order)
    df_bray_o = order_square(df_bray, sample_order)

    (outdir / "samples").mkdir(exist_ok=True, parents=True)
    df_sor_o.to_csv(outdir / "samples" / "beta_sorensen_samples.csv")
    df_sim_o.to_csv(outdir / "samples" / "beta_simpson_turnover_samples.csv")
    df_sne_o.to_csv(outdir / "samples" / "beta_nestedness_samples.csv")
    df_bray_o.to_csv(outdir / "samples" / "bray_curtis_samples.csv")

    # sample heatmaps
    heatmap(df_sor_o,  outdir / "samples" / "heat_beta_sor_samples.png",  "β_Sørensen (samples)")
    heatmap(df_sim_o,  outdir / "samples" / "heat_beta_sim_samples.png",  "β_Simpson turnover (samples)")
    heatmap(df_sne_o,  outdir / "samples" / "heat_beta_sne_samples.png",  "β_Nestedness (samples)")
    heatmap(df_bray_o, outdir / "samples" / "heat_bray_curtis_samples.png","Bray–Curtis (samples)")

    # ---------- RANCH-LEVEL ----------
    s2r = meta.set_index("SampleID")["ranch"].to_dict()

    # presence OR within ranch (taxa present in any sample of that ranch)
    # (avoid deprecated axis=1 groupby)
    pres_by_ranch = (
        pres_bool.T.assign(ranch=[s2r[s] for s in pres_bool.columns])
                   .groupby("ranch")
                   .any()
                   .T
                   .astype(int)
    )
    # hellinger average within ranch
    hell_by_ranch = (
        hell.T.assign(ranch=[s2r[s] for s in hell.columns])
              .groupby("ranch")
              .mean()
              .T
    )

    ranch_cols = list(pres_by_ranch.columns)
    if user_ranch_order:
        # keep requested order, then any extras
        ranch_cols = [r for r in user_ranch_order if r in ranch_cols] + [r for r in pres_by_ranch.columns if r not in user_ranch_order]

    # compute Baselga components on ranch-level presence
    B_ranch = pres_by_ranch[ranch_cols].T.values  # ranches x taxa
    sor_r, sim_r, sne_r = compute_baselga_components(B_ranch)

    df_sor_r = pd.DataFrame(sor_r, index=ranch_cols, columns=ranch_cols)
    df_sim_r = pd.DataFrame(sim_r, index=ranch_cols, columns=ranch_cols)
    df_sne_r = pd.DataFrame(sne_r, index=ranch_cols, columns=ranch_cols)

    # Bray–Curtis on ranch-level hellinger means
    df_bray_r = bray_curtis_from_hellinger(hell_by_ranch, order=ranch_cols)

    (outdir / "ranches").mkdir(exist_ok=True, parents=True)
    df_sor_r.to_csv(outdir / "ranches" / "beta_sorensen_ranches.csv")
    df_sim_r.to_csv(outdir / "ranches" / "beta_simpson_turnover_ranches.csv")
    df_sne_r.to_csv(outdir / "ranches" / "beta_nestedness_ranches.csv")
    df_bray_r.to_csv(outdir / "ranches" / "bray_curtis_ranches.csv")

    # ranch heatmaps (ordered La_Venada -> ... -> El_Tule, etc.)
    heatmap(df_sor_r,  outdir / "ranches" / "heat_beta_sor_ranches.png",  "β_Sørensen (ranches)")
    heatmap(df_sim_r,  outdir / "ranches" / "heat_beta_sim_ranches.png",  "β_Simpson turnover (ranches)")
    heatmap(df_sne_r,  outdir / "ranches" / "heat_beta_sne_ranches.png",  "β_Nestedness (ranches)")
    heatmap(df_bray_r, outdir / "ranches" / "heat_bray_curtis_ranches.png","Bray–Curtis (ranches)")

    print("[OK] Wrote beta diversity matrices + heatmaps with preserved ranch order to:", outdir)

if __name__ == "__main__":
    main()