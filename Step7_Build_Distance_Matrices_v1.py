#!/usr/bin/env python3
"""
Step7_Build_Distance_Matrices_v1.py

Builds distance matrices among ranch centroids:

- D_bray.csv: Bray–Curtis on Hellinger (aggregated to ranch centroids)
- D_geo_km.csv: Geodesic distance (km) between ranch centroids
- D_env_scaled.csv: Euclidean distance in z-scored env space (from env_table_sites.csv)
- D_cost_<name>.csv: Least-cost distances on each resistance raster (if scikit-image is available)

USAGE (example with your paths)
-------------------------------
ROOT="/Users/estiner/ranch_meta_2025/bin/baja_ranch_meta_pipeline_v2.0/ericsmetapipe/pipeline bin/MB Master Files/Nogales_METACOMMUNITY"

python Step7_Build_Distance_Matrices_v1.py \
  --hellinger  "$ROOT/USE_Analysis_Pseudo_4_ranches/Step2_Normalized/hellinger_matrix.csv" \
  --samples    "$ROOT/samples.csv" \
  --env_table  "$ROOT/USE_Analysis_Pseudo_4_ranches/Step5_Env/env_table_sites.csv" \
  --resistance "$ROOT/USE_Analysis_Pseudo_4_ranches/Step6_Resistance/combined_global.tif" \
               "$ROOT/USE_Analysis_Pseudo_4_ranches/Step6_Resistance/res_slope_global.tif" \
               "$ROOT/USE_Analysis_Pseudo_4_ranches/Step6_Resistance/res_rough_global.tif" \
  --outdir     "$ROOT/USE_Analysis_Pseudo_4_ranches/Step7_Dist"
"""

import os, argparse, math, warnings
import numpy as np
import pandas as pd
from itertools import combinations
from pyproj import Geod, CRS, Transformer
import rasterio
from rasterio.warp import reproject, Resampling
from rasterio.transform import rowcol
from scipy.spatial.distance import pdist, squareform

def info(s): print(f"[INFO] {s}")
def ok(s):   print(f"[OK]  {s}")
def warn(s): print(f"[WARN] {s}")

def norm_id(x):
    return pd.Series(x, dtype=str).str.strip().str.replace(r"\s+","_", regex=True)

def load_hellinger(path):
    """
    Load a Hellinger matrix in one of two formats.

    1) Your current Step 2 output (taxa-as-rows):
       - First column = Taxon
       - Remaining columns = sample IDs (e.g., La_Venada_1, El_Mesquite_F, ...)
       This is what hellinger_matrix.csv looks like right now.
       We:
         * take all non-Taxon columns as samples
         * force them to numeric
         * transpose so rows = samples, cols = taxa
         * return sample IDs in 'sample_id' and taxa labels

    2) Fallback / alternate: samples-as-rows format
       - Has a 'sample_id' column already, rest = taxa.
       We keep the old behavior for this case.

    Returns
    -------
    hel_ids : DataFrame with one column 'sample_id'
    taxa    : list/Index of taxon names
    X       : numpy array of shape (n_samples x n_taxa)
    """
    df = pd.read_csv(path)
    cols_lower = {c.lower(): c for c in df.columns}

    # ---- Case 1: taxa-as-rows (your Step 2 hellinger)
    if "taxon" in cols_lower:
        tax_col = cols_lower["taxon"]          # usually 'Taxon'
        taxa = df[tax_col].astype(str)

        # all other columns are samples
        mat = df.drop(columns=[tax_col])

        # force numeric, bad values -> NaN -> 0.0
        mat = mat.apply(pd.to_numeric, errors="coerce").fillna(0.0)

        # sample IDs come from column names; normalize them
        sample_ids = norm_id(mat.columns)

        # transpose so rows = samples, cols = taxa
        X = mat.values.T   # shape (n_samples, n_taxa)

        hel_ids = pd.DataFrame({"sample_id": sample_ids})
        return hel_ids, taxa, X

    # ---- Case 2: already samples-as-rows with sample_id
    if "sample_id" in df.columns:
        df["sample_id"] = norm_id(df["sample_id"])
        taxa_cols = [c for c in df.columns if c != "sample_id"]
        X = df[taxa_cols].astype(float).values
        return df[["sample_id"]], taxa_cols, X

    # ---- Fallback: treat first column as sample_id (legacy behavior)
    df = df.rename(columns={df.columns[0]: "sample_id"})
    df["sample_id"] = norm_id(df["sample_id"])
    taxa_cols = [c for c in df.columns if c != "sample_id"]
    X = df[taxa_cols].astype(float).values
    return df[["sample_id"]], taxa_cols, X

def load_samples(path):
    meta = pd.read_csv(path)
    # flexible column names
    cols = {c.lower(): c for c in meta.columns}
    def get(*cands):
        for c in cands:
            if c in cols: return cols[c]
        return None
    c_sid = get("sample_id","sample","id")
    c_ran = get("ranch")
    c_lon = get("lon","longitude","x")
    c_lat = get("lat","latitude","y")
    if not all([c_sid,c_ran,c_lon,c_lat]):
        raise ValueError("samples.csv must have sample_id, ranch, lon, lat")
    meta = meta.rename(columns={c_sid:"sample_id", c_ran:"ranch", c_lon:"lon", c_lat:"lat"})
    meta["sample_id"] = norm_id(meta["sample_id"])
    meta["ranch"]     = norm_id(meta["ranch"])
    meta["lon"] = pd.to_numeric(meta["lon"], errors="coerce")
    meta["lat"] = pd.to_numeric(meta["lat"], errors="coerce")
    meta = meta.dropna(subset=["lon","lat"])
    return meta

def load_env_table(path):
    env = pd.read_csv(path)
    if "ranch" not in env.columns:
        env = env.rename(columns={env.columns[0]: "ranch"})
    env["ranch"] = norm_id(env["ranch"])
    # keep only numeric columns for env distance
    num = env.select_dtypes(include=[np.number]).copy()
    # z-score cols with nonzero std
    for c in list(num.columns):
        sd = float(num[c].std(ddof=0))
        if not np.isfinite(sd) or sd == 0:
            num.drop(columns=[c], inplace=True)
        else:
            num[c] = (num[c] - float(num[c].mean()))/sd
    env_z = pd.concat([env[["ranch"]], num], axis=1)
    return env, env_z

def ensure_dir(d):
    os.makedirs(d, exist_ok=True)

def to_square(df_labels, dist_vec):
    labs = list(df_labels)
    D = squareform(dist_vec)
    return pd.DataFrame(D, index=labs, columns=labs)

def ranch_centroids(meta):
    c = meta.groupby("ranch", as_index=False).agg(
        lon=("lon","mean"), lat=("lat","mean")
    )
    return c

def bray_from_hellinger_by_ranch(hel_ids, X, samples):
    # Attach sample->ranch, average Hellinger rows within ranch
    sid_to_ranch = samples.set_index("sample_id")["ranch"]
    df = pd.DataFrame(X, index=hel_ids["sample_id"])
    df["ranch"] = sid_to_ranch.reindex(df.index).values
    df = df.dropna(subset=["ranch"])
    ranch_means = df.groupby("ranch").mean(numeric_only=True)
    # Bray–Curtis on ranch-level Hellinger means
    bc = pdist(ranch_means.values, metric="braycurtis")
    D = to_square(ranch_means.index, bc)
    return D

def geo_distance_km(centroids):
    geod = Geod(ellps="WGS84")
    ranches = centroids["ranch"].tolist()
    coords = centroids.set_index("ranch")[["lon","lat"]].to_dict("index")
    labs = ranches
    D = np.zeros((len(labs), len(labs)), dtype=float)
    for i,a in enumerate(labs):
        for j,b in enumerate(labs):
            if j<=i: continue
            la,fa = coords[a]["lon"], coords[a]["lat"]
            lb,fb = coords[b]["lon"], coords[b]["lat"]
            _,_,d = geod.inv(la,fa,lb,fb)
            D[i,j] = D[j,i] = d/1000.0
    return pd.DataFrame(D, index=labs, columns=labs)

def env_distance(env_z):
    # env_z: ranch + z-scored numeric columns
    env_z = env_z.dropna(axis=1)  # drop any all-NA vars
    labs = env_z["ranch"]
    X = env_z.drop(columns=["ranch"]).values
    D = pd.DataFrame(squareform(pdist(X, metric="euclidean")),
                     index=labs, columns=labs)
    return D

def read_cost_align(fp, ref_meta):
    with rasterio.open(fp) as src:
        arr = src.read(1, masked=True)
        arr = np.array(np.ma.filled(arr, np.nan), dtype=np.float32)
        # Reproject/resample to ref grid if needed
        if (src.crs != ref_meta["crs"]) or (src.transform != ref_meta["transform"]) or (src.width!=ref_meta["width"]) or (src.height!=ref_meta["height"]):
            dst = np.full((ref_meta["height"], ref_meta["width"]), np.nan, dtype=np.float32)
            reproject(
                source=arr, destination=dst,
                src_transform=src.transform, src_crs=src.crs,
                dst_transform=ref_meta["transform"], dst_crs=ref_meta["crs"],
                src_nodata=np.nan, dst_nodata=np.nan,
                resampling=Resampling.bilinear
            )
            arr = dst
        # ensure strictly positive costs where finite
        finite = np.isfinite(arr)
        if not finite.any():
            raise ValueError("All NaN in cost raster")
        mn = float(np.nanmin(arr))
        if not np.isfinite(mn):
            mn = 0.0
        arr = np.where(finite, arr, np.nan)
        # normalize 0..1 then add floor
        mx = float(np.nanmax(arr))
        if np.isfinite(mx) and mx>mn:
            arr = (arr - mn)/(mx-mn)
        arr = np.where(np.isfinite(arr), arr + 0.05, np.nan).astype(np.float32)
        return arr

def build_ref_grid_from_one(res_fp):
    with rasterio.open(res_fp) as src:
        meta = src.meta.copy()
        meta.update(count=1, dtype="float32", nodata=np.nan)
    return meta

def cost_distance_matrix(cost_fp, centroids):
    # Reference grid taken from this cost raster itself
    ref_meta = build_ref_grid_from_one(cost_fp)
    cost = read_cost_align(cost_fp, ref_meta)
    # Map centroids (lon/lat) -> row/col in cost grid (assumes cost grid CRS is projected)
    crs_cost = ref_meta["crs"]
    transformer = Transformer.from_crs(CRS.from_epsg(4326), crs_cost, always_xy=True)
    xs, ys = transformer.transform(centroids["lon"].values, centroids["lat"].values)
    # Convert to indices
    rr = []; cc=[]
    for x,y in zip(xs,ys):
        r,c = rowcol(ref_meta["transform"], x, y, op=float)
        r = int(round(r)); c = int(round(c))
        rr.append(r); cc.append(c)
    # Check bounds
    H,W = cost.shape
    inb = [(0<=r<H and 0<=c<W and np.isfinite(cost[r,c])) for r,c in zip(rr,cc)]
    if not all(inb):
        warn("Some centroids fall outside/masked cost grid; will skip those pairs.")
    # Try least-cost via scikit-image
    try:
        from skimage.graph import route_through_array
        labs = centroids["ranch"].tolist()
        D = np.full((len(labs), len(labs)), np.nan, dtype=float)
        px = abs(ref_meta["transform"].a)
        for i in range(len(labs)):
            for j in range(i+1, len(labs)):
                if not (inb[i] and inb[j]): 
                    continue
                start = (rr[i], cc[i]); end=(rr[j], cc[j])
                _, cost_sum = route_through_array(cost, start, end,
                                                  fully_connected=True, geometric=True)
                # cost_sum is unitless “cost”; multiply by pixel size to express "cost-meters"
                D[i,j] = D[j,i] = float(cost_sum) * px
        return pd.DataFrame(D, index=labs, columns=labs)
    except Exception as e:
        warn(f"Least-cost failed for {os.path.basename(cost_fp)} ({e}); skipping.")
        return None

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--hellinger", required=True)
    ap.add_argument("--samples",   required=True)
    ap.add_argument("--env_table", required=True)
    ap.add_argument("--resistance", nargs="*", default=[], help="one or more .tif resistance rasters")
    ap.add_argument("--outdir",    required=True)
    args = ap.parse_args()

    ensure_dir(args.outdir)

    # Load data
    hel_ids, taxa, X = load_hellinger(args.hellinger)
    samples = load_samples(args.samples)
    env, env_z = load_env_table(args.env_table)

    # Centroids
    cents = ranch_centroids(samples)  # ranch, lon, lat
    ok(f"Ranches: {', '.join(cents['ranch'])}")

    # D_bray (Hellinger -> Bray at ranch)
    D_bray = bray_from_hellinger_by_ranch(hel_ids, X, samples)
    D_bray.to_csv(os.path.join(args.outdir, "D_bray.csv"))
    ok("Wrote D_bray.csv")

    # D_geo (km)
    cents_geo = cents.copy()
    D_geo = geo_distance_km(cents_geo)
    D_geo.to_csv(os.path.join(args.outdir, "D_geo_km.csv"))
    ok("Wrote D_geo_km.csv")

    # D_env (z-scored env space)
    # Keep rows that match ranch centroids order
    env_z_use = env_z[env_z["ranch"].isin(cents["ranch"])].copy()
    # Reorder
    env_z_use = env_z_use.set_index("ranch").loc[cents["ranch"]].reset_index()
    D_env = env_distance(env_z_use)
    D_env.to_csv(os.path.join(args.outdir, "D_env_scaled.csv"))
    ok("Wrote D_env_scaled.csv")

    # D_cost_* for each resistance raster
    for fp in args.resistance:
        name = os.path.splitext(os.path.basename(fp))[0]
        info(f"Cost distance: {name}")
        Dc = cost_distance_matrix(fp, cents)
        if Dc is not None:
            out = os.path.join(args.outdir, f"D_cost_{name}.csv")
            Dc.to_csv(out)
            ok(f"Wrote {os.path.basename(out)}")

    # Summary
    with open(os.path.join(args.outdir, "INDEX.txt"), "w") as f:
        f.write("Step7 distance outputs\n")
        f.write("- D_bray.csv          (Bray–Curtis on ranch-mean Hellinger)\n")
        f.write("- D_geo_km.csv        (geodesic km)\n")
        f.write("- D_env_scaled.csv    (Euclidean in z-scored env space)\n")
        if len(args.resistance):
            f.write("- D_cost_*.csv        (least-cost distances per resistance raster)\n")

    ok(f"All outputs written to: {os.path.abspath(args.outdir)}")

if __name__ == "__main__":
    # quieter numpy warnings on degenerate slices
    np.seterr(all="ignore")
    main()