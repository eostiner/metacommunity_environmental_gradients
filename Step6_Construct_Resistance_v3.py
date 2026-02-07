#!/usr/bin/env python3
"""
Step6_Construct_Resistance_v3.py

Build both:
  (A) combined_global.tif                     — physical baseline (slope + roughness + climate, globally scaled 0..1)
  (B) combined_sites_climate_only.tif         — site-referenced climate mismatch blend (climate-only, 0..1)

Also writes individual layers:
  - res_slope_global.tif, res_rough_global.tif
  - res_climate_<KEY>_global.tif
  - res_climate_<KEY>_sites.tif   (if --sites given and sampling succeeds)

Quicklook .png previews are optional (on by default) and lightweight.

USAGE (example with your paths)
-------------------------------
ROOT="/Users/estiner/ranch_meta_2025/bin/baja_ranch_meta_pipeline_v2.0/ericsmetapipe/pipeline bin/MB Master Files/Nogales_METACOMMUNITY"

python Step6_Construct_Resistance_v3.py \
  --dem "$ROOT/GEO_TIFFS/DEM_SRTM30_RANCH_AOI.tif" \
  --climate \
    "$ROOT/GEO_TIFFS/ERA5_TmaxK_RANCH_AOI_1979_2024.tif" \
    "$ROOT/GEO_TIFFS/ERA5_TmeanK_RANCH_AOI_1979_2024.tif" \
    "$ROOT/GEO_TIFFS/MODIS_EVI_Amplitude_RANCH_AOI_2001_2024.tif" \
    "$ROOT/GEO_TIFFS/MODIS_EVI_Area_RANCH_AOI_2001_2020.tif" \
    "$ROOT/GEO_TIFFS/MODIS_Greenup_RANCH_AOI_2001_2020.tif" \
  --sites "$ROOT/samples.csv" \
  --weights_global "SlopeDeg_SRTM30_RANCH_AOI=0.4,RoughStdM_SRTM30_RANCH_AOI=0.3,ERA5_TmeanK_RANCH_AOI_1979_2024=0.2,MODIS_EVI_Amplitude_RANCH_AOI_2001_2024=0.1" \
  --weights_sites  "ERA5_TmeanK_RANCH_AOI_1979_2024=0.7,MODIS_EVI_Amplitude_RANCH_AOI_2001_2024=0.3" \
  --outdir "$ROOT/USE_Analysis_Pseudo_4_ranches/Step6_Resistance"

Notes:
- weights_global can include slope/roughness + any climate keys.
- weights_sites should name only climate keys that have *_sites mismatch layers.
- If --weights_sites is omitted, the script uses equal weights across all available *_sites layers.
"""

import os, argparse, math, warnings
import numpy as np
import pandas as pd
import rasterio
from rasterio.enums import Resampling
from rasterio.transform import Affine
from rasterio.warp import reproject
from pyproj import CRS, Transformer
from shapely.geometry import Point
from shapely.ops import transform as shp_transform
from scipy.ndimage import generic_filter

def info(s): print(f"[INFO] {s}")
def ok(s):   print(f"[OK]  {s}")
def warn(s): print(f"[WARN] {s}")
def err(s):  print(f"[ERR]  {s}")

def parse_weights(s):
    if not s: return {}
    out = {}
    for chunk in s.split(","):
        k, v = chunk.split("=")
        out[k.strip()] = float(v)
    # normalize to sum=1
    total = sum(out.values())
    if total > 0:
        for k in out: out[k] /= total
    return out

def read_to_match(fp, ref_meta, dst_crs=None, resampling=Resampling.bilinear):
    with rasterio.open(fp) as src:
        src_arr = src.read(1, masked=True)
        src_arr = np.array(np.ma.filled(src_arr, np.nan), dtype=np.float32)
        dst = np.full((ref_meta['height'], ref_meta['width']), np.nan, dtype=np.float32)
        reproject(
            source=src_arr,
            destination=dst,
            src_transform=src.transform,
            src_crs=src.crs,
            dst_transform=ref_meta['transform'],
            dst_crs=dst_crs or ref_meta['crs'],
            src_nodata=np.nan, dst_nodata=np.nan,
            resampling=resampling
        )
    return dst

def write_raster(fp, arr, meta):
    m = meta.copy()
    m.update(count=1, dtype="float32", nodata=np.nan, compress="deflate", predictor=2, zlevel=4)
    with rasterio.open(fp, "w", **m) as dst:
        dst.write(arr.astype(np.float32), 1)

def quicklook_png(fp_png, arr):
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        v = np.nanpercentile(arr, [2, 98])
        plt.figure(figsize=(5,4))
        plt.imshow(arr, vmin=v[0], vmax=v[1])
        plt.axis("off")
        plt.tight_layout()
        plt.savefig(fp_png, dpi=150)
        plt.close()
    except Exception as e:
        warn(f"quicklook failed: {e}")

def to_utm_meta(dem_fp):
    with rasterio.open(dem_fp) as src:
        if CRS(src.crs).to_epsg() == 4326:
            # build a local UTM from center
            h, w = src.height, src.width
            cx = src.bounds.left + (src.bounds.right - src.bounds.left)/2
            cy = src.bounds.bottom + (src.bounds.top - src.bounds.bottom)/2
            zone = int((cx + 180)//6) + 1
            south = (cy < 0)
            utm = CRS.from_epsg((32700 if south else 32600) + zone)
            # reproject DEM to UTM at ~native meter resolution
            res_x = (src.bounds.right - src.bounds.left) / w
            res_y = (src.bounds.top - src.bounds.bottom) / h
            # pick a reasonable ~30m target if degrees; otherwise keep scale
            target_res = 30.0
            # compute target shape from bounds in UTM
            # first reproject corners to UTM to get span
            transformer = Transformer.from_crs(src.crs, utm, always_xy=True)
            (x0,y0) = transformer.transform(src.bounds.left,  src.bounds.bottom)
            (x1,y1) = transformer.transform(src.bounds.right, src.bounds.top)
            width  = max(1, int(round((x1-x0)/target_res)))
            height = max(1, int(round((y1-y0)/target_res)))
            transform = Affine(target_res, 0, min(x0,x1), 0, -target_res, max(y0,y1))
            meta = dict(driver="GTiff", width=width, height=height, count=1, crs=utm, transform=transform, dtype="float32")
            # reproject dem band
            src_arr = np.array(np.ma.filled(src.read(1, masked=True), np.nan), dtype=np.float32)
            dst = np.full((height, width), np.nan, dtype=np.float32)
            reproject(
                source=src_arr, destination=dst,
                src_transform=src.transform, src_crs=src.crs,
                dst_transform=transform, dst_crs=utm,
                src_nodata=np.nan, dst_nodata=np.nan,
                resampling=Resampling.bilinear
            )
            ok(f"DEM reprojected to {utm.to_string()} with shape ({height}, {width})")
            return dst, meta
        else:
            # already projected
            arr = np.array(np.ma.filled(src.read(1, masked=True), np.nan), dtype=np.float32)
            meta = src.meta.copy()
            meta.update(dtype="float32")
            ok(f"DEM already projected ({src.crs}); shape {arr.shape}")
            return arr, meta

def slope_from_dem_m(arr, res):
    # simple finite-diff slope (deg); res = pixel size (m)
    gy, gx = np.gradient(arr, res, edge_order=1)
    slope_rad = np.arctan(np.sqrt(gx*gx + gy*gy))
    return np.degrees(slope_rad)

def roughness_std(arr, size=3):
    def _std9(x): return np.nanstd(x)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)
        out = generic_filter(arr, _std9, size=size, mode="nearest")
    return out

def scale01(a):
    mn = np.nanmin(a); mx = np.nanmax(a)
    if not np.isfinite(mn) or not np.isfinite(mx) or mx <= mn:
        return np.zeros_like(a, dtype=np.float32)
    out = (a - mn) / (mx - mn)
    return out.astype(np.float32)

def z_mismatch01(a, mean_site, std_site, clip=3.0):
    # absolute z-score -> clip to ±clip -> 0..1 by /clip
    z = np.abs((a - mean_site) / (std_site if std_site > 0 else np.nan))
    z = np.where(np.isfinite(z), np.minimum(z, clip), np.nan)
    return (z / clip).astype(np.float32)

def read_sites_csv(fp, crs_src, crs_dst):
    df = pd.read_csv(fp)
    # flexible lon/lat column names
    cols = {c.lower(): c for c in df.columns}
    def pick(*names):
        for n in names:
            if n in cols: return cols[n]
        return None
    c_lon = pick("lon","longitude","x"); c_lat = pick("lat","latitude","y")
    c_ranch = pick("ranch")
    if not (c_lon and c_lat):
        raise ValueError("sites CSV must have lon/lat columns")
    lon = pd.to_numeric(df[c_lon], errors="coerce").values
    lat = pd.to_numeric(df[c_lat], errors="coerce").values
    mask = np.isfinite(lon) & np.isfinite(lat)
    lon, lat = lon[mask], lat[mask]
    transformer = Transformer.from_crs(crs_src, crs_dst, always_xy=True)
    xs, ys = transformer.transform(lon, lat)
    return np.array(xs), np.array(ys)

def sample_at_points(raster_arr, transform, xs, ys):
    # map x,y to row,col (col=x, row=y)
    cols = ((xs - transform.c) / transform.a).astype(int)
    rows = ((ys - transform.f) / transform.e).astype(int)
    mask = (rows>=0)&(rows<raster_arr.shape[0])&(cols>=0)&(cols<raster_arr.shape[1])
    vals = raster_arr[rows[mask], cols[mask]]
    vals = vals[np.isfinite(vals)]
    return vals

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dem", required=True)
    ap.add_argument("--climate", nargs="*", default=[])
    ap.add_argument("--sites", type=str, default=None)
    ap.add_argument("--weights_global", type=str, default="")  # keys: slope/roughness/climate
    ap.add_argument("--weights_sites",  type=str, default="")  # keys: climate ONLY (for *_sites)
    ap.add_argument("--quicklook", action="store_true", default=True)
    ap.add_argument("--outdir", required=True)
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    # --- DEM to UTM grid ---
    info("Preparing DEM grid…")
    dem_arr_utm, dem_meta = to_utm_meta(args.dem)
    res_m = abs(dem_meta['transform'].a)

    # --- Slope & Roughness (global) ---
    info("Computing slope & roughness…")
    slope_deg = slope_from_dem_m(dem_arr_utm, res=res_m)
    rough_std = roughness_std(dem_arr_utm, size=3)
    res_slope = scale01(slope_deg)           # higher = more resistance
    res_rough = scale01(rough_std)

    fp_slope = os.path.join(args.outdir, "res_slope_global.tif")
    fp_rough = os.path.join(args.outdir, "res_rough_global.tif")
    write_raster(fp_slope, res_slope, dem_meta)
    write_raster(fp_rough, res_rough, dem_meta)
    ok("Wrote slope/roughness (global).")
    if args.quicklook:
        quicklook_png(fp_slope.replace(".tif",".png"), res_slope)
        quicklook_png(fp_rough.replace(".tif",".png"), res_rough)

    # --- Climate rasters (global + sites mismatch) ---
    xs = ys = None
    if args.sites:
        try:
            xs, ys = read_sites_csv(args.sites, CRS.from_epsg(4326), dem_meta['crs'])
            ok(f"Loaded {len(xs)} sites.")
        except Exception as e:
            warn(f"Sites failed; proceeding without site mismatch ({e})")

    res_clim_global = {}  # key -> array
    res_clim_sites  = {}  # key -> array

    for fp in args.climate:
        key = os.path.splitext(os.path.basename(fp))[0]
        info(f"Climate: {fp} -> key={key}")
        arr = read_to_match(fp, dem_meta)
        # GLOBAL scaling 0..1 (more = more resistance)
        res_g = scale01(arr)
        res_clim_global[key] = res_g
        out_g = os.path.join(args.outdir, f"res_climate_{key}_global.tif")
        write_raster(out_g, res_g, dem_meta)
        if args.quicklook:
            quicklook_png(out_g.replace(".tif",".png"), res_g)

        # SITE mismatch 0..1 via absolute z-score against site mean/std
        if xs is not None and ys is not None and len(xs) >= 2:
            vals = sample_at_points(arr, dem_meta['transform'], xs, ys)
            if np.size(vals) >= 2 and np.isfinite(vals).any():
                m = float(np.nanmean(vals)); s = float(np.nanstd(vals))
                if s > 0:
                    res_s = z_mismatch01(arr, m, s, clip=3.0)
                    res_clim_sites[key] = res_s
                    out_s = os.path.join(args.outdir, f"res_climate_{key}_sites.tif")
                    write_raster(out_s, res_s, dem_meta)
                    if args.quicklook:
                        quicklook_png(out_s.replace(".tif",".png"), res_s)
                else:
                    warn(f"{key}: site std = 0; skipping sites mismatch.")
            else:
                warn(f"{key}: no finite site samples; skipping sites mismatch.")

    # --- Combined GLOBAL (slope + roughness + climate global) ---
    wg = parse_weights(args.weights_global)
    # If user provided nothing, default blend: slope 0.5, rough 0.3, average of all climate globals 0.2
    if not wg:
        wg = {}
        if res_slope is not None: wg["SlopeDeg_SRTM30_RANCH_AOI"] = 0.5
        if res_rough is not None: wg["RoughStdM_SRTM30_RANCH_AOI"] = 0.3
        if len(res_clim_global):
            # split the remaining 0.2 equally among climate layers
            rem = 0.2 / max(1, len(res_clim_global))
            for k in res_clim_global.keys():
                wg[k] = rem
        # renormalize
        s = sum(wg.values())
        if s > 0:
            for k in wg: wg[k] /= s

    info(f"Combining GLOBAL with weights: {wg}")
    combo = np.zeros_like(dem_arr_utm, dtype=np.float32)
    for k, w in wg.items():
        if k == "SlopeDeg_SRTM30_RANCH_AOI":
            combo += w * res_slope
        elif k == "RoughStdM_SRTM30_RANCH_AOI":
            combo += w * res_rough
        elif k in res_clim_global:
            combo += w * res_clim_global[k]
        else:
            warn(f"Weight key '{k}' not found; skipping.")
    out_combo_g = os.path.join(args.outdir, "combined_global.tif")
    write_raster(out_combo_g, combo, dem_meta)
    if args.quicklook:
        quicklook_png(out_combo_g.replace(".tif",".png"), combo)
    ok("Wrote combined_global.tif")

    # --- Combined SITES (climate-only mismatch), always if any *_sites exist ---
    if len(res_clim_sites):
        ws = parse_weights(args.weights_sites)
        if not ws:
            # equal weights across all *_sites layers
            eq = 1.0 / len(res_clim_sites)
            ws = {k: eq for k in res_clim_sites.keys()}
        info(f"Combining SITES (climate-only) with weights: {ws}")
        combo_s = np.zeros_like(dem_arr_utm, dtype=np.float32)
        for k, w in ws.items():
            if k in res_clim_sites:
                combo_s += w * res_clim_sites[k]
            else:
                warn(f"Sites key '{k}' not found among climate mismatch layers; skipping.")
        out_combo_s = os.path.join(args.outdir, "combined_sites_climate_only.tif")
        write_raster(out_combo_s, combo_s, dem_meta)
        if args.quicklook:
            quicklook_png(out_combo_s.replace(".tif",".png"), combo_s)
        ok("Wrote combined_sites_climate_only.tif")
    else:
        warn("No climate mismatch (sites) layers produced; combined_sites_climate_only not written.")

    # --- SUMMARY ---
    with open(os.path.join(args.outdir, "SUMMARY.txt"), "w") as f:
        f.write("Step6 v3 outputs\n")
        f.write(f"- Grid CRS: {dem_meta['crs']}\n")
        f.write(f"- Global: slope + roughness + {len(res_clim_global)} climate layers\n")
        f.write(f"- Sites mismatch layers: {len(res_clim_sites)}\n")
        f.write(f"- combined_global.tif written with weights={wg}\n")
        if len(res_clim_sites):
            f.write(f"- combined_sites_climate_only.tif written\n")

    ok(f"All outputs in: {args.outdir}")

if __name__ == "__main__":
    main()