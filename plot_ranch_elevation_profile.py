#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
========================================
STUDENT HEADER: What this script does
========================================
Given:
  1) ranch_segment_results.csv  -> per-segment summary (from your segment script)
  2) ranch_vertices.csv         -> all vertices with columns [lat, lon, elev_ft, segment]

This script:
  • Reconstructs a single continuous transect following YOUR ranch order.
  • Computes cumulative distance along the transect (miles).
  • Plots elevation (feet) vs distance with vertical markers and labels at ranch boundaries.
  • Saves a high-resolution PNG and PDF for publication.

Usage (run from the folder with the CSVs):
  python3 plot_ranch_elevation_profile.py \
      --segments ranch_segment_results.csv \
      --vertices ranch_vertices.csv \
      --out-prefix ranch_elevation_profile

Output files:
  ranch_elevation_profile.png
  ranch_elevation_profile.pdf

Requirements:
  pandas, numpy, matplotlib
"""

# --------- Imports (you don't need to change these) ---------
import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from math import radians, sin, cos, sqrt, atan2
from pathlib import Path

# --------- Simple distance function (meters) ---------
def haversine_m(lat1, lon1, lat2, lon2):
    """Great-circle distance between two points (lat/lon in degrees)."""
    R = 6371000.0
    p1, p2 = radians(lat1), radians(lat2)
    dphi  = radians(lat2 - lat1)
    dlmb  = radians(lon2 - lon1)
    a = sin(dphi/2)**2 + cos(p1)*cos(p2)*sin(dlmb/2)**2
    return 2 * R * atan2(sqrt(a), sqrt(1-a))

# --------- CLI args ---------
ap = argparse.ArgumentParser(description="Plot elevation profile with ranch markers.")
ap.add_argument("--segments", default="ranch_segment_results.csv", help="Per-segment CSV.")
ap.add_argument("--vertices", default="ranch_vertices.csv", help="Per-vertex CSV with elev_ft and segment label.")
ap.add_argument("--out-prefix", default="ranch_elevation_profile", help="Prefix for PNG/PDF outputs.")
args = ap.parse_args()

# --------- Read inputs ---------
segs  = pd.read_csv(args.segments)
verts = pd.read_csv(args.vertices)

# Basic sanity on columns
need = {"lat","lon","elev_ft","segment"}
missing = need - set(verts.columns)
if missing:
    raise SystemExit(f"[ERROR] {args.vertices} missing columns: {missing}. Found: {list(verts.columns)}")

if segs.empty:
    raise SystemExit(f"[ERROR] {args.segments} is empty.")

# --------- Reconstruct your ranch order from the segment table ----------
# Example rows:
#   from: La Venada  to: El Mesquite
#   from: El Mesquite  to: Casa Viejas
#   from: Casa Viejas  to: El Tule
order = [segs.loc[0, "from"]]
for _, row in segs.iterrows():
    order.append(row["to"])

# De-duplicate in sequence -> [La Venada, El Mesquite, Casa Viejas, El Tule]
ranch_order = []
for r in order:
    if not ranch_order or ranch_order[-1] != r:
        ranch_order.append(r)

# Segment labels must match verts['segment'] exactly (e.g., "La Venada -> El Mesquite")
segment_labels = [f"{ranch_order[i]} -> {ranch_order[i+1]}" for i in range(len(ranch_order)-1)]

# --------- Build one continuous profile ----------
profile_pieces = []
x_offset_m = 0.0
boundary_marks = []  # list of tuples: (x_m, ranch_name, elev_ft_at_boundary)

for i, lab in enumerate(segment_labels):
    seg_df = verts[verts["segment"] == lab].copy()
    if seg_df.empty:
        # Sometimes non-breaking spaces sneak in; normalize as a fallback
        seg_df = verts[verts["segment"].str.replace("\u00a0"," ", regex=False) == lab].copy()
    if seg_df.empty:
        raise SystemExit(f"[ERROR] No vertices found for segment label: '{lab}'. "
                         f"Check spellings in both CSVs.")

    # Compute cumulative distance along THIS segment (meters)
    pts = seg_df[["lat","lon"]].to_numpy()
    d = np.zeros(len(seg_df), dtype=float)
    if len(seg_df) > 1:
        step = [haversine_m(pts[j,0], pts[j,1], pts[j+1,0], pts[j+1,1]) for j in range(len(pts)-1)]
        d[1:] = np.cumsum(step)

    # Shift by the current total offset so segments stitch correctly
    seg_df["x_m"] = d + x_offset_m
    profile_pieces.append(seg_df[["x_m","elev_ft"]])

    # Record boundary ticks:
    # - At start of the very first segment -> first ranch
    if i == 0:
        start_elev = float(seg_df["elev_ft"].iloc[0]) if seg_df["elev_ft"].notna().any() else float("nan")
        boundary_marks.append((x_offset_m, ranch_order[0], start_elev))
    # - At end of this segment -> next ranch
    x_offset_m = float(seg_df["x_m"].iloc[-1])
    end_elev = float(seg_df["elev_ft"].iloc[-1]) if seg_df["elev_ft"].notna().any() else float("nan")
    boundary_marks.append((x_offset_m, ranch_order[i+1], end_elev))

# Concatenate to one DataFrame
profile = pd.concat(profile_pieces, ignore_index=True)

# --------- Plot (one chart; neutral style) ----------
plt.figure(figsize=(11, 4), dpi=200)
plt.plot(profile["x_m"] / 1609.34, profile["elev_ft"])
plt.xlabel("Distance (miles)")
plt.ylabel("Elevation (ft)")
plt.title("Elevation Profile Across Baja Ranch Transect")

# Vertical dashed lines + label at each ranch
for x_m, name, elev in boundary_marks:
    x_mi = x_m / 1609.34
    plt.axvline(x=x_mi, linestyle="--", linewidth=1)
    # place label slightly above boundary elevation if available; otherwise near bottom
    y = elev + 30 if np.isfinite(elev) else np.nan
    if np.isfinite(y):
        plt.text(x_mi, y, name, rotation=90, va="bottom", ha="center")
    else:
        ymin, ymax = plt.ylim()
        plt.text(x_mi, ymin + 0.05*(ymax - ymin), name, rotation=90, va="bottom", ha="center")

plt.tight_layout()

# --------- Save high-res outputs ----------
out_png = Path(f"{args.out_prefix}.png")
out_pdf = Path(f"{args.out_prefix}.pdf")
plt.savefig(out_png, bbox_inches="tight")
plt.savefig(out_pdf, bbox_inches="tight")
print(f"[OK] Wrote: {out_png.resolve()}")
print(f"[OK] Wrote: {out_pdf.resolve()}")

# --------- Optional: print basic stats for the caption ----------
total_miles = (profile['x_m'].iloc[-1] - profile['x_m'].iloc[0]) / 1609.34
gain_pos = np.maximum(np.diff(profile['elev_ft'].to_numpy(dtype=float)), 0).sum()
loss_neg = np.minimum(np.diff(profile['elev_ft'].to_numpy(dtype=float)), 0).sum()
print(f"[INFO] Total distance ≈ {total_miles:.1f} mi | Elevation gain ≈ {gain_pos:.0f} ft | loss ≈ {abs(loss_neg):.0f} ft")