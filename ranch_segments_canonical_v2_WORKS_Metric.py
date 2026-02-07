#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ranch_segments_canonical_v2.py

Student notes & goals
---------------------
1) Input is a single KML track exported from your GPS app.
   We handle BOTH KML styles:
     - <gx:Track><gx:coord>lon lat alt</gx:coord>...</gx:Track>
     - <LineString><coordinates>lon,lat,alt ...</coordinates></LineString>
2) We also read <Placemark><Point> for named ranches and snap each name to the
   nearest vertex along the track.
3) We automatically pick the track orientation that makes your ranch names
   appear in the order you requested. If needed, we reverse the track.
4) We then TRIM the data to just the span from the first ranch to the last
   ranch in that canonical direction, so *only one direction* remains.
5) We compute:
     - per-vertex cumulative distance and elevation
     - per-leg distance, elevation gain, elevation loss
6) We save:
     - vertices CSV (`--out-vertices`)
     - segment summary CSV (`--out-summary`)
     - labeled elevation profile PNG (`--plot`)

CLI
---
python3 ranch_segments_canonical_v2.py \
  --kml onRanch_route_Xmaps-09_11_25-140607.kml \
  --names "La Venada,El Mesquite,Casa Viejas,El Tule" \
  --out-vertices ranch_vertices.csv \
  --out-summary  ranch_segment_results.csv \
  --plot         ranch_elevation_profile.png \
  --unit miles \
  --debug

Dependencies: only stdlib + numpy + pandas + matplotlib
(Parses KML with xml.etree; no extra KML packages needed.)
"""

import argparse
import math
import sys
import xml.etree.ElementTree as ET
from typing import List, Tuple, Dict, Optional

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# ----------------------------- helpers ---------------------------------

def haversine_m(lat1, lon1, lat2, lon2) -> float:
    """Great-circle distance in meters."""
    R = 6371000.0
    p1, p2 = math.radians(lat1), math.radians(lat2)
    dphi = math.radians(lat2 - lat1)
    dlmb = math.radians(lon2 - lon1)
    a = math.sin(dphi/2)**2 + math.cos(p1)*math.cos(p2)*math.sin(dlmb/2)**2
    return 2*R*math.asin(math.sqrt(a))

def meters_to_unit(m, unit: str) -> float:
    unit = unit.lower()
    if unit in ("mile", "miles"):
        return m / 1609.344
    if unit in ("meter", "meters"):
        return m
    if unit in ("foot", "feet"):
        return m * 3.28084
    raise ValueError(f"Unknown unit: {unit}")

def parse_coords_text(coords_text: str) -> List[Tuple[float,float,float]]:
    """Parse KML coordinates string 'lon,lat,alt lon,lat,alt ...' into [(lat,lon,alt_m)]."""
    pts = []
    for tok in coords_text.strip().replace("\n", " ").split():
        parts = tok.split(",")
        if len(parts) < 2:
            continue
        lon = float(parts[0])
        lat = float(parts[1])
        alt = float(parts[2]) if len(parts) >= 3 and parts[2] != "" else 0.0
        pts.append((lat, lon, alt))
    return pts

# ------------------------- KML parsing ---------------------------------

def parse_kml_track(kml_path: str) -> Tuple[List[Tuple[float,float,float]], Dict[str, Tuple[float,float,float]]]:
    """Return (track_vertices, placemarks) from KML.
       track_vertices: list of (lat, lon, alt_m)
       placemarks: {name: (lat, lon, alt_m)}
    """
    ns = {
        "kml": "http://www.opengis.net/kml/2.2",
        "gx":  "http://www.google.com/kml/ext/2.2",
    }
    tree = ET.parse(kml_path)
    root = tree.getroot()

    # 1) Try gx:Track
    coords = []
    for tr in root.findall(".//gx:Track", ns):
        for c in tr.findall(".//gx:coord", ns):
            vals = c.text.strip().split()
            if len(vals) >= 2:
                lon = float(vals[0]); lat = float(vals[1])
                alt = float(vals[2]) if len(vals) >= 3 else 0.0
                coords.append((lat, lon, alt))

    # 2) Fallback: LineString coordinates
    if not coords:
        for ls in root.findall(".//kml:LineString", ns):
            cnode = ls.find("kml:coordinates", ns)
            if cnode is not None and cnode.text:
                coords.extend(parse_coords_text(cnode.text))

    if not coords:
        raise RuntimeError("No Track/LineString coordinates found in the KML.")

    # Placemark Points (ranch names)
    placemarks: Dict[str, Tuple[float,float,float]] = {}
    for pm in root.findall(".//kml:Placemark", ns):
        name_node = pm.find("kml:name", ns)
        pnode = pm.find(".//kml:Point/kml:coordinates", ns)
        if name_node is not None and pnode is not None and pnode.text:
            name = name_node.text.strip()
            pts = parse_coords_text(pnode.text)
            if pts:
                placemarks[name] = pts[0]

    return coords, placemarks

# ---------------------- snapping & orientation --------------------------

def snap_names_to_track(track: List[Tuple[float,float,float]],
                        wanted_names: List[str],
                        placemarks: Dict[str, Tuple[float,float,float]]) -> Dict[str, int]:
    """Map each wanted name to the nearest vertex index along track."""
    # If the KML didn’t include named Point placemarks, we’ll try to infer
    # from the names list by ignoring snapping (we’ll handle in caller).
    idxs = {}
    if not placemarks:
        return idxs

    latlon = np.array([[lat, lon] for (lat,lon,alt) in track])
    for nm in wanted_names:
        if nm not in placemarks:
            continue
        plat, plon, _ = placemarks[nm]
        # brute-force nearest
        dists = np.hypot(latlon[:,0] - plat, latlon[:,1] - plon)
        idxs[nm] = int(np.argmin(dists))
    return idxs

def choose_canonical_orientation(n: int, order: List[str], name2idx: Dict[str,int]) -> Tuple[bool, Dict[str,int]]:
    """Decide whether to keep as-is or reverse the track so that indices increase
       in the provided order. Returns (reversed?, adjusted_mapping).
    """
    # Gather indices in the order list (skip names we failed to snap)
    seq = [name2idx[nm] for nm in order if nm in name2idx]
    if len(seq) < 2:
        return False, name2idx  # nothing to decide

    increasing = all(b > a for a,b in zip(seq, seq[1:]))
    if increasing:
        return False, name2idx

    # Try the reversed view: idx' = (n-1) - idx
    rev_map = {nm: (n-1 - idx) for nm, idx in name2idx.items()}
    rev_seq = [rev_map[nm] for nm in order if nm in rev_map]
    rev_increasing = all(b > a for a,b in zip(rev_seq, rev_seq[1:]))

    if rev_increasing:
        return True, rev_map

    # If neither strictly increasing, pick the orientation with the larger positive sum of forward deltas
    def forward_score(s): return sum(max(0, b-a) for a,b in zip(s,s[1:]))
    return (forward_score(rev_seq) > forward_score(seq)), (rev_map if forward_score(rev_seq) > forward_score(seq) else name2idx)

# --------------------------- main logic ---------------------------------

def build_dataframe(track: List[Tuple[float,float,float]]) -> pd.DataFrame:
    lat = np.array([p[0] for p in track])
    lon = np.array([p[1] for p in track])
    alt_m = np.array([p[2] for p in track], dtype=float)

    # distances
    seg_m = np.zeros_like(lat, dtype=float)
    for i in range(1, len(lat)):
        seg_m[i] = haversine_m(lat[i-1], lon[i-1], lat[i], lon[i])
    cum_m = np.cumsum(seg_m)

    df = pd.DataFrame({
        "lat": lat,
        "lon": lon,
        "elev_m": alt_m,
        "cum_dist_m": cum_m
    })
    df["elev_ft"] = df["elev_m"] * 3.28084
    df["cum_dist_mi"] = df["cum_dist_m"] / 1609.344
    return df

def slice_one_direction(df: pd.DataFrame, first_idx: int, last_idx: int) -> pd.DataFrame:
    """Return trimmed dataframe from first_idx..last_idx inclusive (ensure first<last)."""
    if last_idx <= first_idx:
        raise RuntimeError(f"Indices not increasing after orientation fix (start={first_idx}, end={last_idx}).")
    sub = df.iloc[first_idx:last_idx+1].copy()
    # re-zero distance from the new start
    sub["cum_dist_m"]  = sub["cum_dist_m"]  - float(sub["cum_dist_m"].iloc[0])
    sub["cum_dist_mi"] = sub["cum_dist_mi"] - float(sub["cum_dist_mi"].iloc[0])
    return sub

def leg_stats(df: pd.DataFrame, i0: int, i1: int, unit: str) -> Dict[str, float]:
    sub = df.iloc[i0:i1+1]
    d_m = float(sub["cum_dist_m"].iloc[-1] - sub["cum_dist_m"].iloc[0])
    elev = sub["elev_ft"].to_numpy()
    de = np.diff(elev)
    gain = float(np.clip(de, 0, None).sum())
    loss = float(np.clip(-de, 0, None).sum())
    return {
        "distance": meters_to_unit(d_m, unit),
        "elev_gain_ft": gain,
        "elev_loss_ft": loss
    }

def plot_profile(df, ticks_at, out_png, unit):
    # Use meters → km for X; meters for Y
    x = (df["cum_dist_m"].to_numpy() / 1000.0)   # km
    y = df["elev_m"].to_numpy()                  # m

    import matplotlib.pyplot as plt
    plt.figure(figsize=(14, 3.5), dpi=160)
    plt.plot(x, y, linewidth=1.5)

    # ticks_at currently in miles; convert to km for plotting
    for name, xmile in ticks_at:
        xkm = xmile * 1.609344
        plt.axvline(x=xkm, linestyle="--", linewidth=1, color="tab:blue", alpha=0.6)
        plt.text(xkm, y.max()*0.95, name, rotation=90, ha="right", va="top", fontsize=10)

    plt.title("Elevation Profile Across Baja Ranch Transect", fontsize=16, pad=10)
    plt.xlabel("Distance (km)")
    plt.ylabel("Elevation (m)")
    plt.tight_layout()
    plt.savefig(out_png, bbox_inches="tight")
    plt.close()
    
# ------------------------------ CLI ------------------------------------

def main():
    ap = argparse.ArgumentParser(description="Compute one-direction ranch segments from KML (canonicalized orientation).")
    ap.add_argument("--kml", required=True, help="Input KML file with the track (gx:Track or LineString).")
    ap.add_argument("--names", required=True,
                    help="Comma-separated ranch names in your desired order, e.g. 'La Venada,El Mesquite,Casa Viejas,El Tule'.")
    ap.add_argument("--out-vertices", required=True, help="CSV of per-vertex data (trimmed one-direction).")
    ap.add_argument("--out-summary", required=True, help="CSV of segment summary per leg.")
    ap.add_argument("--plot", required=False, default=None, help="PNG elevation profile (optional).")
    ap.add_argument("--unit", default="miles", choices=("miles", "meters", "feet"), help="Distance unit for outputs.")
    ap.add_argument("--debug", action="store_true")
    args = ap.parse_args()

    wanted_names = [s.strip() for s in args.names.split(",") if s.strip()]
    if len(wanted_names) < 2:
        sys.exit("[ERROR] Please provide at least two names in --names (start and end).")

    # Parse KML
    track, pmarks = parse_kml_track(args.kml)
    n = len(track)
    if args.debug:
        print(f"[DEBUG] Track vertices: {n}")
        print(f"[DEBUG] Found placemarks: {list(pmarks.keys())}")

    # Snap names to nearest vertices
    name2idx = snap_names_to_track(track, wanted_names, pmarks)
    # If KML lacks Point placemarks, we cannot snap → fail clearly.
    missing = [nm for nm in wanted_names if nm not in name2idx]
    if missing:
        sys.exit(f"[ERROR] Could not find Point placemarks for: {', '.join(missing)}.\n"
                 f"Add Point markers named exactly as in --names to the KML, or we cannot place the vertical ticks and segments.")

    # Orientation selection
    reversed_needed, adj_map = choose_canonical_orientation(n, wanted_names, name2idx)
    if reversed_needed:
        track = list(reversed(track))
        # update mapping accordingly
        adj_map = {nm: (n-1 - idx) for nm, idx in name2idx.items()}
        name2idx = adj_map
        if args.debug:
            print("[DEBUG] Reversed track to honor requested ranch order.")
    else:
        if args.debug:
            print("[DEBUG] Kept original track orientation.")

    # Build DF on (possibly reversed) full track
    df_all = build_dataframe(track)

    # Trim to first → last name
    first_idx = name2idx[wanted_names[0]]
    last_idx  = name2idx[wanted_names[-1]]
    df = slice_one_direction(df_all, first_idx, last_idx)

    # Recompute/translate ranch tick positions into trimmed frame (in miles)
    # Use the original indices remapped to the trimmed domain
    ticks = []
    for nm in wanted_names:
        i = name2idx[nm]
        if i < first_idx or i > last_idx:
            continue
        xm = float(df_all.loc[i, "cum_dist_mi"] - df_all.loc[first_idx, "cum_dist_mi"])
        ticks.append((nm, xm))

    # Segment summaries (between consecutive snapped names)
    rows = []
    # make a mapping from global index → trimmed index
    # trimmed range is [first_idx, last_idx]
    for a, b in zip(wanted_names[:-1], wanted_names[1:]):
        i0 = name2idx[a]
        i1 = name2idx[b]
        if not (first_idx <= i0 <= last_idx and first_idx <= i1 <= last_idx):
            continue
        # convert to trimmed indices
        ti0 = i0 - first_idx
        ti1 = i1 - first_idx
        if ti1 <= ti0:
            # defensive; shouldn’t happen after orientation fix
            continue
        stats = leg_stats(df, ti0, ti1, args.unit)
        rows.append({
            "from": a,
            "to": b,
            f"distance_{args.unit}": round(stats["distance"], 2),
            "elev_gain_ft": round(stats["elev_gain_ft"], 1),
            "elev_loss_ft": round(stats["elev_loss_ft"], 1),
            "num_points": int(ti1 - ti0 + 1)
        })

    # Save CSVs
    df_out = df[["lat","lon","elev_ft","cum_dist_m","cum_dist_mi"]].copy()
    # Keep both m and mi so you have flexibility
    df_out.rename(columns={"cum_dist_m":"cum_dist_meters","cum_dist_mi":"cum_dist_miles"}, inplace=True)
    df_out.to_csv(args.out_vertices, index=False)
    pd.DataFrame(rows).to_csv(args.out_summary, index=False)

    if args.debug:
        print("[DEBUG] Ranch indices (canonical):")
        for nm in wanted_names:
            if nm in name2idx:
                gi = name2idx[nm]
                xm = float(df_all.loc[gi, "cum_dist_mi"])
                print(f"   {nm:12s}  idx={gi}  mile={xm:.2f}")
        print(f"[OK] Wrote vertices: {args.out_vertices}")
        print(f"[OK] Wrote summary : {args.out_summary}")

    # Plot
    if args.plot:
        plot_profile(df, ticks, args.plot, unit=args.unit)
        if args.debug:
            print(f"[OK] Wrote plot: {args.plot}")

if __name__ == "__main__":
    main()