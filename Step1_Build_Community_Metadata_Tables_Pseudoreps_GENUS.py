#!/usr/bin/env python3

"""
Step 1 — Build Community + Metadata Tables (GENUS version, pseudoreps preserved)

WHAT THIS DOES
• Prefer GENUS for the Taxon column:
  - If a 'Genus' (case-insensitive) column exists, we copy it into a new
    'Taxon' column and use that going forward (outputs will report genus).
  - Else we fall back to 'Taxon' if already present, or 'Species' / 'ASV' / 'OTU'.
• Preserve pseudoreps (no collapsing), just align community ↔ metadata.
• Safer sample column ordering:
  - Only keep community columns that appear in metadata SampleID.
  - Warn if metadata SampleIDs are missing from the community matrix.

INPUTS (required)
1) --community  CSV/XLSX of community data with taxa in the first column and sample IDs as columns.
   • This GENUS version prefers a 'Genus' column, but also accepts 'Taxon' / 'Species' / 'ASV' / 'OTU'.
   • Numeric sample columns may contain counts or relative abundances (zeros are OK).
2) --metadata   CSV/XLSX with at least:
   • SampleID  (string; must match column names in community)
   • Ranch/Site (called 'Ranch', 'Site', or variant)
   Optional: Latitude, Longitude, Elevation (used only if present).

KEY FLAGS
• --drop-negs              Drop negative controls by regex pattern.
• --neg-name-pattern REGEX  Pattern to identify negs in SampleID   [default: "blank|neg|control"]
• --require-geo             Error if ALL of Latitude/Longitude/Elevation are missing.
• --readme                  Save a README.txt with provenance and counts.

OUTPUTS (written to --outdir)
• community_aligned.csv   Taxon × SampleID matrix aligned to metadata (reports genus if available).
• metadata_aligned.csv    Filtered metadata aligned to the kept sample IDs.
• long_joined.csv         Long format: SampleID, Taxon, Abundance + metadata columns.
• README.txt              Optional provenance summary if --readme is provided.

TYPICAL USAGE
python Step1_Build_Community_Metadata_Tables_Pseudoreps_GENUS.py \
  --community /path/USE_GENUS_ONLY_matrix_by_trap_MEAN_REL_WeighteAvg_copy.csv \
  --metadata  /path/samples_env.csv \
  --outdir    STEP1_BY_TRAP_GENUS_MEANREL \
  --drop-negs --readme

NOTES
• This script does not collapse pseudoreps; it only aligns and cleans.
• If you need species→genus collapsing, do that upstream (you already did).
• If your SampleIDs include negatives (e.g., "NC", "R1", "blank"), tune --neg-name-pattern.
• If metadata has fewer IDs than community, only the intersection is kept and a warning is printed.

Author: Eric Olaf Stiner + Arlo
Version: 1.2 (2025-11-03)
"""

from __future__ import annotations
import argparse, os, sys, re
from datetime import datetime
from typing import List, Optional
import pandas as pd

# ---------- small utils ----------
def eprint(*a, **k): print(*a, **k, file=sys.stderr)

def load_table_auto(path: str) -> pd.DataFrame:
    """Load CSV/XLSX/XLS with sane fallbacks."""
    lower = path.lower()
    if lower.endswith(".xlsx") or lower.endswith(".xls"):
        return pd.read_excel(path)
    try:
        return pd.read_csv(path)
    except UnicodeDecodeError:
        return pd.read_csv(path, encoding="latin-1")

def find_col(df: pd.DataFrame, candidates: List[str], required: bool = True, nice_name: str = "") -> Optional[str]:
    """Find a column in df by candidate names (case-insensitive)."""
    lut = {c.lower(): c for c in df.columns}
    for cand in candidates:
        if cand.lower() in lut:
            return lut[cand.lower()]
    if required:
        want = nice_name or "/".join(candidates)
        raise ValueError(f"Required column not found: one of [{', '.join(candidates)}] for {want}")
    return None

def ensure_outdir(path: str):
    if path and not os.path.isdir(path):
        os.makedirs(path, exist_ok=True)

# ---------- main ----------
def main():
    ap = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="Assemble community + metadata tables with pseudoreps preserved. Prefers Genus for Taxon."
    )
    ap.add_argument("--community", required=True, help="Community matrix (CSV/XLSX/XLS). Contains taxa rows and sample columns.")
    ap.add_argument("--metadata",  required=True, help="Metadata (CSV/XLSX/XLS). Must contain SampleID and Ranch/Site.")
    ap.add_argument("--outdir",    required=True, help="Output directory.")
    ap.add_argument("--drop-negs", action="store_true", help="Drop negative controls by name pattern from both sides.")
    ap.add_argument("--neg-name-pattern", default=r"blank|neg|control",
                    help="Regex (case-insensitive) to match negatives in SampleID.")
    ap.add_argument("--require-geo", action="store_true",
                    help="Error if ALL of {Latitude, Longitude, Elevation} are missing.")
    ap.add_argument("--readme", action="store_true", help="Write README.txt with provenance and counts.")
    args = ap.parse_args()

    # Load
    comm = load_table_auto(args.community)
    meta = load_table_auto(args.metadata)

    # --- Choose Taxon column (PREFER GENUS) ---
    # 1) If Genus exists → use as Taxon (outputs will report genus).
    # 2) Else fall back to an existing Taxon column or Species/ASV/OTU.
    genus_col = find_col(comm, ["Genus", "genus"], required=False)
    if genus_col:
        comm = comm.rename(columns={genus_col: "Taxon"})
    else:
        taxon_col = find_col(comm, ["Taxon", "taxon", "Species", "species", "ASV", "OTU"], required=True, nice_name="Taxon/Species/ASV/OTU")
        if taxon_col != "Taxon":
            comm = comm.rename(columns={taxon_col: "Taxon"})

    # Make sure 'Taxon' is first column for clarity
    other_cols = [c for c in comm.columns if c != "Taxon"]
    comm = comm[["Taxon"] + other_cols]

    # Fill numeric NAs with 0 (common for abundance/rel-abundance)
    for c in other_cols:
        if pd.api.types.is_numeric_dtype(comm[c]):
            comm[c] = comm[c].fillna(0)

    # --- Metadata columns ---
    sample_col = find_col(meta, ["SampleID", "sampleid", "Sample", "sample", "Sample_ID"], required=True, nice_name="SampleID")
    ranch_col  = find_col(meta, ["Ranch", "ranch", "Site", "site"], required=True, nice_name="Ranch/Site")
    lat_col    = find_col(meta, ["Latitude", "latitude", "Lat", "lat"], required=False)
    lon_col    = find_col(meta, ["Longitude", "longitude", "Lon", "lon", "Long", "long"], required=False)
    elev_col   = find_col(meta, ["Elevation", "elevation", "Elev", "elev", "Altitude", "altitude"], required=False)

    if args.require_geo and (lat_col is None and lon_col is None and elev_col is None):
        raise ValueError("`--require-geo` set but none of Latitude/Longitude/Elevation are present in metadata.")

    # Normalize SampleID to string
    meta[sample_col] = meta[sample_col].astype(str)

    # Drop negatives if requested
    missing_in_comm = []
    not_merged = []
    if args.drop_negs:
        rx = re.compile(args.neg_name_pattern, flags=re.IGNORECASE)
        keep_mask = ~meta[sample_col].astype(str).str.contains(rx)
        dropped = (~keep_mask).sum()
        if dropped > 0:
            eprint(f"[INFO] Dropping {dropped} negative control samples by pattern: {args.neg_name_pattern}")
        meta = meta.loc[keep_mask].copy()

    # Align: order community sample columns to MATCH metadata order
    comm_sample_cols = [c for c in comm.columns if c != "Taxon"]
    meta_order = meta[sample_col].astype(str).tolist()
    keep_comm_cols = [sid for sid in meta_order if sid in comm_sample_cols]

    # warn if any metadata SampleIDs are missing in the community
    missing_in_comm = [sid for sid in meta_order if sid not in comm_sample_cols]
    if missing_in_comm:
        eprint(f"[WARN] {len(missing_in_comm)} SampleIDs in metadata not found in community. e.g., {missing_in_comm[:5]}")

    if not keep_comm_cols:
        raise ValueError("No overlapping SampleIDs between metadata and community matrix after filtering/alignment.")

    comm_aligned = comm[["Taxon"] + keep_comm_cols].copy()
    meta_aligned = meta.loc[meta[sample_col].astype(str).isin(keep_comm_cols)].copy()

    # Long table join
    long = comm_aligned.melt(id_vars=["Taxon"], var_name="SampleID", value_name="Abundance")
    long["SampleID"] = long["SampleID"].astype(str)
    long = long.merge(meta_aligned, left_on="SampleID", right_on=sample_col, how="left")

    # Any SampleIDs that failed to merge?
    not_merged = long[long[ranch_col].isna()]["SampleID"].unique().tolist()
    if not_merged:
        eprint(f"[WARN] {len(not_merged)} SampleIDs failed metadata merge. e.g., {not_merged[:5]}")

    # Summaries
    n_taxa = (comm_aligned["Taxon"].astype(str).str.strip() != "").sum()
    n_samp = len(keep_comm_cols)
    n_meta = meta_aligned.shape[0]

    # Write
    ensure_outdir(args.outdir)
    out_comm = os.path.join(args.outdir, "community_aligned.csv")
    out_meta = os.path.join(args.outdir, "metadata_aligned.csv")
    out_long = os.path.join(args.outdir, "long_joined.csv")
    comm_aligned.to_csv(out_comm, index=False)
    meta_aligned.to_csv(out_meta, index=False)
    long.to_csv(out_long, index=False)

    eprint(f"[OK] Wrote: {out_comm}")
    eprint(f"[OK] Wrote: {out_meta}")
    eprint(f"[OK] Wrote: {out_long}")
    eprint(f"[SUMMARY] Taxa (non-empty): {n_taxa} | Samples kept: {n_samp} | Metadata rows: {n_meta}")

    # README (optional)
    if args.readme:
        readme_path = os.path.join(args.outdir, "README.txt")
        lines = [
            "Elevational Metacommunity — Step 1 Provenance",
            f"Timestamp: {datetime.now().isoformat(timespec='seconds')}",
            "",
            "CLI:",
            "  " + " ".join(map(str, sys.argv)),
            "",
            "Inputs:",
            f"  community: {os.path.abspath(args.community)}",
            f"  metadata : {os.path.abspath(args.metadata)}",
            "",
            "Settings:",
            f"  drop_negs       : {args.drop_negs}",
            f"  neg_name_pattern: {args.neg_name_pattern}",
            f"  require_geo     : {args.require_geo}",
            "",
            "Outputs:",
            f"  community_aligned.csv  -> {out_comm}",
            f"  metadata_aligned.csv   -> {out_meta}",
            f"  long_joined.csv        -> {out_long}",
            "",
            "Counts:",
            f"  taxa (non-empty Taxon) : {n_taxa}",
            f"  samples kept           : {n_samp}",
            f"  metadata rows kept     : {n_meta}",
            "",
        ]
        if missing_in_comm:
            lines += [f"Metadata SampleIDs missing from community: {len(missing_in_comm)}",
                      "  e.g., " + ", ".join(missing_in_comm[:10]), ""]
        if not_merged:
            lines += [f"SampleIDs in community not found in metadata after filtering: {len(not_merged)}",
                      "  e.g., " + ", ".join(not_merged[:10]), ""]
        if (lat_col is None and lon_col is None and elev_col is None):
            lines.append("NOTE: No geo columns detected (Latitude/Longitude/Elevation).")
        with open(readme_path, "w") as fh:
            fh.write("\n".join(lines))
        eprint(f"[OK] Wrote: {readme_path}")

if __name__ == "__main__":
    main()