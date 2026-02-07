#!/usr/bin/env python3
"""
Step_9_Figures_pub.py
Publication-style Step 9 figure that matches the Step 8 look/flow.

Inputs (in --in_dir):
  - single_block_adjR2_<method>.csv
  - pairwise_varpart_<method>.csv
  - RUN_METADATA.txt  (optional; used for footer info)

Outputs:
  - <out_base>.png / .pdf

Example:
  python Step_9_Figures_pub.py \
    --in_dir "/.../Step9_VarPart_MIN_bray" \
    --method bray \
    --out_base "Step9_VarPart_bray" \
    --title "Drivers of Community Dissimilarity Across the Nogales Metacommunity — Step 9"
"""

import argparse, os, sys, textwrap
import pandas as pd
import matplotlib.pyplot as plt


# ---------- small utils ----------
def die(msg):
    print(f"[ERR] {msg}", file=sys.stderr); sys.exit(1)

def read_csv_required(p):
    if not os.path.exists(p): die(f"Missing file: {p}")
    try:
        return pd.read_csv(p)
    except Exception as e:
        die(f"Could not read {p}: {e}")

def read_meta(p):
    meta = {}
    if os.path.exists(p):
        with open(p, "r", encoding="utf-8") as fh:
            for line in fh:
                if "," in line:
                    k, v = line.strip().split(",", 1)
                    meta[k.strip()] = v.strip()
    return meta

def collapse_pairs(df_pairs):
    """Normalize/collapse pairwise varpart rows (E-G, E-R, G-R)."""
    if "pair" not in df_pairs:
        die("pairwise_varpart CSV must have a 'pair' column.")
    df = df_pairs.copy()
    df["pair"] = df["pair"].astype(str).str.strip()

    # required numeric columns
    for c in ["total_adjR2", "shared", "unique_X", "unique_Y"]:
        if c not in df.columns:
            df[c] = pd.NA

    # collapse duplicates by mean
    if df["pair"].duplicated().any():
        df = df.groupby("pair", as_index=False).agg({
            "total_adjR2":"mean", "shared":"mean", "unique_X":"mean", "unique_Y":"mean"
        })

    # keep order E-G, E-R, G-R if present
    wanted = ["E-G", "E-R", "G-R"]
    present = [p for p in wanted if p in set(df["pair"])]
    df = df.set_index("pair").loc[present].reset_index()
    return df


# ---------- plotting helpers ----------
def hbar_block(ax, labels, totals, best_individual, panel_title):
    """Single panel of horizontal bars with total + best-individual annotations."""
    # order by totals descending
    order = sorted(range(len(labels)), key=lambda i: (totals[i] if totals[i] is not None else -1), reverse=True)
    labels = [labels[i] for i in order]
    totals = [totals[i] for i in order]
    best_individual = [best_individual[i] for i in order]

    y = range(len(labels))
    heights = [t if t is not None else 0 for t in totals]
    ax.barh(y, heights)

    ax.set_yticks(y)
    ax.set_yticklabels(labels, fontsize=9)
    ax.invert_yaxis()
    ax.set_xlabel(r"$R^2$")
    ax.set_title(panel_title, fontsize=11)

    max_v = max([v for v in heights] + [0])
    ax.set_xlim(0, max_v*1.25 if max_v>0 else 1)

    # annotate
    for i, (t, b) in enumerate(zip(totals, best_individual)):
        if t is not None:
            ax.text(t + max_v*0.02 if max_v>0 else 0.02, i,
                    f"total R²={t:.2f}   best individual R²={b:.2f}" if b is not None else f"total R²={t:.2f}",
                    va="center", fontsize=8)


def main():
    ap = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=textwrap.dedent(__doc__)
    )
    ap.add_argument("--in_dir", required=True, help="Directory with Step 9 SAFE outputs")
    ap.add_argument("--method", default="bray", help="Method label (e.g., bray)")
    ap.add_argument("--out_base", required=True, help="Base name for outputs (no extension)")
    ap.add_argument("--title", default=None, help="Figure title")
    ap.add_argument("--width", type=float, default=8.5)
    ap.add_argument("--height", type=float, default=11.0)
    args = ap.parse_args()

    single_path = os.path.join(args.in_dir, f"single_block_adjR2_{args.method}.csv")
    pair_path   = os.path.join(args.in_dir, f"pairwise_varpart_{args.method}.csv")
    meta_path   = os.path.join(args.in_dir, "RUN_METADATA.txt")

    single = read_csv_required(single_path)
    pairs  = read_csv_required(pair_path)
    meta   = read_meta(meta_path)

    # ---- Single-block (Panel A) ----
    # Expect rows: baseline_null (0), E_only, G_only, R_only (maybe NA)
    s = single.copy()
    s["model"] = s["model"].astype(str).str.strip()
    get = lambda k: (float(s.loc[s["model"]==k, "adj_R2"].iloc[0])
                     if (k in set(s["model"]) and pd.notna(s.loc[s["model"]==k, "adj_R2"].iloc[0]))
                     else None)
    E = get("E_only"); G = get("G_only"); R = get("R_only")

    panelA_labels, panelA_totals, panelA_best = [], [], []
    if E is not None:
        panelA_labels.append("Environmental distance")
        panelA_totals.append(E)
        panelA_best.append(E)  # best individual for a single is itself
    if R is not None:
        panelA_labels.append("Combined resistance" if E is not None else "Resistance")
        panelA_totals.append(R)
        panelA_best.append(R)
    if G is not None:
        panelA_labels.append("Geographic distance")
        panelA_totals.append(G)
        panelA_best.append(G)

    # ---- Pairwise (Panel B) ----
    p = collapse_pairs(pairs)

    def best_of_singles_for_pair(tag):
        # E-G -> best among E and G (ignore missing)
        want = []
        if "E" in tag: want.append(E)
        if "G" in tag: want.append(G)
        if "R" in tag: want.append(R)
        want = [w for w in want if (w is not None)]
        return max(want) if want else None

    panelB_labels, panelB_totals, panelB_best = [], [], []
    for pair in p["pair"]:
        row = p.set_index("pair").loc[pair]
        tot = float(row["total_adjR2"]) if pd.notna(row["total_adjR2"]) else None
        panelB_labels.append({
            "E-G":"Geographic distance\nEnvironmental distance",
            "E-R":"Combined resistance\nEnvironmental distance",
            "G-R":"Geographic distance\nCombined resistance"
        }[pair])
        panelB_totals.append(tot)
        panelB_best.append(best_of_singles_for_pair(pair))

    # ---- Plot (two stacked panels to mimic Step 8 layout) ----
    fig = plt.figure(figsize=(args.width, args.height))
    gs = fig.add_gridspec(nrows=2, ncols=1, hspace=0.45)
    axA = fig.add_subplot(gs[0, 0])
    axB = fig.add_subplot(gs[1, 0])

    hbar_block(axA, panelA_labels, panelA_totals, panelA_best,
               panel_title="Panel A – Single-Predictor Models:")
    hbar_block(axB, panelB_labels, panelB_totals, panelB_best,
               panel_title="Panel B – Multi-Predictor Models:")

    # big title
    if args.title:
        fig.suptitle(args.title, y=0.98, fontsize=15, fontweight="bold")

    # footer like Step 8
    footer = []
    if meta:
        if "env_threshold" in meta or "kmax" in meta:
            footer.append(f"PCoA: threshold={meta.get('env_threshold','?')}, kmax={meta.get('kmax','?')}")
        if "env_file" in meta or "geo_file" in meta or "res_file" in meta:
            footer.append(f"E={os.path.basename(meta.get('env_file','?'))} | "
                          f"G={os.path.basename(meta.get('geo_file','?'))} | "
                          f"R={os.path.basename(meta.get('res_file','(none)'))}")
    if footer:
        fig.text(0.5, 0.02, "   ".join(footer), ha="center", fontsize=8)

    # save
    png = f"{args.out_base}.png"
    pdf = f"{args.out_base}.pdf"
    fig.savefig(png, dpi=300, bbox_inches="tight")
    fig.savefig(pdf, bbox_inches="tight")
    plt.close(fig)
    print(f"[OK] Wrote {png} and {pdf}")


if __name__ == "__main__":
    main()