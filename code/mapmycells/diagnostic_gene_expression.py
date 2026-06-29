"""
Diagnostic: Gene Expression Completeness Check
===============================================
Investigates missing / zero gene-expression counts in the pairwise-unmixing
cell-by-gene tables *before* any taxonomy mapping runs.

Root cause discovered
---------------------
For some imaging rounds (notably R2 and R3 in several HCR datasets) the
upstream segmentation pipeline writes cell label images as uint16, capping
cell IDs at 65 535 (0xFFFF).  The combined CSV is built from a reference
segmentation that uses uint32, yielding cell IDs up to ~127 000.  When the
multi-round tables are merged on ``cell_id``, all cells with cell_id > 65 535
receive zeros for the genes measured in the uint16-capped rounds.

What this script checks
-----------------------
Per HCR dataset (any folder matching ``HCR_*_pairwise-unmixing_*``):

  1. Per-round per-cell-format files (``*_R*/unmixed_cell_by_gene.csv`` and
     ``unmixed_cell_by_gene_filtered.csv``):
       - cell count, cell_id range, all-zero rows
  2. The merged table (``all_cells_unmixed/unmixed_all_cells.csv``):
       - per-round all-zero cell counts and fraction
       - cell_id range of the zero vs non-zero cells per round
       - cells with ≥ N rounds completely missing (default N=2)
  3. Identification of the uint16-cap boundary (65 535) per round
  4. Overall summary table across all datasets

Usage
-----
    python diagnostic_gene_expression.py [--data-root /root/capsule/data]
                                         [--multi-zero-thresh 2]
                                         [--output-csv /root/capsule/results/gene_expr_diagnostic.csv]
"""

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

UINT16_MAX = np.iinfo(np.uint16).max  # 65 535

# ── helpers ──────────────────────────────────────────────────────────────────

def _round_label(col: str) -> str:
    """Return the round prefix from a column name like 'R2-514-Hpse'."""
    return col.split("-")[0]


def _parse_rounds(gene_cols) -> dict:
    """Map round-prefix → list of column names."""
    rounds: dict = {}
    for col in gene_cols:
        r = _round_label(col)
        rounds.setdefault(r, []).append(col)
    return rounds


def check_per_round_files(pw_dir: Path) -> list[dict]:
    """
    Inspect the per-round unmixed_cell_by_gene[_filtered].csv files.
    Returns a list of dicts (one per round found).
    """
    results = []
    for rnd_dir in sorted(pw_dir.glob("*_R*")):
        for suffix, label in [
            ("unmixed_cell_by_gene.csv", "unfiltered"),
            ("unmixed_cell_by_gene_filtered.csv", "filtered"),
        ]:
            fpath = rnd_dir / suffix
            if not fpath.exists():
                continue

            df = pd.read_csv(fpath)
            id_col = "cell_id" if "cell_id" in df.columns else df.columns[0]
            cell_ids = df[id_col]
            cell_id_max = int(cell_ids.max())
            unique_cells = int(cell_ids.nunique())

            # All-zero check only applies to wide-format files
            gene_cols = [c for c in df.columns if c not in (id_col, "gene",
                         "round_chan_gene", "spot_id", "volume", "centroid",
                         "z_centroid", "y_centroid", "x_centroid", "z", "y",
                         "x", "z_spots", "y_spots", "x_spots", "spot_count",
                         "Unnamed: 0")]
            if gene_cols:
                all_zero_rows = int((df[gene_cols] == 0).all(axis=1).sum())
            else:
                all_zero_rows = None  # long-format file, skip

            cap_flag = cell_id_max <= UINT16_MAX

            results.append(dict(
                round=rnd_dir.name,
                file=label,
                n_rows=len(df),
                unique_cells=unique_cells,
                cell_id_min=int(cell_ids.min()),
                cell_id_max=cell_id_max,
                capped_at_uint16=cap_flag,
                all_zero_rows=all_zero_rows,
            ))
    return results


def check_combined_csv(combined_path: Path) -> dict:
    """
    Analyse the merged all_cells_unmixed/unmixed_all_cells.csv.
    Returns a summary dict plus a per-round breakdown DataFrame.
    """
    df = pd.read_csv(combined_path)
    gene_cols = [c for c in df.columns if c != "cell_id"]
    rounds = _parse_rounds(gene_cols)

    per_round = []
    round_zero_flags = {}
    for r, cols in sorted(rounds.items()):
        zero_mask = (df[cols] == 0).all(axis=1)
        round_zero_flags[r] = zero_mask
        zero_ids = df.loc[zero_mask, "cell_id"]
        nonzero_ids = df.loc[~zero_mask, "cell_id"]
        per_round.append(dict(
            round=r,
            n_genes=len(cols),
            all_zero_cells=int(zero_mask.sum()),
            pct_zero=round(100.0 * zero_mask.mean(), 1),
            zero_cell_id_max=int(zero_ids.max()) if len(zero_ids) else None,
            nonzero_cell_id_max=int(nonzero_ids.max()) if len(nonzero_ids) else None,
            nonzero_cell_id_min=int(nonzero_ids.min()) if len(nonzero_ids) else None,
            capped_at_uint16_boundary=(
                int(nonzero_ids.max()) <= UINT16_MAX if len(nonzero_ids) else False
            ),
        ))

    round_zero_df = pd.DataFrame(round_zero_flags)
    multi_zero = round_zero_df.sum(axis=1)

    summary = dict(
        path=str(combined_path),
        n_cells=len(df),
        n_genes=len(gene_cols),
        cell_id_min=int(df["cell_id"].min()),
        cell_id_max=int(df["cell_id"].max()),
        cells_zero_in_0_rounds=int((multi_zero == 0).sum()),
        cells_zero_in_1_round=int((multi_zero == 1).sum()),
        cells_zero_in_2plus_rounds=int((multi_zero >= 2).sum()),
        pct_zero_in_2plus_rounds=round(100.0 * (multi_zero >= 2).mean(), 1),
        per_round=pd.DataFrame(per_round),
    )
    return summary


# ── main diagnostic ───────────────────────────────────────────────────────────

def run_diagnostic(data_root: Path, multi_zero_thresh: int = 2) -> list[dict]:
    """
    Run the full diagnostic over all HCR datasets under *data_root*.
    Returns a list of per-dataset result dicts.
    """
    hcr_dirs = sorted(data_root.glob("HCR_*_pairwise-unmixing_*"))
    if not hcr_dirs:
        print(f"[WARN] No HCR_*_pairwise-unmixing_* directories found under {data_root}")
        return []

    all_results = []

    for hcr_dir in hcr_dirs:
        mouse_id = hcr_dir.name.split("_")[1]
        pw_dir = hcr_dir / "pairwise_unmixing"
        print(f"\n{'='*70}")
        print(f"Dataset : {hcr_dir.name}")
        print(f"Mouse ID: {mouse_id}")
        print(f"{'='*70}")

        if not pw_dir.exists():
            print("  [SKIP] No pairwise_unmixing sub-folder found.")
            all_results.append({"mouse_id": mouse_id, "status": "no_pairwise_unmixing"})
            continue

        # ── 1. Per-round file inspection ──────────────────────────────────────
        per_round_files = check_per_round_files(pw_dir)
        if per_round_files:
            print("\n[1] Per-round unmixed_cell_by_gene files")
            print(f"    {'Round':<20} {'File':<12} {'Rows':>8} {'UniqCells':>10} "
                  f"{'ID min':>8} {'ID max':>8} {'≤uint16?':>9} {'AllZeroRows':>12}")
            print("    " + "-"*95)
            for r in per_round_files:
                zero_str = str(r["all_zero_rows"]) if r["all_zero_rows"] is not None else "N/A"
                cap_str  = "YES ⚠" if r["capped_at_uint16"] else "no"
                print(f"    {r['round']:<20} {r['file']:<12} {r['n_rows']:>8,} "
                      f"{r['unique_cells']:>10,} {r['cell_id_min']:>8,} "
                      f"{r['cell_id_max']:>8,} {cap_str:>9} {zero_str:>12}")

        # ── 2. Combined CSV ───────────────────────────────────────────────────
        combined_path = pw_dir / "all_cells_unmixed" / "unmixed_all_cells.csv"
        if not combined_path.exists():
            print("\n  [SKIP] No all_cells_unmixed/unmixed_all_cells.csv found.")
            all_results.append({
                "mouse_id": mouse_id,
                "status": "no_combined_csv",
                "per_round_files": per_round_files,
            })
            continue

        combined = check_combined_csv(combined_path)
        pr_df = combined["per_round"]

        print(f"\n[2] Combined CSV: {combined_path.name}")
        print(f"    Total cells : {combined['n_cells']:,}")
        print(f"    Total genes : {combined['n_genes']}")
        print(f"    cell_id range: {combined['cell_id_min']:,} – {combined['cell_id_max']:,}")
        print(f"\n    Per-round all-zero cell analysis:")
        print(f"    {'Round':<6} {'Genes':>6} {'AllZeroCells':>14} {'%Zero':>7} "
              f"{'NonZero ID max':>16} {'Capped?':>9}")
        print("    " + "-"*65)
        for _, row in pr_df.iterrows():
            cap_str = "YES ⚠" if row["capped_at_uint16_boundary"] else "no"
            print(f"    {row['round']:<6} {row['n_genes']:>6} {row['all_zero_cells']:>14,} "
                  f"{row['pct_zero']:>6.1f}% {str(row['nonzero_cell_id_max']):>16} "
                  f"{cap_str:>9}")

        print(f"\n    Cross-round zero summary (cells with ≥ {multi_zero_thresh} rounds all-zero):")
        print(f"      No rounds zero      : {combined['cells_zero_in_0_rounds']:,}")
        print(f"      Exactly 1 round zero: {combined['cells_zero_in_1_round']:,}")
        print(f"      ≥ {multi_zero_thresh} rounds zero    : {combined['cells_zero_in_2plus_rounds']:,} "
              f"({combined['pct_zero_in_2plus_rounds']:.1f}%)")

        # ── 3. uint16 cap assessment ──────────────────────────────────────────
        capped_rounds = pr_df[pr_df["capped_at_uint16_boundary"] == True]["round"].tolist()
        print(f"\n[3] uint16 cap assessment (UINT16_MAX = {UINT16_MAX:,})")
        if capped_rounds:
            for _, row in pr_df[pr_df["capped_at_uint16_boundary"]].iterrows():
                n_affected = row["all_zero_cells"]
                pct = row["pct_zero"]
                print(f"    ⚠  Round {row['round']}: non-zero cell_id maxes at "
                      f"{row['nonzero_cell_id_max']:,} (≤ {UINT16_MAX:,}). "
                      f"{n_affected:,} cells ({pct:.1f}%) have ALL zeros for "
                      f"this round's {row['n_genes']} genes.")
            print(f"\n    DIAGNOSIS: The segmentation for round(s) {capped_rounds} used a "
                  f"uint16 label image.\n"
                  f"    Cells with cell_id > {UINT16_MAX:,} have no entry in those per-round\n"
                  f"    files, so they receive zeros when the combined table is built.\n"
                  f"    This is an UPSTREAM DATA ISSUE in the pairwise-unmixing pipeline,\n"
                  f"    not in this mapping capsule.")
        else:
            print(f"    No rounds capped at uint16 max — no obvious ID overflow detected.")

        result = {
            "mouse_id": mouse_id,
            "status": "ok",
            "n_cells": combined["n_cells"],
            "n_genes": combined["n_genes"],
            "cell_id_max": combined["cell_id_max"],
            "cells_zero_2plus_rounds": combined["cells_zero_in_2plus_rounds"],
            "pct_zero_2plus_rounds": combined["pct_zero_in_2plus_rounds"],
            "capped_rounds": capped_rounds,
            "per_round_df": pr_df,
            "per_round_files": per_round_files,
        }
        all_results.append(result)

    return all_results


def print_cross_dataset_summary(results: list[dict]) -> None:
    """Print a compact summary table across all datasets."""
    print(f"\n{'='*70}")
    print("CROSS-DATASET SUMMARY")
    print(f"{'='*70}")
    header = (f"{'Mouse ID':<12} {'Cells':>8} {'Affected':>10} {'%Affected':>10} "
              f"{'Capped rounds'}")
    print(header)
    print("-" * 70)
    for r in results:
        if r.get("status") not in ("ok",):
            print(f"  {r['mouse_id']:<12}  [skipped – {r['status']}]")
            continue
        capped = ", ".join(r["capped_rounds"]) or "none"
        print(f"  {r['mouse_id']:<12} {r['n_cells']:>8,} {r['cells_zero_2plus_rounds']:>10,} "
              f"{r['pct_zero_2plus_rounds']:>9.1f}%  {capped}")


def save_summary_csv(results: list[dict], output_path: Path) -> None:
    """Save a flat per-dataset + per-round CSV to *output_path*."""
    rows = []
    for r in results:
        if r.get("status") != "ok":
            continue
        for _, row in r["per_round_df"].iterrows():
            rows.append({
                "mouse_id": r["mouse_id"],
                "round": row["round"],
                "n_genes_in_round": row["n_genes"],
                "total_cells": r["n_cells"],
                "all_zero_cells": row["all_zero_cells"],
                "pct_zero": row["pct_zero"],
                "nonzero_cell_id_max": row["nonzero_cell_id_max"],
                "uint16_cap_detected": row["capped_at_uint16_boundary"],
                "cells_zero_in_2plus_rounds": r["cells_zero_2plus_rounds"],
                "pct_zero_2plus_rounds": r["pct_zero_2plus_rounds"],
            })
    if not rows:
        print("[WARN] No data to save.")
        return
    df_out = pd.DataFrame(rows)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    df_out.to_csv(output_path, index=False)
    print(f"\nSummary CSV saved → {output_path}")


# ── entry point ───────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description="Diagnostic: hunt for missing gene expression counts in HCR unmixing tables."
    )
    parser.add_argument(
        "--data-root",
        default="/root/capsule/data",
        help="Root folder that contains HCR_*_pairwise-unmixing_* subdirectories.",
    )
    parser.add_argument(
        "--multi-zero-thresh",
        type=int,
        default=2,
        help="Flag cells that are all-zero in this many or more rounds (default: 2).",
    )
    parser.add_argument(
        "--output-csv",
        default="/root/capsule/results/gene_expr_diagnostic.csv",
        help="Path to write the per-dataset/per-round summary CSV.",
    )
    args = parser.parse_args()

    data_root = Path(args.data_root)
    print(f"Data root : {data_root}")
    print(f"uint16 max: {UINT16_MAX:,}")

    results = run_diagnostic(data_root, multi_zero_thresh=args.multi_zero_thresh)
    print_cross_dataset_summary(results)
    save_summary_csv(results, Path(args.output_csv))


if __name__ == "__main__":
    main()
