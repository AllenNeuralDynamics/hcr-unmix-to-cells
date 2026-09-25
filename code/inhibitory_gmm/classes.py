"""All-cell class labels.

* ``Inhibitory`` — the GMM-selected inhibitory set (after the Slc17a7 cutoff); takes precedence,
  so a Slc17a7+ inhibitory cell stays ``Inhibitory``.
* ``Excitatory`` — Slc17a7+ (GMM threshold fitted on all cells) and not inhibitory.
* ``Unassigned`` — neither (not ``None``: pandas reads that string back as missing).

Outputs in ``output_dir`` (``results/inhibitory_gmm/classes``): ``slc17a7_gmm_threshold.csv``
and its QC plot.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

CLASSES = ("Inhibitory", "Excitatory", "Unassigned")


def assign_classes(cell_ids: pd.Index, inhibitory: pd.Index,
                   slc17a7_positive: pd.Index) -> pd.DataFrame:
    inh = cell_ids.isin(inhibitory)
    slc = cell_ids.isin(slc17a7_positive)
    cls = np.select([inh, slc], ["Inhibitory", "Excitatory"], "Unassigned")
    return pd.DataFrame({"class": cls, "slc17a7_positive": slc}, index=cell_ids)


def fit_slc17a7_positive(genes: pd.DataFrame, output_dir: Path) -> tuple[pd.Index, pd.DataFrame]:
    """Slc17a7+ cells from a GMM threshold fitted on all cells (same routine as the GMM gate)."""
    import aind_hcr_qc.viz as viz

    save = dict(save=True, output_dir=str(output_dir), show=False)
    thresholds, floored = viz.run_inhibitory_gmm_thresholding(
        genes, inh_genes=["Slc17a7"], marker_gene=None,
        grid_save_kwargs={**save, "filename": "slc17a7_gmm_threshold"},
    )
    positive = viz.filter_cells_by_gmm_thresholds(floored, thresholds, logic="any").index
    thresholds.to_csv(output_dir / "slc17a7_gmm_threshold.csv")
    return positive, thresholds


def run_classes(genes: pd.DataFrame, inhibitory: pd.Index, output_dir: Path) -> pd.DataFrame:
    """Per-cell class for every cell in *genes*."""
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    slc_pos, _ = fit_slc17a7_positive(genes, output_dir)
    classes = assign_classes(genes.index, inhibitory, slc_pos)
    print(f"\n--- Cell classes ({len(classes):,} cells) ---")
    print(classes["class"].value_counts().reindex(CLASSES, fill_value=0).to_string())
    n_both = int((classes["class"].eq("Inhibitory") & classes["slc17a7_positive"]).sum())
    print(f"  Inhibitory and Slc17a7+ (kept Inhibitory): {n_both:,}")
    return classes
