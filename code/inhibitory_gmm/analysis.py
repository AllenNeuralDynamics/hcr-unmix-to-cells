"""Inhibitory cells by per-gene GMM thresholds, Slc17a7 leakage QC, then k-means clusters.

Ported from ``run_inhibitory_cell_analysis`` in hcr-pairwise-spot-unmixing
(``code/run_capsule.py``). Folder and file names are unchanged, so outputs line up with the
pairwise asset's ``inhibitory_cells_<spot type>_<subset>/`` folders. Differences:

* No ``HCRDataset``: columns are ordered from their ``R<N>-<channel>-<gene>`` labels.
* ``*_cluster_labels*.csv`` carries a ``cell_id`` column. The pairwise file is indexed by
  heatmap row, which only lines up with ``*_sorted_cell_ids*.csv``.
* The per-gene GMM thresholds are saved as ``*_gmm_thresholds*.csv``.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

import aind_hcr_qc.viz as viz
from aind_hcr_qc.utils.utils import saveable_plot

from .spot_tables import plain_gene_table, sort_columns_by_round_channel

DEFAULT_GENES = ["GFP", "Pvalb", "Sst", "Vip", "Npy", "Gad2"]


@saveable_plot()
def plot_slc17a7_qc(inh_cells_pivot, xlim=600, ylim=600, figsize=(18, 5)):
    """Slc17a7 histogram plus Gad2 and GFP against Slc17a7 for GMM-selected cells."""

    fig, axes = plt.subplots(1, 3, figsize=figsize)
    sns.histplot(data=inh_cells_pivot, x="Slc17a7", bins=1000, kde=True, ax=axes[0])
    axes[0].set_xlim(0, xlim)
    axes[0].set_title("Slc17a7 Distribution")
    for ax, gene in zip(axes[1:], ("Gad2", "GFP")):
        if gene in inh_cells_pivot.columns:
            ax.scatter(inh_cells_pivot["Slc17a7"], inh_cells_pivot[gene], alpha=0.5, s=5)
        ax.set_xlabel("Slc17a7")
        ax.set_ylabel(gene)
        ax.set_title(f"{gene} vs Slc17a7")
        ax.set_xlim(0, xlim)
        ax.set_ylim(0, ylim)
        ax.set_aspect("equal")
    plt.tight_layout()
    return fig


def run_inhibitory_cell_analysis(cxg_pivot: pd.DataFrame, output_dir: Path, mouse_id: str,
                                 table_type: str = "mixed_all_spots",
                                 genes: list[str] | None = None, slc17a7_max: int = 150,
                                 k: int = 20, clip_max: int = 200):
    """Select, QC, cluster and save inhibitory cells from an all-rounds cell x gene pivot.

    ``table_type`` is ``<spot type>_<subset>`` (e.g. ``mixed_all_spots``) and names the
    output folder ``inhibitory_cells_<table_type>`` and the file prefixes/suffixes.
    Returns ``(inhibitory pivot, labels frame, thresholds)``.
    """

    spot_type, _, filter_tag = table_type.partition("_")
    suffix = f"_{filter_tag}" if filter_tag else ""
    output_dir = Path(output_dir) / f"inhibitory_cells_{table_type}"
    output_dir.mkdir(parents=True, exist_ok=True)

    gene_cxg = plain_gene_table(cxg_pivot)
    requested = genes or DEFAULT_GENES
    genes = [g for g in requested if g in gene_cxg.columns]
    missing = [g for g in requested if g not in gene_cxg.columns]
    print(f"\n--- Inhibitory GMM analysis ({table_type}) ---")
    print(f"  Input: {gene_cxg.shape[0]:,} cells x {gene_cxg.shape[1]} genes; GMM genes {genes}")
    if missing:
        print(f"  WARNING: GMM genes not in this panel, skipped: {missing}")

    save = dict(save=True, output_dir=str(output_dir), show=False)
    thresholds, pivot_filtered = viz.run_inhibitory_gmm_thresholding(
        gene_cxg,
        inh_genes=genes,
        grid_save_kwargs={**save, "filename": f"{spot_type}_gmm_threshold_grid{suffix}"},
        refit_save_kwargs={**save, "filename": f"{spot_type}_gmm_refit{suffix}"},
    )
    inh_cells_plain = viz.filter_cells_by_gmm_thresholds(pivot_filtered, thresholds, logic="any")
    print(f"  Inhibitory cells identified: {len(inh_cells_plain):,}")
    inh_cells_pivot = cxg_pivot.loc[cxg_pivot.index.isin(inh_cells_plain.index)]

    if "Slc17a7" in inh_cells_plain.columns:
        plot_slc17a7_qc(inh_cells_plain, filename=f"{spot_type}_slc17a7_qc{suffix}", **save)
        keep = inh_cells_plain.index[inh_cells_plain["Slc17a7"] <= slc17a7_max]
        n_before = len(inh_cells_pivot)
        inh_cells_pivot = inh_cells_pivot.loc[inh_cells_pivot.index.isin(keep)].copy()
        print(f"  Slc17a7 > {slc17a7_max} cutoff removed {n_before - len(inh_cells_pivot):,} cells.")

    inh_cells_pivot = sort_columns_by_round_channel(inh_cells_pivot)
    # gene_sort=None: columns are already in round/channel order and carry R<N>-<chan>- labels.
    fig, cluster_labels, sorted_cell_ids = viz.plot_cell_x_gene_clustered(
        inh_cells_pivot,
        clip_range=(0, clip_max),
        fig_size=(8, 10),
        add_cluster_labels=True,
        title=f"{mouse_id} Inhibitory cells (GMM identified) \n k-means clusters k={k}",
        cluster_sort_gene="Gad2",
        gene_sort=None,
        k=k,
        filename=f"{table_type}_inhibitory_cells_clustered",
        **save,
    )
    plt.close(fig)

    labels = pd.DataFrame({"cell_id": list(sorted_cell_ids), "cluster": list(cluster_labels)})
    inh_cells_pivot.to_csv(output_dir / f"{spot_type}_inhibitory_cells{suffix}.csv")
    labels.to_csv(output_dir / f"{spot_type}_cluster_labels{suffix}.csv", index=False)
    labels[["cell_id"]].to_csv(output_dir / f"{spot_type}_sorted_cell_ids{suffix}.csv")
    if isinstance(thresholds, pd.DataFrame):
        thresholds.to_csv(output_dir / f"{spot_type}_gmm_thresholds{suffix}.csv")
    print(f"  Saved inhibitory cell outputs to {output_dir}")
    return inh_cells_pivot, labels, thresholds
