"""All-cell class labels and an excitatory-only clustering.

* ``Inhibitory`` — the GMM-selected inhibitory set (after the Slc17a7 cutoff), not Slc17a7+.
* ``Excitatory`` — Slc17a7+ (GMM threshold fitted on all cells) and not in the inhibitory set.
* ``Ambiguous`` — in the inhibitory set *and* Slc17a7+.
* ``Unassigned`` — neither (not ``None``: pandas reads that string back as missing).

``inhibitory_gate`` also flags cells that passed the inhibitory GMM gate but were cut for
Slc17a7 > the cutoff; they are ``Excitatory`` when Slc17a7+, so they can be separated.

Outputs in ``output_dir`` (``results/excitatory``): ``slc17a7_gmm_threshold.csv`` and its QC
plots, ``cell_excitatory_k<k>.csv``, ``excitatory_k<k>_summary.csv`` and
``excitatory_k<k>_cxg.png``.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

CLASSES = ("Inhibitory", "Excitatory", "Ambiguous", "Unassigned")
SILHOUETTE_SAMPLE = 10_000  # silhouette on a sample: excitatory sets reach ~60k cells


def assign_classes(cell_ids: pd.Index, inhibitory: pd.Index, gate: pd.Index,
                   slc17a7_positive: pd.Index) -> pd.DataFrame:
    inh = cell_ids.isin(inhibitory)
    slc = cell_ids.isin(slc17a7_positive)
    cls = np.select([inh & slc, inh, slc], ["Ambiguous", "Inhibitory", "Excitatory"], "Unassigned")
    return pd.DataFrame({"class": cls, "inhibitory_gate": cell_ids.isin(gate),
                         "slc17a7_positive": slc}, index=cell_ids)


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


def cluster_excitatory(genes: pd.DataFrame, k: int, exclude=("GFP",)):
    """k-means on log1p counts; clusters named ``Exc_1..k`` by median total spots, highest first."""
    from sklearn.cluster import KMeans
    from sklearn.metrics import silhouette_score

    features = [g for g in genes.columns if g not in set(exclude)]
    X = np.log1p(genes[features].to_numpy(float))
    model = KMeans(n_clusters=k, n_init=10, random_state=0).fit(X)
    lab = model.labels_
    sample = min(SILHOUETTE_SAMPLE, len(X))
    score = silhouette_score(X, lab, sample_size=sample, random_state=0) if k > 1 else float("nan")
    # Named by brightness, not by a marker: excitatory k-means splits mostly on detection depth.
    totals = genes[features].sum(axis=1).groupby(lab).median()
    order = list(totals.sort_values(ascending=False).index)
    names = {c: f"Exc_{i + 1}" for i, c in enumerate(order)}
    rank = np.array([order.index(c) for c in lab])
    distance = np.linalg.norm(X - model.cluster_centers_[lab], axis=1)
    return (pd.Series([names[c] for c in lab], index=genes.index), rank, distance,
            [names[c] for c in order], features, score)


def run_classes(genes: pd.DataFrame, inhibitory: pd.Index, gate: pd.Index, output_dir: Path,
                mouse_id: str, excitatory_k: int = 4, excitatory_exclude=("GFP",),
                plots: bool = True) -> pd.DataFrame:
    """Per-cell class (+ excitatory cluster for Excitatory cells) for every cell in *genes*."""
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    slc_pos, _ = fit_slc17a7_positive(genes, output_dir)
    classes = assign_classes(genes.index, inhibitory, gate, slc_pos)
    print(f"\n--- Cell classes ({len(classes):,} cells) ---")
    print(classes["class"].value_counts().reindex(CLASSES, fill_value=0).to_string())

    exc = genes.loc[classes["class"] == "Excitatory"]
    classes["excitatory_cluster"] = pd.NA
    if len(exc) > excitatory_k:
        names, rank, distance, order, features, score = cluster_excitatory(
            exc, excitatory_k, excitatory_exclude)
        classes.loc[exc.index, "excitatory_cluster"] = names
        print(f"  excitatory k={excitatory_k} on {features}: silhouette {score:.3f} "
              f"(sample of {min(SILHOUETTE_SAMPLE, len(exc)):,})")
        rows = np.lexsort((distance, rank))
        cells = exc.index[rows]
        col = "excitatory_cluster"
        table = pd.DataFrame({col: names}).loc[cells].join(exc.loc[cells])
        table.index.name = "cell_id"
        table.to_csv(output_dir / f"cell_excitatory_k{excitatory_k}.csv")
        summary = exc.groupby(names).median().add_suffix("_median_count")
        summary.columns.name = None
        summary.insert(0, "n_cells", names.value_counts())
        summary.rename_axis(col).to_csv(output_dir / f"excitatory_k{excitatory_k}_summary.csv")
        print(summary.to_string())
        if plots:
            import matplotlib.pyplot as plt

            from .subtypes import LEAD_GENES, _plot_cxg

            tab20 = plt.get_cmap("tab20")
            codes = rank[rows]
            lead = [g for g in LEAD_GENES if g in exc.columns]
            columns = lead + [g for g in exc.columns if g not in lead]
            _plot_cxg(exc, cells, [(codes, [tab20(i % 20) for i in range(excitatory_k)])], codes,
                      order, columns, set(features), output_dir / f"excitatory_k{excitatory_k}_cxg.png",
                      f"{mouse_id} excitatory cells, k-means k={excitatory_k} (n={len(cells):,})\n"
                      f"Clusters from bold genes only")
    return classes
