"""Inhibitory subclasses and subtypes from the GMM-selected cells.

1. **Subclass**: each cell is ``Pvalb`` / ``Sst`` / ``Vip`` when it passes that gene's GMM
   threshold (counts below the noise floor zeroed, as in the GMM step). A cell positive for
   several goes to the gene it exceeds its own threshold by the most (log2 fold) and is
   flagged ``conflict``; a cell positive for none is ``Other``.
2. **Subtype**: k-means within each subclass on log1p counts of the clustering genes, with a
   fixed k per subclass (or, where none is given, the best silhouette if it reaches
   ``min_silhouette``). A sub-cluster is named ``<subclass>_<gene>`` after the expressed gene
   that most separates it from its sibling sub-clusters, or ``<subclass>_only`` when it is
   above its siblings on no gene.
3. **All-gene k-means**: one k-means on log1p counts of every gene, for comparison.

Outputs (in ``output_dir``): ``cell_subtypes_clustering_genes.csv``,
``cell_subtypes_all_genes.csv``, ``cell_all_gene_k<k>.csv``, the matching three raw-count
heatmaps, ``subtype_summary.csv``, ``subtype_vs_all_gene_k<k>.csv`` and
``silhouette_by_subclass.csv``.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_samples, silhouette_score
from sklearn.preprocessing import StandardScaler

DEFAULT_SUBCLASS_GENES = ("Pvalb", "Sst", "Vip")
DEFAULT_EXCLUDE = ("GFP", "Slc17a7", "Gad2")
LEAD_GENES = ("GFP", "Slc17a7", "Gad2")
NOISE_FLOOR = 5
K_RANGE = range(2, 9)
MIN_SILHOUETTE = 0.25  # weak-structure cutoff (Kaufman & Rousseeuw) for automatic k
MIN_MARKER_Z = 0.25  # a gene names a sub-cluster only if expressed there above the all-cell mean
VMIN, VMAX = 10, 200
COLORS = {"Pvalb": "#D95F5F", "Sst": "#3F78BF", "Vip": "#2E8B57", "Other": "#8c959f"}


def gmm_thresholds(thresholds: pd.DataFrame) -> pd.Series:
    """Raw-count threshold per gene; the last row wins (Gad2's sub-gene refit)."""
    return thresholds.drop_duplicates("gene", keep="last").set_index("gene")["threshold_raw"]


def assign_subclasses(genes: pd.DataFrame, thresholds: pd.Series,
                      subclass_genes=DEFAULT_SUBCLASS_GENES,
                      noise_floor: int = NOISE_FLOOR) -> pd.DataFrame:
    subclass_genes = [g for g in subclass_genes if g in genes.columns and g in thresholds.index]
    floored = genes[subclass_genes].where(genes[subclass_genes] >= noise_floor, 0)
    positive = floored.ge(thresholds[subclass_genes])
    fold = np.log2((floored + 1) / thresholds[subclass_genes])
    n_pos = positive.sum(axis=1)
    subclass = fold.where(positive, -np.inf).idxmax(axis=1).where(n_pos > 0, "Other")
    labels = pd.DataFrame({"subclass": subclass, "n_subclass_genes_positive": n_pos,
                           "conflict": n_pos > 1}, index=genes.index)
    for gene in subclass_genes:
        labels[f"{gene}_positive"] = positive[gene]
    return labels


def _subtype_names(name: str, subcluster: pd.Series, display: pd.DataFrame) -> dict[int, str]:
    if subcluster.nunique() == 1:
        return {int(subcluster.iloc[0]): name}
    # All-cell z-scores: within-subclass scaling would inflate near-zero genes.
    means = display.loc[subcluster.index].groupby(subcluster).mean()
    out = {}
    for c in means.index:
        diff = means.loc[c] - means.drop(index=c).mean()
        if diff.max() <= 0:
            out[c] = f"{name}_only"
            continue
        expressed = means.loc[c].drop(index=name, errors="ignore")
        candidates = diff[expressed[expressed >= MIN_MARKER_Z].index]
        out[c] = f"{name}_{candidates.idxmax()}" if len(candidates) else name
    for tag in {v for v in out.values() if list(out.values()).count(v) > 1}:
        gene = tag.split("_", 1)[1] if "_" in tag else tag
        tied = sorted((c for c in out if out[c] == tag),
                      key=lambda c: -means.loc[c, gene] if gene in means else 0)
        for rank, c in enumerate(tied):
            out[c] = f"{tag}{rank + 1}"
    return out


def cluster_subtypes(genes: pd.DataFrame, labels: pd.DataFrame, features: list[str],
                     fixed_k: dict[str, int] | None = None, k_range=K_RANGE,
                     min_silhouette: float = MIN_SILHOUETTE,
                     subclass_genes=DEFAULT_SUBCLASS_GENES):
    """Add ``subcluster``, ``silhouette`` and ``subtype`` to *labels*; return the silhouette table."""
    fixed_k = fixed_k or {}
    labels = labels.copy()
    labels["subcluster"] = 0
    labels["silhouette"] = np.nan
    display = pd.DataFrame(StandardScaler().fit_transform(np.log1p(genes[features].to_numpy(float))),
                           index=genes.index, columns=features)
    scores, chosen = {}, {}
    present = set(labels["subclass"])
    order = [s for s in [*subclass_genes, "Other"] if s in present]
    for name in order:
        idx = labels.index[labels["subclass"] == name]
        X = np.log1p(genes.loc[idx, features].to_numpy(float))
        fits = {k: KMeans(n_clusters=k, n_init=10, random_state=0).fit_predict(X)
                for k in k_range if k < len(idx)}
        scores[name] = {k: silhouette_score(X, lab) for k, lab in fits.items()}
        best = max(scores[name], key=scores[name].get) if scores[name] else 1
        auto = best if scores[name].get(best, 0) >= min_silhouette else 1
        k = int(fixed_k.get(name, auto))
        chosen[name] = k
        if k > 1:
            lab = fits[k] if k in fits else KMeans(n_clusters=k, n_init=10, random_state=0).fit_predict(X)
            labels.loc[idx, "subcluster"] = lab
            labels.loc[idx, "silhouette"] = silhouette_samples(X, lab)
        print(f"  {name:6s} n={len(idx):5d} best k={best} -> k={k}")
    names = {s: _subtype_names(s, labels.loc[labels["subclass"] == s, "subcluster"], display)
             for s in order}
    labels["subtype"] = [names[s][c] for s, c in zip(labels["subclass"], labels["subcluster"])]
    table = pd.DataFrame(scores).rename_axis("k")
    return labels, table, chosen, order


def all_gene_kmeans(genes: pd.DataFrame, columns: list[str], k: int):
    """k-means on log1p counts of *columns*; clusters named by their most elevated gene."""
    X = np.log1p(genes[columns].to_numpy(float))
    lab = KMeans(n_clusters=k, n_init=10, random_state=0).fit_predict(X)
    sil = silhouette_samples(X, lab)
    z = pd.DataFrame(StandardScaler().fit_transform(X), index=genes.index, columns=columns)
    means = z.groupby(lab).mean()
    tops = {c: (means.loc[c] - means.drop(index=c).mean()).idxmax() for c in means.index}
    order = sorted(means.index, key=lambda c: (tops[c], c))
    names = {c: f"k{i}_{tops[c]}" for i, c in enumerate(order)}
    rank = np.array([order.index(c) for c in lab])
    return pd.Series([names[c] for c in lab], index=genes.index), pd.Series(sil, index=genes.index), \
        rank, [names[c] for c in order]


def _plot_cxg(genes, order, bars, group_codes, group_names, columns, clustered, path, title):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import ListedColormap

    width = 1.2 + 0.85 * len(columns)
    fig, axes = plt.subplots(1, len(bars) + 1, figsize=(width, 10),
                             gridspec_kw={"width_ratios": [0.35] * len(bars) + [0.85 * len(columns)],
                                          "wspace": 0.06})
    for ax, (codes, colors) in zip(axes, bars):
        ax.imshow(codes[:, None], aspect="auto", interpolation="none",
                  cmap=ListedColormap(colors), vmin=0, vmax=len(colors) - 1)
        ax.set_xticks([])
        ax.set_yticks([])
    bounds = np.flatnonzero(np.diff(group_codes)) + 0.5
    centers = [(lo + hi) / 2 for lo, hi in zip(np.r_[-0.5, bounds], np.r_[bounds, len(order) - 0.5])]
    axes[0].set_yticks(centers, [group_names[group_codes[int(np.ceil(c))]] for c in centers],
                       fontsize=8)
    ax = axes[-1]
    counts = genes.loc[order, columns].to_numpy(float)
    im = ax.imshow(counts.clip(VMIN, VMAX), aspect="auto", interpolation="none", cmap="gray_r",
                   vmin=VMIN, vmax=VMAX)
    for b in bounds:
        ax.axhline(b, color="#D95F5F", lw=0.6)
    ax.set_xticks(range(len(columns)), columns, rotation=90)
    for tick, gene in zip(ax.get_xticklabels(), columns):
        tick.set_fontweight("bold" if gene in clustered else "normal")
    ax.set_yticks([])
    ax.set_title(title, fontsize=10)
    fig.colorbar(im, ax=ax, fraction=0.05, pad=0.02, label=f"Spot count ({VMIN}-{VMAX})")
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)


def run_subtypes(inhibitory_pivot: pd.DataFrame, thresholds: pd.DataFrame, output_dir: Path,
                 mouse_id: str, subclass_genes=DEFAULT_SUBCLASS_GENES, exclude=DEFAULT_EXCLUDE,
                 fixed_k: dict[str, int] | None = None, all_gene_k: int = 10,
                 plots: bool = True) -> pd.DataFrame:
    """Subclass + subtype labels, all-gene k-means, tables and heatmaps. Returns per-cell labels."""
    from .spot_tables import plain_gene_table

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    genes = plain_gene_table(inhibitory_pivot)
    features = [g for g in genes.columns if g not in set(exclude)]
    print(f"\n--- Inhibitory subtypes: {len(genes):,} cells, clustering genes {features} ---")

    labels = assign_subclasses(genes, gmm_thresholds(thresholds), subclass_genes)
    labels, silhouettes, chosen, order_sub = cluster_subtypes(genes, labels, features, fixed_k,
                                                              subclass_genes=subclass_genes)
    silhouettes.to_csv(output_dir / "silhouette_by_subclass.csv")

    sub_rank = {s: i for i, s in enumerate(order_sub)}
    subtype_order = sorted(labels["subtype"].unique(),
                           key=lambda t: (sub_rank[labels.loc[labels["subtype"] == t, "subclass"].iloc[0]], t))
    rank = labels["subtype"].map({t: i for i, t in enumerate(subtype_order)})
    rows = np.lexsort((-labels["silhouette"].fillna(0).to_numpy(), rank.to_numpy()))
    cells = labels.index[rows]

    lead = [g for g in LEAD_GENES if g in genes.columns]
    all_genes = lead + [g for g in features if g not in lead] + sorted(
        g for g in genes.columns if g not in lead and g not in features)
    ids = labels.loc[cells, ["subclass", "subtype"]]
    ids.index.name = "cell_id"
    ids.join(genes.loc[cells, features]).to_csv(output_dir / "cell_subtypes_clustering_genes.csv")
    ids.join(genes.loc[cells, all_genes]).to_csv(output_dir / "cell_subtypes_all_genes.csv")

    summary = labels.groupby(["subclass", "subtype"]).agg(
        n_cells=("conflict", "size"), n_conflict=("conflict", "sum"),
        mean_silhouette=("silhouette", "mean")).round(3)
    summary = summary.join(genes[all_genes].groupby([labels["subclass"], labels["subtype"]])
                           .median().add_suffix("_median_count"))
    summary.to_csv(output_dir / "subtype_summary.csv")
    print(summary[["n_cells", "n_conflict", "mean_silhouette"]].to_string())

    k_label, k_sil, k_rank, k_names = all_gene_kmeans(genes, all_genes, all_gene_k)
    k_rows = np.lexsort((-k_sil.to_numpy(), k_rank))
    k_cells = genes.index[k_rows]
    k_col = f"all_gene_k{all_gene_k}"
    k_table = pd.DataFrame({"subclass": labels["subclass"], "subtype": labels["subtype"],
                            k_col: k_label, "all_gene_silhouette": k_sil.round(3)})
    k_table.index.name = "cell_id"
    k_table.loc[k_cells].join(genes.loc[k_cells, all_genes]).to_csv(
        output_dir / f"cell_{k_col}.csv")
    pd.crosstab(k_table["subtype"], k_table[k_col]).to_csv(output_dir / f"subtype_vs_{k_col}.csv")
    print(f"  all-gene k={all_gene_k}: mean silhouette {k_sil.mean():.3f}")

    if plots:
        import matplotlib.pyplot as plt

        tab20 = plt.get_cmap("tab20")
        sub_colors = [COLORS.get(s, "#636c76") for s in order_sub]
        subtype_colors = [tab20(i % 20) for i in range(len(subtype_order))]
        sub_codes = labels.loc[cells, "subclass"].map(sub_rank).to_numpy()
        type_codes = rank.loc[cells].to_numpy()
        bars = [(sub_codes, sub_colors), (type_codes, subtype_colors)]
        n = f"n={len(cells):,}"
        _plot_cxg(genes, cells, bars, type_codes, subtype_order, features, set(features),
                  output_dir / "subtypes_cxg.png", f"{mouse_id} inhibitory subtypes ({n})")
        _plot_cxg(genes, cells, bars, type_codes, subtype_order, all_genes, set(features),
                  output_dir / "subtypes_cxg_all_genes.png",
                  f"{mouse_id} inhibitory subtypes, all genes ({n})\n"
                  f"Clusters from bold genes only: {', '.join(features)}")
        k_codes = k_rank[k_rows]
        k_bars = [(k_codes, [tab20(i % 20) for i in range(all_gene_k)]),
                  (labels.loc[k_cells, "subclass"].map(sub_rank).to_numpy(), sub_colors)]
        _plot_cxg(genes, k_cells, k_bars, k_codes, k_names, all_genes, set(all_genes),
                  output_dir / f"{k_col}_cxg.png",
                  f"{mouse_id} inhibitory cells, k-means on all genes, k={all_gene_k} ({n})")
    labels.attrs.update(chosen_k=chosen, features=features)
    return labels
