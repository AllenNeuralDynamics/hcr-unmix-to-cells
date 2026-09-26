"""Entry point for the mixed-spot inhibitory GMM cell-typing method."""

from __future__ import annotations

import json
from pathlib import Path

from .spot_tables import build_cell_by_gene, discover_rounds


def main(mouse_id: str, data_root: Path, output_dir: Path, spots: str = "all_spots",
         source: str = "auto", genes: list[str] | None = None, slc17a7_max: int = 150,
         k: int = 20, clip_max: int = 200, subtypes: dict | None = None,
         classes: bool = False) -> Path:
    """Build the mixed all-rounds cell x gene table, then run the GMM analysis on it.

    Outputs go under ``output_dir`` (``results/inhibitory_gmm``):
    ``all_cells_mixed_<spots>/mixed_all_cells_<spots>.csv``, the
    ``inhibitory_cells_mixed_<spots>/`` folder, ``subtypes/`` when *subtypes* is given
    (keyword arguments for :func:`subtypes.run_subtypes`), and ``inputs.json`` (resolved
    inputs and parameters). With *classes* every cell gets a class (:func:`classes.run_classes`,
    Slc17a7 threshold in ``classes/``) and ``cell_classes.csv`` holds all per-cell labels.
    Returns the inhibitory-cells folder.
    """
    from .analysis import run_inhibitory_cell_analysis  # lazy: plotting / GMM deps

    output_dir = Path(output_dir)
    kind, rounds = discover_rounds(Path(data_root), mouse_id, spots=spots, source=source)
    print(f"\n{'=' * 60}")
    print(f"Inhibitory GMM : mouse {mouse_id}")
    print(f"Source         : {kind} ({len(rounds)} rounds: {[r.key for r in rounds]})")
    print(f"Spots          : mixed, {spots}")
    print(f"Output dir     : {output_dir}")
    print(f"{'=' * 60}\n")

    cxg = build_cell_by_gene(rounds, mouse_id, spots)
    table_type = f"mixed_{spots}"
    all_cells_dir = output_dir / f"all_cells_{table_type}"
    all_cells_dir.mkdir(parents=True, exist_ok=True)
    cxg.to_csv(all_cells_dir / f"mixed_all_cells_{spots}.csv")
    print(f"  All-rounds mixed table: {cxg.shape[0]:,} cells x {cxg.shape[1]} columns")

    inhibitory, labels, thresholds = run_inhibitory_cell_analysis(
        cxg, output_dir, mouse_id, table_type=table_type, genes=genes,
        slc17a7_max=slc17a7_max, k=k, clip_max=clip_max,
    )
    subtype_counts = None
    cell_types = None
    if subtypes is not None:
        from .subtypes import run_subtypes

        cell_types = run_subtypes(inhibitory, thresholds, output_dir / "subtypes", mouse_id,
                                  **subtypes)
        subtype_counts = {k: int(v) for k, v in cell_types["subtype"].value_counts().items()}
    class_counts = None
    if classes:
        from .classes import run_classes
        from .spot_tables import plain_gene_table
        from .subtypes import LEAD_GENES

        all_genes = plain_gene_table(cxg)
        per_cell = run_classes(all_genes, inhibitory.index, output_dir / "classes")
        if cell_types is not None:
            sub = cell_types[["subclass", "subtype", "silhouette"]]
            per_cell = per_cell.join(sub.rename(columns={"silhouette": "subtype_silhouette"}))
            per_cell["subtype_silhouette"] = per_cell["subtype_silhouette"].round(3)
        lead = [g for g in LEAD_GENES if g in all_genes.columns]
        per_cell = per_cell.join(all_genes[lead + [g for g in all_genes.columns if g not in lead]])
        per_cell.index.name = "cell_id"
        per_cell.to_csv(output_dir / "cell_classes.csv")
        class_counts = {k: int(v) for k, v in per_cell["class"].value_counts().items()}
    record = {
        "mouse_id": mouse_id,
        "source": kind,
        "spots": spots,
        "rounds": [{"round": r.key, "kind": r.kind, "source_asset": r.source_asset,
                    "path": str(r.path), "genes": r.genes} for r in rounds],
        "parameters": {"genes": genes, "slc17a7_max": slc17a7_max, "k": k,
                       "clip_max": clip_max, "subtypes": subtypes, "classes": classes},
        "n_cells": int(len(cxg)),
        "n_inhibitory": int(len(inhibitory)),
        "n_clusters": int(labels["cluster"].nunique()),
        "subtype_counts": subtype_counts,
        "class_counts": class_counts,
    }
    (output_dir / "inputs.json").write_text(json.dumps(record, indent=2))
    return output_dir / f"inhibitory_cells_{table_type}"
