"""
Build the final **cell typing table** at the root of the results folder.

This capsule can run two cell-typing strategies (TASIC superclusters and
MapMyCells / cell_type_mapper). Whichever ran, we emit one consolidated table
of every HCR cell with its assignment(s) at ``{results}/cell_typing_table.csv``:

* both ran   → the two per-method tables merged on the (mouse-stripped) cell id
* only one   → a cleaned copy of that method's table

Per-method column handling (see the module constants below)
-----------------------------------------------------------
TASIC (``hcr_assignments_leiden_named.csv``)
    cell_id, branch→leiden_subclass, leiden_cluster (drop),
    assignment→leiden_assignment, confidence→leiden_confidence, mouse_id

MapMyCells (``basic_results.csv``)
    cell_id kept as the join key; the label/consistency columns are dropped and
    every remaining column is prefixed ``mapmycells_``.

Cell-id alignment
-----------------
TASIC cell ids are ``{mouse_id}_{cell}`` (e.g. ``790322_6868``); MapMyCells cell
ids are the bare cell number (e.g. ``6868``). We strip the ``{mouse_id}_`` prefix
from the TASIC ids so both tables merge on the same key, and keep ``mouse_id`` as
its own column.
"""

from pathlib import Path

import pandas as pd

# --- output location ---------------------------------------------------------
CELL_TYPING_TABLE_NAME = "cell_typing_table.csv"

# --- TASIC ("hcr_assignments_leiden_named.csv") ------------------------------
TASIC_RESULT_SUBPATH = Path("tasic_superclusters") / "stage5" / "hcr_assignments_leiden_named.csv"
TASIC_RENAME = {
    "branch": "leiden_subclass",
    "assignment": "leiden_assignment",
    "confidence": "leiden_confidence",
}
TASIC_DROP = ["leiden_cluster"]

# --- MapMyCells ("basic_results.csv") ----------------------------------------
# The inhibitory-cells run is the one that aligns with TASIC's inhibitory branches.
MMC_RESULT_SUBPATH_TEMPLATE = "inhibitory_cells_{spots}/mapped_data/basic_results.csv"
MMC_PREFIX = "mapmycells_"
MMC_DROP = [
    "hierarchy_consistent",
    "class_label",
    "subclass_label",
    "supertype_label",
    "cluster_label",
]

MERGE_KEY = "cell_id"


def _strip_mouse_prefix(cell_id: str, mouse_id: str) -> str:
    """Return *cell_id* with a leading ``{mouse_id}_`` removed, if present.

    Falls back to the substring after the first underscore so a stray/renamed
    mouse id still yields the bare cell number.
    """
    cell_id = str(cell_id)
    prefix = f"{mouse_id}_"
    if cell_id.startswith(prefix):
        return cell_id[len(prefix):]
    if "_" in cell_id:
        return cell_id.split("_", 1)[1]
    return cell_id


def load_tasic_table(results_root: Path) -> pd.DataFrame | None:
    """Load + clean the TASIC assignments table, or return ``None`` if absent."""
    path = results_root / TASIC_RESULT_SUBPATH
    if not path.exists():
        print(f"[cell_typing_table] TASIC table not found: {path}")
        return None

    df = pd.read_csv(path)
    print(f"[cell_typing_table] Loaded TASIC table ({len(df)} rows): {path}")

    df = df.drop(columns=[c for c in TASIC_DROP if c in df.columns])
    df = df.rename(columns={k: v for k, v in TASIC_RENAME.items() if k in df.columns})

    # Strip the mouse-id prefix so the cell id merges with MapMyCells. Keep
    # mouse_id as a string so an outer merge's missing rows don't promote it to
    # float (which would render 790322 as "790322.0").
    if "mouse_id" in df.columns:
        df["mouse_id"] = df["mouse_id"].astype(str)
        df[MERGE_KEY] = [
            _strip_mouse_prefix(cid, mid)
            for cid, mid in zip(df[MERGE_KEY], df["mouse_id"])
        ]
    else:
        df[MERGE_KEY] = df[MERGE_KEY].map(lambda c: _strip_mouse_prefix(c, ""))
    df[MERGE_KEY] = df[MERGE_KEY].astype(str)
    return df


def load_mapmycells_table(results_root: Path, spots: str) -> pd.DataFrame | None:
    """Load + clean the MapMyCells basic_results table, or ``None`` if absent."""
    subpath = MMC_RESULT_SUBPATH_TEMPLATE.format(spots=spots)
    path = results_root / "mapmycells" / subpath
    if not path.exists():
        print(f"[cell_typing_table] MapMyCells table not found: {path}")
        return None

    df = pd.read_csv(path)
    print(f"[cell_typing_table] Loaded MapMyCells table ({len(df)} rows): {path}")

    df = df.drop(columns=[c for c in MMC_DROP if c in df.columns])
    # Prefix every column except the merge key.
    df = df.rename(
        columns={c: f"{MMC_PREFIX}{c}" for c in df.columns if c != MERGE_KEY}
    )
    df[MERGE_KEY] = df[MERGE_KEY].astype(str)
    return df


def build_cell_typing_table(
    results_root: Path,
    mouse_id: str,
    spots: str = "filtered",
) -> Path | None:
    """Write the consolidated cell typing table to *results_root*.

    Combines whichever of the TASIC / MapMyCells outputs exist. Returns the path
    to the written table, or ``None`` if neither method produced a table.
    """
    tasic = load_tasic_table(results_root)
    mmc = load_mapmycells_table(results_root, spots)

    if tasic is None and mmc is None:
        print(
            "[cell_typing_table] No TASIC or MapMyCells table found; "
            "skipping cell typing table."
        )
        return None

    if tasic is not None and mmc is not None:
        combined = tasic.merge(mmc, on=MERGE_KEY, how="outer")
        print(
            f"[cell_typing_table] Merged TASIC + MapMyCells on {MERGE_KEY!r}: "
            f"{len(combined)} rows."
        )
    else:
        combined = tasic if tasic is not None else mmc

    # Ensure a mouse_id column exists and is populated for every row (MapMyCells
    # rows have no mouse id of their own).
    if "mouse_id" not in combined.columns:
        combined["mouse_id"] = str(mouse_id)
    else:
        combined["mouse_id"] = combined["mouse_id"].fillna(str(mouse_id)).astype(str)

    # Column order: identity first, then TASIC (leiden_*), then MapMyCells.
    lead = [c for c in (MERGE_KEY, "mouse_id") if c in combined.columns]
    tasic_cols = [c for c in combined.columns if c.startswith("leiden_")]
    mmc_cols = [c for c in combined.columns if c.startswith(MMC_PREFIX)]
    ordered = lead + tasic_cols + mmc_cols
    ordered += [c for c in combined.columns if c not in ordered]
    combined = combined[ordered]

    out_path = results_root / CELL_TYPING_TABLE_NAME
    out_path.parent.mkdir(parents=True, exist_ok=True)
    combined.to_csv(out_path, index=False)
    print(
        f"[cell_typing_table] Wrote cell typing table ({len(combined)} rows, "
        f"{len(combined.columns)} cols): {out_path}"
    )
    return out_path
