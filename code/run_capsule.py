"""
Top-level run script – thin wrapper around run_taxonomy_mapper.main().

Pass --mouse-id <ID> and this script will locate the pairwise-unmixing asset
that was mounted under /root/capsule/data, derive the input CSV path, and
forward everything to run_taxonomy_mapper with sensible defaults.

Example
-------
python run_capsule.py --mouse-id 767018
# or with overrides:
python run_capsule.py --mouse-id 767018 --num-workers 8 --no-generate-plots
"""

import argparse
import sys
from pathlib import Path

from run_taxonomy_mapper import main as _mapper_main

DATA_ROOT = Path("/root/capsule/data")
UNMIXED_CSV_SUBPATH = "inhibitory_cells_unmixed/unmixed_inhibitory_cells.csv"
UNMIXED_ALL_CELLS_SUBPATH = "all_cells_unmixed/unmixed_all_cells.csv"


def find_pairwise_unmixing_asset(mouse_id: str, data_root: Path = DATA_ROOT) -> Path:
    """Return the pairwise-unmixing folder for *mouse_id*.

    Looks for a directory matching ``HCR_{mouse_id}_pairwise-unmixing_*``
    inside *data_root* and returns the first match.

    Raises
    ------
    FileNotFoundError
        If no matching folder is found.
    """
    pattern = f"HCR_{mouse_id}_pairwise-unmixing_*"
    matches = sorted(data_root.glob(pattern))
    if not matches:
        raise FileNotFoundError(
            f"No pairwise-unmixing asset found for mouse_id={mouse_id!r} "
            f"(searched {data_root / pattern})"
        )
    # Use the most-recent folder if multiple timestamps exist
    return matches[-1]


def find_all_cells_csv(asset_folder: Path) -> Path | None:
    """Return the all-cells CSV, checking two candidate locations.

    Candidates (in order):
    1. ``<asset_folder>/all_cells_unmixed/unmixed_all_cells.csv``
    2. ``<asset_folder>/unmixed_cell_by_gene_all_rounds.csv``

    Returns ``None`` if neither exists.
    """
    candidates = [
        asset_folder / UNMIXED_ALL_CELLS_SUBPATH,
        asset_folder / "unmixed_cell_by_gene_all_rounds.csv",
    ]
    for path in candidates:
        if path.exists():
            return path
    return None


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Run capsule – taxonomy mapping wrapper",
        add_help=True,
    )
    parser.add_argument(
        "--mouse-id",
        type=str,
        required=True,
        help="Mouse ID used to locate the pairwise-unmixing data asset "
             "(e.g. 767018 → HCR_767018_pairwise-unmixing_*/)",
    )
    parser.add_argument(
        "--output-dir",
        type=str,
        default="/root/capsule/results",
        help="Base directory for all outputs (default: /root/capsule/results)",
    )
    # Consume only known args; pass everything else straight to the mapper
    args, remaining = parser.parse_known_args()

    asset_folder = find_pairwise_unmixing_asset(args.mouse_id)
    output_name = asset_folder.name  # e.g. HCR_767018_pairwise-unmixing_2026-03-06_12-00-00

    print(f"Found asset : {asset_folder}")

    # --- resolve input CSVs ---------------------------------------------------
    inhibitory_csv = asset_folder / UNMIXED_CSV_SUBPATH
    all_cells_csv  = find_all_cells_csv(asset_folder)

    runs = [
        {
            "label":     "inhibitory_cells",
            "input_csv": inhibitory_csv,
            "output_dir": str(Path(args.output_dir) / "inhibitory_cells"),
        },
        {
            "label":     "all_cells",
            "input_csv": all_cells_csv,
            "output_dir": str(Path(args.output_dir) / "all_cells"),
        },
    ]

    for run in runs:
        input_csv = run["input_csv"]
        if input_csv is None or not input_csv.exists():
            print(f"WARNING: CSV not found, skipping {run['label']}: {input_csv}")
            continue

        print(f"\n{'='*60}")
        print(f"Run         : {run['label']}")
        print(f"Input CSV   : {input_csv}")
        print(f"Output dir  : {run['output_dir']}")
        print(f"Output name : {output_name}")
        print(f"{'='*60}\n")

        # Build argv for run_taxonomy_mapper, applying defaults then any user overrides
        defaults = [
            "--config", "/root/capsule/code/params.json",
            "--input-csv", str(input_csv),  # already resolved above
            "--output-name", output_name,
            "--output-dir", run["output_dir"],
            #"--log-norm-data", # using raw counts, dont need
            "--drop-layers", "VISp6a", "VISp6b",
            "--bootstrap-iteration", "100",
            "--bootstrap-factor", "1.0",
            "--n-runners-up", "2",
            "--num-workers", "4",
            "--overwrite-input",
            "--overwrite-mapping",
            "--overwrite-formatted",
            "--generate-plots",
        ]

        # remaining args from the user override / extend the defaults
        sys.argv = [sys.argv[0]] + defaults + remaining
        _mapper_main()
