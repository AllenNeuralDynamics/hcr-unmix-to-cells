"""
Top-level run script – thin wrapper around run_taxonomy_mapper.main().

Two run modes are supported:

Pipeline mode (automatic)
    If the folder ``/data/pipeline_data`` exists the script reads the
    single ``.txt`` file inside it and uses its stem (filename without the
    ``.txt`` extension) as the mouse ID.

    Expected layout::

        /data/pipeline_data/
            767018.txt          ← stem becomes mouse_id

Standalone mode
    When the pipeline_data folder is absent, pass ``--mouse-id <ID>`` on the
    command line.  The script will locate the pairwise-unmixing asset mounted
    under /root/capsule/data and forward everything to run_taxonomy_mapper
    with sensible defaults.

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

DATA_ROOT          = Path("/root/capsule/data")
PIPELINE_DATA_ROOT = Path("/data/pipeline_data")
UNMIXED_CSV_SUBPATH = "inhibitory_cells_unmixed/unmixed_inhibitory_cells.csv"
UNMIXED_ALL_CELLS_SUBPATH = "all_cells_unmixed/unmixed_all_cells.csv"


def resolve_mouse_id_from_pipeline(pipeline_root: Path = PIPELINE_DATA_ROOT) -> str:
    """Return the mouse ID encoded in the single .txt file inside *pipeline_root*.

    The convention is that exactly one file named ``<mouse_id>.txt`` lives in
    the pipeline_data folder.  Its stem is returned as the mouse ID string.

    Raises
    ------
    FileNotFoundError
        If no .txt file is found in *pipeline_root*.
    ValueError
        If more than one .txt file is found (ambiguous).
    """
    txt_files = sorted(pipeline_root.glob("*.txt"))
    if not txt_files:
        raise FileNotFoundError(
            f"Pipeline mode: no .txt file found in {pipeline_root}"
        )
    if len(txt_files) > 1:
        raise ValueError(
            f"Pipeline mode: expected exactly one .txt file in {pipeline_root}, "
            f"found {[f.name for f in txt_files]}"
        )
    mouse_id = txt_files[0].stem
    print(f"Pipeline mode: resolved mouse_id={mouse_id!r} from {txt_files[0]}")
    return mouse_id


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
        required=False,
        default=None,
        help="Mouse ID used to locate the pairwise-unmixing data asset "
             "(e.g. 767018 → HCR_767018_pairwise-unmixing_*/). "
             "Not required when running in pipeline mode "
             "(i.e. when /root/capsule/pipeline_data exists).",
    )
    parser.add_argument(
        "--output-dir",
        type=str,
        default="/root/capsule/results",
        help="Base directory for all outputs (default: /root/capsule/results)",
    )
    # Consume only known args; pass everything else straight to the mapper
    args, remaining = parser.parse_known_args()

    # --- resolve mouse_id (pipeline mode takes priority) ---------------------
    if PIPELINE_DATA_ROOT.exists():
        print(f"Pipeline mode detected: reading mouse_id from {PIPELINE_DATA_ROOT}")
        mouse_id = resolve_mouse_id_from_pipeline(PIPELINE_DATA_ROOT)
    elif args.mouse_id:
        mouse_id = args.mouse_id
        print(f"Standalone mode: using mouse_id={mouse_id!r} from --mouse-id argument")
    else:
        parser.error(
            "--mouse-id is required when not running in pipeline mode "
            "(i.e. when /root/capsule/pipeline_data does not exist)"
        )

    asset_folder = find_pairwise_unmixing_asset(mouse_id)
    output_name = asset_folder.name  # e.g. HCR_767018_pairwise-unmixing_2026-03-06_12-00-00

    print(f"Found asset : {asset_folder}")

    # Some assets nest data under a pairwise_unmixing/ subfolder
    pairwise_subfolder = asset_folder / "pairwise_unmixing"
    csv_root = pairwise_subfolder if pairwise_subfolder.is_dir() else asset_folder
    print(f"CSV root    : {csv_root}")

    # --- resolve input CSVs ---------------------------------------------------
    inhibitory_csv = csv_root / UNMIXED_CSV_SUBPATH
    all_cells_csv  = find_all_cells_csv(csv_root)

    runs = [
        {
            "label":     "inhibitory_cells",
            "input_csv": inhibitory_csv,
            "output_dir": str(Path(args.output_dir) / "inhibitory_cells"),
        },
        # {
        #     "label":     "all_cells",
        #     "input_csv": all_cells_csv,
        #     "output_dir": str(Path(args.output_dir) / "all_cells"),
        # },
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
        # output_name is "." so results land directly in output_dir (no extra subfolder)
        defaults = [
            "--config", "/root/capsule/code/params.json",
            "--input-csv", str(input_csv),  # already resolved above
            "--output-name", ".",
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
