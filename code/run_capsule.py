"""
Top-level run script – orchestrates the two cell-typing strategies in this
capsule, selected by boolean flags:

    --run-mapmycells            Taxonomy mapping (MapMyCells / cell_type_mapper).
                                Outputs → results/mapmycells/
    --run-tasic-superclusters   Tasic supercluster matching pipeline.
                                Outputs → results/tasic_superclusters/

Default (no flag): run MapMyCells only (preserves the original behavior).
Pass both flags to run both strategies in one invocation.

Mouse-id resolution (shared by both strategies)
-----------------------------------------------
Pipeline mode (automatic)
    If the folder ``/data/pipeline_data`` exists the script reads the single
    ``.txt`` file inside it and uses its stem as the mouse ID.

Standalone mode
    Otherwise pass ``--mouse-id <ID>`` on the command line.

The capsule processes one mouse at a time.

Examples
--------
python run_capsule.py --mouse-id 767018
python run_capsule.py --mouse-id 767018 --run-tasic-superclusters
python run_capsule.py --mouse-id 767018 --run-mapmycells --run-tasic-superclusters
python run_capsule.py --mouse-id 767018 --spots all_spots
"""

import argparse
import sys
from pathlib import Path

# The two strategies live in sibling subfolders. Add them to sys.path so their
# bare intra-package imports (``import config``, ``import taxonomy_mapper``,
# ``import cluster_validation_utils``) resolve. Heavy modules are imported
# lazily inside each branch so running one strategy does not pay the import cost
# of the other.
CODE_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(CODE_DIR / "mapmycells"))
sys.path.insert(0, str(CODE_DIR / "tasic_superclusters"))

DATA_ROOT          = Path("/root/capsule/data")
PIPELINE_DATA_ROOT = Path("/data/pipeline_data")
RESULTS_ROOT       = Path("/root/capsule/results")
SCRATCH_ROOT       = Path("/root/capsule/scratch")
PARAMS_PATH        = "/root/capsule/code/params.json"
SPOT_CHOICES = ("filtered", "all_spots", "nnls")

_TRUE_STRINGS = {"1", "true", "t", "yes", "y", "on"}
_FALSE_STRINGS = {"0", "false", "f", "no", "n", "off"}


def str2bool(value):
    """Parse a boolean from a string (for Code Ocean app boolean parameters).

    Accepts true/false/yes/no/1/0/on/off (case-insensitive). Lets the run flags
    be passed as ``--run-tasic-superclusters true`` rather than as bare presence
    flags, which the Code Ocean parameter UI cannot toggle cleanly.
    """
    if isinstance(value, bool):
        return value
    v = str(value).strip().lower()
    if v in _TRUE_STRINGS:
        return True
    if v in _FALSE_STRINGS:
        return False
    raise argparse.ArgumentTypeError(f"expected a boolean value, got {value!r}")

INHIBITORY_CSV_CANDIDATES_BY_SPOTS = {
    "filtered": [
        "inhibitory_cells_unmixed_filtered/unmixed_inhibitory_cells_filtered.csv",
        "inhibitory_cells_unmixed_filtered/unmixed_inhibitory_cells.csv",
        "inhibitory_cells_unmixed/unmixed_inhibitory_cells.csv",  # legacy
    ],
    "all_spots": [
        "inhibitory_cells_unmixed_all_spots/unmixed_inhibitory_cells_all_spots.csv",
        "inhibitory_cells_unmixed_all_spots/unmixed_inhibitory_cells.csv",
        "inhibitory_cells_unmixed/unmixed_inhibitory_cells.csv",  # legacy
    ],
    "nnls": [
        "inhibitory_cells_nnls/nnls_inhibitory_cells.csv",
    ],
}

ALL_CELLS_CSV_CANDIDATES_BY_SPOTS = {
    "filtered": [
        "all_cells_unmixed_filtered/unmixed_all_cells_filtered.csv",
        "all_cells_unmixed_filtered/unmixed_all_cells.csv",
        "all_cells_unmixed/unmixed_all_cells.csv",                # legacy
        "unmixed_cell_by_gene_all_rounds.csv",                    # legacy fallback
    ],
    "all_spots": [
        "all_cells_unmixed_all_spots/unmixed_all_cells_all_spots.csv",
        "all_cells_unmixed_all_spots/unmixed_all_cells.csv",
        "all_cells_unmixed/unmixed_all_cells.csv",                # legacy
        "unmixed_cell_by_gene_all_rounds.csv",                    # legacy fallback
    ],
    "nnls": [
        "all_cells_nnls/nnls_all_cells.csv",
    ],
}


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


# Which mounted derived asset holds the cell-by-gene tables for each spot mode.
# filtered/all_spots come from the pairwise-unmixing asset; nnls comes from the
# NNLS-only run's own derived asset (HCR_<mouse>_nnls-unmixing_*).
UNMIXING_ASSET_PATTERN_BY_SPOTS = {
    "filtered":  "HCR_{mouse}_pairwise-unmixing_*",
    "all_spots": "HCR_{mouse}_pairwise-unmixing_*",
    "nnls":      "HCR_{mouse}_nnls-unmixing_*",
}


def find_unmixing_asset(mouse_id: str, spots: str = "filtered",
                        data_root: Path = DATA_ROOT) -> Path:
    """Return the mounted unmixing asset folder for *mouse_id* and *spots* mode.

    Looks for the pattern in UNMIXING_ASSET_PATTERN_BY_SPOTS (pairwise-unmixing for
    filtered/all_spots, nnls-unmixing for nnls) inside *data_root*; returns the
    most-recent match.

    Raises
    ------
    FileNotFoundError
        If no matching folder is found.
    """
    pattern = UNMIXING_ASSET_PATTERN_BY_SPOTS.get(spots, UNMIXING_ASSET_PATTERN_BY_SPOTS["filtered"]).format(mouse=mouse_id)
    matches = sorted(data_root.glob(pattern))
    if not matches:
        raise FileNotFoundError(
            f"No unmixing asset found for mouse_id={mouse_id!r} spots={spots!r} "
            f"(searched {data_root / pattern})"
        )
    # Use the most-recent folder if multiple timestamps exist
    return matches[-1]


def find_pairwise_unmixing_asset(mouse_id: str, data_root: Path = DATA_ROOT) -> Path:
    """Backwards-compatible alias: the pairwise-unmixing (filtered/all_spots) asset."""
    return find_unmixing_asset(mouse_id, spots="filtered", data_root=data_root)


def find_first_existing_csv(asset_folder: Path, candidate_subpaths: list[str]) -> Path | None:
    """Return the first existing CSV from candidate relative subpaths.

    Returns ``None`` if none of the candidates exist.
    """
    for subpath in candidate_subpaths:
        path = asset_folder / subpath
        if path.exists():
            return path
    return None


def run_mapmycells(mouse_id: str, spots: str, output_root: Path, remaining: list[str]) -> None:
    """Run the MapMyCells (taxonomy mapping) strategy for one mouse.

    Outputs land under ``{output_root}/mapmycells/{inhibitory,all}_cells_{spots}``.
    """
    from run_taxonomy_mapper import main as _mapper_main  # lazy: heavy deps

    asset_folder = find_unmixing_asset(mouse_id, spots)
    output_name = asset_folder.name  # e.g. HCR_767018_pairwise-unmixing_2026-03-06_12-00-00
    mmc_root = output_root / "mapmycells"

    print(f"Found asset : {asset_folder}")
    print(f"Spot mode   : {spots}")

    # Some assets nest data under a pairwise_unmixing/ subfolder
    pairwise_subfolder = asset_folder / "pairwise_unmixing"
    csv_root = pairwise_subfolder if pairwise_subfolder.is_dir() else asset_folder
    print(f"CSV root    : {csv_root}")

    inhibitory_csv = find_first_existing_csv(
        csv_root, INHIBITORY_CSV_CANDIDATES_BY_SPOTS[spots],
    )
    all_cells_csv = find_first_existing_csv(
        csv_root, ALL_CELLS_CSV_CANDIDATES_BY_SPOTS[spots],
    )

    runs = [
        {
            "label":     "inhibitory_cells",
            "input_csv": inhibitory_csv,
            "output_dir": str(mmc_root / f"inhibitory_cells_{spots}"),
        },
        {
            "label":     "all_cells",
            "input_csv": all_cells_csv,
            "output_dir": str(mmc_root / f"all_cells_{spots}"),
        },
    ]

    for run in runs:
        input_csv = run["input_csv"]
        if input_csv is None or not input_csv.exists():
            print(f"WARNING: CSV not found, skipping {run['label']}: {input_csv}")
            continue

        print(f"\n{'='*60}")
        print(f"MapMyCells run : {run['label']}")
        print(f"Input CSV      : {input_csv}")
        print(f"Output dir     : {run['output_dir']}")
        print(f"Output name    : {output_name}")
        print(f"{'='*60}\n")

        # Build argv for run_taxonomy_mapper, applying defaults then user overrides.
        # output_name is "." so results land directly in output_dir (no extra subfolder)
        defaults = [
            "--config", PARAMS_PATH,
            "--input-csv", str(input_csv),
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


def run_tasic_superclusters(
    mouse_id: str,
    output_root: Path,
    scratch_root: Path,
    normalization: str = "log_zscore",
    hcr_apply_pf: bool = True,
    spots: str = "filtered",
) -> None:
    """Run the Tasic supercluster matching strategy for one mouse.

    HCR query data is loaded internally by the pipeline via
    aind_hcr_data_loader.get_hcr_dataset_pairwise(mouse_id, data_dir=DATA_ROOT).
    Outputs land directly under ``{output_root}/tasic_superclusters``; heavy
    intermediate .h5ad files go under ``{scratch_root}/tasic_superclusters``.
    (No per-mouse subfolder — the capsule runs one mouse at a time and the mouse
    id is already carried on the captured results data asset name.)

    normalization : per-cell normalization method passed to the pipeline
        ("log_zscore" default, "clr_shift", or "pflogpf"). hcr_apply_pf toggles
        the HCR depth-normalizing PF step for pflogpf.
    spots : which HCR inhibitory cell-by-gene table to load — "filtered"
        (default), "all_spots" (unfiltered), or "nnls" (per-spot NNLS
        optical-ghost table). Same axis as the MapMyCells --spots flag; Tasic
        is always inhibitory-only.

    NOTE(data-asset): requires the Tasic 2018 Smart-seq reference to be mounted
    and pointed at via the TASIC_SMARTSEQ_DIR env var (see SS_PATH in
    run_tasic_superclusters.py). Find/attach that asset before running.
    """
    from run_tasic_superclusters import main as _tasic_main  # lazy: heavy deps

    out_dir = output_root / "tasic_superclusters"
    scratch_dir = scratch_root / "tasic_superclusters"

    print(f"\n{'='*60}")
    print(f"Tasic superclusters : mouse {mouse_id}")
    print(f"Output dir          : {out_dir}")
    print(f"Scratch dir         : {scratch_dir}")
    print(f"Normalization       : {normalization} (hcr_apply_pf={hcr_apply_pf})")
    print(f"Spots               : {spots}")
    print(f"{'='*60}\n")

    # One mouse at a time → cross-mouse batch correction is a no-op; use "none".
    _tasic_main(
        mouse_ids=[mouse_id],
        output_dir=out_dir,
        scratch_dir=scratch_dir,
        batch_mode="none",
        normalization=normalization,
        hcr_apply_pf=hcr_apply_pf,
        all_spots=(spots == "all_spots"),
        spots=spots,
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Run capsule – cell-typing strategy orchestrator",
        add_help=True,
    )
    parser.add_argument(
        "--mouse-id",
        type=str,
        required=False,
        default=None,
        help="Mouse ID used to locate the pairwise-unmixing data asset "
             "(e.g. 767018 → HCR_767018_pairwise-unmixing_*/). "
             "Not required in pipeline mode (when /data/pipeline_data exists).",
    )
    parser.add_argument(
        "--output-dir",
        type=str,
        default=str(RESULTS_ROOT),
        help="Base directory for all outputs (default: /root/capsule/results)",
    )
    parser.add_argument(
        "--spots",
        type=str,
        choices=SPOT_CHOICES,
        default="filtered",
        help="Which HCR spot subset to use (applies to both strategies): "
             "'filtered' (default) or 'all_spots' (unfiltered). MapMyCells maps "
             "both the inhibitory and all-cells tables for this subset; Tasic "
             "uses the inhibitory table for this subset.",
    )
    parser.add_argument(
        "--run-mapmycells",
        type=str2bool,
        nargs="?",
        const=True,
        default=False,
        help="Run taxonomy mapping (MapMyCells) → results/mapmycells/ "
             "(true/false). Default strategy when no --run-* flag is true.",
    )
    parser.add_argument(
        "--run-tasic-superclusters",
        type=str2bool,
        nargs="?",
        const=True,
        default=False,
        help="Run Tasic supercluster matching → results/tasic_superclusters/ "
             "(true/false).",
    )
    parser.add_argument(
        "--normalization",
        type=str,
        default="log_zscore",
        choices=["log_zscore", "clr_shift", "pflogpf"],
        help="(Tasic) per-cell normalization before gene z-scoring: 'log_zscore' "
             "(default, original), 'clr_shift' (base log-norm + per-cell centering), "
             "or 'pflogpf' (PF -> log1p -> centering).",
    )
    parser.add_argument(
        "--hcr-apply-pf",
        type=str2bool,
        nargs="?",
        const=True,
        default=False,
        help="(Tasic, pflogpf only) apply the depth-normalizing PF step to HCR "
             "(true/false). False keeps HCR depth-free but still centers cells.",
    )
    # Consume only known args; pass everything else straight to the mapper
    args, remaining = parser.parse_known_args()

    # --- select strategies (default: MapMyCells only) ------------------------
    run_mmc = args.run_mapmycells
    run_tasic = args.run_tasic_superclusters
    if not run_mmc and not run_tasic:
        run_mmc = True
        print("No strategy flag given; defaulting to --run-mapmycells.")

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
            "(i.e. when /data/pipeline_data does not exist)"
        )

    output_root = Path(args.output_dir)

    if run_mmc:
        run_mapmycells(mouse_id, args.spots, output_root, remaining)

    if run_tasic:
        run_tasic_superclusters(
            mouse_id, output_root, SCRATCH_ROOT,
            normalization=args.normalization,
            hcr_apply_pf=args.hcr_apply_pf,
            spots=args.spots,
        )

    # --- consolidated cell typing table --------------------------------------
    # Merge whichever method(s) ran into one root-level table of every HCR cell
    # with its assignment(s): results/cell_typing_table.csv.
    from cell_typing_table import build_cell_typing_table

    build_cell_typing_table(output_root, mouse_id, spots=args.spots)
