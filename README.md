# HCR Pairwise-Unmixing Cell Typing

Assigns cell types to pairwise-unmixed HCR data using **two interchangeable strategies**,
one mouse at a time:

| Strategy | Flag | Output | Code |
|---|---|---|---|
| **MapMyCells** taxonomy mapping (ABC Atlas) | `--run-mapmycells` (default) | `results/mapmycells/` | [code/mapmycells/](code/mapmycells/) |
| **Tasic supercluster** matching (Tasic 2018 Smart-seq) | `--run-tasic-superclusters` | `results/tasic_superclusters/` | [code/tasic_superclusters/](code/tasic_superclusters/) |
| **Reference compare** (collect per-mouse atlas-comparison CSVs) | `--run-reference-compare` | `results/reference_compare/` | [code/tasic_superclusters/collect_results.py](code/tasic_superclusters/collect_results.py) |

Pass no flag → MapMyCells only (original behavior). Flags combine; `--run-reference-compare`
is an aggregation step and does **not** require `--mouse-id`.

> `--run-reference-compare` consumes outputs of the reference-atlas compare flow
> (`run_atlas_compare.py` / `atlas_compare.py`), which is **not yet in this capsule** — until
> that upstream flow is added it will report "no comparison CSVs found".

---

## Change log

See [CHANGELOG.md](CHANGELOG.md) for full details.

| Date | Summary |
|---|---|
| 2026-06-29 | Added the **Tasic supercluster matching** strategy (extracted from `hcr-integrated-qc-capsule`). `run_capsule.py` now selects strategies via `--run-mapmycells` / `--run-tasic-superclusters`; MapMyCells code moved to `code/mapmycells/`, new pipeline in `code/tasic_superclusters/`. |
| 2026-05-28 | Added spot-mode input selection to `run_capsule.py` via `--spots {filtered\|all_spots}` (default: `filtered`), with backward-compatible legacy CSV fallbacks and spot-specific output folders |

---

## Requirements

### Data asset
Attach a **pairwise-unmixing output asset** to the capsule before running. The asset must be
mounted under `/root/capsule/data/` and follow the naming convention:

```
HCR_{mouse_id}_pairwise-unmixing_{YYYY-MM-DD_HH-MM-SS}/
├── inhibitory_cells_unmixed_filtered/
│   └── unmixed_inhibitory_cells_filtered.csv
├── all_cells_unmixed_filtered/
│   └── unmixed_all_cells_filtered.csv
├── inhibitory_cells_unmixed_all_spots/
│   └── unmixed_inhibitory_cells_all_spots.csv
└── all_cells_unmixed_all_spots/
    └── unmixed_all_cells_all_spots.csv
```

Example: `HCR_767018_pairwise-unmixing_2026-03-06_12-00-00`

Legacy layouts are still supported as fallbacks (for example,
`inhibitory_cells_unmixed/unmixed_inhibitory_cells.csv` and
`all_cells_unmixed/unmixed_all_cells.csv`).

### ABC Atlas asset (MapMyCells)
The ABC Atlas reference data must also be attached and available at `/root/capsule/data/abc_atlas/`.

### Tasic Smart-seq reference (Tasic superclusters)
> **TODO (data asset):** the Tasic supercluster strategy needs the **Tasic 2018 Smart-seq VISp
> reference** (the `mouse_VISp_2018-06-14_*` matrices). Attach it as a data asset and point the
> `TASIC_SMARTSEQ_DIR` env var at the mounted folder (default location:
> `/root/capsule/data/mouse_VISp_gene_expression_matrices_2018-06-14`). See `SS_PATH` in
> [code/tasic_superclusters/run_tasic_superclusters.py](code/tasic_superclusters/run_tasic_superclusters.py).

---

## Usage

Run the capsule by passing the **mouse ID**:

```bash
# MapMyCells only (default)
python run_capsule.py --mouse-id 767018

# all-spots subset (MapMyCells)
python run_capsule.py --mouse-id 767018 --spots all_spots

# Tasic supercluster matching only
python run_capsule.py --mouse-id 767018 --run-tasic-superclusters

# both strategies in one run
python run_capsule.py --mouse-id 767018 --run-mapmycells --run-tasic-superclusters
```

**MapMyCells** writes under `/root/capsule/results/mapmycells/` into spot-specific folders:
- `inhibitory_cells_filtered` / `all_cells_filtered`
- `inhibitory_cells_all_spots` / `all_cells_all_spots`

**Tasic superclusters** writes figures/tables under
`/root/capsule/results/tasic_superclusters/HCR_{mouse_id}/`, with heavy intermediate `.h5ad`
files under `/root/capsule/scratch/tasic_superclusters/HCR_{mouse_id}/`. The pipeline loads its
HCR query internally from the pairwise-unmixing asset and runs `batch_mode="none"` for a single
mouse. For the full option set (effect threshold, min cells, optional 10x-HMB matching), run
`code/tasic_superclusters/run_tasic_superclusters.py` directly — see its `--help`.

---

## Tuning parameters

Mapping parameters are currently set as defaults inside `run_capsule.py` and can be adjusted
manually in the `defaults` list (a future release will expose these as CLI flags):

| Parameter | Default | Description |
|---|---|---|
| `--bootstrap-iteration` | `100` | Number of bootstrap iterations |
| `--bootstrap-factor` | `1.0` | Fraction of markers sampled per bootstrap round |
| `--n-runners-up` | `2` | Number of runner-up cell types reported |
| `--num-workers` | `4` | Parallel workers for mapping |
| `--drop-layers` | `VISp6a VISp6b` | Taxonomy layers excluded from mapping |
| `--log-norm-data` | `True` | Apply expm1 to log-normalised input |
| `--generate-plots` | `True` | Produce QC and mapping plots |

Additional parameters (normalization, bootstrap seed, chunk size, etc.) are set in
`/root/capsule/code/params.json`.

---

## Output structure

```
scratch/{output_name}/
├── input_data/
│   └── input_cellxgene.h5ad
└── mapped_data/
    ├── basic_results.csv
    ├── extended_results.json
    └── mapped_cellxgene.h5ad
plots/
    └── *.png
```