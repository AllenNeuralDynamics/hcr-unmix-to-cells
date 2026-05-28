# HCR Pairwise-Unmixing Taxonomy Mapper

Maps inhibitory cell types from pairwise-unmixed HCR data against the ABC Atlas using MapMyCells.

---

## Change log

See [CHANGELOG.md](CHANGELOG.md) for full details.

| Date | Summary |
|---|---|
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

### ABC Atlas asset
The ABC Atlas reference data must also be attached and available at `/root/capsule/data/abc_atlas/`.

---

## Usage

Run the capsule by passing the **mouse ID** as the only required argument:

```bash
python run_capsule.py --mouse-id 767018
```

By default, the capsule uses the **filtered** spot subset. To run with all spots:

```bash
python run_capsule.py --mouse-id 767018 --spots all_spots
```

The script will:
1. Locate the matching pairwise-unmixing folder in `/root/capsule/data/`
2. Resolve input CSVs for both inhibitory and all-cells runs based on `--spots` (`filtered` by default)
3. Name the output after the asset folder (e.g. `HCR_767018_pairwise-unmixing_2026-03-06_12-00-00`)
4. Run all mapping steps and generate plots

Results are written under `/root/capsule/results/` into spot-specific folders:
- `inhibitory_cells_filtered` / `all_cells_filtered`
- `inhibitory_cells_all_spots` / `all_cells_all_spots`

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