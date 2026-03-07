# HCR Pairwise-Unmixing Taxonomy Mapper

Maps inhibitory cell types from pairwise-unmixed HCR data against the ABC Atlas using MapMyCells.

---

## Requirements

### Data asset
Attach a **pairwise-unmixing output asset** to the capsule before running. The asset must be
mounted under `/root/capsule/data/` and follow the naming convention:

```
HCR_{mouse_id}_pairwise-unmixing_{YYYY-MM-DD_HH-MM-SS}/
└── inhibitory_cells_unmixed/
    └── unmixed_inhibitory_cells.csv   ← input to the mapper
```

Example: `HCR_767018_pairwise-unmixing_2026-03-06_12-00-00`

### ABC Atlas asset
The ABC Atlas reference data must also be attached and available at `/root/capsule/data/abc_atlas/`.

---

## Usage

Run the capsule by passing the **mouse ID** as the only required argument:

```bash
python run_capsule.py --mouse-id 767018
```

The script will:
1. Locate the matching pairwise-unmixing folder in `/root/capsule/data/`
2. Set the input CSV to `inhibitory_cells_unmixed/unmixed_inhibitory_cells.csv` inside that folder
3. Name the output after the asset folder (e.g. `HCR_767018_pairwise-unmixing_2026-03-06_12-00-00`)
4. Run all mapping steps and generate plots

Results are written to `/root/capsule/results/{output_name}/`.

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