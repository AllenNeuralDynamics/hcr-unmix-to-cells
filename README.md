# HCR Pairwise-Unmixing Cell Typing

Assigns cell types to HCR data using **interchangeable strategies**, one mouse at a time:

| Strategy | Flag | Output | Code |
|---|---|---|---|
| **MapMyCells** taxonomy mapping (ABC Atlas) | `--run-mapmycells` (default) | `results/mapmycells/` | [code/mapmycells/](code/mapmycells/) |
| **Tasic supercluster** matching (Tasic 2018 Smart-seq) | `--run-tasic-superclusters` | `results/tasic_superclusters/` | [code/tasic_superclusters/](code/tasic_superclusters/) |
| **Inhibitory GMM** on the mixed spot table (per-gene GMM gating, Slc17a7 QC, k-means clusters) | `--run-inhibitory-gmm` | `results/inhibitory_gmm/` | [code/inhibitory_gmm/](code/inhibitory_gmm/) |

Pass no flag → MapMyCells only (original behavior). Pass several → run several.

## Configurations

`--config <name>` loads a preset from [code/configs/](code/configs/); pick it from the app
panel's **config** list. An explicit argument always overrides the preset, so a preset can be
tweaked for one run (e.g. `--config p3_mixed_inhibitory_gmm --gmm-source processed`).

| Preset | Runs |
|---|---|
| `mapmycells` | MapMyCells on the pairwise-unmixing `filtered` tables |
| `tasic_superclusters` | Tasic supercluster matching, `log_zscore` normalization |
| `p3_mixed_inhibitory_gmm` | Inhibitory GMM on the mixed `all_spots` table (P3 data) |
| `p3_mixed_inhibitory_gmm_no_gfp` | Same, without GFP in the gate — use when GFP is pan-neuronal (e.g. 839909) — plus subclasses and fixed-k subtypes |

A preset is a JSON object of `run_capsule.py` options (`run_*`, `spots`, `normalization`,
`hcr_apply_pf`, `gmm_*`, plus an optional `description`); unknown keys are rejected. The app
panel's run flags have no default value so that an untouched panel does not override a preset.

---

## Change log

See [CHANGELOG.md](CHANGELOG.md) for full details.

| Date | Summary |
|---|---|
| 2026-09-25 | Added the **inhibitory GMM** strategy (ported from `hcr-pairwise-spot-unmixing`) and `--config` presets in `code/configs/`. The GMM strategy reads processed round assets, a spot-parquet asset or a pairwise-unmixing asset. |
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
The Tasic supercluster strategy needs the **Tasic 2018 Smart-seq VISp reference** (the
`mouse_VISp_2018-06-14_*` matrices), mounted at
`/root/capsule/data/tasic2018_VISp_gene_expression_matrices`. Override the location with the
`TASIC_SMARTSEQ_DIR` env var. See `SS_PATH` in
[code/tasic_superclusters/run_tasic_superclusters.py](code/tasic_superclusters/run_tasic_superclusters.py).

### Mixed spot tables (inhibitory GMM)
The inhibitory GMM strategy needs no pairwise-unmixing asset. It builds the all-rounds mixed
cell × gene table from whichever of these is mounted (`--gmm-source auto` tries them in order):

| `--gmm-source` | Mount | Reads |
|---|---|---|
| `pairwise` | `HCR_{mouse}_pairwise-unmixing_*` | `{mouse}_R{N}/mixed_{spots}_cell_by_gene.csv` (already built) |
| `spot_parquet` | an hcr-cache-spot-table asset (matched on `subject.json`) | `spots_R{N}.parquet` + `meta_R{N}.json` |
| `processed` | one `HCR_{mouse}_..._processed_...` asset per round | `processing_manifest.json` + `image_spot_spectral_unmixing/mixed_spots_R{N}.pkl` |

Spot-level sources use the pairwise capsule's rules: per-mouse gene deletions and gene-dict
overrides, and the Slc17a7 ≥ 200 intensity floor. `--gmm-spots filtered` also gates on the
pairwise `valid_spot` QC, so it needs the pairwise asset; processed and spot-parquet sources
support `all_spots` only.

---

## Usage

Run the capsule by passing the **mouse ID**:

The strategy flags take explicit boolean values (`true`/`false`) so they map cleanly to Code
Ocean app parameters; a bare flag (e.g. `--run-tasic-superclusters`) still means `true`.

```bash
# MapMyCells only (default when no flag is true)
python run_capsule.py --mouse-id 767018

# all-spots subset (MapMyCells)
python run_capsule.py --mouse-id 767018 --spots all_spots

# Tasic supercluster matching only
python run_capsule.py --mouse-id 767018 --run-tasic-superclusters true

# both strategies in one run
python run_capsule.py --mouse-id 767018 --run-mapmycells true --run-tasic-superclusters true

# P3: inhibitory GMM on the mixed spot table from processed round assets
python run_capsule.py --mouse-id 839909 --config p3_mixed_inhibitory_gmm
```

**Inhibitory GMM** writes under `/root/capsule/results/inhibitory_gmm/`. The folder and file
names match the pairwise capsule's `inhibitory_cells_mixed_*` output:

```
inhibitory_gmm/
├── inputs.json                                  source, per-round files, genes, parameters, counts
├── all_cells_mixed_{spots}/mixed_all_cells_{spots}.csv   all-rounds mixed cell x gene
└── inhibitory_cells_mixed_{spots}/
    ├── mixed_inhibitory_cells_{spots}.csv       GMM-selected, Slc17a7-QC'd cells
    ├── mixed_cluster_labels_{spots}.csv         cell_id, cluster (k-means)
    ├── mixed_sorted_cell_ids_{spots}.csv        heatmap row order
    ├── mixed_gmm_thresholds_{spots}.csv         per-gene GMM thresholds
    └── *.png                                    GMM grid / refit, Slc17a7 QC, clustered heatmap
```

Unlike the pairwise capsule's file, `mixed_cluster_labels_*.csv` has a `cell_id` column. The
pairwise file is indexed by heatmap row and only lines up with `*_sorted_cell_ids*.csv`.
GMM genes, the Slc17a7 cutoff, `k` and the heatmap clip are preset keys (`gmm_genes`,
`gmm_slc17a7_max`, `gmm_k`, `gmm_clip_max`).

**Subclasses and subtypes** (`gmm_subtypes: true`, on in `p3_mixed_inhibitory_gmm_no_gfp`) run on
the GMM-selected inhibitory cells and write `inhibitory_gmm/subtypes/`:

1. *Subclass* — `Pvalb` / `Sst` / `Vip` when the cell passes that gene's GMM threshold (counts
   under 5 zeroed). A cell positive for several goes to the gene it exceeds its own threshold by
   most (log2 fold) and is flagged `conflict`; a cell positive for none is `Other`.
2. *Subtype* — k-means within each subclass on log1p counts of every gene except
   `gmm_subtype_exclude` (GFP, Slc17a7, Gad2), with k from `gmm_subtype_k`
   (Pvalb 1, Sst 3, Vip 3, Other 2). Named `<subclass>_<gene>` after the expressed gene that
   most separates the sub-cluster from its siblings, or `<subclass>_only`.
3. *All-gene k-means* — k = `gmm_all_gene_k` (10) on log1p counts of every gene, for comparison.

```
inhibitory_gmm/subtypes/
├── cell_subtypes_clustering_genes.csv    cell_id, subclass, subtype, clustering-gene counts
├── cell_subtypes_all_genes.csv           same, all genes (GFP, Slc17a7, Gad2 first)
├── cell_all_gene_k10.csv                 + all-gene k-means label and silhouette
├── subtypes_cxg.png                      raw counts (10-200), clustering genes
├── subtypes_cxg_all_genes.png            raw counts, all genes (clustering genes in bold)
├── all_gene_k10_cxg.png                  raw counts, all-gene k-means
├── subtype_summary.csv                   cells, conflicts, silhouette, median counts per subtype
├── subtype_vs_all_gene_k10.csv           subtype x all-gene cluster counts
└── silhouette_by_subclass.csv            silhouette vs k (2-8) per subclass
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

---

## Cell typing table

After the strategies run, a single consolidated table of every HCR cell with its
assignment(s) is written to `results/cell_typing_table.csv`. The per-method tables
are merged (outer join) on the mouse-stripped cell id; if only one method ran, it
is a cleaned copy of that method's table. Columns from a method that did not type
a given cell are left blank for that row. When the inhibitory GMM ran, `gmm_cluster`
(k-means cluster) and `gmm_inhibitory` (`True` for GMM-selected cells) are appended, plus
`gmm_subclass` / `gmm_subtype` when the subtype step ran.

| Column | Source | Description |
|---|---|---|
| `cell_id` | both | Cell id with the `{mouse_id}_` prefix stripped; the merge key. |
| `mouse_id` | both | Mouse (subject) id. |
| `leiden_subclass` | TASIC | Inhibitory branch the cell was gated into (`Pvalb`, `Sst`, `Vip`, `Lamp5`). |
| `leiden_assignment` | TASIC | Leiden-named within-branch cluster assignment. |
| `leiden_confidence` | TASIC | Match confidence (Pearson correlation to the assigned cluster). |
| `mapmycells_class_name` | MapMyCells | ABC-atlas **class** label (coarsest taxonomy level). |
| `mapmycells_class_bootstrapping_probability` | MapMyCells | Fraction of bootstrap iterations supporting the class call. |
| `mapmycells_subclass_name` | MapMyCells | ABC-atlas **subclass** label. |
| `mapmycells_subclass_bootstrapping_probability` | MapMyCells | Bootstrap support for the subclass call. |
| `mapmycells_supertype_name` | MapMyCells | ABC-atlas **supertype** label. |
| `mapmycells_supertype_bootstrapping_probability` | MapMyCells | Bootstrap support for the supertype call. |
| `mapmycells_cluster_name` | MapMyCells | ABC-atlas **cluster** label (finest taxonomy level). |
| `mapmycells_cluster_alias` | MapMyCells | Numeric alias id for the assigned cluster. |
| `mapmycells_cluster_bootstrapping_probability` | MapMyCells | Bootstrap support for the cluster call. |