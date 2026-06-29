# Refactor plan: move Tasic cluster-matching code into `hcr-unmix-to-cells`

**Goal** — Move the Tasic supercluster matching pipeline out of
`hcr-integrated-qc-capsule/code/cluster_validation_tasic/` and into the
`hcr-unmix-to-cells` capsule, alongside the existing MapMyCells (taxonomy
mapping) workflow. The target capsule should be able to run **either** strategy,
selected by argparse boolean flags, **one mouse at a time**, writing to two
separate results subfolders:

| Strategy | Flag | Results subfolder |
|----------|------|-------------------|
| Taxonomy mapping (current behavior) | `--run-mapmycells` | `results/mapmycells/` |
| Tasic supercluster matching (new)   | `--run-tasic-superclusters` | `results/tasic_superclusters/` |

**Branches (created):**
- target `hcr-unmix-to-cells` → `feat/tasic-superclusters`
- source `hcr-integrated-qc-capsule` → `feat/extract-tasic-cluster-validation`

---

## Execution status (2026-06-29) — DONE (decisions confirmed)

Decisions: **D1-B** (full reorg into `mapmycells/`), **D2** intermediates → `scratch/`,
**D3** default → MapMyCells. Optional pieces (batch/aggregation/reference-atlas/old) **skipped**.

Completed in the target repo:
- Moved MapMyCells files (`config.py`, `taxonomy_mapper.py`, `diagnostic_gene_expression.py`,
  `run_taxonomy_mapper.py`) → `code/mapmycells/` (`params.json` kept at `code/` root).
- Copied Tasic core → `code/tasic_superclusters/`: `cluster_validation_utils.py`,
  `run_tasic_superclusters.py` (renamed from `run_tasic_cluster_uplevel.py`), 4 notes, all 7
  notebooks → `code/notebooks/tasic_superclusters/`.
- Refactored `run_tasic_superclusters.py`: `main()` accepts `mouse_ids` / `output_dir` /
  `scratch_dir`; heavy `.h5ad` → scratch; `SS_PATH` overridable via `TASIC_SMARTSEQ_DIR` env var.
- Rewrote `run_capsule.py`: `--run-mapmycells` / `--run-tasic-superclusters` flags, sys.path
  bootstraps for the two subfolders, MapMyCells → `results/mapmycells/`, Tasic →
  `results/tasic_superclusters/HCR_<id>/` (`batch_mode="none"` for single mouse).
- Updated README. All edited/new Python files pass `py_compile`.

**Remaining (cannot complete locally):**
- ⚠️ **Tasic Smart-seq reference data asset** — find/attach it; set `TASIC_SMARTSEQ_DIR`.
- ⚠️ **Dockerfile deps** — the target `environment/Dockerfile` is GUI-managed (has `# hash:` line),
  so it was **not** hand-edited. Verify these resolve in the base image, add via the CO GUI if
  missing: `polars`, `scanpy`, `anndata`, `seaborn`, `leidenalg`, `igraph`
  (`aind-hcr-data-loader` + `aind-hcr-qc` are already present).
- Source repo cleanup: decide whether to delete `cluster_validation_tasic/` after verification.

---

## 1. Source inventory (`cluster_validation_tasic/`, 19 files)

### Python (12)

| File | Lines | Internal deps | External deps | Verdict |
|------|------:|---------------|---------------|---------|
| `cluster_validation_utils.py` | 2095 | — | numpy, polars, matplotlib, anndata, scanpy, seaborn, sklearn, scipy | **MOVE — core lib.** Required by every pipeline script. |
| `run_tasic_cluster_uplevel.py` | 4804 | cluster_validation_utils | + aind_hcr_data_loader, aind_hcr_qc.viz | **MOVE — primary pipeline** (Stages 1–5). Rename → `run_tasic_superclusters.py`. Needs refactor (§4). |
| `atlas_compare.py` | ~ | — | numpy, pandas, matplotlib | **MOVE (optional)** — only used by the reference-atlas prep script. |
| `run_collect_reference_atlas_cellxgene.py` | ~900 | atlas_compare, cluster_validation_utils | + aind_hcr_data_loader, aind_hcr_qc.viz | **MOVE (optional)** — builds the 10x-HMB / MERFISH / Tasic reference matrices. Only needed for the optional 10x-HMB matching path; the core HCR↔Tasic path does not require it. |
| `run_tasic_cluster_coarsen_old.py` | 3950 | cluster_validation_utils | same as uplevel | **SKIP / archive** — explicitly the superseded "old" version of the uplevel pipeline. Git history preserves it; optionally drop in `tasic_superclusters/old/`. |
| `run_tasic_subclass_batch.py` | ~730 | cluster_validation_utils | + aind_hcr_data_loader, aind_hcr_qc.viz | **SKIP** — multi-mouse batch driver. You intend to run one mouse at a time, so this is not needed. |
| `compare_mice_results.py` | ~770 | — | pandas, seaborn, matplotlib, argparse | **SKIP for now (optional)** — cross-mouse Stage 5 aggregation. Useful later if you want to combine per-mouse runs; not needed for a single-mouse run. |
| `collect_results.py` | ~110 | — | pandas, argparse | **SKIP for now (optional)** — concatenates per-mouse atlas-comparison CSVs (pairs with the reference/atlas prep flow). |

### Notebooks (7) — **bring ALL over** (per request) → `code/notebooks/tasic_superclusters/`

| Notebook | Lines |
|----------|------:|
| `20260506-smartseq-tasic2018.ipynb` | 5023 |
| `20260603-post-tasic-matching-hcr-cellxgene.ipynb` | 3051 |
| `Explore_10x_VIS_Inhibitory_Depth.ipynb` | 1916 |
| `20260608-CompareAtlases-with-recluster.ipynb` | 1054 |
| `20260608-CompareAtlases.ipynb` | 420 |
| `2026-06-11-highly-variable-genes-and-unmixing.ipynb` | 136 |
| `ideal-data.ipynb` | 99 |

> Notebooks reference modules by bare name (e.g. `import cluster_validation_utils`).
> They will work if launched from inside `code/tasic_superclusters/`, or we add a
> short path-bootstrap cell. Note as a follow-up, not a blocker.

### Notes (4) → `code/tasic_superclusters/notes/`
`BATCH_CORRECTION_NOTES.md`, `IMPLEMENTATION_PLAN.md`,
`cross_platform_expression_comparison_pipeline.md`,
`hcr_tasic_matching_protocol.md` — **MOVE**, valuable design docs.

---

## 2. Dependency check (target environment)

The target `environment/Dockerfile` **already installs** the two AIND deps the
Tasic code needs:
- `aind-hcr-data-loader` ✓ (provides `get_hcr_dataset_pairwise`)
- `aind-hcr-qc` ✓ (provides `aind_hcr_qc.viz`)

(Pinned to different commits than the source capsule — fine, but worth a sanity
check that `aind_hcr_qc.viz` still exposes the symbols the Tasic code calls.)

**Not explicitly pinned in the target Dockerfile** (the Tasic code imports these):
`polars`, `scanpy`, `anndata`, `seaborn`, `leidenalg`/`igraph` (scanpy Leiden
`flavor="igraph"`), `scipy`, `scikit-learn`.

> **ACTION:** confirm these resolve in the target base image
> (`code-server-python-extensions-pack:code-server4.108.2python3.12.4...`). The
> source QC Dockerfile only had to add `polars` + `numpy` explicitly, implying
> scanpy/anndata/seaborn came from its base image — but the target uses a
> **different** base tag, so this must be verified and any missing packages added
> to the target Dockerfile (unlock if GUI-managed; see CLAUDE.md).

**Data assets** (critical, currently sourced from `scratch/`):
- `SS_PATH = scratch/mouse_VISp_gene_expression_matrices_2018-06-14` — the Tasic
  2018 Smart-seq VISp reference matrices. **Must be available in the target
  capsule** (attach as a Code Ocean data asset and point a config var at it).
- `TENX_HMB_DIR = scratch/reference_atlas_cellxgene/10x-hmb` — only for the
  optional 10x-HMB path.
- HCR query data via `get_hcr_dataset_pairwise(mouse_id, data_dir=/root/capsule/data)`
  — **same pairwise-unmixing asset pattern the existing MapMyCells flow already
  uses**, so this aligns cleanly with `run_capsule.py`'s asset discovery.

---

## 3. Proposed target layout (`code/`)

```
code/
  run                              # unchanged thin driver: python -u run_capsule.py "$@"
  run_capsule.py                   # REFACTORED orchestrator (adds the two strategy flags)
  params.json                      # mapmycells config (unchanged; optional tasic section)

  mapmycells/                      # Strategy 1 (existing taxonomy-mapping code, relocated)
    run_taxonomy_mapper.py
    taxonomy_mapper.py
    config.py
    diagnostic_gene_expression.py

  tasic_superclusters/             # Strategy 2 (new, from cluster_validation_tasic/)
    run_tasic_superclusters.py     # renamed from run_tasic_cluster_uplevel.py (refactored)
    cluster_validation_utils.py
    atlas_compare.py               # optional
    reference/                     # optional reference-atlas prep
      run_collect_reference_atlas_cellxgene.py
    old/                           # optional archive
      run_tasic_cluster_coarsen_old.py
    notes/                         # 4 design docs

  notebooks/
    demo-plot-taxonomy-mapping.ipynb   # existing
    tasic_superclusters/               # all 7 Tasic notebooks
  examples/                        # existing, unchanged
  old/                             # existing, unchanged
```

> **Decision point (D1):** Moving the existing MapMyCells files into a
> `mapmycells/` subfolder is cleaner but **breaks current imports/paths**
> (`from run_taxonomy_mapper import main`, `--config /root/capsule/code/params.json`).
> Two options:
> - **D1-A (recommended, lower risk):** leave the 4 MapMyCells files at `code/`
>   root; only the **new** Tasic code goes in `tasic_superclusters/`. Minimal
>   disruption to a working pipeline.
> - **D1-B:** fully reorganize both strategies into subfolders (matches the
>   "subfolders for two strategies" request) — requires fixing imports, the
>   `--config` path, and adding `sys.path` bootstraps. More work, more risk.
>
> The layout above shows D1-B. If you prefer D1-A, drop the `mapmycells/` move.

---

## 4. Refactor of the Tasic pipeline for single-mouse + results output

`run_tasic_cluster_uplevel.py` currently hard-codes (module level):
```python
DATA_DIR  = /root/capsule/data
SS_PATH   = /root/capsule/scratch/mouse_VISp_gene_expression_matrices_2018-06-14
OUT_ROOT  = /root/capsule/scratch/hcr_tasic_matching
TENX_HMB_DIR = /root/capsule/scratch/reference_atlas_cellxgene/10x-hmb
MOUSE_IDS = ["790322", "788406"]          # ← multi-mouse, hard-coded
```
and `main()` reads `MOUSE_IDS`/`OUT_ROOT` directly.

**Changes:**
1. **Parameterize `main()`** to accept `mouse_ids: list[str]` (or a single
   `mouse_id`) and `output_dir: Path`, instead of reading module globals.
   `run_stage1(...)` already takes `mouse_ids` — thread the value through.
2. **Single-mouse default:** orchestrator passes `[mouse_id]`. The batch-correction
   stage (Stage 2) is cross-*mouse*; with one mouse it's effectively a no-op /
   `batch_mode="none"` — verify Stage 2 handles a single mouse gracefully (it
   should: per-mouse centering of one group). Note in code.
3. **Output to results:** `OUT_ROOT` → `results/tasic_superclusters/<asset_or_mouse_name>/`.
   Heavy intermediate `.h5ad` files (`tasic_z`, `hcr_corrected`, `hcr_log`,
   `tasic_log`) are large — **decision (D2):** keep these in `scratch/` and write
   only final tables/figures to `results/`, or write everything to `results/`.
   Recommend: intermediates → scratch, summaries/figures/assignments → results.
4. **Tasic reference path** `SS_PATH` → a configurable arg / env var pointing at
   the attached Tasic Smart-seq data asset (not scratch).
5. Keep the existing rich argparse (`--batch-mode`, `--effect-threshold`,
   `--drop-minor-subclasses`, `--min-cells-per-cluster`, etc.) as pass-through
   options from the orchestrator.

---

## 5. Orchestrator (`run_capsule.py`) changes

Add two boolean flags and branch:

```python
parser.add_argument("--run-mapmycells", action="store_true",
                    help="Run taxonomy mapping (MapMyCells) → results/mapmycells/")
parser.add_argument("--run-tasic-superclusters", action="store_true",
                    help="Run Tasic supercluster matching → results/tasic_superclusters/")
```

- **mouse_id resolution stays as-is** (pipeline_data `.txt` stem, or `--mouse-id`).
- **MapMyCells branch:** existing two-pass logic (inhibitory_cells + all_cells),
  but output base becomes `results/mapmycells/...` instead of `results/...`.
- **Tasic branch:** call `run_tasic_superclusters.main(mouse_ids=[mouse_id],
  output_dir=results/tasic_superclusters/..., **passthrough_opts)`.

> **Decision point (D3): default behavior when neither flag is passed.**
> Options: (a) default to MapMyCells only (preserves today's behavior); (b) run
> both; (c) require at least one flag (error if none). Recommend **(a)** so
> existing pipeline triggers keep working unchanged.

---

## 6. Step-by-step execution checklist

1. ✅ Create feature branches in both repos.
2. **Target (`feat/tasic-superclusters`):**
   1. Create `code/tasic_superclusters/` (+ `notes/`, optional `reference/`, `old/`)
      and `code/notebooks/tasic_superclusters/`.
   2. Copy: `cluster_validation_utils.py`, `run_tasic_cluster_uplevel.py`
      (→ `run_tasic_superclusters.py`), the 4 notes, all 7 notebooks. Optional:
      `atlas_compare.py`, `run_collect_reference_atlas_cellxgene.py`,
      `run_tasic_cluster_coarsen_old.py`.
   3. Refactor `run_tasic_superclusters.py` per §4.
   4. (If D1-B) move MapMyCells files into `code/mapmycells/`, fix imports + paths.
   5. Refactor `run_capsule.py` per §5 (two flags, two results subfolders).
   6. Update `environment/Dockerfile` per §2 (add any missing libs; unlock if needed).
   7. Update `params.json` / add a tasic config block if useful.
   8. Smoke-test imports locally where possible (heavy deps may be cloud-only).
   9. Update README.
3. **Source (`feat/extract-tasic-cluster-validation`):**
   - Decide whether to **delete** `cluster_validation_tasic/` from the QC capsule
     (it now lives in `hcr-unmix-to-cells`) or leave it. Recommend a follow-up
     commit that removes it after the target is verified, with a pointer note.
4. Push both branches; sync + run in Code Ocean (attach Tasic Smart-seq + HCR
   pairwise-unmixing data assets); iterate.

---

## 7. Open decisions to confirm before I start editing

- **D1** — Reorganize MapMyCells into `mapmycells/` (D1-B) or leave at `code/` root (D1-A, recommended)?
- **D2** — Tasic intermediates (`.h5ad`) to `scratch/` (recommended) or `results/`?
- **D3** — Default when no flag passed: MapMyCells only (recommended), both, or require a flag?
- **Optional pieces** — include reference-atlas prep (`atlas_compare.py` +
  `run_collect_reference_atlas_cellxgene.py`), the `_old` coarsen script, and/or
  the multi-mouse aggregation scripts (`compare_mice_results.py`,
  `collect_results.py`)? Default plan: **skip all four**, bring only the core path.
- **Tasic reference data asset** — what is the Code Ocean data asset name/id for
  `mouse_VISp_gene_expression_matrices_2018-06-14` so we can wire `SS_PATH`?
