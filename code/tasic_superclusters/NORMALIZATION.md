# Normalization exploration: PFlog1pPF / shifted-CLR

Status: **Phases 0–3 implemented** (opt-in, default behavior unchanged). Phases 4–5
(evaluation + decision) are not yet done.

## Motivation

A thread summarizing Booeshaghi & Pachter (bioRxiv 2022.05.06.490859v3) argues the
field-standard depth normalization — `log(x/s·K + 1)` with `K=10⁴` (Scanpy/Seurat
`normalize_total → log1p`) — is a poor depth normalizer, and that **PFlog1pPF**
(a.k.a. shifted CLR) is provably better: it is the unique transform satisfying rank
monotonicity, perturbation additivity (PCA compatibility), relabeling equivariance,
depth invariance, and calibration. Reported ~35/50 vs ~5/50 nearest-neighbor
preservation over 437 datasets.

We want to know whether it improves the **TASIC ↔ HCR** cross-platform matching in
the tasic-superclusters pipeline, so we implemented it behind a flag and can A/B it.

## The algorithm (confirmed from the paper's reference code)

```python
def pf(mtx):                              # proportional fitting to mean depth
    d = mtx.sum(1); return diags(d.mean() / d) @ mtx
x = pf(raw); x.data = np.log1p(x.data)   # PF → log1p
pflog1ppf = x - x.mean(axis=1)           # "second PF" = per-cell centering (CLR)
```

So **PFlog1pPF = PF(each cell → mean depth) → log1p → subtract per-cell mean**. The
pseudocount is just `log1p` (+1); the centering densifies the matrix.

## Current normalization in the repo (baseline)

In `run_stage1`, genes are first subset to the shared panel, then per platform:

| Matrix | Steps |
|--------|-------|
| TASIC (Smart-seq) | `normalize_total(target_sum=1e4)` → `log1p` |
| HCR (spot counts) | `log1p` only — **no** library-size normalization (targeted panel) |
| 10x-HMB | `log2` as-is (or `log1p` if raw) |

Everything is then **gene-wise z-scored** (`zscore_genes`), which is the
cross-platform comparison currency feeding batch correction (stage 2), centroid
matching (stage 3), and assignment (stage 5).

## Key considerations / caveats

- **Cell-wise vs gene-wise are orthogonal.** The thread's fix centers *cells*
  (across genes); the current pipeline standardizes *genes* (across cells). CLR
  slots in as `log → CLR(center cells) → [gene z-score]`; it does not replace the
  gene z-score. We keep gene z-scoring after every method so the **only variable is
  the per-cell normalization** (clean A/B; batch-correction currency held fixed).
- **Small panel.** CLR's guarantees are shown on whole-transcriptome data
  (thousands of genes). Here the shared panel is a *targeted* set (tens of genes),
  where the per-cell mean is taken over few genes — exactly the regime where CLR
  may help *or* hurt. This must be settled empirically (Phase 4).
- **HCR has no depth normalization by design.** Full PFlog1pPF would impose
  depth-invariance on HCR spot counts, whose per-cell totals partly reflect cell
  size / detection efficiency (possibly signal, possibly artifact). So HCR's PF step
  is made toggleable (`hcr_apply_pf`).
- **MapMyCells is out of scope.** That path uses `cell_type_mapper`'s own
  normalization tied to how its reference taxonomy was built; changing it would
  require matching the taxonomy's expectations.

## Implemented methods (`--normalization`)

All keep raw counts in `layers["raw"]` and apply gene z-scoring afterward.

| method | TASIC | HCR | 10x-HMB |
|--------|-------|-----|---------|
| `log_zscore` (default) | CP10k + log1p | log1p | log2 as-is / log1p |
| `clr_shift` | base log-norm + **center** | base log-norm + **center** | + center |
| `pflogpf` | PF → log1p → **center** | PF\* → log1p → **center** | center-only if log2, else PF→log1p→center |

\* HCR PF controlled by `hcr_apply_pf` (default `True` = faithful; `False` = keep HCR
depth-free but still center).

Code: `_proportional_fit`, `_clr_center`, `normalize_platform` in
`run_tasic_superclusters.py`; `NORMALIZATION_METHODS` enumerates the options.

## How to run

Capsule (via `run_capsule.py`): `--normalization {log_zscore,clr_shift,pflogpf}` and
`--hcr-apply-pf {true,false}`. Standalone tasic CLI: `--normalization ...` and
`--no-hcr-apply-pf`. The chosen method + `hcr_apply_pf` are recorded in the output
`run_params.json`. (Exposing these in the Code Ocean app panel is a separate GUI
step; the CLI/kwargs path works regardless.)

## Phase plan

- **Phase 0 — confirm algorithm.** ✅ Done (recipe above, from the paper's code).
- **Phase 1 — normalization switch.** ✅ `normalization`/`hcr_apply_pf` threaded
  through `run_stage1` → `main` → CLI → `run_capsule.py`; default unchanged.
- **Phase 2 — transforms.** ✅ `_proportional_fit`, `_clr_center`,
  `normalize_platform`.
- **Phase 3 — platform-specific handling.** ✅ HCR PF toggle; 10x log2 center-only;
  applied across the HCR+TASIC, 10x-only, and extra-10x paths.
- **Phase 4 — evaluation harness (TODO).** Run 790322 (and ≥1 more mouse) under each
  method into a `normalization_comparison/` dir and compare:
  - kNN preservation (neighbors retained vs raw),
  - cross-platform mixing (how well HCR/TASIC interleave pre/post batch correction),
  - stage-3 separability metrics (already computed),
  - cluster ARI vs the current pipeline,
  - stage-5 assignment confidence + marker separability.
- **Phase 5 — decide & document (TODO).** Keep opt-in, flip the default, or drop —
  based on Phase 4, written up here.

## Open questions for Phase 4+

- Does gene z-scoring *after* CLR blunt the benefit? Consider a "pure pflogpf"
  variant that skips gene z-score.
- Right PF target for HCR (mean depth vs none) given the targeted panel.
- Interaction with stage-2 batch correction — does better normalization reduce how
  much correction is needed (or make `batch_mode="none"` sufficient)?
