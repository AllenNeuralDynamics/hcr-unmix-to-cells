# `cell_typing_table.csv` — inhibitory GMM cell typing

One row per segmented HCR cell of one mouse. Written to the results root by the
`hcr-unmix-to-cells` capsule when the inhibitory GMM strategy runs with classes
(preset `p3_mixed_inhibitory_gmm_no_gfp`). Parameters below are that preset's; the values
used in a run are in `inhibitory_gmm/inputs.json`.

## Procedure

1. **Cell × gene counts.** Mixed (pre-unmixing) spots from every round, no per-spot QC gate
   (`all_spots`), counted per cell and gene. Slc17a7 spots with intensity below 200 are
   dropped, as in the pairwise-unmixing capsule. Table:
   `inhibitory_gmm/all_cells_mixed_all_spots/mixed_all_cells_all_spots.csv`.
2. **GMM thresholds.** Counts below 5 are set to 0. For each gate gene (Pvalb, Sst, Vip,
   Npy, Gad2), a 2-component Gaussian mixture (tied covariance) is fitted to log2(count + 1)
   of the cells with a non-zero count. The threshold is the count where the posterior of
   the higher component reaches 0.5. Gad2 is refitted on only the cells positive for Pvalb,
   Sst, Vip or Npy, and that threshold is used. Thresholds:
   `inhibitory_gmm/inhibitory_cells_mixed_all_spots/mixed_gmm_thresholds_all_spots.csv`.
3. **Inhibitory set.** A cell passes the gate when it exceeds the threshold of *any* gate
   gene. Gate-passing cells with more than 150 Slc17a7 spots are removed. The rest form the
   inhibitory set.
4. **Slc17a7 threshold.** The same GMM fit (step 2) on Slc17a7 alone, over all cells. A cell
   above it is Slc17a7-positive. Threshold and fit plot: `inhibitory_gmm/classes/`.
5. **Class.** First matching rule wins:
   `Inhibitory` (in the inhibitory set, even when Slc17a7-positive) →
   `Excitatory` (Slc17a7-positive) → `Unassigned`.
6. **Subclass** (inhibitory cells only). `Pvalb`, `Sst` or `Vip` when the cell reaches that
   gene's GMM threshold (counts below 5 set to 0). A cell positive for several gets the gene
   it exceeds by the largest log2 fold over the threshold. A cell positive for none (it
   passed the gate on Npy or Gad2) is `Other`.
7. **Subtype** (inhibitory cells only). k-means within each subclass on log1p counts of every
   gene except GFP, Slc17a7 and Gad2 (here Calb2, Npy, Pvalb, Sst, Vip). Fixed k: Pvalb 1,
   Sst 3, Vip 3, Other 2. Each sub-cluster is named `<subclass>_<gene>` after the gene
   whose z-scored mean most exceeds that of its sibling sub-clusters, counting only genes
   with a z-scored mean of at least 0.25 in the sub-cluster. A sub-cluster above its
   siblings on no gene is `<subclass>_only`. If no gene qualifies, it keeps the plain
   subclass name. Repeated names get a numeric suffix. With k = 1 the subtype is the
   subclass (`Pvalb`).
8. **Silhouette.** For each cell with a subtype, the silhouette of its subtype assignment
   within its subclass, on the same log1p features (step 7).

Per-subtype summaries, the silhouette-vs-k scan and heatmaps are in `inhibitory_gmm/subtypes/`.

## Columns

| Column | Type | Meaning |
|---|---|---|
| `cell_id` | int | HCR segmentation cell id (the label value in the segmentation mask). |
| `mouse_id` | str | Subject id. |
| `class` | str | `Inhibitory`, `Excitatory` or `Unassigned` (step 5). |
| `subclass` | str | `Pvalb`, `Sst`, `Vip` or `Other` (step 6); blank for non-inhibitory cells. |
| `subtype` | str | Subtype within the subclass (step 7), e.g. `Sst_Npy`, `Vip_Calb2`, `Vip_only`; blank for non-inhibitory cells. |
| `subtype_silhouette` | float | Silhouette of the subtype assignment, −1 to 1 (step 8). Near 1: close to its own subtype's centroid and far from its sibling subtypes. Near 0: on a boundary. Below 0: closer to a sibling subtype. Blank when the subclass has one subtype (Pvalb) or the cell is not inhibitory. |
| `gmm_inhibitory_positive` | bool | In the inhibitory set (step 3); identical to `class == "Inhibitory"`. |
| `slc17a7_positive` | bool | Above the Slc17a7 GMM threshold (step 4). Can be true for `Inhibitory` cells. |
| `GFP`, `Slc17a7`, `Gad2`, … | int | Raw mixed spot count per gene (step 1, before the noise floor), GFP/Slc17a7/Gad2 first. |

## Caveats

- Counts are **mixed** (not spectrally unmixed) spots, so some channel crosstalk remains.
- `subtype_silhouette` measures separation from the *other subtypes of the same subclass*
  only; it says nothing about the subclass call itself. The k values are fixed, not chosen
  by silhouette; `inhibitory_gmm/subtypes/silhouette_by_subclass.csv` shows how well each k
  separates each subclass.
- `Unassigned` cells are below both the inhibitory gate and the Slc17a7 threshold. Many are
  dim or small cells, not necessarily non-neuronal.
