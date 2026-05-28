# Changelog

All notable changes to this capsule are documented in this file.

## 2026-05-28

- Added spot-mode input selection to [code/run_capsule.py](code/run_capsule.py) with `--spots {filtered|all_spots}`.
- Default behavior now runs taxonomy mapping on `filtered` spots.
- Added `all_spots` mode for taxonomy mapping runs.
- Implemented CSV resolution for both inhibitory and all-cells inputs using the new folder/file naming convention:
  - `inhibitory_cells_unmixed_filtered/unmixed_inhibitory_cells_filtered.csv`
  - `all_cells_unmixed_filtered/unmixed_all_cells_filtered.csv`
  - `inhibitory_cells_unmixed_all_spots/unmixed_inhibitory_cells_all_spots.csv`
  - `all_cells_unmixed_all_spots/unmixed_all_cells_all_spots.csv`
- Kept backward-compatible fallback support for legacy paths:
  - `inhibitory_cells_unmixed/unmixed_inhibitory_cells.csv`
  - `all_cells_unmixed/unmixed_all_cells.csv`
  - `unmixed_cell_by_gene_all_rounds.csv` (all-cells fallback)
- Updated output folder naming to include spot mode to avoid collisions:
  - `inhibitory_cells_filtered`, `all_cells_filtered`
  - `inhibitory_cells_all_spots`, `all_cells_all_spots`

## 2026-05-14

- Spectral ratio metrics and nn5 density computed on full pre-split spot population.
- Both unmixed and removed spots now carry matching QC columns.

## 2026-05-13

- Added neighbor ratio QC metrics.
- Added all-rounds parquet export.
- Added AnnData cell-by-gene output.
