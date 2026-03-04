"""
Example script: Process a single mouse round through the spot unmixing pipeline.
"""

from pathlib import Path
from spot_analysis.spot_pipeline import (
    SpotPipelineConfig,
    run_spot_pipeline_v2,
)
from aind_hcr_data_loader.hcr_dataset import HCRDataset

# ============================================================================
# 1. CONFIGURE
# ============================================================================

MOUSE_ID = "mouse123"
ROUND_KEY = "R3"
DATA_DIR = Path("/root/capsule/data")

# ============================================================================
# 2. LOAD DATASET
# ============================================================================

ds = HCRDataset(DATA_DIR / MOUSE_ID)
datasets = {MOUSE_ID: ds}

# ============================================================================
# 3. SET UP PIPELINE CONFIG
# ============================================================================

config = SpotPipelineConfig(
    expt_key="test_run",
    min_distances=[1],
    dist_cutoff=1.0,
)

# ============================================================================
# 4. RUN
# ============================================================================

results = run_spot_pipeline_v2(datasets, MOUSE_ID, ROUND_KEY, config)

# ============================================================================
# 5. INSPECT OUTPUTS
# ============================================================================

print(f"\nOutput folder: {results['output_folder']}")
print(f"Total spots loaded:   {len(results['pipeline_data'].spots_df)}")
print(f"Spots after cleaning: {len(results['spots_df_cleaned'])}")
print(f"Spots after unmixing: {len(results['spots_df_full'])}")
print(f"\nMixed total counts:   {results['mixed_results']['spot_count'].sum()}")
print(f"Unmixed total counts: {results['unmixed_results']['spot_count'].sum()}")

print("\nRatios matrix:")
print(results['ratios_matrix'])

print("\nReassignment matrix:")
print(results['reassignment_matrix'])

print("\nSpot fate matrix:")
print(results['spot_fate_matrix'])