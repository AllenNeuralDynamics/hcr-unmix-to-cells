#!/bin/bash
# Example usage script for taxonomy mapper

# This script demonstrates how to run the refactored taxonomy mapper
# with the same parameters as the original map_taxonomy.py script

python run_taxonomy_mapper.py \
    --config /root/capsule/code/params.json \
    --input-csv /root/capsule/data/HCR_767018_Oregano_251104/767018_Oregano_251104_inhibitory_clustered_cellxgene_lognorm.csv \
    --output-name 767018_Oregano_251104_inhibitory_clustered_cellxgene_lognorm_runners_up \
    --log-norm-data \
    --drop-layers VISp6a VISp6b \
    --bootstrap-iteration 100 \
    --bootstrap-factor 0.95 \
    --n-runners-up 2 \
    --num-workers 4 \
    --overwrite-input \
    --overwrite-mapping \
    --overwrite-formatted
