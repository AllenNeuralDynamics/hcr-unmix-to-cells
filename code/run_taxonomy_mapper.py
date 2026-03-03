#!/usr/bin/env python
"""
Command-line interface for running taxonomy mapping pipeline.

This script provides a CLI for mapping cell types using the ABC atlas
and MapMyCells pipeline.
"""

import argparse
import sys
from pathlib import Path

from config import TaxonomyMapperConfig
from taxonomy_mapper import (
    create_input_adata,
    format_and_save_results,
    load_abc_cache,
    prepare_taxonomy_filters,
    run_mapping,
    setup_environment
)
from xenium_analysis_tools.map_xenium.map_sections import (
    get_abc_paths,
    validate_input_adata
)


def parse_args():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description='Run taxonomy mapping pipeline for spatial transcriptomics data'
    )
    
    # Required arguments
    parser.add_argument(
        '--config',
        type=str,
        default='/root/capsule/code/params.json',
        help='Path to configuration JSON file'
    )
    parser.add_argument(
        '--input-csv',
        type=str,
        required=True,
        help='Path to input CSV file (cellxgene data)'
    )
    
    # Optional arguments
    parser.add_argument(
        '--output-name',
        type=str,
        help='Custom output folder name (default: derived from input CSV name)'
    )
    parser.add_argument(
        '--log-norm-data',
        action='store_true',
        default=True,
        help='Input data is log-normalized (will apply expm1)'
    )
    parser.add_argument(
        '--no-log-norm-data',
        dest='log_norm_data',
        action='store_false',
        help='Input data is NOT log-normalized'
    )
    parser.add_argument(
        '--drop-layers',
        type=str,
        nargs='+',
        help='Layer names to exclude (e.g., VISp6a VISp6b)'
    )
    
    # Overwrite flags
    parser.add_argument(
        '--overwrite-input',
        action='store_true',
        help='Overwrite existing input h5ad file'
    )
    parser.add_argument(
        '--overwrite-mapping',
        action='store_true',
        help='Overwrite existing mapping results'
    )
    parser.add_argument(
        '--overwrite-formatted',
        action='store_true',
        help='Overwrite existing formatted output'
    )
    parser.add_argument(
        '--overwrite-all',
        action='store_true',
        help='Overwrite all existing outputs'
    )
    
    # Mapping parameters (overrides)
    parser.add_argument(
        '--bootstrap-iteration',
        type=int,
        help='Number of bootstrap iterations'
    )
    parser.add_argument(
        '--bootstrap-factor',
        type=float,
        help='Bootstrap sampling factor'
    )
    parser.add_argument(
        '--n-runners-up',
        type=int,
        help='Number of runner-up cell types to report'
    )
    parser.add_argument(
        '--num-workers',
        type=int,
        help='Number of parallel workers for mapping'
    )
    
    return parser.parse_args()


def main():
    """Main execution function."""
    args = parse_args()
    
    # Set up environment
    setup_environment()
    
    # Load configuration
    print(f"Loading configuration from {args.config}")
    config = TaxonomyMapperConfig.from_json(args.config)
    
    # Parse input CSV path
    cellxgene_path = Path(args.input_csv)
    if not cellxgene_path.exists():
        print(f"Error: Input file not found: {cellxgene_path}")
        sys.exit(1)
    
    # Update config with command line arguments
    config.dataset_folder = cellxgene_path.parent
    config.data_csv = cellxgene_path.name
    config.log_norm_data = args.log_norm_data
    config.drop_layers = args.drop_layers
    
    if args.output_name:
        config.output_folder_name = args.output_name
    
    # Set overwrite flags
    if args.overwrite_all:
        config.overwrite_input_adata = True
        config.overwrite_mapping_results = True
        config.overwrite_formatted_outputs = True
    else:
        config.overwrite_input_adata = args.overwrite_input or config.overwrite_input_adata
        config.overwrite_mapping_results = args.overwrite_mapping or config.overwrite_mapping_results
        config.overwrite_formatted_outputs = args.overwrite_formatted or config.overwrite_formatted_outputs
    
    # Override mapping parameters if provided
    if args.bootstrap_iteration:
        config.mapping_params.bootstrap_iteration = args.bootstrap_iteration
    if args.bootstrap_factor:
        config.mapping_params.bootstrap_factor = args.bootstrap_factor
    if args.n_runners_up is not None:
        config.mapping_params.n_runners_up = args.n_runners_up
    if args.num_workers:
        config.mapping_params.num_workers = args.num_workers
    
    # Validate input file exists
    if not cellxgene_path.exists():
        print(f"Error: Input file not found: {cellxgene_path}")
        sys.exit(1)
    
    input_folder, output_folder, dataset_output_folder = config.get_output_paths()
    input_folder.mkdir(parents=True, exist_ok=True)
    output_folder.mkdir(parents=True, exist_ok=True)
    
    print(f"\nInput CSV: {cellxgene_path}")
    print(f"Output folder: {dataset_output_folder}")
    
    # Define output file paths
    input_adata_path = input_folder / config.input_h5ad_name
    extended_results_path = output_folder / config.extended_results_name
    basic_results_path = output_folder / config.basic_results_name
    mapped_adata_path = output_folder / config.mapped_data_h5ad_name
    
    # Load ABC atlas cache
    print(f"\nLoading ABC atlas from {config.paths.abc_path}")
    abc_cache = load_abc_cache(config.paths.abc_path)
    precomputed_stats_path, mouse_markers_path, gene_mapper_db_path = get_abc_paths(abc_cache)
    
    # Step 1: Create input AnnData
    print("\n" + "="*80)
    print("STEP 1: Creating input AnnData")
    print("="*80)
    
    if input_adata_path.exists() and not config.overwrite_input_adata:
        print(f"Input adata already exists at {input_adata_path}.")
    else:
        create_input_adata(
            cellxgene_path,
            input_adata_path,
            abc_cache,
            log_norm_data=config.log_norm_data
        )
    
    # Step 2: Prepare taxonomy filters
    print("\n" + "="*80)
    print("STEP 2: Preparing taxonomy filters")
    print("="*80)
    
    filter_v1_config = {
        'enabled': config.filter_config.enabled,
        'h_level': config.filter_config.h_level,
        'min_cells': config.filter_config.min_cells,
        'saved_df_name': config.filter_config.saved_df_name
    }
    
    nodes_to_drop = prepare_taxonomy_filters(
        abc_cache,
        config.paths.data_root,
        drop_nodes_dict=config.drop_nodes_dict,
        filter_v1_config=filter_v1_config,
        drop_layers=config.drop_layers
    )
    
    # Step 3: Validate input
    print("\n" + "="*80)
    print("STEP 3: Validating input AnnData")
    print("="*80)
    
    query_path = validate_input_adata(
        input_adata_path,
        input_adata_path.parent,
        mouse_markers_path,
        gene_mapper_db_path
    )
    
    # Step 4: Run mapping
    print("\n" + "="*80)
    print("STEP 4: Running cell type mapping")
    print("="*80)
    
    if (basic_results_path.exists() and extended_results_path.exists() 
        and not config.overwrite_mapping_results):
        print(f"Mapping results already exist. Skipping mapping.")
    else:
        mapping_params_dict = {
            'normalization': config.mapping_params.normalization,
            'bootstrap_iteration': config.mapping_params.bootstrap_iteration,
            'bootstrap_factor': config.mapping_params.bootstrap_factor,
            'n_runners_up': config.mapping_params.n_runners_up,
        }
        
        run_mapping(
            query_path,
            extended_results_path,
            basic_results_path,
            precomputed_stats_path,
            mouse_markers_path,
            mapping_params_dict,
            nodes_to_drop=nodes_to_drop
        )
    
    # Step 5: Format and save results
    print("\n" + "="*80)
    print("STEP 5: Formatting and saving results")
    print("="*80)
    
    if mapped_adata_path.exists() and not config.overwrite_formatted_outputs:
        print(f"Formatted mapping output already exists at {mapped_adata_path}.")
    else:
        format_and_save_results(
            extended_results_path,
            mapped_adata_path,
            query_path,
            {'nodes_to_drop': nodes_to_drop}
        )
    
    print("\n" + "="*80)
    print("PIPELINE COMPLETE")
    print("="*80)
    print(f"Results saved to: {output_folder}")
    print(f"  - Basic results: {basic_results_path}")
    print(f"  - Extended results: {extended_results_path}")
    print(f"  - Mapped AnnData: {mapped_adata_path}")


if __name__ == '__main__':
    main()
