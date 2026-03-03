"""
Example: Programmatic usage of the taxonomy mapper

This script demonstrates how to use the refactored taxonomy mapper
functions programmatically in your own Python code.
"""

from pathlib import Path
from config import TaxonomyMapperConfig
from taxonomy_mapper import (
    setup_environment,
    load_abc_cache,
    create_input_adata,
    prepare_taxonomy_filters,
    run_mapping,
    format_and_save_results
)
from xenium_analysis_tools.map_xenium.map_sections import (
    get_abc_paths,
    validate_input_adata
)


def run_custom_mapping():
    """Example of running the mapping pipeline programmatically."""
    
    # Set up environment
    setup_environment()
    
    # Load configuration
    config = TaxonomyMapperConfig.from_json('/root/capsule/code/params.json')
    
    # Customize configuration
    config.dataset_folder = Path('/root/capsule/data/HCR_767018_Oregano_251104')
    config.data_csv = '767018_Oregano_251104_inhibitory_clustered_cellxgene_lognorm.csv'
    config.log_norm_data = True
    config.drop_layers = ['VISp6a', 'VISp6b']
    config.output_folder_name = 'custom_oregano_mapping'
    
    # Customize mapping parameters
    config.mapping_params.bootstrap_iteration = 150
    config.mapping_params.n_runners_up = 3
    
    # Set up paths
    cellxgene_path = config.dataset_folder / config.data_csv
    input_folder, output_folder, _ = config.get_output_paths()
    input_folder.mkdir(parents=True, exist_ok=True)
    output_folder.mkdir(parents=True, exist_ok=True)
    
    input_adata_path = input_folder / config.input_h5ad_name
    extended_results_path = output_folder / config.extended_results_name
    basic_results_path = output_folder / config.basic_results_name
    mapped_adata_path = output_folder / config.mapped_data_h5ad_name
    
    # Load ABC cache
    abc_cache = load_abc_cache(config.paths.abc_path)
    precomputed_stats_path, mouse_markers_path, gene_mapper_db_path = get_abc_paths(abc_cache)
    
    # Step 1: Create input AnnData
    if not input_adata_path.exists():
        create_input_adata(
            cellxgene_path,
            input_adata_path,
            abc_cache,
            log_norm_data=config.log_norm_data
        )
    
    # Step 2: Prepare filters
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
    query_path = validate_input_adata(
        input_adata_path,
        input_adata_path.parent,
        mouse_markers_path,
        gene_mapper_db_path
    )
    
    # Step 4: Run mapping
    if not (basic_results_path.exists() and extended_results_path.exists()):
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
    
    # Step 5: Format results
    if not mapped_adata_path.exists():
        format_and_save_results(
            extended_results_path,
            mapped_adata_path,
            query_path,
            {'nodes_to_drop': nodes_to_drop}
        )
    
    print(f"\nMapping complete! Results saved to: {output_folder}")
    return mapped_adata_path


if __name__ == '__main__':
    result_path = run_custom_mapping()
    print(f"Final output: {result_path}")
