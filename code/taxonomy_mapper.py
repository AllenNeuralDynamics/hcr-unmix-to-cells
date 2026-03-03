"""Core functions for taxonomy mapping pipeline."""

import os
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import anndata as ad
import numpy as np
import pandas as pd
from cell_type_mapper.cli.from_specified_markers import FromSpecifiedMarkersRunner
from cell_type_mapper.test_utils.cache_wrapper import AbcCacheWrapper

from xenium_analysis_tools.map_xenium.map_sections import (
    format_mapping_outputs,
    get_abc_paths,
    get_nodes_to_drop,
    get_v1_merfish_cells,
    validate_input_adata
)


def setup_environment():
    """Set up environment variables for threading control."""
    os.environ['NUMEXPR_NUM_THREADS'] = '1'
    os.environ['MKL_NUM_THREADS'] = '1'
    os.environ['OMP_NUM_THREADS'] = '1'


def load_abc_cache(abc_atlas_path: Path) -> AbcCacheWrapper:
    """Load ABC atlas cache and return paths to required resources."""
    return AbcCacheWrapper.from_local_cache(abc_atlas_path)


def create_input_adata(
    cellxgene_path: Path,
    output_path: Path,
    abc_cache: AbcCacheWrapper,
    log_norm_data: bool = True
) -> Path:
    """
    Create and save AnnData object from cellxgene CSV data.
    
    Args:
        cellxgene_path: Path to cellxgene CSV file
        output_path: Path to save the h5ad file
        abc_cache: ABC atlas cache wrapper
        log_norm_data: Whether the data is log-normalized (will apply expm1 if True)
    
    Returns:
        Path to saved h5ad file
    """
    print(f"Creating input adata from cellxgene data at {cellxgene_path}")
    print(f"Saving to {output_path}...")
    
    # Load cellxgene data
    cellxgene = pd.read_csv(cellxgene_path)
    
    # Identify cell ID and cluster columns
    cell_id_col = [col for col in cellxgene.columns if 'cell_' in col.lower()]
    cluster_cols = [col for col in cellxgene.columns if 'cluster' in col.lower()]
    print(f"Identified cluster columns in cellxgene data: {cluster_cols}")
    
    # Create cell metadata
    obs = cellxgene[cell_id_col + cluster_cols]
    obs.rename(columns={col: f'hcr_{col}' for col in cluster_cols}, inplace=True)
    cellxgene.set_index(cell_id_col, inplace=True)
    cellxgene.drop(columns=cluster_cols, inplace=True)
    
    # Process gene/variable metadata
    if any(cellxgene.columns.str.match(r'^R\d+-')):
        # Parse round-channel-gene format
        var = pd.DataFrame(
            [gn.split('-') for gn in cellxgene.columns],
            columns=['round', 'chan', 'gene_symbol']
        )
        var['round'] = var['round'].str.replace('R', '').astype(int)
        cellxgene.columns = var['gene_symbol'].values
    else:
        var = pd.DataFrame({'gene_symbol': cellxgene.columns})
    
    # Map gene symbols to identifiers
    gene_metadata = abc_cache.get_metadata_dataframe(directory='WMB-10X', file_name='gene')
    var['gene_identifier'] = var['gene_symbol'].map(
        dict(zip(gene_metadata['gene_symbol'], gene_metadata['gene_identifier']))
    ).fillna(var['gene_symbol'])
    
    # Replace column names with ensembl IDs
    cellxgene.columns = cellxgene.columns.map(
        dict(zip(var['gene_symbol'], var['gene_identifier']))
    )
    var.set_index('gene_identifier', inplace=True)
    
    # Ensure same ordering
    cellxgene = cellxgene.loc[:, var.index]
    
    # Convert to array and handle log-normalization
    cxg_array = np.array(cellxgene)
    if log_norm_data:
        cxg_array = np.expm1(cxg_array)
    
    # Validation output
    print(f"Obs columns: \n{obs.columns.values}")
    print(f"\nVar columns: \n{var.columns.values}")
    print(f"\nVar indices (should be ensembl IDs): \n\t{var.index[:4].values}")
    if not np.any(var.index.isna()):
        print("All var indices have gene identifiers (no NaNs).")
    else:
        print("Warning: Some var indices are NaNs")
    
    # Create and save AnnData
    adata = ad.AnnData(X=cxg_array, obs=obs, var=var)
    print(f"Saving AnnData with shape: {adata.shape}")
    adata.write_h5ad(output_path)
    
    return output_path


def prepare_taxonomy_filters(
    abc_cache: AbcCacheWrapper,
    data_root: Path,
    drop_nodes_dict: Optional[Dict[str, List[str]]] = None,
    filter_v1_config: Optional[Dict] = None,
    drop_layers: Optional[List[str]] = None
) -> List[Tuple[str, str]]:
    """
    Prepare taxonomy node filters for mapping.
    
    Args:
        abc_cache: ABC atlas cache wrapper
        data_root: Root data directory
        drop_nodes_dict: Dictionary mapping hierarchy levels to node names to drop
        filter_v1_config: Configuration for V1 MERFISH filtering
        drop_layers: List of layer names to exclude
    
    Returns:
        List of (hierarchy_level, node_name) tuples to drop
    """
    nodes_to_drop = []
    
    # Add nodes from drop_nodes_dict
    if drop_nodes_dict:
        for h_level in drop_nodes_dict:
            nodes_to_drop.extend([(h_level, cl) for cl in drop_nodes_dict[h_level]])
        print(f"Dropping {len(nodes_to_drop)} nodes based on drop_nodes_dict.")
    
    # Filter to only V1 cell types
    if filter_v1_config and filter_v1_config.get('enabled', False):
        h_level = filter_v1_config.get('h_level', 'subclass')
        min_cells = filter_v1_config.get('min_cells', 0)
        v1_types_df_name = filter_v1_config.get('saved_df_name', 'v1_merfish_cells.csv')
        
        v1_types_path = data_root / v1_types_df_name if v1_types_df_name else None
        v1_merfish_cells = get_v1_merfish_cells(abc_cache, df_path=v1_types_path)
        
        if drop_layers:
            v1_merfish_cells = v1_merfish_cells.loc[
                ~v1_merfish_cells['parcellation_substructure'].isin(drop_layers)
            ]
        
        v1_nodes_to_drop = get_nodes_to_drop(
            v1_merfish_cells, abc_cache, h_level=h_level, min_cells=min_cells
        )
        print(f"Dropping {len(v1_nodes_to_drop)} {h_level} nodes not present in V1 "
              f"MERFISH data with at least {min_cells if min_cells > 0 else 1} cell(s).")
        nodes_to_drop.extend(v1_nodes_to_drop)
    
    return nodes_to_drop


def run_mapping(
    query_path: Path,
    extended_result_path: Path,
    basic_result_path: Path,
    precomputed_stats_path: Path,
    mouse_markers_path: Path,
    mapping_params: Dict,
    nodes_to_drop: Optional[List[Tuple[str, str]]] = None
) -> None:
    """
    Run the cell type mapping pipeline.
    
    Args:
        query_path: Path to validated query h5ad file
        extended_result_path: Path to save extended results JSON
        basic_result_path: Path to save basic results CSV
        precomputed_stats_path: Path to precomputed statistics
        mouse_markers_path: Path to mouse marker genes
        mapping_params: Dictionary of mapping parameters
        nodes_to_drop: List of taxonomy nodes to exclude
    """
    # Build type assignment configuration
    type_assignment = {
        'normalization': mapping_params.get('normalization', 'raw'),
        'bootstrap_iteration': mapping_params.get('bootstrap_iteration', 100),
        'bootstrap_factor': mapping_params.get('bootstrap_factor', 0.95),
        'n_runners_up': mapping_params.get('n_runners_up', 2),
    }
    
    print("Type assignment parameters:")
    for key, val in type_assignment.items():
        print(f"  {key}: {val}")
    
    # Build mapper configuration
    mapper_config = {
        'query_path': query_path,
        'extended_result_path': str(extended_result_path),
        'csv_result_path': str(basic_result_path),
        'flatten': False,
        'precomputed_stats': {'path': precomputed_stats_path},
        'query_markers': {'serialized_lookup': mouse_markers_path},
        'type_assignment': type_assignment,
        'nodes_to_drop': nodes_to_drop,
        'verbose_stdout': True,
        'tmp_dir': '/tmp'
    }
    
    print("\nMapper configuration:")
    for key, val in mapper_config.items():
        if key != 'type_assignment':
            print(f"  {key}: {val}")
    
    # Run the mapper
    runner = FromSpecifiedMarkersRunner(args=[], input_data=mapper_config)
    runner.run()


def format_and_save_results(
    extended_results_path: Path,
    mapped_adata_path: Path,
    query_path: Path,
    mapping_params: Dict
) -> None:
    """
    Format mapping results and save as AnnData.
    
    Args:
        extended_results_path: Path to extended results JSON
        mapped_adata_path: Path to save formatted results h5ad
        query_path: Path to original query h5ad file
        mapping_params: Dictionary of mapping parameters
    """
    print(f"Formatting mapping outputs and saving to {mapped_adata_path}...")
    format_mapping_outputs(
        extended_results_path,
        mapped_adata_path,
        mapping_params,
        h5ad_path=query_path
    )
