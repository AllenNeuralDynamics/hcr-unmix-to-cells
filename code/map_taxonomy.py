# %%


from pathlib import Path
import pandas as pd
import os
import numpy as np
import anndata as ad

# MapMyCells mapping pipeline
from cell_type_mapper.test_utils.cache_wrapper import AbcCacheWrapper
from cell_type_mapper.cli.from_specified_markers import FromSpecifiedMarkersRunner

# Local imports
from xenium_analysis_tools.map_xenium.map_sections import ( 
    get_abc_paths, 
    get_v1_merfish_cells,
    get_nodes_to_drop,
    validate_input_adata,
    format_mapping_outputs
)

from xenium_analysis_tools.utils.io_utils import (
    load_config, 
)

# Environment setup (limit threads for numpy operations)
os.environ['NUMEXPR_NUM_THREADS'] = '1'
os.environ['MKL_NUM_THREADS'] = '1'
os.environ['OMP_NUM_THREADS'] = '1'

# %%
# Paths/typical configs
config = load_config('/root/capsule/code/params.json')
paths = config['paths']
mapping_config = config['mapping_config']

# Overwrite any existing mapping outputs
overwrite_input_adata = True
overwrite_mapping_results = True
overwrite_formatted_outputs = True


# %%

# Dataset folder
dataset_folder = Path(paths['data_root']) / 'HCR_767018_Oregano_251104'

# Cellxgene data
data_csv = '767018_Oregano_251104_inhibitory_clustered_cellxgene_lognorm.csv'
log_norm_data = True
cellxgene_path = dataset_folder / data_csv

# Output folders
folder_name = f'{cellxgene_path.stem}_runners_up'
dataset_output_folder = Path(paths['scratch_root']) / 'oregano_mapping' / folder_name
input_folder = dataset_output_folder / mapping_config['input_data_folder_name']
output_folder = dataset_output_folder / mapping_config['mapped_data_folder_name']
input_folder.mkdir(parents=True, exist_ok=True)
output_folder.mkdir(parents=True, exist_ok=True)

# Layer selection
drop_layers = ['VISp6a', 'VISp6b'] # None

# Output files
extended_results_path = output_folder / mapping_config['extended_results_name']
basic_results_path = output_folder / mapping_config['basic_results_name']


# %%

abc_atlas_path = Path(paths['abc_path'])
abc_cache = AbcCacheWrapper.from_local_cache(abc_atlas_path)
precomputed_stats_path, mouse_markers_path, gene_mapper_db_path = get_abc_paths(abc_cache)

# %%
input_adata_path = input_folder / mapping_config['input_h5ad_name']
if input_adata_path.exists() and not overwrite_input_adata:
    print(f"Input adata already exists at {input_adata_path}.")
else:
    print(f"Creating input adata from cellxgene data at {cellxgene_path}\nSaving to {input_adata_path}...")

    # Load cellxgene data
    cellxgene = pd.read_csv(cellxgene_path)

    ## Format adata
    # Cells
    cell_id_col = [col for col in cellxgene.columns if 'cell_' in col.lower()]
    cluster_cols = [col for col in cellxgene.columns if 'cluster' in col.lower()]
    print(f"Identified cluster columns in cellxgene data: {cluster_cols}")

    # Cell metadata
    obs = cellxgene[cell_id_col + cluster_cols]
    obs.rename(columns={col: f'hcr_{col}' for col in cluster_cols}, inplace=True)
    cellxgene.set_index(cell_id_col, inplace=True)
    cellxgene.drop(columns=cluster_cols, inplace=True)

    # Organized var
    if any(cellxgene.columns.str.match(r'^R\d+-')):
        var = pd.DataFrame([gn.split('-') for gn in cellxgene.columns], columns=['round','chan','gene_symbol'])
        var['round'] = var['round'].str.replace('R', '').astype(int)
        cellxgene.columns = var['gene_symbol'].values
    else:
        var = pd.DataFrame({'gene_symbol': cellxgene.columns})
    gene_metadata = abc_cache.get_metadata_dataframe(directory='WMB-10X', file_name='gene')
    var['gene_identifier'] = var['gene_symbol'].map(dict(zip(gene_metadata['gene_symbol'], gene_metadata['gene_identifier']))).fillna(var['gene_symbol'])

    # Replace cellxgene cols with ensembl ids
    cellxgene.columns = cellxgene.columns.map(dict(zip(var['gene_symbol'], var['gene_identifier'])))
    var.set_index('gene_identifier', inplace=True)

    # Order the same way
    cellxgene = cellxgene.loc[:,var.index]

    # X
    cxg_array = np.array(cellxgene)
    if log_norm_data:
        cxg_array = np.expm1(cxg_array)

    # Print out values of columns and indices to validate formatting
    print(f"Obs columns: \n{obs.columns.values}")
    print(f"\nVar columns: \n{var.columns.values}")
    print(f"\nVar indices (should be ensembl IDs): \n\t{var.index[:4].values}")
    if not np.any(var.index.isna()):
        print("All var indices have gene identifiers (no NaNs).")
    else:
        print("Warning: Some var indices are NaNs")

    # Create and save adata
    adata = ad.AnnData(X=cxg_array, obs=obs, var=var)
    print(f"Saving AnnData with shape: {adata.shape}")
    adata.write_h5ad(input_adata_path)


# %%

# Taxonomy filters - nodes to drop for mapping
nodes_to_drop=[]
# If specified any specific nodes to drop in config
drop_nodes_dict = mapping_config.get('drop_nodes_dict', None)
if drop_nodes_dict:
    for h_level in drop_nodes_dict:
        nodes_to_drop.extend([(h_level, cl) for cl in drop_nodes_dict[h_level]])
    print(f"Dropping {len(nodes_to_drop)} nodes based on drop_nodes_dict.")
# Filter to only V1 cells
filter_v1_types_config = mapping_config.get('filter_mapping_v1_types', None)
if filter_v1_types_config and filter_v1_types_config.get('enabled', False):
    h_level = filter_v1_types_config.get('h_level', 'subclass')
    min_cells = filter_v1_types_config.get('min_cells', 0)
    v1_types_df_name = filter_v1_types_config.get('saved_df_name', 'v1_merfish_cells.csv')
    if v1_types_df_name:
        v1_types_path = Path(paths['data_root']) / v1_types_df_name
    else:
        v1_types_path = None
    v1_merfish_cells = get_v1_merfish_cells(abc_cache, df_path=v1_types_path)
    if drop_layers:
        v1_merfish_cells = v1_merfish_cells.loc[~v1_merfish_cells['parcellation_substructure'].isin(drop_layers)]
    v1_nodes_to_drop = get_nodes_to_drop(v1_merfish_cells, abc_cache, h_level=h_level, min_cells=min_cells)
    print(f"Dropping {len(v1_nodes_to_drop)} {h_level} nodes not present in V1 MERFISH data with at least {min_cells if min_cells>0 else 1} cell(s).")
    nodes_to_drop.extend(v1_nodes_to_drop)

# %%
query_path = validate_input_adata(input_adata_path, input_adata_path.parent, mouse_markers_path, gene_mapper_db_path)

# %%

# ----- Mapper parameters -----
mapping_params = mapping_config.get('mapping_params', {})
mapping_params['nodes_to_drop'] = nodes_to_drop

# n_processors for mapper
num_workers = 4

# Type assignment parameters for mapper
type_assignment = {
    'normalization': 'raw',
    'bootstrap_iteration': 100,  
    'bootstrap_factor': 0.95,      
    'n_runners_up': 2,                                      
}

for key, val in type_assignment.items():
    print(f"{key}: {val}")

mapper_config = {
    'query_path': query_path,
    'extended_result_path': str(extended_results_path),
    'csv_result_path': str(basic_results_path),
    'flatten': False,
    'precomputed_stats': {'path': precomputed_stats_path},
    'query_markers': {'serialized_lookup': mouse_markers_path},
    'type_assignment': type_assignment,
    # 'gene_mapping': {'db_path': str(gene_mapper_db_path)}, # we know the genes
    'nodes_to_drop': mapping_params.get('nodes_to_drop', None),
    'verbose_stdout': True,
    'tmp_dir': '/tmp'
}
for key, val in mapper_config.items():
    if key != 'type_assignment':
        print(f"{key}: {val}")
    else:
        print(f"{key}:")
        for sub_key, sub_val in val.items():
            print(f"\t{sub_key}: {sub_val}")


# %%

if basic_results_path.exists() and extended_results_path.exists() and not overwrite_mapping_results:
    print(f"Mapping results already exist at {basic_results_path} and {extended_results_path}. Skipping mapping.")
else:
    runner = FromSpecifiedMarkersRunner(args=[], input_data=mapper_config)
    runner.run()

# %%

mapped_adata_path = output_folder / mapping_config['mapped_data_h5ad_name']
if mapped_adata_path.exists() and not overwrite_formatted_outputs:
    print(f"Formatted mapping output adata already exists at {mapped_adata_path}. Skipping formatting.")
else:
    format_mapping_outputs(extended_results_path, mapped_adata_path, mapping_params, h5ad_path=query_path)