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


# ----
# Plotting functions
# -----

from pathlib import Path
import pandas as pd
import matplotlib
import scanpy as sc
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from matplotlib.patches import Patch
import plotly.graph_objects as go
import matplotlib.colors as mcolors
matplotlib.rcParams['svg.fonttype'] = 'none'

def map_to_broad_subclass_name(subclass_name, patterns):
    if pd.isna(subclass_name):
        return 'Other'
    
    subclass_str = str(subclass_name).strip()
    
    for broad_name, pattern_list in patterns.items():
        for pattern in pattern_list:
            # Exact match (case insensitive)
            if subclass_str.lower() == pattern.lower():
                return broad_name
            
            # Contains match (case insensitive)
            if pattern.lower() in subclass_str.lower():
                return broad_name
    return 'Other'

def add_broad_types(adata):
    broad_subclass_name_patterns = {
        'L2/3 IT': ['L2/3 IT', 'L2/3IT', 'L23 IT', 'L2-3 IT'],
        'L4/5 IT': ['L4/5 IT', 'L4/5IT', 'L4 IT', 'L45 IT', 'L4-5 IT'],
        'L5 IT': ['L5 IT', 'L5IT'],
        'L5 ET': ['L5 ET', 'L5ET'],
        'L6 IT': ['L6 IT', 'L6IT'],
        'L6 CT': ['L6 CT', 'L6CT'],
        'Lamp5': ['Lamp5', 'LAMP5'],
        'Vip': ['Vip', 'VIP'],
        'Sncg': ['Sncg', 'SNCG'],
        'Sst': ['Sst', 'SST'],
        'Pvalb': ['Pvalb', 'PVALB', 'PV'],
        'NN': ['NN', 'Neuronal Non-Neuronal', 'Neuronal-Non-Neuronal'],
        'L6b': ['L6b', 'L6-b', 'L6 b'],
        'L2 IT': ['L2 IT', 'L2IT', 'L2-IT'],
        'Ob-in': ['Ob-in', 'OB-in', 'OB in', 'Ob in'],
        'STR': ['STR Prox1 Lhx6', 'STR Prox1 LHX6', 'STR Prox1-Lhx6', 'STR Prox1-LHX6'],
        'CB Gran': ['CB Gran', 'CB Granule', 'CB Granular', 'CBGran'],
    }

    adata.obs['broad_class_name'] = 'Other'
    adata.obs.loc[adata.obs['class_name'].str.contains('Glut', case=False, na=False), 'broad_class_name'] = 'Glut'
    adata.obs.loc[adata.obs['class_name'].str.contains('GABA', case=False, na=False), 'broad_class_name'] = 'GABA'
    
    def parse_numeric_class(class_name):
        try:
            if pd.isna(class_name):
                return False
            parts = str(class_name).strip().split(' ')
            if len(parts) > 0 and parts[0].isdigit():
                return int(parts[0]) >= 30
            return False
        except (ValueError, IndexError):
            return False
    
    adata.obs.loc[adata.obs['class_name'].apply(parse_numeric_class), 'broad_class_name'] = 'NN'

    adata.obs['broad_subclass_name'] = adata.obs['subclass_name'].apply(
        lambda x: map_to_broad_subclass_name(x, broad_subclass_name_patterns)
    )
    adata.obs['broad_class_name'] = adata.obs['broad_class_name'].astype('category')
    adata.obs['broad_subclass_name'] = adata.obs['broad_subclass_name'].astype('category')

    return adata

def get_types_breakdown(adata, col_name, sub_col_name=None, col_val=None, col_width=[30, 30, 10], print_output=False):
    if col_name in ['broad_class_name','broad_subclass_name']:
        col_width[0] = 10

    if sub_col_name is None:
        overall_counts = adata.obs[col_name].value_counts()
        if print_output:
            print(f"{col_name} counts:")
            print("-" * 70)
            for class_name, count in overall_counts.items():
                print(f"{class_name}: {count}")
            print("-" * 70)
        print(f"Total: {len(adata.obs)}")
    else:
        if col_val is not None:
            adata = adata[adata.obs[col_name] == col_val]
        sub_counts = adata.obs.groupby([col_name, sub_col_name]).size().sort_values(ascending=False)
        filtered_counts = sub_counts.loc[sub_counts > 1].sort_index()

        if print_output:
            print(f"{f'\n{col_name}':<{col_width[0]}} {sub_col_name:<{col_width[1]}} {'count':<{col_width[2]}}")
            print("-" * 70)
            for (class_name, subclass_name), count in filtered_counts.items():
                print(f"{class_name:<{col_width[0]}} {subclass_name:<{col_width[1]}} {count:<{col_width[2]}}")
            print("-" * 70)
            print(f"Total: {filtered_counts.sum()}")
            print()

        return pd.DataFrame(filtered_counts, columns=['count']).reset_index()

def get_shared_colormap(data):
    joined_colormap = {}
    gaba_palette = 'plasma'
    glut_palette = 'cividis'
    nn_palette = 'copper'

    def _add_color_for_subclasses_supertypes(subclass_data, broad_subclass_names, cmap):
        for i, sbcl in enumerate(broad_subclass_names):
            broad_subclass_name_col = cmap[i]
            joined_colormap[sbcl] = broad_subclass_name_col
            # Get subclasses
            subclass_mask = subclass_data['broad_subclass_name'] == sbcl
            subclasses = sorted(subclass_data.loc[subclass_mask, 'subclass_name'].unique())

            # All subclasses of the same broad_subclass_name assigned the same color
            for subclass in subclasses:
                joined_colormap[subclass] = broad_subclass_name_col

            # Supertypes as assigned progressively desaturated versions of the broad_subclass_name color
            supertypes = sorted(subclass_data.loc[subclass_mask, 'supertype_name'].unique())
            saturation_factors = np.linspace(1.0, 0.4, num=len(supertypes))
            for supertype, saturation_factor in zip(supertypes, saturation_factors):
                joined_colormap[supertype] = sns.desaturate(broad_subclass_name_col, saturation_factor)

    # GABA types
    joined_colormap['GABA'] = 'tab:red'
    gaba_data = data.obs[data.obs['broad_class_name'] == 'GABA']
    gaba_subclasses = sorted(gaba_data['broad_subclass_name'].unique())
    n_gaba_subclasses = len(gaba_subclasses)
    # GABA subclasses/supertypes palette
    cmap = sns.color_palette(gaba_palette, n_colors=n_gaba_subclasses)
    _add_color_for_subclasses_supertypes(gaba_data, gaba_subclasses, cmap)

    # Glut types
    joined_colormap['Glut'] = 'tab:blue'
    glut_data = data.obs[data.obs['broad_class_name'] == 'Glut']
    glut_subclasses = sorted(glut_data['broad_subclass_name'].unique())
    n_glut_subclasses = len(glut_subclasses)
    # Glut subclasses/supertypes palette
    cmap = sns.color_palette(glut_palette, n_colors=n_glut_subclasses)
    _add_color_for_subclasses_supertypes(glut_data, glut_subclasses, cmap)

    # NN types
    joined_colormap['NN'] = 'tab:brown'
    nn_data = data.obs[data.obs['broad_class_name'] == 'NN']
    nn_subclasses = sorted(nn_data['broad_subclass_name'].unique())
    n_nn_subclasses = len(nn_subclasses)
    # NN subclasses/supertypes palette
    cmap = sns.color_palette(nn_palette, n_colors=n_nn_subclasses)
    _add_color_for_subclasses_supertypes(nn_data, nn_subclasses, cmap)

    return joined_colormap

def add_colormap_adata(adata, colormap):
    broad_class_name_colors = []
    for tp in sorted(adata.obs['broad_class_name'].unique()):
        color = colormap.get(tp, 'gray')
        broad_class_name_colors.append(color)

    broad_subclass_name_colors = []
    for tp in sorted(adata.obs['broad_subclass_name'].unique()):
        color = colormap.get(tp, 'gray')
        broad_subclass_name_colors.append(color)

    subclass_colors = []
    for tp in sorted(adata.obs['subclass_name'].unique()):
        color = colormap.get(tp, 'gray')
        subclass_colors.append(color)

    supertype_colors = []
    for tp in sorted(adata.obs['supertype_name'].unique()):
        color = colormap.get(tp, 'gray')
        supertype_colors.append(color)

    # Store as numpy arrays
    adata.uns['broad_class_name_colors'] = np.array(broad_class_name_colors, dtype=object)
    adata.uns['broad_subclass_name_colors'] = np.array(broad_subclass_name_colors, dtype=object)
    adata.uns['supertype_name_colors'] = np.array(supertype_colors, dtype=object)
    adata.uns['subclass_name_colors'] = np.array(subclass_colors, dtype=object)

    return adata

def create_sankey_diagram(data, columns, colormap=None, title=None, 
                         height=900, width=1200, sort_columns=None,
                         default_node_color='#87CEEB'):
    
    def to_hex(color):
        try: return mcolors.to_hex(color)
        except: return default_node_color

    if isinstance(data, sc.AnnData):
        obs_data = data.obs[columns].astype(str).copy()
    else:
        obs_data = data[columns].astype(str).copy()

    nodes, node_colors, x_pos, y_pos = [], [], [], []
    node_to_idx = {}
    
    # Spacing
    n_cols = len(columns)
    x_coords = np.linspace(0.01, 0.99, n_cols)

    for i, col in enumerate(columns):
        unique_vals = obs_data[col].unique()
        
        if sort_columns and col in sort_columns:
            sort_type = sort_columns[col]
            if sort_type == 'reverse':
                unique_vals = sorted(unique_vals, reverse=True)
            elif sort_type == 'normal':
                unique_vals = sorted(unique_vals)
            elif callable(sort_type):
                unique_vals = sorted(unique_vals, key=sort_type)
        else:
            unique_vals = sorted(unique_vals)
            
        n_nodes = len(unique_vals)
        y_coords = np.linspace(0.01, 0.99, n_nodes)

        for j, val in enumerate(unique_vals):
            node_key = (col, val)
            node_to_idx[node_key] = len(nodes)
            nodes.append(val)
            x_pos.append(x_coords[i])
            y_pos.append(y_coords[j])
            
            raw_color = colormap.get(val, default_node_color) if colormap else default_node_color
            node_colors.append(to_hex(raw_color))

    sources, targets, values = [], [], []
    for i in range(len(columns) - 1):
        s_col, t_col = columns[i], columns[i+1]
        counts = obs_data.groupby([s_col, t_col]).size().reset_index(name='count')
        for _, row in counts.iterrows():
            sources.append(node_to_idx[(s_col, row[s_col])])
            targets.append(node_to_idx[(t_col, row[t_col])])
            values.append(row['count'])

    fig = go.Figure(data=[go.Sankey(
        arrangement="fixed", 
        node=dict(
            pad=15, thickness=20,
            line=dict(color="black", width=0.5),
            label=nodes,
            color=node_colors,
            x=x_pos,
            y=y_pos,
            hovertemplate='%{label}<br>Total: %{value}<extra></extra>'
        ),
        link=dict(
            source=sources, target=targets, value=values,
            color="rgba(211, 211, 211, 0.4)",
            hovercolor="rgba(80, 80, 80, 0.8)",
            hovertemplate='<b>%{source.label}</b> → <b>%{target.label}</b><br>Count: %{value}<extra></extra>'
        )
    )])

    fig.update_layout(title_text=title or "River Plot", 
                height=height, width=width, 
                font_size=12)
    return fig

def plot_mapping_comparison(mapping_1, 
                            mapping_2, 
                            level='subclass_name',
                            mapping_1_name='mapping_1',
                            mapping_2_name='mapping_2',
                            ax=None,
                            colormap='viridis',
                            figsize=(10, 10)):
    # Load if necessary
    if isinstance(mapping_1, str) or isinstance(mapping_1, Path):
        mapping_1 = sc.read_h5ad(mapping_1)
        mapping_1 = mapping_1
    if isinstance(mapping_2, str) or isinstance(mapping_2, Path):
        mapping_2 = sc.read_h5ad(mapping_2)
        mapping_2 = mapping_2
        # Format columns                       
    mapping_2.obs.columns = mapping_2.obs.columns.str.replace('CDM_', '', regex=False)
    mapping_2 = add_broad_types(mapping_2)

    # Merge the two mappings on cell_id and format
    mapping_cols = [col for col in mapping_1.obs.columns if col.endswith('_name')]
    mapping_cols = ['cell_id'] + mapping_cols
    if isinstance(mapping_1, sc.AnnData):
        mapping_1 = mapping_1.obs
    if isinstance(mapping_2, sc.AnnData):
        mapping_2 = mapping_2.obs

    # Merge
    plot_data = mapping_1[mapping_cols].merge(
            mapping_2[mapping_cols].copy(),
            on='cell_id',
            suffixes=(f'_{mapping_1_name}', f'_{mapping_2_name}')
        )

    # Calculate proportion of cells with consistent mapping at the specified level
    plot_data[plot_data.select_dtypes(include=['category']).columns] = plot_data[plot_data.select_dtypes(include=['category']).columns].astype(str)

    # Plot
    plot_groups = len(np.unique(plot_data[[f'{level}_{mapping_1_name}', f'{level}_{mapping_2_name}']]))
    cross_tabs = pd.crosstab(plot_data[f'{level}_{mapping_1_name}'], plot_data[f'{level}_{mapping_2_name}'], normalize='index')
    if ax is None:
        fig, ax = plt.subplots(1,1,figsize=figsize)
    sns.heatmap(cross_tabs, 
                cmap=colormap, 
                ax=ax, 
                square=True, 
                annot=True if plot_groups < 40 else False,
                fmt='.1f',
                linewidths=0.1 if plot_groups < 40 else 0,
                linecolor='gray',
                annot_kws={'size': 6},
                cbar_kws={'label': 'Proportion of cells', 'shrink': 0.5})
    ax.set_title(level)
    ax.set_xlabel(f'{mapping_1_name}')
    ax.set_ylabel(f'{mapping_2_name}')

    return plot_data, fig

def plot_stacked_categories(data, main_cat, second_cat, colormap, ax, kwargs={}, figsize=(12, 8)):
    counts = get_types_breakdown(data, main_cat, second_cat, print_output=False)
    counts['proportion'] = counts.groupby(main_cat)['count'].transform(lambda x: x / x.sum())
    counts = counts.sort_values(by=['proportion'], ascending=False)
    pivot_data = counts.pivot(index=main_cat, columns=second_cat, values='proportion').fillna(0)

    # Plot stacked bars
    pivot_data.plot(kind='bar', stacked=True, ax=ax, 
                    color=[colormap.get(col, 'gray') for col in pivot_data.columns],
                    figsize=figsize)

    for i, (idx, row) in enumerate(pivot_data.iterrows()):
        bottom = 0
        for j, (col, value) in enumerate(row.items()):
            if value > 0.05:
                actual_count = counts[(counts[main_cat] == idx) & 
                                    (counts[second_cat] == col)]['count'].values
                if len(actual_count) > 0:
                    count_text = f'{col}\n{actual_count[0]} ({value*100:.1f}%)'
                    ax.text(i, bottom + value/2, count_text, 
                        ha='center', va='center', fontweight=kwargs.get('font_weight', 'bold'),
                        fontsize=kwargs.get('font_size', 8), color=kwargs.get('font_color', 'white'))
            bottom += value
            
    ax.set_title(f'Proportion of {second_cat.replace("_", " ").title()} within each {main_cat.replace("_", " ").title()}')
    ax.set_xlabel(main_cat.replace('_', ' ').title())
    ax.set_ylabel('Proportion')
    ax.legend().set_visible(False)
    plt.xticks(rotation=0)

def mapping_quality_pairplot(plot_data, 
                            hierarchy_level, 
                            color_col, 
                            joined_colormap, 
                            avg_corr_thresh=0.4, 
                            agg_prob_thresh=0.4):
    plot_vars = [f'{hierarchy_level}_avg_correlation', f'{hierarchy_level}_bootstrapping_probability', f'{hierarchy_level}_aggregate_probability', color_col]
    g = sns.pairplot(
        plot_data[plot_vars],
        corner=True,
        hue=color_col,
        palette=joined_colormap,
        plot_kws={'s': 2, 'alpha': 0.85}, 
        diag_kind='hist',
        height=2.5
    )
    # figure variables 
    vars_on_grid = getattr(g, 'x_vars', None)
    if vars_on_grid is None:
        vars_on_grid = [v for v in plot_vars if v != color_col]
    # avg correlation
    avg_var = f'{hierarchy_level}_avg_correlation'
    if avg_var in vars_on_grid:
        idx = vars_on_grid.index(avg_var)
        n_grid = len(vars_on_grid)
        # diagonal (hist) - vertical line
        ax_diag = g.axes[idx, idx]
        if ax_diag is not None:
            ax_diag.axvline(avg_corr_thresh, color='red', linestyle='--', alpha=0.7)
        # horizontal lines where this variable is the y-axis (cols < idx)
        for col in range(0, idx):
            ax_row = g.axes[idx, col]
            if ax_row is None:
                continue
            ax_row.axhline(avg_corr_thresh, color='red', linestyle='--', alpha=0.7)
        # vertical lines where this variable is the x-axis (rows > idx)
        for row in range(idx + 1, n_grid):
            ax_col = g.axes[row, idx]
            if ax_col is None:
                continue
            ax_col.axvline(avg_corr_thresh, color='red', linestyle='--', alpha=0.7)

    # aggregate probability
    agg_var = f'{hierarchy_level}_aggregate_probability'
    if agg_var in vars_on_grid:
        idx = vars_on_grid.index(agg_var)
        n_grid = len(vars_on_grid)
        # diagonal (hist) - vertical line
        ax_diag = g.axes[idx, idx]
        if ax_diag is not None:
            ax_diag.axvline(agg_prob_thresh, color='red', linestyle='--', alpha=0.7)
        # horizontal lines where this variable is the y-axis (cols < idx)
        for col in range(0, idx):
            ax_row = g.axes[idx, col]
            if ax_row is None:
                continue
            ax_row.axhline(agg_prob_thresh, color='red', linestyle='--', alpha=0.7)
        # vertical lines where this variable is the x-axis (rows > idx)
        for row in range(idx + 1, n_grid):
            ax_col = g.axes[row, idx]
            if ax_col is None:
                continue
            ax_col.axvline(agg_prob_thresh, color='red', linestyle='--', alpha=0.7)

    g.fig.tight_layout()

def mapping_quality_boxplots(data, 
                            x=['broad_class_name', 'broad_subclass_name', 'supertype_name'],
                            y='cluster_avg_correlation',
                            hue=['broad_class_name', 'broad_subclass_name', 'broad_subclass_name'],
                            colormap='magma',
                            figsize=(20, 10),
                            threshold_val=0.4):
    # Create figure and gridspec properly
    fig = plt.figure(figsize=figsize)
    gs = fig.add_gridspec(2, 4)
    
    # First subplot - broad_class_name
    ax_class = fig.add_subplot(gs[0, 0])
    sns.boxplot(data=data, 
                x=x[0], 
                y=y, 
                hue=hue[0],
                palette=colormap,
                ax=ax_class)
    ax_class.set_title(f'{x[0]}')
    ax_class.tick_params(axis='x', which='both', bottom=False, top=False)
    ax_class.axhline(threshold_val, color='red', linestyle='--', label=f'{y} Threshold')

    # Second subplot - broad_subclass_name
    ax_subclass = fig.add_subplot(gs[0, 1:3])  # spans columns 1-2
    sns.boxplot(data=data, 
                x=x[1], 
                y=y, 
                hue=hue[1],
                palette=colormap,
                ax=ax_subclass)
    ax_subclass.set_title(f'{x[1]}')
    ax_subclass.tick_params(axis='x', which='both', bottom=False, top=False)
    plt.setp(ax_subclass.xaxis.get_majorticklabels(), rotation=45, ha='right')
    ax_subclass.axhline(threshold_val, color='red', linestyle='--', label=f'{y} Threshold')

    # Third subplot - supertype_name
    ax_supertype = fig.add_subplot(gs[1, :])  # spans all columns in row 1
    sns.boxplot(data=data, 
                x=x[2], 
                y=y, 
                hue=hue[2],
                palette=colormap,
                ax=ax_supertype)
    ax_supertype.set_title(f'{x[2]}')
    ax_supertype.tick_params(axis='x', which='both', bottom=False, top=False)
    plt.setp(ax_supertype.xaxis.get_majorticklabels(), rotation=90, ha='right')
    ax_supertype.axhline(threshold_val, color='red', linestyle='--', label=f'{y} Threshold')

    # Remove legends
    for ax in [ax_class, ax_subclass, ax_supertype]:
        legend = ax.get_legend()
        if legend is not None:
            legend.remove()
    
    plt.tight_layout()
    return fig

def plot_mapping_quality(adata, avg_corr_thresh=0.4, agg_prob_thresh=0.25, plot_vars = ['avg_correlation', 'bootstrapping_probability', 'aggregate_probability']):
    fig, axes = plt.subplots(1, 3, figsize=(12, 5))
    for i, var in enumerate(plot_vars):
        ax = axes[i]
        plot_cols = [col for col in adata.obs.columns if col.endswith(var)]
        
        # Melt the data to long format for seaborn boxplot
        melted_data = pd.melt(adata.obs[plot_cols], 
            value_vars=plot_cols,
            var_name='value_type', 
            value_name='value')

        sns.boxplot(data=melted_data, x='value_type', y='value', ax=ax)
        
        if var == 'avg_correlation':
            ax.axhline(avg_corr_thresh, color='red', linestyle='--', label='Correlation Threshold')
            ax.legend(loc='lower left')
        
        if not np.all(adata.obs[plot_cols]==1) and var=='aggregate_probability':
            ax.axhline(agg_prob_thresh, color='red', linestyle='--', label='Correlation Threshold')
            ax.legend(loc='lower left')
        
        # Clean up x-axis labels
        ax.set_xticklabels([label.get_text().split('_')[0] for label in ax.get_xticklabels()], 
                        rotation=45, ha='right')
        ax.set_title(var.replace('_', ' ').title())
    plt.tight_layout()
    return fig

def plot_mapping_quality_comparison(datasets, avg_corr_thresh=0.4, agg_prob_thresh=0.25, plot_vars=['avg_correlation', 'bootstrapping_probability', 'aggregate_probability']): 
    fig, axes = plt.subplots(1,3,figsize=(12, 5))
    for i, var in enumerate(plot_vars):
        ax=axes[i]
        plot_cols = [col for col in datasets[0].columns if col.endswith(var)]
        plot_data = pd.concat([ds[plot_cols + ['dataset']] for ds in datasets], axis=0)

        # Melt the data to long format for seaborn boxplot
        melted_data = pd.melt(plot_data[plot_cols + ['dataset']], 
            id_vars=['dataset'], 
            value_vars=plot_cols,
            var_name='value_type', 
            value_name='value')

        sns.boxplot(data=melted_data, x='value_type', y='value', hue='dataset', ax=ax)
        if var=='avg_correlation':
            ax.axhline(avg_corr_thresh, color='red', linestyle='--', label='Correlation Threshold')  # Fixed variable name
        elif var=='aggregate_probability':
            ax.axhline(agg_prob_thresh, color='red', linestyle='--', label='Aggregate Probability Threshold')  # Fixed variable name
        ax.set_xticklabels([label.get_text().split('_')[0] for label in ax.get_xticklabels()], rotation=45, ha='right')
        ax.set_title(var.replace('_', ' ').title())
        ax.legend(loc='lower left')

    plt.tight_layout()
    return fig

def save_plot(fig, plots_folder, plot_name, save_plots, overwrite_plots, save_format, save_plot_params, close_plots=False):
    if save_plots:
        save_path = plots_folder / f'{plot_name}.{save_format}'
        if overwrite_plots or not save_path.exists():
            fig.savefig(save_path, **save_plot_params)
            print(f"\nSaved to: {'/'.join(save_path.parts[-3:])}")
    if close_plots:
        plt.close()
    else:
        plt.show()