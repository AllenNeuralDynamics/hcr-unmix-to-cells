"""Configuration management for taxonomy mapping."""

import json
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Tuple


@dataclass
class PathConfig:
    """Path configuration for data, scratch, and results directories."""
    data_root: Path
    scratch_root: Path
    results_root: Path
    abc_path: Path
    mmc_path: Path

    @classmethod
    def from_dict(cls, config: Dict) -> 'PathConfig':
        """Create PathConfig from dictionary."""
        return cls(
            data_root=Path(config['data_root']),
            scratch_root=Path(config['scratch_root']),
            results_root=Path(config['results_root']),
            abc_path=Path(config['abc_path']),
            mmc_path=Path(config['mmc_path'])
        )


@dataclass
class MappingParams:
    """Parameters for cell type mapping."""
    normalization: str = 'raw'
    bootstrap_iteration: int = 100
    bootstrap_factor: float = 0.95
    n_runners_up: int = 2
    nodes_to_drop: Optional[List[Tuple[str, str]]] = None
    num_workers: int = 4
    
    # Additional params needed by format_mapping_outputs
    drop_level: str = 'supertype'
    mapping_type: str = 'hrc'
    clobber: str = '1'
    chunk_size: int = 15000
    rng_seed: int = 42
    
    @classmethod
    def from_dict(cls, config: Dict) -> 'MappingParams':
        """Create MappingParams from dictionary."""
        return cls(
            normalization=config.get('normalization', 'raw'),
            bootstrap_iteration=config.get('bootstrap_iteration', 100),
            bootstrap_factor=config.get('bootstrap_factor', 0.95),
            n_runners_up=config.get('n_runners_up', config.get('n_runner_ups', 2)),
            num_workers=config.get('num_workers', 4),
            drop_level=config.get('drop_level', 'supertype'),
            mapping_type=config.get('mapping_type', 'hrc'),
            clobber=str(config.get('clobber', '1')),
            chunk_size=config.get('chunk_size', 15000),
            rng_seed=config.get('rng_seed', 42)
        )
    
    def to_dict(self) -> Dict:
        """Convert to dictionary for passing to mapping functions."""
        return {
            'normalization': self.normalization,
            'bootstrap_iteration': self.bootstrap_iteration,
            'bootstrap_factor': self.bootstrap_factor,
            'n_runners_up': self.n_runners_up,
            'num_workers': self.num_workers,
            'drop_level': self.drop_level,
            'mapping_type': self.mapping_type,
            'clobber': self.clobber,
            'chunk_size': self.chunk_size,
            'rng_seed': self.rng_seed
        }


@dataclass
class FilterConfig:
    """Configuration for V1 MERFISH filtering."""
    enabled: bool = True
    h_level: str = 'subclass'
    min_cells: int = 0
    saved_df_name: Optional[str] = 'v1_merfish_cells.csv'
    
    @classmethod
    def from_dict(cls, config: Dict) -> 'FilterConfig':
        """Create FilterConfig from dictionary."""
        return cls(
            enabled=config.get('enabled', True),
            h_level=config.get('h_level', 'subclass'),
            min_cells=config.get('min_cells', 0),
            saved_df_name=config.get('saved_df_name', 'v1_merfish_cells.csv')
        )


@dataclass
class TaxonomyMapperConfig:
    """Complete configuration for taxonomy mapping pipeline."""
    paths: PathConfig
    mapping_params: MappingParams
    filter_config: FilterConfig
    
    # Input/output settings
    dataset_folder: Path = None
    data_csv: str = None
    log_norm_data: bool = True
    output_folder_name: str = None
    drop_layers: Optional[List[str]] = None
    drop_nodes_dict: Optional[Dict[str, List[str]]] = None
    
    # Overwrite flags
    overwrite_input_adata: bool = True
    overwrite_mapping_results: bool = True
    overwrite_formatted_outputs: bool = True
    
    # Folder names
    input_data_folder_name: str = 'input_data'
    mapped_data_folder_name: str = 'mapped_data'
    input_h5ad_name: str = 'input_cellxgene.h5ad'
    mapped_data_h5ad_name: str = 'mapped_cellxgene.h5ad'
    basic_results_name: str = 'basic_results.csv'
    extended_results_name: str = 'extended_results.json'

    @classmethod
    def from_json(cls, config_path: str) -> 'TaxonomyMapperConfig':
        """Load configuration from JSON file."""
        with open(config_path, 'r') as f:
            config = json.load(f)
        
        paths = PathConfig.from_dict(config['paths'])
        
        mapping_config = config.get('mapping_config', {})
        mapping_params = MappingParams.from_dict(mapping_config.get('mapping_params', {}))
        
        filter_config_dict = mapping_config.get('filter_mapping_v1_types', {})
        filter_config = FilterConfig.from_dict(filter_config_dict)
        
        return cls(
            paths=paths,
            mapping_params=mapping_params,
            filter_config=filter_config,
            input_data_folder_name=mapping_config.get('input_data_folder_name', 'input_data'),
            mapped_data_folder_name=mapping_config.get('mapped_data_folder_name', 'mapped_data'),
            input_h5ad_name=mapping_config.get('input_h5ad_name', 'input_cellxgene.h5ad'),
            mapped_data_h5ad_name=mapping_config.get('mapped_data_h5ad_name', 'mapped_cellxgene.h5ad'),
            basic_results_name=mapping_config.get('basic_results_name', 'basic_results.csv'),
            extended_results_name=mapping_config.get('extended_results_name', 'extended_results.json'),
            drop_nodes_dict=mapping_config.get('drop_nodes_dict')
        )

    def get_output_paths(self) -> Tuple[Path, Path, Path]:
        """Get input folder, output folder, and dataset output folder paths."""
        if self.output_folder_name is None:
            if self.data_csv:
                stem = Path(self.data_csv).stem
            else:
                stem = 'default_output'
            self.output_folder_name = f'{stem}_mapping'
        
        dataset_output_folder = self.paths.scratch_root / self.output_folder_name
        input_folder = dataset_output_folder / self.input_data_folder_name
        output_folder = dataset_output_folder / self.mapped_data_folder_name
        
        return input_folder, output_folder, dataset_output_folder
