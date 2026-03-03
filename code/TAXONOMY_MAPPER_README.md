# Taxonomy Mapper - Refactored Pipeline

This directory contains a refactored version of the taxonomy mapping pipeline, organized into modular components with a command-line interface.

## Architecture

The pipeline is now organized into three main modules:

### 1. `config.py` - Configuration Management
- **`PathConfig`**: Manages data, scratch, results, and ABC atlas paths
- **`MappingParams`**: Cell type mapping parameters (normalization, bootstrap settings, etc.)
- **`FilterConfig`**: V1 MERFISH filtering configuration
- **`TaxonomyMapperConfig`**: Complete pipeline configuration with JSON loading

### 2. `taxonomy_mapper.py` - Core Functions
- **`setup_environment()`**: Configure threading for numpy operations
- **`load_abc_cache()`**: Load ABC atlas cache
- **`create_input_adata()`**: Convert cellxgene CSV to AnnData format
- **`prepare_taxonomy_filters()`**: Set up taxonomy node filters
- **`run_mapping()`**: Execute cell type mapping
- **`format_and_save_results()`**: Format and save mapping outputs

### 3. `run_taxonomy_mapper.py` - CLI Entry Point
- Command-line interface with argparse
- Orchestrates the complete pipeline workflow
- Provides flexible parameter override options

## Usage

### Basic Usage

```bash
python run_taxonomy_mapper.py \
    --input-csv /path/to/input_data.csv
```

### Advanced Usage

```bash
python run_taxonomy_mapper.py \
    --config /path/to/params.json \
    --input-csv /root/capsule/data/HCR_767018_Oregano_251104/767018_Oregano_251104_inhibitory_clustered_cellxgene_lognorm.csv \
    --output-name oregano_inhibitory_mapping \
    --drop-layers VISp6a VISp6b \
    --bootstrap-iteration 100 \
    --bootstrap-factor 0.95 \
    --n-runners-up 2 \
    --num-workers 4 \
    --overwrite-all
```

### Command-Line Arguments

#### Required
- `--input-csv`: Path to input CSV file (cellxgene data)

#### Optional
- `--config`: Path to configuration JSON file (default: `/root/capsule/code/params.json`)
- `--output-name`: Custom output folder name
- `--log-norm-data` / `--no-log-norm-data`: Toggle log-normalization handling
- `--drop-layers`: Space-separated list of layers to exclude

#### Overwrite Flags
- `--overwrite-input`: Overwrite existing input h5ad
- `--overwrite-mapping`: Overwrite mapping results
- `--overwrite-formatted`: Overwrite formatted outputs
- `--overwrite-all`: Overwrite everything

#### Mapping Parameters (overrides config file)
- `--bootstrap-iteration`: Number of bootstrap iterations
- `--bootstrap-factor`: Bootstrap sampling factor
- `--n-runners-up`: Number of runner-up cell types
- `--num-workers`: Number of parallel workers

## Pipeline Steps

The pipeline executes five main steps:

1. **Create Input AnnData**: Convert cellxgene CSV to h5ad format
2. **Prepare Taxonomy Filters**: Set up node filtering based on V1 MERFISH data
3. **Validate Input**: Ensure input data matches ABC atlas markers
4. **Run Mapping**: Execute MapMyCells cell type assignment
5. **Format Results**: Save formatted outputs as AnnData

## Configuration File

The pipeline uses `params.json` for default configuration. Example structure:

```json
{
    "paths": {
        "data_root": "/root/capsule/data",
        "scratch_root": "/root/capsule/scratch",
        "abc_path": "/root/capsule/data/abc_atlas"
    },
    "mapping_config": {
        "filter_mapping_v1_types": {
            "enabled": true,
            "h_level": "subclass",
            "min_cells": 0
        },
        "mapping_params": {
            "normalization": "raw",
            "bootstrap_iteration": 100,
            "bootstrap_factor": 0.95,
            "n_runners_up": 2,
            "num_workers": 4
        }
    }
}
```

## Output Structure

```
scratch_root/
└── {output_folder_name}/
    ├── input_data/
    │   └── input_cellxgene.h5ad
    └── mapped_data/
        ├── basic_results.csv
        ├── extended_results.json
        └── mapped_cellxgene.h5ad
```
