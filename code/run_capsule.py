""" top level run script """
import argparse
import os
from pathlib import Path

from codeocean import CodeOcean
from aind_spot_spectral_unmixing.spot_analysis.spot_pipeline import SpotPipelineConfig, run_spot_pipeline_v2
from aind_hcr_data_loader.hcr_dataset import create_hcr_dataset_from_config


DATA_DIR = Path("/root/capsule/data")
CONFIG_PATH = "/root/capsule/code/MOUSE_HCR_CONFIG.json"

os.environ['CODEOCEAN_TOKEN'] = os.environ.get('API_SECRET')
os.environ['CODEOCEAN_DOMAIN'] = "https://codeocean.allenneuraldynamics.org/"

def attach_dataset_assets_from_HCRdataset(ds):
    """
    For each round in the dataset, search for a data asset by the dataset name
    and attach it to the capsule.

    Args:
        ds: HCRDataset object (ds.name is used to query asset titles)
    """

    co = CodeOcean(
        api_token=os.environ["CODEOCEAN_TOKEN"],
        base_url=os.environ.get("CODEOCEAN_DOMAIN"),
    )
    capsule_id = os.environ["CO_CAPSULE_ID"]

    asset_attach_params = []

    for round_key in ds.rounds.keys():
        title = ds.name  # use dataset name as the asset title to search for
        search_results = co.data_assets.search(query=f"name:{title}")

        if search_results and len(search_results) > 0:
            asset_id = search_results[0]["id"]
            asset_attach_params.append({
                "id": asset_id,
                "mount": f"/data/{title}",
            })
        else:
            print(f"Warning: Asset '{title}' not found for round '{round_key}'")

    if asset_attach_params:
        co.capsules.attach_data_assets(
            capsule_id=capsule_id,
            attach_params=asset_attach_params,
        )
        print(f"Attached {len(asset_attach_params)} data assets to capsule {capsule_id}")
    else:
        print("No assets found to attach.")

    return asset_attach_params


def run_pairwise_unmix(mouse_id, input_dir, output_dir):
    """Run the pairwise unmixing pipeline on all datasets found in input_dir."""
    config = SpotPipelineConfig(
        expt_key="pairwise_unmxing",
        min_distances=[1],
        dist_cutoff=1.0,
        do_clean_spots = False,  # default is True, maybe I used this in my test?
        unmixing_method ="reassignment",
        #channel_pairs = # default works well for 6 and 3 channel
        spatial_scale =  (1.0, 0.24, 0.24),
        ratio_spot_filter_method = "95percentile",
        output_base_folder = Path("/root/capsule/scratch")
    )

    ds = create_hcr_dataset_from_config(
        mouse_id,
        data_dir=DATA_DIR,
        config_path=CONFIG_PATH,
    )
    attach_dataset_assets_from_HCRdataset(ds)

    for round_key in ["R2"]:
    #for round_key in ds.rounds.keys():
        results = run_spot_pipeline_v2({mouse_id: ds}, mouse_id, round_key, config)

        #print(f"\nOutput folder: {results['output_folder']}")
        #print(f"Mixed total counts:   {results['mixed_results']['spot_count'].sum()}")
        #print(f"Unmixed total counts: {results['unmixed_results']['spot_count'].sum()}")

    return results

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run capsule pipeline")
    parser.add_argument(
        "--do-pairwise-unmix",
        action="store_true",
        help="Whether to perform new pairwise unmixing analysis",
    )
    parser.add_argument(
        "--do-taxonomy-mapping",
        action="store_true",
        help="Whether to perform taxonomy mapping",
    )
    parser.add_argument(
        "--cell-by-gene-path",
        type=str,
        default=None,
        help="Path to the cell-by-gene file",
    )
    parser.add_argument(
        "--output-dir",
        type=str,
        default="../results/",
        help="Directory to save the output",
    )
    parser.add_argument(
        "--input-dir",
        type=str,
        default="../data/",
        help="Input directory",
    )
    parser.add_argument(
        "--mouse-id",
        type=str,
        default=None,
        help="Mouse ID (required when --do-pairwise-unmix is set)",
    )
    args = parser.parse_args()

    if args.do_pairwise_unmix and not args.mouse_id:
        parser.error("--mouse-id is required when --do-pairwise-unmix is set")
    

    input_dir = Path(args.input_dir).resolve()
    output_dir = Path(args.output_dir).resolve()

    if args.do_pairwise_unmix:
        run_pairwise_unmix(args.mouse_id, input_dir, output_dir)

    if args.do_taxonomy_mapping:
        pass  # TODO