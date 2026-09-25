import json

import pandas as pd
import pytest

from inhibitory_gmm.spot_tables import (
    build_cell_by_gene,
    discover_rounds,
    sort_columns_by_round_channel,
    spots_to_cell_by_gene,
)

GENES_R1 = {"405": "Rn28s", "488": "GFP", "561": "Slc17a7"}
GENES_R2 = {"488": "Calb2", "561": "Npy", "638": "Pvalb"}


def spots(rows):
    return pd.DataFrame(rows, columns=["spot_id", "chan", "cell_id", "chan_561_intensity"])


def write_processed(root, name, rnd, genes, frame, pkl_name=None):
    folder = root / name
    (folder / "image_spot_spectral_unmixing").mkdir(parents=True)
    manifest = {"round": rnd, "gene_dict": {c: {"gene": g} for c, g in genes.items()}}
    (folder / "processing_manifest.json").write_text(json.dumps(manifest))
    frame.to_pickle(folder / "image_spot_spectral_unmixing" / (pkl_name or f"mixed_spots_R{rnd}.pkl"))


def test_spots_to_cell_by_gene_maps_genes_and_applies_slc17a7_floor():
    frame = spots([
        (1, "488", 10, 0.0), (2, "488", 10, 0.0),
        (3, "561", 10, 250.0), (4, "561", 10, 150.0),  # second Slc17a7 spot is below the floor
        (5, "594", 11, 0.0),  # channel with no gene -> dropped
    ])
    long = spots_to_cell_by_gene(frame, "R1", GENES_R1, "839909", "all_spots")
    counts = dict(zip(long["round_chan_gene"], long["spot_count"]))
    assert counts == {"R1-488-GFP": 2, "R1-561-Slc17a7": 1}


def test_filtered_subset_needs_valid_spot_column():
    with pytest.raises(ValueError, match="valid_spot"):
        spots_to_cell_by_gene(spots([(1, "488", 10, 0.0)]), "R1", GENES_R1, "839909", "filtered")
    frame = spots([(1, "488", 10, 0.0), (2, "488", 10, 0.0)]).assign(valid_spot=[True, False])
    long = spots_to_cell_by_gene(frame, "R1", GENES_R1, "839909", "filtered")
    assert long["spot_count"].tolist() == [1]


def test_gene_deletions_follow_the_pairwise_capsule():
    frame = spots([(1, "488", 10, 0.0), (2, "561", 10, 300.0)])
    long = spots_to_cell_by_gene(frame, "R3", {"488": "Tac1", "561": "Gad2"}, "767022", "all_spots")
    assert long["gene"].tolist() == ["Gad2"]


def test_processed_rounds_use_manifest_round_and_stale_pickle_fallback(tmp_path):
    write_processed(tmp_path, "HCR_839909_2026-07-23_13-00-00_processed_2026-07-30_21-17-42", 1,
                    GENES_R1, spots([(1, "488", 10, 0.0)]))
    write_processed(tmp_path, "HCR_839909_2026-07-30_13-00-00_processed_2026-08-03_22-48-02", 2,
                    GENES_R2, spots([(1, "638", 10, 0.0)]), pkl_name="mixed_spots_R-1.pkl")
    write_processed(tmp_path, "HCR_111111_2026-07-30_13-00-00_processed_2026-08-03_22-48-02", 1,
                    GENES_R1, spots([(1, "488", 10, 0.0)]))
    kind, rounds = discover_rounds(tmp_path, "839909")
    assert kind == "processed"
    assert [(r.key, r.path.name) for r in rounds] == [("R1", "mixed_spots_R1.pkl"),
                                                      ("R2", "mixed_spots_R-1.pkl")]

    table = build_cell_by_gene(rounds, "839909", "all_spots")
    assert table.loc[10].to_dict() == {"R1-488-GFP": 1, "R2-638-Pvalb": 1}


def test_auto_prefers_the_pairwise_asset_and_reads_its_tables(tmp_path):
    write_processed(tmp_path, "HCR_837568-01_2026-07-23_13-00-00_processed_2026-07-30_21-17-42", 1,
                    GENES_R1, spots([(1, "488", 10, 0.0)]))
    round_dir = tmp_path / "837568-01_pairwise-unmixing_2026-07-22_00-00-01" / "837568-01_R1"
    round_dir.mkdir(parents=True)
    pd.DataFrame({"cell_id": [7, 7], "gene": ["GFP", "Slac17a7"],
                  "round_chan_gene": ["R1-488-GFP", "R1-561-Slac17a7"], "spot_count": [3, 4]}
                 ).to_csv(round_dir / "mixed_all_spots_cell_by_gene.csv", index=False)
    kind, rounds = discover_rounds(tmp_path, "837568-01", spots="all_spots")
    assert kind == "pairwise" and rounds[0].kind == "table"
    table = build_cell_by_gene(rounds, "837568-01", "all_spots")
    assert table.loc[7].to_dict() == {"R1-488-GFP": 3, "R1-561-Slc17a7": 4}

    kind, _ = discover_rounds(tmp_path, "837568-01", source="processed")
    assert kind == "processed"


def test_spot_parquet_asset_is_matched_by_subject(tmp_path):
    for subject in ("839909", "800995"):
        asset = tmp_path / f"cell-types-and-learning_spot-parquet_{subject}"
        asset.mkdir()
        (asset / "subject.json").write_text(json.dumps({"subject_id": subject}))
        (asset / "meta_R1.json").write_text(json.dumps({"round": 1, "key": "R1", "genes": GENES_R1}))
        spots([(1, "488", 10, 0.0)]).to_parquet(asset / "spots_R1.parquet")
    kind, rounds = discover_rounds(tmp_path, "839909")
    assert kind == "spot_parquet"
    assert [r.source_asset for r in rounds] == ["cell-types-and-learning_spot-parquet_839909"]
    table = build_cell_by_gene(rounds, "839909", "all_spots")
    assert table.loc[10].to_dict() == {"R1-488-GFP": 1}


def test_missing_inputs_raise(tmp_path):
    with pytest.raises(FileNotFoundError):
        discover_rounds(tmp_path, "839909")


def test_columns_sort_by_round_then_channel():
    cxg = pd.DataFrame(columns=["R2-488-Calb2", "R1-561-Slc17a7", "other", "R1-488-GFP"])
    assert list(sort_columns_by_round_channel(cxg).columns) == [
        "R1-488-GFP", "R1-561-Slc17a7", "R2-488-Calb2", "other"
    ]
