import json

import pandas as pd
import pytest

import run_capsule
from cell_typing_table import build_cell_typing_table


def test_shipped_presets_are_valid():
    for path in sorted(run_capsule.CONFIG_DIR.glob("*.json")):
        run_capsule.resolve_settings({}, path.stem)


def test_explicit_args_override_preset_which_overrides_defaults():
    settings = run_capsule.resolve_settings(
        {"gmm_k": 8, "spots": None, "mouse_id": "839909"}, "p3_mixed_inhibitory_gmm"
    )
    assert settings["run_inhibitory_gmm"] is True
    assert settings["gmm_spots"] == "all_spots"
    assert settings["gmm_k"] == 8
    assert settings["spots"] == "filtered"
    assert settings["run_mapmycells"] is False


def test_unknown_preset_keys_and_bad_choices_are_rejected(tmp_path):
    bad_key = tmp_path / "bad_key.json"
    bad_key.write_text(json.dumps({"run_inhibitory_gmn": True}))
    with pytest.raises(ValueError, match="unknown keys"):
        run_capsule.resolve_settings({}, str(bad_key))
    with pytest.raises(ValueError, match="gmm_source"):
        run_capsule.resolve_settings({"gmm_source": "s3"}, None)
    with pytest.raises(FileNotFoundError):
        run_capsule.resolve_settings({}, "no_such_preset")


def test_cell_typing_table_includes_gmm_labels(tmp_path):
    folder = tmp_path / "inhibitory_gmm" / "inhibitory_cells_mixed_all_spots"
    folder.mkdir(parents=True)
    pd.DataFrame({"cell_id": [5, 9], "cluster": [3, 0]}).to_csv(
        folder / "mixed_cluster_labels_all_spots.csv", index=False
    )
    out = build_cell_typing_table(tmp_path, "839909", gmm_spots="all_spots")
    table = pd.read_csv(out)
    assert list(table.columns) == ["cell_id", "mouse_id", "gmm_cluster", "gmm_inhibitory"]
    assert table.set_index("cell_id")["gmm_cluster"].to_dict() == {5: 3, 9: 0}
    assert build_cell_typing_table(tmp_path / "empty", "839909") is None

    subtypes = tmp_path / "inhibitory_gmm" / "subtypes"
    subtypes.mkdir()
    pd.DataFrame({"cell_id": [5, 9], "subclass": ["Vip", "Sst"], "subtype": ["Vip_Calb2", "Sst"],
                  "Vip": [300, 0]}).to_csv(subtypes / "cell_subtypes_clustering_genes.csv", index=False)
    table = pd.read_csv(build_cell_typing_table(tmp_path, "839909", gmm_spots="all_spots"))
    assert table.set_index("cell_id")["gmm_subtype"].to_dict() == {5: "Vip_Calb2", 9: "Sst"}


def test_cell_typing_table_covers_all_cells_when_classes_ran(tmp_path):
    folder = tmp_path / "inhibitory_gmm"
    folder.mkdir()
    pd.DataFrame({
        "cell_id": [5, 9, 11], "class": ["Inhibitory", "Excitatory", "Unassigned"],
        "inhibitory_gate": [True, False, False], "slc17a7_positive": [False, True, False],
        "excitatory_cluster": [None, "e0_Slc17a7", None], "subclass": ["Vip", None, None],
        "subtype": ["Vip_Calb2", None, None], "gmm_cluster": [3, None, None],
        "GFP": [1, 2, 3], "Slc17a7": [4, 300, 0],
    }).to_csv(folder / "cell_classes.csv", index=False)
    table = pd.read_csv(build_cell_typing_table(tmp_path, "839909", gmm_spots="all_spots"))
    assert len(table) == 3
    assert list(table.columns)[:8] == ["cell_id", "mouse_id", "class", "subclass", "subtype",
                                       "excitatory_cluster", "inhibitory_gate", "slc17a7_positive"]
    assert list(table.columns)[-2:] == ["GFP", "Slc17a7"]
    assert table["gmm_inhibitory"].tolist() == [True, False, False]
