import numpy as np
import pandas as pd
import pytest

pytest.importorskip("sklearn")
pytest.importorskip("matplotlib")

from inhibitory_gmm.subtypes import assign_subclasses, cluster_subtypes, run_subtypes  # noqa: E402

THRESHOLDS = pd.Series({"Pvalb": 100.0, "Sst": 100.0, "Vip": 100.0})


def test_subclass_from_positivity_with_conflicts_to_the_larger_fold():
    genes = pd.DataFrame(
        {"Pvalb": [300, 0, 150, 0, 4], "Sst": [0, 250, 600, 0, 0], "Vip": [0, 0, 0, 90, 0]},
        index=[1, 2, 3, 4, 5],
    )
    labels = assign_subclasses(genes, THRESHOLDS)
    assert labels["subclass"].tolist() == ["Pvalb", "Sst", "Sst", "Other", "Other"]
    assert labels["conflict"].tolist() == [False, False, True, False, False]


def test_fixed_k_splits_by_co_marker_and_names_the_rest_only():
    rng = np.random.default_rng(0)
    n = 60
    calb2 = np.r_[rng.poisson(150, n), rng.poisson(2, n)]
    genes = pd.DataFrame({"Calb2": calb2, "Npy": rng.poisson(3, 2 * n),
                          "Vip": rng.poisson(300, 2 * n)}, index=range(2 * n))
    labels = pd.DataFrame({"subclass": "Vip", "conflict": False}, index=genes.index)
    out, _, chosen, _ = cluster_subtypes(genes, labels, ["Calb2", "Npy", "Vip"], {"Vip": 2})
    assert chosen == {"Vip": 2}
    assert set(out.loc[:n - 1, "subtype"]) == {"Vip_Calb2"}
    assert set(out.loc[n:, "subtype"]) == {"Vip_only"}


def test_run_subtypes_writes_tables_with_gene_order(tmp_path):
    rng = np.random.default_rng(1)
    n = 40
    pivot = pd.DataFrame({
        "R1-488-GFP": rng.poisson(300, 3 * n), "R1-561-Slc17a7": rng.poisson(20, 3 * n),
        "R2-638-Pvalb": np.r_[rng.poisson(400, n), rng.poisson(2, 2 * n)],
        "R3-488-Gad2": rng.poisson(100, 3 * n),
        "R3-561-Vip": np.r_[rng.poisson(2, n), rng.poisson(300, n), rng.poisson(2, n)],
        "R3-638-Sst": np.r_[rng.poisson(2, 2 * n), rng.poisson(300, n)],
        "R2-488-Calb2": rng.poisson(20, 3 * n),
    }, index=pd.Index(range(3 * n), name="cell_id"))
    thresholds = pd.DataFrame({"gene": ["Pvalb", "Sst", "Vip"], "threshold_raw": [100, 100, 100]})
    labels = run_subtypes(pivot, thresholds, tmp_path, "000001", all_gene_k=3, plots=False,
                          fixed_k={"Pvalb": 1, "Sst": 1, "Vip": 1})
    assert labels["subtype"].value_counts().to_dict() == {"Pvalb": n, "Vip": n, "Sst": n}
    table = pd.read_csv(tmp_path / "cell_subtypes_all_genes.csv")
    assert list(table.columns) == ["cell_id", "subclass", "subtype", "GFP", "Slc17a7", "Gad2",
                                   "Calb2", "Pvalb", "Sst", "Vip"]
    clustering = pd.read_csv(tmp_path / "cell_subtypes_clustering_genes.csv")
    assert list(clustering.columns)[3:] == ["Calb2", "Pvalb", "Sst", "Vip"]
    assert (tmp_path / "cell_all_gene_k3.csv").exists()
