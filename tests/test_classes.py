import numpy as np
import pandas as pd
import pytest

pytest.importorskip("sklearn")

from inhibitory_gmm.classes import assign_classes, cluster_excitatory  # noqa: E402


def test_classes_mark_collisions_ambiguous():
    cells = pd.Index([1, 2, 3, 4, 5, 6])
    out = assign_classes(cells, inhibitory=pd.Index([1, 2]), gate=pd.Index([1, 2, 3]),
                         slc17a7_positive=pd.Index([2, 3, 4]))
    assert out["class"].tolist() == ["Inhibitory", "Ambiguous", "Excitatory", "Excitatory",
                                     "Unassigned", "Unassigned"]
    assert out["inhibitory_gate"].tolist() == [True, True, True, False, False, False]


def test_excitatory_clusters_are_named_by_brightness():
    rng = np.random.default_rng(0)
    n = 80
    genes = pd.DataFrame({"GFP": rng.poisson(300, 2 * n),
                          "Slc17a7": np.r_[rng.poisson(100, n), rng.poisson(500, n)],
                          "Calb2": rng.poisson(10, 2 * n)})
    names, rank, distance, order, features, score = cluster_excitatory(genes, 2, exclude=("GFP",))
    assert features == ["Slc17a7", "Calb2"]
    assert set(names.iloc[:n]) == {"Exc_2"} and set(names.iloc[n:]) == {"Exc_1"}
    assert order == ["Exc_1", "Exc_2"] and score > 0.5
    assert len(distance) == 2 * n
