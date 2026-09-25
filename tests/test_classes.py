import pandas as pd

from inhibitory_gmm.classes import assign_classes


def test_inhibitory_takes_precedence_over_slc17a7():
    cells = pd.Index([1, 2, 3, 4, 5])
    out = assign_classes(cells, inhibitory=pd.Index([1, 2]), slc17a7_positive=pd.Index([2, 3]))
    assert out["class"].tolist() == ["Inhibitory", "Inhibitory", "Excitatory", "Unassigned",
                                     "Unassigned"]
    assert out["slc17a7_positive"].tolist() == [False, True, True, False, False]
