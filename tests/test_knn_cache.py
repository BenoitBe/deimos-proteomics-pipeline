"""Imputation mixte : la composante kNN pré-calculée (score de robustesse) donne
exactement la même imputation que son recalcul à chaque tirage."""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import deimos  # noqa: E402


def _data(n=300, reps=4, conds=("A", "B", "C"), seed=0):
    rng = np.random.default_rng(seed)
    labels = [f"{c}{r}" for c in conds for r in range(reps)]
    Y = rng.normal(14, 2, (n, 1)) + rng.normal(0, 0.4, (n, len(labels)))
    Y[rng.random(Y.shape) > 1 / (1 + np.exp(-(Y - 12.5) * 1.5))] = np.nan
    Y[:20, :reps] = np.nan                                  # absences totales (MNAR)
    design = pd.DataFrame({"label": labels, "condition": [c for c in conds for _ in range(reps)]})
    return pd.DataFrame(Y, columns=labels), design


def test_cached_knn_gives_identical_imputation():
    mat, design = _data()
    knn = deimos.knn_component(mat)
    for method in ("mixed", "mixed_legacy"):
        for seed in (1, 7, 123):
            a = deimos.impute(mat, method, design, rng=np.random.default_rng(seed))
            b = deimos.impute(mat, method, design, rng=np.random.default_rng(seed), knn_full=knn)
            np.testing.assert_array_equal(a.values, b.values)
            assert not np.isnan(b.values).any()


def test_knn_is_deterministic():
    mat, _ = _data(seed=3)
    np.testing.assert_array_equal(deimos.knn_component(mat), deimos.knn_component(mat))
