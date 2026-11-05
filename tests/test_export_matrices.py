"""Matrices pour d'autres analyses (Mfuzz, MOFA+/MEFISTO, multi-omique)."""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import deimos  # noqa: E402


def _inputs(n=50, seed=0):
    rng = np.random.default_rng(seed)
    labels = [f"S{i}" for i in range(9)]
    design = pd.DataFrame({"label": labels, "condition": np.repeat(["T0", "T2", "T8"], 3),
                           "replicate": np.tile([1, 2, 3], 3)})
    X = rng.normal(20, 1, (n, 9))
    X[rng.random(X.shape) < 0.1] = np.nan
    mat_filt = pd.DataFrame(X, columns=labels)
    mat_imp = mat_filt.fillna(15.0)
    meta = pd.DataFrame({"name": [f"P{i}" for i in range(n)],
                         "Protein.Group": [f"P{i};Q{i}" for i in range(n)],
                         "Genes": [np.nan if i % 7 == 0 else f"G{i}" for i in range(n)],
                         "First.Protein.Description": [f"prot {i}" for i in range(n)]})
    meta.loc[3, "name"] = "P1"                                   # doublon d'identifiant
    E = np.where(np.isnan(X), 16.0, X)
    dpc = {"E": E, "SE": np.full(X.shape, 0.2), "imputed": np.isnan(X)}
    return mat_filt, mat_imp, meta, design, dpc


def test_files_shapes_and_alignment(tmp_path):
    mat_filt, mat_imp, meta, design, dpc = _inputs()
    d = deimos.export_matrices(str(tmp_path), mat_filt, mat_imp, meta,
                               design.iloc[::-1].reset_index(drop=True),   # ordre différent
                               dpc_primary=dpc, params={"impute_method": "mixed", "seed": 7})
    d = Path(d)
    for f in ("log2_filtered_NA.tsv", "log2_imputed.tsv", "dpc_abundance.tsv",
              "dpc_standard_error.tsv", "dpc_observed.tsv", "condition_means_dpc.tsv",
              "samples.tsv", "README_matrices.txt"):
        assert (d / f).exists(), f
    na = pd.read_csv(d / "log2_filtered_NA.tsv", sep="\t")
    s = pd.read_csv(d / "samples.tsv", sep="\t")
    assert list(s["label"]) == list(mat_filt.columns)            # ordre des colonnes
    assert na["protein"].is_unique and len(na) == len(mat_filt)
    np.testing.assert_allclose(na[list(mat_filt.columns)].to_numpy(),
                               mat_filt.round(4).to_numpy(), equal_nan=True)
    assert (na["num_NA"] == mat_filt.isna().sum(axis=1)).all()
    obs = pd.read_csv(d / "dpc_observed.tsv", sep="\t")
    assert obs[mat_filt.columns[0]].dtype.kind == "i"            # 0/1 entiers
    assert (obs[list(mat_filt.columns)].to_numpy() == (~mat_filt.isna()).astype(int).to_numpy()).all()
    cm = pd.read_csv(d / "condition_means_dpc.tsv", sep="\t")
    np.testing.assert_allclose(cm["T2"], dpc["E"][:, 3:6].mean(axis=1).round(4))
    assert list(cm.columns[-3:]) == ["T0", "T2", "T8"]           # ordre de samples.tsv
    assert "MIXED, seed 7" in (d / "README_matrices.txt").read_text()


def test_without_dpc(tmp_path):
    mat_filt, mat_imp, meta, design, _ = _inputs()
    d = Path(deimos.export_matrices(str(tmp_path), mat_filt, mat_imp, meta, design))
    assert not (d / "dpc_abundance.tsv").exists() and (d / "log2_imputed.tsv").exists()
    assert "dpc_abundance" not in (d / "README_matrices.txt").read_text()
