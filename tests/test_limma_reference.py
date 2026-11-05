"""limma_ebayes.py contre les références externes.

1. R limma 3.58.1 / DEqMS 1.31.0 : sorties stockées dans tests/fixtures/limma_r
   (données simulées ; R n'est pas nécessaire pour lancer ces tests).
2. InMoose (portage Python de limma), si installé : comparaison directe.
"""
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

import limma_ebayes as le                                        # noqa: E402

FIX = ROOT / "tests" / "fixtures" / "limma_r"
CASES = sorted(p.name for p in FIX.glob("case_*"))


def _read(case, name):
    p = FIX / case / name
    if not p.exists():
        return None
    return pd.read_csv(p, header=None, na_values=["NA"]).to_numpy(float)


def _inputs(case):
    E = np.loadtxt(FIX / case / "expr.csv.gz", delimiter=",", ndmin=2)
    X = np.loadtxt(FIX / case / "design.csv", delimiter=",", ndmin=2)
    C = np.loadtxt(FIX / case / "contrasts.csv", delimiter=",", ndmin=2)
    return E, X, C


def _close_logp(a, b, tol=1e-9):
    la = -np.log10(np.maximum(np.asarray(a, float), 1e-300))
    lb = -np.log10(np.maximum(np.asarray(b, float), 1e-300))
    np.testing.assert_allclose(la, lb, rtol=tol, atol=tol)


@pytest.mark.parametrize("case", CASES)
def test_ebayes_matches_r_limma(case):
    E, X, C = _inputs(case)
    f = le.ebayes(le.contrasts_fit(le.lm_fit(E, X), C))
    hyper = _read(case, "R/hyper.csv.gz").ravel()          # df.prior, s2.prior, df.total, df.res
    np.testing.assert_allclose(f["coefficients"], _read(case, "R/coef.csv.gz"), atol=1e-9)
    np.testing.assert_allclose(f["sigma"], _read(case, "R/sigma.csv.gz").ravel(), rtol=1e-9,
                               atol=1e-12)
    assert (np.isinf(f["d0"]) and np.isinf(hyper[0])) or np.isclose(f["d0"], hyper[0], rtol=1e-9)
    assert np.isclose(f["s0_sq"], hyper[1], rtol=1e-9)
    assert np.isclose(f["df_total"], hyper[2], rtol=1e-9)
    np.testing.assert_allclose(f["t_stat"], _read(case, "R/t.csv.gz"), rtol=1e-8, atol=1e-9)
    _close_logp(f["p_value"], _read(case, "R/p.csv.gz"))
    np.testing.assert_allclose(f["p_adj"], _read(case, "R/padj.csv.gz"), rtol=1e-8, atol=1e-12)
    g = le.ebayes(le.contrasts_fit(le.lm_fit(E, X), C), fdr_global=True)
    np.testing.assert_allclose(g["p_adj"].ravel(order="F"),
                               _read(case, "R/padj_global.csv.gz").ravel(), rtol=1e-8, atol=1e-12)


def test_deqms_matches_r():
    pytest.importorskip("skmisc.loess")
    case = "case_deqms"
    E, X, C = _inputs(case)
    k = np.loadtxt(FIX / case / "counts.csv")
    q = le.spectra_count_ebayes(le.contrasts_fit(le.lm_fit(E, X), C), k)
    assert q["loess_engine"].startswith("R loess")
    np.testing.assert_allclose(q["loess_fit"], _read(case, "R/sca_loess.csv.gz").ravel(),
                               rtol=1e-9, atol=1e-10)
    assert q["df_prior"] == pytest.approx(_read(case, "R/sca_dfprior.csv.gz").item())
    np.testing.assert_allclose(q["t_stat"], _read(case, "R/sca_t.csv.gz"), rtol=1e-8, atol=1e-9)
    _close_logp(q["p_value"], _read(case, "R/sca_p.csv.gz"))
    np.testing.assert_allclose(q["p_adj"], _read(case, "R/sca_padj.csv.gz"), rtol=1e-8, atol=1e-12)


def test_deqms_without_scikit_misc_stays_close(monkeypatch):
    """Sans scikit-misc : loess directe, écart faible à R et mêmes appels ou presque."""
    import builtins
    real_import = builtins.__import__

    def no_skmisc(name, *a, **kw):
        if name.startswith("skmisc"):
            raise ImportError(name)
        return real_import(name, *a, **kw)

    monkeypatch.setattr(builtins, "__import__", no_skmisc)
    case = "case_deqms"
    E, X, C = _inputs(case)
    k = np.loadtxt(FIX / case / "counts.csv")
    q = le.spectra_count_ebayes(le.contrasts_fit(le.lm_fit(E, X), C), k)
    assert q["loess_engine"].startswith("direct loess")
    assert np.max(np.abs(q["loess_fit"] - _read(case, "R/sca_loess.csv.gz").ravel())) < 0.05
    sig_r = _read(case, "R/sca_padj.csv.gz") < 0.05
    assert ((q["p_adj"] < 0.05) != sig_r).sum() <= 2


def test_deqms_prior_df_grid():
    """Grille séquentielle de DEqMS : arrêt dès que |m - trigamma(d0/2)| remonte."""
    from scipy.special import polygamma
    for target in (1.7, 4.2, 12.0):
        m = float(polygamma(1, target / 2))
        assert le._deqms_prior_df(m, 1000) == pytest.approx(target, abs=0.05)
    assert le._deqms_prior_df(-1.0, 3) == pytest.approx(3.0)   # croît sans fin : bout de grille


def test_bh_matches_p_adjust_with_missing_values():
    p = np.array([0.01, np.nan, 0.04, 0.03, 0.005, np.nan, 0.5])
    # R : p.adjust(c(0.01, NA, 0.04, 0.03, 0.005, NA, 0.5), "BH")
    ref = np.array([0.025, np.nan, 0.05, 0.05, 0.025, np.nan, 0.5])
    np.testing.assert_allclose(le._bh_correction(p), ref, equal_nan=True)


def test_rank_deficient_design_df():
    rng = np.random.default_rng(0)
    E = rng.normal(size=(50, 6))
    X = np.c_[np.r_[np.ones(3), np.zeros(3)], np.r_[np.zeros(3), np.ones(3)], np.ones(6)]
    assert le.lm_fit(E, X)["df_residual"] == 4                  # rang 2, pas 3 colonnes


@pytest.mark.filterwarnings("ignore::RuntimeWarning")
def test_inmoose_agrees():
    L = pytest.importorskip("inmoose.limma")
    import contextlib
    import io
    for case in ("case_2g", "case_4g", "case_paired", "case_zero"):
        E, X, C = _inputs(case)
        f = le.ebayes(le.contrasts_fit(le.lm_fit(E, X), C))
        Ed = pd.DataFrame(E, columns=[f"s{i}" for i in range(E.shape[1])])
        Xd = pd.DataFrame(X, index=Ed.columns, columns=[f"c{i}" for i in range(X.shape[1])])
        with contextlib.redirect_stdout(io.StringIO()):
            fi = L.lmFit(Ed, Xd)
            Cd = pd.DataFrame(C, index=fi.coefficients.columns,
                              columns=[f"k{i}" for i in range(C.shape[1])])
            fi = L.eBayes(L.contrasts_fit(fi, contrasts=Cd))
        np.testing.assert_allclose(f["t_stat"], np.asarray(fi.t), rtol=1e-8, atol=1e-9)
        _close_logp(f["p_value"], np.asarray(fi.p_value))
