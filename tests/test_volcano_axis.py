"""Volcanos : axe Y en p brute, ligne de seuil à la p brute équivalente au seuil p.adj,
même règle DEP pour volcanos, scatters et dashboard."""
import importlib.util
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

import deimos                                                     # noqa: E402
from limma_ebayes import _bh_correction                           # noqa: E402

PARAMS = dict(volcano_use_padj=True, volcano_p_thresh=0.05, volcano_ratio_min=2.0,
              volcano_lfc_min=1.0)


def _df(seed=0, n=2000, n_sig=60):
    rng = np.random.default_rng(seed)
    p = rng.uniform(size=n)
    p[:n_sig] = rng.uniform(0, 1e-4, n_sig)
    lfc = rng.normal(0, 0.5, n)
    lfc[:n_sig] = rng.choice([-2, 2], n_sig) + rng.normal(0, 0.2, n_sig)
    return pd.DataFrame({"name": [f"P{i}" for i in range(n)], "A_vs_B_diff": lfc,
                         "A_vs_B_p.val": p, "A_vs_B_p.adj": _bh_correction(p)})


def test_axis_is_raw_p_and_line_separates_padj_calls():
    df = _df()
    p_axis, cut, name = deimos.volcano_axis(df, "A_vs_B", PARAMS)
    assert name == "p.value"
    np.testing.assert_array_equal(p_axis, df["A_vs_B_p.val"].to_numpy())
    passed = df["A_vs_B_p.adj"].to_numpy() < 0.05
    p = df["A_vs_B_p.val"].to_numpy()
    assert p[passed].max() <= cut < p[~passed].min()   # au-dessus de la ligne <=> p.adj < 0.05


def test_no_significant_protein_means_no_line():
    df = _df(n_sig=0)
    assert (df["A_vs_B_p.adj"] >= 0.05).all()
    _, cut, _ = deimos.volcano_axis(df, "A_vs_B", PARAMS)
    assert cut is None


def test_padj_axis_option_keeps_previous_behaviour():
    df = _df()
    p_axis, cut, name = deimos.volcano_axis(df, "A_vs_B", {**PARAMS, "volcano_y_axis": "padj"})
    assert name == "p.adj" and cut == 0.05
    np.testing.assert_array_equal(p_axis, df["A_vs_B_p.adj"].to_numpy())


def test_raw_p_threshold_line_at_threshold():
    df = _df()
    _, cut, name = deimos.volcano_axis(df, "A_vs_B", {**PARAMS, "volcano_use_padj": False})
    assert name == "p.value" and cut == 0.05


def test_status_uses_chosen_statistic_and_upset_rule():
    df = _df()
    lfc, _, _, up, down = deimos.volcano_status(df, "A_vs_B", PARAMS)
    q = df["A_vs_B_p.adj"].to_numpy()
    assert np.array_equal(up, (q < 0.05) & (lfc >= 1.0))
    assert np.array_equal(down, (q < 0.05) & (lfc <= -1.0))


def test_plots_render(tmp_path):
    df = _df()
    files, facet = deimos.plot_volcanoes(df, ["A_vs_B"], PARAMS, str(tmp_path))
    assert Path(files[0]).exists() and Path(facet).exists()
    df0 = _df(n_sig=0)
    df0.loc[0, "A_vs_B_p.val"] = 0.0                                  # p sous-dépassée
    files, _ = deimos.plot_volcanoes(df0, ["A_vs_B"], PARAMS, str(tmp_path))
    assert Path(files[0]).exists()


def test_dashboard_line_matches_pipeline():
    spec = importlib.util.spec_from_file_location("bd", ROOT / "build_dashboardv7.py")
    bd = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(bd)
    df = _df()
    cuts = bd.build_volcano_bh(df, ["A_vs_B"])["A_vs_B"]
    _, cut, _ = deimos.volcano_axis(df, "A_vs_B", PARAMS)
    assert cuts[bd.BH_GRID.index(0.05)] == pytest.approx(cut, rel=1e-3)
    assert cuts[0] is None or cuts[0] <= cut                           # seuil plus strict
