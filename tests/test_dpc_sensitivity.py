"""Sensitivity to missing values (dpc_sensitivity.py): DPC estimation, imputation-free
analysis, differential detection. Simulated data with known truth; fast."""
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

import dpc_sensitivity as ds                                       # noqa: E402
from limma_ebayes import make_all_contrasts, make_design_matrix   # noqa: E402

B0, B1 = -0.8 * 19.0, 0.8          # 50 % detection at log2 = 19


def simulate(seed=1, n_prot=1500, n_rep=4, de_frac=0.12, onoff_frac=0.03):
    rng = np.random.default_rng(seed)
    groups = ["Ctrl"] * n_rep + ["Treat"] * n_rep
    mu = rng.normal(22.0, 2.5, n_prot)
    sd = np.sqrt(0.15 / rng.chisquare(6, n_prot) * 6)
    lfc = np.zeros(n_prot)
    de = rng.random(n_prot) < de_frac
    lfc[de] = rng.choice([-1, 1], de.sum()) * rng.uniform(0.8, 3.0, de.sum())
    onoff = (~de) & (rng.random(n_prot) < onoff_frac)
    lfc[onoff] = -8.0
    Y = np.empty((n_prot, 2 * n_rep))
    Y[:, :n_rep] = mu[:, None] + rng.normal(0, 1, (n_prot, n_rep)) * sd[:, None]
    Y[:, n_rep:] = (mu + lfc)[:, None] + rng.normal(0, 1, (n_prot, n_rep)) * sd[:, None]
    obs = rng.random(Y.shape) < ds.expit(B0 + B1 * Y)
    Yobs = np.where(obs, Y, np.nan)
    keep = (~np.isnan(Yobs)).any(axis=1)
    return {"Y": Yobs[keep], "truth": (de | onoff)[keep], "onoff": onoff[keep],
            "groups": groups, "n_rep": n_rep}


def filtered(sim):
    g = np.array(sim["groups"])
    keep = np.zeros(len(sim["Y"]), bool)
    for gg in np.unique(g):
        keep |= (~np.isnan(sim["Y"][:, g == gg])).sum(axis=1) >= sim["n_rep"] - 1
    return keep


@pytest.fixture(scope="module")
def sim():
    return simulate()


@pytest.fixture(scope="module")
def dpc(sim):
    return ds.estimate_dpc(sim["Y"])


class TestDPC:
    def test_recovers_the_true_curve(self, dpc):
        assert abs(-dpc["beta0"] / dpc["beta1"] - 19.0) < 0.3      # 50 % detection point
        assert 0.6 < dpc["beta1"] < 0.95

    def test_intercept_matches_overall_detection_rate(self, sim):
        Y = sim["Y"][:1000]                       # < 10 000 cells: exact equation
        b0 = ds.estimate_dpc_intercept(Y, 0.7)
        yc = ds.exp_tilt_complete(Y, 0.7)
        assert ds.expit(b0 + 0.7 * yc).mean() == pytest.approx((~np.isnan(Y)).mean(), rel=1e-6)
        # > 10 000 cells: binned like limpa, close to the exact solution
        b_big = ds.estimate_dpc_intercept(sim["Y"], 0.7)
        yc = ds.exp_tilt_complete(sim["Y"], 0.7)
        assert ds.expit(b_big + 0.7 * yc).mean() == pytest.approx(
            (~np.isnan(sim["Y"])).mean(), rel=1e-3)

    def test_exp_tilt_fills_below_the_observed_mean(self, sim):
        yc = ds.exp_tilt_complete(sim["Y"], 0.8)
        miss = np.isnan(sim["Y"])
        assert not np.isnan(yc).any()
        row_means = np.nanmean(sim["Y"], axis=1)
        assert np.mean(yc[miss] < np.broadcast_to(row_means[:, None], yc.shape)[miss]) > 0.9

    def test_complete_data_is_refused(self):
        with pytest.raises(ValueError, match="no missing"):
            ds.estimate_dpc(np.random.default_rng(0).normal(20, 1, (50, 6)))


class TestCompletion:
    def test_observed_kept_missing_below_and_uncertain(self, sim, dpc):
        Y = sim["Y"][filtered(sim)]
        c = ds.dpc_complete(Y, dpc["beta0"], dpc["beta1"])
        obs = ~np.isnan(Y)
        assert np.allclose(c["E"][obs], Y[obs], atol=1e-3)
        assert np.median(c["SE"][~obs]) > 10 * np.median(c["SE"][obs])
        below = c["E"][~obs] < np.broadcast_to(c["row_mean"][:, None], Y.shape)[~obs]
        assert below.all()


class TestImputationFreeDE:
    def test_controls_fdr_and_finds_partial_missing(self, sim, dpc):
        keep = filtered(sim)
        Y, truth = sim["Y"][keep], sim["truth"][keep]
        X, g = make_design_matrix(sim["groups"], control_condition="Ctrl")
        C, names = make_all_contrasts(g)
        c = ds.dpc_complete(Y, dpc["beta0"], dpc["beta1"])
        r = ds.dpc_de(c["E"], c["SE"], c["imputed"], X, C)
        sig = r["p_adj"][:, 0] < 0.05
        fdp = (sig & ~truth).sum() / max(sig.sum(), 1)
        assert fdp < 0.10
        assert (sig & truth).sum() / truth.sum() > 0.5
        assert r["n_df_corrected"] > 0          # whole-group absences detected

    def test_no_df_correction_without_absent_group(self):
        rng = np.random.default_rng(3)
        Y = rng.normal(22, 1, (300, 6))
        Y[rng.random(Y.shape) < 0.05] = np.nan
        Y[:, :3][np.isnan(Y[:, :3]).all(axis=1)] = 22.0
        X, g = make_design_matrix(["A"] * 3 + ["B"] * 3)
        C, _ = make_all_contrasts(g)
        c = ds.dpc_complete(Y, -15.0, 0.8)
        all_missing_group = np.isnan(Y[:, :3]).all(1) | np.isnan(Y[:, 3:]).all(1)
        r = ds.dpc_de(c["E"], c["SE"], c["imputed"], X, C)
        assert r["n_df_corrected"] == int(all_missing_group.sum())

    def test_lowess_span_matches_limma(self):
        assert ds.choose_lowess_span(2500) == pytest.approx(0.4900092, abs=1e-6)   # limma 3.58
        assert ds.choose_lowess_span(10) == 1.0


class TestDifferentialDetection:
    def test_onoff_found_without_false_positive(self, sim, dpc):
        keep = filtered(sim)
        Y, truth, onoff = sim["Y"][keep], sim["truth"][keep], sim["onoff"][keep]
        groups = np.array([0] * 4 + [1] * 4)
        dd = ds.differential_detection(Y, groups, [("Treat_vs_Ctrl", 1, 0)],
                                       dpc["beta0"], dpc["beta1"])
        r = dd["by_contrast"]["Treat_vs_Ctrl"]
        sig = r["p_adj"] < 0.05
        assert (sig & ~truth).sum() == 0
        assert (sig & onoff).sum() >= 0.3 * onoff.sum()
        assert np.all(r["higher_in"][sig & onoff] == 0)      # present in Ctrl

    def test_conservative_under_the_null(self, dpc):
        sim0 = simulate(seed=11, de_frac=0.0, onoff_frac=0.0)
        Y = sim0["Y"][filtered(sim0)]
        dd = ds.differential_detection(Y, np.array([0] * 4 + [1] * 4),
                                       [("B_vs_A", 1, 0)], dpc["beta0"], dpc["beta1"])
        p = dd["by_contrast"]["B_vs_A"]["p"]
        has_missing = np.isnan(Y).any(axis=1)
        assert np.mean(p[has_missing] < 0.05) <= 0.05
        assert (dd["by_contrast"]["B_vs_A"]["p_adj"] < 0.05).sum() == 0

    def test_absence_near_the_detection_limit_is_not_evidence(self):
        """4/4 detected at the 50 % point vs 0/4: compatible with chance."""
        Y = np.array([[19.0, 19.2, 18.9, 19.1, np.nan, np.nan, np.nan, np.nan],
                      [26.0, 26.2, 25.9, 26.1, np.nan, np.nan, np.nan, np.nan]])
        Y = np.vstack([Y, np.random.default_rng(1).normal(23, 0.3, (40, 8))])
        dd = ds.differential_detection(Y, np.array([0] * 4 + [1] * 4),
                                       [("B_vs_A", 1, 0)], B0, B1)
        p = dd["by_contrast"]["B_vs_A"]["p"]
        assert p[0] > 0.05             # near the LOD: not significant
        assert p[1] < 1e-4             # far above it: strong evidence

    def test_sporadic_dropout_is_not_evidence(self):
        """Real data: some abundant proteins miss in 30-50 % of the runs whatever the
        group (identification failures, presence in some individuals only). Without
        the sporadic term they come out as confident false 'off' calls."""
        rng = np.random.default_rng(4)
        n_prot, n = 3000, 8
        mu = rng.normal(24, 1.5, n_prot)                     # far above the 50 % point (19)
        Y = mu[:, None] + rng.normal(0, 0.3, (n_prot, 2 * n))
        pi = np.where(rng.random(n_prot) < 0.08, rng.uniform(0.3, 0.5, n_prot), 0.0)
        Y[rng.random(Y.shape) < pi[:, None]] = np.nan        # same rate in both groups
        off = np.arange(20)
        Y[off, :n] = np.nan                                  # true 'off' in group 0
        Y[off, n:] = mu[off, None] + rng.normal(0, 0.3, (len(off), n))
        Y = Y[(~np.isnan(Y)).sum(axis=1) >= n - 1]
        groups = np.array([0] * n + [1] * n)
        null_rows = np.arange(len(off), len(Y))
        res = {dr: ds.differential_detection(Y, groups, [("B_vs_A", 1, 0)], B0, B1,
                                             dropout=dr) for dr in (None, "auto")}
        fp_old = (res[None]["by_contrast"]["B_vs_A"]["p_adj"][null_rows] < 0.05).sum()
        r = res["auto"]["by_contrast"]["B_vs_A"]
        assert fp_old > 10                                   # the DPC-only model fails here
        assert (r["p_adj"][null_rows] < 0.05).sum() <= 1
        assert (r["p_adj"][:len(off)] < 0.05).mean() >= 0.9  # true on/off still found
        prior = res["auto"]["dropout"]
        assert 0.01 < prior["mean"] < 0.06 and prior["frac_above_0.1"] > 0.04


class TestOrchestration:
    def test_tables_and_main_results_untouched(self, sim):
        keep = filtered(sim)
        Y = sim["Y"][keep]
        X, g = make_design_matrix(sim["groups"], control_condition="Ctrl")
        C, names = make_all_contrasts(g)
        c = names[0]
        rng = np.random.default_rng(0)
        df = pd.DataFrame({"name": [f"P{i}" for i in range(len(Y))],
                           "Genes": [f"G{i}" for i in range(len(Y))],
                           f"{c}_diff": rng.normal(0, 1, len(Y)),
                           f"{c}_p.val": rng.uniform(0, 1, len(Y)),
                           f"{c}_p.adj": rng.uniform(0, 1, len(Y)),
                           f"Robustness_Score_{c}": rng.integers(0, 100, len(Y))})
        before = df.copy()
        Ycopy = Y.copy()
        params = {"volcano_use_padj": False, "volcano_p_thresh": 0.05,
                  "volcano_ratio_min": 1.5, "volcano_lfc_min": float(np.log2(1.5))}
        out = ds.run_missingness_sensitivity(sim["Y"], Y, df, X, C, names, len(g), params,
                                             verbose=False)
        pd.testing.assert_frame_equal(df, before)
        assert np.array_equal(Y, Ycopy, equal_nan=True)
        s = out["dpc_sheet"]
        assert {f"{c}_logFC_dpc", f"{c}_p.val_dpc", f"{c}_p.adj_dpc",
                f"{c}_concordance"} <= set(s.columns)
        assert set(s[f"{c}_concordance"]) <= {"both", "imputation only", "DPC only", ""}
        d = out["detection_sheet"]
        assert {"contrast", "pattern", "higher_in", "p.adj_detection", "verdict"} <= set(d.columns)
        summ = out["summary"].iloc[0]
        assert summ["count_both"] + summ["count_imputation_only"] == summ["count_sig_imputed"]
        assert "robust_median_both" in out["summary"].columns

    def test_numeric_columns_never_look_like_samples(self, sim):
        """Pathfinder picks its PCA matrix among numeric columns that are not
        statistics: every numeric column of the new sheets must carry a token it
        excludes ('_vs_', p.val, p.adj, logfc, rank, count, imputed, mean_intensity...)."""
        keep = filtered(sim)
        Y = sim["Y"][keep]
        X, g = make_design_matrix(sim["groups"], control_condition="Ctrl")
        C, names = make_all_contrasts(g)
        c = names[0]
        df = pd.DataFrame({"name": [f"P{i}" for i in range(len(Y))],
                           f"{c}_diff": np.zeros(len(Y)), f"{c}_p.val": np.ones(len(Y)),
                           f"{c}_p.adj": np.ones(len(Y))})
        params = {"volcano_use_padj": False, "volcano_p_thresh": 0.05,
                  "volcano_ratio_min": 1.5, "volcano_lfc_min": float(np.log2(1.5))}
        out = ds.run_missingness_sensitivity(sim["Y"], Y, df, X, C, names, len(g), params,
                                             verbose=False)
        tokens = ("_vs_", "p.val", "p.adj", "logfc", "rank", "count", "imputed",
                  "mean_intensity", "robust", "score")
        for key in ("dpc_sheet", "detection_sheet", "summary"):
            t = out[key]
            for col in t.columns:
                if pd.api.types.is_numeric_dtype(t[col]):
                    assert any(tok in str(col).lower() for tok in tokens), (key, col)


class TestDeimosStep:
    """run_sensitivity_step / export_sensitivity_sheets in deimos.py."""

    @pytest.fixture
    def deimos(self):
        import deimos
        return deimos

    @pytest.fixture
    def inputs(self, sim):
        keep = filtered(sim)
        Y = sim["Y"][keep]
        X, g = make_design_matrix(sim["groups"], control_condition="Ctrl")
        C, names = make_all_contrasts(g)
        c = names[0]
        df = pd.DataFrame({"name": [f"P{i}" for i in range(len(Y))],
                           f"{c}_diff": np.zeros(len(Y)), f"{c}_p.val": np.ones(len(Y)),
                           f"{c}_p.adj": np.ones(len(Y))})
        model = {"design_mat": X, "contrast_mat": C, "contrast_names": names,
                 "n_groups": len(g)}
        params = {"volcano_use_padj": False, "volcano_p_thresh": 0.05,
                  "volcano_ratio_min": 1.5, "volcano_lfc_min": float(np.log2(1.5))}
        return pd.DataFrame(sim["Y"]), pd.DataFrame(Y), df, model, params

    def test_disabled_does_nothing(self, deimos, inputs, tmp_path):
        all_, filt, df, model, params = inputs
        p = {**params, "dpc_sensitivity": False, "differential_detection": False}
        assert deimos.run_sensitivity_step(all_, filt, df, model, p, str(tmp_path)) is None

    def test_no_missing_value_skips(self, deimos, inputs, tmp_path, capsys):
        all_, filt, df, model, params = inputs
        full = filt.fillna(20.0)
        assert deimos.run_sensitivity_step(all_, full, df, model, params, str(tmp_path)) is None
        assert "[SKIP]" in capsys.readouterr().out

    def test_failure_never_raises(self, deimos, inputs, tmp_path, capsys):
        all_, filt, df, model, params = inputs
        bad = df.iloc[:10]                                   # misaligned on purpose
        assert deimos.run_sensitivity_step(all_, filt, bad, model, params, str(tmp_path)) is None
        assert "[WARN]" in capsys.readouterr().out

    def test_sheets_and_manifest_rows(self, deimos, inputs, tmp_path):
        import openpyxl
        all_, filt, df, model, params = inputs
        res = deimos.run_sensitivity_step(all_, filt, df, model, params, str(tmp_path))
        assert res is not None and res["methods"]["dpc"] and res["methods"]["detection"]
        wb = openpyxl.Workbook()

        def write_df(ws, d):
            ws.append(list(d.columns))
            for r in d.itertuples(index=False):
                ws.append(list(r))
        sheets = deimos.export_sensitivity_sheets(wb, res, write_df, lambda *a, **k: None)
        assert sheets == ["Sensitivity_Summary", "Sensitivity_DPC", "Differential_Detection"]
        assert all(s in wb.sheetnames for s in sheets)
        for s in sheets:       # never captured by the dashboard's GO / WGCNA / UMAP filters
            assert not any(t in s.lower() for t in ("go", "wgcna", "umap"))
        rows = deimos._sensitivity_manifest_rows(res, sheets)
        keys = {k for s, k, v in rows if s == "sensitivity"}
        assert {"summary_sheet", "dpc_sheet", "detection_sheet", "dpc_beta0", "dpc_beta1"} <= keys


def simulate_precursors(seed=5, n_prot=300, n_rep=4, b0=-12.0, b1=0.8):
    """Proteins with 1-8 precursors (additive precursor effect), logistic MNAR."""
    rng = np.random.default_rng(seed)
    n = 2 * n_rep
    K = rng.integers(1, 9, n_prot)
    mu = rng.normal(20, 2, n_prot)
    lfc = np.where(rng.random(n_prot) < 0.2, rng.choice([-1.5, 1.5], n_prot), 0.0)
    rows, pid, true_z = [], [], []
    for p in range(n_prot):
        z = mu[p] + np.r_[np.zeros(n_rep), np.full(n_rep, lfc[p])] + rng.normal(0, 0.3, n)
        true_z.append(z)
        for k in range(K[p]):
            rows.append(z + rng.normal(0, 1.0) + rng.normal(0, 0.2, n))
            pid.append(f"PG{p:04d}")
    Y = np.array(rows)
    Y = np.where(rng.random(Y.shape) < ds.expit(b0 + b1 * Y), Y, np.nan)
    return Y, np.array(pid), np.array(true_z), lfc


class TestPrecursorLevel:
    def test_quantification_tracks_truth_with_honest_se(self):
        Y, pid, true_z, _ = simulate_precursors()
        d = ds.estimate_dpc(Y)
        q = ds.dpc_quant_precursors(Y, pid, d["beta0"], d["beta1"])
        order = np.array([int(p[2:]) for p in q["proteins"]])
        z = true_z[order]
        # abundances are identifiable up to a protein-wise constant: compare
        # centred profiles
        est = q["E"] - q["E"].mean(axis=1, keepdims=True)
        tru = z - z.mean(axis=1, keepdims=True)
        assert np.corrcoef(est.ravel(), tru.ravel())[0, 1] > 0.85
        # more precursors observed -> smaller standard error
        few, many = q["n_obs"] <= 1, q["n_obs"] >= 4
        assert np.median(q["SE"][few]) > np.median(q["SE"][many])

    def test_orchestration_uses_precursors(self):
        Y, pid, _, lfc = simulate_precursors(seed=6)
        # protein-level matrix: mean of observed precursors (stand-in for MaxLFQ)
        prots = np.unique(pid)
        with np.errstate(all="ignore"), warnings.catch_warnings():
            warnings.simplefilter("ignore", RuntimeWarning)
            P = np.array([np.nanmean(Y[pid == p], axis=0) for p in prots])
        keep = (~np.isnan(P)).sum(axis=1) >= 3
        P, prots = P[keep], prots[keep]
        X, g = make_design_matrix(["A"] * 4 + ["B"] * 4)
        C, names = make_all_contrasts(g)
        c = names[0]
        df = pd.DataFrame({"name": prots, "Protein.Group": prots,
                           f"{c}_diff": np.zeros(len(prots)), f"{c}_p.val": np.ones(len(prots)),
                           f"{c}_p.adj": np.ones(len(prots))})
        params = {"volcano_use_padj": True, "volcano_p_thresh": 0.05,
                  "volcano_ratio_min": 1.5, "volcano_lfc_min": float(np.log2(1.5))}
        out = ds.run_missingness_sensitivity(P, P, df, X, C, names, len(g), params,
                                             precursors={"Y": Y, "protein_ids": pid},
                                             verbose=False)
        assert out["dpc_level"] == "precursor"
        assert out["dpc_precursor"]["n_proteins_quantified"] == len(prots)
        sig = out["dpc_sheet"][f"{c}_p.adj_dpc"].to_numpy() < 0.05
        truth = np.isin(prots, np.unique(pid)[lfc != 0][np.isin(np.unique(pid)[lfc != 0], prots)])
        assert (sig & truth).sum() > 0.5 * truth.sum()
        assert (sig & ~truth).sum() <= max(2, 0.1 * sig.sum())
