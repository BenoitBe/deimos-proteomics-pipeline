"""primary_statistic: dpc (statistique principale sans imputation, méthode limpa),
robustesse leave-one-sample-out, règle MIXED corrigée, appariement des runs du
pr_matrix. Données simulées, rapide."""
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

PARAMS = {"volcano_use_padj": True, "volcano_p_thresh": 0.05,
          "volcano_ratio_min": 1.5, "volcano_lfc_min": float(np.log2(1.5)),
          "n_iter_robustness": 10, "fdr_global": False}


def simulate_precursors(seed=5, n_prot=240, n_rep=4, b0=-12.0, b1=0.8):
    """Protéines à 1-8 précurseurs (effet précurseur additif), MNAR logistique."""
    rng = np.random.default_rng(seed)
    n = 2 * n_rep
    K = rng.integers(1, 9, n_prot)
    mu = rng.normal(20, 2, n_prot)
    lfc = np.where(rng.random(n_prot) < 0.2, rng.choice([-1.5, 1.5], n_prot), 0.0)
    rows, pid = [], []
    for p in range(n_prot):
        z = mu[p] + np.r_[np.zeros(n_rep), np.full(n_rep, lfc[p])] + rng.normal(0, 0.3, n)
        for _ in range(K[p]):
            rows.append(z + rng.normal(0, 1.0) + rng.normal(0, 0.2, n))
            pid.append(f"PG{p:04d}")
    Y = np.array(rows)
    Y = np.where(rng.random(Y.shape) < ds.expit(b0 + b1 * Y), Y, np.nan)
    return Y, np.array(pid), lfc


@pytest.fixture(scope="module")
def data():
    Y, pid, lfc = simulate_precursors()
    prots = np.unique(pid)
    with np.errstate(all="ignore"), warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        P = np.array([np.nanmean(Y[pid == p], axis=0) for p in prots])
    keep = (~np.isnan(P)).sum(axis=1) >= 3
    X, g = make_design_matrix(["A"] * 4 + ["B"] * 4)
    C, names = make_all_contrasts(g)
    truth = lfc[keep] != 0
    return {"Y": Y, "pid": pid, "P_all": P, "P": P[keep], "prots": prots[keep],
            "X": X, "C": C, "names": names, "n_groups": len(g), "truth": truth}


@pytest.fixture(scope="module")
def primary(data):
    # deux protéines retirées du pr_matrix : repli au niveau protéine
    drop = np.isin(data["pid"], data["prots"][:2])
    return ds.dpc_primary_analysis(
        data["P_all"], data["P"], data["prots"], data["X"], data["C"], verbose=False,
        precursors={"Y": data["Y"][~drop], "protein_ids": data["pid"][~drop]})


class TestPrimaryAnalysis:
    def test_sources_and_fallback(self, data, primary):
        assert primary["level"] == "precursor"
        assert list(primary["source"][:2]) == ["protein", "protein"]
        assert (primary["source"][2:] == "precursor").all()
        assert np.isfinite(primary["p_value"]).all()          # aucune protéine perdue
        assert np.isnan(primary["n_precursors"][:2]).all()

    def test_bh_is_common_to_all_proteins(self, primary):
        p = primary["p_value"][:, 0]
        assert np.allclose(primary["p_adj"][:, 0], ds._bh_correction(p))

    def test_precursor_block_identical_to_dpc_de(self, primary):
        b = primary["_blocks"]["precursor"]
        r = ds.dpc_de(b["E"], b["SE"], b["imputed"], primary["_design"], primary["_contrast"])
        assert np.allclose(primary["p_value"][b["rows"], 0], r["p_value"][:, 0])
        assert np.allclose(primary["logFC"][b["rows"], 0], r["logFC"][:, 0])

    def test_calibrated_and_powerful(self, data, primary):
        sig = primary["p_adj"][:, 0] < 0.05
        assert (sig & data["truth"]).sum() > 0.5 * data["truth"].sum()
        assert (sig & ~data["truth"]).sum() <= max(2, 0.1 * sig.sum())

    def test_protein_level_without_precursors(self, data):
        r = ds.dpc_primary_analysis(data["P_all"], data["P"], None, data["X"], data["C"],
                                    verbose=False)
        assert r["level"] == "protein" and r["dpc_precursor"] is None
        assert (r["source"] == "protein").all()


class TestLooRobustness:
    def test_scores(self, data, primary):
        rob = ds.dpc_loo_robustness(primary, data["n_groups"], PARAMS, verbose=False)
        s = rob["score"][:, 0]
        assert rob["n_fits"][0] == 8 and rob["n_skipped"] == 0
        assert np.nanmin(s) >= 0 and np.nanmax(s) <= 100
        sig = (primary["p_adj"][:, 0] < 0.05) & (np.abs(primary["logFC"][:, 0]) >= np.log2(1.5))
        assert np.median(s[sig]) >= 75            # appels solides : robustes au retrait
        assert np.median(s[~sig]) == 0

    def test_single_sample_group_is_skipped(self, data):
        X, g = make_design_matrix(["A"] * 4 + ["B"] * 3 + ["C"])
        C, _ = make_all_contrasts(g)
        r = ds.dpc_primary_analysis(data["P_all"], data["P"], None, X, C, verbose=False)
        rob = ds.dpc_loo_robustness(r, len(g), PARAMS, verbose=False)
        assert rob["n_skipped"] == 1               # retirer le seul C rend le design non estimable


@pytest.fixture(scope="module")
def deimos():
    import deimos
    return deimos


class TestApplyPrimary:
    def test_columns_switched_nothing_lost(self, deimos, data, primary):
        c = data["names"][0]
        n = len(data["P"])
        df = pd.DataFrame({"name": data["prots"], "Protein.Group": data["prots"],
                           f"{c}_diff": np.full(n, 9.0), f"{c}_p.val": np.full(n, 0.5),
                           f"{c}_p.adj": np.full(n, 0.6), f"{c}_p.val_limma": np.full(n, 0.7),
                           f"{c}_p.adj_limma": np.full(n, 0.8), f"Pi_Score_{c}": np.full(n, 1.0),
                           f"Robustness_Score_{c}": np.full(n, 20)})
        rob = {"score": np.full((n, 1), 50.0)}
        out = deimos.apply_dpc_primary(df, {**primary, "robustness": rob}, data["names"])
        assert np.allclose(out[f"{c}_p.val"], primary["p_value"][:, 0])
        assert np.allclose(out[f"{c}_diff"], primary["logFC"][:, 0])
        assert (out[f"Robustness_Score_{c}"] == 50).all()
        assert (out[f"{c}_diff_imputed"] == 9).all() and (out[f"{c}_p.adj_imputed"] == 0.6).all()
        assert (out[f"{c}_p.adj_imputed_limma"] == 0.8).all()
        assert (out[f"Robustness_Score_imputed_{c}"] == 20).all()
        assert not any(col.endswith("_p.adj_limma") for col in out.columns)
        assert len(out.columns) == len(df.columns) + 2 + 5      # 2 méta + 5 colonnes DPC
        # le dashboard ne détecte que le vrai contraste
        assert [x for x in out.columns if x.endswith("_p.val")] == [f"{c}_p.val"]


class TestMixedRule:
    def _mat(self):
        rng = np.random.default_rng(0)
        X = rng.normal(20, 2, (200, 6))
        X[rng.random(X.shape) < 0.1] = np.nan
        X[0, :3] = np.nan
        X[0, 5] = np.nan
        X[0, 3:5] = [20.2, 18.9]
        cols = [f"S{i}" for i in range(6)]
        return pd.DataFrame(X, columns=cols), pd.DataFrame(
            {"label": cols, "condition": ["A"] * 3 + ["B"] * 3})

    def test_whole_absence_row_imputed_low(self, deimos):
        mat, dsg = self._mat()
        new = deimos.impute(mat, "mixed", design=dsg, rng=np.random.default_rng(1))
        old = deimos.impute(mat, "mixed_legacy", design=dsg, rng=np.random.default_rng(1))
        q = deimos.impute_qrilc(mat, rng=np.random.default_rng(1))
        assert np.isclose(new.values[0, 5], q.values[0, 5])          # QRILC, comme le groupe A
        assert not np.isclose(old.values[0, 5], q.values[0, 5])      # kNN dans l'ancienne règle
        na = mat.isna().values
        whole = na[:, :3].all(1) | na[:, 3:].all(1)
        assert np.allclose(new.values[~whole], old.values[~whole])
        assert not new.isna().values.any()


class TestPrecursorColumns:
    def test_injection_match_and_ambiguity(self, deimos, tmp_path):
        pg = tmp_path / "pg.tsv"
        meta = ["Protein.Group", "Protein.Names", "Genes", "First.Protein.Description",
                "N.Sequences", "N.Proteotypic.Sequences"]
        runs_pg = [r"E:\d\x-Sple-1_Slot1-1_1_100.d", r"E:\d\x-Sple-2_Slot1-2_1_101.d",
                   r"E:\d\x-Sple-3_Slot1-3_1_102.d"]
        pg.write_text("\t".join(meta + runs_pg) + "\n")
        # dans le pr_matrix, l'injection 102 porte encore l'ancien label (Sple-2, doublon)
        runs_pr = [r"E:\d\x-Sple-1_Slot1-1_1_100.d", r"E:\d\x-Sple-2_Slot1-2_1_101.d",
                   r"E:\d\x-Sple-2_Slot1-3_1_102.d"]
        m, notes = deimos._precursor_columns(runs_pr, ["Sple-1", "Sple-2", "Sple-3"], str(pg))
        assert m == {"Sple-1": runs_pr[0], "Sple-2": runs_pr[1], "Sple-3": runs_pr[2]}
        assert len(notes) == 1 and notes[0].endswith("-> Sple-3")
        m2, notes2 = deimos._precursor_columns(runs_pr, ["Sple-1", "Sple-2", "Sple-3"], None)
        assert "Sple-2" not in m2 and "Sple-3" not in m2          # jamais choisi au hasard
        assert any(n.startswith("AMBIGUOUS Sple-2") for n in notes2)


class TestPrimaryStep:
    def test_fallback_without_pr_matrix(self, deimos, data, capsys):
        c = data["names"][0]
        df = pd.DataFrame({"name": data["prots"], "Protein.Group": data["prots"],
                           f"{c}_diff": 0.0, f"{c}_p.val": 1.0, f"{c}_p.adj": 1.0})
        model = {"design_mat": data["X"], "contrast_mat": data["C"],
                 "contrast_names": data["names"], "n_groups": data["n_groups"]}
        params = dict(PARAMS, primary_statistic="dpc")
        r = deimos.run_dpc_primary_step(pd.DataFrame(data["P_all"]), pd.DataFrame(data["P"]),
                                        df, model, params, "/nonexistent/pr.tsv")
        assert r is None and "no usable pr_matrix" in params["_dpc_primary_fallback"]
        assert "stays the imputed analysis" in capsys.readouterr().out

    def test_imputation_requested(self, deimos, data):
        params = dict(PARAMS, primary_statistic="imputation")
        assert deimos.run_dpc_primary_step(None, None, None, {}, params, None) is None

    def test_runs_with_precursor_file(self, deimos, data, tmp_path):
        labels = [f"Sple-{i + 1}" for i in range(8)]
        runs = [f"x-{l}_Slot1-{i}_1_{500 + i}.d" for i, l in enumerate(labels)]
        pr = pd.DataFrame(np.where(np.isnan(data["Y"]), 0.0, 2.0 ** data["Y"]), columns=runs)
        pr.insert(0, "Protein.Group", data["pid"])
        pr.insert(1, "Protein.Names", data["pid"])
        pr.to_csv(tmp_path / "pr.tsv", sep="\t", index=False)
        c = data["names"][0]
        df = pd.DataFrame({"name": data["prots"], "Protein.Group": data["prots"],
                           f"{c}_diff": 0.0, f"{c}_p.val": 1.0, f"{c}_p.adj": 1.0})
        model = {"design_mat": data["X"], "contrast_mat": data["C"],
                 "contrast_names": data["names"], "n_groups": data["n_groups"]}
        filt = pd.DataFrame(data["P"], columns=labels)
        r = deimos.run_dpc_primary_step(pd.DataFrame(data["P_all"]), filt, df, model,
                                        dict(PARAMS, primary_statistic="dpc"),
                                        str(tmp_path / "pr.tsv"))
        assert r is not None and r["level"] == "precursor"
        assert r["robustness"]["score"].shape == (len(data["P"]), 1)


class TestConfig:
    def test_keys_validated(self):
        import config
        base = dict(config._DEFAULTS, tsv_path="x", design_path="y")
        assert not [e for e in config._validate(dict(base, primary_statistic="imputation",
                                                     impute_method="mixed_legacy"))
                    if "primary_statistic" in e or "impute_method" in e]
        errs = config._validate(dict(base, primary_statistic="limpa"))
        assert any("primary_statistic" in e for e in errs)
        assert config._DEFAULTS["primary_statistic"] == "dpc"
