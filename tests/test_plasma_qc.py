"""Qualité des échantillons de plasma (plasma_qc.py, clé plasma_qc de Deimos).

Trois niveaux, comme pour batch_column : fonctions unitaires sur des cas calculés à
la main ; scénarios synthétiques à contamination connue (tests/plasma_fixture.py) ;
un run Deimos complet qui doit produire les feuilles, la section Methods et le
manifeste — un test qui n'exécute pas la clé ne prouve pas qu'elle agit.
"""
import os
import pathlib
import subprocess
import sys

import numpy as np
import pandas as pd
import pytest

ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "tests"))

import plasma_qc as pq  # noqa: E402
from plasma_fixture import make_plasma_dataset, write_dataset  # noqa: E402


def _log2_matrix(pg, design):
    cols = [c for c in pg.columns if str(c).startswith("C:")]
    m = pg[cols].replace(0, np.nan).astype(float).apply(np.log2)
    m.columns = design["label"].tolist()
    return m


# ─────────────────────────── panels ─────────────────────────────────────────

def test_panels_follow_the_article():
    p = pq.load_panels()
    assert {k: len(v) for k, v in p.items()} == {"erythrocyte": 29, "platelet": 29,
                                                  "coagulation": 31}
    alias = {k: set().union(*v) for k, v in p.items()}
    # exclusions énoncées dans le texte de Geyer et al. 2019
    assert "NIF3L1" not in alias["erythrocyte"]
    assert "F13A1" not in alias["platelet"] and "F13A1" in alias["coagulation"]
    # seuls ACTB et GAPDH sont communs aux panels érythrocytaire et plaquettaire
    assert alias["erythrocyte"] & alias["platelet"] == {"ACTB", "GAPDH"}
    assert {"PF4", "PF4V1"} <= alias["coagulation"]


def test_shared_proteins_are_left_out_of_both_index_sums():
    p = pq.load_panels()
    for panel in ("erythrocyte", "platelet"):
        used = set().union(*pq.index_entries(p, panel))
        assert not used & {"ACTB", "GAPDH"}
        assert len(pq.index_entries(p, panel)) == 27
    assert set().union(*pq.index_entries(p, "coagulation")) == {"FGA", "FGB", "FGG"}


def test_gene_matching_is_case_insensitive_and_knows_renamed_symbols():
    p = pq.load_panels()
    genes = ["Hbb-bs", "CAVIN2", "UBB;UBC", "ALB", ""]
    assert pq.membership(genes, p["erythrocyte"]).tolist() == [True, False, True, False, False]
    assert pq.membership(genes, p["platelet"]).tolist() == [False, True, False, False, False]


@pytest.mark.parametrize("symbol", ["Hba-a1", "Hbb-bs",          # souris
                                    "Hba-a3", "Hba1", "Hbb",       # rat
                                    "HBA", "HBB"])                 # chien (UniProt)
def test_hemoglobins_of_mouse_rat_and_dog_are_recognised(symbol):
    assert pq.membership([symbol], pq.load_panels()["erythrocyte"]).tolist() == [True]


def test_user_panels_file(tmp_path):
    f = tmp_path / "panels.tsv"
    f.write_text("gene\tpanel\nHBB\terythrocyte\nHBA1;HBA2\terythrocyte\nPF4\tplatelet\n")
    p = pq.load_panels(str(f))
    assert len(p["erythrocyte"]) == 2 and frozenset({"HBA1", "HBA2"}) in p["erythrocyte"]
    bad = tmp_path / "bad.tsv"
    bad.write_text("gene\tpanel\nHBB\tliver\n")
    with pytest.raises(ValueError, match="unknown panel"):
        pq.load_panels(str(bad))


# ─────────────────────────── indices et drapeaux ────────────────────────────

def test_indices_are_the_ratios_of_the_article():
    genes = ["HBB", "CA1", "PRDX2", "TLN1", "MYH9", "FLNA", "FGA", "FGB", "FGG", "ALB"]
    lin = np.array([[8, 2], [4, 1], [4, 1], [1, 1], [1, 1], [2, 2],
                    [10, 5], [10, 5], [10, 5], [950, 977]], float)
    m = pd.DataFrame(np.log2(lin), columns=["s1", "s2"])
    idx, info = pq.compute_indices(m, genes, pq.load_panels())
    tot = lin.sum(0)
    np.testing.assert_allclose(idx["erythrocyte"], lin[:3].sum(0) / tot)
    np.testing.assert_allclose(idx["platelet"], lin[3:6].sum(0) / tot)
    np.testing.assert_allclose(idx["coagulation"], tot / lin[6:9].sum(0))
    assert info["erythrocyte"]["n_in_index"] == 3


def test_empty_gene_names_say_nothing_can_be_assessed(capsys):
    # FASTA sans noms de gènes : sous pandas 3 la colonne Genes est toute NaN
    pg, des, _ = make_plasma_dataset("confounded")
    res = pq.run_sample_qc(_log2_matrix(pg, des), [np.nan] * len(pg), des["condition"])
    out = capsys.readouterr().out
    assert "No plasma index could be computed" in out and "[OK]" not in out
    assert not res["flags"]["flagged"].any()


def test_too_few_panel_proteins_gives_no_index():
    m = pd.DataFrame(np.log2([[5.0, 6.0], [100.0, 90.0]]), columns=["a", "b"])
    idx, _ = pq.compute_indices(m, ["HBB", "ALB"], pq.load_panels())
    assert idx["erythrocyte"].isna().all()


def test_mean_plus_2sd_cannot_flag_five_samples_but_the_robust_rule_can():
    assert pq.max_z_reachable(5) < 2 < pq.max_z_reachable(7)
    idx = pd.DataFrame({"erythrocyte": [1e-3, 1.1e-3, 0.9e-3, 1.05e-3, 5e-2],
                        "platelet": 1e-3, "coagulation": 10.0},
                       index=list("abcde"))
    fl = pq.flag_samples(idx)
    assert not fl["erythrocyte_flag_sd"].any()
    assert fl["erythrocyte_flag_robust"].tolist() == [False] * 4 + [True]
    assert fl["flagged"].tolist() == [False] * 4 + [True]


# ─────────────────────────── corrélation intra-condition ────────────────────

def test_within_condition_r_ignores_the_condition_effect():
    rng = np.random.default_rng(1)
    g = np.array(["A"] * 6 + ["B"] * 6)
    x = np.where(g == "B", 3.0, 0.0) + rng.normal(0, 1, 12)   # indice confondu
    y_cond = np.where(g == "B", 2.0, 0.0) + rng.normal(0, 1, 12)  # effet de condition seul
    y_cov = x + rng.normal(0, 0.1, 12)                        # suit l'indice
    raw = np.corrcoef(x, y_cond)[0, 1]
    r, n, df = pq.within_condition_r(np.vstack([y_cond, y_cov]), x, g)
    assert raw > 0.6                    # la corrélation brute se laisse piéger
    assert abs(r[0]) < 0.6 and r[1] > 0.95
    assert n.tolist() == [12, 12] and df.tolist() == [9, 9]


def test_panel_members_are_correlated_with_the_index_computed_without_them():
    # HBB varie seul (bruit propre), les autres membres du panel sont stables : avec
    # HBB dans la somme, l'indice est HBB lui-même (r ~ 1) ; sans lui, r ~ 0.
    rng = np.random.default_rng(3)
    g = np.array(["A"] * 6 + ["B"] * 6)
    genes = ["HBB", "CA1", "PRDX2", "CAT", "SOD1", "ALB"]
    lin = np.vstack([2 ** rng.normal(24, 1.0, 12),
                     *[2 ** (18 + rng.normal(0, 0.05, 12)) for _ in range(4)],
                     np.full(12, 2.0 ** 34)])
    m = pd.DataFrame(np.log2(lin), columns=[f"s{i}" for i in range(12)])
    res = pq.run_sample_qc(m, genes, g, verbose=False)
    corr = pq.protein_correlations(m, genes, res["flags"], g, res["panels"])
    r_full, _, _ = pq.within_condition_r(m.to_numpy()[:1],
                                         res["flags"]["erythrocyte_log2"].to_numpy(), g)
    assert r_full[0] > 0.99
    assert abs(corr.at[0, "r_within_erythrocyte"]) < 0.6


def test_bh_matches_the_textbook_values():
    # m = 4 : 0.01*4/1 = 0.04 ; 0.03*4/2 = 0.06 -> min cumulé 0.0533 ; 0.04*4/3 = 0.0533 ; 0.5
    q = pq.bh([0.01, 0.04, 0.03, np.nan, 0.5])
    np.testing.assert_allclose(q[[0, 1, 2, 4]], [0.04, 0.16 / 3, 0.16 / 3, 0.5], rtol=1e-9)
    assert np.isnan(q[3])


# ─────────────────────────── scénarios à vérité connue ──────────────────────

def test_outlier_scenario_flags_only_the_hemolysed_sample():
    pg, des, truth = make_plasma_dataset("outlier")
    res = pq.run_sample_qc(_log2_matrix(pg, des), pg["Genes"], des["condition"],
                           contrast_names=["A_vs_B"], verbose=False)
    fl = res["flags"]
    assert fl.index[fl["flagged"]].tolist() == ["Sple-A-01"]
    assert "erythrocyte" in fl.at["Sple-A-01", "flag_reason"]
    assert "platelet" not in fl.at["Sple-A-01", "flag_reason"]
    assert res["confounded"] == []
    r = np.corrcoef(fl["erythrocyte_log2"], np.log2(truth["h"]))[0, 1]
    assert r > 0.95


def test_confounded_scenario_flags_no_sample_but_reports_the_confounding():
    pg, des, _ = make_plasma_dataset("confounded")
    res = pq.run_sample_qc(_log2_matrix(pg, des), pg["Genes"], des["condition"],
                           contrast_names=["A_vs_B"], verbose=False)
    assert not res["flags"]["flagged"].any()
    assert res["confounded"] == ["erythrocyte"]
    row = res["confounding"].query("index == 'erythrocyte' and comparison == 'A_vs_B'")
    assert row["delta_log2"].iloc[0] < -2 and row["p_value"].iloc[0] < 0.01


def _welch_dep(m, des, contrast="A_vs_B"):
    import warnings
    from scipy.stats import ttest_ind
    a, b = contrast.split("_vs_")
    A = m.loc[:, (des["condition"] == a).to_numpy()].to_numpy()
    B = m.loc[:, (des["condition"] == b).to_numpy()].to_numpy()
    with warnings.catch_warnings():          # protéines absentes d'un groupe : NaN voulus
        warnings.simplefilter("ignore")
        lfc = np.nanmean(A, 1) - np.nanmean(B, 1)
        p = ttest_ind(A, B, axis=1, equal_var=False, nan_policy="omit").pvalue
    p = np.asarray(p, float)
    res = pd.DataFrame({"name": m.index.astype(str), f"{contrast}_diff": lfc,
                        f"{contrast}_p.val": p})
    return res, {contrast: (p < 0.05) & (np.abs(lfc) >= np.log2(1.5))}


@pytest.mark.parametrize("scenario", ["clean", "confounded"])
def test_true_biological_effects_are_not_flagged(scenario):
    pg, des, truth = make_plasma_dataset(scenario)
    m = _log2_matrix(pg, des)
    genes = pg["Genes"].tolist()
    res = pq.run_sample_qc(m, genes, des["condition"], verbose=False)
    dfr, masks = _welch_dep(m, des)
    res = pq.run_dep_qc(res, m, genes, dfr, ["A_vs_B"], masks, des["condition"].tolist(),
                        verbose=False)
    d = res["dep"]
    eff = d[d["Genes"].isin(truth["effect_genes"])]
    assert len(eff) >= 15
    assert (eff["flag"] == "").mean() >= 0.9
    if scenario == "confounded":
        hi = set(d.loc[d["flag"] == "high", "Genes"])
        assert {"HBB", "HBA1", "CA1"} <= hi
        assert "RBCTRK" in hi        # hors panel : attrapé par la seule corrélation
    else:
        assert not (d["flag"] == "high").any()


# ─────────────────────────── configuration ──────────────────────────────────

def test_keys_are_declared_with_the_same_default_everywhere():
    import config
    import yaml
    assert config._DEFAULTS["plasma_qc"] == "__auto__"
    assert config._DEFAULTS["plasma_qc_panels"] is None
    ex = yaml.safe_load((ROOT / "config_example.yaml").read_text(encoding="utf-8"))
    assert ex["plasma_qc"] == "__auto__" and ex["plasma_qc_panels"] is None


def test_invalid_values_are_refused(tmp_path):
    import config
    errs = config._validate({**config._DEFAULTS, "plasma_qc": "yes"})
    assert any("plasma_qc doit" in e for e in errs)
    for ok in (True, False, "__auto__"):
        assert not any("plasma_qc doit" in e
                       for e in config._validate({**config._DEFAULTS, "plasma_qc": ok}))
    errs = config._validate({**config._DEFAULTS, "plasma_qc_panels": str(tmp_path / "x.tsv")})
    assert any("plasma_qc_panels" in e for e in errs)


def test_step_is_inert_when_switched_off():
    import deimos
    pg, des, _ = make_plasma_dataset("outlier")
    m = _log2_matrix(pg, des)
    meta = pd.DataFrame({"Genes": pg["Genes"]})
    assert deimos.run_plasma_qc_samples(m, meta, des, {"plasma_qc": False}, None) is None


# ─────────────────────────── choix au lancement ─────────────────────────────

def test_plasma_signature_tells_plasma_from_a_lysate():
    pg, des, _ = make_plasma_dataset("clean")
    assert pq.plasma_signature(_log2_matrix(pg, des), pg["Genes"])["looks_like"]
    ex = pd.read_csv(ROOT / "example_data" / "report.pg_matrix.tsv", sep="\t")
    cols = [c for c in ex.columns if str(c).startswith(("C:", "/")) or str(c).endswith(".d")]
    mex = ex[cols].replace(0, np.nan).astype(float).apply(np.log2)
    assert not pq.plasma_signature(mex, ex["Genes"].fillna(""))["looks_like"]
    # lysat dont l'albumine n'est pas en tête
    genes = ["ALB", "APOA1", "TF", "ACTB", "GAPDH", "HSPA8", "ENO1", "PKM", "TUBB", "VIM"]
    m = pd.DataFrame(np.log2(np.array([[30.0], [20], [20], [40], [38], [36], [35], [34],
                                       [33], [32]]) * np.ones((1, 4))))
    assert not pq.plasma_signature(m, genes)["looks_like"]


class _Stdin:
    def __init__(self, tty):
        self._tty = tty

    def isatty(self):
        return self._tty


@pytest.fixture
def plasma_inputs():
    pg, des, _ = make_plasma_dataset("clean")
    return _log2_matrix(pg, des), pd.DataFrame({"Genes": pg["Genes"]})


def test_auto_asks_on_plasma_in_an_interactive_run(monkeypatch, plasma_inputs):
    import deimos
    m, meta = plasma_inputs
    monkeypatch.setattr(deimos.sys, "stdin", _Stdin(True))
    asked = []
    yes = lambda prompt: asked.append(prompt) or ""            # Entrée = oui
    no = lambda prompt: "n"
    assert deimos.resolve_plasma_qc({"plasma_qc": "__auto__"}, m, meta, ask=yes) is True
    assert asked and "plasma panels" in asked[0]
    assert deimos.resolve_plasma_qc({"plasma_qc": "__auto__"}, m, meta, ask=no) is False


def test_auto_never_blocks_a_non_interactive_run(monkeypatch, capsys, plasma_inputs):
    import deimos
    m, meta = plasma_inputs
    monkeypatch.setattr(deimos.sys, "stdin", _Stdin(False))

    def boom(prompt):
        raise AssertionError("asked in a non-interactive run")
    assert deimos.resolve_plasma_qc({"plasma_qc": "__auto__"}, m, meta, ask=boom) is False
    assert "--plasma-qc" in capsys.readouterr().out


def test_auto_is_silent_on_a_lysate(monkeypatch):
    import deimos
    monkeypatch.setattr(deimos.sys, "stdin", _Stdin(True))
    m = pd.DataFrame(np.log2(np.full((3, 4), 1e6)))
    meta = pd.DataFrame({"Genes": ["ACTB", "GAPDH", "ALB"]})

    def boom(prompt):
        raise AssertionError("asked on a lysate")
    assert deimos.resolve_plasma_qc({}, m, meta, ask=boom) is False


def test_explicit_values_win(plasma_inputs):
    import deimos
    m, meta = plasma_inputs
    assert deimos.resolve_plasma_qc({"plasma_qc": True}, m, meta, ask=None) is True
    assert deimos.resolve_plasma_qc({"plasma_qc": False}, m, meta, ask=None) is False


def test_command_line_flag_sets_the_key(tmp_path, monkeypatch):
    import config
    cfg = tmp_path / "c.yaml"
    cfg.write_text(f"out_dir: {tmp_path / 'out'}\ngo_organism: null\n")
    monkeypatch.setattr(sys, "argv", ["deimos.py", "--config", str(cfg), "--plasma-qc"])
    params = config.resolve_config(lambda *a: {}, lambda: None, False, False, "none.tsv")
    assert params["plasma_qc"] is True
    monkeypatch.setattr(sys, "argv", ["deimos.py", "--config", str(cfg)])
    params = config.resolve_config(lambda *a: {}, lambda: None, False, False, "none.tsv")
    assert params["plasma_qc"] == "__auto__"


# ─────────────────────────── run Deimos complet ─────────────────────────────

@pytest.fixture(scope="module")
def deimos_run(tmp_path_factory):
    d = tmp_path_factory.mktemp("plasma_run")
    write_dataset(str(d / "data"), "confounded")
    cfg = d / "cfg.yaml"
    cfg.write_text(f"""tsv_path: {d / 'data' / 'report.pg_matrix.tsv'}
design_path: {d / 'data' / 'ExperimentalDesign.csv'}
out_dir: {d / 'out'}
go_organism: null
run_gsea: false
run_perseverance: false
make_wgcna: false
make_pathfinder: false
make_insight: false
make_dashboard: false
n_iter_robustness: 0
plasma_qc: true
""")
    p = subprocess.run([sys.executable, str(ROOT / "deimos.py"), "--config", str(cfg)],
                       cwd=str(ROOT), capture_output=True, text=True, timeout=900)
    assert p.returncode == 0, p.stdout[-3000:] + p.stderr[-3000:]
    return d / "out", p.stdout


def test_full_run_writes_the_sheets_methods_and_manifest(deimos_run):
    out, log = deimos_run
    wb = out / "ProteomicAnalysis_Results.xlsx"
    xl = pd.ExcelFile(wb)
    assert {"Plasma_QC", "Plasma_QC_DEP"} <= set(xl.sheet_names)
    assert "Erythrocyte lysis index differs between conditions" in log
    d = pd.read_excel(xl, "Plasma_QC_DEP")
    hi = set(d.loc[d["flag"] == "high", "Genes"])
    assert {"HBB", "HBA1", "RBCTRK"} <= hi
    man = pd.read_excel(xl, "_manifest")
    rows = dict(zip(man.loc[man["section"] == "plasma_qc", "key"],
                    man.loc[man["section"] == "plasma_qc", "value"]))
    assert rows["dep_sheet"] == "Plasma_QC_DEP" and rows["confounded_indices"] == "erythrocyte"
    meth = pd.read_excel(xl, "Methods", header=None).astype(str)
    assert meth.apply(lambda c: c.str.contains("PLASMA SAMPLE QUALITY")).any().any()
    from openpyxl import load_workbook
    assert len(load_workbook(wb)["Plasma_QC"]._images) == 1   # figure intégrée
    assert not (out / "plasma_qc_indices.png").exists()        # puis balayée


def test_command_line_on_the_finished_run(deimos_run):
    out, _ = deimos_run
    assert pq.main([str(out)]) == 0
    xl = pd.ExcelFile(out / "plasma_qc" / "Plasma_QC.xlsx")
    assert xl.sheet_names == ["Plasma_QC", "Plasma_QC_DEP"]
    d = pd.read_excel(xl, "Plasma_QC_DEP")
    assert "HBB" in set(d.loc[d["flag"] == "high", "Genes"])
