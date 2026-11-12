"""ANOVA par strate : une ANOVA et ses heatmaps par sous-ensemble de conditions."""
import numpy as np
import pandas as pd
import pytest

from deimos import parse_anova_strata, run_anova_strata

CONDS = ["T0_Dark_Col_0", "T0_Dark_OMP_like_Col6D", "T0_Dark_OMP_like_Col5F",
         "T3_Dark_Col_0", "T3_Dark_OMP_like_Col6D", "T3_Dark_OMP_like_Col5F"]


def _design(conds=CONDS, n_rep=4):
    rows = [dict(label=f"S{i}", condition=c, replicate=r + 1)
            for i, (c, r) in enumerate((c, r) for c in conds for r in range(n_rep))]
    return pd.DataFrame(rows)


def _data(design, n_prot=60, seed=0):
    rng = np.random.default_rng(seed)
    meta = pd.DataFrame({"name": [f"P{i}" for i in range(n_prot)]})
    x = rng.normal(12, 0.3, (n_prot, len(design)))
    # un effet génotype franc, présent À T3 SEULEMENT
    t3 = (design.condition.str.startswith("T3") & design.condition.str.contains("Col6D")).values
    x[:12, t3] += 4.0
    return pd.DataFrame(x, columns=design.label), meta


def test_regex_and_mapping_give_the_same_strata():
    d = _design()
    by_rx = parse_anova_strata(d, r"^(T\d+)_")
    assert set(by_rx) == {"T0", "T3"}
    assert by_rx["T0"] == CONDS[:3] and by_rx["T3"] == CONDS[3:]
    named = parse_anova_strata(d, r"^(?P<stratum>T\d+)_")
    assert named == by_rx
    assert parse_anova_strata(d, {"T0": CONDS[:3], "T3": CONDS[3:]}) == by_rx
    for off in (False, None, "false", ""):
        assert parse_anova_strata(d, off) == {}


def test_raw_and_normalised_condition_names_both_work(capsys):
    """design garde les noms bruts du CSV (T0-Dark-Col-0) alors que les contrastes
    portent les noms normalisés (T0_Dark_Col_0) : les deux écritures doivent marcher,
    sinon la regex échoue en silence."""
    raw = ["T0-Dark-Col-0", "T0-Dark-OMP-like-Col6D", "T3-Dark-Col-0",
           "T3-Dark-OMP-like-Col6D"]
    d = _design(conds=raw, n_rep=2)
    for rx in (r"^(T\d+)_", r"^(T\d+)-"):                # écriture normalisée OU brute
        out = parse_anova_strata(d, rx)
        assert set(out) == {"T0", "T3"}, rx
        assert out["T0"] == raw[:2]                      # rendus dans la forme de design
    # un mapping écrit avec les noms normalisés retrouve les conditions brutes
    m = parse_anova_strata(d, {"T0": ["T0_Dark_Col_0", "T0_Dark_OMP_like_Col6D"]})
    assert m == {"T0": raw[:2]}
    # une regex qui ne colle à rien le dit, au lieu de ne rien faire en silence
    assert parse_anova_strata(d, r"^ZZZ(\d+)") == {}
    log = capsys.readouterr().out
    assert "matches none of the 4 condition names" in log


def test_bad_strata_are_dropped_with_a_warning(capsys):
    d = _design()
    out = parse_anova_strata(d, {"ok": CONDS[:3], "seul": [CONDS[0]],
                                 "fantome": ["NEXISTE_PAS"], "partiel": [CONDS[0], "X"]})
    assert set(out) == {"ok"}                       # 1 condition ou 0 -> écartée
    log = capsys.readouterr().out
    assert "unknown condition" in log and "skipped" in log
    assert parse_anova_strata(d, "([") == {}        # regex invalide, pas d'exception


def test_one_anova_per_stratum_finds_the_local_effect(tmp_path):
    d = _design()
    mat, meta = _data(d)
    params = {"anova_use_padj": True, "anova_p_thresh": 0.05, "n_heatmap_clusters": 3,
              "anova_strata": r"^(T\d+)_"}
    res = run_anova_strata(mat, meta, d, params, str(tmp_path))
    assert set(res) == {"T0", "T3"}
    # 3 groupes x 4 réplicats = 12 échantillons, 9 df résiduels dans chaque strate
    for r in res.values():
        assert r["n_samples"] == 12 and r["df_resid"] == 9
        assert list(r["table"]["stratum"].unique()) == [r["table"]["stratum"].iloc[0]]
    # l'effet n'existe qu'à T3 : il doit ressortir là, et pas à T0
    assert res["T3"]["n_sig"] >= 12 > res["T0"]["n_sig"]
    sig_t3 = set(res["T3"]["table"].query("significant")["name"])
    assert {f"P{i}" for i in range(12)} <= sig_t3
    # BH propre à chaque strate : la correction ne voit que ses propres tests
    t3 = res["T3"]["table"]
    assert (t3["p.adj"] >= t3["p.value"] - 1e-12).all()
    assert len(t3) == 60                            # toutes les protéines rapportées
    # les figures sont écrites dans leur sous-dossier
    assert (tmp_path / "anova_strata" / "T3" / "heatmap_annotated.png").exists()


def test_a_stratum_without_residual_df_is_skipped(tmp_path, capsys):
    d = _design(conds=CONDS[:3], n_rep=1)           # 3 groupes, 3 échantillons, 0 df
    mat, meta = _data(d)
    params = {"anova_use_padj": False, "anova_p_thresh": 0.05, "n_heatmap_clusters": 3,
              "anova_strata": r"^(T\d+)_"}
    assert run_anova_strata(mat, meta, d, params, str(tmp_path)) == {}
    assert "residual degree" in capsys.readouterr().out


def test_disabled_by_default_and_global_anova_untouched(tmp_path):
    d = _design()
    mat, meta = _data(d)
    base = {"anova_use_padj": True, "anova_p_thresh": 0.05, "n_heatmap_clusters": 3}
    assert run_anova_strata(mat, meta, d, base, str(tmp_path)) == {}
    assert not (tmp_path / "anova_strata").exists()


def test_dashboard_counts_imputed_whatever_the_dtype():
    """pandas >= 2.1 (infer_string) et pandas 3 donnent le dtype `str`, pas `object`,
    à une colonne de texte : le test `dtype == object` du dashboard tombait dans la
    branche numérique et affichait 0 % d'imputation sur des VRAI/FAUX."""
    import build_dashboardv7 as B

    def n_imputed(col):
        df = pd.DataFrame({"imputed": col, "num_NAs": [1, 0, 2]})
        rep = B.build_qc_report(None, df, ["s1", "s2", "s3"], {})
        sec = next(s for s in rep["sections"] if s["title"] == "Imputation")
        m = next(m for m in sec["metrics"] if m["label"] == "Proteines avec imputation")
        return int(m["detail"].split("/")[0])

    assert n_imputed(pd.Series(["VRAI", "FAUX", "VRAI"], dtype="str")) == 2
    assert n_imputed(pd.Series(["VRAI", "FAUX", "VRAI"], dtype=object)) == 2
    assert n_imputed(pd.Series(["TRUE", "FALSE", "true"], dtype=object)) == 2
    assert n_imputed(pd.Series([True, False, True])) == 2
    assert n_imputed(pd.Series([1, 0, 1])) == 2


def test_missing_value_status_scales_with_the_design():
    """2,83 NA/protéine, c'est 6 % sur 47 échantillons (bon) et 47 % sur 6 (mauvais) :
    le statut suit la part de manquants, le nombre absolu reste affiché."""
    import build_dashboardv7 as B

    def metric(n_samples, mean_na):
        df = pd.DataFrame({"imputed": ["FAUX"] * 4, "num_NAs": [mean_na] * 4})
        rep = B.build_qc_report(None, df, [f"s{i}" for i in range(n_samples)], {})
        sec = next(s for s in rep["sections"] if s["title"] == "Imputation")
        return next(m for m in sec["metrics"] if m["label"] == "NA moyen par proteine")

    big = metric(47, 2.83)
    assert big["status"] == "OK" and big["value"] == 2.83      # 6 % : correct
    assert "6.0% des valeurs" in big["detail"] and "47 echantillons" in big["detail"]
    small = metric(6, 2.83)
    assert small["status"] == "FAIL" and small["value"] == 2.83  # 47 % : mauvais
    assert metric(20, 2.83)["status"] == "WARN"                  # 14 %


def test_pathfinder_skipped_steps_are_relayed(tmp_path):
    """Pathfinder continue son rapport quand une étape échoue et l'écrit dans son log.
    Sans relais, un onglet entier disparaît sans un mot dans la sortie de Deimos."""
    from pathfinder_bridge import _skipped_steps
    log = tmp_path / "pathfinder_log.txt"
    log.write_text(
        "$ python pathfinder.py --stats x.xlsx\n"
        "[Pathfinder] reading stats\n"
        "[Pathfinder] ⚠ Biological reading skipped (ImportError: no module 'reading')\n"
        "[Pathfinder] ⚠ figure 'volcano' failed\n"
        "[Pathfinder] ⚠ Biological reading skipped (ImportError: no module 'reading')\n"
        "[Pathfinder] report written\n", encoding="utf-8")
    out = _skipped_steps(str(log))
    assert len(out) == 2                                  # doublon fusionné
    assert any("Biological reading skipped" in l and "ImportError" in l for l in out)
    assert not any(l.startswith("[Pathfinder]") for l in out)
    assert _skipped_steps(str(tmp_path / "absent.txt")) == []
    quiet = tmp_path / "ok.txt"
    quiet.write_text("[Pathfinder] all good\n", encoding="utf-8")
    assert _skipped_steps(str(quiet)) == []


def test_strata_are_suggested_from_the_condition_names():
    """La question interactive doit marcher sur n'importe quel projet : les découpages
    sont déduits des noms de conditions, sans rien coder en dur."""
    from deimos import suggest_strata

    def labels(conds):
        return {p["label"]: p["values"] for p in suggest_strata(conds)}

    # plan réel : temps x génotype, nombre de jetons variable selon la condition
    ara = [f"T{t}-Dark-{g}" for t in (0, 3, 5, 10)
           for g in ("Col-0", "OMP-like-Col6D", "OMP-like-Col5F")]
    out = labels(ara)
    assert out["1st segment from the left"] == ["T0", "T3", "T5", "T10"]
    assert len(out["1st segment from the right"]) == 3           # les 3 génotypes
    # les découpages déséquilibrés (OMP/Col, nom coupé au mauvais endroit) sont écartés
    assert all(len({len(v) for v in p["strata"].values()}) == 1 for p in suggest_strata(ara))

    # plan 2 x 3 avec un autre séparateur
    out = labels(["A_1h", "A_6h", "A_24h", "B_1h", "B_6h", "B_24h"])
    assert ["1h", "6h", "24h"] in out.values() and ["A", "B"] in out.values()

    # pas de structure factorielle : on ne propose rien plutôt que n'importe quoi
    assert suggest_strata(["WT", "KO", "WT2", "KO2"]) == []
    assert suggest_strata(["A", "B"]) == []                      # trop peu de conditions
    # un niveau présent une seule fois ne fait pas une strate
    assert all("uniq" not in p["values"] for p in suggest_strata(ara + ["uniq-Dark-Col-0"]))


def test_design_conditions_are_read_whatever_the_separator(tmp_path):
    from deimos import read_design_conditions
    for sep in (";", ",", "\t"):
        f = tmp_path / f"d{ord(sep)}.csv"
        f.write_text(f"label{sep}condition{sep}replicate\n"
                     f"S1{sep}T0-Col{sep}1\nS2{sep}T0-Col{sep}2\nS3{sep}T3-Col{sep}1\n")
        assert read_design_conditions(str(f)) == ["T0-Col", "T3-Col"]
    assert read_design_conditions(str(tmp_path / "absent.csv")) == []
    bad = tmp_path / "bad.csv"
    bad.write_text("a;b\n1;2\n")                                 # pas de colonne condition
    assert read_design_conditions(str(bad)) == []
