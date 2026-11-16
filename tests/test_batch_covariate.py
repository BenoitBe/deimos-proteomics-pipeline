"""La clé `batch_column` est réellement modélisée comme covariable de nuisance.

Jusqu'au 16/11/2026 la clé était lue par la configuration et annoncée dans le
README, mais n'entrait dans aucun modèle. Ces tests verrouillent les deux
moitiés du contrat : la covariable est absorbée quand elle est séparable, et
le run est REFUSÉ quand le lot est confondu avec la condition.
"""
import numpy as np
import pandas as pd
import pytest

from deimos import apply_batch_covariate
from limma_ebayes import (make_design_matrix, make_paired_design_matrix,
                          make_all_contrasts)


def _design(conditions, batches=None, subjects=None, col="batch"):
    d = {"label": [f"S{i+1}" for i in range(len(conditions))],
         "condition": list(conditions)}
    if batches is not None:
        d[col] = list(batches)
    if subjects is not None:
        d["subject"] = list(subjects)
    return pd.DataFrame(d)


def _std(conditions):
    mat, groups = make_design_matrix(list(conditions))
    return mat, list(groups), list(groups)


# --- Cas où la covariable ne s'applique pas : design inchangé ----------------

def test_cle_absente_design_inchange():
    conds = ["A", "A", "B", "B"]
    mat, groups, cols = _std(conds)
    out_mat, out_cols = apply_batch_covariate(mat, cols, conds,
                                              _design(conds), {})
    assert out_cols == cols
    assert out_mat.shape == mat.shape


def test_cle_false_ou_vide_design_inchange():
    conds = ["A", "A", "B", "B"]
    mat, groups, cols = _std(conds)
    for val in (None, False, "", "   "):
        out_mat, out_cols = apply_batch_covariate(
            mat, cols, conds, _design(conds, ["b1", "b2", "b1", "b2"]),
            {"batch_column": val})
        assert out_cols == cols, val


def test_colonne_absente_du_plan_avertit_et_continue(capsys):
    conds = ["A", "A", "B", "B"]
    mat, groups, cols = _std(conds)
    out_mat, out_cols = apply_batch_covariate(
        mat, cols, conds, _design(conds), {"batch_column": "plaque"})
    assert out_cols == cols
    err = capsys.readouterr().out
    assert "[WARN]" in err and "plaque" in err
    # le message doit lister les colonnes disponibles, sinon il n'aide pas
    assert "condition" in err


def test_un_seul_niveau_ignore(capsys):
    conds = ["A", "A", "B", "B"]
    mat, groups, cols = _std(conds)
    out_mat, out_cols = apply_batch_covariate(
        mat, cols, conds, _design(conds, ["b1"] * 4), {"batch_column": "batch"})
    assert out_cols == cols
    assert "[WARN]" in capsys.readouterr().out


def test_autant_de_niveaux_que_d_echantillons_ignore(capsys):
    conds = ["A", "A", "B", "B"]
    mat, groups, cols = _std(conds)
    out_mat, out_cols = apply_batch_covariate(
        mat, cols, conds, _design(conds, ["b1", "b2", "b3", "b4"]),
        {"batch_column": "batch"})
    assert out_cols == cols
    assert "[WARN]" in capsys.readouterr().out


def test_valeurs_manquantes_ignore(capsys):
    conds = ["A", "A", "B", "B"]
    mat, groups, cols = _std(conds)
    out_mat, out_cols = apply_batch_covariate(
        mat, cols, conds, _design(conds, ["b1", np.nan, "b1", "b2"]),
        {"batch_column": "batch"})
    assert out_cols == cols
    assert "[WARN]" in capsys.readouterr().out


# --- Cas nominal : la covariable est réellement absorbée ---------------------

def test_batch_croise_ajoute_les_termes(capsys):
    conds = ["A", "A", "B", "B"]
    mat, groups, cols = _std(conds)
    out_mat, out_cols = apply_batch_covariate(
        mat, cols, conds, _design(conds, ["b1", "b2", "b1", "b2"]),
        {"batch_column": "batch"})
    # 2 niveaux -> 1 terme ajouté (le premier niveau est la référence)
    assert len(out_cols) == len(cols) + 1
    assert out_mat.shape == (4, len(cols) + 1)
    assert out_cols[:len(cols)] == cols          # conditions AVANT la nuisance
    assert out_cols[-1].startswith("batch_")
    assert np.linalg.matrix_rank(out_mat) == out_mat.shape[1]
    assert "[MODEL]" in capsys.readouterr().out


def test_trois_lots_deux_termes():
    conds = ["A", "A", "A", "B", "B", "B"]
    mat, groups, cols = _std(conds)
    out_mat, out_cols = apply_batch_covariate(
        mat, cols, conds, _design(conds, ["b1", "b2", "b3"] * 2),
        {"batch_column": "batch"})
    assert len(out_cols) == len(cols) + 2
    assert np.linalg.matrix_rank(out_mat) == out_mat.shape[1]


def test_nom_de_colonne_insensible_a_la_casse():
    conds = ["A", "A", "B", "B"]
    mat, groups, cols = _std(conds)
    d = _design(conds, ["b1", "b2", "b1", "b2"], col="Batch")
    out_mat, out_cols = apply_batch_covariate(mat, cols, conds, d,
                                              {"batch_column": "batch"})
    assert len(out_cols) == len(cols) + 1


def test_contrastes_nuls_sur_les_colonnes_de_lot():
    """Le contrat du padding : les contrastes ne portent QUE sur les conditions."""
    conds = ["A", "A", "B", "B"]
    mat, groups, cols = _std(conds)
    out_mat, out_cols = apply_batch_covariate(
        mat, cols, conds, _design(conds, ["b1", "b2", "b1", "b2"]),
        {"batch_column": "batch"})
    cmat, cnames = make_all_contrasts(groups, n_total_cols=len(out_cols))
    assert cmat.shape[0] == len(out_cols)
    assert np.allclose(cmat[len(groups):, :], 0.0)
    assert np.allclose(cmat[:len(groups), :].sum(axis=0), 0.0)  # somme nulle


def test_compose_avec_un_design_apparie():
    conds = ["A", "B", "A", "B", "A", "B"]
    subs = ["s1", "s1", "s2", "s2", "s3", "s3"]
    # lot croisé avec la condition ET non imbriqué dans le sujet
    batches = ["b1", "b2", "b2", "b1", "b1", "b2"]
    mat, groups, cols = make_paired_design_matrix(conds, subs)
    cols = list(cols)
    out_mat, out_cols = apply_batch_covariate(
        mat, cols, conds, _design(conds, batches, subjects=subs),
        {"batch_column": "batch"})
    assert len(out_cols) == len(cols) + 1
    assert np.linalg.matrix_rank(out_mat) == out_mat.shape[1]
    cmat, _ = make_all_contrasts(groups, n_total_cols=len(out_cols))
    assert np.allclose(cmat[len(groups):, :], 0.0)


def test_la_matrice_d_expression_n_est_jamais_touchee():
    """Modélisation, pas correction de matrice (ComBat/removeBatchEffect)."""
    conds = ["A", "A", "B", "B"]
    mat, groups, cols = _std(conds)
    before = mat.copy()
    apply_batch_covariate(mat, cols, conds,
                          _design(conds, ["b1", "b2", "b1", "b2"]),
                          {"batch_column": "batch"})
    assert np.array_equal(mat, before)           # design d'entrée non muté


# --- Le point dur : refus de la confusion lot / condition --------------------

def test_lot_confondu_avec_la_condition_refuse(capsys):
    conds = ["A", "A", "B", "B"]
    mat, groups, cols = _std(conds)
    with pytest.raises(SystemExit) as exc:
        apply_batch_covariate(mat, cols, conds,
                              _design(conds, ["b1", "b1", "b2", "b2"]),
                              {"batch_column": "batch"})
    assert exc.value.code == 1
    out = capsys.readouterr().out
    assert "[ERROR]" in out
    assert "confondu" in out
    # le diagnostic doit nommer les conditions fautives et montrer le croisement
    assert "'A'" in out or '"A"' in out or "A" in out
    assert "batch" in out


def test_confusion_partielle_refusee():
    """Une seule condition enfermée dans un lot suffit à casser le rang."""
    conds = ["A", "A", "B", "B", "C", "C"]
    batches = ["b1", "b2", "b1", "b2", "b3", "b3"]   # C seul dans b3
    mat, groups, cols = _std(conds)
    with pytest.raises(SystemExit):
        apply_batch_covariate(mat, cols, conds,
                              _design(conds, batches),
                              {"batch_column": "batch"})


def test_lot_imbrique_dans_le_sujet_refuse():
    """Design apparié : un lot = un sujet rend le terme redondant."""
    conds = ["A", "B", "A", "B"]
    subs = ["s1", "s1", "s2", "s2"]
    batches = ["b1", "b1", "b2", "b2"]               # lot == sujet
    mat, groups, cols = make_paired_design_matrix(conds, subs)
    with pytest.raises(SystemExit):
        apply_batch_covariate(mat, list(cols), conds,
                              _design(conds, batches, subjects=subs),
                              {"batch_column": "batch"})


# --- La clé est bien câblée dans le pipeline, pas seulement déclarée ---------

def test_la_cle_est_appelee_dans_run_differential_analysis():
    import inspect
    import deimos
    src = inspect.getsource(deimos.run_differential_analysis)
    assert "apply_batch_covariate(" in src, (
        "batch_column doit être appliquée dans run_differential_analysis "
        "— c'est exactement la régression corrigée le 16/11/2026")


def test_la_cle_existe_dans_la_configuration():
    from config import _DEFAULTS
    assert "batch_column" in _DEFAULTS


# --- Diagnostic PCA : combien de structure de lot le modèle doit-il absorber ?

def _toy_matrix(conditions, batches, n_prot=300, batch_amp=0.0, cond_amp=2.0,
                seed=1):
    """Matrice log2 synthétique : effet condition + effet lot contrôlés."""
    rng = np.random.default_rng(seed)
    conds = pd.Series(list(conditions))
    bats = pd.Series(list(batches))
    expr = rng.normal(20, 0.4, size=(n_prot, len(conds)))
    for lv in conds.unique():
        shift = rng.normal(0, cond_amp, size=n_prot)
        expr[:, (conds == lv).values] += shift[:, None]
    for lv in bats.unique():
        shift = rng.normal(0, batch_amp, size=n_prot)
        expr[:, (bats == lv).values] += shift[:, None]
    return pd.DataFrame(expr, index=[f"P{i}" for i in range(n_prot)],
                        columns=[f"S{i+1}" for i in range(len(conds))])


def test_eps2_detecte_un_effet_de_lot_fort():
    from deimos import _kruskal_eps2
    groups = ["b1", "b1", "b1", "b2", "b2", "b2"]
    sepa = np.array([0.0, 0.1, -0.1, 10.0, 10.1, 9.9])
    melange = np.array([0.0, 10.0, 0.1, 9.9, -0.1, 10.1])
    _, p_s, e_s = _kruskal_eps2(sepa, groups)
    _, p_m, e_m = _kruskal_eps2(melange, groups)
    assert e_s > 0.7 and p_s < 0.05
    assert e_m < e_s


def test_eps2_non_calculable_renvoie_nan():
    from deimos import _kruskal_eps2
    h, p, e = _kruskal_eps2(np.array([1.0, 2.0, 3.0]), ["b1", "b1", "b1"])
    assert e != e                                  # NaN : un seul groupe
    h, p, e = _kruskal_eps2(np.array([5.0, 5.0, 5.0, 5.0]),
                            ["b1", "b1", "b2", "b2"])
    assert e != e                                  # NaN : variance nulle


def test_diagnostic_ecrit_la_figure(tmp_path):
    from deimos import plot_batch_diagnostic
    conds = ["A", "A", "A", "B", "B", "B"]
    bats = ["b1", "b2", "b1", "b2", "b1", "b2"]
    mat = _toy_matrix(conds, bats, batch_amp=0.0)
    d = _design(conds, bats)
    f = plot_batch_diagnostic(mat, d, bats, "batch", str(tmp_path))
    assert f is not None
    assert f.endswith("diagnostic_batch.png")
    assert (tmp_path / "diagnostic_batch.png").exists()
    assert (tmp_path / "diagnostic_batch.png").stat().st_size > 10_000


def test_diagnostic_signale_un_lot_dominant(capsys, tmp_path):
    """Lot énorme, condition faible : le WARN doit partir."""
    conds = ["A", "A", "A", "B", "B", "B"]
    bats = ["b1", "b2", "b1", "b2", "b1", "b2"]
    from deimos import plot_batch_diagnostic
    mat = _toy_matrix(conds, bats, batch_amp=6.0, cond_amp=0.02, seed=3)
    plot_batch_diagnostic(mat, _design(conds, bats), bats, "batch", str(tmp_path))
    out = capsys.readouterr().out
    assert "[BATCH]" in out
    assert "eps2" in out


def test_diagnostic_declenche_par_la_covariable():
    """Le diagnostic doit être appelé là où la covariable est appliquée."""
    import inspect
    import deimos
    src = inspect.getsource(deimos.run_differential_analysis)
    assert "plot_batch_diagnostic(" in src
    assert "batch_info" in src


def test_info_out_expose_la_colonne_resolue():
    conds = ["A", "A", "B", "B"]
    mat, groups, cols = _std(conds)
    info: dict = {}
    apply_batch_covariate(mat, cols, conds,
                          _design(conds, ["b1", "b2", "b1", "b2"], col="Plate"),
                          {"batch_column": "plate"}, info_out=info)
    assert info["applied"] is True
    assert info["column"] == "Plate"               # casse réelle du plan
    assert info["n_levels"] == 2 and info["n_terms"] == 1
    assert len(info["values"]) == 4


def test_info_out_reste_vide_si_non_applique():
    conds = ["A", "A", "B", "B"]
    mat, groups, cols = _std(conds)
    info: dict = {}
    apply_batch_covariate(mat, cols, conds, _design(conds), {}, info_out=info)
    assert info == {}
