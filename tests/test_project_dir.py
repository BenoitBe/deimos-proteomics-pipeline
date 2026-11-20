"""`project_dir` : un dossier par projet, au lieu de tout empiler chez Deimos.

Sans la clé, un chemin relatif part du répertoire courant : les matrices, le
plan et le FASTA finissent à côté de deimos.py et `proteomics_output/` est
écrasé au projet suivant. Avec elle, les chemins relatifs du projet sont
résolus dans son dossier — et les chemins absolus, qui désignent des
références partagées entre projets, ne bougent pas.
"""
import os

import pytest

from config import PROJECT_PATH_KEYS, apply_project_dir, _DEFAULTS


# ── Sans la clé : rien ne change ────────────────────────────────────────────

def test_sans_cle_les_chemins_sont_intacts():
    p = {"tsv_path": "report.pg_matrix.tsv", "out_dir": "proteomics_output"}
    avant = dict(p)
    assert apply_project_dir(p) == avant


@pytest.mark.parametrize("vide", [None, "", "   ", False])
def test_valeur_vide_equivaut_a_pas_de_cle(vide):
    p = {"project_dir": vide, "tsv_path": "a.tsv", "out_dir": "out"}
    apply_project_dir(p)
    assert p["tsv_path"] == "a.tsv" and p["out_dir"] == "out"


def test_la_cle_existe_dans_les_defauts_et_vaut_None():
    assert "project_dir" in _DEFAULTS
    assert _DEFAULTS["project_dir"] is None      # comportement historique par défaut


# ── Résolution ──────────────────────────────────────────────────────────────

def test_les_chemins_relatifs_sont_resolus_dans_le_projet(tmp_path, capsys):
    p = {"project_dir": str(tmp_path),
         "tsv_path": "report.pg_matrix.tsv",
         "design_path": "ExperimentalDesign.csv",
         "out_dir": "output"}
    apply_project_dir(p)
    assert p["tsv_path"] == str(tmp_path / "report.pg_matrix.tsv")
    assert p["design_path"] == str(tmp_path / "ExperimentalDesign.csv")
    assert p["out_dir"] == str(tmp_path / "output")
    out = capsys.readouterr().out
    assert "[PROJECT]" in out and str(tmp_path) in out


def test_un_chemin_absolu_n_est_jamais_touche(tmp_path):
    """Une référence partagée entre projets reste hors du dossier de projet."""
    partage = tmp_path / "references" / "uniprot_all.fasta"
    partage.parent.mkdir()
    partage.write_text(">x", encoding="utf-8")
    p = {"project_dir": str(tmp_path),
         "tsv_path": "report.pg_matrix.tsv",
         "string_reference_fasta": str(partage),
         "out_dir": "/var/resultats"}
    apply_project_dir(p)
    assert p["string_reference_fasta"] == str(partage)
    assert p["out_dir"] == "/var/resultats"          # absolu : il gagne
    assert p["tsv_path"].startswith(str(tmp_path))   # relatif : résolu


def test_auto_reste_auto(tmp_path):
    """`annotation_fasta: __auto__` est un marqueur, pas un chemin."""
    p = {"project_dir": str(tmp_path), "annotation_fasta": "__auto__"}
    apply_project_dir(p)
    assert p["annotation_fasta"] == "__auto__"


def test_les_cles_non_renseignees_restent_nulles(tmp_path):
    p = {"project_dir": str(tmp_path), "perseverance_fasta_a": None,
         "gsea_go_obo": "", "tsv_path": "a.tsv"}
    apply_project_dir(p)
    assert p["perseverance_fasta_a"] is None
    assert p["gsea_go_obo"] == ""


def test_tilde_et_variables_d_environnement(tmp_path, monkeypatch):
    monkeypatch.setenv("MON_PROJET", str(tmp_path))
    p = {"project_dir": "$MON_PROJET", "tsv_path": "a.tsv"}
    apply_project_dir(p)
    assert p["project_dir"] == str(tmp_path)
    assert p["tsv_path"] == str(tmp_path / "a.tsv")


def test_sous_dossier_relatif(tmp_path):
    (tmp_path / "donnees").mkdir()
    p = {"project_dir": str(tmp_path), "tsv_path": "donnees/report.pg_matrix.tsv"}
    apply_project_dir(p)
    assert p["tsv_path"] == str(tmp_path / "donnees" / "report.pg_matrix.tsv")


def test_project_dir_est_rendu_absolu(tmp_path, monkeypatch):
    (tmp_path / "projet_x").mkdir()
    monkeypatch.chdir(tmp_path)
    p = {"project_dir": "projet_x", "tsv_path": "a.tsv"}
    apply_project_dir(p)
    assert os.path.isabs(p["project_dir"])
    assert p["project_dir"] == str(tmp_path / "projet_x")


def test_appliquer_deux_fois_ne_double_pas_le_chemin(tmp_path):
    """Les chemins sont absolus après le premier passage, donc intouchables."""
    p = {"project_dir": str(tmp_path), "tsv_path": "a.tsv", "out_dir": "out"}
    apply_project_dir(p)
    attendu = dict(p)
    apply_project_dir(p)
    assert p == attendu


# ── Refus net plutôt qu'un « fichier introuvable » obscur ───────────────────

def test_un_dossier_de_projet_inexistant_arrete_le_run(tmp_path, capsys):
    p = {"project_dir": str(tmp_path / "nexiste_pas"), "tsv_path": "a.tsv"}
    with pytest.raises(SystemExit) as exc:
        apply_project_dir(p)
    assert exc.value.code == 1
    out = capsys.readouterr().out
    assert "[ERROR]" in out and "project_dir" in out


def test_un_fichier_n_est_pas_un_dossier_de_projet(tmp_path):
    f = tmp_path / "pas_un_dossier.txt"
    f.write_text("x", encoding="utf-8")
    with pytest.raises(SystemExit):
        apply_project_dir({"project_dir": str(f), "tsv_path": "a.tsv"})


# ── Le périmètre des clés ───────────────────────────────────────────────────

def test_les_clones_d_outils_et_le_cache_sont_hors_perimetre():
    """Les y rattacher obligerait à dupliquer un clone de Pathfinder par projet."""
    for k in ("pathfinder_path", "insight_path", "gsea_cache_dir"):
        assert k not in PROJECT_PATH_KEYS, k


def test_toutes_les_cles_du_perimetre_existent_dans_les_defauts():
    """Une clé mal orthographiée ici serait silencieusement sans effet."""
    for k in PROJECT_PATH_KEYS:
        assert k in _DEFAULTS, k


def test_les_entrees_et_la_sortie_du_projet_sont_couvertes():
    for k in ("tsv_path", "design_path", "pr_path", "out_dir"):
        assert k in PROJECT_PATH_KEYS, k


def test_les_trois_modes_de_configuration_resolvent(tmp_path):
    """--config, rechargement de last_config, et questions interactives."""
    import inspect
    import config
    src = inspect.getsource(config.resolve_config)
    assert src.count("apply_project_dir(") == 3, (
        "les trois chemins de resolve_config doivent résoudre project_dir")


# ── Toute clé de chemin doit être classée, sinon le build tombe ─────────────

def test_toute_cle_de_chemin_est_classee():
    """Une clé de chemin ajoutée plus tard et oubliée resterait relative au
    répertoire courant pendant que le reste du projet a déménagé — une erreur
    muette, et exactement celle que `plasma_qc_panels` a failli introduire le
    20/11/2026 en arrivant sur main pendant ce développement."""
    import re
    from config import SHARED_PATH_KEYS

    motif = re.compile(r"(path|_dir|fasta|obo|panels|gmt|_file)$")
    non_classees = [k for k in _DEFAULTS
                    if motif.search(k)
                    and k not in PROJECT_PATH_KEYS
                    and k not in SHARED_PATH_KEYS]
    assert not non_classees, (
        "clé(s) de chemin non classée(s) — ajoutez-les à PROJECT_PATH_KEYS "
        "(propres au projet) ou à SHARED_PATH_KEYS (partagées entre projets) :\n"
        + "\n".join(f"  {k}" for k in sorted(non_classees)))


def test_les_deux_listes_sont_disjointes():
    from config import SHARED_PATH_KEYS
    assert not set(PROJECT_PATH_KEYS) & set(SHARED_PATH_KEYS)


def test_un_panel_plasma_personnalise_suit_le_projet(tmp_path):
    p = {"project_dir": str(tmp_path), "plasma_qc_panels": "mes_panels.tsv"}
    apply_project_dir(p)
    assert p["plasma_qc_panels"] == str(tmp_path / "mes_panels.tsv")
