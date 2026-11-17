"""Nom du dashboard HTML configurable, sans pouvoir sortir de out_dir."""
import pytest

from deimos import (DASHBOARD_DEFAULT_NAME, DASHBOARD_LEGACY_NAME,
                    dashboard_filename)


def test_the_default_is_the_new_name():
    """Renommage demandé le 12/11/2026 : le défaut n'est plus proteogen_*."""
    for params in ({}, None, {"dashboard_filename": None}, {"dashboard_filename": ""}):
        assert dashboard_filename(params) == "Proteomic_dashboard.html"
    assert DASHBOARD_DEFAULT_NAME == "Proteomic_dashboard.html"


def test_the_old_name_is_still_reachable():
    """Un projet qui tient à l'ancien nom le remet dans sa config."""
    assert DASHBOARD_LEGACY_NAME == "proteogen_dashboard.html"
    assert dashboard_filename({"dashboard_filename": DASHBOARD_LEGACY_NAME}) \
        == "proteogen_dashboard.html"


def test_a_plain_name_is_taken_as_is():
    assert dashboard_filename({"dashboard_filename": "Proteomic_dashboard.html"}) \
        == "Proteomic_dashboard.html"
    assert dashboard_filename({"dashboard_filename": "resultats.htm"}) == "resultats.htm"


def test_the_html_extension_is_added_when_missing():
    assert dashboard_filename({"dashboard_filename": "Proteomic_dashboard"}) \
        == "Proteomic_dashboard.html"
    assert dashboard_filename({"dashboard_filename": "ostrea.v2"}) == "ostrea.v2.html"


def test_a_path_is_reduced_to_its_file_name():
    """Le dashboard doit rester dans out_dir : une valeur de config ne peut pas
    le faire écrire ailleurs, ni remonter l'arborescence."""
    for value in ("C:\\\\Windows\\\\System32\\\\evil.html", "/etc/evil.html",
                  "../../evil.html", "sous/dossier/evil.html",
                  "..\\\\..\\\\evil.html"):
        out = dashboard_filename({"dashboard_filename": value})
        assert out == "evil.html", value
        assert "/" not in out and "\\\\" not in out and ".." not in out


def test_characters_windows_refuses_are_replaced():
    out = dashboard_filename({"dashboard_filename": 'tab:1?<x>|y*.html'})
    assert out == "tab_1__x__y_.html"
    assert dashboard_filename({"dashboard_filename": "   "}) == "Proteomic_dashboard.html"
    assert dashboard_filename({"dashboard_filename": "..."}) == "Proteomic_dashboard.html"


def test_the_config_default_matches_the_code_default():
    """config.py et deimos.py doivent annoncer le même nom : c'est l'écart entre
    les deux qui a fait sortir un proteogen_dashboard.html alors que la mise à
    jour promettait l'autre nom."""
    from config import _DEFAULTS
    assert _DEFAULTS["dashboard_filename"] == DASHBOARD_DEFAULT_NAME


def test_the_cleanup_keeps_the_renamed_dashboard():
    """Le ménage de fin de run ne doit pas supprimer un dashboard renommé,
    ni celui qu'un run antérieur avait écrit sous l'ancien nom.

    Vérifié en EXÉCUTANT la règle de conservation, pas en lisant le source :
    la version précédente de ce test cherchait `keep_exact = {` dans
    `deimos.py` et serait passée sur n'importe quel bloc contenant les bons
    mots, même inatteignable.
    """
    from deimos import figures_to_keep, DASHBOARD_LEGACY_NAME

    _, exact = figures_to_keep({})
    assert "Proteomic_dashboard.html" in exact        # défaut courant
    assert DASHBOARD_LEGACY_NAME in exact             # écrit par un run antérieur

    # un nom choisi par l'utilisateur doit être épargné lui aussi
    _, exact = figures_to_keep({"dashboard_filename": "mon_rapport.html"})
    assert "mon_rapport.html" in exact
    assert "Proteomic_dashboard.html" in exact
    assert DASHBOARD_LEGACY_NAME in exact
