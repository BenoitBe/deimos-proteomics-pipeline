"""Recherche du volcano : tous les gènes, et des sélections qui s'accumulent.

Deux limites corrigées le 25/11/2026.

1. `build_proteome_data` sous-échantillonne les non significatifs à
   NS_VOLCANO_MAX par contraste — sans quoi le nuage fige le navigateur. Une
   protéine NS absente du tirage était donc introuvable : elle n'était pas dans
   le fichier. Sur un protéome de 6 000 entrées, cela représentait les deux
   tiers du jeu. `build_volcano_index` porte désormais nom, gène et
   coordonnées de TOUTES les protéines.
2. Chaque recherche effaçait la précédente. Les protéines épinglées
   s'accumulent jusqu'au bouton Reset (logique côté navigateur, vérifiée ici
   sur la présence des points d'entrée).
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

import build_dashboardv7 as bd                                    # noqa: E402

TEMPLATE = (ROOT / "dashboard_template.html").read_text(encoding="utf-8")


def _df(n=4000, n_sig=80, seed=0):
    """Protéome assez grand pour que le plafond NS morde réellement."""
    rng = np.random.default_rng(seed)
    p = rng.uniform(0.2, 1.0, n)
    p[:n_sig] = rng.uniform(0, 1e-5, n_sig)
    lfc = rng.normal(0, 0.1, n)
    lfc[:n_sig] = rng.choice([-2.0, 2.0], n_sig)
    return pd.DataFrame({
        "name":  [f"P{i:05d}" for i in range(n)],
        "Genes": [f"GN{i:05d}" for i in range(n)],
        "A_vs_B_diff":  lfc,
        "A_vs_B_p.val": p,
        "A_vs_B_p.adj": np.clip(p * 10, 0, 1),
    })


# ── L'index couvre tout le protéome ─────────────────────────────────────────

def test_l_index_contient_toutes_les_proteines():
    df = _df()
    idx = bd.build_volcano_index(df, ["A_vs_B"])
    assert len(idx["prot"]) == len(df)
    assert len(idx["xy"]["A_vs_B"]) == 2 * len(df)


def test_l_index_couvre_ce_que_le_nuage_omet():
    """Le cœur du défaut : la part du protéome qui était introuvable."""
    df = _df()
    idx = bd.build_volcano_index(df, ["A_vs_B"])
    dessine = {p["name"] for p in
               bd.build_proteome_data(df, ["A_vs_B"], {})["A_vs_B"]}
    tout = {p[0] for p in idx["prot"]}
    manquantes = tout - dessine
    assert len(manquantes) > 0.5 * len(tout), (
        "le jeu de test doit dépasser le plafond NS, sinon le test ne prouve rien")
    assert manquantes <= tout


def test_les_coordonnees_sont_alignees_sur_la_table():
    df = _df(n=50, n_sig=5)
    idx = bd.build_volcano_index(df, ["A_vs_B"])
    xy = idx["xy"]["A_vs_B"]
    for i, (nom, _gene) in enumerate(idx["prot"]):
        ligne = df.loc[df["name"] == nom].iloc[0]
        assert xy[2 * i] == round(float(ligne["A_vs_B_diff"]), 3)
        attendu = round(float(-np.log10(max(ligne["A_vs_B_p.val"], 1e-300))), 2)
        assert xy[2 * i + 1] == attendu


def test_le_gene_accompagne_le_nom():
    """La recherche porte sur le gène ; sans lui l'index ne sert à rien."""
    df = _df(n=20, n_sig=2)
    idx = bd.build_volcano_index(df, ["A_vs_B"])
    assert all(len(e) == 2 for e in idx["prot"])
    assert idx["prot"][0][1].startswith("GN")


def test_une_valeur_manquante_devient_null_et_non_NaN():
    """NaN n'est pas du JSON : json.dumps l'écrirait et le navigateur casserait."""
    import json
    df = _df(n=30, n_sig=3)
    df.loc[0, "A_vs_B_diff"] = np.nan
    df.loc[1, "A_vs_B_p.val"] = np.nan
    idx = bd.build_volcano_index(df, ["A_vs_B"])
    assert idx["xy"]["A_vs_B"][0] is None
    assert idx["xy"]["A_vs_B"][3] is None
    assert "NaN" not in json.dumps(idx)


def test_un_contraste_sans_colonnes_est_ignore_sans_echec():
    df = _df(n=20, n_sig=2)
    idx = bd.build_volcano_index(df, ["A_vs_B", "C_vs_D"])
    assert "A_vs_B" in idx["xy"] and "C_vs_D" not in idx["xy"]


def test_un_tableau_vide_ne_fait_pas_echouer_le_dashboard():
    idx = bd.build_volcano_index(pd.DataFrame(), ["A_vs_B"])
    assert idx == {"prot": [], "xy": {}}


def test_l_index_est_branche_dans_l_assemblage_et_l_injection():
    import inspect
    src = inspect.getsource(bd)
    assert "data['volc_index'] = build_volcano_index(" in src
    assert "const VOLC_INDEX=" in src


# ── L'interface : épinglage cumulatif et Reset ──────────────────────────────

def test_le_bouton_reset_existe():
    assert 'id="btnVolcReset"' in TEMPLATE
    assert 'onclick="volcResetPins()"' in TEMPLATE


def test_les_fonctions_d_epinglage_existent():
    for fn in ("volcPinCurrent", "volcResetPins", "volcUpdatePinUI",
               "volcIndexHits", "volcIndexXY"):
        assert f"function {fn}(" in TEMPLATE, fn


def test_entree_epingle_et_le_champ_se_vide():
    assert "onkeydown" in TEMPLATE and "volcPinCurrent();" in TEMPLATE
    pin = TEMPLATE[TEMPLATE.index("function volcPinCurrent("):]
    pin = pin[:pin.index("\nfunction ")]
    assert "box.value = '';" in pin, "le champ doit se vider pour enchaîner"


def test_une_epingle_reste_jaune_apres_une_autre_recherche():
    """Le point dur : le surlignage ne doit plus dépendre du seul terme courant."""
    assert "||!!volcPins[p.name]" in TEMPLATE


def test_le_dessin_utilise_l_index():
    assert "volcIndexHits(search).forEach(addIdx)" in TEMPLATE
    assert "canvas._pts=pts.concat(extra)" in TEMPLATE


def test_le_reset_vide_les_epingles_et_redessine():
    reset = TEMPLATE[TEMPLATE.index("function volcResetPins("):]
    reset = reset[:reset.index("\nfunction ")]
    assert "volcPins = {}" in reset
    assert "drawVolcano()" in reset


def test_la_legende_ne_peut_plus_diverger():
    """Une surcharge au rendu réécrivait la légende : les deux doivent coïncider."""
    legende = ("non-significant points subsampled for display; "
               "search, counts, DEP table and exports cover every protein")
    assert legende in TEMPLATE
    assert legende in (ROOT / "build_dashboardv7.py").read_text(encoding="utf-8")
