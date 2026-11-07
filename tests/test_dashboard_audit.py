"""Audit du dashboard (06/11/2026) : le dashboard montre et exporte ce que contient le
classeur. GO appariés par le _manifest, moyennes sans manquants mis à 0, p non
arrondies à 0, table complète FULL_DE, q-value de la heatmap = p.adj de l'ANOVA, PCA de
Deimos, réécriture p.Genes -> p.g et script de correction syntaxiquement valide."""
import importlib.util
import re
import shutil
import subprocess
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[1]


@pytest.fixture()
def bd():
    spec = importlib.util.spec_from_file_location("bd_audit", ROOT / "build_dashboardv7.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


CONTRASTS = ["KO_vs_WT_stress", "KO_vs_KO_stress", "KO_vs_WT"]


def test_go_sheet_follows_manifest(bd, tmp_path):
    # sans manifeste, le préfixe « KOvWT » colle au premier contraste qui commence ainsi
    assert bd.map_go_sheet("GO_KOvWT", CONTRASTS) == "KO_vs_WT_stress"
    xl = tmp_path / "w.xlsx"
    man = pd.DataFrame([["section", "key", "value"], ["go_sheet", "GO_KOvWT", "KO_vs_WT"],
                        ["go_sheet", "GO_KOvKO_stres", "KO_vs_KO_stress"]])
    with pd.ExcelWriter(xl) as w:
        man.to_excel(w, sheet_name="_manifest", header=False, index=False)
    bd._GO_MANIFEST = bd.read_manifest_go(str(xl))
    assert bd._GO_MANIFEST == {"GO_KOvWT": "KO_vs_WT", "GO_KOvKO_stres": "KO_vs_KO_stress"}
    assert bd.map_go_sheet("GO_KOvWT", CONTRASTS) == "KO_vs_WT"
    assert bd.read_manifest_go(str(tmp_path / "absent.xlsx")) == {}


def _comp():
    return pd.DataFrame({
        "name": ["P1", "P2", "P3"], "Protein.Group": ["P1", "P2", "P3"],
        "Genes": ["G1", "G2", ""], "First.Protein.Description": ["d1", "d2", "d3"],
        "KO_1": [28.0, 20.0, np.nan], "KO_2": [np.nan, 20.0, np.nan], "WT_1": [27.0, 21.0, 25.0],
        "KO_vs_WT_diff": [1.0, -1.0, 0.1], "KO_vs_WT_p.val": [3.2e-9, 0.01, 0.7],
        "KO_vs_WT_p.adj": [1.0e-8, 0.02, 0.7], "Pi_Score_KO_vs_WT": [8.0, 2.0, 0.0],
        "Robustness_Score_KO_vs_WT": [20, 12, 0]})


def test_means_ignore_missing_and_prefer_analysis_matrix(bd):
    df = _comp()
    cond_map = {"KO": ["KO_1", "KO_2"], "WT": ["WT_1"]}
    pts = {p["name"]: p for p in bd.build_proteome_data(df, ["KO_vs_WT"], cond_map)["KO_vs_WT"]}
    assert pts["P1"]["x"] == 28.0            # (28 + manquant) / 1, et non (28 + 0) / 2 = 14
    assert pts["P3"]["x"] is None            # jamais vue dans KO
    assert pts["P1"]["p"] == pytest.approx(3.2e-9) and pts["P1"]["p"] > 0   # plus d'arrondi à 0
    imp = pd.DataFrame({"KO_1": [28.0, 20.0, 18.0], "KO_2": [26.0, 20.0, 18.4], "WT_1": [27.0, 21.0, 25.0]},
                       index=["P1", "P2", "P3"])
    pts = {p["name"]: p for p in bd.build_proteome_data(df, ["KO_vs_WT"], cond_map, imp=imp)["KO_vs_WT"]}
    assert pts["P1"]["x"] == 27.0 and pts["P3"]["x"] == 18.2


def test_full_table_has_every_protein(bd):
    full = bd.build_full_de(_comp(), ["KO_vs_WT"])
    assert full["n"] == ["P1", "P2", "P3"] and full["g"][:2] == ["G1", "G2"]
    c = full["c"]["KO_vs_WT"]
    assert c["d"] == [1.0, -1.0, 0.1] and c["r"] == [20, 12, 0]
    assert c["p"][0] == pytest.approx(3.2e-9) and c["q"][1] == 0.02


def test_heatmap_q_is_anova_padj(bd):
    z = pd.DataFrame({"Protein.Group": ["P1", "P2"], "S_1": [1.0, -1.0], "S_2": [-1.0, 1.0]})
    an = pd.DataFrame({"Protein.Group": ["P1", "P2"], "p.adj": [0.001, 0.04],
                       "significant": [True, True], "Cluster_ID": ["Cluster_1", "Cluster_2"]})
    rows = bd.build_heatmap(z, an, ["S_1", "S_2"])
    assert [r["q"] for r in rows] == [0.001, 0.04]


def test_pca_panel_is_the_deimos_pca(bd, tmp_path):
    xl = tmp_path / "s.xlsx"
    sc = pd.DataFrame({"label": ["Sple-KO-01", "Sple-WT-01"], "condition": ["KO", "WT"],
                       "PC1": [12.5, -12.5], "PC2": [1.0, -1.0], "PC3": [0.2, -0.2],
                       "PC1_var_pct": 40.0, "PC2_var_pct": 20.0, "PC3_var_pct": 10.0})
    with pd.ExcelWriter(xl) as w:
        sc.to_excel(w, sheet_name="PCA_Scores", index=False)
    bd.SAMPLE_LABELS = {"KO_1": "Sple-KO-01", "WT_1": "Sple-WT-01"}
    out = bd.build_pca(pd.DataFrame(), ["KO_1", "WT_1"], {"KO_1": "KO", "WT_1": "WT"},
                       stats_path=str(xl))
    assert out["source"].startswith("Deimos") and out["var_exp"] == [40.0, 20.0, 10.0]
    p = {q["sample"]: q for q in out["points"]}
    assert p["KO_1"]["pc1"] == 12.5 and p["KO_1"]["cond"] == "KO" and p["KO_1"]["label"] == "Sple-KO-01"
    # repli sans PCA_Scores : protéines incomplètes écartées (plus de NaN -> 0)
    df = pd.DataFrame({"A_1": [1.0, 2.0, np.nan, 4.0], "B_1": [2.0, 1.0, 5.0, 3.0]})
    out = bd.build_pca(df, ["A_1", "B_1"], {"A_1": "A", "B_1": "B"})
    assert len(out["points"]) == 2 and "complete" in out["source"]


def test_audit_patch_rewrites_genes_and_injects_valid_js(bd, tmp_path):
    html = ("<script>function drawVolcano(){const isU=sig&&p.d>fc,isD=sig&&p.d<-fc;"
            "[-3,-2,-1,0,1,2,3].filter(v=>v>=xmin&&v<=xmax).forEach(v=>{});"
            "ctx.font='bold 12px sans-serif';ctx.textAlign='left';}\n"
            "function f(){var a=(p.Genes||'');var b=best.p.Genes||best.p.name;}</script>"
            "<!--PROTEOGEN_BODY_END--></body>")
    out = bd.patch_audit_2026_11(html)
    assert "((p.g||p.Genes)||'')" in out and "(best.p.g||best.p.Genes)" in out
    assert "const isU=sig&&p.d>0,isD=sig&&p.d<0;" in out and "fullCounts(activeContrast" in out
    assert "[-3,-2,-1,0,1,2,3]" not in out
    assert out.index("Audit du 06/11/2026") < out.index("<!--PROTEOGEN_BODY_END-->")
    # page GO : un contraste sans terme affiche un message, pas la vue précédente
    for fn in ("window.renderGO = function", "window.drawGOBubble = function",
               "window.renderGOTable = function", "window.drawGONetwork = function"):
        assert fn in bd.AUDIT_FIX_JS
    # bouton « Inverser » : les exports du volcano suivent l'affichage
    assert "volcOrient(fullRows(c))" in bd.AUDIT_FIX_JS and "volcName()" in bd.AUDIT_FIX_JS
    node = shutil.which("node")
    if node is None:
        pytest.skip("node absent : syntaxe JS non vérifiée")
    js = re.search(r"<script>(.*?)</script>", bd.AUDIT_FIX_JS, re.S).group(1)
    f = tmp_path / "fix.js"
    f.write_text(js, encoding="utf-8")
    assert subprocess.run([node, "--check", str(f)], capture_output=True).returncode == 0
