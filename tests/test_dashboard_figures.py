"""Figures exploratoires du dashboard (09/11/2026) : raincloud observé / imputé, motifs de
valeurs manquantes (UpSet), identifiants iPath3, panneaux et script injectés une fois,
ancres présentes dans le modèle et les autres patchs, JavaScript syntaxiquement valide."""
import importlib.util
import json
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
    spec = importlib.util.spec_from_file_location("bd_figs", ROOT / "build_dashboardv7.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def _design(n_prot=300, seed=1):
    """2 conditions x 3 réplicats ; intensités brutes (log2, NaN = absente) et matrice des
    tests (NaN remplacés par des valeurs basses, comme une imputation à gauche)."""
    rng = np.random.default_rng(seed)
    cols = ["A_1", "A_2", "A_3", "B_1", "B_2", "B_3"]
    raw = pd.DataFrame(rng.normal(20, 2, (n_prot, 6)), columns=cols)
    raw.iloc[:20, 0] = np.nan                  # 20 manquantes dans A_1
    raw.iloc[20:30, 3:6] = np.nan              # 10 absentes de toute la condition B
    pg = [f"P{i}" for i in range(n_prot)]
    comp = pd.concat([pd.DataFrame({"name": pg, "Protein.Group": pg}), raw], axis=1)
    imp = raw.copy()
    imp[imp.isna()] = 12.0
    imp.index = pg
    return comp, imp, cols, {c: c.rsplit("_", 1)[0] for c in cols}


def test_raincloud_flags_imputed_values(bd):
    comp, imp, cols, cmap = _design()
    out = bd.build_raincloud(comp, imp, cols, cmap)
    assert out["mask"] is True and out["n_prot"] == 300 and out["order"] == cols
    a1, b1, a2 = out["samples"]["A_1"], out["samples"]["B_1"], out["samples"]["A_2"]
    assert (a1["n_obs"], a1["n_imp"]) == (280, 20) and (b1["n_obs"], b1["n_imp"]) == (290, 10)
    assert a2["n_imp"] == 0 and a2["di"] == [] and a2["pi"] == []
    assert a1["mi"] == 12.0 and a1["cond"] == "A"
    # densités en nombres : l'aire des observées vaut leur nombre (à la troncature près)
    g = np.array(out["grid"])
    assert abs(np.trapezoid(a1["do"], g) - 280) / 280 < 0.03
    # pluie : tirage conjoint, valeurs imputées = 12 et part réelle préservée
    assert set(a1["pi"]) <= {12.0} and len(a1["po"]) + len(a1["pi"]) <= 220
    assert len(a1["pi"]) <= 0.2 * (len(a1["po"]) + len(a1["pi"]))
    assert a1["q"][0] < a1["q"][2] < a1["q"][4]


def test_raincloud_without_raw_intensities_or_matrix(bd):
    comp, imp, cols, cmap = _design()
    out = bd.build_raincloud(comp[["name", "Protein.Group"]], imp, cols, cmap)
    assert out["mask"] is False and out["samples"]["A_1"]["n_imp"] == 0
    assert bd.build_raincloud(comp, None, cols, cmap) is None


def test_missing_patterns_classes(bd, tmp_path):
    rng = np.random.default_rng(2)
    labels = ["s1", "s2", "s3", "s4", "s5", "s6"]          # libellés d'origine
    v = pd.DataFrame(rng.uniform(1e5, 1e7, (100, 6)), columns=[f"LFQ.intensity.{l}" for l in labels])
    v.iloc[0:8, 3:6] = np.nan          # condition B entière -> 'cond'
    v.iloc[8:20, 1] = np.nan           # un échantillon -> 'single'
    v.iloc[20:25, [0, 1]] = np.nan     # deux échantillons de A -> 'within'
    v.iloc[25:28, [0, 4]] = np.nan     # A et B -> 'mixed'
    v.iloc[28, :] = 0                  # 0 = absent aussi -> 'all'
    raw = pd.concat([pd.DataFrame({"Protein.Group": [f"P{i}" for i in range(100)],
                                   "Genes": [f"G{i}" for i in range(100)]}), v], axis=1)
    path = tmp_path / "s.xlsx"
    raw.to_excel(path, sheet_name="raw_data", index=False)
    cols = ["A_1", "A_2", "A_3", "B_1", "B_2", "B_3"]
    bd.SAMPLE_LABELS = dict(zip(cols, labels))
    out = bd.build_missing_patterns(str(path), cols, {c: c[0] for c in cols})
    assert out["samples"] == cols and out["labels"] == labels and out["cond"] == list("AAABBB")
    assert out["class_tot"] == {"cond": 8, "single": 12, "within": 5, "mixed": 3, "all": 1}
    assert out["n_with_na"] == 29 and out["n_prot"] == 100
    assert out["set_size"] == [9, 18, 1, 9, 12, 9]
    top = out["patterns"][0]
    assert top["k"] == "single" and top["m"] == [1] and top["n"] == 12 and top["g"][0] == "G8"
    b = next(p for p in out["patterns"] if p["k"] == "cond")
    assert b["m"] == [3, 4, 5] and b["n"] == 8
    # sans libellés : repli sur l'ordre du plan
    bd.SAMPLE_LABELS = {}
    assert bd.build_missing_patterns(str(path), cols, {c: c[0] for c in cols})["class_tot"]["cond"] == 8
    assert bd.build_missing_patterns(str(tmp_path / "absent.xlsx"), cols, {}) is None
    # plans larges : motifs par condition entièrement absente, classes toujours par échantillon
    bd.MVP_MAX_SAMPLES = 4
    out = bd.build_missing_patterns(str(path), cols, {c: c[0] for c in cols})
    assert out["level"] == "condition" and out["samples"] == ["A", "B"] and out["n_samples"] == 6
    assert [(p["m"], p["n"], p["k"]) for p in out["patterns"]] == [([1], 8, "cond"), ([0, 1], 1, "all")]
    assert out["set_size"] == [1, 9] and out["class_tot"]["single"] == 12


def test_ipath_ids_prefer_reference_ortholog(bd, tmp_path):
    path = tmp_path / "o.xlsx"
    pd.DataFrame({"protein_id": ["A0X1"], "ortholog_id": ["P45952"], "pident": [83.0],
                  "coverage": [100], "confidence": ["high"]}).to_excel(
        path, sheet_name="Orthologues_Perseverance", index=False)
    comp = pd.DataFrame({"name": ["A0X1", "A0X2", "A0X3"],
                         "Protein.Group": ["A0X1", "A0X2;A0X9", "A0X3-2;A0X4"]})
    out = bd.build_ipath_ids(str(path), comp)
    assert out["ids"] == {"A0X1": "P45952", "A0X3": "A0X3-2"} and out["n_orth"] == 1 and out["n"] == 3
    assert bd.build_ipath_ids(str(tmp_path / "absent.xlsx"), comp)["ids"] == {"A0X3": "A0X3-2"}


def test_anchors_exist_in_template_and_patches(bd):
    tpl = (ROOT / "dashboard_template.html").read_text(encoding="utf-8")
    assert '<canvas id="missingCanvas" style="display:block;width:100%;"></canvas>\n  </div>\n</div>' in bd.QC_TAB_HTML
    assert ('<div id="gseaNesMsg" style="display:none;color:var(--tx2);font-size:12px;padding:18px 4px;"></div>\n  </div>'
            in bd.GSEA_PAGE_HTML)
    assert '<span style="font-size:10px;color:var(--tx2);" id="pcaLoadInfo"></span>' in bd.PCA_LOAD_HTML
    for a in ('<button class="exp-btn" onclick="exportVolcAll()">📋 All proteins</button>',
              "<a onclick=\"showPage('volcano')\">Volcano</a>", '<div id="heatmap" class="page">'):
        assert a in tpl
    # la page GSEA expose son état au ridge plot et l'appelle à chaque rendu
    assert "window.gseaState = function" in bd.GSEA_JS and "window.drawGseaRidges()" in bd.GSEA_JS


def test_patch_injects_once_and_js_is_valid(bd, tmp_path):
    html = ("<nav><a onclick=\"showPage('volcano')\">Volcano</a></nav>"
            '<button class="exp-btn" onclick="exportVolcAll()">📋 All proteins</button>'
            + bd.QC_TAB_HTML + bd.PCA_LOAD_HTML + bd.GSEA_PAGE_HTML
            + '<div id="heatmap" class="page"></div><!--PROTEOGEN_BODY_END--></body>')
    data = {"raincloud": {"grid": [1, 2], "samples": {}, "order": [], "mask": False, "n_prot": 0},
            "mv_patterns": {"labels": ["a</script>b"]}, "ipath_ids": {"ids": {}, "n_orth": 0, "n": 0}}
    out = bd.patch_figures_2026_11(html, data)
    for k in ('id="rainPanel"', 'id="mvpPanel"', 'id="gseaRidgePanel"', 'id="pcaBiplot"',
              'id="btnIPath"', 'id="fcpage"', "showPage('fcpage')"):
        assert out.count(k) == 1, k
    assert out.index('id="missingCanvas"') < out.index('id="rainPanel"')
    assert out.index('id="gseaNesCanvas"') < out.index('id="gseaRidgePanel"') < out.index('id="gseaTable"')
    assert out.index('id="pcaBiplot"') < out.index('id="pcaLoadInfo"')
    assert out.index('id="fcpage"') < out.index('id="heatmap"')
    assert out.index("Figures exploratoires (09/11/2026)") < out.index("<!--PROTEOGEN_BODY_END-->")
    assert "a</script>b" not in out                               # chaîne neutralisée
    node = shutil.which("node")
    if node is None:
        pytest.skip("node absent : syntaxe JS non vérifiée")
    i = out.index("const RAIN_DATA=")
    data_js = out[i:out.index("</script>", i)]
    js = re.search(r"<script>(.*?)</script>", bd.FIG2611_JS, re.S).group(1)
    for name, src in (("data.js", data_js), ("fig.js", js)):
        f = tmp_path / name
        f.write_text(src, encoding="utf-8")
        r = subprocess.run([node, "--check", str(f)], capture_output=True, text=True)
        assert r.returncode == 0, r.stderr
    assert json.loads(data_js.split("const MVP_DATA=")[1].split(";\n")[0].replace("<\\/", "</")) == data["mv_patterns"]
