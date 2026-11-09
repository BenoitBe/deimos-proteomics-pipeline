"""Page GSEA du dashboard (07/11/2026) : onglets GSEA_* lus d'après le _manifest (section
gsea_sheet, schéma 1.3), table de gseapy repérée par sa colonne Term, leading edge rendu
en protéines (colonne Lead_proteins, sinon symbole -> accession), noms d'onglets uniques
à l'export, page et script injectés une fois et syntaxiquement valides."""
import importlib.util
import re
import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np
import openpyxl
import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
import gsea_enrichment as gs  # noqa: E402


@pytest.fixture()
def bd():
    spec = importlib.util.spec_from_file_location("bd_gsea", ROOT / "build_dashboardv7.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def _table(n=4, lead_prot=True):
    t = pd.DataFrame({
        "Name": "prerank",
        "Term": ["GO:BP__GO:0009414 response to water deprivation",
                 "KEGG__KEGG:00010 Glycolysis",
                 "GO_Biological_Process_2023__translation (GO:0006412)",
                 "WP__custom set"][:n],
        "ES": [0.7, -0.6, 0.5, -0.4][:n], "NES": [2.1, -1.9, 1.5, -1.2][:n],
        "NOM p-val": [0.0, 0.001, 0.01, 0.02][:n], "FDR q-val": [0.0, 0.01, 0.05, 0.2][:n],
        "FWER p-val": [0.0, 0.02, 0.1, 0.4][:n], "Tag %": ["5/20"] * n, "Gene %": ["10%"] * n,
        "Lead_genes": ["G1;G2", "G3", "G4;G9", "G5"][:n]})
    if lead_prot:
        t["Lead_proteins"] = ["P1;P2", "P3", "P4", "P5"][:n]
    return t


def _workbook(tmp_path, contrasts, manifest=True, lead_prot=True):
    wb = openpyxl.Workbook()
    wb.active.title = "Differential_Expression"
    res = {c: {"table": _table(lead_prot=lead_prot),
               "provenance": "GSEA pre-ranked (gseapy), 500 permutations, FDR <= 0.25"}
           for c in contrasts}
    sm = gs.export_gsea_sheets(wb, res, lambda ws, img, cell: None)
    if manifest:
        ws = wb.create_sheet("_manifest")
        ws.append(("section", "key", "value"))
        for sh, c in sm.items():
            ws.append(("gsea_sheet", sh, c))
    p = tmp_path / "w.xlsx"
    wb.save(p)
    return str(p), sm


def test_export_names_are_unique_and_returned(tmp_path):
    # deux contrastes dont le nom d'onglet tronqué à 31 caractères est identique
    a = "Genotype_Tron_Light_T3_vs_Genotype_Tron_Dark_T3"
    b = "Genotype_Tron_Light_T3_vs_Genotype_Tron_Dark_T7"
    path, sm = _workbook(tmp_path, [a, b])
    assert sorted(sm.values()) == [a, b] and len(set(k.lower() for k in sm)) == 2
    assert all(len(k) <= 31 for k in sm)
    assert openpyxl.load_workbook(path)[next(iter(sm))]["A1"].value.startswith("GSEA pre-ranked")


def test_build_gsea_reads_manifest_table_and_lead_proteins(bd, tmp_path):
    a = "Genotype_Tron_Light_T3_vs_Genotype_Tron_Dark_T3"
    b = "Genotype_Tron_Light_T3_vs_Genotype_Tron_Dark_T7"
    path, sm = _workbook(tmp_path, [a, b])
    bd._GSEA_MANIFEST = bd.read_manifest_go(path, "gsea_sheet")
    assert bd._GSEA_MANIFEST == sm
    out = bd.build_gsea(path, [a, b], pd.DataFrame({"name": ["P1"], "Genes": ["G1"]}))
    assert set(out) == {a, b}
    d = out[a]
    assert d["fdr_max"] == 0.25 and d["n_perm"] == 500 and d["sheet"] in sm
    t = {x["id"]: x for x in d["terms"]}
    assert t["GO:0009414"]["n"] == "response to water deprivation" and t["GO:0009414"]["s"] == "GO:BP"
    assert t["GO:0006412"]["n"] == "translation"                    # identifiant en fin de nom
    assert t["KEGG:00010"]["nes"] == -1.9 and t["KEGG:00010"]["q"] == 0.01
    assert t["GO:0009414"]["lp"] == ["P1", "P2"] and t["GO:0009414"]["q"] == 0.0
    assert d["terms"][3]["n"] == "custom set" and d["terms"][3]["id"] == ""


def test_lead_proteins_fallback_and_name_mapping_without_manifest(bd, tmp_path):
    path, sm = _workbook(tmp_path, ["KO_vs_WT"], manifest=False, lead_prot=False)
    bd._GSEA_MANIFEST = {}
    comp = pd.DataFrame({"name": ["P1", "P2", "P9"], "Genes": ["g1", "G2;G2b", "G9"]})
    out = bd.build_gsea(path, ["KO_vs_WT", "KO_vs_WT_stress"], comp)
    assert list(out) == ["KO_vs_WT"]                                 # GSEA_KOvWT -> KO_vs_WT
    t = out["KO_vs_WT"]["terms"]
    assert t[0]["lp"] == ["P1", "P2"] and t[2]["lp"] == ["P9"]
    assert bd.build_gsea(None, ["KO_vs_WT"], comp) == {}


def test_run_gsea_adds_lead_proteins(tmp_path, monkeypatch):
    gseapy = pytest.importorskip("gseapy")
    df = pd.DataFrame({"name": ["P1", "P2", "P3"], "Genes": ["a", "B;b2", "c"]})
    k2p = gs._key_proteins(df)
    assert k2p == {"A": ["P1"], "B": ["P2"], "C": ["P3"]}
    assert gs._lead_proteins("A;C;ZZ", k2p) == "P1;P3"
    k2p = gs._key_proteins(df, {"P1": "AT1G01010", "P2": "AT1G01010"})
    assert gs._lead_proteins("AT1G01010", k2p) == "P1;P2"            # deux protéines, une clé


def test_page_is_injected_once_and_js_is_valid(bd, tmp_path):
    html = ("<nav><a onclick=\"showPage('gopage')\">GO</a>\n<button class=\"theme-toggle\"></button></nav>"
            "<div id=\"heatmap\" class=\"page\"></div><!--PROTEOGEN_BODY_END--></body>")
    out = bd.patch_gsea_page(html)
    assert out.count("showPage('gseapage')") == 1 and out.count('id="gseapage"') == 1
    assert out.index("showPage('gopage')") < out.index("showPage('gseapage')")
    assert out.index('id="gseapage"') < out.index('id="heatmap"')
    assert out.index("Page GSEA (07/11/2026)") < out.index("<!--PROTEOGEN_BODY_END-->")
    node = shutil.which("node")
    if node is None:
        pytest.skip("node absent : syntaxe JS non vérifiée")
    js = re.search(r"<script>(.*?)</script>", bd.GSEA_JS, re.S).group(1)
    f = tmp_path / "gsea.js"
    f.write_text(js, encoding="utf-8")
    assert subprocess.run([node, "--check", str(f)], capture_output=True).returncode == 0


def test_manifest_schema_announces_gsea_sheets():
    import pathfinder_bridge as pb
    assert pb.SCHEMA_VERSION >= "1.3"
