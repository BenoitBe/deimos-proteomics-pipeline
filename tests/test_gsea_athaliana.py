"""GSEA sur A. thaliana (et toute espèce g:Profiler) : jeux de gènes du GMT g:Profiler
(téléchargé une fois, en cache), protéines reliées par g:Convert ou par les
orthologues Perseverance (B. napus -> athaliana). Sans réseau : téléchargement,
g:Convert et prerank remplacés par des doublures ; un prerank réel en fin de fichier."""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import gsea_enrichment as gs  # noqa: E402
import config  # noqa: E402

gseapy = pytest.importorskip("gseapy")

N = 120
AGI = [f"AT{1 + i % 5}G{10000 + 10 * i:05d}" for i in range(N)]
ACC = [f"Q{i:05d}" for i in range(N)]                      # UniProt Arabidopsis
NAME = [f"GENE{i}" if i % 3 else a for i, a in enumerate(AGI)]   # 1/3 sans symbole


def _gmt(path):
    up = "\t".join(NAME[:30])                                # 30 premières : en hausse
    down = "\t".join(NAME[-30:])
    path.write_text(
        f"GO:0009414\tresponse to water deprivation\t{up}\n"
        f"GO:0005739\tmitochondrion\t{down}\n"
        f"KEGG:00010\tGlycolysis\t{NAME[40]}\t{NAME[41]}\n"
        f"WP:WP1\tsmall\tNOTINDATA1\tNOTINDATA2\n")
    return str(path)


def _df(prefix=""):
    rng = np.random.default_rng(0)
    lfc = np.r_[rng.uniform(1, 3, 30), rng.normal(0, 0.3, N - 60), rng.uniform(-3, -1, 30)]
    p = np.r_[rng.uniform(1e-6, 1e-3, 30), rng.uniform(0.05, 1, N - 60), rng.uniform(1e-6, 1e-3, 30)]
    return pd.DataFrame({"name": [prefix + a for a in ACC],
                         "Genes": [f"BnaSym{i}" for i in range(N)],
                         "C_vs_A_diff": lfc, "C_vs_A_p.val": p, "C_vs_A_p.adj": np.minimum(1, p * 10)})


def _convert(ids, code):
    assert code == "athaliana"
    m = dict(zip(ACC, zip(AGI, NAME)))
    return [(i, *(m[i] if i in m else ("None", "None"))) for i in ids]


def test_organism_codes_and_settings():
    assert gs.gprofiler_code("athaliana") == "athaliana" == gs.gprofiler_code("Arabidopsis")
    assert gs.gprofiler_code("gprofiler:hvulgare") == "hvulgare"
    assert gs.gprofiler_code("human") is None and gs.gprofiler_code("gprofiler:Bad!") is None
    for o in ("human", "mouse", "athaliana", "gprofiler:osativa"):
        assert gs.valid_gsea_organism(o)
    assert not gs.valid_gsea_organism("rat")
    s = gs.gsea_settings({"gsea_organism": "athaliana"}, {})
    assert s["gprofiler"] == "athaliana" and s["gene_sets"] == ["gprofiler:athaliana"]
    assert gs.gsea_settings({"gsea_organism": "gprofiler:hvulgare"}, {})["gene_sets"] == ["gprofiler:hvulgare"]
    assert gs.gsea_settings({}, {})["gprofiler"] is None          # human inchangé


def test_read_gmt_with_and_without_go_namespaces(tmp_path):
    p = _gmt(tmp_path / "g.gmt")
    lib = gs.read_gprofiler_gmt(p)
    assert ("GO", "GO:0009414 response to water deprivation") in lib
    assert ("KEGG", "KEGG:00010 Glycolysis") in lib
    lib = gs.read_gprofiler_gmt(p, {"GO:0009414": "biological_process",
                                    "GO:0005739": "cellular_component"})
    assert {k[0] for k in lib} == {"GO:BP", "GO:CC", "KEGG", "WP"}


def test_fetch_downloads_once_then_uses_cache(tmp_path, monkeypatch):
    import requests
    calls = []

    class R:
        def __init__(self, code, text):
            self.status_code, self.text = code, text

    def get(url, **kw):
        calls.append(url)
        if "athaliana" in url:
            return R(200, "GO:1\tterm\tA\tB\n")
        if "html" in url:
            return R(200, "<html>error</html>")
        return R(404, "")
    monkeypatch.setattr(requests, "get", get)
    p1 = gs.fetch_gprofiler_gmt("athaliana", str(tmp_path), verbose=False)
    p2 = gs.fetch_gprofiler_gmt("athaliana", str(tmp_path), verbose=False)
    assert p1 == p2 and len(calls) == 1 and calls[0].endswith("gprofiler_full_athaliana.name.gmt")
    with pytest.raises(ValueError, match="unknown to g:Profiler"):
        gs.fetch_gprofiler_gmt("nosuchorg", str(tmp_path), verbose=False)
    monkeypatch.setattr(gs, "GPROFILER_GMT_URL", "https://x/{org}html")
    with pytest.raises(ValueError, match="not a GMT"):
        gs.fetch_gprofiler_gmt("zzz", str(tmp_path), verbose=False)


def test_gene_keys_direct_orthologs_and_fallback():
    km, alias, summ = gs.gprofiler_gene_keys(_df(), "athaliana", convert_fn=_convert)
    assert len(km) == N and km["Q00000"] == AGI[0] and "g:Convert" in summ
    assert alias["GENE1"] == AGI[1] and alias[AGI[0]] == AGI[0] and alias["Q00002"] == AGI[2]
    # B. napus : accessions propres, orthologues Perseverance vers les accessions Arabidopsis
    napus = _df("BNA_")
    orth = {"BNA_" + a: a for a in ACC}
    km, _, summ = gs.gprofiler_gene_keys(napus, "athaliana", orthologs=orth, convert_fn=_convert)
    assert len(km) == N and km["BNA_Q00005"] == AGI[5] and "orthologues" in summ
    # sans orthologues, les accessions B. napus ne se convertissent pas -> repli Genes + conseil
    km, _, summ = gs.gprofiler_gene_keys(napus, "athaliana", convert_fn=_convert)
    assert summ.startswith("[WARN]") and "Perseverance" in summ and km["BNA_Q00000"] == "BNASYM0"

    def offline(ids, code):
        raise ConnectionError("no network")
    _, _, summ = gs.gprofiler_gene_keys(_df(), "athaliana", convert_fn=offline)
    assert "ConnectionError" in summ


def test_translate_drops_members_absent_from_data():
    sets = {"GO__x": ["GENE1", "AT1G10000", "NOTINDATA"], "WP__y": ["NOTINDATA"]}
    out = gs.translate_gene_sets(sets, {"GENE1": "AT2G10010", "AT1G10000": "AT1G10000"})
    assert out == {"GO__x": ["AT2G10010", "AT1G10000"]}


def _patch_network(monkeypatch, tmp_path, seen):
    gmt = _gmt(tmp_path / "gprofiler_full_athaliana.name.gmt")
    monkeypatch.setattr(gs, "fetch_gprofiler_gmt", lambda code, cache_dir=None, verbose=True: gmt)
    monkeypatch.setattr(gs, "_convert_ids", _convert)

    class _Pre:
        def __init__(self, gene_sets, rnk):
            seen["sets"], seen["rnk"] = gene_sets, rnk
            self.res2d = pd.DataFrame({"Term": list(gene_sets)[:1], "NES": [2.0], "FDR q-val": [0.01]})
    monkeypatch.setattr(gseapy, "prerank", lambda rnk, gene_sets, **kw: _Pre(gene_sets, rnk))


def test_run_gsea_athaliana_direct(tmp_path, monkeypatch):
    seen = {}
    _patch_network(monkeypatch, tmp_path, seen)
    out = gs.run_gsea(_df(), ["C_vs_A"], {"gsea_organism": "athaliana"}, {"run_gsea": True},
                      str(tmp_path))
    assert out and "C_vs_A" in out
    assert set(seen["rnk"].index) <= set(AGI) and len(seen["rnk"]) == N
    s = seen["sets"]
    assert set(s["GO__GO:0009414 response to water deprivation"]) == set(AGI[:30])
    assert not any(k.startswith("WP__") for k in s)               # membres absents des données


def test_run_gsea_uses_perseverance_orthologs_only_for_same_species(tmp_path, monkeypatch):
    seen = {}
    _patch_network(monkeypatch, tmp_path, seen)
    napus = _df("BNA_")
    orth = {"organism": "athaliana", "map": {"BNA_" + a: a for a in ACC}}
    gs.run_gsea(napus, ["C_vs_A"], {"gsea_organism": "athaliana"}, {"run_gsea": True},
                str(tmp_path), orthologs=orth)
    assert set(seen["rnk"].index) <= set(AGI)
    seen.clear()
    gs.run_gsea(napus, ["C_vs_A"], {"gsea_organism": "athaliana"}, {"run_gsea": True},
                str(tmp_path), orthologs={"organism": "osativa", "map": orth["map"]})
    assert "rnk" not in seen or not set(seen["rnk"].index) & set(AGI)   # orthologues d'une autre espèce ignorés


def test_config_and_interactive_choice(monkeypatch):
    ok = {**config._DEFAULTS}
    for o in ("athaliana", "gprofiler:hvulgare"):
        assert not [e for e in config._validate({**ok, "gsea_organism": o}) if "gsea" in e]
    assert [e for e in config._validate({**ok, "gsea_organism": "gprofiler:Bad!"}) if "gsea" in e]
    import builtins
    import go_enrichment as ge
    for answers, exp in ((["y", "1", "athaliana", "1", "y", "3"], "athaliana"),
                         (["y", "1", "hvulgare", "1", "y", "4", "hvulgare"], "gprofiler:hvulgare")):
        it = iter(answers)
        monkeypatch.setattr(builtins, "input", lambda prompt="", it=it: next(it, ""))
        assert ge.ask_go_params()["gsea_organism"] == exp


def test_real_prerank_on_translated_sets(tmp_path, monkeypatch):
    """gseapy.prerank réel (100 permutations) sur les jeux traduits : la réponse au
    manque d'eau (30 gènes en tête) sort enrichie vers le haut, la mitochondrie vers le bas."""
    gmt = _gmt(tmp_path / "gprofiler_full_athaliana.name.gmt")
    monkeypatch.setattr(gs, "fetch_gprofiler_gmt", lambda code, cache_dir=None, verbose=True: gmt)
    monkeypatch.setattr(gs, "_convert_ids", _convert)
    out = gs.run_gsea(_df(), ["C_vs_A"], {"gsea_organism": "athaliana", "gsea_permutations": 100},
                      {"run_gsea": True}, str(tmp_path))
    t = out["C_vs_A"]["table"].set_index("Term")
    assert t.loc["GO__GO:0009414 response to water deprivation", "NES"] > 0
    assert t.loc["GO__GO:0005739 mitochondrion", "NES"] < 0
    prov = out["C_vs_A"]["provenance"]
    assert "g:Profiler athaliana GMT" in prov and "g:Convert" in prov
    from openpyxl import Workbook
    wb = Workbook()
    gs.export_gsea_sheets(wb, out, lambda ws, img, cell: None)
    assert wb["GSEA_CvA"]["A1"].value == prov
