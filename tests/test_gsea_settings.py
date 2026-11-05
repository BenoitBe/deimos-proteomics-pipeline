"""GSEA : choix des jeux de gènes par organisme, fichiers GMT locaux, classement
sur la p brute (sans réseau : gseapy.prerank / Msigdb remplacés par des doublures)."""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import gsea_enrichment as gs  # noqa: E402
import config  # noqa: E402

gseapy = pytest.importorskip("gseapy")


def _df(n=200, seed=0):
    rng = np.random.default_rng(seed)
    genes = [f"Gene{i}" for i in range(n)]          # casse murine
    p = rng.uniform(1e-6, 1, n)
    padj = np.minimum(1, np.sort(p)[np.argsort(np.argsort(p))] * n / (np.argsort(np.argsort(p)) + 1))
    return pd.DataFrame({"Genes": genes, "C_vs_A_diff": rng.normal(0, 1, n),
                         "C_vs_A_p.val": p, "C_vs_A_p.adj": np.clip(padj, 0, 1)})


def test_settings_defaults_by_organism():
    h = gs.gsea_settings({}, {})
    m = gs.gsea_settings({"gsea_organism": "mouse"}, {})
    assert h["gene_sets"] == gs.DEFAULT_GENE_SETS and h["rank_p"] == "pvalue"
    assert m["gene_sets"] == gs.DEFAULT_GENE_SETS_MOUSE and m["msigdb_version"].endswith(".Mm")
    assert gs.gsea_settings({"gsea_gene_sets": "x.gmt"}, {})["gene_sets"] == ["x.gmt"]
    assert gs.gsea_settings({}, {"gsea_gene_sets": ["old"]})["gene_sets"] == ["old"]


def test_ranking_uses_raw_p_by_default():
    df = _df()
    r_p = gs._ranking_metric(df, "C_vs_A", {"volcano_use_padj": True})
    r_v = gs._ranking_metric(df, "C_vs_A", {"volcano_use_padj": True, "gsea_rank_p": "volcano"})
    exp = (-np.log10(df["C_vs_A_p.val"]) * np.sign(df["C_vs_A_diff"]))
    exp.index = df["Genes"].str.upper()
    exp = exp[exp != 0]
    assert dict(zip(map(str, r_p.index), r_p.values)) == pytest.approx(dict(zip(map(str, exp.index), exp.values)))
    assert r_p.nunique() > r_v.nunique()            # moins d'ex aequo qu'avec la p.adj


def test_local_gmt_and_msigdb_are_upper_cased(tmp_path, monkeypatch):
    gmt = tmp_path / "my_sets.gmt"
    gmt.write_text("SET_A\tdesc\tGene1\tGene2\tGene3\nSET_B\tdesc\tgene4\tGene5\n")
    monkeypatch.setattr(gseapy.Msigdb, "get_gmt",
                        classmethod(lambda cls, category, dbver, entrez=False:
                                    {"GOBP_X": ["Actb", "Gapdh"]}
                                    if (category, dbver) == ("m5.go.bp", "2025.1.Mm") else None))
    sets = gs.load_gene_sets([str(gmt), "msigdb:m5.go.bp", "msigdb:nope"], "mouse",
                             "2025.1.Mm", verbose=False)
    assert sets["my_sets__SET_A"] == ["GENE1", "GENE2", "GENE3"]
    assert sets["m5.go.bp__GOBP_X"] == ["ACTB", "GAPDH"]
    assert len(sets) == 3                           # 'nope' ignoré sans planter (A, B, X)


def test_run_gsea_passes_loaded_sets_to_prerank(tmp_path, monkeypatch):
    gmt = tmp_path / "s.gmt"
    gmt.write_text("UP\td\t" + "\t".join(f"Gene{i}" for i in range(20)) + "\n")
    seen = {}

    class _Pre:
        def __init__(self, gene_sets, rnk):
            seen["sets"], seen["rnk"] = gene_sets, rnk
            self.res2d = pd.DataFrame({"Term": ["s__UP"], "NES": [2.0], "FDR q-val": [0.01]})

    monkeypatch.setattr(gseapy, "prerank",
                        lambda rnk, gene_sets, **kw: _Pre(gene_sets, rnk))
    out = gs.run_gsea(_df(), ["C_vs_A"], {"gsea_organism": "mouse", "gsea_gene_sets": [str(gmt)]},
                      {"run_gsea": True}, str(tmp_path))
    assert list(seen["sets"]) == ["s__UP"] and "GENE0" in seen["sets"]["s__UP"]
    assert out and "C_vs_A" in out


def test_config_validation():
    ok = {**config._DEFAULTS}
    assert not [e for e in config._validate(ok) if "gsea" in e]
    bad = {**ok, "gsea_organism": "rat", "gsea_rank_p": "t", "gsea_fdr": 2,
           "gsea_permutations": 10, "gsea_gene_sets": [1]}
    errs = [e for e in config._validate(bad) if "gsea" in e]
    assert len(errs) == 5


def test_interactive_gsea_asks_for_mouse(monkeypatch):
    """Mode interactif : la question GSEA est suivie du choix des jeux de gènes,
    reporté dans params (et donc dans last_config.yaml)."""
    import builtins
    import go_enrichment as ge
    answers = iter(["y", "1", "hsapiens", "1", "y", "2"])   # GO, gProfiler, organisme, fond, GSEA, souris
    asked = []

    def fake_input(prompt=""):
        asked.append(prompt)
        return next(answers, "")
    monkeypatch.setattr(builtins, "input", fake_input)
    go = ge.ask_go_params()
    assert go["run_gsea"] is True and go["gsea_organism"] == "mouse"
    assert any("Gene sets [1/2]" in a for a in asked)


def test_resolve_config_keeps_interactive_gsea_organism(monkeypatch, tmp_path):
    import builtins
    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr(sys, "argv", ["deimos.py", "--out-dir", str(tmp_path / "out")])
    monkeypatch.setattr(builtins, "input", lambda prompt="": "")
    params = config.resolve_config(
        ask_params_fn=lambda: {},
        ask_go_params_fn=lambda: {"organism": "hsapiens", "run_gsea": True,
                                  "gsea_organism": "mouse"},
        go_available=True, dash_available=False, pr_path="/nonexistent/pr.tsv")
    assert params["run_gsea"] is True and params["gsea_organism"] == "mouse"
