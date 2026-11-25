"""Backend STRING : couverture directe faible -> orthologie vers le protéome STRING
de la même espèce (si DIAMOND est disponible)."""
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import go_enrichment as ge  # noqa: E402


def _det(share):
    return {"taxid": 10181, "name": "Heterocephalus glaber", "method": "uniprot",
            "share": share, "others": {}, "reliable": True, "unreachable": False,
            "strain_taxid": 10181, "recommendation": "x"}


@pytest.mark.parametrize("share, diamond, mode", [
    (0.32, True, "ortholog"),        # cas du rat-taupe nu (16/50)
    (0.32, False, "direct"),         # sans DIAMOND : direct + avertissement
    (0.80, True, "direct"),          # couverture suffisante
])
def test_auto_detection(monkeypatch, share, diamond, mode, capsys):
    monkeypatch.setattr(ge, "detect_string_taxid", lambda ids, caller: _det(share))
    monkeypatch.setattr(ge, "_diamond_available", lambda: diamond)
    plan = ge.plan_string_mapping(["A0A1"], "auto", "x")
    assert plan["mode"] == mode and plan["taxid"] == 10181
    if share < ge.STRING_DIRECT_MIN_COVERAGE and not diamond:
        assert "Install DIAMOND" in capsys.readouterr().out


def test_explicit_taxid_low_coverage(monkeypatch):
    monkeypatch.setattr(ge, "diagnose_string_coverage",
                        lambda ids, T, caller, **kw: {"reliable": True, "coverage": 0.4,
                                                       "recommendation": "x"})
    monkeypatch.setattr(ge, "_diamond_available", lambda: True)
    assert ge.plan_string_mapping(["A0A1"], 10181, "x")["mode"] == "ortholog"
    monkeypatch.setattr(ge, "diagnose_string_coverage",
                        lambda ids, T, caller, **kw: {"reliable": True, "coverage": 0.9,
                                                       "recommendation": "x"})
    assert ge.plan_string_mapping(["A0A1"], 10181, "x")["mode"] == "direct"


def test_perseverance_not_run_with_string_backend():
    """Le pipeline n'appelle pas Perseverance quand le backend GO est STRING."""
    src = (Path(__file__).resolve().parents[1] / "deimos.py").read_text(encoding="utf-8")
    i = src.index('elif local_go_params and local_go_params.get("backend") == "string":')
    j = src.index("run_perseverance_step(df_results, params, out_dir)", i)
    assert "correspondence_df, persev_override = None, None" in src[i:j]


def test_uniparc_ids_without_gene_names_go_to_orthology(monkeypatch):
    """Recherche DIA-NN sur un FASTA UniParc (en-têtes '>UPI… status=active') : ni
    accession indexée dans STRING ni nom de gène (colonne Genes vide, lue NaN puis
    'nan' par astype(str)). 'nan' ressemble à un symbole bactérien (3 minuscules) :
    il ne doit pas servir de sonde, sinon 0/1 gène reconnu désactive le GO STRING
    au lieu de passer en orthologie vers le taxid demandé."""
    calls = []
    monkeypatch.setattr(ge, "diagnose_string_coverage",
                        lambda ids, T, caller, **kw: {"reliable": False, "unreachable": False,
                                                       "coverage": 0.0, "recommendation": "x"})
    monkeypatch.setattr(ge, "detect_string_taxid",
                        lambda ids, caller: {"reliable": False, "unreachable": False, "taxid": None,
                                             "name": None, "strain_taxid": None,
                                             "recommendation": "x"})

    def fake_map(names, taxid, caller, **kw):
        calls.append(list(names))
        import pandas as pd
        return pd.DataFrame(columns=["queryItem"])
    monkeypatch.setattr(ge, "_string_id_map", fake_map)
    monkeypatch.setattr(ge, "_diamond_available", lambda: True)
    ids = [f"UPI0000{i:06d}" for i in range(50)]
    plan = ge.plan_string_mapping(ids, 93061, "x", gene_names=["nan"] * 50, reference_taxid=93061)
    assert plan["mode"] == "ortholog" and plan["taxid"] == 93061
    assert not any("nan" in c for c in calls)        # 'nan' never sent as a gene name


def test_probe_ignores_missing_gene_tokens(monkeypatch):
    seen = []
    monkeypatch.setattr(ge, "_string_id_map",
                        lambda names, taxid, caller, **kw: seen.append(list(names)))
    assert ge.probe_string_organism(["nan", "None", "", "NaN", "null"], 93061) == (0, 0)
    assert seen == []
