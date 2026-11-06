"""Pathfinder et InSight sélectionnent les protéines de Deimos : la règle volcano est
écrite dans le _manifest (texte + clés structurées, schéma 1.2), relue par les deux,
et la feuille Methods décrit la règle réellement appliquée par Pathfinder."""
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

import pathfinder_bridge as pb                         # noqa: E402

DESIGN = None


def _m(params):
    out = {}
    for s, k, v in pb.manifest_rows(params, ["KO_vs_WT"], {}, DESIGN, deimos_commit="x"):
        out.setdefault(s, {})[k] = v
    return out


@pytest.mark.parametrize("padj,thr,ratio,col", [(True, 0.05, 2.0, "p.adj"),
                                                (False, 0.01, 1.5, "p.val")])
def test_manifest_structured_rule(padj, thr, ratio, col):
    st = _m({"volcano_use_padj": padj, "volcano_p_thresh": thr, "volcano_ratio_min": ratio})["stats"]
    assert st["volcano_p_column"] == col
    assert float(st["volcano_p_thresh"]) == thr and float(st["volcano_ratio_min"]) == ratio
    # le texte historique est conservé (lecteurs schéma 1.0/1.1)
    assert st["volcano_rule"].startswith("p.adj" if padj else "p.value")
    assert pb.SCHEMA_VERSION.startswith("1.") and pb.SCHEMA_VERSION >= "1.2"


def _clone(tmp_path, version):
    d = tmp_path / f"pf{version}"
    d.mkdir(exist_ok=True)
    (d / "pathfinder.py").write_text(f'__version__ = "{version}"\n')
    return str(d)


def test_methods_text_follows_the_installed_pathfinder(tmp_path):
    new = pb.pathfinder_selection_text({}, _clone(tmp_path, "1.2.0"))
    assert "volcano rule above" in new and "same proteins" in new and "q < 0.05" in new
    old = pb.pathfinder_selection_text({}, _clone(tmp_path, "1.1.1"))
    assert "p.adj < 0.05 without fold-change criterion" in old
    over = pb.pathfinder_selection_text({"pathfinder_options": {"min_ratio": 1}},
                                        _clone(tmp_path, "1.2.0"))
    assert "overridden" in over and "min_ratio=1" in over


def test_min_ratio_option(tmp_path):
    p = {"pathfinder_options": {"min_ratio": 1.5}}
    argv = pb.build_pathfinder_argv("/pf", str(tmp_path / "r.xlsx"), str(tmp_path / "d.csv"),
                                    str(tmp_path), p)
    assert argv[argv.index("--min-ratio") + 1] == "1.5"
    assert "--pthr" not in pb.build_pathfinder_argv("/pf", "r.xlsx", "d.csv", ".", {})
    assert pb.validate_pathfinder_params(p) == []
    assert pb.validate_pathfinder_params({"pathfinder_options": {"min_ratio": 0.5}})
