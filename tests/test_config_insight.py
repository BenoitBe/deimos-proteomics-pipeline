"""Configuration de l'étape InSight : clés YAML, --insight, question interactive
(posée seulement si un clone est trouvé, réponse par défaut = non)."""
import builtins
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

import config                                          # noqa: E402
import insight_bridge as ib                            # noqa: E402
import pathfinder_bridge as pb                         # noqa: E402


@pytest.fixture(autouse=True)
def _no_pathfinder_clone(monkeypatch):
    monkeypatch.delenv("PATHFINDER_HOME", raising=False)
    monkeypatch.setattr(pb, "pathfinder_search_paths", lambda params: [])


def _resolve(monkeypatch, argv, answers=()):
    monkeypatch.setattr(sys, "argv", ["deimos.py", *argv])
    it = iter(answers)
    asked = []

    def fake_input(prompt=""):
        asked.append(prompt)
        return next(it)
    monkeypatch.setattr(builtins, "input", fake_input)
    params = config.resolve_config(ask_params_fn=lambda: {}, ask_go_params_fn=lambda: None,
                                   go_available=False, dash_available=True,
                                   pr_path="/nonexistent/report.pr_matrix.tsv")
    return params, asked


def _clone(tmp_path):
    d = tmp_path / "InSight"
    (d / "insight").mkdir(parents=True)
    (d / "insight.py").write_text("")
    (d / "insight" / "__init__.py").write_text('__version__ = "0.1.0"\n')
    return d


def test_question_when_a_clone_is_found(tmp_path, monkeypatch):
    monkeypatch.setenv("INSIGHT_HOME", str(_clone(tmp_path)))
    monkeypatch.chdir(tmp_path)
    params, asked = _resolve(monkeypatch, [], answers=["n", "y", "y"])
    assert params["make_insight"] is True
    assert any("InSight" in a for a in asked)


def test_default_answer_is_no(tmp_path, monkeypatch):
    monkeypatch.setenv("INSIGHT_HOME", str(_clone(tmp_path)))
    monkeypatch.chdir(tmp_path)
    params, _ = _resolve(monkeypatch, [], answers=["n", "y", ""])
    assert params["make_insight"] is False


def test_no_question_without_a_clone(tmp_path, monkeypatch):
    monkeypatch.delenv("INSIGHT_HOME", raising=False)
    monkeypatch.setattr(ib, "insight_search_paths", lambda params: [])
    monkeypatch.chdir(tmp_path)
    params, asked = _resolve(monkeypatch, [], answers=["n", "y"])
    assert params["make_insight"] is False
    assert not any("InSight" in a for a in asked)


def test_flag_forces_insight_without_question(tmp_path, monkeypatch):
    monkeypatch.setenv("INSIGHT_HOME", str(_clone(tmp_path)))
    monkeypatch.chdir(tmp_path)
    params, asked = _resolve(monkeypatch, ["--insight"], answers=["n", "y"])
    assert params["make_insight"] is True
    assert not any("InSight" in a for a in asked)


def test_bad_options_are_reported_by_config_validation():
    p = config._resolve(dict(config._DEFAULTS, insight_options={"kk": 3}))
    assert any("insight_options" in e for e in config._validate(p))
    p = config._resolve(dict(config._DEFAULTS))
    assert not any("insight" in e for e in config._validate(p))
