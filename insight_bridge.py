# ==============================================================================
# insight_bridge.py — Deimos -> InSight
# Deimos — DIA Expression Integrated Multi-Omics Suite
# ==============================================================================
# Lance InSight (lecture biologique : figures + page récit) dans un PROCESSUS
# SÉPARÉ sur le dossier de sortie qui vient d'être écrit (classeur + matrices/).
# Aucune question, aucun recalcul, et un échec n'interrompt jamais Deimos.
# Même contrat que pathfinder_bridge.py ; bibliothèque standard uniquement.
# ==============================================================================

import os
import subprocess
import sys
import time

_HERE = os.path.dirname(os.path.abspath(__file__))

# Options transmissibles via le YAML (insight_options)
_INSIGHT_OPTIONS = {"k", "max_volcanos", "volcanos", "categorical", "order_regex",
                    "reference_fasta", "title"}


def insight_search_paths(params):
    """Emplacements essayés, dans l'ordre : clé YAML, $INSIGHT_HOME, ../InSight."""
    return [p for p in (params.get("insight_path"),
                        os.environ.get("INSIGHT_HOME"),
                        os.path.join(_HERE, "..", "InSight")) if p]


def locate_insight(params):
    """Dossier d'un clone InSight (insight.py + paquet insight/), ou None."""
    for cand in insight_search_paths(params):
        cand = os.path.expanduser(str(cand))
        if (os.path.isfile(os.path.join(cand, "insight.py"))
                and os.path.isdir(os.path.join(cand, "insight"))):
            return os.path.abspath(cand)
    return None


def insight_version(ins_dir):
    try:
        with open(os.path.join(ins_dir, "insight", "__init__.py"), encoding="utf-8") as fh:
            for line in fh:
                if line.startswith("__version__"):
                    return line.split("=", 1)[1].strip().strip("\"'")
    except OSError:
        pass
    return "unknown"


def build_insight_argv(ins_dir, out_dir, xlsx_path, params):
    opts = params.get("insight_options") or {}
    argv = [sys.executable, os.path.join(ins_dir, "insight.py"), os.path.abspath(out_dir),
            "--out", os.path.join(os.path.abspath(out_dir), "insight")]
    if xlsx_path:
        argv += ["--xlsx", os.path.abspath(xlsx_path)]
    if opts.get("reference_fasta"):
        argv += ["--reference-fasta", os.path.abspath(os.path.expanduser(
            str(opts["reference_fasta"])))]
    if opts.get("k"):
        argv += ["--k", str(int(opts["k"]))]
    if opts.get("max_volcanos"):
        argv += ["--max-volcanos", str(int(opts["max_volcanos"]))]
    for c in opts.get("volcanos") or []:
        argv += ["--volcano", str(c)]
    if opts.get("categorical"):
        argv.append("--categorical")
    if opts.get("order_regex"):
        argv += ["--order-regex", str(opts["order_regex"])]
    if opts.get("title"):
        argv += ["--title", str(opts["title"])]
    return argv


def _tail(path, n=5):
    try:
        with open(path, encoding="utf-8", errors="replace") as fh:
            return [l.rstrip() for l in fh.readlines()[-n:]]
    except OSError:
        return []


def run_insight_step(xlsx_path, params, out_dir):
    """Lance InSight sur le dossier de sortie. Ne lève jamais : retourne le
    chemin de insight_report.html, ou None (raison affichée). Sorties dans
    <out_dir>/insight/, sortie complète dans insight_log.txt."""
    try:
        ins_dir = locate_insight(params)
        if ins_dir is None:
            print("  [SKIP] InSight not found. Looked in: "
                  + ", ".join(os.path.abspath(os.path.expanduser(str(p)))
                              for p in insight_search_paths(params))
                  + " (set insight_path or INSIGHT_HOME).")
            return None
        ins_out = os.path.join(out_dir, "insight")
        os.makedirs(ins_out, exist_ok=True)
        argv = build_insight_argv(ins_dir, out_dir, xlsx_path, params)
        log_path = os.path.join(ins_out, "insight_log.txt")
        env = dict(os.environ, PYTHONIOENCODING="utf-8", PYTHONUTF8="1", MPLBACKEND="Agg")
        timeout = int(params.get("insight_timeout") or 900)
        t0 = time.time()
        with open(log_path, "w", encoding="utf-8") as log:
            log.write("$ " + " ".join(f'"{a}"' if " " in a else a for a in argv) + "\n\n")
            log.flush()
            try:
                # stdin = tube vide (cf. pathfinder_bridge : isatty() False partout)
                r = subprocess.run(argv, cwd=ins_dir, input=b"", stdout=log,
                                   stderr=subprocess.STDOUT, env=env, timeout=timeout)
            except subprocess.TimeoutExpired:
                print(f"  [WARN] InSight stopped (timeout of {timeout} s exceeded) "
                      f"— see {log_path}")
                return None
        html = os.path.join(ins_out, "insight_report.html")
        if r.returncode != 0 or not os.path.exists(html):
            print(f"  [WARN] InSight failed (exit code {r.returncode}) — see {log_path}")
            for line in _tail(log_path):
                print(f"         | {line}")
            return None
        skipped = [l.strip() for l in _tail(log_path, 200) if "[WARN]" in l]
        for line in skipped[:5]:
            print(f"         | {line}")
        print(f"  [OK] InSight {insight_version(ins_dir)} finished in {time.time() - t0:.0f} s")
        return html
    except Exception as e:      # garde-fou : InSight ne doit jamais interrompre Deimos
        print(f"  [WARN] InSight not run ({type(e).__name__}: {e})")
        return None


def validate_insight_params(params):
    """Erreurs de configuration des clés insight_* (liste vide = OK)."""
    errors = []
    opts = params.get("insight_options")
    if opts is not None and not isinstance(opts, dict):
        errors.append("  insight_options doit être un dictionnaire, "
                      "ex: {k: 6, max_volcanos: 4, categorical: false}")
    elif opts:
        unknown = sorted(set(opts) - _INSIGHT_OPTIONS)
        if unknown:
            errors.append(f"  insight_options : clé(s) inconnue(s) {unknown} "
                          f"(attendu : {sorted(_INSIGHT_OPTIONS)})")
        for key, lo in (("k", 2), ("max_volcanos", 0)):
            if opts.get(key) is not None:
                try:
                    if int(opts[key]) < lo:
                        raise ValueError
                except (TypeError, ValueError):
                    errors.append(f"  insight_options.{key} doit être un entier >= {lo}")
        if opts.get("volcanos") is not None and not isinstance(opts["volcanos"], list):
            errors.append("  insight_options.volcanos doit être une liste de contrastes")
    try:
        if int(params.get("insight_timeout") or 900) <= 0:
            raise ValueError
    except (TypeError, ValueError):
        errors.append("  insight_timeout doit être un entier > 0 (secondes)")
    return errors
