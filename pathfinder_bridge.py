# ==============================================================================
# pathfinder_bridge.py — Deimos -> Pathfinder
# Deimos — DIA Expression Integrated Multi-Omics Suite
# ==============================================================================
# Deux responsabilités, rien d'autre :
#   1. write_manifest()      : onglet caché '_manifest' = contrat explicite et
#                              versionné entre le classeur Deimos et ses lecteurs
#                              (Pathfinder en premier). Une ligne par fait :
#                              section | key | value.
#   2. run_pathfinder_step() : lance Pathfinder dans un PROCESSUS SÉPARÉ sur le
#                              classeur qui vient d'être écrit. Aucune question,
#                              aucun recalcul, et un échec n'interrompt jamais Deimos.
#
# Module léger (bibliothèque standard uniquement) : importé aussi par config.py
# pour la question interactive, sans charger deimos.py.
# ==============================================================================

import os
import re
import subprocess
import sys
import time

MANIFEST_SHEET = "_manifest"
SCHEMA_NAME = "deimos-results"
SCHEMA_VERSION = "1.3"   # majeure +1 si un onglet/une colonne lu en aval est renommé ;
                         # mineure +1 pour un ajout (1.1 : section 'sensitivity' ;
                         # 1.2 : règle volcano structurée, volcano_p_column/_thresh/_ratio_min;
                         # 1.3 : section gsea_sheet (onglet GSEA_* -> contraste), ajout)

STATS_SHEET = "Differential_Expression"
MATRIX_SHEET = "Log2_Impute"

# Index de réplicat (1, 2, 3 dans chaque condition) : pas un facteur du design.
_NOT_FACTORS = {"label", "condition", "replicate", "rep", "replicat"}

# Options Pathfinder transmissibles via le YAML (pathfinder_options)
_PF_OPTIONS = {"pthr", "min_ratio", "jaccard", "no_umap",
               # lecture biologique (ex-InSight, Pathfinder >= 1.3.0)
               "no_reading", "k", "volcanos", "max_volcanos", "categorical", "order_regex",
               "reference_fasta", "groups_from", "prog_boot", "factor_regex"}
# Options d'InSight reprises par la lecture biologique de Pathfinder >= 1.3.0, quand
# pathfinder_options ne les donne pas (make_insight absorbé par make_pathfinder)
_READING_FROM_INSIGHT = ("k", "volcanos", "max_volcanos", "categorical", "order_regex",
                         "reference_fasta")

_HERE = os.path.dirname(os.path.abspath(__file__))


# ── Versions ──────────────────────────────────────────────────────────────────

def _git_describe(path):
    """Commit court du dépôt (suffixe '-dirty' si des fichiers suivis sont modifiés),
    ou None hors dépôt git / git absent."""
    try:
        r = subprocess.run(["git", "-C", path, "describe", "--always", "--dirty",
                            "--abbrev=7"], capture_output=True, text=True, timeout=10,
                           stdin=subprocess.DEVNULL)
    except (OSError, subprocess.SubprocessError):
        return None
    out = r.stdout.strip()
    return out if r.returncode == 0 and out else None


def deimos_version():
    return _git_describe(_HERE) or "unknown"


def pathfinder_version(pf_dir):
    """'1.1.0 (abc1234)' : __version__ lu dans pathfinder.py (sans l'importer) + commit."""
    ver = None
    try:
        with open(os.path.join(pf_dir, "pathfinder.py"), encoding="utf-8") as f:
            m = re.search(r'^__version__\s*=\s*["\']([^"\']+)["\']', f.read(), re.M)
            ver = m.group(1) if m else None
    except OSError:
        pass
    commit = _git_describe(pf_dir)
    if ver and commit:
        return f"{ver} ({commit})"
    return ver or commit or "unknown version"


# ── Organisme réellement analysé ──────────────────────────────────────────────

def pathfinder_selection_text(params, pf_dir):
    """Règle de sélection des protéines appliquée par Pathfinder, pour la feuille Methods.
    Pathfinder >= 1.2 applique la règle volcano du manifeste (mêmes protéines que les
    volcanos) sauf surcharge pthr / min_ratio ; les versions antérieures, p.adj < pthr
    sans critère de fold-change."""
    opts = params.get("pathfinder_options") or {}
    ver = pathfinder_version(pf_dir) if pf_dir else ""
    m = re.match(r"(\d+)\.(\d+)", ver or "")
    modern = bool(m) and (int(m.group(1)), int(m.group(2))) >= (1, 2)
    if not modern:
        return (f"Protein signature = p.adj < {opts.get('pthr', 0.05)} without fold-change "
                "criterion (independent of the volcano thresholds above)")
    if "pthr" in opts or "min_ratio" in opts:
        return ("Protein selection = volcano rule above, overridden by pathfinder_options "
                f"({', '.join(f'{k}={opts[k]}' for k in ('pthr', 'min_ratio') if k in opts)}): "
                "counts differ from the volcanos. GO / pathway terms: q < 0.05 as computed "
                "above")
    return ("Protein selection = the volcano rule above (same proteins as the volcanos, "
            "UpSet and dashboard); GO / pathway terms: q < 0.05 as computed above")


def analysed_organism(params, perseverance=False):
    """('taxid', '1304') en backend STRING, ('gprofiler', 'bnapus') sinon,
    (None, None) si inconnu.

    En mode Perseverance le GO porte sur l'espèce de RÉFÉRENCE, pas sur
    l'organisme mesuré : aucune espèce n'est déclarée (pas d'axe Publications).
    En backend STRING, go_species_taxid a déjà été remplacé par main() par le
    taxid réellement utilisé (organisme indexé, ou espèce parente en orthologie)."""
    if perseverance:
        return None, None
    tx = params.get("go_species_taxid")
    if params.get("go_backend", "gprofiler") == "string":
        if tx not in (None, "", "auto", "__auto__"):
            return "taxid", str(tx)
        return None, None
    if params.get("go_organism"):
        return "gprofiler", str(params["go_organism"])
    return None, None


def pathfinder_org(params, perseverance=False):
    """Valeur passée à 'pathfinder --org' (taxid STRING ou code gProfiler), ou None."""
    return analysed_organism(params, perseverance)[1]


# ── 1. Manifeste ──────────────────────────────────────────────────────────────

def manifest_rows(params, contrast_names, go_sheet_map, design, *,
                  perseverance=False, primary_statistic="limma eBayes",
                  deimos_commit=None, extra_rows=None):
    """Lignes (section, key, value) du manifeste — séparé de l'écriture pour être testé."""
    p_col = "p.adj" if params.get("volcano_use_padj") else "p.value"
    rows = [("schema", "name", SCHEMA_NAME),
            ("schema", "version", SCHEMA_VERSION),
            ("run", "deimos_commit", deimos_commit or deimos_version()),
            ("run", "date", time.strftime("%Y-%m-%d %H:%M:%S")),
            ("stats", "sheet", STATS_SHEET),
            ("stats", "primary_statistic", primary_statistic),
            ("stats", "p_adjust", "BH global" if params.get("fdr_global") else "BH per contrast"),
            ("stats", "volcano_rule",
             f"{p_col} < {params.get('volcano_p_thresh', 0.05)} & ratio >= "
             f"{params.get('volcano_ratio_min', 1.5)}"),
            # même règle, structurée (schéma 1.2) : les lecteurs (Pathfinder, InSight)
            # sélectionnent exactement les protéines des volcanos de Deimos
            ("stats", "volcano_p_column", "p.adj" if params.get("volcano_use_padj") else "p.val"),
            ("stats", "volcano_p_thresh", str(params.get("volcano_p_thresh", 0.05))),
            ("stats", "volcano_ratio_min", str(params.get("volcano_ratio_min", 1.5))),
            ("matrix", "sheet", MATRIX_SHEET),
            ("matrix", "imputed", "true"),
            ("matrix", "imputation",
             f"{str(params.get('impute_method', 'auto')).upper()}, seed "
             f"{params.get('seed', 42)}")]
    for c in contrast_names:
        num, _, den = str(c).partition("_vs_")
        rows.append(("contrast", str(c), f"{num}|{den}"))
    for c, sheet in (go_sheet_map or {}).items():
        rows.append(("go_sheet", str(sheet), str(c)))
    if go_sheet_map:
        backend = "perseverance" if perseverance else params.get("go_backend", "gprofiler")
        rows.append(("go", "backend", backend))
        if perseverance and params.get("perseverance_reference_organism"):
            rows.append(("go", "reference_organism",
                         str(params["perseverance_reference_organism"])))
    kind, value = analysed_organism(params, perseverance)
    if kind == "taxid":
        rows.append(("organism", "taxid", value))
        rows.append(("organism", "string_mode", params.get("go_string_mode") or ""))
    elif kind == "gprofiler":
        rows.append(("organism", "gprofiler", value))
    cols = [] if design is None else list(design.columns)
    factors = [str(c) for c in cols if str(c).lower() not in _NOT_FACTORS]
    rows.append(("design", "condition_column", "condition"))
    rows.append(("design", "factors", ", ".join(factors)))
    rows.extend(tuple(map(str, r)) for r in (extra_rows or []))
    return rows


def write_manifest(wb, params, contrast_names, go_sheet_map, design, **kw):
    """Onglet caché '_manifest', dernier du classeur. À appeler juste avant wb.save().
    Lisible par pandas/openpyxl, invisible pour le client. Retourne la feuille."""
    if MANIFEST_SHEET in wb.sheetnames:
        del wb[MANIFEST_SHEET]
    ws = wb.create_sheet(MANIFEST_SHEET)
    ws.append(("section", "key", "value"))
    for r in manifest_rows(params, contrast_names, go_sheet_map, design, **kw):
        ws.append(r)
    ws.column_dimensions["A"].width = 12
    ws.column_dimensions["B"].width = 34
    ws.column_dimensions["C"].width = 60
    ws.sheet_state = "hidden"
    return ws


# ── 2. Lancement de Pathfinder ────────────────────────────────────────────────

def pathfinder_search_paths(params):
    """Emplacements essayés, dans l'ordre : clé YAML, $PATHFINDER_HOME, ../Pathfinder."""
    return [p for p in (params.get("pathfinder_path"),
                        os.environ.get("PATHFINDER_HOME"),
                        os.path.join(_HERE, "..", "Pathfinder")) if p]


def locate_pathfinder(params):
    """Dossier d'un clone Pathfinder (contenant pathfinder.py), ou None."""
    for cand in pathfinder_search_paths(params):
        cand = os.path.expanduser(str(cand))
        if os.path.isfile(os.path.join(cand, "pathfinder.py")):
            return os.path.abspath(cand)
    return None


def _tail(path, n=5):
    try:
        with open(path, encoding="utf-8", errors="replace") as f:
            lines = [l.rstrip() for l in f if l.strip()]
        return lines[-n:]
    except OSError:
        return []


def pathfinder_has_reading(pf_dir):
    """True si ce Pathfinder porte la lecture biologique (>= 1.3.0, option --no-reading)."""
    try:
        with open(os.path.join(pf_dir, "pathfinder.py"), encoding="utf-8") as f:
            return "--no-reading" in f.read()
    except OSError:
        return False


def resolve_reading(params, verbose=True):
    """make_insight absorbé par make_pathfinder : un seul rapport HTML.

    Si un Pathfinder >= 1.3.0 (lecture biologique intégrée) est trouvé, demander InSight
    revient à demander Pathfinder : make_pathfinder passe à True, make_insight à False
    (pas de second rapport insight_report.html) et les insight_options sont reprises par
    Pathfinder (reading_argv). Sans ce Pathfinder, InSight reste lancé à part, comme avant.
    Retourne 'pathfinder', 'insight' ou None (lecture non demandée)."""
    if not params.get("make_insight"):
        return None
    pf = locate_pathfinder(params)
    if pf and pathfinder_has_reading(pf):
        was_pf = bool(params.get("make_pathfinder"))
        params["make_pathfinder"] = True
        params["make_insight"] = False
        params["_reading_in_pathfinder"] = True
        if verbose:
            print("  [INFO] InSight is merged into Pathfinder >= 1.3.0: the biological reading is "
                  "the 'Biological reading' tab of pathfinder_report.html (no separate "
                  "insight_report.html)"
                  + ("." if was_pf else "; Pathfinder is run for it."))
        return "pathfinder"
    return "insight"


def reading_argv(params):
    """Options de la lecture biologique : pathfinder_options, puis insight_options."""
    opts = dict(params.get("pathfinder_options") or {})
    for k, v in (params.get("insight_options") or {}).items():
        if k in _READING_FROM_INSIGHT and k not in opts:
            opts[k] = v
    argv = []
    if opts.get("no_reading"):
        return ["--no-reading"]
    if opts.get("k"):
        argv += ["--k", str(int(opts["k"]))]
    for c in opts.get("volcanos") or []:
        argv += ["--volcano", str(c)]
    if opts.get("max_volcanos") is not None:
        argv += ["--max-volcanos", str(int(opts["max_volcanos"]))]
    if opts.get("categorical"):
        argv.append("--categorical")
    if opts.get("order_regex"):
        argv += ["--order-regex", str(opts["order_regex"])]
    if opts.get("reference_fasta"):
        argv += ["--reference-fasta", os.path.abspath(str(opts["reference_fasta"]))]
    if opts.get("groups_from"):
        argv += ["--groups-from", str(opts["groups_from"])]
    if opts.get("prog_boot") is not None:
        argv += ["--prog-boot", str(int(opts["prog_boot"]))]
    if opts.get("factor_regex"):
        argv += ["--factor-regex", str(opts["factor_regex"])]
    return argv


def build_pathfinder_argv(pf_dir, xlsx_path, design_csv, pf_out, params, *,
                          go_ran=True, perseverance=False, reading=None):
    argv = [sys.executable, os.path.join(pf_dir, "pathfinder.py"),
            "--stats", os.path.abspath(xlsx_path),
            # jamais d'auto-détection : un classeur GO d'un autre projet pourrait être pris
            "--go", os.path.abspath(xlsx_path) if go_ran else "none",
            "--design", os.path.abspath(design_csv),
            "--outdir", os.path.abspath(pf_out)]
    org = pathfinder_org(params, perseverance)
    if params.get("pathfinder_publications") and org:
        argv += ["--org", org]
    else:
        argv.append("--no-online")
    opts = params.get("pathfinder_options") or {}
    if "pthr" in opts:
        argv += ["--pthr", str(float(opts["pthr"]))]
    if "min_ratio" in opts:
        argv += ["--min-ratio", str(float(opts["min_ratio"]))]
    if "jaccard" in opts:
        argv += ["--jaccard", str(float(opts["jaccard"]))]
    if opts.get("no_umap"):
        argv.append("--no-umap")
    # lecture biologique : seulement si ce Pathfinder la connaît (un 1.2 refuserait l'option)
    if pathfinder_has_reading(pf_dir) if reading is None else reading:
        argv += reading_argv(params)
    return argv


def run_pathfinder_step(xlsx_path, design, params, out_dir, *, go_ran=True,
                        perseverance=False):
    """Lance Pathfinder sur le classeur qui vient d'être écrit.

    Ne lève jamais : retourne le chemin du rapport HTML, ou None (raison affichée).
    Sorties dans <out_dir>/pathfinder/, sortie complète dans pathfinder_log.txt."""
    try:
        pf_dir = locate_pathfinder(params)
        if pf_dir is None:
            print("  [SKIP] Pathfinder not found. Looked in: "
                  + ", ".join(os.path.abspath(os.path.expanduser(str(p)))
                              for p in pathfinder_search_paths(params))
                  + " (set pathfinder_path or PATHFINDER_HOME).")
            return None
        if params.get("pathfinder_publications") and not pathfinder_org(params, perseverance):
            print("  [INFO] Publications axis requested but no species is known "
                  "(go_organism or STRING taxid; never in Perseverance mode): skipped.")
        pf_out = os.path.join(out_dir, "pathfinder")
        os.makedirs(pf_out, exist_ok=True)

        # Design réellement utilisé (après filtrage), sans l'index de réplicat
        keep = [c for c in design.columns
                if c == "label" or str(c).lower() not in (_NOT_FACTORS - {"condition"})]
        design_csv = os.path.join(pf_out, "design_pathfinder.csv")
        design[keep].to_csv(design_csv, index=False)

        argv = build_pathfinder_argv(pf_dir, xlsx_path, design_csv, pf_out, params,
                                     go_ran=go_ran, perseverance=perseverance)
        log_path = os.path.join(pf_out, "pathfinder_log.txt")
        # Windows : sortie redirigée en cp1252 sinon (⚠, ε², — font planter Pathfinder)
        env = dict(os.environ, PYTHONIOENCODING="utf-8", PYTHONUTF8="1", MPLBACKEND="Agg")
        timeout = int(params.get("pathfinder_timeout") or 1800)
        t0 = time.time()
        with open(log_path, "w", encoding="utf-8") as log:
            log.write("$ " + " ".join(f'"{a}"' if " " in a else a for a in argv) + "\n\n")
            log.flush()
            try:
                # stdin = tube vide et fermé (et non DEVNULL) : sous Windows, NUL est
                # un périphérique caractère, isatty() y répond True et Pathfinder
                # croirait pouvoir poser une question. Avec un tube : isatty()
                # False et lecture = fin de fichier immédiate, sur tous les systèmes.
                r = subprocess.run(argv, cwd=pf_out, input=b"",
                                   stdout=log, stderr=subprocess.STDOUT, env=env,
                                   timeout=timeout)
            except subprocess.TimeoutExpired:
                print(f"  [WARN] Pathfinder stopped (timeout of {timeout} s exceeded) "
                      f"— see {log_path}")
                return None
        html = os.path.join(pf_out, "pathfinder_report.html")
        if r.returncode != 0 or not os.path.exists(html):
            print(f"  [WARN] Pathfinder failed (exit code {r.returncode}) — see {log_path}")
            for line in _tail(log_path):
                print(f"         | {line}")
            return None
        print(f"  [OK] Pathfinder {pathfinder_version(pf_dir)} finished in {time.time() - t0:.0f} s")
        return html
    except Exception as e:      # garde-fou : Pathfinder ne doit jamais interrompre Deimos
        print(f"  [WARN] Pathfinder not run ({type(e).__name__}: {e})")
        return None


def validate_pathfinder_params(params):
    """Erreurs de configuration des clés pathfinder_* (liste vide = OK)."""
    errors = []
    opts = params.get("pathfinder_options")
    if opts is not None and not isinstance(opts, dict):
        errors.append("  pathfinder_options doit être un dictionnaire, "
                      "ex: {pthr: 0.05, jaccard: 0.5, no_umap: true}")
    elif opts:
        unknown = sorted(set(opts) - _PF_OPTIONS)
        if unknown:
            errors.append(f"  pathfinder_options : clé(s) inconnue(s) {unknown} "
                          f"(attendu : {sorted(_PF_OPTIONS)})")
        if "min_ratio" in opts:
            try:
                if float(opts["min_ratio"]) < 1:
                    raise ValueError
            except (TypeError, ValueError):
                errors.append("  pathfinder_options.min_ratio doit être >= 1 (1 = sans critère "
                              "de fold-change)")
        if "groups_from" in opts and opts["groups_from"] not in ("auto", "ora", "gsea"):
            errors.append("  pathfinder_options.groups_from : auto, ora ou gsea")
        for k, lo in (("k", 2), ("max_volcanos", 0), ("prog_boot", 0)):
            if k in opts:
                try:
                    if int(opts[k]) < lo or float(opts[k]) != int(opts[k]):
                        raise ValueError
                except (TypeError, ValueError):
                    errors.append(f"  pathfinder_options.{k} doit être un entier >= {lo}")
        if opts.get("factor_regex"):
            try:
                if not re.compile(str(opts["factor_regex"])).groupindex:
                    errors.append("  pathfinder_options.factor_regex : groupes nommés requis, "
                                  "ex: '(?P<genotype>[^_]+)_T(?P<time>\\d+)'")
            except re.error as e:
                errors.append(f"  pathfinder_options.factor_regex invalide ({e})")
        if "volcanos" in opts and not isinstance(opts["volcanos"], (list, tuple)):
            errors.append("  pathfinder_options.volcanos : liste de contrastes, ex: [A_vs_B]")
        for k in ("pthr", "jaccard"):
            if k in opts:
                try:
                    v = float(opts[k])
                    if not 0 < v <= 1:
                        raise ValueError
                except (TypeError, ValueError):
                    errors.append(f"  pathfinder_options.{k} doit être dans ]0, 1]")
    try:
        if int(params.get("pathfinder_timeout") or 1800) <= 0:
            raise ValueError
    except (TypeError, ValueError):
        errors.append("  pathfinder_timeout doit être un entier > 0 (secondes)")
    return errors
