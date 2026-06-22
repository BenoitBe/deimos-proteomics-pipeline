# ==============================================================================
# config.py — Gestion de la configuration YAML réutilisable
# Deimos — DIA Expression Integrated Multi-Omics Suite
# ==============================================================================
# Usage CLI :
#   python deimos.py                        → questions interactives
#   python deimos.py --config mon_projet.yaml  → chargement direct
#   python deimos.py --save-config          → questions + sauvegarde
#
# Après chaque run, last_config.yaml est toujours écrit dans out_dir.
# ==============================================================================

import os
import sys
import json
import argparse
from datetime import datetime
from typing import Optional

# PyYAML optionnel → fallback JSON transparent
try:
    import yaml
    _YAML_OK = True
except ImportError:
    _YAML_OK = False


# ── Clés obligatoires et leurs valeurs par défaut ──────────────────────────────
_DEFAULTS = {
    # Volcanos / Robustness
    "volcano_use_padj":    False,
    "volcano_p_thresh":    0.05,
    "volcano_ratio_min":   1.5,
    # ANOVA / Heatmap
    "anova_use_padj":      False,
    "anova_p_thresh":      0.05,
    # Heatmap clusters
    "n_heatmap_clusters":  3,
    # Robustness
    "n_iter_robustness":   0,
    # Correction FDR
    "fdr_global":          False,
    # Imputation
    "impute_method":       "qrilc",
    # Modules optionnels
    "go_organism":         None,    # None = désactivé
    "make_wgcna":          False,
    "make_dashboard":      True,
    "use_deqms":           False,
    # Chemins (peuvent être surchargés)
    "tsv_path":            "report.pg_matrix.tsv",
    "design_path":         "ExperimentalDesign.csv",
    "pr_path":             "report.pr_matrix.tsv",
    "out_dir":             "proteomics_output",
}

# Clés dérivées calculées automatiquement (non stockées dans le YAML)
_COMPUTED = {"volcano_lfc_min"}


# ── Helpers de sérialisation ────────────────────────────────────────────────────

def _dumps(obj: dict) -> str:
    """Sérialise en YAML si disponible, sinon JSON indenté."""
    if _YAML_OK:
        return yaml.dump(obj, allow_unicode=True, sort_keys=False,
                         default_flow_style=False)
    return json.dumps(obj, ensure_ascii=False, indent=2)


def _loads(text: str) -> dict:
    """Désérialise YAML ou JSON (auto-détection)."""
    stripped = text.lstrip()
    if stripped.startswith("{") or stripped.startswith("["):
        return json.loads(text)
    if _YAML_OK:
        return yaml.safe_load(text)
    # Pas de yaml installé et pas du JSON → erreur claire
    raise RuntimeError(
        "PyYAML non installé et le fichier ne semble pas être du JSON.\n"
        "  → pip install pyyaml"
    )


# ── Résolution des clés dérivées ────────────────────────────────────────────────

def _resolve(params: dict) -> dict:
    """Calcule les clés dérivées à partir des clés brutes."""
    import numpy as np
    params["volcano_lfc_min"] = float(np.log2(params["volcano_ratio_min"]))
    return params


# ── Validation ──────────────────────────────────────────────────────────────────

def _validate(params: dict) -> list[str]:
    """Retourne la liste des erreurs (vide = OK)."""
    errors = []
    for p_key in ("volcano_p_thresh", "anova_p_thresh"):
        v = params.get(p_key)
        if not (isinstance(v, (int, float)) and 0 < v <= 1):
            errors.append(f"  {p_key} doit être dans ]0, 1] (valeur : {v!r})")
    if params.get("volcano_ratio_min", 0) < 1:
        errors.append(f"  volcano_ratio_min doit être ≥ 1 "
                      f"(valeur : {params.get('volcano_ratio_min')!r})")
    if params.get("n_heatmap_clusters", 0) < 2:
        errors.append(f"  n_heatmap_clusters doit être ≥ 2")
    if params.get("n_iter_robustness", 0) < 0:
        errors.append(f"  n_iter_robustness doit être ≥ 0 (0 = désactivé)")
    if params.get("impute_method") not in ("qrilc", "mixed"):
        errors.append(f"  impute_method doit être 'qrilc' ou 'mixed'")
    return errors


# ── Chargement depuis fichier ───────────────────────────────────────────────────

def load_config(path: str) -> dict:
    """
    Charge un fichier YAML ou JSON de configuration.
    Complète les clés manquantes avec _DEFAULTS.
    Valide et lève ValueError si incohérent.
    """
    if not os.path.exists(path):
        raise FileNotFoundError(f"Fichier de config introuvable : {path}")

    with open(path, "r", encoding="utf-8") as f:
        raw = _loads(f.read())

    if not isinstance(raw, dict):
        raise ValueError(f"Le fichier de config doit être un dictionnaire (trouvé : {type(raw).__name__})")

    # Fusion avec les defaults (les clés du fichier ont priorité)
    params = {**_DEFAULTS, **raw}

    errors = _validate(params)
    if errors:
        raise ValueError("Config invalide :\n" + "\n".join(errors))

    return _resolve(params)


# ── Sauvegarde ──────────────────────────────────────────────────────────────────

def save_config(params: dict, path: str, comment: str = "") -> None:
    """
    Sauvegarde les paramètres dans un fichier YAML (ou JSON si PyYAML absent).
    Exclut les clés calculées (_COMPUTED) et ajoute un en-tête horodaté.
    """
    # On ne stocke que les clés sérialisables (pas les clés dérivées)
    out = {k: v for k, v in params.items()
           if k not in _COMPUTED and k in _DEFAULTS}

    ext = ".yaml" if _YAML_OK else ".json"
    # Forcer l'extension correcte si l'utilisateur n'a pas mis d'extension
    base, file_ext = os.path.splitext(path)
    if file_ext.lower() not in (".yaml", ".yml", ".json"):
        path = base + ext

    header_lines = [
        f"# Pipeline Protéomique — Config générée le {datetime.now().strftime('%Y-%m-%d %H:%M')}",
    ]
    if comment:
        header_lines.append(f"# {comment}")
    header_lines.append("")
    header = "\n".join(header_lines)

    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    with open(path, "w", encoding="utf-8") as f:
        if _YAML_OK:
            f.write(header)
        f.write(_dumps(out))

    print(f"  [SAVE] Config saved: {path}")


# ── Affichage récapitulatif ─────────────────────────────────────────────────────

def print_config_summary(params: dict) -> None:
    """Affiche un récapitulatif lisible des paramètres chargés."""
    import numpy as np
    v_stat  = "p.adj (FDR)" if params["volcano_use_padj"] else "raw p.value"
    a_stat  = "p.adj (FDR)" if params["anova_use_padj"]   else "raw p.value"
    fdr_txt = "global (whole study)" if params["fdr_global"] else "per contrast"
    rob_txt = ("disabled" if params["n_iter_robustness"] == 0
               else f"{params['n_iter_robustness']} iterations")
    imp_txt = ("QRILC (pure MNAR)" if params["impute_method"] == "qrilc"
               else "Mixed (QRILC + kNN)")
    go_txt  = params.get("go_organism") or "no"
    wgcna_txt = "yes" if params.get("make_wgcna") else "no"
    dash_txt  = "yes" if params.get("make_dashboard") else "no"
    deqms_txt = "yes" if params.get("use_deqms") else "no"

    print("\n  [OK] Parameters:")
    print(f"     Volcanos/Robustness -> {v_stat} < {params['volcano_p_thresh']}"
          f"  |  ratio >= {params['volcano_ratio_min']}"
          f"  (log2FC >= {params['volcano_lfc_min']:.3f})")
    print(f"     ANOVA/Heatmap       -> {a_stat} < {params['anova_p_thresh']}")
    print(f"     Heatmap clusters    -> {params['n_heatmap_clusters']}")
    print(f"     Robustness          -> {rob_txt}")
    print(f"     FDR correction      -> {fdr_txt}")
    print(f"     Imputation          -> {imp_txt}")
    print(f"     GO (organism)       -> {go_txt}")
    print(f"     WGCNA               -> {wgcna_txt}")
    print(f"     Dashboard           -> {dash_txt}")
    print(f"     DEqMS               -> {deqms_txt}")


# ── Parsing des arguments CLI ───────────────────────────────────────────────────

def parse_cli_args() -> argparse.Namespace:
    """
    Parse les arguments de ligne de commande.
    Appelé UNE seule fois au début de main().
    """
    parser = argparse.ArgumentParser(
        description="Pipeline protéomique LFQ DIA — Proteogen",
        formatter_class=argparse.RawTextHelpFormatter,
    )
    parser.add_argument(
        "--config", "-c",
        metavar="CONFIG.yaml",
        help="Charger une configuration YAML/JSON existante (court-circuite les questions).",
    )
    parser.add_argument(
        "--save-config",
        metavar="NOM.yaml",
        nargs="?",
        const="__AUTO__",  # présent sans valeur → nom auto
        help=(
            "Sauvegarder la config après les questions interactives.\n"
            "  Sans argument : nom automatique (projet_YYYYMMDD.yaml)\n"
            "  Avec argument : chemin spécifié."
        ),
    )
    parser.add_argument(
        "--tsv",
        metavar="report.pg_matrix.tsv",
        default=None,
        help="Chemin du fichier pg_matrix.tsv (surcharge la config).",
    )
    parser.add_argument(
        "--design",
        metavar="ExperimentalDesign.csv",
        default=None,
        help="Chemin du fichier ExperimentalDesign.csv (surcharge la config).",
    )
    parser.add_argument(
        "--out-dir",
        metavar="proteomics_output",
        default=None,
        help="Dossier de sortie (surcharge la config).",
    )
    return parser.parse_args()


# ── Proposition de rechargement de la config précédente ────────────────────────

def maybe_reload_last_config(out_dir: str) -> Optional[dict]:
    """
    Si un last_config.yaml existe dans out_dir, propose à l'utilisateur
    de le recharger pour éviter de rerépondre aux questions.
    Retourne le dict params si accepté, None sinon.
    """
    # Cherche dans out_dir ET dans le répertoire courant
    candidates = [
        os.path.join(out_dir, "last_config.yaml"),
        os.path.join(out_dir, "last_config.json"),
        "last_config.yaml",
        "last_config.json",
    ]
    found = next((p for p in candidates if os.path.exists(p)), None)
    if not found:
        return None

    mtime = datetime.fromtimestamp(os.path.getmtime(found))
    print(f"\n  [RELOAD] Previous config found: {found}")
    print(f"     (generated on {mtime.strftime('%Y-%m-%d %H:%M')})")

    rep = input("  Reload this config? (Y/n) -> ").strip().lower()
    if rep in ("", "o", "oui", "y", "yes"):
        try:
            params = load_config(found)
            print("  [OK] Config reloaded.")
            print_config_summary(params)
            confirm = input("\n  Continue with these parameters? (Y/n) -> ").strip().lower()
            if confirm in ("", "o", "oui", "y", "yes"):
                return params
            else:
                print("  -> Interactive reconfiguration.")
                return None
        except Exception as e:
            print(f"  [WARN] Could not load config ({e}). "
                  f"Interactive reconfiguration.")
            return None
    return None


# ── Point d'entrée principal (appelé depuis main()) ────────────────────────────

def resolve_config(ask_params_fn, ask_go_params_fn,
                   go_available: bool, dash_available: bool,
                   pr_path: str) -> dict:
    """
    Logique centrale de résolution de la configuration :

    1. Parse les args CLI.
    2. Si --config → charger le fichier, compléter les flags optionnels manquants.
    3. Sinon, proposer de recharger last_config.yaml.
    4. Sinon, lancer les questions interactives (ask_params_fn + ask_go_params_fn).
    5. Si --save-config → sauvegarder.
    6. Toujours sauvegarder last_config.yaml dans out_dir.

    Retourne le dict params complet (avec clés dérivées).
    """
    args = parse_cli_args()

    print("\n" + "="*60)
    print("  PROTEOMICS PIPELINE — Statistical threshold configuration")
    print("="*60)

    # -- A. Load from --config --------------------------------------------------
    if args.config:
        print(f"\n  [LOAD] Loading config: {args.config}")
        try:
            params = load_config(args.config)
        except (FileNotFoundError, ValueError, RuntimeError) as e:
            print(f"  [ERROR] {e}")
            sys.exit(1)

        # Surchargez les chemins CLI si fournis
        if args.tsv:     params["tsv_path"]    = args.tsv
        if args.design:  params["design_path"] = args.design
        if args.out_dir: params["out_dir"]     = args.out_dir

        # Les flags optionnels liés à la présence de fichiers/modules
        # peuvent ne pas être dans la config — on les complète
        if "use_deqms" not in params:
            params["use_deqms"] = False
        if "make_wgcna" not in params:
            params["make_wgcna"] = False
        if "make_dashboard" not in params:
            params["make_dashboard"] = dash_available
        if "go_organism" not in params:
            params["go_organism"] = None

        print_config_summary(params)
        return params

    # ── B. Proposition de rechargement last_config ─────────────────────────────
    out_dir = args.out_dir or _DEFAULTS["out_dir"]
    reloaded = maybe_reload_last_config(out_dir)
    if reloaded is not None:
        # Surchargez les chemins CLI si fournis même en cas de rechargement
        if args.tsv:     reloaded["tsv_path"]    = args.tsv
        if args.design:  reloaded["design_path"] = args.design
        if args.out_dir: reloaded["out_dir"]      = args.out_dir
        return reloaded

    # ── C. Questions interactives ──────────────────────────────────────────────
    # ask_params_fn() retourne uniquement les clés stats. On fusionne _DEFAULTS
    # en premier pour garantir tsv_path, design_path, pr_path, out_dir, etc.
    # Les réponses interactives ont priorité via le second spread.
    params = {**_DEFAULTS, **ask_params_fn()}

    # Surchargez les chemins CLI si fournis (priorité maximale)
    if args.tsv:     params["tsv_path"]    = args.tsv
    if args.design:  params["design_path"] = args.design
    if args.out_dir: params["out_dir"]     = args.out_dir

    # GO
    go_organism = None
    if go_available:
        go_result = ask_go_params_fn()
        go_organism = go_result.get("organism") if go_result else None
    params["go_organism"] = go_organism

    # WGCNA
    rep = input("\n  Run the WGCNA co-expression analysis? "
                "(longest step) (y/N) -> ").strip().lower()
    params["make_wgcna"] = rep in ("o", "oui", "y", "yes")

    # Dashboard
    params["make_dashboard"] = False
    if dash_available:
        rep = input("\n  Generate the interactive HTML dashboard at the end? (Y/n) -> "
                    ).strip().lower()
        params["make_dashboard"] = rep in ("", "o", "oui", "y", "yes")

    # DEqMS
    params["use_deqms"] = False
    if os.path.exists(pr_path):
        rep = input(
            f"\n  pr_matrix detected ({pr_path}). ALSO compute DEqMS "
            f"(peptide-count weighting) alongside limma? (Y/n) -> "
        ).strip().lower()
        params["use_deqms"] = rep in ("", "o", "oui", "y", "yes")
    else:
        print(f"\n  [INFO] {pr_path} not found — DEqMS unavailable, limma only.")

    # ── D. Sauvegarde optionnelle --save-config ────────────────────────────────
    if args.save_config is not None:
        if args.save_config == "__AUTO__":
            ts   = datetime.now().strftime("%Y%m%d_%H%M")
            name = f"config_{ts}.yaml"
            save_path = os.path.join(params["out_dir"], name)
        else:
            save_path = args.save_config
        save_config(params, save_path,
                    comment="Généré par --save-config")

    # ── E. Sauvegarde automatique last_config (toujours) ──────────────────────
    _save_last_config(params)

    return params


def _save_last_config(params: dict) -> None:
    """Sauvegarde silencieuse de last_config.yaml dans out_dir."""
    out_dir = params.get("out_dir", _DEFAULTS["out_dir"])
    os.makedirs(out_dir, exist_ok=True)
    ext  = ".yaml" if _YAML_OK else ".json"
    path = os.path.join(out_dir, f"last_config{ext}")
    try:
        save_config(params, path,
                    comment="Sauvegarde automatique — dernier run")
    except Exception as e:
        # Non bloquant : la sauvegarde last_config ne doit jamais planter le pipeline
        print(f"  [WARN] Could not save last_config ({e})")


# ── Exemple de fichier config (affiché avec --help ou --example) ───────────────

CONFIG_EXAMPLE = """\
# Pipeline Protéomique — Exemple de configuration YAML
# Utilisation : python deimos.py --config config_example.yaml

# ── Seuils statistiques ──────────────────────────────────────────────────────
volcano_use_padj:    false        # true = FDR (p.adj), false = p-value brute
volcano_p_thresh:    0.05
volcano_ratio_min:   1.5          # ratio linéaire (ex: 1.5 = 50% de changement)

anova_use_padj:      false
anova_p_thresh:      0.05

n_heatmap_clusters:  3

# ── Options d'analyse ────────────────────────────────────────────────────────
n_iter_robustness:   100          # 0 = désactivé
fdr_global:          false        # true = correction BH sur toutes les p-values
impute_method:       qrilc        # 'qrilc' ou 'mixed'

# ── Modules optionnels ───────────────────────────────────────────────────────
go_organism:         null          # ex: 'bnapus', 'hsapiens', null = désactivé
make_wgcna:          false
make_dashboard:      true
use_deqms:           false

# ── Chemins (optionnel — surchargeables aussi via CLI) ───────────────────────
tsv_path:            report.pg_matrix.tsv
design_path:         ExperimentalDesign.csv
pr_path:             report.pr_matrix.tsv
out_dir:             proteomics_output
"""
