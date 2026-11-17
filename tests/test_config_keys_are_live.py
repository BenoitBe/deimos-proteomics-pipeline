"""Aucune clé de configuration ne doit être lue puis rester sans effet.

`batch_column` était analysée par config.py, exposée en `--batch-column` et
documentée dans le README, et n'entrait dans aucun modèle : un run qui
déclarait un lot produisait exactement le même résultat qu'un run sans lot,
sans un mot. Compter les occurrences ne l'aurait pas trouvée — la clé
apparaissait bien dans le code, dans la plomberie de la ligne de commande.

Ce module suit la valeur sur un pas depuis son point de lecture :
  - affectée à une variable  -> cette variable doit être relue ;
  - passée en argument nommé -> la fonction appelée doit lire ce paramètre ;
  - utilisée dans une expression, un test, un retour -> consommée.

Un pas suffit pour ce défaut : ce qui est mort l'est au point de lecture. Les
clés dont la valeur voyage plus loin avant de mourir ne sont attrapables que
par un run complet — d'où le job `smoke` de la CI.

Le détecteur se valide lui-même (test_le_detecteur_attrape_une_cle_inerte) :
un détecteur silencieusement cassé est pire que pas de détecteur.
"""
import ast
import pathlib
import textwrap

import pytest

ROOT = pathlib.Path(__file__).resolve().parents[1]


# ─────────────────────────── le détecteur ───────────────────────────────────

def _key_read_by(node):
    """La clé de configuration lue par ce nœud, ou None."""
    if (isinstance(node, ast.Call) and isinstance(node.func, ast.Attribute)
            and node.func.attr == "get" and node.args
            and isinstance(node.args[0], ast.Constant)
            and isinstance(node.args[0].value, str)):
        return node.args[0].value
    if (isinstance(node, ast.Subscript) and isinstance(node.slice, ast.Constant)
            and isinstance(node.slice.value, str)):
        return node.slice.value
    return None


def _parents(tree):
    out = {}
    for n in ast.walk(tree):
        for ch in ast.iter_child_nodes(n):
            out[id(ch)] = n
    return out


def _enclosing_function(tree, node):
    best = None
    for f in ast.walk(tree):
        if (isinstance(f, (ast.FunctionDef, ast.AsyncFunctionDef))
                and f.lineno <= node.lineno <= (f.end_lineno or f.lineno)):
            if best is None or f.lineno > best.lineno:
                best = f
    return best


def _names_loaded(scope):
    return {n.id for n in ast.walk(scope)
            if isinstance(n, ast.Name) and isinstance(n.ctx, ast.Load)}


def verdict(source: str, key: str):
    """('consommee'|'inerte'|'absente', détail) pour une clé dans un source."""
    tree = ast.parse(source)
    par = _parents(tree)
    funcs = {n.name: n for n in ast.walk(tree)
             if isinstance(n, (ast.FunctionDef, ast.AsyncFunctionDef))}
    reads = [n for n in ast.walk(tree) if _key_read_by(n) == key]
    if not reads:
        return "absente", "aucune lecture"

    for node in reads:
        parent = par.get(id(node))
        # affectation à une variable : la variable doit être relue
        if (isinstance(parent, ast.Assign) and len(parent.targets) == 1
                and isinstance(parent.targets[0], ast.Name)):
            var = parent.targets[0].id
            scope = _enclosing_function(tree, node) or tree
            later = [n for n in ast.walk(scope)
                     if isinstance(n, ast.Name) and n.id == var
                     and isinstance(n.ctx, ast.Load) and n.lineno > node.lineno]
            if later:
                return "consommee", f"{var} relue L{later[0].lineno}"
            continue
        # argument nommé : la fonction appelée doit lire ce paramètre
        if isinstance(parent, ast.keyword):
            call = par.get(id(parent))
            fn_name = (getattr(getattr(call, "func", None), "id", None)
                       or getattr(getattr(call, "func", None), "attr", None))
            fn = funcs.get(fn_name)
            if fn is None:
                return "consommee", f"passée à {fn_name}() (hors module)"
            if parent.arg in _names_loaded(fn):
                return "consommee", f"{fn_name}() lit {parent.arg}"
            return "inerte", f"{fn_name}() ne lit jamais {parent.arg}"
        # tout autre contexte : la valeur participe à une expression
        if parent is not None and not isinstance(parent, ast.Expr):
            return "consommee", type(parent).__name__
    return "inerte", "lue puis jamais relue"


# ─────────────────── validation du détecteur lui-même ───────────────────────

INERT_SAMPLES = {
    "cle_affectee_jamais_relue": textwrap.dedent("""
        def run(params):
            mort = params.get("cle_affectee_jamais_relue")
    """),
    "cle_passee_a_callee_qui_ignore": textwrap.dedent("""
        def callee(seuil=None, ignore=None):
            print(seuil)
        def run(params):
            callee(ignore=params.get("cle_passee_a_callee_qui_ignore"))
    """),
    "cle_jamais_lue": "def run(params):\n    pass\n",
}


@pytest.mark.parametrize("key,src", sorted(INERT_SAMPLES.items()))
def test_le_detecteur_attrape_une_cle_inerte(key, src):
    """Sans ceci, un détecteur cassé rendrait un quitus silencieux."""
    assert verdict(src, key)[0] in ("inerte", "absente")


def test_le_detecteur_ne_crie_pas_sur_une_cle_vivante():
    src = textwrap.dedent("""
        def run(params):
            seuil = params.get("vivante")
            if seuil > 1:
                print(seuil)
    """)
    assert verdict(src, "vivante")[0] == "consommee"


def test_le_detecteur_reconnait_le_defaut_historique():
    """La forme exacte de batch_column après correctif."""
    src = textwrap.dedent("""
        def apply_batch_covariate(design_mat, all_cols, conditions, design, params):
            raw = params.get("batch_column")
            if not raw:
                return design_mat, all_cols
    """)
    assert verdict(src, "batch_column")[0] == "consommee"


# ─────────────────────────── l'audit du dépôt ───────────────────────────────

def _sources():
    return [p for p in ROOT.rglob("*.py")
            if ".git" not in p.parts and "__pycache__" not in p.parts
            and "tests" not in p.parts and p.name != "config.py"]


def test_aucune_cle_de_configuration_n_est_inerte():
    from config import _DEFAULTS

    sources = {p: p.read_text(encoding="utf-8", errors="replace") for p in _sources()}
    assert sources, "aucun module à auditer : le test ne prouverait rien"

    suspectes = {}
    for key in sorted(_DEFAULTS):
        verdicts = []
        for p, src in sources.items():
            try:
                verdicts.append((p.name,) + verdict(src, key))
            except SyntaxError:
                continue
        if any(v[1] == "consommee" for v in verdicts):
            continue
        inerte = next((v for v in verdicts if v[1] == "inerte"), None)
        suspectes[key] = (f"{inerte[0]}: {inerte[2]}" if inerte
                          else "lue nulle part hors config.py")

    assert not suspectes, (
        "clé(s) de configuration déclarée(s) mais sans effet — "
        "c'est le défaut de batch_column :\n"
        + "\n".join(f"  {k}: {d}" for k, d in sorted(suspectes.items())))


DEAD_FLAG_SAMPLE = textwrap.dedent("""
    import argparse
    def main():
        ap = argparse.ArgumentParser()
        ap.add_argument("--utilise")
        ap.add_argument("--par-getattr", dest="par_getattr", action="store_true")
        ap.add_argument("--mort")
        a = ap.parse_args()
        print(a.utilise, getattr(a, "par_getattr", False))
""")


def test_le_detecteur_de_drapeaux_distingue_vivant_et_mort():
    """Les deux formes de lecture doivent compter ; seul --mort doit sortir."""
    flags, reads = _flags_and_reads({pathlib.Path("echantillon.py"): DEAD_FLAG_SAMPLE})
    mortes = {d for d in flags if d not in reads}
    assert mortes == {"mort"}, mortes


DEAD_FLAG_SAMPLE = textwrap.dedent("""
    import argparse
    def main():
        ap = argparse.ArgumentParser()
        ap.add_argument("--utilise")
        ap.add_argument("--par-getattr", dest="par_getattr", action="store_true")
        ap.add_argument("--mort")
        a = ap.parse_args()
        print(a.utilise, getattr(a, "par_getattr", False))
""")


def test_le_detecteur_de_drapeaux_distingue_vivant_et_mort():
    """Les deux formes de lecture doivent compter ; seul --mort doit sortir."""
    flags, reads = _flags_and_reads({"echantillon.py": DEAD_FLAG_SAMPLE})
    assert {d for d in flags if d not in reads} == {"mort"}


def _flags_and_reads(sources: dict):
    """(drapeaux déclarés -> emplacement, noms lus quelque part)."""
    flags, reads = {}, set()
    for p, src in sources.items():
        try:
            tree = ast.parse(src)
        except SyntaxError:
            continue
        for n in ast.walk(tree):
            if (isinstance(n, ast.Call) and isinstance(n.func, ast.Attribute)
                    and n.func.attr == "add_argument"):
                dest = next((k.value.value for k in n.keywords if k.arg == "dest"
                             and isinstance(k.value, ast.Constant)), None)
                if dest is None:
                    longs = [a.value for a in n.args if isinstance(a, ast.Constant)
                             and str(a.value).startswith("--")]
                    if not longs:
                        continue
                    dest = longs[0][2:].replace("-", "_")
                flags[dest] = f"{getattr(p, 'name', p)}:{n.lineno}"
        # Lectures sous leurs DEUX formes. Deimos lit la plupart de ses drapeaux
        # par getattr(args, "qc_only", False) et non par args.qc_only : ne
        # regarder que les attributs déclarait trois drapeaux vivants morts.
        for n in ast.walk(tree):
            if isinstance(n, ast.Attribute) and isinstance(n.ctx, ast.Load):
                reads.add(n.attr)
            if (isinstance(n, ast.Call) and isinstance(n.func, ast.Name)
                    and n.func.id == "getattr" and len(n.args) >= 2
                    and isinstance(n.args[1], ast.Constant)
                    and isinstance(n.args[1].value, str)):
                reads.add(n.args[1].value)
    return flags, reads


def test_aucune_option_de_ligne_de_commande_n_est_morte():
    """Un drapeau annoncé dans --help et jamais relu est le même défaut."""
    sources = {p: p.read_text(encoding="utf-8", errors="replace")
               for p in ROOT.rglob("*.py")
               if ".git" not in p.parts and "__pycache__" not in p.parts
               and "tests" not in p.parts}
    flags, reads = _flags_and_reads(sources)
    assert flags, "aucune option trouvée : le test ne prouverait rien"
    mortes = {d: loc for d, loc in flags.items() if d not in reads}
    assert not mortes, (
        "option(s) de ligne de commande jamais relue(s) :\n"
        + "\n".join(f"  --{d.replace('_', '-')} ({loc})"
                    for d, loc in sorted(mortes.items())))
