# -*- coding: utf-8 -*-
"""
check_gsea.py — Diagnostic de l'intégration GSEA dans Deimos.
Lance : python check_gsea.py
Vérifie que toutes les pièces sont en place et cohérentes.
"""
import os, sys, importlib, inspect

def ok(m):   print(f"  [OK]   {m}")
def bad(m):  print(f"  [FAIL] {m}")
def info(m): print(f"  [..]   {m}")

print("=" * 60)
print("  DIAGNOSTIC GSEA — Deimos")
print("=" * 60)

here = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, here)

# 1. Fichier présent
print("\n[1] Fichier gsea_enrichment.py")
p = os.path.join(here, "gsea_enrichment.py")
if os.path.isfile(p):
    ok(f"présent ({os.path.getsize(p)} octets)")
else:
    bad("ABSENT du dossier — c'est la cause.")
    sys.exit(1)

# 2. gseapy importable
print("\n[2] Dépendance gseapy")
try:
    import gseapy
    ok(f"gseapy {gseapy.__version__} importable")
except ImportError:
    bad("gseapy non installé -> pip install gseapy")

# 3. Fonctions du module
print("\n[3] Fonctions de gsea_enrichment")
try:
    import gsea_enrichment as gs
    for fn in ("run_gsea", "export_gsea_sheets", "_ranking_metric"):
        if hasattr(gs, fn):
            ok(f"{fn} présente")
        else:
            bad(f"{fn} MANQUANTE")
except Exception as e:
    bad(f"import échoué : {e}")

# 4. ask_go_params pose la question
print("\n[4] Question GSEA dans go_enrichment.ask_go_params")
try:
    import go_enrichment as ge
    src = inspect.getsource(ge.ask_go_params)
    if "run_gsea" in src and "Also run GSEA" in src:
        ok("la question GSEA est présente")
    else:
        bad("ask_go_params ne pose PAS la question GSEA -> go_enrichment.py obsolète")
except Exception as e:
    bad(f"lecture impossible : {e}")

# 5. deimos.py appelle run_gsea  <-- LE POINT CRITIQUE
print("\n[5] Appel de run_gsea dans deimos.py  (POINT CRITIQUE)")
dp = os.path.join(here, "deimos.py")
if os.path.isfile(dp):
    with open(dp, encoding="utf-8") as f:
        dsrc = f.read()
    checks = {
        "import gsea_enrichment": "from gsea_enrichment import run_gsea" in dsrc,
        "condition run_gsea":     'go_params.get("run_gsea")' in dsrc,
        "appel run_gsea(":        "run_gsea(" in dsrc,
        "export_gsea_sheets(":    "export_gsea_sheets(" in dsrc,
    }
    for k, v in checks.items():
        (ok if v else bad)(k)
    if not all(checks.values()):
        print("\n  >>> CAUSE PROBABLE : ton deimos.py ne contient pas le bloc")
        print("      d'appel GSEA. Remplace-le par la version à jour.")
    else:
        ok("deimos.py intègre bien l'appel GSEA")
else:
    bad("deimos.py introuvable")

print("\n" + "=" * 60)
print("  Fin du diagnostic")
print("=" * 60)
