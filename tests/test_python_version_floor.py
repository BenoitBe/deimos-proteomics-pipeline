"""Compatibilite de version : tout doit compiler avec le Python le plus ancien."""
import pytest


def test_every_module_compiles_on_the_oldest_supported_python():
    """Deimos tourne chez l'utilisateur en 3.11. Depuis 3.12 (PEP 701) une f-string peut contenir un
    saut de ligne DANS son expression, et des guillemets identiques imbriqués : écrit
    sur une machine récente, le fichier ne compile plus chez quelqu'un en 3.11, et
    l'erreur ne se voit qu'à l'exécution (un onglet entier disparaît du rapport).
    On compile donc tout avec l'interpréteur le plus ancien disponible."""
    import pathlib
    import shutil
    import subprocess
    import sys

    older = next((p for v in ("3.9", "3.10", "3.11")
                  if (p := shutil.which(f"python{v}"))), None)
    if older is None or sys.version_info[:2] <= (3, 11):
        pytest.skip("no older interpreter than this one available")
    root = pathlib.Path(__file__).resolve().parent.parent
    files = sorted(str(p) for p in root.rglob("*.py")
                   if ".git" not in p.parts and "__pycache__" not in p.parts)
    assert files
    script = ("import sys\n"
              "bad = []\n"
              "for f in sys.argv[1:]:\n"
              "    src = open(f, encoding='utf-8', errors='replace').read()\n"
              "    try:\n"
              "        compile(src, f, 'exec')\n"
              "    except SyntaxError as e:\n"
              "        bad.append(f'{f}:{e.lineno}: {e.msg}')\n"
              "print('\\n'.join(bad))\n"
              "sys.exit(1 if bad else 0)\n")
    r = subprocess.run([older, "-c", script, *files], capture_output=True, text=True)
    assert r.returncode == 0, f"ne compile pas avec {older} :\n{r.stdout}{r.stderr}"
