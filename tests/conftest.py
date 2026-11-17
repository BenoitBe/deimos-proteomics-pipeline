"""Rend les modules de Deimos importables quel que soit le mode d'invocation.

Les tests font `from deimos import ...`, et deimos.py est à la racine du dépôt.
`python -m pytest` depuis la racine ajoute le répertoire courant à sys.path et
cela fonctionne ; le script `pytest` ne le fait pas, et toute la suite échoue
sur ModuleNotFoundError. C'est ainsi que la première exécution de la CI est
tombée alors que la suite passait en local.

Même mécanisme que tests/conftest.py de Pathfinder.
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
