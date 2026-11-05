"""Identifiants d'orthologues : en-têtes UniProt (sp|ACC|NAME) ramenés à l'accession,
seule forme reconnue par gProfiler."""
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
rbh = pytest.importorskip("rbh")


@pytest.mark.parametrize("raw, clean", [
    ("sp|Q9C5R8|TEST1_ARATH", "Q9C5R8"),
    ("tr|F4HXX1|F4HXX1_ARATH", "F4HXX1"),
    ("A0ABQ8DMI4", "A0ABQ8DMI4"),
    ("AT1G01010.1", "AT1G01010.1"),
    ("BnaA01g00010D", "BnaA01g00010D"),
    ("gi|12345|ref|XP_1|", "gi|12345|ref|XP_1|"),
])
def test_clean_seqid(raw, clean):
    assert rbh.clean_seqid(raw) == clean
