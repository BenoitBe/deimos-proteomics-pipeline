"""id_mapping : conversion des identifiants des matrices DIA-NN (pg + pr) par le
FASTA de recherche ou une table, et branchement dans deimos.main()."""
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from id_remap import apply_id_mapping  # noqa: E402

PG = ("Protein.Group\tProtein.Ids\tProtein.Names\tGenes\tFirst.Protein.Description\t"
      "N.Sequences\tN.Proteotypic.Sequences\tSple-1\tSple-2\n"
      "UPI0000000001\tUPI0000000001\tUPI0000000001\t\t\t3\t3\t1000\t2000\n"
      "UPI0000000002;UPI0000000003\tUPI0000000002;UPI0000000003\t\t\t\t2\t1\t500\t600\n"
      "cRAP-001\tcRAP-001\tcRAP-001\t\t\t1\t1\t10\t20\n")
PR = ("Protein.Group\tProtein.Ids\tProtein.Names\tGenes\tFirst.Protein.Description\t"
      "Proteotypic\tStripped.Sequence\tModified.Sequence\tPrecursor.Charge\tPrecursor.Id\t"
      "Sple-1\tSple-2\n"
      "UPI0000000001\tUPI0000000001\t\t\t\t1\tPEPTIDEK\tPEPTIDEK\t2\tPEPTIDEK2\t100\t200\n")
FASTA = (">UPI0000000001 Protein one OS=Foo OX=1 GN=abcA AC=P11111\nMAAA\n"
         ">UPI0000000002 Protein two OS=Foo OX=1 GN=abcB AC=P22222\nMBBB\n"
         ">UPI0000000003 Protein three OS=Foo OX=1 AC=P33333\nMCCC\n")


def test_fasta_mapping_pg_and_pr(tmp_path):
    (tmp_path / "pg.tsv").write_text(PG)
    (tmp_path / "pr.tsv").write_text(PR)
    (tmp_path / "s.fasta").write_text(FASTA)
    params = {"id_mapping": {"path": str(tmp_path / "s.fasta"), "tag": "FOO"}}
    t, r = apply_id_mapping(params, str(tmp_path / "pg.tsv"), str(tmp_path / "pr.tsv"),
                            str(tmp_path / "out"))
    pg, pr = pd.read_csv(t, sep="\t"), pd.read_csv(r, sep="\t")
    assert list(pg["Protein.Group"]) == ["P11111", "P22222;P33333", "cRAP-001"]
    assert pr["Protein.Group"].iloc[0] == "P11111"             # pg et pr appariables
    assert pg["Sple-1"].tolist() == [1000, 500, 10]             # intensités intactes
    assert "2/3" in params["id_mapping_summary"]
    # idempotent : relancé sur les fichiers convertis, ne change rien
    t2, _ = apply_id_mapping(params, t, r, str(tmp_path / "out2"))
    assert pd.read_csv(t2, sep="\t")["Protein.Group"].tolist() == pg["Protein.Group"].tolist()


def test_deimos_calls_id_mapping():
    """Garde-fou : le branchement dans deimos.py (perdu une fois lors d'un
    remplacement de fichier) est présent, avec la ligne Methods."""
    src = (ROOT / "deimos.py").read_text(encoding="utf-8")
    assert "from id_remap import apply_id_mapping" in src
    assert "tsv_path, pr_path = apply_id_mapping(params, tsv_path, pr_path, out_dir)" in src
    assert 'params.get("id_mapping_summary")' in src
