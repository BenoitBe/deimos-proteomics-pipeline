"""Genes rempli depuis le FASTA de recherche quand DIA-NN l'a laissé vide."""
import pandas as pd
import pytest

from deimos import (annotate_from_fasta, fill_annotation_from_fasta,
                    find_annotation_fasta, parse_fasta_annotation)

ARAPORT = """>AT1G01010.1 | Symbols: NAC001, ANAC001, NTL10 | NAC domain containing protein 1 | chr1:3760-5630
MEDQVGFGFRPND
>AT1G01020.1 | Symbols:  | ARV1 family protein | chr1:6915-8666 REVERSE
MAASEHRCVGCGF
>AT1G01050.1 | Symbols: PPa1 | pyrophosphorylase 1 | chr1:31170-33153
MSEEGKLNVGKL
"""

UNIPROT = """>sp|P12345|AATM_ARATH Aspartate aminotransferase OS=Arabidopsis thaliana OX=3702 GN=GOT2 PE=1 SV=3
MALLHSARVLS
>tr|Q9ZZZ9|Q9ZZZ9_ARATH Uncharacterized protein OS=Arabidopsis thaliana GN=At1g99999 PE=4 SV=1
MKTAYIAKQR
"""

PLAIN = """>contig_17_g2 hypothetical protein, no symbol
MKTAYIAKQRQ
"""


def _fasta(tmp_path, text, name="query.fasta"):
    p = tmp_path / name
    p.write_text(text)
    return str(p)


def test_parse_the_three_header_conventions(tmp_path):
    a = parse_fasta_annotation(_fasta(tmp_path, ARAPORT))
    assert a["AT1G01010.1"] == ("NAC001", "NAC domain containing protein 1")
    assert a["AT1G01020.1"] == ("", "ARV1 family protein")      # Symbols vide, pas de gene
    assert a["AT1G01010"][0] == "NAC001"                        # indexé sans l'isoforme

    u = parse_fasta_annotation(_fasta(tmp_path, UNIPROT, "u.fasta"))
    assert u["P12345"] == ("GOT2", "Aspartate aminotransferase")
    assert u["Q9ZZZ9"][0] == "At1g99999"

    p = parse_fasta_annotation(_fasta(tmp_path, PLAIN, "p.fasta"))
    assert p["contig_17_g2"] == ("", "hypothetical protein, no symbol")
    assert parse_fasta_annotation(str(tmp_path / "absent.fasta")) == {}


def test_only_empty_cells_are_filled(tmp_path):
    tsv = pd.DataFrame({
        "Protein.Group": ["AT1G01010.1", "AT1G01020.1", "AT1G01050.1;AT1G01010.1", "AT9G99999.1"],
        "Genes": ["", "|", "DEJA_LA", ""],
        "First.Protein.Description": ["", "", "", ""]})
    done = fill_annotation_from_fasta(tsv, _fasta(tmp_path, ARAPORT))
    assert done["genes"] == 1                                   # seule la 1re gagne un gène
    assert done["matched"] == 3                                 # 3 lignes sur 4 reconnues
    assert done["descriptions"] == 3                            # les 3 descriptions vides
    assert list(tsv["Genes"]) == ["NAC001", "|", "DEJA_LA", ""]  # DIA-NN jamais écrasé
    assert tsv["First.Protein.Description"][1] == "ARV1 family protein"
    assert tsv["First.Protein.Description"][3] == ""            # inconnu du FASTA : vide


def test_auto_waits_for_a_mostly_empty_column(tmp_path, capsys):
    def capsys_text():
        t = capsys.readouterr().out
        print(t, end="")
        return t

    f = _fasta(tmp_path, ARAPORT)
    full = pd.DataFrame({"Protein.Group": ["AT1G01010.1", "AT1G01050.1"],
                         "Genes": ["X", "Y"], "First.Protein.Description": ["", ""]})
    assert annotate_from_fasta(full, f, "__auto__") is None     # DIA-NN a déjà les noms
    assert list(full["Genes"]) == ["X", "Y"]
    assert annotate_from_fasta(full, f, f) == f                 # un chemin force la lecture
    assert list(full["Genes"]) == ["X", "Y"]                    # mais ne remplit que le vide
    assert "0 gene names and 2 descriptions" in capsys_text()   # dit exactement ce qu'il a fait
    assert full["First.Protein.Description"][0] == "NAC domain containing protein 1"

    empty = pd.DataFrame({"Protein.Group": ["AT1G01010.1", "AT1G01050.1"],
                          "Genes": ["", ""], "First.Protein.Description": ["", ""]})
    assert annotate_from_fasta(empty, f, "__auto__") == f
    assert list(empty["Genes"]) == ["NAC001", "PPa1"]

    for off in (False, None, "false"):
        e2 = pd.DataFrame({"Protein.Group": ["AT1G01010.1"], "Genes": [""],
                           "First.Protein.Description": [""]})
        assert annotate_from_fasta(e2, f, off) is None and e2["Genes"][0] == ""


def test_the_fasta_is_found_next_to_the_pg_matrix(tmp_path):
    pg = tmp_path / "report.pg_matrix.tsv"
    pg.write_text("x")
    assert find_annotation_fasta(str(pg)) is None
    notfasta = tmp_path / "notes.fa"
    notfasta.write_text("ceci n'est pas un FASTA\n")             # pas de '>' : écarté
    assert find_annotation_fasta(str(pg)) is None
    small = tmp_path / "panel.fa"
    small.write_text(ARAPORT)                                   # 300 octets, mais valide
    assert find_annotation_fasta(str(pg)) == str(small)
    bigger = tmp_path / "reference.fasta"
    bigger.write_text(ARAPORT + "\n>x | Symbols: Z | z |\n" + "M" * 5000)
    assert find_annotation_fasta(str(pg)) == str(bigger)        # le plus gros gagne
    assert find_annotation_fasta(str(pg), str(small)) == str(small)
    assert find_annotation_fasta(str(pg), str(notfasta)) is None
    assert find_annotation_fasta(str(tmp_path / "x.tsv"), "/nope.fasta") is None
    # le FASTA n'est que lu : rien n'est réécrit
    assert bigger.read_text().startswith(">AT1G01010.1")


def test_no_identifier_matches_leaves_the_column_alone(tmp_path, capsys):
    """Le FASTA est lisible mais ne parle pas des mêmes protéines : on prévient et on
    ne touche à rien (cas d'un FASTA d'une autre espèce laissé dans le dossier)."""
    pg = tmp_path / "report.pg_matrix.tsv"
    pg.write_text("x")
    _fasta(tmp_path, ARAPORT)
    tsv = pd.DataFrame({"Protein.Group": ["FOO1", "BAR2"], "Genes": ["", ""],
                        "First.Protein.Description": ["", ""]})
    assert annotate_from_fasta(tsv, str(pg), "__auto__") is None
    assert list(tsv["Genes"]) == ["", ""]
    assert "No identifier" in capsys.readouterr().out


def test_protein_ids_column_is_preferred(tmp_path):
    """Protein.Group peut être tronqué ; Protein.Ids porte tous les accessions."""
    tsv = pd.DataFrame({"Protein.Group": ["AT9G99999.1"],
                        "Protein.Ids": ["AT9G99999.1;AT1G01050.1"],
                        "Genes": [""], "First.Protein.Description": [""]})
    assert fill_annotation_from_fasta(tsv, _fasta(tmp_path, ARAPORT))["genes"] == 1
    assert tsv["Genes"][0] == "PPa1"
