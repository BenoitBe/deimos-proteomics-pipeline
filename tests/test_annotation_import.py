"""Import d'une annotation de collaborateur (Excel Trinotate, background topGO) vers
le format Perseverance / Deimos, et ORA local nommé (BP/CC/MF) sur les têtes de
groupes protéiques."""
import sys
import textwrap
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

import annotation_import as ai                       # noqa: E402
from go_enrichment import (annotation_to_group_leaders, load_go_terms,  # noqa: E402
                           run_local_ora)

OBO = textwrap.dedent("""\
    format-version: 1.2

    [Term]
    id: GO:0008150
    name: biological_process
    namespace: biological_process

    [Term]
    id: GO:0003674
    name: molecular_function
    namespace: molecular_function

    [Term]
    id: GO:0005575
    name: cellular_component
    namespace: cellular_component

    [Term]
    id: GO:0000001
    name: parent process
    namespace: biological_process
    is_a: GO:0008150 ! biological_process

    [Term]
    id: GO:0000002
    name: child process
    namespace: biological_process
    alt_id: GO:0000099
    is_a: GO:0000001 ! parent process
    relationship: part_of GO:0000003 ! other process

    [Term]
    id: GO:0000003
    name: other process
    namespace: biological_process
    is_a: GO:0008150

    [Term]
    id: GO:0000004
    name: regulation of child process
    namespace: biological_process
    intersection_of: GO:0000001 ! parent process
    intersection_of: regulates GO:0000002

    [Term]
    id: GO:0000005
    name: some activity
    namespace: molecular_function
    is_a: GO:0003674
    relationship: part_of GO:0000001 ! cross-ontology, not propagated

    [Term]
    id: GO:0000006
    name: old term
    namespace: biological_process
    is_obsolete: true
    """)


def _trinotate_xlsx(path):
    rows = [
        ["TRINITY_scaffold1_g1_i1", "Q922K7", 64.2, 428, 145, 4, 694, 1974, 190, 610, 0.0, 592.0,
         "x [Gene Symbol: Nop2]", "GO:0000002 [Name: child process]",
         "GO:0005575 [Name: cellular_component]", "GO:0000005 [Name: some activity]",
         "hsa00270; dme01100", "PF01189; PF17125"],
        ["TRINITY_scaffold1_g2_i1", "P19623", 40.0, 100, 50, 2, 1, 300, 1, 300, 1e-30, 200.0,
         "-", "GO:0000004 [Name: regulation of child process]", "-", "-", "-", "PF00087"],
        ["TRINITY_scaffold1_g3_i1", "A7X3R6", 25.0, 37, 24, 0, 343, 233, 32, 68, 3.4, 30.4,
         "-", "GO:0000003 [Name: other process]", "-", "-", "-", "-"],       # e = 3.4
        ["TRINITY_scaffold2_g9_i2", "B1", 50.0, 90, 10, 1, 1, 90, 1, 90, 1e-20, 150.0,
         "-", "GO:0000006 [Name: old term]; GO:0000099 [Name: child process]", "-", "-",
         "-", "-"],
    ]
    cols = ["", "UniProt"] + [f"Unnamed: {i}" for i in range(2, 12)] + [
        "Gene Info", "GO - Biological Process", "GO - Cellular Component",
        "GO - Molecular Function", "KEGG Pathway ID", "Pfam ID"]
    pd.DataFrame(rows, columns=cols).rename(columns={"": "Unnamed: 0"}).to_excel(path, index=False)


def _files(tmp_path):
    x = tmp_path / "Annotation.xlsx"
    _trinotate_xlsx(x)
    bg = tmp_path / "background.txt"
    bg.write_text('ID\tGos\nscaffold1.g1.t1\t"GO:0008150,GO:0000001,GO:0000003"\n'
                  'scaffold1..g7.t1\t"GO:0008150,GO:0000003"\n')
    obo = tmp_path / "go.obo"
    obo.write_text(OBO)
    fa = tmp_path / "collab.faa"
    fa.write_text(">scaffold1.g1.t1 x\nMA\n>scaffold1.g2.t1\nMA\n>scaffold1.g7.t1\nMA\n"
                  ">scaffold2.g9.t2\nMA\n")
    return x, bg, obo, fa


def test_ids_are_normalised():
    assert ai.canonical_id("TRINITY_scaffold1_g1_i1") == ai.canonical_id("scaffold1.g1.t1")
    assert ai.canonical_id("scaffold1..g449.t1") == ai.canonical_id("scaffold1.g449.t1")
    assert ai.canonical_id(" scaffold3.g5.t2") == ai.canonical_id("TRINITY_scaffold3_g5_i2")
    assert ai.default_display_id("TRINITY_scaffold1_g1_i1") == "scaffold1.g1.t1"


def test_obo_parents_follow_go_basic_rules(tmp_path):
    p = tmp_path / "go.obo"
    p.write_text(OBO)
    names, nss, parents, alt = ai.read_obo(str(p))
    assert "GO:0000006" not in names                         # obsolète
    assert alt["GO:0000099"] == "GO:0000002"
    assert parents["GO:0000002"] == {"GO:0000001", "GO:0000003"}    # is_a + part_of
    assert parents["GO:0000004"] == {"GO:0000001"}                  # genre logique
    assert parents["GO:0000005"] == {"GO:0003674"}       # part_of vers BP : pas propagé
    anc = ai.propagate({"GO:0000002"}, parents, alt, {})
    assert anc == {"GO:0000002", "GO:0000001", "GO:0000003", "GO:0008150"}


def test_build_with_evalue_fasta_and_obo(tmp_path):
    x, bg, obo, fa = _files(tmp_path)
    ann, terms, rep = ai.build([str(x), str(bg)], evalue_max=1e-5, fasta=str(fa), obo=str(obo))
    a = ann.set_index("protein_id")
    # g3 : e-value 3.4 -> retiré ; identifiants = ceux du FASTA
    assert set(a.index) == {"scaffold1.g1.t1", "scaffold1.g2.t1", "scaffold1.g7.t1",
                            "scaffold2.g9.t2"}
    g1 = set(a.at["scaffold1.g1.t1", "GO"].split(";"))
    assert {"GO:0000002", "GO:0000001", "GO:0000003", "GO:0008150", "GO:0005575",
            "GO:0000005", "GO:0003674"} == g1
    assert a.at["scaffold1.g1.t1", "KEGG"] == "map00270;map01100"
    assert a.at["scaffold1.g1.t1", "Pfam"] == "PF01189;PF17125"
    assert "background" in a.at["scaffold1.g1.t1", "source"]
    assert set(a.at["scaffold2.g9.t2", "GO"].split(";")) >= {"GO:0000002"}   # alt_id
    assert "GO:0000006" not in a.at["scaffold2.g9.t2", "GO"]                 # obsolète
    assert rep["obsolete_dropped"] == 1
    t = terms.set_index("GO")
    assert t.at["GO:0000004", "name"] == "regulation of child process"
    assert t.at["GO:0000005", "namespace"] == "molecular_function"
    assert (terms.name != "").all()


def test_mixed_depth_warning_without_obo(tmp_path):
    x, bg, _, _ = _files(tmp_path)
    _, terms, rep = ai.build([str(x), str(bg)], evalue_max=1e-5)
    assert "depth_warning" in rep
    # sans OBO : noms et espaces lus dans le tableau Excel
    t = terms.set_index("GO")
    assert t.at["GO:0000002", "name"] == "child process"
    assert t.at["GO:0000002", "namespace"] == "biological_process"


def test_cli(tmp_path):
    x, bg, obo, fa = _files(tmp_path)
    out, tout = tmp_path / "a.tsv", tmp_path / "t.tsv"
    assert ai.main(["--table", str(x), "--table", str(bg), "--obo", str(obo),
                    "--fasta", str(fa), "--out", str(out), "--terms-out", str(tout)]) == 0
    assert list(pd.read_csv(out, sep="\t").columns) == ["protein_id", "GO", "KEGG", "Pfam",
                                                        "source"]
    assert list(pd.read_csv(tout, sep="\t").columns) == ["GO", "name", "namespace"]


def test_local_ora_named_and_split(tmp_path):
    terms = tmp_path / "t.tsv"
    pd.DataFrame({"GO": ["GO:1", "GO:2"], "name": ["alpha process", "beta compartment"],
                  "namespace": ["biological_process", "cellular_component"]}) \
        .to_csv(terms, sep="\t", index=False)
    info = load_go_terms(str(terms))
    bgp = [f"P{i}" for i in range(40)]
    ann = pd.DataFrame({"protein_id": bgp,
                        "GO": ["GO:1;GO:2" if i < 8 else "GO:3" for i in range(40)]})
    res = run_local_ora(bgp[:8], bgp, ann, term_info=info)
    r = res.set_index("native")
    assert r.at["GO:1", "name"] == "alpha process" and r.at["GO:1", "source"] == "GO:BP"
    assert r.at["GO:2", "source"] == "GO:CC"
    res0 = run_local_ora(bgp[:8], bgp, ann)                 # sans table : comportement inchangé
    assert set(res0.source) == {"GO:PERSEVERANCE"} and set(res0.name) == set(res0.native)


def _go_kegg_annotation():
    bgp = [f"P{i}" for i in range(60)]
    ann = pd.DataFrame({
        "protein_id": bgp,
        "GO": ["GO:1;GO:2" if i < 10 else ("GO:3" if i < 50 else None) for i in range(60)],
        # map00900 sur-représenté dans la requête (P0-P9) ; map00010 non
        "KEGG": ["map00900" if i < 9 else ("map00010" if i % 3 == 0 else None)
                 for i in range(60)]})
    return bgp, ann


def test_local_ora_tests_kegg_as_its_own_family(tmp_path):
    terms = tmp_path / "t.tsv"
    pd.DataFrame({"GO": ["GO:1", "map00900"],
                  "name": ["alpha process", "Terpenoid backbone biosynthesis"],
                  "namespace": ["biological_process", "KEGG"]}) \
        .to_csv(terms, sep="\t", index=False)
    info = load_go_terms(str(terms))
    assert info["map00900"] == ("Terpenoid backbone biosynthesis", "KEGG")
    bgp, ann = _go_kegg_annotation()
    res = run_local_ora(bgp[:10], bgp, ann, term_info=info)
    r = res.set_index("native")
    assert r.at["map00900", "source"] == "KEGG"
    assert r.at["map00900", "name"] == "Terpenoid backbone biosynthesis"
    # requête KEGG = P0-P9 annotées en KEGG (P9 porte map00010) ; univers = 26, pas 60
    assert r.at["map00900", "query_size"] == 10 and r.at["map00900", "term_size"] == 9
    assert r.at["GO:1", "source"] == "GO:BP"
    # sans table de noms : source KEGG quand même (vient de la colonne)
    r0 = run_local_ora(bgp[:10], bgp, ann).set_index("native")
    assert r0.at["map00900", "source"] == "KEGG" and r0.at["map00900", "name"] == "map00900"


def test_kegg_orthology_numbers_are_not_tested_as_pathways():
    """annotation_import.py garde les KO (K00001) dans la colonne KEGG : ce ne sont
    pas des voies, ils ne doivent pas former de termes KEGG."""
    bgp, ann = _go_kegg_annotation()
    ann["KEGG"] = ann["KEGG"].map(lambda v: f"K00001;{v}" if isinstance(v, str) else "K00002")
    res = run_local_ora(bgp[:10], bgp, ann)
    kegg = set(res.loc[res.source == "KEGG", "native"])
    assert "map00900" in kegg and not any(k.startswith("K") for k in kegg)
    only_ko = ann.assign(KEGG="K00001")
    assert not (run_local_ora(bgp[:10], bgp, only_ko).source == "KEGG").any()


def test_kegg_column_leaves_go_results_unchanged():
    """Familles BH séparées : la présence du KEGG ne doit changer aucune ligne GO."""
    bgp, ann = _go_kegg_annotation()
    with_kegg = run_local_ora(bgp[:10], bgp, ann)
    go_only = run_local_ora(bgp[:10], bgp, ann.drop(columns="KEGG"))
    pd.testing.assert_frame_equal(
        with_kegg[with_kegg.source != "KEGG"].reset_index(drop=True), go_only)


def test_group_leaders_get_member_kegg():
    df = pd.DataFrame({"name": ["A1", "B1"], "Protein.Group": ["A1;C1", "B1"]})
    ann = pd.DataFrame({"protein_id": ["A1", "C1"], "GO": ["GO:1", "GO:2"],
                        "KEGG": [None, "map00900"]})
    out = annotation_to_group_leaders(ann, df)
    kegg = out.dropna(subset=["KEGG"])
    assert set(zip(kegg.protein_id, kegg.KEGG)) == {("C1", "map00900"), ("A1", "map00900")}
    assert out.dropna(subset=["GO"]).query("protein_id == 'A1'").GO.tolist() == ["GO:1"]


def test_group_leaders_get_member_annotation():
    df = pd.DataFrame({"name": ["XP_1", "XP_3", "XP_5"],
                       "Protein.Group": ["XP_1;XP_2", "XP_3;XP_4", "XP_5"]})
    ann = pd.DataFrame({"protein_id": ["XP_2", "XP_3", "XP_4"],
                        "GO": ["GO:1;GO:2", "GO:9", "GO:8"]})
    out = annotation_to_group_leaders(ann, df)
    by = dict(zip(out.protein_id, out.GO))
    assert by["XP_1"] == "GO:1;GO:2"            # tête sans annotation <- membre annoté
    assert by["XP_3"] == "GO:9"                 # tête annotée : inchangée
    assert "XP_5" not in by
    assert len(out) == len(ann) + 1             # rien n'est retiré


def test_gene_level_annotation_reaches_every_isoform(tmp_path):
    """Une annotation au niveau du gène (scaffold1.g1) doit atteindre toutes les
    isoformes du FASTA (g1.t1, g1.t2), malgré le préfixe de scaffold absent du FASTA."""
    from annotation_import import match_to_fasta
    fa = tmp_path / "prot.fa"
    fa.write_text(">g1.t1\nMK\n>g1.t2\nMKK\n>g2.t1\nMA\n>g7.t1\nMW\n")
    out, rep = match_to_fasta(["scaffold1_g1", "scaffold1_g2"], str(fa))
    assert out["scaffold1_g1"] == ["g1.t1", "g1.t2"]        # les deux isoformes
    assert out["scaffold1_g2"] == ["g2.t1"]
    assert rep["matched_suffix"] == 2 and rep["proteins_annotated"] == 3
    assert rep["matched_exact"] == 0 and rep["ambiguous_dropped"] == 0

    # voie exacte et voie gène, quand les identifiants se correspondent déjà
    out, rep = match_to_fasta(["g1_t1", "g2"], str(fa))
    assert out["g1_t1"] == ["g1.t1"] and rep["matched_exact"] == 1
    assert out["g2"] == ["g2.t1"] and rep["matched_gene"] == 1


def test_an_ambiguous_suffix_is_dropped_not_guessed(tmp_path):
    """Si deux scaffolds portent chacun un g1, le suffixe ne tranche pas : mieux vaut
    ne pas annoter que d'attribuer l'annotation du mauvais gène."""
    from annotation_import import match_to_fasta
    fa = tmp_path / "prot.fa"
    fa.write_text(">g1.t1\nMK\n>g5.t1\nMA\n")
    out, rep = match_to_fasta(["scaffold1_g1", "scaffold2_g1", "scaffold1_g5"], str(fa))
    assert "scaffold1_g1" not in out and "scaffold2_g1" not in out
    assert rep["ambiguous_dropped"] == 2
    assert out["scaffold1_g5"] == ["g5.t1"] and rep["matched_suffix"] == 1


def test_an_evalue_column_is_recognised_by_its_shape(tmp_path):
    """Trinotate nomme son e-value 'p-value'. Une colonne qui descend à 1e-125 et
    plafonne à 10 ne peut être qu'une e-value : ni une identité (0-100), ni un
    bitscore. On la reconnaît à sa forme, pas à son nom."""
    import numpy as np
    import pandas as pd
    from annotation_import import find_evalue_column, looks_like_evalue

    ev = pd.Series(list(np.logspace(-125, 1, 60)) + [0.0] * 5)
    assert looks_like_evalue(ev)
    assert not looks_like_evalue(pd.Series(np.linspace(20, 99, 60)))      # % identité
    assert not looks_like_evalue(pd.Series(np.linspace(50, 2000, 60)))    # bitscore
    assert not looks_like_evalue(pd.Series([0.04, 0.2, 0.9] * 30))        # vraie p-value
    assert not looks_like_evalue(pd.Series([1e-30, 2.0]))                 # trop peu de lignes

    df = pd.DataFrame({"Seq.Name": [f"g{i}" for i in range(65)], "p-value": ev,
                       "Biological": ["GO:0008150 [Name: biological_process]"] * 65})
    col, maybe = find_evalue_column(df)
    assert col == "p-value" and maybe == []                 # retenue d'office
    df2 = df.assign(**{"p-value": np.linspace(20, 99, 65)})
    col2, maybe2 = find_evalue_column(df2)
    assert col2 is None and maybe2 == ["p-value"]           # signalée, pas devinée
