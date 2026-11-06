#!/usr/bin/env python3
# ==============================================================================
# annotation_import.py — annotation fournie par un collaborateur -> format Deimos
# Deimos — DIA Expression Integrated Multi-Omics Suite
# ==============================================================================
# Convertit une annotation fonctionnelle livrée sous forme de tableau (Excel
# Trinotate/Blast2GO, CSV, TSV, fichier « background » topGO : ID <tab> GO,GO,...)
# en deux tables :
#   - <out>          : protein_id, GO, KEGG, Pfam, source  (GO ';'-séparés)
#                      = format --ref-annotations de Perseverance et
#                        perseverance_annotation_path de Deimos ;
#   - <terms-out>    : GO, name, namespace  (noms lus dans le tableau, '[Name: ...]',
#                      ou dans go-basic.obo) = perseverance_go_terms_path de Deimos,
#                      pour des résultats GO nommés et séparés en BP / CC / MF.
#
# Les GO sont repérés partout (motif GO:0000000), quel que soit le format des
# cellules. L'espace (BP/CC/MF) vient de l'en-tête de colonne quand il le dit, ou
# de l'OBO. Les identifiants sont normalisés (TRINITY_scaffold1_g1_i1 ~
# scaffold1.g1.t1) et, avec --fasta, réécrits exactement comme dans le FASTA.
#
# Usage (Ostrea edulis : annotation sur les modèles de gènes du collaborateur) :
#   python annotation_import.py --table AnnotationOstreaedulis.xlsx \
#       --table "background file for GOostrea.txt" --evalue-max 1e-5 \
#       --fasta proteines_collaborateur.fasta --obo go-basic.obo \
#       --out ostrea_ref_annotation.tsv --terms-out ostrea_go_terms.tsv
# puis Perseverance (identifiants du FASTA de recherche <- modèles annotés) :
#   python perseverance.py --query GCF_947568905.1_protein.faa \
#       --reference proteines_collaborateur.fasta \
#       --ref-annotations ostrea_ref_annotation.tsv --out ostrea_annotation.tsv \
#       --min-identity 90 --min-coverage 70
# ==============================================================================
from __future__ import annotations

import argparse
import os
import re
import sys

import pandas as pd

GO_RE = re.compile(r"GO:\d{7}")
NAMED_GO_RE = re.compile(r"(GO:\d{7})\s*\[\s*Name:\s*([^\]]+?)\s*\]")
KEGG_RE = re.compile(r"\b(?:ko:)?K\d{5}\b|\b[a-z]{2,4}\d{5}\b")


def _kegg(cell) -> set:
    """KO (K00001) gardés ; voies propres à une espèce (hsa00270, dme01100) ramenées
    à la carte générique (map00270), comparable quelle que soit l'espèce du hit."""
    out = set()
    for k in KEGG_RE.findall(str(cell)):
        k = k.replace("ko:", "")
        out.add(k if k.startswith("K") else "map" + k[-5:])
    return out
PFAM_RE = re.compile(r"\bPF\d{5}\b")
NS = {"biological_process": "biological_process", "cellular_component": "cellular_component",
      "molecular_function": "molecular_function"}
_NS_FROM_HEADER = [(re.compile(r"biolog|\bBP\b|process", re.I), "biological_process"),
                   (re.compile(r"cellul|\bCC\b|component|compartment", re.I), "cellular_component"),
                   (re.compile(r"molec|\bMF\b|function", re.I), "molecular_function")]
_EVALUE_HEADERS = re.compile(r"^(e[-_ ]?value|evalue|e)$", re.I)
_BLAST12 = ["qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart",
            "qend", "sstart", "send", "evalue", "bitscore"]


# ── Identifiants ──────────────────────────────────────────────────────────────

def canonical_id(x) -> str:
    """Clé de rapprochement : TRINITY_scaffold1_g1_i1, scaffold1.g1.t1,
    scaffold1_g1_t1 -> 'scaffold1_g1_t1' (casse, séparateurs, préfixe Trinity,
    suffixe d'isoforme _iN / .tN / -RA... ignorés dans leur forme)."""
    s = str(x).strip().split()[0] if str(x).strip() else ""
    s = s.split("|")[1] if s.count("|") >= 2 and s.split("|")[0] in ("sp", "tr") else s
    s = re.sub(r"^TRINITY_", "", s, flags=re.I)
    s = re.sub(r"[.\-_]+", "_", s).lower().strip("_")   # 'scaffold1..g449.t1' (coquille)
    s = re.sub(r"_i(\d+)$", r"_t\1", s)
    return s


def default_display_id(x) -> str:
    """Forme AUGUSTUS/BRAKER (scaffold1.g1.t1) quand l'ID en a la structure,
    sinon l'ID tel quel (sans préfixe TRINITY_)."""
    s = re.sub(r"^TRINITY_", "", str(x).strip(), flags=re.I)
    m = re.match(r"^(.+?)[_.]+(g\d+)[_.]+(?:i|t)(\d+)$", s)
    return f"{m.group(1)}.{m.group(2)}.t{m.group(3)}" if m else s


def fasta_ids(path: str) -> list:
    out = []
    with open(path, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line.startswith(">"):
                tok = line[1:].split()[0] if line[1:].split() else ""
                parts = tok.split("|")
                out.append(parts[1] if len(parts) >= 3 and parts[0] in ("sp", "tr") else tok)
    return out


# ── Lecture des tableaux ──────────────────────────────────────────────────────

def read_table(path: str) -> pd.DataFrame:
    low = path.lower()
    if low.endswith((".xlsx", ".xlsm", ".xls")):
        return pd.read_excel(path, dtype=str)
    with open(path, encoding="utf-8", errors="replace") as fh:
        head = fh.readline()
    sep = "\t" if "\t" in head else ("," if head.count(",") > head.count(";") else ";")
    return pd.read_csv(path, sep=sep, dtype=str, quotechar='"')


def _blast_like(df: pd.DataFrame) -> pd.DataFrame:
    """Tableau Trinotate/Excel sans en-têtes pour les 12 colonnes BLAST
    ('Unnamed: 2'...) : les renomme d'après leur position quand le motif colle."""
    cols = list(df.columns)
    unnamed = [i for i, c in enumerate(cols) if str(c).startswith("Unnamed")]
    if len(cols) >= 12 and len(unnamed) >= 8:
        ren = {cols[i]: _BLAST12[i] for i in range(12) if str(cols[i]).startswith("Unnamed")
               or i == 0}
        return df.rename(columns=ren)
    return df


def parse_table(df: pd.DataFrame, evalue_max=None, id_col=None):
    """-> (records {id_brut: {GO:set, KEGG:set, Pfam:set}}, names {GO: name},
    namespaces {GO: ns}, stats)."""
    df = _blast_like(df)
    id_col = id_col or df.columns[0]
    stats = {"rows": len(df), "dropped_evalue": 0}
    ecol = next((c for c in df.columns if _EVALUE_HEADERS.match(str(c).strip())), None)
    if evalue_max is not None and ecol is not None:
        e = pd.to_numeric(df[ecol], errors="coerce")
        keep = e.le(float(evalue_max)).fillna(False)       # pas de hit -> pas de transfert
        stats["dropped_evalue"] = int((~keep).sum())
        df = df[keep]
    go_cols = [c for c in df.columns if c != id_col
               and df[c].astype(str).str.contains("GO:", regex=False).any()]
    kegg_cols = [c for c in df.columns if re.search(r"kegg|\bko\b", str(c), re.I)]
    pfam_cols = [c for c in df.columns if re.search(r"pfam", str(c), re.I)]
    names, nss, rec = {}, {}, {}
    col_ns = {}
    for c in go_cols:
        for rx, ns in _NS_FROM_HEADER:
            if rx.search(str(c)) and not re.search(r"^gos?$", str(c).strip(), re.I):
                col_ns[c] = ns
                break
    for _, r in df.iterrows():
        pid = str(r[id_col]).strip()
        if not pid or pid.lower() == "nan":
            continue
        e = rec.setdefault(pid, {"GO": set(), "KEGG": set(), "Pfam": set()})
        for c in go_cols:
            cell = str(r[c])
            ids = GO_RE.findall(cell)
            e["GO"].update(ids)
            for g, nm in NAMED_GO_RE.findall(cell):
                names.setdefault(g, nm)
            if c in col_ns:
                for g in ids:
                    nss.setdefault(g, col_ns[c])
        for c in kegg_cols:
            e["KEGG"].update(_kegg(r[c]))
        for c in pfam_cols:
            e["Pfam"].update(PFAM_RE.findall(str(r[c])))
    roots = {"GO:0008150", "GO:0003674", "GO:0005575"}
    with_go = [v["GO"] for v in rec.values() if v["GO"]]
    stats.update(ids=len(rec), with_go=len(with_go), go_columns=go_cols, evalue_column=ecol,
                 # part des protéines portant une racine GO : ~1 si déjà propagé
                 root_share=(sum(1 for g in with_go if g & roots) / len(with_go))
                 if with_go else 0.0)
    return rec, names, nss, stats


# ── Gene Ontology (optionnel) : noms, espaces, propagation ────────────────────

def read_obo(path: str):
    """go-basic.obo (ou go-edit.obo) -> (names, namespaces, parents, alt).

    parents{GO: set} = is_a + part_of de même espace (règle de go-basic : pas de
    propagation entre BP, CC et MF) + genre et part_of des définitions logiques
    (intersection_of), qui portent le parent des termes que go-edit.obo ne relie
    plus par is_a (celui-ci est alors inféré par le raisonneur à la publication).
    Les termes obsolètes sont ignorés ; alt_id -> id principal."""
    names, nss, isa, partof, alt = {}, {}, {}, {}, {}
    cur = None

    def flush(c):
        if not c or c.get("obsolete") or not c.get("id"):
            return
        i = c["id"]
        names[i], nss[i] = c.get("name", i), c.get("namespace", "")
        isa[i], partof[i] = c.get("isa", set()), c.get("partof", set())
        for a in c.get("alt", []):
            alt[a] = i

    with open(path, encoding="utf-8", errors="replace") as fh:
        for raw in fh:
            line = raw.strip()
            if line.startswith("[") and line.endswith("]"):
                flush(cur)
                cur = {} if line == "[Term]" else None
                continue
            if cur is None or ":" not in line:
                continue
            k, v = line.split(":", 1)
            v = v.split("!")[0].split("{")[0].strip()
            tok = v.split()
            if k == "id":
                cur["id"] = v
            elif k == "name":
                cur["name"] = raw.split(":", 1)[1].strip()
            elif k == "namespace":
                cur["namespace"] = v
            elif k == "alt_id":
                cur.setdefault("alt", []).append(v)
            elif k == "is_obsolete" and v == "true":
                cur["obsolete"] = True
            elif k == "is_a" and tok:
                cur.setdefault("isa", set()).add(tok[0])
            elif k == "relationship" and len(tok) >= 2 and tok[0] == "part_of":
                cur.setdefault("partof", set()).add(tok[1])
            elif k == "intersection_of" and tok:
                if len(tok) == 1:                       # genre de la définition logique
                    cur.setdefault("isa", set()).add(tok[0])
                elif tok[0] == "part_of":
                    cur.setdefault("partof", set()).add(tok[1])
        flush(cur)
    parents = {}
    for i in names:
        same_ns = {p for p in partof.get(i, ()) if nss.get(p) == nss.get(i)}
        parents[i] = {alt.get(p, p) for p in isa.get(i, set()) | same_ns}
    return names, nss, parents, alt


def propagate(gos: set, parents: dict, alt: dict, cache: dict) -> set:
    """Ajoute tous les ancêtres (is_a, part_of) : règle du « vrai chemin »."""
    out = set()
    for g in gos:
        g = alt.get(g, g)
        if g not in cache:
            seen, stack = set(), [g]
            while stack:
                x = stack.pop()
                if x in seen:
                    continue
                seen.add(x)
                stack.extend(parents.get(x, ()))
            cache[g] = seen
        out |= cache[g]
    return out


# ── Assemblage ────────────────────────────────────────────────────────────────

def build(tables, evalue_max=1e-5, fasta=None, obo=None, do_propagate=True, verbose=True):
    """-> (annotation DataFrame, terms DataFrame, report dict)."""
    merged, names, nss, report = {}, {}, {}, {"tables": []}
    for path in tables:
        rec, nm, ns, st = parse_table(read_table(path), evalue_max=evalue_max)
        st["file"] = os.path.basename(path)
        report["tables"].append(st)
        for k, v in nm.items():
            names.setdefault(k, v)
        for k, v in ns.items():
            nss.setdefault(k, v)
        for pid, e in rec.items():
            key = canonical_id(pid)
            m = merged.setdefault(key, {"raw": pid, "GO": set(), "KEGG": set(), "Pfam": set(),
                                        "src": []})
            for f in ("GO", "KEGG", "Pfam"):
                m[f] |= e[f]
            if e["GO"] or e["KEGG"] or e["Pfam"]:
                m["src"].append(os.path.basename(path))
    shares = [t["root_share"] for t in report["tables"] if t["with_go"]]
    report["depth_warning"] = (not obo and len(shares) > 1
                               and max(shares) > 0.5 and min(shares) < 0.2)
    parents = alt = None
    if obo:
        onames, onss, parents, alt = read_obo(obo)
        names = {**names, **{k: v for k, v in onames.items()}}
        nss = {**nss, **{k: v for k, v in onss.items() if v}}
        report["obo"] = os.path.basename(obo)
    cache = {}
    if obo:
        # identifiants secondaires -> principal ; termes obsolètes / inconnus retirés
        dropped = set()
        for m in merged.values():
            kept = {alt.get(g, g) for g in m["GO"]}
            dropped |= {g for g in kept if g not in parents}
            m["GO"] = {g for g in kept if g in parents}
            if do_propagate and m["GO"]:
                m["GO"] = {g for g in propagate(m["GO"], parents, alt, cache)
                           if g in parents and g.startswith("GO:")}
        report["obsolete_dropped"] = len(dropped)
    # identifiants de sortie : ceux du FASTA de référence si fourni
    if fasta:
        fmap = {}
        for fid in fasta_ids(fasta):
            fmap.setdefault(canonical_id(fid), fid)
        out_id = {k: fmap.get(k) for k in merged}
        report["fasta_ids"] = len(fmap)
        report["matched_to_fasta"] = sum(1 for v in out_id.values() if v)
    else:
        out_id = {k: default_display_id(m["raw"]) for k, m in merged.items()}
    rows = []
    for k, m in merged.items():
        if not out_id.get(k) or not (m["GO"] or m["KEGG"] or m["Pfam"]):
            continue
        rows.append({"protein_id": out_id[k], "GO": ";".join(sorted(m["GO"])),
                     "KEGG": ";".join(sorted(m["KEGG"])), "Pfam": ";".join(sorted(m["Pfam"])),
                     "source": "+".join(dict.fromkeys(m["src"])) or "table"})
    ann = pd.DataFrame(rows, columns=["protein_id", "GO", "KEGG", "Pfam", "source"])
    ann = ann.drop_duplicates("protein_id").sort_values("protein_id").reset_index(drop=True)
    all_go = sorted({g for gos in ann.GO for g in gos.split(";") if g})
    terms = pd.DataFrame({"GO": all_go,
                          "name": [names.get(g, "") for g in all_go],
                          "namespace": [nss.get(g, "") for g in all_go]})
    report.update(proteins=len(ann), with_go=int((ann.GO != "").sum()),
                  go_terms=len(all_go), named=int((terms.name != "").sum()),
                  with_namespace=int((terms.namespace != "").sum()),
                  propagated=bool(obo and do_propagate))
    return ann, terms, report


def print_report(rep):
    for t in rep["tables"]:
        ev = (f", {t['dropped_evalue']} rows dropped (e-value > threshold, column "
              f"'{t['evalue_column']}')" if t.get("evalue_column") else ", no e-value column")
        print(f"  [READ] {t['file']}: {t['rows']} rows, {t['ids']} IDs, {t['with_go']} with GO{ev}")
    if "fasta_ids" in rep:
        print(f"  [IDS]  {rep['matched_to_fasta']} annotated IDs found in the FASTA "
              f"({rep['fasta_ids']} sequences)")
    if rep.get("depth_warning"):
        print("  [WARN] one table lists GO terms with all their ancestors (propagated), another "
              "only the direct terms: mixing them biases the enrichment towards the propagated "
              "proteins. Give --obo go-basic.obo (propagates every protein and names every "
              "term), or import a single table.")
    if rep.get("obsolete_dropped"):
        print(f"  [OBO]  {rep['obsolete_dropped']} obsolete or unknown GO IDs removed "
              f"({rep['obo']})")
    print(f"  [OUT]  {rep['proteins']} proteins, {rep['with_go']} with GO; {rep['go_terms']} GO "
          f"terms ({rep['named']} named, {rep['with_namespace']} with BP/CC/MF)"
          + ("; annotations propagated to ancestors (OBO)" if rep["propagated"] else ""))


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__ or "annotation import",
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--table", action="append", required=True,
                    help="annotation table (xlsx / csv / tsv / topGO background); repeatable, "
                         "annotations of the same protein are merged")
    ap.add_argument("--out", required=True, help="annotation TSV (protein_id, GO, KEGG, Pfam)")
    ap.add_argument("--terms-out", default=None, help="GO names TSV (GO, name, namespace)")
    ap.add_argument("--evalue-max", type=float, default=1e-5,
                    help="drop rows whose e-value column is above (default 1e-5; tables "
                         "without e-value are kept)")
    ap.add_argument("--fasta", default=None,
                    help="FASTA whose headers the IDs must match (the reference of Perseverance)")
    ap.add_argument("--obo", default=None,
                    help="go-basic.obo: names and namespaces of every term + propagation")
    ap.add_argument("--no-propagate", action="store_true",
                    help="with --obo: keep names/namespaces but do not add ancestors")
    a = ap.parse_args(argv)
    ann, terms, rep = build(a.table, a.evalue_max, a.fasta, a.obo, not a.no_propagate)
    ann.to_csv(a.out, sep="\t", index=False)
    print_report(rep)
    print(f"  [OK] {a.out}")
    if a.terms_out:
        terms.to_csv(a.terms_out, sep="\t", index=False)
        print(f"  [OK] {a.terms_out}")
    if rep["proteins"] == 0:
        print("  [WARN] no annotation written (check --fasta IDs / --evalue-max).")
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
