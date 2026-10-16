"""Diagnostic du GO via orthologues (Perseverance -> gProfiler), sans relancer Deimos.

Usage, depuis le dossier Deimos, après un run avec Perseverance :
    python diag_go.py [contraste] [organisme]
    (défauts : embMajor24_vs_embMajor0, athaliana)

Vérifie dans l'ordre :
  1. le format des identifiants d'orthologues (correctif rbh.clean_seqid installé ?) ;
  2. la part des orthologues reconnue par gProfiler pour l'organisme ;
  3. l'enrichissement selon le fond statistique : génome annoté, 'custom'
     (fond = protéines quantifiées), 'custom_annotated' (méthode Deimos).
"""
import os
import sys

import pandas as pd

sys.path.insert(0, os.getcwd())
import rbh                                                         # noqa: E402
from deimos_to_perseverance import build_ortholog_map              # noqa: E402

CONTRAST = sys.argv[1] if len(sys.argv) > 1 else "embMajor24_vs_embMajor0"
ORGANISM = sys.argv[2] if len(sys.argv) > 2 else "athaliana"
OUT = "proteomics_output"
PV = os.path.join(OUT, "perseverance")
XLSX = os.path.join(OUT, "ProteomicAnalysis_Results.xlsx")


def _clean(s):
    if hasattr(rbh, "clean_seqid"):
        return rbh.clean_seqid(s)
    return s


# --- 1. Orthologues (relus depuis les sorties DIAMOND du run) ---
print(f"1. Correctif v6.3 (rbh.clean_seqid) installé : {hasattr(rbh, 'clean_seqid')}")
fwd = pd.read_csv(os.path.join(PV, "fwd.tsv"), sep="\t", header=None, names=rbh.BLAST_COLUMNS)
rev = pd.read_csv(os.path.join(PV, "rev.tsv"), sep="\t", header=None, names=rbh.BLAST_COLUMNS)
fb, rb = rbh._best_hit_per_query(fwd), rbh._best_hit_per_query(rev)
rmap = dict(zip(rb["qseqid"], rb["sseqid"]))
rows = []
for r in fb.itertuples(index=False):
    if rmap.get(r.sseqid) == r.qseqid:
        cov = 100.0 * r.length / min(r.qlen, r.slen)
        if r.pident >= 25 and cov >= 50:
            rows.append({"query_id": _clean(r.qseqid), "reference_id": _clean(r.sseqid),
                         "bitscore": r.bitscore})
omap = build_ortholog_map(pd.DataFrame(rows), os.path.join(PV, "mapping.tsv"))
ex = list(omap.values())[:3]
print(f"   {len(omap)} orthologues ; exemples envoyés à gProfiler : {ex}")
if any("|" in str(x) for x in ex):
    print("   -> PROBLÈME : identifiants au format sp|ACC|NOM, non reconnus par gProfiler. "
          "Installer le correctif v6.3 (rbh.py).")

# --- 2. Requête et fond du contraste ---
df = pd.read_excel(XLSX, sheet_name="Differential_Expression")
names = df["name"].astype(str).str.split(";").str[0].str.strip()
sig = (df[f"{CONTRAST}_p.adj"] < 0.05) & (df[f"{CONTRAST}_diff"].abs() >= 1)
query = [omap[p] for p in names[sig] if p in omap]
bg = list(dict.fromkeys([omap[p] for p in names if p in omap] + query))
print(f"\n2. {CONTRAST} : {int(sig.sum())} protéines significatives, "
      f"{len(query)} avec orthologue ; fond = {len(bg)} orthologues")

from gprofiler import GProfiler                                    # noqa: E402
gp = GProfiler(return_dataframe=True)
conv = gp.convert(organism=ORGANISM, query=bg, target_namespace="ENSG")
ok = conv[conv["converted"].astype(str) != "None"]["incoming"].nunique()
print(f"   gProfiler '{ORGANISM}' reconnaît {ok}/{len(bg)} identifiants du fond "
      f"({100 * ok / max(len(bg), 1):.0f} %)")

# --- 3. Enrichissement selon le fond ---
import go_enrichment as ge                                         # noqa: E402
kw = dict(organism=ORGANISM, query=query, sources=ge.GO_SOURCES, user_threshold=0.05,
          no_evidences=True, no_iea=False, significance_threshold_method="fdr")
print("\n3. Termes significatifs (FDR < 0.05) selon le fond :")
tests = [("génome annoté (domain_scope=annotated)", lambda: gp.profile(**kw)),
         ("fond quantifié, 'custom' (client gProfiler)",
          lambda: gp.profile(background=bg, **kw)),
         ("fond quantifié, 'custom_annotated' (Deimos)",
          lambda: ge.run_gost(query, ORGANISM, background=bg))]
for label, fn in tests:
    try:
        res = fn()
        n = 0 if res is None else len(res)
        top = "" if not n else " | " + "; ".join(res.sort_values("p_value")["name"].head(4))
        print(f"   - {label} : {n}{top}")
    except Exception as e:                                        # noqa: BLE001
        print(f"   - {label} : ERREUR {type(e).__name__}: {str(e)[:150]}")
