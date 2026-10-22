"""
id_remap.py — Conversion des identifiants protéiques des matrices DIA-NN
(report.pg_matrix.tsv + report.pr_matrix.tsv) par une table de correspondance.

Cas d'usage : recherche DIA-NN faite sur un FASTA dont les identifiants ne sont
reconnus ni par STRING ni par UniProt (ex. UniParc 'UPI000000305D') -> on remplace
le Protein.Group par un identifiant exploitable (locus tag, accession UniProtKB),
sans relancer DIA-NN.

Principes
  - Seules les colonnes d'identifiants/annotation sont réécrites (Protein.Group,
    Protein.Ids, Protein.Names, Genes, First.Protein.Description) ; les intensités
    et les fins de ligne restent identiques à l'octet près.
  - Correspondance EXACTE (espaces retirés), membre par membre des groupes 'A;B'.
    Un identifiant absent de la table (contaminant cRAP, entrée non couverte) est
    laissé tel quel.
  - pg et pr_matrix sont convertis avec la MÊME table résolue -> les Protein.Group
    restent appariables entre les deux fichiers (DPC, DEqMS).
  - Collision (deux Protein.Group différents -> même nouvel identifiant) : les
    groupes concernés gardent leur identifiant d'origine (jamais de fusion).
  - Idempotent : sur des fichiers déjà convertis, ne fait rien.

Usage Deimos (YAML) — source = le FASTA donné à DIA-NN (recommandé, exact) :
    id_mapping:
      path: JIM8777_UniParc.fasta   # en-têtes '>UPI… desc OS=… OX=… GN=… AC=…'
      tag: SSAL            # optionnel -> Protein.Names = <ancien ID>_<tag>
    Par défaut pour un FASTA : from = id (1er mot de l'en-tête), to = AC,
    gene = GN, description = texte avant le 1er champ KEY=. Tout champ KEY= de
    l'en-tête est utilisable (ex. to: GN). En-tête UniProt 'sp|ACC|NAME' :
    id = ACC, champ 'entry' = NAME.

Usage Deimos (YAML) — source = table (.tsv/.csv/.txt/.xlsx) :
    id_mapping:
      path: JIM8777_UPI_locus_mapping.tsv
      from: upi            # colonne = identifiants actuels du Protein.Group
      to: locus_tag        # colonne = nouvel identifiant
      gene: gene           # optionnel -> colonne Genes
      description: description   # optionnel -> First.Protein.Description
      tag: SSAL            # optionnel -> Protein.Names = <ancien ID>_<tag>
      sheet: 0             # optionnel (Excel)

Usage autonome :
    python id_remap.py --map table.tsv --from upi --to locus_tag \
        --gene gene --description description --tag SSAL \
        --pg report.pg_matrix.tsv --pr report.pr_matrix.tsv --out-dir converted/
"""
from __future__ import annotations

import os
import re
import pandas as pd

ID_COLS = ("Protein.Group", "Protein.Ids")
_ALL_COLS = ID_COLS + ("Protein.Names", "Genes", "First.Protein.Description")
FASTA_EXT = (".fasta", ".fa", ".faa", ".fas")
FASTA_DEFAULTS = {"from": "id", "to": "AC", "gene": "GN", "description": "description"}
_KV = re.compile(r"\s([A-Z]{2})=(.*?)(?=\s[A-Z]{2}=|$)")


# ── Table de correspondance ───────────────────────────────────────────────────

def _fasta_headers(path: str) -> pd.DataFrame:
    """Une ligne par en-tête : id, entry (UniProt), description, champs KEY=VALUE
    (OS, OX, GN, AC, SS, PC...). Formats UniParc ('UPI… desc AC=…') et UniProt
    ('sp|ACC|NAME desc OS=…')."""
    rows = []
    with open(path, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if not line.startswith(">"):
                continue
            h = line[1:].strip()
            first, _, rest = h.partition(" ")
            rec = {"id": first, "entry": None}
            if first.startswith(("sp|", "tr|")) and first.count("|") >= 2:
                _, rec["id"], rec["entry"] = first.split("|", 2)
            m = _KV.search(" " + rest)
            rec["description"] = (rest[:m.start()] if m else rest).strip() or None
            rec.update({k: v.strip() for k, v in _KV.findall(" " + rest)})
            rows.append(rec)
    if not rows:
        raise ValueError(f"id_mapping : aucun en-tête FASTA dans {path}")
    return pd.DataFrame(rows).astype(object)


def _read_table(path: str, id_from: str, sheet=0) -> pd.DataFrame:
    ext = os.path.splitext(path)[1].lower()
    if ext in FASTA_EXT:
        return _fasta_headers(path)
    if ext in (".xlsx", ".xlsm", ".xls"):
        raw = pd.read_excel(path, sheet_name=sheet, header=None, dtype=str)
        hit = raw.apply(lambda r: r.astype(str).str.strip().eq(id_from).any(), axis=1)
        if not hit.any():
            raise ValueError(f"id_mapping : colonne '{id_from}' introuvable dans {path}")
        h = int(hit.idxmax())                     # 1re ligne contenant le nom de colonne
        df = raw.iloc[h + 1:].copy()
        df.columns = raw.iloc[h].astype(str).str.strip()
        return df
    sep = "\t" if ext in (".tsv", ".txt", ".tab") else None
    return pd.read_csv(path, sep=sep, engine="python", dtype=str)


def load_id_mapping(path: str, id_from: str, id_to: str, gene: str | None = None,
                    description: str | None = None, sheet=0) -> pd.DataFrame:
    """Table indexée par l'identifiant d'origine : colonnes to / gene / description.
    Cellules nettoyées (espaces). Lignes sans from/to ignorées ; un 'from' associé
    à plusieurs 'to' différents est exclu (ambigu)."""
    df = _read_table(path, id_from, sheet)
    df.columns = [str(c).strip() for c in df.columns]
    need = [c for c in (id_from, id_to, gene, description) if c]
    miss = [c for c in need if c not in df.columns]
    if miss:
        raise ValueError(f"id_mapping : colonne(s) {miss} absente(s) de {path} "
                         f"(colonnes : {list(df.columns)})")
    df = df[need].apply(lambda s: s.str.strip()).replace({"": None, "nan": None})
    df = df.dropna(subset=[id_from, id_to])
    amb = df.groupby(id_from)[id_to].nunique()
    amb = amb[amb > 1].index
    if len(amb):
        print(f"  [WARN] id_mapping : {len(amb)} identifiant(s) associé(s) à plusieurs "
              f"cibles, exclus : {list(amb[:5])}")
    df = df[~df[id_from].isin(amb)].drop_duplicates(id_from)
    out = pd.DataFrame({"to": df[id_to].values}, index=df[id_from].values)
    out["gene"] = df[gene].values if gene else None
    out["description"] = df[description].values if description else None
    return out.astype(object).where(out.notna(), None)      # NaN -> None


# ── Conversion ────────────────────────────────────────────────────────────────

def _split_eol(line: str) -> tuple[str, str]:
    if line.endswith("\r\n"):
        return line[:-2], "\r\n"
    if line.endswith("\n"):
        return line[:-1], "\n"
    return line, ""


def _read_lines(path: str) -> tuple[list[str], dict]:
    with open(path, "r", encoding="utf-8", newline="") as f:
        lines = f.readlines()
    head, _ = _split_eol(lines[0])
    cols = head.split("\t")
    return lines, {c: cols.index(c) for c in _ALL_COLS if c in cols}


def _resolve_groups(groups: set[str], mp: pd.DataFrame) -> tuple[dict, dict]:
    """{Protein.Group d'origine -> nouveau}, et statut par groupe. Les collisions
    (plusieurs groupes -> même nouvel ID) sont annulées."""
    to = mp["to"].to_dict()
    new, status = {}, {}
    for g in groups:
        mem = [m.strip() for m in g.split(";")]
        hit = [m in to for m in mem]
        if not any(hit):
            status[g] = "unmapped"
            continue
        new[g] = ";".join(to.get(m, m) for m in mem)
        status[g] = "mapped" if all(hit) else "partial"
    inv = {}
    for g, n in new.items():
        inv.setdefault(n, []).append(g)
    for n, gs in inv.items():
        if len(gs) > 1:
            for g in gs:
                new.pop(g)
                status[g] = f"collision -> {n} (kept original)"
    return new, status


def _rewrite(lines: list[str], idx: dict, new: dict, mp: pd.DataFrame,
             tag: str | None) -> list[str]:
    to, gene, desc = (mp["to"].to_dict(), mp["gene"].to_dict(),
                      mp["description"].to_dict())
    has_gene, has_desc = mp["gene"].notna().any(), mp["description"].notna().any()
    ipg = idx["Protein.Group"]
    out = [lines[0]]
    for line in lines[1:]:
        body, eol = _split_eol(line)
        f = body.split("\t")
        pg = f[ipg]
        if pg not in new:
            out.append(line)
            continue
        mem = [m.strip() for m in pg.split(";")]
        f[ipg] = new[pg]
        if "Protein.Ids" in idx:
            ids = [m.strip() for m in f[idx["Protein.Ids"]].split(";")]
            f[idx["Protein.Ids"]] = ";".join(to.get(m, m) for m in ids)
        if tag and "Protein.Names" in idx:
            f[idx["Protein.Names"]] = ";".join(f"{m}_{tag}" for m in mem)
        if has_gene and "Genes" in idx:
            gs = [gene.get(m) for m in mem if m in to]
            f[idx["Genes"]] = ";".join(dict.fromkeys(g for g in gs if isinstance(g, str) and g))
        if has_desc and "First.Protein.Description" in idx:
            d = next((desc.get(m) for m in mem if desc.get(m)), None)
            if d:
                f[idx["First.Protein.Description"]] = d
        out.append("\t".join(f) + eol)
    return out


def remap_diann_files(files: dict, mp: pd.DataFrame, out_dir: str,
                      tag: str | None = None) -> tuple[dict, pd.DataFrame]:
    """files : {'pg': chemin, 'pr': chemin ou None}. Retourne ({clé: nouveau
    chemin}, rapport par Protein.Group). Fichier déjà converti -> renvoyé tel quel."""
    loaded = {}
    for k, p in files.items():
        if p and os.path.exists(str(p)):
            lines, idx = _read_lines(p)
            if "Protein.Group" not in idx:
                raise ValueError(f"id_mapping : pas de colonne Protein.Group dans {p}")
            loaded[k] = (p, lines, idx)
    groups = set()
    for _, lines, idx in loaded.values():
        i = idx["Protein.Group"]
        groups.update(_split_eol(l)[0].split("\t")[i] for l in lines[1:])
    members = {m.strip() for g in groups for m in g.split(";")}
    n_from = len(members & set(mp.index))
    n_to = len(members & set(mp["to"]))
    if n_from == 0:
        if n_to >= 0.5 * len(members):
            print("  [ID-MAP] Identifiants déjà convertis — fichiers utilisés tels quels.")
        else:
            print(f"  [WARN] id_mapping : aucun Protein.Group trouvé dans la colonne "
                  f"'from' de la table — fichiers utilisés tels quels.")
        return {k: v[0] for k, v in loaded.items()}, pd.DataFrame()

    new, status = _resolve_groups(groups, mp)
    os.makedirs(out_dir, exist_ok=True)
    paths = {}
    for k, (p, lines, idx) in loaded.items():
        dst = os.path.join(out_dir, os.path.basename(p))
        if os.path.abspath(dst) == os.path.abspath(p):
            raise ValueError("id_mapping : le dossier de sortie est celui des fichiers d'origine.")
        with open(dst, "w", encoding="utf-8", newline="") as f:
            f.writelines(_rewrite(lines, idx, new, mp, tag))
        paths[k] = dst
    rep = (pd.DataFrame({"original": sorted(groups)})
             .assign(new=lambda d: d["original"].map(new).fillna(d["original"]),
                     status=lambda d: d["original"].map(status)))
    return paths, rep


# ── Point d'entrée Deimos ─────────────────────────────────────────────────────

def apply_id_mapping(params: dict, tsv_path: str, pr_path: str | None,
                     out_dir: str) -> tuple[str, str | None]:
    """Hook appelé par deimos.main() quand params['id_mapping'] est renseigné.
    Écrit les matrices converties dans <out_dir>/id_mapping/ et renvoie leurs
    chemins. params['tsv_path'] n'est PAS modifié (labels des runs, dossier des
    FASTA, last_config.yaml restent ceux d'origine)."""
    cfg = params.get("id_mapping") or {}
    if isinstance(cfg, str):
        raise ValueError("id_mapping doit être un bloc {path, from, to, ...}, pas une chaîne.")
    path = cfg.get("path")
    if path and not os.path.exists(path):
        alt = os.path.join(os.path.dirname(os.path.abspath(tsv_path)), path)
        path = alt if os.path.exists(alt) else path
    if not path or not os.path.exists(path):
        raise FileNotFoundError(f"id_mapping.path introuvable : {cfg.get('path')}")
    if path.lower().endswith(FASTA_EXT):          # FASTA : colonnes par défaut
        cfg = {**FASTA_DEFAULTS, **{k: v for k, v in cfg.items() if v is not None}}
    for k in ("from", "to"):
        if not cfg.get(k):
            raise ValueError(f"id_mapping.{k} manquant (obligatoire pour une table).")
    print(f"\n[STEP] Conversion des identifiants ({os.path.basename(path)} : "
          f"'{cfg.get('from')}' -> '{cfg.get('to')}')...")
    mp = load_id_mapping(path, cfg["from"], cfg["to"], cfg.get("gene"),
                         cfg.get("description"), cfg.get("sheet", 0))
    work = os.path.join(out_dir, "id_mapping")
    paths, rep = remap_diann_files({"pg": tsv_path, "pr": pr_path}, mp, work,
                                   tag=cfg.get("tag"))
    if rep.empty:
        params["id_mapping_summary"] = "not applied (identifiers already converted or not found)"
        return paths.get("pg", tsv_path), paths.get("pr", pr_path)
    rep.to_csv(os.path.join(work, "id_mapping_report.tsv"), sep="\t", index=False)
    st = rep["status"].str.split(" ").str[0].value_counts().to_dict()
    n = len(rep)
    print(f"  [OK] {st.get('mapped', 0)}/{n} Protein.Group convertis"
          + (f", {st['partial']} partiellement" if st.get("partial") else "")
          + f", {st.get('unmapped', 0)} inchangés (contaminants / absents de la table)"
          + (f", {st['collision']} collision(s) laissée(s) d'origine" if st.get("collision") else ""))
    unm = rep.loc[rep["status"] == "unmapped", "original"]
    unm = unm[~unm.str.lower().str.startswith(("crap", "con__", "rev_"))]
    if len(unm):
        print(f"  [WARN] {len(unm)} protéine(s) hors contaminants non converties : "
              f"{unm.head(5).tolist()}{'...' if len(unm) > 5 else ''}")
    print(f"  -> {work}  (rapport : id_mapping_report.tsv)")
    params["id_mapping_summary"] = (
        f"{os.path.basename(path)}: '{cfg['from']}' -> '{cfg['to']}'; "
        f"{st.get('mapped', 0)}/{n} protein groups converted, "
        f"{st.get('unmapped', 0)} unchanged, {st.get('collision', 0)} collisions kept original")
    return paths.get("pg", tsv_path), paths.get("pr", pr_path)


if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--map", required=True)
    ap.add_argument("--from", dest="id_from")   # défaut FASTA : id
    ap.add_argument("--to", dest="id_to")       # défaut FASTA : AC
    ap.add_argument("--gene")
    ap.add_argument("--description")
    ap.add_argument("--tag")
    ap.add_argument("--sheet", default=0)
    ap.add_argument("--pg", required=True)
    ap.add_argument("--pr")
    ap.add_argument("--out-dir", required=True)
    a = ap.parse_args()
    params = {"id_mapping": {"path": a.map, "from": a.id_from, "to": a.id_to,
                             "gene": a.gene, "description": a.description,
                             "tag": a.tag, "sheet": a.sheet}}
    pg, pr = apply_id_mapping(params, a.pg, a.pr, a.out_dir)
    # apply_id_mapping écrit dans <out_dir>/id_mapping/ (convention Deimos)
    print(f"pg -> {pg}\npr -> {pr}")
