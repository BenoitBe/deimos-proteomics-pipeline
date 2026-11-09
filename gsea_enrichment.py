# ==============================================================================
# gsea_enrichment.py — GSEA rank-based (pré-ranked) optionnel pour Deimos
# Version Python — Proteogen
# ==============================================================================
# Complément de l'ORA (go_enrichment.py). Là où l'ORA part de la LISTE des DEP
# (coupée au seuil), le GSEA classe TOUTES les protéines par un score signé et
# détecte les gene sets concentrés en haut (up) ou en bas (down) du classement —
# captant les effets coordonnés que l'ORA rate (beaucoup de protéines qui
# bougent un peu dans le même sens, sans passer le seuil).
#
# Métrique de classement : -log10(p) * sign(LFC)
#   -> sépare mieux les extrêmes que le LFC seul (combine sens ET significativité)
#
# Moteur : gseapy.prerank (algorithme Broad, permutations). L'algorithme tourne
# en local ; seul le téléchargement des gene sets (Enrichr) nécessite le réseau.
#
# ENTIÈREMENT NON BLOQUANT : toute erreur (gseapy absent, réseau, mapping vide,
# pas assez de protéines) est capturée et n'interrompt jamais le pipeline.
# ==============================================================================

import os
import re
import time
import warnings
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

warnings.filterwarnings("ignore")

# Gene sets Enrichr par défaut (couvrent GO + pathways, cohérent avec l'ORA).
# Jeux humains : défaut de gsea_organism: human (comportement d'origine).
DEFAULT_GENE_SETS = [
    "GO_Biological_Process_2021",
    "GO_Molecular_Function_2021",
    "GO_Cellular_Component_2021",
    "KEGG_2021_Human",
    "Reactome_2022",
]

# gsea_organism: mouse — collections MSigDB souris (symboles murins). Pas de KEGG
# dans MSigDB souris : Reactome + WikiPathways pour les voies.
DEFAULT_GENE_SETS_MOUSE = [
    "msigdb:mh.all",
    "msigdb:m2.cp.reactome",
    "msigdb:m2.cp.wikipathways",
    "msigdb:m5.go.bp",
    "msigdb:m5.go.mf",
    "msigdb:m5.go.cc",
]
DEFAULT_MSIGDB_VERSION = {"human": "2025.1.Hs", "mouse": "2025.1.Mm"}
GSEA_ORGANISMS = ("human", "mouse")

# gsea_organism: athaliana (ou gprofiler:<code> pour toute espèce de g:Profiler,
# ex. gprofiler:hvulgare, gprofiler:osativa) — jeux de gènes de g:Profiler (même
# source que l'ORA : GO, et les voies que g:Profiler distribue pour l'espèce),
# fichier GMT téléchargé une fois puis gardé en cache. Les protéines sont reliées
# aux gènes par g:Convert (accession UniProt -> gène Ensembl, AGI chez Arabidopsis),
# directement ou via les orthologues Perseverance vers la même espèce (B. napus ->
# A. thaliana). Une espèce sans symboles humains/murins n'a plus besoin de la
# colonne Genes.
GPROFILER_GMT_URL = "https://biit.cs.ut.ee/gprofiler/static/gprofiler_full_{org}.name.gmt"
GPROFILER_NAMED = {"athaliana": "athaliana", "arabidopsis": "athaliana",
                   "a.thaliana": "athaliana", "arabidopsis_thaliana": "athaliana"}
GMT_CACHE_DAYS = 30            # GMT g:Profiler retéléchargé au-delà (ou supprimer le fichier)
_GP_CODE = re.compile(r"^[a-z][a-z0-9_]+$")
_NO_HIT = {"", "NONE", "N/A", "NAN"}
LAST_SOURCES = []              # provenance des jeux chargés (écrite dans les onglets GSEA)


def gprofiler_code(organism) -> "str | None":
    """Code g:Profiler de gsea_organism ('athaliana', 'gprofiler:hvulgare'), None pour
    human / mouse (Enrichr / MSigDB) ou une valeur invalide."""
    o = str(organism or "").strip().lower()
    if o in GPROFILER_NAMED:
        return GPROFILER_NAMED[o]
    if o.startswith("gprofiler:"):
        code = o.split(":", 1)[1].strip()
        return code if _GP_CODE.match(code) else None
    return None


def valid_gsea_organism(organism) -> bool:
    return str(organism or "").strip().lower() in GSEA_ORGANISMS or gprofiler_code(organism) is not None


def gsea_settings(params, go_params=None) -> dict:
    """Réglages GSEA lus dans la configuration (clés gsea_*), avec repli sur
    go_params (ancien passage par go_params['gsea_gene_sets'])."""
    gp = go_params or {}
    org = str(params.get("gsea_organism") or gp.get("gsea_organism") or "human").lower()
    code = gprofiler_code(org)
    if code:
        org = code
    elif org not in GSEA_ORGANISMS:
        org = "human"
    sets = params.get("gsea_gene_sets") or gp.get("gsea_gene_sets")
    if isinstance(sets, str):
        sets = [sets]
    if not sets:
        sets = ([f"gprofiler:{code}"] if code else
                list(DEFAULT_GENE_SETS_MOUSE if org == "mouse" else DEFAULT_GENE_SETS))
    return {
        "organism": org,
        "gprofiler": code,
        "gene_sets": list(sets),
        "go_obo": params.get("gsea_go_obo") or None,
        "cache_dir": params.get("gsea_cache_dir") or None,
        "msigdb_version": str(params.get("gsea_msigdb_version")
                              or DEFAULT_MSIGDB_VERSION.get(org, DEFAULT_MSIGDB_VERSION["human"])),
        "fdr": float(params.get("gsea_fdr", gp.get("gsea_fdr", 0.25))),
        "permutations": int(params.get("gsea_permutations",
                                       gp.get("gsea_permutations", 1000))),
        "rank_p": str(params.get("gsea_rank_p", "pvalue")).lower(),
    }


def read_gmt(path: str) -> dict:
    """Fichier GMT (nom <tab> description <tab> gènes...) -> {nom: [gènes]}."""
    out = {}
    with open(path, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            row = line.rstrip("\n\r").split("\t")
            if len(row) >= 3 and row[0]:
                out[row[0]] = [g for g in row[2:] if g]
    return out


def default_cache_dir() -> str:
    return os.path.join(os.path.expanduser("~"), ".deimos_cache", "gsea")


def fetch_gprofiler_gmt(code: str, cache_dir: "str | None" = None, verbose: bool = True,
                        max_age_days: float = GMT_CACHE_DAYS) -> str:
    """Chemin local du GMT complet de g:Profiler pour l'espèce `code` (gène = nom
    Ensembl, colonne 1 = identifiant du terme, colonne 2 = son nom). Téléchargé dans
    cache_dir (~/.deimos_cache/gsea) et réutilisé tant qu'il a moins de max_age_days.
    Lève une erreur claire si l'espèce est inconnue de g:Profiler ou hors ligne."""
    import requests
    cache_dir = cache_dir or default_cache_dir()
    os.makedirs(cache_dir, exist_ok=True)
    path = os.path.join(cache_dir, f"gprofiler_full_{code}.name.gmt")
    if os.path.isfile(path) and (time.time() - os.path.getmtime(path)) < max_age_days * 86400:
        return path
    url = GPROFILER_GMT_URL.format(org=code)
    if verbose:
        print(f"     g:Profiler gene sets for '{code}': downloading {url}")
    r = requests.get(url, timeout=300, headers={"User-Agent": "Deimos-proteomics-pipeline"})
    if r.status_code != 200:
        if os.path.isfile(path):                    # cache périmé mais utilisable
            print(f"     [WARN] download failed (HTTP {r.status_code}) — cached copy kept.")
            return path
        raise ValueError(f"HTTP {r.status_code} for {url} (organism code unknown to "
                         f"g:Profiler, or offline)")
    first = r.text.split("\n", 1)[0]
    if first.count("\t") < 2:
        raise ValueError(f"{url} is not a GMT file")
    tmp = path + ".part"
    with open(tmp, "w", encoding="utf-8") as fh:
        fh.write(r.text)
    os.replace(tmp, path)
    return path


def read_gprofiler_gmt(path: str, go_namespaces: "dict | None" = None) -> dict:
    """GMT g:Profiler -> {(source, 'ID nom'): [gènes]}. Source = préfixe de l'ID
    (GO, KEGG, REAC, WP…) ; GO séparé en GO:BP / GO:MF / GO:CC si go_namespaces
    ({GO:id: biological_process|molecular_function|cellular_component}) est fourni."""
    ns_short = {"biological_process": "GO:BP", "molecular_function": "GO:MF",
                "cellular_component": "GO:CC"}
    out = {}
    with open(path, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            row = line.rstrip("\n\r").split("\t")
            if len(row) < 3 or not row[0]:
                continue
            tid, name = row[0].strip(), row[1].strip()
            src = tid.split(":", 1)[0] if ":" in tid else "gprofiler"
            if src == "GO" and go_namespaces:
                src = ns_short.get(go_namespaces.get(tid, ""), "GO")
            label = f"{tid} {name}" if name and name != tid else tid
            out[(src, label)] = [g for g in row[2:] if g]
    return out


def load_gene_sets(specs, organism: str = "human",
                   msigdb_version: "str | None" = None, verbose: bool = True,
                   cache_dir: "str | None" = None, go_obo: "str | None" = None) -> dict:
    """Charge les jeux de gènes UNE fois pour tous les contrastes.

    specs : liste de
      - chemins de fichiers .gmt locaux (version figée, sans réseau) ;
      - "msigdb:<collection>" (ex. msigdb:m5.go.bp, msigdb:mh.all), version
        msigdb_version (ex. 2025.1.Mm) ;
      - noms de bibliothèques Enrichr (ex. GO_Biological_Process_2021).
    Retourne {"<source>__<terme>": [SYMBOLES EN MAJUSCULES]} : même convention de
    nommage que gseapy avec plusieurs bibliothèques ; symboles en majuscules comme
    le classement (_ranking_metric), donc Actb (souris) = ACTB (humain).
    Une source introuvable est signalée et ignorée.
      - "gprofiler:<code>" (ex. gprofiler:athaliana) : GMT complet de g:Profiler,
        une entrée par source (GO, ou GO:BP/MF/CC si go_obo est donné ; KEGG, REAC,
        WP… selon l'espèce), terme nommé « ID nom »."""
    import gseapy
    version = msigdb_version or DEFAULT_MSIGDB_VERSION.get(organism, "2025.1.Hs")
    sets = {}
    del LAST_SOURCES[:]
    for spec in specs:
        spec = str(spec).strip()
        try:
            if spec.lower().startswith("gprofiler:"):
                code = spec.split(":", 1)[1].strip().lower()
                gmt = fetch_gprofiler_gmt(code, cache_dir, verbose)
                nss = None
                if go_obo and os.path.isfile(go_obo):
                    from annotation_import import read_obo
                    nss = read_obo(go_obo)[1]
                lib = read_gprofiler_gmt(gmt, nss)
                n0 = len(sets)
                by_src = {}
                for (src, term), genes in lib.items():
                    g = list(dict.fromkeys(str(x).strip().upper() for x in genes if str(x).strip()))
                    if g:
                        sets[f"{src}__{term}"] = g
                        by_src[src] = by_src.get(src, 0) + 1
                LAST_SOURCES.append(
                    f"g:Profiler {code} GMT ({time.strftime('%Y-%m-%d', time.localtime(os.path.getmtime(gmt)))}; "
                    + ", ".join(f"{k} {v}" for k, v in sorted(by_src.items())) + ")")
                if verbose:
                    age = (time.time() - os.path.getmtime(gmt)) / 86400
                    print(f"     g:Profiler {code} ({os.path.basename(gmt)}, "
                          f"{'downloaded today' if age < 1 else f'cached {age:.0f} d ago'}): "
                          + ", ".join(f"{k} {v}" for k, v in sorted(by_src.items()))
                          + f" — {len(sets) - n0} gene sets")
                continue
            if os.path.isfile(spec):
                src = os.path.splitext(os.path.basename(spec))[0]
                lib = read_gmt(spec)
            elif spec.lower().startswith("msigdb:"):
                cat = spec.split(":", 1)[1]
                src = cat
                lib = gseapy.Msigdb.get_gmt(category=cat, dbver=version)
                if not lib:
                    raise ValueError(f"MSigDB {cat} v{version} not found")
            else:
                src = spec
                lib = gseapy.get_library(name=spec, organism=organism.capitalize())
        except Exception as e:                                        # noqa: BLE001
            if verbose:
                print(f"  [WARN] GSEA gene sets '{spec}' unavailable "
                      f"({type(e).__name__}: {str(e)[:80]}) — skipped. Offline: give "
                      f"a local .gmt file in gsea_gene_sets.")
            continue
        n0 = len(sets)
        for term, genes in lib.items():
            g = list(dict.fromkeys(str(x).strip().upper() for x in genes if str(x).strip()))
            if g:
                sets[f"{src}__{term}"] = g
        LAST_SOURCES.append(f"{src}" + (f" (MSigDB {version})" if spec.lower().startswith("msigdb:") else ""))
        if verbose:
            print(f"     {src}: {len(sets) - n0} gene sets")
    return sets


def _leaders(df):
    """Accession de tête (convention Deimos) de chaque ligne."""
    col = "name" if "name" in df.columns else ("Protein.Group" if "Protein.Group" in df.columns
                                               else df.columns[0])
    return df[col].astype(str).str.split(";").str[0].str.strip()


def _key_proteins(df, key_map=None):
    """Clé de classement (majuscules) -> accessions de tête des protéines qui la portent :
    inverse de la clé utilisée par _ranking_metric (g:Profiler ou colonne Genes). Sert à
    rendre le leading edge en protéines (colonne Lead_proteins, dashboard)."""
    pids = _leaders(df)
    if key_map is not None:
        keys = pids.map(key_map)
    elif "Genes" in df.columns:
        keys = df["Genes"]
    else:
        return {}
    keys = keys.astype(str).str.split(";").str[0].str.strip().str.upper()
    out = {}
    for k, p in zip(keys, pids):
        if k and k not in ("NAN", "NONE"):
            out.setdefault(k, []).append(p)
    return out


def _lead_proteins(lead, key2prot):
    """'G1;G2' (Lead_genes de gseapy) -> 'P1;P2;P3'."""
    out = []
    for g in str(lead or "").split(";"):
        out += key2prot.get(g.strip().upper(), [])
    return ";".join(dict.fromkeys(out))


def _convert_ids(ids, code, chunk=2500):
    """g:Convert (g:Profiler) : identifiants -> [(incoming, gène Ensembl, nom)]."""
    from gprofiler import GProfiler
    gp = GProfiler(user_agent="Deimos-proteomics-pipeline", return_dataframe=False)
    out = []
    ids = list(dict.fromkeys(i for i in ids if i))
    for k in range(0, len(ids), chunk):
        res = gp.convert(organism=code, query=ids[k:k + chunk], target_namespace="ENSG")
        out += [(str(r.get("incoming", "")), str(r.get("converted", "")), str(r.get("name", "")))
                for r in res]
    return out


def gprofiler_gene_keys(df, code, orthologs=None, convert_fn=None, verbose=True):
    """Clé de classement par protéine pour les jeux g:Profiler de l'espèce `code`.

    Requête g:Convert = accession de tête, ou son orthologue Perseverance si
    `orthologs` ({protein_id: accession de l'espèce de référence}) est donné.
    Clé = gène Ensembl (AGI chez Arabidopsis). Retourne (key_map {protein_id: clé},
    alias {NOM / AGI / ACCESSION en majuscules: clé}, résumé). Hors ligne ou si
    g:Convert ne reconnaît presque rien : repli sur la colonne Genes (symboles)."""
    convert_fn = convert_fn or _convert_ids
    pids = _leaders(df)
    query = {p: (orthologs.get(p) if orthologs else p) for p in pids}
    query = {p: q for p, q in query.items() if q}
    try:
        rows = convert_fn(list(query.values()), code)
    except Exception as e:                                            # noqa: BLE001
        rows, err = [], f"{type(e).__name__}: {str(e)[:80]}"
    else:
        err = None
    hit = {}
    for inc, conv, name in rows:
        if conv.upper() in _NO_HIT or inc in hit:
            continue
        hit[inc] = (conv.upper(), name.upper())
    key_map, alias = {}, {}
    for p, q in query.items():
        if q in hit:
            key, name = hit[q]
            key_map[p] = key
            alias[key] = key
            alias[q.upper()] = key
            if name not in _NO_HIT:
                alias.setdefault(name, key)
    n_q = len(query)
    if len(key_map) >= max(15, 0.1 * max(n_q, 1)):
        summ = (f"{len(key_map)}/{len(pids)} proteins linked to {code} genes by g:Convert"
                + (f" (through {len(query)} Perseverance orthologues)" if orthologs else ""))
        return key_map, alias, summ
    why = err or f"g:Convert linked only {len(key_map)}/{n_q} identifiers to {code} genes"
    if "Genes" in df.columns:
        g = df["Genes"].astype(str).str.split(";").str[0].str.strip().str.upper()
        key_map = {p: v for p, v in zip(pids, g) if v and v not in _NO_HIT}
        alias = {v: v for v in key_map.values()}
        summ = (f"[WARN] {why} — fallback on the Genes column ({len(key_map)} symbols)"
                + ("" if orthologs else "; for another species, run Perseverance towards "
                   f"'{code}' (perseverance_reference_organism: {code})"))
        return key_map, alias, summ
    return {}, {}, f"[WARN] {why} and no Genes column"


def translate_gene_sets(gene_sets: dict, alias: dict) -> dict:
    """Membres des jeux (noms, AGI, accessions) -> clés de classement ; les membres
    absents des données sont retirés (la GSEA ne tient compte que des gènes classés)."""
    out = {}
    for term, genes in gene_sets.items():
        g = list(dict.fromkeys(alias[x] for x in genes if x in alias))
        if g:
            out[term] = g
    return out


def _ranking_metric(df_comparaison, contrast, params, key_map=None):
    """Construit le vecteur de ranking signé pour un contraste.

    Métrique = -log10(p) * sign(LFC), sur TOUTES les protéines (pas seulement
    les DEP). Mapping accession -> symbole de gène via la colonne 'Genes'.
    p : p brute par défaut (gsea_rank_p: pvalue). La p.adj BH forme des paliers
    (beaucoup de protéines à la même valeur) qui écrasent l'ordre du milieu du
    classement, là où la GSEA détecte les effets coordonnés ; gsea_rank_p: volcano
    reprend l'ancien comportement (p.adj si volcano_use_padj).

    Retourne une pd.Series {gene_symbol: score} triée décroissante, ou None.
    """
    if str(params.get("gsea_rank_p", "pvalue")).lower() == "volcano":
        suffixe_p = "_p.adj" if params.get("volcano_use_padj") else "_p.val"
    else:
        suffixe_p = "_p.val"
    col_p    = f"{contrast}{suffixe_p}"
    col_diff = f"{contrast}_diff"
    if col_p not in df_comparaison.columns or col_diff not in df_comparaison.columns:
        return None

    if key_map is not None:                       # clés g:Profiler (gène Ensembl / AGI)
        df = pd.DataFrame({"gene": _leaders(df_comparaison).map(key_map).values,
                           "p": df_comparaison[col_p].values,
                           "lfc": df_comparaison[col_diff].values})
    else:
        gene_col = "Genes" if "Genes" in df_comparaison.columns else None
        if gene_col is None:
            return None
        df = df_comparaison[[gene_col, col_p, col_diff]].copy()
        df.columns = ["gene", "p", "lfc"]
    df = df.dropna(subset=["p", "lfc"])
    # Symbole propre : premier gène avant ';', majuscule, sans espaces
    df["gene"] = (df["gene"].astype(str).str.split(";").str[0]
                  .str.strip().str.upper())
    df = df[df["gene"].notna() & (df["gene"] != "") & (df["gene"] != "NAN")]
    if df.empty:
        return None

    # Métrique signée : -log10(p) * sign(LFC)
    p_clip = df["p"].clip(lower=1e-300)
    df["score"] = -np.log10(p_clip) * np.sign(df["lfc"])
    # En cas de doublons de symboles : garder le score de plus grande amplitude
    df["absscore"] = df["score"].abs()
    df = df.sort_values("absscore", ascending=False).drop_duplicates("gene")
    rnk = df.set_index("gene")["score"].sort_values(ascending=False)
    # Retirer les scores nuls (p=1 exactement) qui n'apportent rien
    rnk = rnk[rnk != 0]
    return rnk if len(rnk) >= 15 else None


def _plot_gsea_nes(res2d, contrast, out_dir, top_n=15):
    """Barplot horizontal des top voies par |NES| (rouge=up, bleu=down)."""
    try:
        df = res2d.copy()
        # Colonnes gseapy : 'Term', 'NES', 'FDR q-val'
        df["NES"] = pd.to_numeric(df["NES"], errors="coerce")
        df["FDR"] = pd.to_numeric(df.get("FDR q-val", np.nan), errors="coerce")
        df = df.dropna(subset=["NES"])
        if df.empty:
            return None
        df["absNES"] = df["NES"].abs()
        df = df.sort_values("absNES", ascending=False).head(top_n)
        df = df.sort_values("NES")
        colors = ["#3498db" if v < 0 else "#e74c3c" for v in df["NES"]]
        labels = [t[:45] + ("…" if len(t) > 45 else "") for t in df["Term"]]
        fig, ax = plt.subplots(figsize=(9, max(3, 0.42 * len(df))))
        ax.barh(range(len(df)), df["NES"], color=colors, edgecolor="black",
                linewidth=0.4)
        ax.set_yticks(range(len(df)))
        ax.set_yticklabels(labels, fontsize=8)
        ax.axvline(0, color="black", lw=0.6)
        ax.set_xlabel("NES (Normalized Enrichment Score)")
        ax.set_title(f"GSEA — {contrast}\n(red = enriched in up / blue = in down)")
        fig.tight_layout()
        f = os.path.join(out_dir, f"gsea_nes_{contrast}.png")
        fig.savefig(f, dpi=150)
        plt.close(fig)
        return f
    except Exception:
        return None


def run_gsea(df_comparaison, contrast_names, params, go_params, out_dir, orthologs=None):
    """
    GSEA pré-ranked par contraste.

    orthologs : {"organism": code g:Profiler, "map": {protein_id: orthologue}} venu de
    Perseverance ; utilisé si l'espèce de référence est celle des jeux g:Profiler.
    Retourne {contrast: {"table": df, "nes_plot": path}} ou None.
    Entièrement protégé : n'interrompt jamais le pipeline.
    """
    if go_params is None or not go_params.get("run_gsea", False):
        return None

    try:
        import gseapy
    except ImportError:
        print("  [WARN] GSEA skipped: gseapy not installed "
              "(pip install gseapy).")
        return None

    cfg       = gsea_settings(params, go_params)
    fdr_max   = cfg["fdr"]
    n_perm    = cfg["permutations"]
    p_lbl = ("p.adj" if cfg["rank_p"] == "volcano" and params.get("volcano_use_padj")
             else "p.value")

    print(f"\n[GSEA] Pre-ranked GSEA — metric: -log10({p_lbl}) x sign(LFC) — "
          f"{cfg['organism']} gene sets"
          + (f" (MSigDB {cfg['msigdb_version']})"
             if any(str(s).lower().startswith("msigdb:") for s in cfg["gene_sets"]) else ""))
    gene_sets = load_gene_sets(cfg["gene_sets"], cfg["organism"], cfg["msigdb_version"],
                               cache_dir=cfg["cache_dir"], go_obo=cfg["go_obo"])
    if not gene_sets:
        print("  [WARN] GSEA skipped: no gene set could be loaded.")
        return None
    key_map, mapping_txt = None, "proteins matched by gene symbol (Genes column)"
    if cfg["gprofiler"]:
        code = cfg["gprofiler"]
        orth = (orthologs or {}).get("map") if (orthologs or {}).get("organism") == code else None
        key_map, alias, summ = gprofiler_gene_keys(df_comparaison, code, orth)
        print(f"  [GSEA] {summ}")
        mapping_txt = summ.replace("[WARN] ", "")
        if not key_map:
            print("  [WARN] GSEA skipped: proteins could not be linked to gene sets.")
            return None
        gene_sets = translate_gene_sets(gene_sets, alias)
        if not gene_sets:
            print("  [WARN] GSEA skipped: no gene set member matches the data.")
            return None
    universe = {g for v in gene_sets.values() for g in v}
    key2prot = _key_proteins(df_comparaison, key_map)
    results = {}
    coverage_done = False
    for contrast in contrast_names:
        try:
            rnk = _ranking_metric(df_comparaison, contrast, params, key_map)
            if rnk is None:
                print(f"  [SKIP] {contrast}: not enough ranked genes.")
                continue
            if not coverage_done:
                n_in = int(rnk.index.isin(list(universe)).sum())
                print(f"  [GSEA] {n_in}/{len(rnk)} ranked genes found in the gene sets "
                      f"({100 * n_in / len(rnk):.0f} %)"
                      + ("  [WARN] low (expected only for a small custom collection): "
                         "check the Genes column / gsea_organism"
                         if n_in < 0.4 * len(rnk) else ""))
                coverage_done = True

            print(f"  [GSEA] {contrast}: {len(rnk)} ranked genes -> prerank...")
            pre = gseapy.prerank(
                rnk=rnk, gene_sets=gene_sets,
                min_size=10, max_size=1000,
                permutation_num=n_perm, seed=42,
                no_plot=True, verbose=False, threads=4,
                outdir=None)

            res = pre.res2d.copy()
            # Filtrer par FDR et trier
            res["FDR q-val"] = pd.to_numeric(res["FDR q-val"], errors="coerce")
            res["NES"] = pd.to_numeric(res["NES"], errors="coerce")
            res_sig = res[res["FDR q-val"] <= fdr_max].copy()
            res_sig = res_sig.sort_values("NES", ascending=False)

            if res_sig.empty:
                print(f"     No pathway below FDR {fdr_max}.")
                continue

            f_nes = _plot_gsea_nes(res_sig, contrast, out_dir)
            if "Lead_genes" in res_sig.columns:       # leading edge en protéines
                res_sig["Lead_proteins"] = res_sig["Lead_genes"].map(
                    lambda x: _lead_proteins(";".join(x) if isinstance(x, (list, tuple)) else x,
                                             key2prot))
            # Nettoyer les colonnes pour l'export (listes -> texte)
            for col in res_sig.columns:
                if res_sig[col].apply(lambda x: isinstance(x, (list, tuple))).any():
                    res_sig[col] = res_sig[col].apply(
                        lambda x: ";".join(map(str, x))
                        if isinstance(x, (list, tuple)) else x)

            results[contrast] = {
                "table": res_sig, "nes_plot": f_nes,
                "provenance": (f"GSEA pre-ranked (gseapy), metric -log10({p_lbl}) x sign(LFC), "
                               f"{n_perm} permutations, FDR <= {fdr_max}, set size 10-1000 · "
                               f"gene sets: {'; '.join(LAST_SOURCES) or ', '.join(cfg['gene_sets'])}"
                               f" · {mapping_txt}")}
            n_up = (res_sig["NES"] > 0).sum()
            n_dn = (res_sig["NES"] < 0).sum()
            print(f"     [OK] {len(res_sig)} pathways (FDR<={fdr_max}): "
                  f"{n_up} up, {n_dn} down.")

        except Exception as e:
            print(f"     [WARN] GSEA error on {contrast}: "
                  f"{type(e).__name__} ({str(e)[:70]}) — skipped.")
            continue

    if not results:
        print("  [INFO] No significant GSEA result across contrasts.")
        return None
    return results


def export_gsea_sheets(wb, gsea_results, ins_img_fn):
    """Ajoute un onglet GSEA_<contrast> par contraste au classeur Excel.

    Nommage sans 'GO' ni 'go' pour ne pas être capté par les filtres GO du
    dashboard (heatmap/réseau GO cross-contrastes). Noms uniques (31 caractères
    au plus, suffixe ~2, ~3... en cas de collision après troncature).
    Retourne {onglet: contraste} pour le _manifest (section gsea_sheet).
    """
    if not gsea_results:
        return {}
    sheet_map = {}
    used = {n.lower() for n in wb.sheetnames}
    for contrast, data in gsea_results.items():
        base = f"GSEA_{contrast.replace('_vs_', 'v').replace('_', '')}"[:31]
        sheet_name, k = base, 2
        while sheet_name.lower() in used:
            suf = f"~{k}"
            sheet_name, k = base[:31 - len(suf)] + suf, k + 1
        used.add(sheet_name.lower())
        sheet_map[sheet_name] = contrast
        ws = wb.create_sheet(sheet_name)
        df = data["table"]
        start_col = 16  # laisse de la place au plot NES à gauche
        for j, col in enumerate(df.columns):
            ws.cell(row=1, column=start_col + j, value=str(col))
        for i, row in enumerate(df.itertuples(index=False), start=2):
            for j, val in enumerate(row):
                if isinstance(val, float) and np.isnan(val):
                    val = None
                elif isinstance(val, (list, tuple, np.ndarray)):
                    val = ";".join(map(str, val))
                ws.cell(row=i, column=start_col + j, value=val)
        if data.get("provenance"):
            ws.cell(row=1, column=1, value=data["provenance"])
        ins_img_fn(ws, data.get("nes_plot"), "A2")
    return sheet_map
