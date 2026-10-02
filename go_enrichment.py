# ==============================================================================
# go_enrichment.py — Enrichissement GO/gProfiler optionnel (portage EnrichGO.R)
# Version Python — Proteogen
# ==============================================================================
# Portage fidèle de 260416-EnrichGO.R :
#   - gProfiler (gost) via le package officiel gprofiler-official
#   - sources : GO:BP, GO:CC, GO:MF, REAC, KEGG
#   - correction FDR, exclusion des termes racines, term_size 5–1000
#   - z-score = moyenne des LFC des protéines du terme
#   - graphiques : Manhattan, lollipop facetté, dotplot facetté, chord (GO:BP top5)
#   - export Excel : un onglet par contraste (données + plots)
#
# SEUILS : identiques à ceux des volcanos du pipeline principal (params).
# ESPÈCE : demandée à l'utilisateur (ex: hsapiens, rnorvegicus, mmusculus...).
# BACKEND : gProfiler (défaut) ou STRING (go_backend='string', espèce = NCBI
#   taxid) — STRING couvre la quasi-totalité des bactéries séquencées. Les
#   catégories STRING sont renommées en labels gProfiler (GO:BP, GO:MF, GO:CC,
#   REAC, KEGG) : plots, export et dashboard fonctionnent à l'identique.
#
# ENTIÈREMENT NON BLOQUANT : toute erreur (réseau, espèce inconnue, pas assez
# de protéines, API indisponible) est capturée et n'interrompt jamais le
# pipeline principal.
# ==============================================================================

import os
import re
import warnings
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

warnings.filterwarnings("ignore")

# Backend STRING (alternative à gProfiler) — utile pour les organismes mal
# couverts par Ensembl/gProfiler, en particulier les bactéries : STRING
# identifie l'espèce par NCBI Taxonomy ID (pas de convention de nommage
# fragile) et couvre >12000 organismes, procaryotes inclus.
STRING_API_URL = "https://string-db.org/api"
DEFAULT_STRING_CALLER = "Proteogen"


def _strip_common_affix(conditions):
    """Retire le préfixe ET le suffixe communs à TOUTES les conditions, pour ne
    garder que la partie discriminante. Identique à la fonction du dashboard
    (build_dashboardv7.py) — gardée synchronisée pour un nommage d'onglets et un
    appariement cohérents de part et d'autre.

    Ex: ['X18_009IG01Ctrl','X18_009IG01di6h','X18_009IG01di24h']
        -> {'...Ctrl':'Ctrl', '...di6h':'di6h', '...di24h':'di24h'}

    Garanties : jamais de label vide, labels uniques (suffixe #i si collision).
    """
    conds = list(dict.fromkeys(conditions))
    if len(conds) <= 1:
        return {c: (c[:10] if len(c) > 10 else c) for c in conds}

    pre = os.path.commonprefix(conds)
    suf = os.path.commonprefix([c[::-1] for c in conds])[::-1]

    def _trim(c):
        core = c
        if pre and core.startswith(pre):
            core = core[len(pre):]
        if suf and core.endswith(suf) and len(core) > len(suf):
            core = core[:len(core) - len(suf)]
        core = core.strip(" _-.")
        return core if core else c

    trimmed = {c: _trim(c) for c in conds}
    vals = list(trimmed.values())
    if any(not v for v in vals) or len(set(vals)) < len(vals):
        def _trim_pre_only(c):
            core = c[len(pre):] if (pre and c.startswith(pre)) else c
            core = core.strip(" _-.")
            return core if core else c
        trimmed2 = {c: _trim_pre_only(c) for c in conds}
        if len(set(trimmed2.values())) == len(trimmed2):
            trimmed = trimmed2

    out, seen = {}, {}
    for c, lbl in trimmed.items():
        s = lbl if len(lbl) <= 8 else lbl[:8]
        if s in seen.values():
            i = 2
            while f"{s[:6]}#{i}" in seen.values():
                i += 1
            s = f"{s[:6]}#{i}"
        seen[c] = s
        out[c] = s
    return out

# Sources d'enrichissement (comme le script R)
GO_SOURCES = ["GO:BP", "GO:CC", "GO:MF", "REAC", "KEGG"]

# Termes racines GO à exclure (molecular_function, biological_process, cellular_component)
GO_ROOT_TERMS = {"GO:0003674", "GO:0008150", "GO:0005575"}

# Palette par source pour le Manhattan plot
SOURCE_COLORS = {
    "GO:BP": "#FF7F0E", "GO:CC": "#2CA02C", "GO:MF": "#1F77B4",
    "REAC":  "#9467BD", "KEGG":  "#D62728", "TF": "#8C564B",
    "MIRNA": "#E377C2", "HPA": "#7F7F7F", "CORUM": "#BCBD22", "HP": "#17BECF",
}


# ==============================================================================
# 0. QUESTION INTERACTIVE (optionnelle)
# ==============================================================================

def ask_go_params() -> dict | None:
    """
    Demande si l'utilisateur veut lancer l'enrichissement GO et avec quelle espèce.
    Les SEUILS (p, ratio) ne sont PAS redemandés : ils proviennent des volcanos.
    Retourne None si refusé → le pipeline saute proprement cette étape.
    """
    print("\n" + "="*60)
    print("  GO ENRICHMENT (optional)")
    print("="*60)
    print("      The thresholds used are those defined for the volcanos.")

    rep = input("\n  Run GO enrichment? (y/N) -> ").strip().lower()
    if rep not in ("o", "oui", "y", "yes"):
        print("  [SKIP] GO enrichment skipped.\n")
        return None

    print("\n  Backend:")
    print("     [1] gProfiler (default) — good coverage for eukaryotes/models,")
    print("         patchy for bacteria. Species = gProfiler code.")
    print("     [2] STRING — species identified by NCBI Taxonomy ID, covers")
    print("         >12000 organisms incl. virtually all sequenced bacteria.")
    backend_choice = input("  Backend [1/2] (default 1) -> ").strip()

    if backend_choice == "2":
        print("\n  [!] Requires internet access to string-db.org.")
        print("      STRING indexes bacteria per STRAIN: use the strain taxid")
        print("      (E. coli K-12 MG1655 = 511145, not the species 562).")
        print("      Leave empty = auto-detect from your accessions (recommended).")
        taxid_raw = input("  NCBI Taxonomy ID [auto] -> ").strip().lower() or "auto"
        if taxid_raw != "auto" and not taxid_raw.isdigit():
            print("  [WARN] Invalid taxid — step skipped.\n")
            return None
        caller = input(
            f"  Caller identity for STRING (project/lab name, not personal "
            f"data) [{DEFAULT_STRING_CALLER}] -> ").strip() or DEFAULT_STRING_CALLER
        print("\n  Statistical background (universe):")
        print("     [1] whole genome (default — same as gProfiler 'known')")
        print("     [2] quantified proteins only (corrects the detectability")
        print("         bias of proteomics: abundant categories such as")
        print("         ribosome/translation are no longer over-called)")
        bg = "quantified" if input("  Background [1/2] (default 1) -> ").strip() == "2" \
            else "genome"
        params = {"backend": "string",
                  "species_taxid": "auto" if taxid_raw == "auto" else int(taxid_raw),
                  "caller_identity": caller, "string_background": bg,
                  "organism": None, "run_gsea": False}
    else:
        print("\n  [!] Requires internet access to g:Profiler (biit.cs.ut.ee).")
        print("  Example gProfiler species codes:")
        print("     hsapiens (human)     | mmusculus (mouse)   | rnorvegicus (rat)")
        print("     drerio (zebrafish)   | scerevisiae (yeast) | ecoli / efaecalis")
        organism = input("  gProfiler species code (e.g. rnorvegicus) -> ").strip()
        if not organism:
            print("  [WARN] No species provided — step skipped.\n")
            return None
        params = {"backend": "gprofiler", "organism": organism, "run_gsea": False}

    # GSEA rank-based (optionnel, complément de l'ORA)
    print("\n  --- GSEA (rank-based, optional) ---")
    print("  Complements ORA: ranks ALL proteins by -log10(p) x sign(LFC) and")
    print("  detects coordinated pathway shifts (even below the DEP threshold).")
    print("  [!] Requires 'gseapy' installed and internet access to Enrichr.")
    print("  [!] Gene sets are human-centric (Enrichr): most informative for")
    print("      human/mouse data mapped by gene symbol.")
    if params.get("backend") == "string":
        print("  [!] Independent of the STRING backend: for bacteria or other")
        print("      non-model species, Enrichr libraries do not apply — answer N.")
    rep_gsea = input("  Also run GSEA? (y/N) -> ").strip().lower()
    if rep_gsea in ("o", "oui", "y", "yes"):
        params["run_gsea"] = True

    return params


# ==============================================================================
# 1. APPEL gProfiler (gost)
# ==============================================================================

def _benjamini_hochberg(pvals: np.ndarray) -> np.ndarray:
    """BH standard, implémentation locale (cohérente avec le reste de
    l'écosystème Deimos qui réimplémente ses propres statistiques plutôt
    que d'ajouter une dépendance pour une fonction de quelques lignes)."""
    n = len(pvals)
    order = np.argsort(pvals)
    ranked = pvals[order] * n / (np.arange(n) + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    adj = np.empty(n)
    adj[order] = np.clip(ranked, 0, 1)
    return adj


def run_local_ora(query: list[str], background: list[str],
                  annotation_long: pd.DataFrame,
                  min_term_size: int = 3, source_label: str = "GO:PERSEVERANCE"
                  ) -> pd.DataFrame | None:
    """
    Remplace run_gost() par un test d'enrichissement local (hypergéométrique
    + BH), SANS appel réseau, à partir d'une annotation fournie par
    l'utilisateur — typiquement le transfert par orthologie de Perseverance pour
    un organisme non-modèle mal couvert par gProfiler.

    Le "background" (univers statistique) est ici l'ensemble des protéines
    RÉELLEMENT QUANTIFIÉES dans l'expérience Deimos en cours (toutes les
    lignes de la matrice, pas le protéome entier de l'espèce) — c'est le
    choix standard et recommandé pour l'ORA en protéomique (contrairement
    au génome entier utilisé par défaut par la plupart des outils GO), car
    seules les protéines quantifiables pouvaient de toute façon apparaître
    comme significatives.

    Parameters
    ----------
    query           : IDs des protéines significatives (sortie de
                      isolate_significant)
    background      : IDs de TOUTES les protéines quantifiées (univers)
    annotation_long : DataFrame avec au moins les colonnes 'protein_id' et
                      'GO' (GO en liste ';'-séparée OU déjà explosée, une
                      ligne par (protéine, terme))
    min_term_size   : taille minimum d'un terme dans le background pour
                      être testé (évite les termes anecdotiques)
    source_label    : valeur de la colonne 'source' en sortie (permet de
                      distinguer visuellement des sources gProfiler dans
                      les mêmes graphiques/exports)

    Returns
    -------
    DataFrame au même format que run_gost() (colonnes natives gProfiler),
    prêt à passer directement à process_enrichment(), ou None si rien de
    testable/significatif.
    """
    from scipy.stats import hypergeom

    df = annotation_long[["protein_id", "GO"]].dropna(subset=["GO"]).copy()
    df["GO"] = df["GO"].astype(str)
    if df["GO"].str.contains(";").any():
        df = df.assign(GO=df["GO"].str.split(";")).explode("GO")
    df["GO"] = df["GO"].str.strip()
    df = df[df["GO"] != ""]

    background_set = set(background)
    df = df[df["protein_id"].isin(background_set)]
    if df.empty:
        print("    [WARN] Aucune protéine du background n'a d'annotation "
              "GO fournie (Perseverance) — ORA local impossible.")
        return None

    annotated_background = set(df["protein_id"].unique())
    N = len(annotated_background)
    query_annotated = [p for p in query if p in annotated_background]
    n = len(query_annotated)
    if n == 0:
        print("    [INFO] Aucune protéine significative n'a d'annotation "
              "GO transférée — pas de test possible pour ce contraste.")
        return None
    query_set = set(query_annotated)

    rows = []
    for go_id, members in df.groupby("GO")["protein_id"].apply(set).items():
        K = len(members)
        if K < min_term_size:
            continue
        inter = members & query_set
        k = len(inter)
        if k == 0:
            continue
        pval = hypergeom.sf(k - 1, N, K, n)
        rows.append({
            "native": go_id, "name": go_id, "p_value": pval,
            "term_size": K, "query_size": n, "intersection_size": k,
            "source": source_label, "intersections": sorted(inter),
        })

    if not rows:
        return None

    res = pd.DataFrame(rows)
    res["p_value"] = _benjamini_hochberg(res["p_value"].values)
    res = res[res["p_value"] < 0.05].sort_values("p_value").reset_index(drop=True)
    return res if not res.empty else None


def diagnose_gprofiler_coverage(protein_ids: list, organism: str,
                                sample_size: int = 200, min_coverage: float = 0.30
                                ) -> dict:
    """
    Diagnostic de couverture gProfiler AVANT de lancer un enrichissement GO
    complet — utilise g:Convert (gp.convert), pas gost(), pour mesurer le
    taux de reconnaissance des accessions sans dépendre des résultats
    différentiels. Peut tourner tôt dans le pipeline, sur la liste brute des
    Protein.Group du run.

    Sert à objectiver la décision "gProfiler direct" vs "passer par une
    espèce de référence via Perseverance" — même philosophie que
    diagnose_deqms_reliability()/diagnose_batch_effect() : informer, pas
    décider silencieusement à la place de l'utilisateur.

    Parameters
    ----------
    protein_ids  : IDs bruts du run (groupes composites "A;B" nettoyés en
                   1er accession automatiquement)
    organism     : code organisme gProfiler à tester (ex: 'bnapus')
    sample_size  : g:Convert accepte de grandes listes, mais un échantillon
                   suffit pour un diagnostic rapide et limite la charge
                   réseau sur un run répété plusieurs fois en dev/debug
    min_coverage : taux de reconnaissance minimum jugé exploitable

    Returns
    -------
    dict : available (bool, False si package/réseau indisponible),
           organism_valid (bool), coverage (float ou None), reliable (bool),
           recommendation (str)
    """
    try:
        from gprofiler import GProfiler
    except ImportError:
        return {"available": False, "organism_valid": None, "coverage": None,
                "reliable": None,
                "recommendation": "gprofiler-official non installé — diagnostic impossible."}

    ids_clean = list(dict.fromkeys(str(p).split(";")[0].strip() for p in protein_ids))
    if len(ids_clean) > sample_size:
        rng = np.random.RandomState(0)
        sample = list(rng.choice(ids_clean, sample_size, replace=False))
    else:
        sample = ids_clean

    try:
        gp = GProfiler(return_dataframe=True)
        res = gp.convert(organism=organism, query=sample, target_namespace="ENSG")
    except Exception as e:
        return {"available": False, "organism_valid": False, "coverage": None,
                "reliable": False,
                "recommendation": (f"Échec gProfiler pour l'organisme '{organism}' "
                                   f"({type(e).__name__}) — organisme probablement "
                                   f"invalide ou inconnu de gProfiler. Perseverance "
                                   f"(orthologie vers une espèce de référence) recommandé.")}

    if res is None or len(res) == 0:
        return {"available": False, "organism_valid": True, "coverage": None,
                "reliable": False,
                "recommendation": "Réponse vide de g:Convert — diagnostic non concluant."}

    n_total = len(res)
    unmapped = res["converted"].isna() | res["converted"].astype(str).isin(["None", "N/A", "nan"])
    n_converted = int((~unmapped).sum())
    coverage = n_converted / n_total if n_total else 0.0
    reliable = coverage >= min_coverage

    reco = (f"{coverage:.0%} des accessions reconnues par gProfiler pour "
           f"'{organism}' (échantillon n={n_total}). ")
    reco += ("Couverture suffisante pour un enrichissement direct." if reliable
            else "Couverture trop faible — Perseverance (orthologie vers une "
                 "espèce mieux annotée) recommandé.")

    return {"available": True, "organism_valid": True, "n_sampled": n_total,
           "n_converted": n_converted, "coverage": round(coverage, 3),
           "reliable": reliable, "recommendation": reco}


# Catégories STRING -> labels gProfiler. Le dashboard (build_dashboardv7.py),
# les plots (SOURCE_COLORS, chord) et Pathfinder filtrent/colorent sur
# 'GO:BP','GO:MF','GO:CC','REAC','KEGG' : sans harmonisation, les résultats
# STRING seraient absents des onglets GO du dashboard.
STRING_CATEGORY_MAP = {"Process": "GO:BP", "Function": "GO:MF",
                       "Component": "GO:CC", "RCTM": "REAC", "KEGG": "KEGG"}
# Par défaut : mêmes sources que gProfiler (GO_SOURCES). Autres catégories
# STRING ajoutables via go_params['string_categories'] (conservent leur nom
# STRING) — utiles pour les bactéries : 'Keyword' (UniProt), 'Pfam',
# 'InterPro', 'SMART', 'WikiPathways'.
DEFAULT_STRING_CATEGORIES = ["Process", "Function", "Component", "KEGG", "RCTM"]

_STRING_CHUNK = 2000          # identifiants par requête get_string_ids
_STRING_NOSPECIES_MAX = 100   # l'API accepte get_string_ids sans espèce jusqu'à 100 IDs
_STRING_LAST_CALL = [0.0]     # horodatage du dernier appel (politesse API)


class StringNoMatch(Exception):
    """HTTP 404 de l'API STRING : aucun identifiant reconnu pour cette
    requête — typiquement un taxid absent de STRING (ex. taxid d'ESPÈCE
    bactérienne alors que STRING indexe les bactéries par SOUCHE), ou un
    format d'accession non indexé."""


_CONTAMINANT_PREFIXES = ("crap", "con_", "con__", "contam", "rev_", "decoy")


def _clean_ids_for_string(protein_ids: list) -> list:
    """1re accession du Protein.Group, dédupliquée, sans contaminants
    (cRAP-, CON__, REV_...) : inutile de les soumettre à STRING, et ils
    diluent la détection d'organisme."""
    out = []
    for p in protein_ids:
        acc = str(p).split(";")[0].strip()
        if not acc or acc.lower() == "nan":
            continue
        if acc.lower().startswith(_CONTAMINANT_PREFIXES):
            continue
        out.append(acc)
    return list(dict.fromkeys(out))


class StringUnreachable(Exception):
    """STRING n'a pas répondu (timeout, connexion) — à ne pas confondre avec
    « aucune correspondance »."""


def _string_post(method: str, data: dict, timeout: int = 120):
    """POST vers l'API STRING, avec au moins 1 s entre deux appels (politique
    d'usage STRING) et remontée explicite des messages d'erreur de l'API.
      - {"Error": ..., "ErrorMessage": ...}  -> ValueError (message STRING)
      - HTTP 404 / 400 sans corps json      -> StringNoMatch (aucune
        correspondance, ou taxid inconnu de STRING)
      - timeout / connexion                 -> StringUnreachable
    """
    import time
    import requests
    wait = 1.0 - (time.time() - _STRING_LAST_CALL[0])
    if wait > 0:
        time.sleep(wait)
    try:
        r = requests.post(f"{STRING_API_URL}/json/{method}", data=data,
                          timeout=timeout)
    except (requests.exceptions.Timeout, requests.exceptions.ConnectionError) as e:
        raise StringUnreachable(f"{type(e).__name__} after {timeout}s") from e
    finally:
        _STRING_LAST_CALL[0] = time.time()
    try:
        payload = r.json()
    except ValueError:
        payload = None
    if isinstance(payload, dict):
        raise ValueError(payload.get("ErrorMessage") or payload.get("Error")
                         or str(payload)[:120])
    if getattr(r, "status_code", 200) in (400, 404):
        raise StringNoMatch(f"HTTP {r.status_code} on '{method}' "
                            f"(no identifier matched / unknown taxid)")
    r.raise_for_status()
    return payload or []


def _as_list(x) -> list:
    """Champs multi-valués STRING : liste en json, chaîne 'a,b,c' en tsv
    (et selon les versions d'API). Gère les deux."""
    if isinstance(x, (list, tuple, np.ndarray)):
        return [str(v).strip() for v in x if str(v).strip()]
    if isinstance(x, str):
        return [v.strip() for v in x.split(",") if v.strip()]
    return []


def _string_id_map(identifiers: list, species_taxid,
                   caller_identity: str = DEFAULT_STRING_CALLER,
                   timeout: int = 120, raise_errors: bool = False
                   ) -> pd.DataFrame | None:
    """
    get_string_ids : résout des accessions vers les identifiants STRING
    (stringId, preferredName), avec écho de l'ID soumis (queryItem).
    Requêtes découpées par blocs de _STRING_CHUNK identifiants.

    species_taxid : taxid STRING (= NCBI taxid ; pour les bactéries, celui de
                    la SOUCHE, ex: 511145 = E. coli K-12 MG1655, pas 562).
                    None = recherche tous organismes (LENT côté STRING :
                    réservé au repli de detect_string_taxid, peu d'IDs).
    raise_errors  : True -> StringUnreachable / autres erreurs remontées à
                    l'appelant (pour distinguer « pas de réponse » de « rien
                    trouvé ») ; False -> None + avertissement.
    Retourne None si aucune correspondance.
    """
    try:
        import requests  # noqa: F401
    except ImportError:
        return None
    ids = [str(i) for i in identifiers]
    step = _STRING_CHUNK if species_taxid is not None else _STRING_NOSPECIES_MAX
    frames = []
    try:
        for k in range(0, len(ids), step):
            payload = {"identifiers": "\r".join(ids[k:k + step]),
                       "limit": 1, "echo_query": 1,
                       "caller_identity": caller_identity}
            if species_taxid is not None:
                payload["species"] = species_taxid
            try:
                data = _string_post("get_string_ids", payload, timeout=timeout)
            except StringNoMatch:
                continue          # aucun ID de ce bloc reconnu : bloc suivant
            if data:
                frames.append(pd.DataFrame(data))
    except Exception as e:
        if raise_errors:
            raise
        print(f"    [WARN] STRING get_string_ids failed "
              f"({type(e).__name__}: {str(e)[:100]})")
        return None
    if not frames:
        return None
    return pd.concat(frames, ignore_index=True)


_UNIPROT_ACC = re.compile(r"^([OPQ][0-9][A-Z0-9]{3}[0-9]|"
                          r"[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2})(-\d+)?$")


def _uniprot_organisms(accessions: list, timeout: int = 30) -> pd.DataFrame | None:
    """Organisme (taxid NCBI de la souche + nom) d'accessions UniProt, via
    l'API REST UniProt — une requête, ~1 s. None si échec ou rien trouvé."""
    try:
        import requests
    except ImportError:
        return None
    accs = [a.split("-")[0] for a in accessions]
    query = " OR ".join(f"accession:{a}" for a in accs)
    try:
        r = requests.get("https://rest.uniprot.org/uniprotkb/search",
                         params={"query": query, "size": len(accs),
                                 "fields": "accession,organism_id,organism_name",
                                 "format": "json"}, timeout=timeout)
        r.raise_for_status()
        results = r.json().get("results", [])
    except Exception as e:
        print(f"    [INFO] UniProt injoignable ({type(e).__name__}) — "
              f"repli sur la recherche STRING tous organismes.")
        return None
    rows = [{"acc": x.get("primaryAccession"),
             "taxid": (x.get("organism") or {}).get("taxonId"),
             "name": (x.get("organism") or {}).get("scientificName")}
            for x in results]
    df = pd.DataFrame(rows).dropna(subset=["taxid"]) if rows else None
    return df if df is not None and not df.empty else None


def detect_string_taxid(protein_ids: list,
                        caller_identity: str = DEFAULT_STRING_CALLER,
                        min_share: float = 0.5) -> dict:
    """
    Identifie l'organisme STRING des accessions du run.

    1. Accessions UniProt (cas DIA-NN sur un protéome UniProt) : UniProt donne
       le taxid de SOUCHE de chaque accession (1 requête, ~1 s) ; ce taxid est
       ensuite vérifié dans STRING (get_string_ids avec espèce, rapide).
       Si STRING ne le connaît pas -> la souche est absente de STRING (message
       explicite, avec le nom de la souche).
    2. Repli (accessions non-UniProt, ou UniProt injoignable) : get_string_ids
       SANS espèce sur 10 accessions — recherche tous organismes, lente côté
       STRING (d'où l'échantillon réduit et le timeout de 90 s).

    Returns
    -------
    dict : taxid (int ou None), name, method, share, others, reliable,
           unreachable (bool), recommendation.
    """
    ids_clean = _clean_ids_for_string(protein_ids)
    base = {"taxid": None, "name": None, "method": None, "share": 0.0,
            "others": {}, "reliable": False, "unreachable": False}
    if not ids_clean:
        return {**base, "recommendation": "Aucune accession exploitable."}
    rng = np.random.RandomState(0)

    # --- 1. Voie UniProt ---
    uni = [a for a in ids_clean if _UNIPROT_ACC.match(a)]
    if len(uni) >= 5:
        sample = list(rng.choice(uni, min(25, len(uni)), replace=False))
        org = _uniprot_organisms(sample)
        if org is not None:
            counts = org["taxid"].astype(int).value_counts()
            tx = int(counts.index[0])
            name = str(org.loc[org["taxid"].astype(int) == tx, "name"].iloc[0])
            share_u = float(counts.iloc[0]) / len(sample)
            others = {int(k): int(v) for k, v in counts.iloc[1:4].items()}
            print(f"  [DIAG] UniProt : {int(counts.iloc[0])}/{len(sample)} accessions "
                  f"-> {name} (taxid {tx}). Vérification dans STRING...")
            check = list(rng.choice(uni, min(50, len(uni)), replace=False))
            try:
                mapped = _string_id_map(check, tx, caller_identity, timeout=60,
                                        raise_errors=True)
            except StringUnreachable:
                return {**base, "taxid": tx, "name": name, "method": "uniprot",
                        "strain_taxid": tx, "unreachable": True, "recommendation": (
                            f"Souche identifiée par UniProt : {name} (taxid {tx}), "
                            f"mais STRING n'a pas répondu — présence dans STRING "
                            f"non vérifiée. Réessayer plus tard ou forcer "
                            f"go_species_taxid: {tx}.")}
            except Exception as e:
                mapped = None
                print(f"    [WARN] STRING : {type(e).__name__} ({str(e)[:80]})")
            n_ok = (mapped["queryItem"].nunique()
                    if mapped is not None and "queryItem" in mapped.columns else 0)
            share = n_ok / len(check)
            if share >= 0.30:
                return {**base, "taxid": tx, "name": name, "method": "uniprot",
                        "strain_taxid": tx,
                        "share": round(share, 3), "others": others,
                        "reliable": share_u >= min_share,
                        "recommendation": (f"{name} (taxid {tx}) : {n_ok}/{len(check)} "
                                           f"accessions reconnues par STRING.")}
            return {**base, "name": name, "method": "uniprot", "others": others,
                    "strain_taxid": tx, "recommendation": (
                        f"Vos accessions appartiennent à {name} (taxid {tx}, "
                        f"UniProt), mais STRING n'en reconnaît que {n_ok}/"
                        f"{len(check)} pour ce taxid : souche absente de STRING "
                        f"(ou accessions non indexées). Il faut passer par une "
                        f"souche présente dans STRING (orthologie).")}

    # --- 2. Repli : STRING tous organismes, petit échantillon ---
    sample = list(rng.choice(ids_clean, min(10, len(ids_clean)), replace=False))
    print(f"  [DIAG] Recherche STRING tous organismes sur {len(sample)} accessions "
          f"(lent côté STRING, jusqu'à ~1-2 min)...")
    try:
        mapped = _string_id_map(sample, None, caller_identity, timeout=90,
                                raise_errors=True)
    except StringUnreachable as e:
        return {**base, "unreachable": True, "recommendation": (
            f"STRING n'a pas répondu ({e}). Organisme non déterminé : renseigner "
            f"go_species_taxid (taxid de SOUCHE, cf. la fiche UniProt/NCBI d'une "
            f"de vos protéines).")}
    except Exception as e:
        return {**base, "recommendation": f"Erreur STRING ({type(e).__name__}: {str(e)[:80]})."}
    if mapped is None or "stringId" not in mapped.columns:
        return {**base, "method": "string", "recommendation": (
            "Aucune accession reconnue par STRING, tous organismes confondus : "
            "format d'accession non indexé (ex. identifiants custom/Prokka) ou "
            "souche absente de STRING.")}
    key = "queryItem" if "queryItem" in mapped.columns else "stringId"
    hits = mapped.drop_duplicates(key).copy()
    hits["_taxid"] = hits["stringId"].astype(str).str.split(".").str[0]
    counts = hits["_taxid"].value_counts()
    best = counts.index[0]
    name = None
    if "taxonName" in hits.columns:
        nm = hits.loc[hits["_taxid"] == best, "taxonName"].dropna()
        name = str(nm.iloc[0]) if len(nm) else None
    share = float(counts.iloc[0]) / len(sample)
    label = f"{best}" + (f" ({name})" if name else "")
    reco = (f"{int(counts.iloc[0])}/{len(sample)} accessions échantillonnées "
            f"appartiennent au taxid STRING {label}.")
    if share < min_share:
        reco += (" Part trop faible pour conclure (accessions de type nom de "
                 "gène, ou organisme absent de STRING).")
    return {**base, "taxid": int(best), "name": name, "method": "string",
            "share": round(share, 3),
            "others": {int(k): int(v) for k, v in counts.iloc[1:4].items()},
            "reliable": share >= min_share, "recommendation": reco}


def resolve_string_taxid(protein_ids: list, species_taxid,
                         caller_identity: str = DEFAULT_STRING_CALLER,
                         interactive: bool = True, min_coverage: float = 0.30
                         ) -> "tuple[int | None, str | None]":
    """
    Pré-vol STRING, à lancer AVANT les statistiques (échec rapide) :
      - species_taxid 'auto'/None : détection depuis les accessions ;
      - taxid explicite : couverture mesurée ; si elle est nulle ou faible,
        détection et proposition du taxid détecté (question, défaut Y).
    Returns (taxid retenu ou None si GO STRING impossible, nom du taxon).
    """
    auto = species_taxid in (None, "", "auto", "AUTO", "Auto")
    if not auto:
        diag = diagnose_string_coverage(protein_ids, int(species_taxid), caller_identity,
                                        min_coverage=min_coverage)
        print(f"  [DIAG] {diag['recommendation']}")
        if diag.get("reliable"):
            return int(species_taxid), None
        if diag.get("unreachable"):
            print("  [WARN] STRING injoignable — taxid conservé tel quel, "
                  "le GO réessaiera en fin de pipeline.")
            return int(species_taxid), None

    print("  [DIAG] Détection de l'organisme depuis les accessions...")
    det = detect_string_taxid(protein_ids, caller_identity)
    print(f"  [DIAG] {det['recommendation']}")
    if det["others"]:
        print(f"         Autres taxids trouvés : {det['others']}")

    if auto:
        if det["reliable"]:
            print(f"  [AUTO] go_species_taxid = {det['taxid']}"
                  + (f" ({det['name']})" if det["name"] else ""))
            return det["taxid"], det["name"]
        if det.get("unreachable") and det.get("taxid"):
            print(f"  [AUTO] go_species_taxid = {det['taxid']} (non vérifié dans "
                  f"STRING) — le GO réessaiera en fin de pipeline.")
            return det["taxid"], det["name"]
        print("  [WARN] Organisme STRING non déterminé — GO STRING désactivé "
              "pour ce run (le reste du pipeline continue).")
        return None, None

    # Taxid explicite défaillant
    if det["reliable"] and det["taxid"] != int(species_taxid):
        print(f"  [!] Le taxid {species_taxid} ne correspond pas à vos accessions "
              f"(pour une bactérie, STRING attend le taxid de SOUCHE).")
        if interactive:
            rep = input(f"  Utiliser le taxid détecté {det['taxid']} à la place ? "
                        f"[Y/n] -> ").strip().lower()
            if rep in ("", "y", "yes", "o", "oui"):
                return det["taxid"], det["name"]
        else:
            print(f"  -> Relancer avec go_species_taxid: {det['taxid']} "
                  f"(ou go_species_taxid: auto).")
    return int(species_taxid), None


# ------------------------------------------------------------------------------
# Mode orthologie STRING : l'organisme existe dans STRING mais pas vos
# accessions (autre souche indexée) -> RBH DIAMOND vers le protéome STRING.
# ------------------------------------------------------------------------------

_GENE_SYMBOL = re.compile(r"^[a-z]{3}[A-Z]?[0-9]?$")   # dnaK, gyrB, prs, rpsP, groL...


def _uniprot_taxon(taxid, timeout: int = 20) -> dict | None:
    """Fiche taxonomie UniProt (rank, scientificName, parent)."""
    try:
        import requests
        r = requests.get(f"https://rest.uniprot.org/taxonomy/{int(taxid)}",
                         params={"format": "json"}, timeout=timeout)
        r.raise_for_status()
        return r.json()
    except Exception as e:
        print(f"    [INFO] Taxonomie UniProt injoignable pour {taxid} ({type(e).__name__}).")
        return None


def _uniprot_parent_species(taxid):
    """Remonte la taxonomie UniProt jusqu'au rang 'species' (souche -> espèce).
    Returns (taxid, nom) ou (None, None)."""
    tx = taxid
    for _ in range(5):
        js = _uniprot_taxon(tx)
        if not js:
            return None, None
        if js.get("rank") == "species":
            return int(js["taxonId"]), js.get("scientificName")
        parent = js.get("parent") or {}
        if not parent.get("taxonId"):
            return None, None
        tx = parent["taxonId"]
    return None, None


def probe_string_organism(gene_names, taxid,
                          caller_identity: str = DEFAULT_STRING_CALLER,
                          n_probe: int = 30):
    """
    Vérifie qu'un organisme existe dans STRING à partir de noms de gènes
    (dnaK, gyrB...) — utile quand les accessions du run n'y sont pas indexées
    (autre souche). Returns (n_reconnus, n_testés) ; n_reconnus None si STRING
    ne répond pas, n_testés 0 si aucun nom de gène exploitable.
    """
    names = list(dict.fromkeys(str(g).split(";")[0].strip() for g in gene_names))
    names = [g for g in names if _GENE_SYMBOL.match(g)]
    if not names:
        return 0, 0
    if len(names) > n_probe:
        names = list(np.random.RandomState(0).choice(names, n_probe, replace=False))
    try:
        mapped = _string_id_map(names, int(taxid), caller_identity, timeout=60,
                                raise_errors=True)
    except StringUnreachable:
        return None, len(names)
    except Exception:
        return 0, len(names)
    n = (mapped["queryItem"].nunique()
         if mapped is not None and "queryItem" in mapped.columns else 0)
    return int(n), len(names)


def _diamond_available() -> bool:
    import shutil
    return (shutil.which("diamond") is not None
            or os.path.exists("diamond.exe") or os.path.exists("diamond"))


def plan_string_mapping(protein_ids: list, species_taxid,
                        caller_identity: str = DEFAULT_STRING_CALLER,
                        gene_names=None, reference_taxid=None,
                        _from_explicit: bool = False) -> dict:
    """
    Pré-vol complet du backend STRING (AVANT les statistiques). Décide
    comment relier les accessions du run à STRING :
      - 'direct'   : accessions indexées dans STRING pour le taxid retenu ;
      - 'ortholog' : l'organisme est dans STRING mais pas vos accessions (STRING
                     indexe une autre souche) -> RBH DIAMOND vers le protéome
                     STRING après les statistiques ;
      - None       : GO STRING impossible (message explicite).

    reference_taxid : organisme STRING à utiliser en mode orthologie (config
                      string_reference_taxid). Par défaut : l'espèce parente
                      (taxonomie UniProt) de la souche de vos accessions.
    Returns dict : mode, taxid, name, strain_taxid, strain_name.
    """
    out = {"mode": None, "taxid": None, "name": None,
           "strain_taxid": None, "strain_name": None}
    genes = [g for g in (gene_names if gene_names is not None else [])]
    auto = species_taxid in (None, "", "auto", "AUTO", "Auto")

    if not auto:
        T = int(species_taxid)
        diag = diagnose_string_coverage(protein_ids, T, caller_identity)
        print(f"  [DIAG] {diag['recommendation']}")
        if diag.get("reliable"):
            return {**out, "mode": "direct", "taxid": T}
        if diag.get("unreachable"):
            print("  [WARN] STRING injoignable — taxid conservé, non vérifié.")
            return {**out, "mode": "direct", "taxid": T}
        cand, cand_name = (int(reference_taxid) if reference_taxid else T), None
    else:
        print("  [DIAG] Détection de l'organisme depuis les accessions...")
        det = detect_string_taxid(protein_ids, caller_identity)
        print(f"  [DIAG] {det['recommendation']}")
        if det["reliable"]:
            print(f"  [AUTO] go_species_taxid = {det['taxid']}"
                  + (f" ({det['name']})" if det["name"] else ""))
            return {**out, "mode": "direct", "taxid": det["taxid"], "name": det["name"],
                    "strain_taxid": det.get("strain_taxid"), "strain_name": det["name"]}
        if det.get("unreachable") and det.get("taxid"):
            print(f"  [AUTO] go_species_taxid = {det['taxid']} (non vérifié dans STRING).")
            return {**out, "mode": "direct", "taxid": det["taxid"], "name": det["name"]}
        out.update(strain_taxid=det.get("strain_taxid"), strain_name=det.get("name"))
        if reference_taxid:
            cand, cand_name = int(reference_taxid), None
        elif det.get("strain_taxid"):
            cand, cand_name = _uniprot_parent_species(det["strain_taxid"])
            if cand:
                print(f"  [DIAG] Espèce parente (taxonomie UniProt) : {cand_name} "
                      f"(taxid {cand}) — recherche dans STRING...")
        else:
            cand, cand_name = None, None
        if cand is None:
            print("  [WARN] Organisme STRING non déterminé — GO STRING désactivé "
                  "(renseigner string_reference_taxid).")
            return out

    # --- L'organisme de référence existe-t-il dans STRING ? ---
    n, m = probe_string_organism(genes, cand, caller_identity)
    if n is None:
        print(f"  [WARN] STRING n'a pas répondu — présence du taxid {cand} non vérifiée.")
    elif m == 0:
        print(f"  [INFO] Pas de noms de gènes exploitables pour vérifier le taxid "
              f"{cand} ; le téléchargement du protéome STRING tranchera.")
    elif n == 0:
        if not auto and not _from_explicit:
            print(f"  [DIAG] Le taxid {cand} ne correspond à aucun organisme STRING "
                  f"(0/{m} gènes reconnus) — détection automatique...")
            return plan_string_mapping(protein_ids, "auto", caller_identity,
                                       gene_names, reference_taxid, _from_explicit=True)
        print(f"  [WARN] Taxid {cand} absent de STRING (0/{m} gènes reconnus). "
              f"GO STRING désactivé : renseigner string_reference_taxid (une "
              f"souche proche présente dans STRING).")
        return out
    else:
        print(f"  [DIAG] Organisme STRING {cand}"
              + (f" ({cand_name})" if cand_name else "")
              + f" présent : {n}/{m} noms de gènes reconnus.")

    if not _diamond_available():
        print("  [WARN] Mode orthologie requis, mais DIAMOND est introuvable "
              "(ni dans le PATH ni dans le dossier courant). GO STRING désactivé "
              "pour ce run. Installation : github.com/bbuchfink/diamond/releases "
              "(Windows : diamond.exe à côté de deimos.py, ou dans le PATH).")
        return out

    src = out.get("strain_name") or "vos accessions"
    print(f"  [PLAN] Accessions non indexées dans STRING -> orthologie (RBH DIAMOND) "
          f"{src} -> protéome STRING {cand}, après les statistiques.")
    return {**out, "mode": "ortholog", "taxid": cand, "name": cand_name}


def _string_version(caller_identity: str = DEFAULT_STRING_CALLER) -> str:
    """Version courante de STRING (endpoint 'version'), pour télécharger le
    protéome de la MÊME version que l'API d'enrichissement."""
    try:
        data = _string_post("version", {"caller_identity": caller_identity}, timeout=30)
        if data and data[0].get("string_version"):
            return str(data[0]["string_version"])
    except Exception:
        pass
    return "12.0"


def download_string_proteome(taxid, out_dir: str,
                             caller_identity: str = DEFAULT_STRING_CALLER) -> str | None:
    """Télécharge (une fois, puis cache) le protéome STRING d'un organisme :
    <out_dir>/<taxid>.protein.sequences.v<version>.fa. Les en-têtes sont les
    identifiants STRING ('<taxid>.<protéine>'), donc directement utilisables
    comme IDs de référence du RBH. None si échec (message avec la marche à
    suivre manuelle)."""
    import gzip
    version = _string_version(caller_identity)
    fname = f"{taxid}.protein.sequences.v{version}.fa"
    out = os.path.join(out_dir, fname)
    if os.path.exists(out) and os.path.getsize(out) > 0:
        print(f"  [OK] Protéome STRING déjà présent : {out}")
        return out
    try:
        import requests
    except ImportError:
        return None
    for base in ("https://stringdb-downloads.org/download",
                 "https://stringdb-static.org/download"):
        url = f"{base}/protein.sequences.v{version}/{fname}.gz"
        try:
            r = requests.get(url, timeout=300)
            if r.status_code != 200:
                continue
            text = gzip.decompress(r.content).decode("utf-8", errors="replace")
        except Exception:
            continue
        if not text.startswith(f">{taxid}."):
            continue
        os.makedirs(out_dir, exist_ok=True)
        with open(out, "w", encoding="utf-8") as fh:
            fh.write(text)
        n = text.count("\n>") + 1
        print(f"  [OK] Protéome STRING {taxid} (v{version}) téléchargé : {n} séquences.")
        return out
    print(f"  [WARN] Téléchargement du protéome STRING {taxid} impossible. "
          f"Manuellement : string-db.org > Download > choisir l'organisme > "
          f"'{fname}.gz', décompresser dans le dossier des données (ou "
          f"string_reference_fasta: chemin).")
    return None


def download_uniprot_fasta(accessions: list, out_path: str, batch: int = 100) -> int:
    """Séquences UniProt des accessions (repli quand le FASTA de recherche
    DIA-NN n'est pas trouvé localement). Returns le nombre de séquences."""
    try:
        import requests
    except ImportError:
        return 0
    accs = [a for a in dict.fromkeys(accessions) if _UNIPROT_ACC.match(str(a))]
    n = 0
    with open(out_path, "w", encoding="utf-8") as fh:
        for k in range(0, len(accs), batch):
            q = " OR ".join(f"accession:{a}" for a in accs[k:k + batch])
            try:
                r = requests.get("https://rest.uniprot.org/uniprotkb/stream",
                                 params={"query": q, "format": "fasta"}, timeout=120)
                r.raise_for_status()
            except Exception as e:
                print(f"  [WARN] UniProt : lot {k // batch + 1} non téléchargé "
                      f"({type(e).__name__}).")
                continue
            fh.write(r.text if r.text.endswith("\n") else r.text + "\n")
            n += r.text.count(">")
    return n


def string_maps_from_orthologs(ortholog_map: dict):
    """Adapte un ortholog_map {accession: stringId} (RBH) au format de
    build_string_maps, pour réutiliser _string_enrich_query."""
    acc2string, string2acc = dict(ortholog_map), {}
    for acc, sid in ortholog_map.items():
        for key in (str(sid), str(sid).split(".", 1)[-1]):
            string2acc.setdefault(key, []).append(acc)
    return acc2string, string2acc, len(ortholog_map)


def build_string_maps(protein_ids: list, species_taxid: int,
                      caller_identity: str = DEFAULT_STRING_CALLER):
    """
    Résout UNE FOIS toutes les accessions du run (convention Deimos : 1re
    accession du Protein.Group) vers STRING.

    Returns
    -------
    (acc2string, string2acc, n_total) ou None si échec.
      acc2string : {accession: stringId}
      string2acc : {clé STRING: [accessions]} — clés = stringId, stringId
                   sans préfixe taxon, et preferredName (repli), pour
                   retraduire les intersections quel que soit le format
                   renvoyé par l'endpoint enrichment.
    """
    ids_clean = _clean_ids_for_string(protein_ids)
    if not ids_clean:
        return None
    mapped = _string_id_map(ids_clean, species_taxid, caller_identity)
    if mapped is None or not {"queryItem", "stringId"} <= set(mapped.columns):
        return None

    acc2string, string2acc = {}, {}
    # pandas >= 3 : astype(str) conserve les NaN (float) -> tests explicites
    names = (mapped["preferredName"] if "preferredName" in mapped.columns
             else pd.Series([None] * len(mapped)))
    for q, s, nm in zip(mapped["queryItem"], mapped["stringId"], names):
        if not isinstance(q, str) or not isinstance(s, str) or not s:
            continue
        if q in acc2string:
            continue
        acc2string[q] = s
        for key in (s, s.split(".", 1)[-1]):
            string2acc.setdefault(key, []).append(q)
        if isinstance(nm, str) and nm:
            string2acc.setdefault(nm, [])
            if q not in string2acc[nm]:
                string2acc[nm].append(q)
    return acc2string, string2acc, len(ids_clean)


def diagnose_string_coverage(protein_ids: list, species_taxid: int,
                             caller_identity: str = DEFAULT_STRING_CALLER,
                             sample_size: int = 200, min_coverage: float = 0.30
                             ) -> dict:
    """
    Équivalent diagnose_gprofiler_coverage() pour le backend STRING — même
    philosophie (informer avant l'enrichissement complet, ne pas décider à
    la place de l'utilisateur), via get_string_ids plutôt que g:Convert.
    """
    try:
        import requests  # noqa: F401
    except ImportError:
        return {"available": False, "reliable": None, "coverage": None,
                "recommendation": "Package 'requests' non installé — "
                                   "diagnostic STRING impossible."}

    ids_clean = _clean_ids_for_string(protein_ids)
    if len(ids_clean) > sample_size:
        rng = np.random.RandomState(0)
        sample = list(rng.choice(ids_clean, sample_size, replace=False))
    else:
        sample = ids_clean

    try:
        mapped = _string_id_map(sample, species_taxid, caller_identity,
                                timeout=60, raise_errors=True)
    except StringUnreachable as e:
        return {"available": False, "reliable": False, "coverage": None,
                "unreachable": True,
                "recommendation": f"STRING n'a pas répondu ({e}) — couverture "
                                  f"du taxid {species_taxid} non vérifiée."}
    except Exception as e:
        mapped = None
        print(f"    [WARN] STRING : {type(e).__name__} ({str(e)[:80]})")
    if mapped is None:
        return {"available": False, "reliable": False, "coverage": 0.0,
                "recommendation": (f"Aucune accession reconnue par STRING pour le "
                                   f"taxid {species_taxid} — taxid absent de STRING "
                                   f"(bactéries : STRING attend le taxid de SOUCHE, "
                                   f"pas d'espèce), ou pas de réseau vers "
                                   f"string-db.org.")}

    n_total = len(sample)
    n_mapped = (int(mapped.loc[mapped["stringId"].notna(), "queryItem"].nunique())
                if {"stringId", "queryItem"} <= set(mapped.columns) else 0)
    coverage = n_mapped / n_total if n_total else 0.0
    reliable = coverage >= min_coverage
    reco = (f"{coverage:.0%} des accessions résolues par STRING pour le "
           f"taxid {species_taxid} (échantillon n={n_total}). ")
    reco += ("Couverture suffisante pour un enrichissement direct." if reliable
            else "Couverture faible — vérifier le taxid (souche exacte ?) ou "
                 "le type d'accession fourni (UniProt vs locus tag).")
    return {"available": True, "n_sampled": n_total, "n_mapped": n_mapped,
           "coverage": round(coverage, 3), "reliable": reliable,
           "recommendation": reco}


def run_string_enrichment(query: list[str], species_taxid: int,
                          caller_identity: str = DEFAULT_STRING_CALLER,
                          categories: "list[str] | None" = None,
                          background: "list[str] | None" = None
                          ) -> pd.DataFrame | None:
    """
    Équivalent run_gost(), backend STRING (endpoint 'enrichment').

    query      : identifiants STRING de préférence (résolus via
                 build_string_maps) — recommandé par STRING, sans ambiguïté.
    categories : catégories STRING à garder (défaut DEFAULT_STRING_CATEGORIES).
    background : identifiants STRING formant l'univers statistique (ex. les
                 protéines quantifiées). None = génome complet de l'espèce
                 (défaut STRING, équivalent domain_scope='known' de gProfiler).

    Retourne un DataFrame au contrat run_gost() (native/name/p_value/
    term_size/query_size/intersection_size/source/intersections) + gene_names.
    'source' est harmonisé sur les labels gProfiler (STRING_CATEGORY_MAP).
    'intersections' contient les identifiants STRING : la retraduction vers
    les accessions d'origine est faite par l'appelant (_string_enrich_query).
    """
    try:
        import requests  # noqa: F401
    except ImportError:
        print("    [WARN] Package 'requests' not installed. STRING GO step skipped.")
        return None

    payload = {"identifiers": "\r".join(str(q) for q in query),
               "species": species_taxid,
               "caller_identity": caller_identity}
    if background:
        payload["background_string_identifiers"] = "\r".join(background)
    try:
        data = _string_post("enrichment", payload, timeout=180)
    except Exception as e:
        print(f"    [WARN] STRING failed ({type(e).__name__}: {str(e)[:100]})")
        return None
    if not data:
        return None
    df = pd.DataFrame(data)
    if "category" not in df.columns:
        return None

    cats = categories or DEFAULT_STRING_CATEGORIES
    df = df[df["category"].isin(cats)].copy()
    if df.empty:
        return None

    # STRING renvoie 'p_value' (brut) ET 'fdr' (BH) : on garde le fdr comme
    # 'p_value' du contrat (process_enrichment attend un p déjà corrigé,
    # comme gProfiler en mode 'fdr' et run_local_ora) — et on retire le brut
    # AVANT le renommage, sinon deux colonnes 'p_value' coexistent.
    df = df.drop(columns=["p_value"], errors="ignore")
    inter_col = "inputGenes" if "inputGenes" in df.columns else "preferredNames"
    df["intersections"] = df[inter_col].apply(_as_list)
    df["gene_names"] = (df["preferredNames"].apply(lambda x: ",".join(_as_list(x)))
                        if "preferredNames" in df.columns else "")
    df["category"] = df["category"].map(lambda c: STRING_CATEGORY_MAP.get(c, c))
    df = df.rename(columns={
        "term": "native", "description": "name", "fdr": "p_value",
        "number_of_genes": "intersection_size",
        "number_of_genes_in_background": "term_size", "category": "source",
    })
    df["p_value"] = pd.to_numeric(df["p_value"], errors="coerce")
    df = df[df["p_value"] < 0.05]
    if df.empty:
        return None
    df["query_size"] = len(query)
    keep = ["native", "name", "p_value", "term_size", "query_size",
           "intersection_size", "source", "intersections", "gene_names"]
    return df[[c for c in keep if c in df.columns]].reset_index(drop=True)


def _string_enrich_query(prots: list, string_maps, species_taxid: int,
                         caller_identity: str, categories=None,
                         background=None, label: str = "") -> pd.DataFrame | None:
    """Traduit accessions -> stringId, lance l'enrichissement STRING, puis
    retraduit les intersections vers les accessions d'origine (restreintes
    à la requête), pour que le z-score de process_enrichment() retrouve les
    LFC. Même logique que la branche ortholog_map (Perseverance)."""
    acc2string, string2acc, _ = string_maps
    translated = list(dict.fromkeys(acc2string[p] for p in prots if p in acc2string))
    n_dropped = sum(1 for p in prots if p not in acc2string)
    if n_dropped:
        print(f"  [GO] {label}: {n_dropped}/{len(prots)} protéine(s) non "
              f"résolue(s) par STRING -> exclue(s) de la requête.")
    if len(translated) <= 5:
        print(f"  [SKIP] {label}: too few STRING-resolved proteins "
              f"({len(translated)}), skipped.")
        return None
    print(f"  [GO] {label}: {len(translated)} proteins -> STRING "
          f"(taxid {species_taxid})...")
    gost_df = run_string_enrichment(translated, species_taxid, caller_identity,
                                    categories=categories, background=background)
    if gost_df is None or gost_df.empty:
        return gost_df

    qset = set(prots)

    def _back(inter):
        out = []
        for x in _as_list(inter):
            accs = [a for a in string2acc.get(x, []) if a in qset]
            out.extend(accs if accs else [x])
        return list(dict.fromkeys(out))

    gost_df = gost_df.copy()
    gost_df["intersections"] = gost_df["intersections"].apply(_back)
    return gost_df


def _prepare_string(protein_ids: list, species_taxid, caller_identity: str,
                    background_mode: str = "genome"):
    """Résolution unique des IDs + background, avec log de couverture.
    Retourne (string_maps, background) ou (None, None) si inutilisable."""
    if species_taxid in (None, "", "None", "auto", "AUTO", "Auto"):
        # Normalement déjà résolu par le pré-vol de deimos.main() ; filet de
        # sécurité si run_go_enrichment est appelé directement.
        species_taxid, _ = resolve_string_taxid(protein_ids, "auto", caller_identity,
                                                interactive=False)
        if species_taxid is None:
            return None, None
    maps = build_string_maps(protein_ids, int(species_taxid), caller_identity)
    if maps is None:
        print(f"  [WARN] STRING : aucune accession reconnue pour le taxid "
              f"{species_taxid} (taxid absent de STRING — pour une bactérie, "
              f"utiliser le taxid de SOUCHE ou go_species_taxid: auto — ou pas "
              f"de réseau vers string-db.org) — GO ignoré.")
        return None, None
    acc2string, _, n_total = maps
    cov = len(acc2string) / n_total if n_total else 0.0
    print(f"   STRING : {len(acc2string)}/{n_total} accessions résolues "
          f"({cov:.0%}) pour le taxid {species_taxid}.")
    if cov < 0.30:
        print("   [WARN] Couverture faible — vérifier le taxid (souche exacte) "
              "ou le type d'accession (UniProt vs locus tag).")
    background = None
    if background_mode == "quantified":
        background = sorted(set(acc2string.values()))
        print(f"   Background : {len(background)} protéines quantifiées "
              f"(résolues STRING).")
    else:
        print("   Background : génome complet (défaut STRING, équivalent "
              "domain_scope='known' de gProfiler).")
    return maps, background


def run_gost(query: list[str], organism: str) -> pd.DataFrame | None:
    """
    Appelle gProfiler (équivalent de gost() en R).
    Retourne un DataFrame de résultats ou None en cas d'échec.

    Paramètres (équivalents du gost() R, adaptés à l'API Python) :
      user_threshold=0.05, no_iea=False (= exclude_iea=FALSE),
      significance_threshold_method='fdr' (= correction_method='fdr'),
      domain_scope='known', no_evidences=False (= evcodes=TRUE).
    """
    try:
        from gprofiler import GProfiler
    except ImportError:
        print("    [WARN] Package 'gprofiler-official' not installed "
              "(pip install gprofiler-official). GO step skipped.")
        return None

    try:
        gp = GProfiler(return_dataframe=True)
        res = gp.profile(
            organism=organism,
            query=query,
            sources=GO_SOURCES,
            user_threshold=0.05,
            no_evidences=False,       # evcodes=TRUE → récupère les intersections
            no_iea=False,             # exclude_iea=FALSE
            domain_scope="known",
            significance_threshold_method="fdr",   # correction_method='fdr'
        )
        if res is None or len(res) == 0:
            return None
        return res
    except Exception as e:
        print(f"    [WARN] gProfiler failed ({type(e).__name__}: {str(e)[:80]})")
        return None


# ==============================================================================
# 2. ISOLATION DES PROTÉINES SIGNIFICATIVES PAR CONTRASTE
# ==============================================================================

def isolate_significant(df_comparaison: pd.DataFrame, contrast: str,
                        params: dict) -> tuple[list[str], pd.DataFrame]:
    """
    Réplique l'isolation du script R (idx_sig) :
    - filtre p < seuil ET |LFC| > lfc_min   (seuils = volcanos)
    - tri par p croissante
    - IDs nettoyés (1er ID avant ';'), dédupliqués
    Retourne (liste protéines, df_local des stats du contraste).

    NB : le script R utilise p.val brute (col_p = _p.val). On respecte ici le
    choix volcano de l'utilisateur (p.val ou p.adj) pour rester cohérent avec
    le pipeline principal.
    """
    suffixe_p = "_p.adj" if params["volcano_use_padj"] else "_p.val"
    col_p    = f"{contrast}{suffixe_p}"
    col_diff = f"{contrast}_diff"

    if col_p not in df_comparaison.columns or col_diff not in df_comparaison.columns:
        return [], pd.DataFrame()

    id_col = "name" if "name" in df_comparaison.columns else df_comparaison.columns[0]
    seuil_p = params["volcano_p_thresh"]
    lfc_min = params["volcano_lfc_min"]

    df_local = df_comparaison[[id_col, col_p, col_diff]].copy()
    df_local.columns = ["name", "p_val", "diff"]
    df_local = df_local.dropna(subset=["p_val", "diff"])
    df_local = df_local[(df_local["p_val"] < seuil_p) &
                        (df_local["diff"].abs() > lfc_min)]
    df_local = df_local.sort_values("p_val")

    # Nettoyage IDs (1er avant ';'), déduplication en gardant l'ordre
    df_local["name_clean"] = (df_local["name"].astype(str)
                              .str.split(";").str[0].str.strip())
    prots_clean = df_local["name_clean"].drop_duplicates().tolist()

    return prots_clean, df_local


# ==============================================================================
# 3. POST-TRAITEMENT DU TABLEAU D'ENRICHISSEMENT (z-score, filtres)
# ==============================================================================

def process_enrichment(gost_df: pd.DataFrame, df_local: pd.DataFrame
                        ) -> pd.DataFrame:
    """
    Réplique le calcul tab_enrich du script R :
    - exclusion des termes racines GO
    - filtre term_size 5–1000
    - gene_ratio = intersection_size / query_size
    - z_score = moyenne des LFC des protéines du terme
    - tri par (source, p_value)
    """
    df = gost_df.copy()

    # Harmonisation des noms de colonnes gProfiler Python → noms du script R
    rename = {
        "native": "term_id", "name": "term_name", "p_value": "p_value",
        "term_size": "term_size", "query_size": "query_size",
        "intersection_size": "intersection_size", "source": "source",
        "intersections": "intersection",   # liste des gènes (evcodes)
    }
    for old, new in rename.items():
        if old in df.columns and old != new:
            df = df.rename(columns={old: new})

    # Exclusion termes racines
    if "term_id" in df.columns:
        df = df[~df["term_id"].isin(GO_ROOT_TERMS)]

    # Filtre term_size
    if "term_size" in df.columns:
        df = df[(df["term_size"] >= 5) & (df["term_size"] <= 1000)]

    if df.empty:
        return df

    # gene_ratio
    if "intersection_size" in df.columns and "query_size" in df.columns:
        df["gene_ratio"] = df["intersection_size"] / df["query_size"].replace(0, np.nan)

    # Map nom → LFC pour le z-score (df_local peut être vide ou sans name_clean)
    if "name_clean" in df_local.columns and "diff" in df_local.columns:
        lfc_map = dict(zip(df_local["name_clean"], df_local["diff"]))
    else:
        lfc_map = {}

    def _zscore(intersection):
        # gProfiler renvoie soit une liste, soit une string séparée par ','
        if isinstance(intersection, (list, tuple, np.ndarray)):
            prots = [str(p) for p in intersection]
        elif isinstance(intersection, str):
            prots = [p.strip() for p in intersection.split(",")]
        else:
            return 0.0
        prots = [p.split(";")[0] for p in prots]
        vals = [lfc_map[p] for p in prots if p in lfc_map]
        return float(np.mean(vals)) if vals else 0.0

    if "intersection" in df.columns:
        df["z_score"] = df["intersection"].apply(_zscore)
    else:
        df["intersection"] = ""
        df["z_score"] = 0.0

    # --- Garantir la présence de TOUTES les colonnes utilisées en aval ---
    # (les plots et l'export y accèdent directement ; on évite tout KeyError
    #  si gProfiler n'a pas renvoyé une colonne attendue).
    defaults = {
        "term_id": "", "term_name": "", "source": "GO:BP",
        "p_value": 1.0, "term_size": 0, "query_size": 0,
        "intersection_size": 0, "gene_ratio": 0.0, "z_score": 0.0,
        "intersection": "",
    }
    for col, val in defaults.items():
        if col not in df.columns:
            df[col] = val

    sort_cols = [c for c in ["source", "p_value"] if c in df.columns]
    df = df.sort_values(sort_cols).reset_index(drop=True)
    return df


# ==============================================================================
# 4. GRAPHIQUES
# ==============================================================================

def plot_manhattan(tab: pd.DataFrame, contrast: str, out_dir: str) -> str | None:
    """
    Manhattan plot façon gostplot : -log10(p) par terme, groupé par source.
    Les 2 termes les plus significatifs par source sont annotés.
    """
    if tab.empty or "source" not in tab.columns:
        return None

    df = tab.copy()
    df["logp"] = -np.log10(df["p_value"].clip(lower=1e-300))
    sources = [s for s in GO_SOURCES if s in df["source"].unique()]
    sources += [s for s in df["source"].unique() if s not in sources]

    total_pts = len(df)
    # Largeur adaptative : ni trop étalé (peu de points) ni trop serré
    fig_w = float(np.clip(total_pts * 0.45 + 3, 7, 16))
    fig, ax = plt.subplots(figsize=(fig_w, 6))
    x_offset = 0
    xticks, xlabels = [], []
    texts = []

    for src in sources:
        sub = df[df["source"] == src].sort_values("p_value")
        n = len(sub)
        if n == 0:
            continue
        xs = np.arange(x_offset, x_offset + n)
        color = SOURCE_COLORS.get(src, "#666666")
        ax.scatter(xs, sub["logp"].values, c=color, s=28, alpha=0.75,
                   label=src, edgecolor="none")
        xticks.append(x_offset + n / 2)
        xlabels.append(src)
        # Top 2 par source à annoter (texte court, repoussé ensuite)
        for rank in range(min(2, n)):
            nm = str(sub["term_name"].values[rank])
            nm = nm if len(nm) <= 32 else nm[:30] + "…"
            texts.append(ax.text(xs[rank], sub["logp"].values[rank], nm,
                                 fontsize=6.5))
        x_offset += n + max(1, n // 8)

    # Anti-chevauchement des labels (adjustText si dispo)
    try:
        from adjustText import adjust_text
        adjust_text(texts, ax=ax, only_move={"text": "y"},
                    arrowprops=dict(arrowstyle="-", color="grey", lw=0.3))
    except ImportError:
        pass

    ax.axhline(-np.log10(0.05), ls="--", color="grey", lw=0.8)
    ax.set_xticks(xticks)
    ax.set_xticklabels(xlabels, fontsize=10, fontweight="bold")
    ax.set_ylabel("-log10(p-value FDR)")
    ax.set_title(f"Manhattan GO/gProfiler — {contrast.replace('_vs_', ' vs ')}")
    # Marges : éviter que points/labels touchent les bords
    ax.set_xlim(-1, x_offset)
    ymax = df["logp"].max()
    ax.set_ylim(0, ymax * 1.15 + 0.5)
    ax.legend(fontsize=8, loc="upper right", ncol=max(1, len(sources)))
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    f = os.path.join(out_dir, f"go_manhattan_{contrast}.png")
    fig.savefig(f, dpi=150, bbox_inches="tight")
    plt.close(fig)
    return f


def _top_terms_per_source(tab: pd.DataFrame, n_per_source: int = 8) -> pd.DataFrame:
    """Top n termes par source, avec term_name tronqué pour l'affichage."""
    df = tab.copy()
    # Garantir les colonnes de base AVANT tout traitement
    if "source" not in df.columns:
        df["source"] = "GO"
    if "p_value" not in df.columns:
        df["p_value"] = 1.0
    df["p_value"] = pd.to_numeric(df["p_value"], errors="coerce").fillna(1.0)

    # Top n par source via boucle (évite les pièges de groupby.apply selon
    # la version de pandas, qui pouvaient faire disparaître la colonne 'source')
    parts = []
    for src in df["source"].dropna().unique():
        sub = df[df["source"] == src].nsmallest(n_per_source, "p_value")
        parts.append(sub)
    top = (pd.concat(parts, ignore_index=True) if parts
           else df.head(0).copy())

    # Garantir toutes les colonnes utilisées par les plots
    top["minus_log10_p"] = -np.log10(top["p_value"].clip(lower=1e-300))
    if "z_score" not in top.columns:
        top["z_score"] = 0.0
    if "gene_ratio" not in top.columns:
        top["gene_ratio"] = 0.0
    if "intersection_size" not in top.columns:
        top["intersection_size"] = 1
    if "term_name" in top.columns:
        top["term_short"] = top["term_name"].astype(str).apply(
            lambda s: s if len(s) <= 45 else s[:43] + "…")
    else:
        top["term_short"] = [f"term_{i}" for i in range(len(top))]
    return top


def plot_lollipop_faceted(tab: pd.DataFrame, contrast: str, out_dir: str) -> str | None:
    """Lollipop facetté par source, couleur = z-score (LFC moyen)."""
    if tab.empty:
        return None
    top = _top_terms_per_source(tab, 8)
    sources = [s for s in GO_SOURCES if s in top["source"].unique()]
    sources += [s for s in top["source"].unique() if s not in sources]
    n_src = len(sources)
    if n_src == 0:
        return None

    n_src = len(sources)
    if n_src == 0:
        return None

    # Layout vertical (1 colonne) : évite tout chevauchement des labels longs.
    fig, axes = plt.subplots(n_src, 1,
                             figsize=(9, max(2.4, 2.2 * n_src)),
                             squeeze=False)
    axes = axes.flatten()

    zmax = max(abs(top["z_score"]).max(), 1e-6)
    for idx, src in enumerate(sources):
        ax = axes[idx]
        sub = top[top["source"] == src].sort_values("minus_log10_p")
        ypos = np.arange(len(sub))
        ax.hlines(ypos, 0, sub["minus_log10_p"], color="grey", lw=1)
        sc = ax.scatter(sub["minus_log10_p"], ypos, c=sub["z_score"],
                        cmap="RdBu_r", vmin=-zmax, vmax=zmax, s=90,
                        edgecolor="black", lw=0.5, zorder=3)
        ax.set_yticks(ypos)
        ax.set_yticklabels(sub["term_short"], fontsize=8)
        ax.set_xlabel("-log10(p)", fontsize=9)
        ax.set_title(src, fontweight="bold", fontsize=11, loc="left")
        ax.set_xlim(0, sub["minus_log10_p"].max() * 1.15 + 0.3)
        ax.margins(y=0.12)
        ax.spines[["top", "right"]].set_visible(False)

    for idx in range(n_src, len(axes)):
        axes[idx].set_visible(False)

    sm = plt.cm.ScalarMappable(cmap="RdBu_r",
                               norm=plt.Normalize(vmin=-zmax, vmax=zmax))
    cbar = fig.colorbar(sm, ax=axes.tolist(), fraction=0.025, pad=0.04)
    cbar.set_label("Z-score (LFC moyen)", fontsize=9)
    fig.suptitle(f"Lollipop GO — {contrast.replace('_vs_', ' vs ')}",
                 fontsize=13, fontweight="bold", y=1.0)
    fig.subplots_adjust(left=0.32, right=0.86, hspace=0.55, top=0.90)
    f = os.path.join(out_dir, f"go_lollipop_{contrast}.png")
    fig.savefig(f, dpi=150, bbox_inches="tight")
    plt.close(fig)
    return f


def plot_dotplot_faceted(tab: pd.DataFrame, contrast: str, out_dir: str) -> str | None:
    """Dotplot facetté : gene_ratio en x, taille = intersection_size, couleur = p."""
    if tab.empty:
        return None
    top = _top_terms_per_source(tab, 8)
    sources = [s for s in GO_SOURCES if s in top["source"].unique()]
    sources += [s for s in top["source"].unique() if s not in sources]
    n_src = len(sources)
    if n_src == 0:
        return None

    n_src = len(sources)
    if n_src == 0:
        return None

    # Layout vertical (1 colonne) : les noms de termes GO sont longs, empiler
    # les sources évite tout chevauchement horizontal entre facettes.
    fig, axes = plt.subplots(n_src, 1,
                             figsize=(9, max(2.4, 2.2 * n_src)),
                             squeeze=False)
    axes = axes.flatten()

    # Échelle de couleur commune à toutes les facettes
    all_logp = -np.log10(top["p_value"].clip(lower=1e-300))
    vmin, vmax = float(all_logp.min()), float(all_logp.max())
    sc = None
    for idx, src in enumerate(sources):
        ax = axes[idx]
        sub = top[top["source"] == src].sort_values("gene_ratio")
        ypos = np.arange(len(sub))
        sizes = sub["intersection_size"].values * 28 if "intersection_size" in sub.columns else 60
        sc = ax.scatter(sub["gene_ratio"], ypos, s=sizes,
                        c=-np.log10(sub["p_value"].clip(lower=1e-300)),
                        cmap="plasma", vmin=vmin, vmax=vmax,
                        edgecolor="black", lw=0.4, zorder=3)
        ax.set_yticks(ypos)
        ax.set_yticklabels(sub["term_short"], fontsize=8)
        ax.set_xlabel("Gene Ratio", fontsize=9)
        ax.set_title(src, fontweight="bold", fontsize=11, loc="left")
        ax.margins(x=0.18, y=0.12)
        ax.grid(axis="x", ls=":", color="grey", alpha=0.4)
        ax.spines[["top", "right"]].set_visible(False)

    for idx in range(n_src, len(axes)):
        axes[idx].set_visible(False)

    # Colorbar unique partagée à droite, bien détachée des panneaux
    if sc is not None:
        cbar = fig.colorbar(sc, ax=axes.tolist(), fraction=0.025, pad=0.04)
        cbar.set_label("-log10(p)", fontsize=9)
    # Légende des tailles (Count) sur la 1re facette
    if "intersection_size" in top.columns:
        for cnt in sorted(top["intersection_size"].unique())[:5]:
            axes[0].scatter([], [], s=cnt * 28, c="grey",
                            edgecolor="black", lw=0.4, label=str(int(cnt)))
        axes[0].legend(title="Count", loc="upper left",
                       bbox_to_anchor=(1.01, 1.0), fontsize=7,
                       title_fontsize=8, framealpha=0.9)

    fig.suptitle(f"Dotplot GO — {contrast.replace('_vs_', ' vs ')}",
                 fontsize=13, fontweight="bold", y=1.0)
    fig.subplots_adjust(left=0.32, right=0.86, hspace=0.55, top=0.90)
    f = os.path.join(out_dir, f"go_dotplot_{contrast}.png")
    fig.savefig(f, dpi=150, bbox_inches="tight")
    plt.close(fig)
    return f


def plot_chord(tab: pd.DataFrame, df_local: pd.DataFrame,
               contrast: str, out_dir: str) -> str | None:
    """
    Chord plot (GO:BP top 5) version matplotlib : matrice protéines × termes,
    protéines ordonnées par LFC, barre LFC à gauche (GOChord-like).
    """
    if tab.empty:
        return None
    bp = tab[tab["source"] == "GO:BP"].nsmallest(5, "p_value")
    if bp.empty:
        return None

    # Récupérer les protéines de chaque terme
    def _prots(intersection):
        if isinstance(intersection, (list, tuple, np.ndarray)):
            return [str(p).split(";")[0] for p in intersection]
        elif isinstance(intersection, str):
            return [p.strip().split(";")[0] for p in intersection.split(",")]
        return []

    term_prots = {row["term_name"]: _prots(row["intersection"])
                  for _, row in bp.iterrows() if "intersection" in bp.columns}
    all_prots = sorted(set().union(*term_prots.values())) if term_prots else []
    lfc_map = dict(zip(df_local["name_clean"], df_local["diff"]))
    all_prots = [p for p in all_prots if p in lfc_map]
    if len(all_prots) <= 2:
        return None

    terms = list(term_prots.keys())
    # Matrice : NaN si gène absent du terme, sinon LFC du gène (pour la couleur)
    mat = np.full((len(all_prots), len(terms)), np.nan)
    for j, t in enumerate(terms):
        for p in term_prots[t]:
            if p in all_prots:
                mat[all_prots.index(p), j] = lfc_map[p]

    lfc_vals = np.array([lfc_map[p] for p in all_prots])
    order = np.argsort(lfc_vals)
    mat = mat[order]
    prots_ordered = [all_prots[i] for i in order]
    lfc_ordered = lfc_vals[order]

    fig, (ax_lfc, ax_mat) = plt.subplots(
        1, 2, figsize=(11, max(6, len(all_prots) * 0.26)),
        gridspec_kw={"width_ratios": [0.16, 1]})

    vmax = max(np.abs(lfc_ordered).max(), 0.5) if len(lfc_ordered) else 1
    # Barre LFC à gauche
    ax_lfc.barh(range(len(prots_ordered)), lfc_ordered,
                color=plt.cm.RdBu_r((lfc_ordered / (2 * vmax)) + 0.5),
                edgecolor="grey", linewidth=0.3)
    ax_lfc.set_yticks(range(len(prots_ordered)))
    ax_lfc.set_yticklabels(prots_ordered, fontsize=6)
    ax_lfc.set_xlabel("logFC", fontsize=8)
    ax_lfc.axvline(0, color="black", lw=0.5)
    ax_lfc.set_ylim(-0.5, len(prots_ordered) - 0.5)
    ax_lfc.invert_yaxis()
    ax_lfc.spines[["top", "right"]].set_visible(False)

    # Matrice colorée par LFC (cases vides = gène hors du terme)
    masked = np.ma.masked_invalid(mat)
    cmap = plt.cm.RdBu_r.copy()
    cmap.set_bad(color="#f0f0f0")   # gris clair pour les cases vides
    im = ax_mat.imshow(masked, aspect="auto", cmap=cmap,
                       vmin=-vmax, vmax=vmax, interpolation="nearest")
    # Grille fine pour séparer les cases
    ax_mat.set_xticks(np.arange(-0.5, len(terms), 1), minor=True)
    ax_mat.set_yticks(np.arange(-0.5, len(prots_ordered), 1), minor=True)
    ax_mat.grid(which="minor", color="white", linewidth=1.2)
    ax_mat.set_xticks(range(len(terms)))
    wrapped = ["\n".join([t[k:k+22] for k in range(0, min(len(t), 44), 22)])
               for t in terms]
    ax_mat.set_xticklabels(wrapped, rotation=40, ha="right", fontsize=7)
    ax_mat.set_yticks([])
    ax_mat.set_title(f"Genes <-> GO:BP terms (top 5) — "
                     f"{contrast.replace('_vs_', ' vs ')}\n"
                     f"color = gene log2FC", fontsize=9)
    cbar = fig.colorbar(im, ax=ax_mat, fraction=0.025, pad=0.02)
    cbar.set_label("log2FC", fontsize=8)

    fig.tight_layout()
    f = os.path.join(out_dir, f"go_chord_{contrast}.png")
    fig.savefig(f, dpi=150, bbox_inches="tight")
    plt.close(fig)
    return f


# ==============================================================================
# 5. POINT D'ENTRÉE PRINCIPAL
# ==============================================================================

def run_go_enrichment(df_comparaison: pd.DataFrame, contrast_names: list[str],
                      params: dict, go_params: dict, out_dir: str) -> dict | None:
    """
    Lance l'enrichissement GO sur tous les contrastes.
    Entièrement protégé : retourne None si tout échoue, sans jamais lever
    d'exception qui interromprait le pipeline principal.

    Seuils de significativité = ceux des volcanos (params).
    Espèce = go_params['organism'].

    Retourne {contrast: {"table": df, "manhattan":.., "lollipop":.., "dotplot":.., "chord":..}}
    """
    if go_params is None:
        return None

    try:
        organism = go_params.get("organism")
        custom_annotation = go_params.get("custom_annotation")   # DataFrame Perseverance (annotation transferee), ou None
        ortholog_map = go_params.get("ortholog_map")              # dict protein_id -> ortholog_id (RBH), ou None
        backend = go_params.get("backend", "gprofiler")
        species_taxid = go_params.get("species_taxid")
        string_caller = go_params.get("caller_identity") or DEFAULT_STRING_CALLER
        string_categories = go_params.get("string_categories")
        string_bg_mode = go_params.get("string_background", "genome")
        background = go_params.get("background")
        if background is None and custom_annotation is not None:
            id_col = "name" if "name" in df_comparaison.columns else df_comparaison.columns[0]
            background = (df_comparaison[id_col].astype(str)
                          .str.split(";").str[0].tolist())

        use_local = custom_annotation is not None
        use_ortholog = (ortholog_map is not None) and not use_local
        use_string = (backend == "string") and not use_local and not use_ortholog
        use_string_orth = use_ortholog and backend == "string"
        if use_local:
            print(f"\n[GO] Enrichissement LOCAL (annotation custom — Perseverance/orthologue), "
                  f"pas d'appel gProfiler.")
        elif use_string_orth:
            print(f"\n[GO] STRING via orthologues (RBH DIAMOND) — protéome STRING "
                  f"taxid {species_taxid}")
        elif use_ortholog:
            print(f"\n[GO] gProfiler via orthologues (Perseverance RBH) — "
                  f"espèce de référence '{organism}'")
        elif use_string:
            print(f"\n[GO] STRING enrichment — NCBI taxid {species_taxid}")
        else:
            print(f"\n[GO] GO/gProfiler enrichment — species '{organism}'")
        print(f"   Seuils (volcanos) : "
              f"{'p.adj' if params['volcano_use_padj'] else 'p.val'} < "
              f"{params['volcano_p_thresh']} | ratio ≥ {params['volcano_ratio_min']}")

        reverse_ortholog_map = ({v: k for k, v in ortholog_map.items()}
                                if use_ortholog else None)

        # STRING : résolution UNIQUE de toutes les accessions du run (et du
        # background éventuel), réutilisée pour tous les contrastes.
        string_maps, string_bg = None, None
        if use_string_orth:
            string_maps = string_maps_from_orthologs(ortholog_map)
            n_q = df_comparaison.shape[0]
            print(f"   {len(ortholog_map)}/{n_q} protéines avec un orthologue STRING.")
            if string_bg_mode == "quantified":
                string_bg = sorted(set(ortholog_map.values()))
                print(f"   Background : {len(string_bg)} orthologues STRING des "
                      f"protéines quantifiées.")
            else:
                print("   Background : génome complet (défaut STRING).")
        elif use_string:
            id_col = "name" if "name" in df_comparaison.columns else df_comparaison.columns[0]
            string_maps, string_bg = _prepare_string(
                df_comparaison[id_col].astype(str).tolist(), species_taxid,
                string_caller, string_bg_mode)
            if string_maps is None:
                return None

        results = {}
        for contrast in contrast_names:
            try:
                prots, df_local = isolate_significant(df_comparaison, contrast, params)

                if len(prots) <= 5:
                    print(f"  [SKIP] {contrast}: too few proteins ({len(prots)}), skipped.")
                    continue

                if use_local:
                    print(f"  [GO] {contrast}: {len(prots)} proteins -> ORA local (Perseverance)...")
                    gost_df = run_local_ora(prots, background, custom_annotation)
                elif use_string_orth:
                    gost_df = _string_enrich_query(
                        prots, string_maps, int(species_taxid), string_caller,
                        categories=string_categories, background=string_bg,
                        label=contrast)
                elif use_ortholog:
                    translated = [ortholog_map[p] for p in prots if p in ortholog_map]
                    n_dropped = len(prots) - len(translated)
                    if n_dropped:
                        print(f"  [GO] {contrast}: {n_dropped}/{len(prots)} protéine(s) "
                              f"sans orthologue -> exclue(s) de la requête gProfiler.")
                    if len(translated) <= 5:
                        print(f"  [SKIP] {contrast}: too few translatable orthologs "
                              f"({len(translated)}), skipped.")
                        continue
                    print(f"  [GO] {contrast}: {len(translated)} orthologue(s) -> "
                          f"gProfiler ({organism})...")
                    gost_df = run_gost(translated, organism)
                    if gost_df is not None and not gost_df.empty and "intersections" in gost_df.columns:
                        # Retraduit les IDs de l'espece de reference vers les IDs
                        # d'origine AVANT process_enrichment, sinon son calcul de
                        # z_score (lookup par ID d'origine) echoue silencieusement
                        # pour tout le monde (0.0 partout).
                        gost_df = gost_df.copy()
                        gost_df["intersections"] = gost_df["intersections"].apply(
                            lambda inter: [reverse_ortholog_map.get(x, x) for x in inter]
                            if isinstance(inter, (list, tuple, np.ndarray)) else inter)
                elif use_string:
                    gost_df = _string_enrich_query(
                        prots, string_maps, int(species_taxid), string_caller,
                        categories=string_categories, background=string_bg,
                        label=contrast)
                else:
                    print(f"  [GO] {contrast}: {len(prots)} proteins -> gProfiler...")
                    gost_df = run_gost(prots, organism)
                if gost_df is None or gost_df.empty:
                    print(f"     Aucun enrichissement significatif.")
                    continue

                tab = process_enrichment(gost_df, df_local)
                if tab.empty:
                    print(f"     No term after filtering (term_size, roots).")
                    continue

                # Chaque plot est isolé : un échec graphique ne doit pas faire
                # perdre le tableau d'enrichissement (l'essentiel pour l'Excel).
                def _safe_plot(fn, *a):
                    try:
                        return fn(*a)
                    except Exception as ex:
                        import traceback as _tb
                        detail = f"clé {ex}" if isinstance(ex, KeyError) else str(ex)[:80]
                        print(f"        (plot {fn.__name__} skipped: "
                              f"{type(ex).__name__} — {detail})")
                        line = _tb.format_exc().strip().splitlines()[-2:]
                        for l in line:
                            print(f"          {l.strip()}")
                        return None
                f_lol = _safe_plot(plot_lollipop_faceted, tab, contrast, out_dir)
                f_dot = _safe_plot(plot_dotplot_faceted, tab, contrast, out_dir)

                # Nettoyer colonnes complexes (listes) avant export Excel
                tab_export = tab.copy()
                if "intersection" in tab_export.columns:
                    tab_export["intersection"] = tab_export["intersection"].apply(
                        lambda x: ",".join(map(str, x)) if isinstance(x, (list, tuple, np.ndarray))
                        else str(x))
                # Retirer colonnes d'affichage internes
                tab_export = tab_export.drop(columns=["term_short", "minus_log10_p"],
                                             errors="ignore")

                results[contrast] = {
                    "table": tab_export,
                    "lollipop": f_lol, "dotplot": f_dot,
                }
                n_bp = (tab["source"] == "GO:BP").sum()   # STRING 'Process' déjà harmonisé en 'GO:BP'
                if use_local:
                    lbl = "Perseverance/orthologue (ORA local)"
                elif use_string_orth:
                    lbl = f"{n_bp} GO:BP (STRING {species_taxid}, via orthologues)"
                elif use_ortholog:
                    lbl = f"{n_bp} GO:BP, via orthologues {organism}"
                elif use_string:
                    lbl = f"{n_bp} GO:BP (STRING, taxid {species_taxid})"
                else:
                    lbl = f"{n_bp} GO:BP"
                print(f"     [OK] {len(tab)} enriched terms ({lbl}).")

            except Exception as e:
                import traceback
                detail = str(e)
                if isinstance(e, KeyError):
                    detail = f"missing column/key: {e}"
                print(f"     [WARN] Error on {contrast}: {type(e).__name__} "
                      f"({detail}) — skipped.")
                # Trace courte pour diagnostic (1ère ligne utile)
                tb = traceback.format_exc().strip().splitlines()
                for line in tb[-3:]:
                    print(f"        {line.strip()}")
                continue

        if not results:
            print("  [INFO] No significant GO enrichment across all contrasts.")
            return None

        return results

    except Exception as e:
        print(f"  [WARN] GO enrichment interrupted ({type(e).__name__}) — "
              f"the main pipeline is not affected.")
        return None


def run_go_enrichment_clusters(cluster_mapping, params: dict, go_params: dict,
                               out_dir: str) -> dict | None:
    """
    Enrichissement GO par cluster de la heatmap ANOVA.

    Contrairement à run_go_enrichment (qui part des DEP par contraste), on prend
    ici les protéines regroupées par profil d'expression (clusters k-means de la
    heatmap). Background = protéome connu de l'espèce (domain_scope='known').

    Parameters
    ----------
    cluster_mapping : liste de tuples (protein_name, "Cluster_N") — telle que
                      produite par run_anova_heatmaps.
    params, go_params, out_dir : identiques à run_go_enrichment.

    Retourne {cluster_id: {"table": df, "lollipop": path, "dotplot": path}}.
    Entièrement protégé : n'interrompt jamais le pipeline.

    Note z-score : les protéines d'un cluster proviennent de contrastes variés,
    il n'existe pas de LFC unique par protéine au niveau cluster. Le z_score des
    termes est donc neutralisé (0) ici — l'information portée par le cluster est
    le profil d'expression, pas un sens de régulation par contraste.
    """
    if go_params is None:
        return None

    try:
        organism = go_params.get("organism")
        use_string = go_params.get("backend", "gprofiler") == "string"
        species_taxid = go_params.get("species_taxid")
        string_caller = go_params.get("caller_identity") or DEFAULT_STRING_CALLER
        ortholog_map = go_params.get("ortholog_map")   # RBH (Perseverance ou STRING)
        # Regrouper les protéines par cluster
        from collections import defaultdict
        clusters = defaultdict(list)
        for prot, cid in cluster_mapping:
            if cid:  # ignorer les non-assignés ("")
                clusters[cid].append(str(prot).split(";")[0].strip())

        if not clusters:
            return None

        if use_string:
            print(f"\n[GO] Cluster GO enrichment — STRING taxid {species_taxid} "
                  f"({len(clusters)} clusters)")
            if ortholog_map:
                string_maps = string_maps_from_orthologs(ortholog_map)
                string_bg = (sorted(set(ortholog_map.values()))
                             if go_params.get("string_background") == "quantified" else None)
            else:
                all_prots = [p for plist in clusters.values() for p in plist]
                string_maps, string_bg = _prepare_string(
                    all_prots, species_taxid, string_caller,
                    go_params.get("string_background", "genome"))
            if string_maps is None:
                return None
        else:
            if not organism:
                return None
            print(f"\n[GO] Cluster GO enrichment — species '{organism}' "
                  f"({len(clusters)} clusters)")

        results = {}
        for cid in sorted(clusters.keys()):
            prots = list(dict.fromkeys(clusters[cid]))  # dédup, ordre conservé
            try:
                if len(prots) <= 5:
                    print(f"  [SKIP] {cid}: too few proteins ({len(prots)}), skipped.")
                    continue

                if use_string:
                    gost_df = _string_enrich_query(
                        prots, string_maps, int(species_taxid), string_caller,
                        categories=go_params.get("string_categories"),
                        background=string_bg, label=cid)
                elif ortholog_map:
                    # Perseverance -> gProfiler : même traduction/retraduction
                    # que run_go_enrichment (avant : IDs non traduits).
                    rev = {v: k for k, v in ortholog_map.items()}
                    tr = [ortholog_map[p] for p in prots if p in ortholog_map]
                    if len(tr) <= 5:
                        print(f"  [SKIP] {cid}: too few orthologs ({len(tr)}), skipped.")
                        continue
                    print(f"  [GO] {cid}: {len(tr)} orthologue(s) -> gProfiler ({organism})...")
                    gost_df = run_gost(tr, organism)
                    if gost_df is not None and not gost_df.empty and "intersections" in gost_df.columns:
                        gost_df = gost_df.copy()
                        gost_df["intersections"] = gost_df["intersections"].apply(
                            lambda inter: [rev.get(x, x) for x in inter]
                            if isinstance(inter, (list, tuple, np.ndarray)) else inter)
                else:
                    print(f"  [GO] {cid}: {len(prots)} proteins -> gProfiler...")
                    gost_df = run_gost(prots, organism)
                if gost_df is None or gost_df.empty:
                    print(f"     No significant enrichment.")
                    continue

                # df_local vide → process_enrichment neutralise le z-score (0.0)
                tab = process_enrichment(gost_df, pd.DataFrame())
                if tab.empty:
                    print(f"     No term after filtering (term_size, roots).")
                    continue

                def _safe_plot(fn, *a):
                    try:
                        return fn(*a)
                    except Exception as ex:
                        print(f"        (plot {fn.__name__} skipped: "
                              f"{type(ex).__name__})")
                        return None
                # On réutilise les plots existants ; l'étiquette = cid
                f_lol = _safe_plot(plot_lollipop_faceted, tab, cid, out_dir)
                f_dot = _safe_plot(plot_dotplot_faceted, tab, cid, out_dir)

                tab_export = tab.copy()
                if "intersection" in tab_export.columns:
                    tab_export["intersection"] = tab_export["intersection"].apply(
                        lambda x: ",".join(map(str, x))
                        if isinstance(x, (list, tuple, np.ndarray)) else str(x))
                tab_export = tab_export.drop(columns=["term_short", "minus_log10_p"],
                                             errors="ignore")

                results[cid] = {"table": tab_export,
                                "lollipop": f_lol, "dotplot": f_dot}
                n_bp = (tab["source"] == "GO:BP").sum()
                print(f"     [OK] {len(tab)} enriched terms ({n_bp} GO:BP).")

            except Exception as e:
                detail = f"missing column/key: {e}" if isinstance(e, KeyError) else str(e)[:80]
                print(f"     [WARN] Error on {cid}: {type(e).__name__} ({detail}) — skipped.")
                continue

        if not results:
            print("  [INFO] No significant GO enrichment across all clusters.")
            return None
        return results

    except Exception as e:
        print(f"  [WARN] Cluster GO enrichment interrupted ({type(e).__name__}) — "
              f"the main pipeline is not affected.")
        return None


def export_go_cluster_sheets(wb, go_cluster_results: dict, write_df_fn, ins_img_fn):
    """Ajoute un onglet Cluster_Enrich_N par cluster au classeur Excel.

    IMPORTANT — nommage : on n'utilise PAS le préfixe 'GO_' ici. Le dashboard
    (build_dashboardv7.py) capte toute feuille dont le nom contient 'go' pour
    construire sa heatmap GO cross-contrastes, puis apparie de force chaque
    feuille à un contraste (sans seuil de rejet). Des onglets de CLUSTER (qui ne
    correspondent à aucun contraste) y seraient appariés à tort et casseraient la
    heatmap. Le préfixe 'Cluster_Enrich_' échappe à tous les filtres du dashboard
    ('go' in name / startswith('GO_')), donc ces onglets restent visibles dans
    l'Excel sans perturber le dashboard.
    """
    if not go_cluster_results:
        return
    for cid, data in go_cluster_results.items():
        # cid = "Cluster_1" → onglet "Cluster_Enrich_1"
        n = cid.replace("Cluster_", "")
        sheet_name = f"Cluster_Enrich_{n}"[:31]
        ws = wb.create_sheet(sheet_name)
        df = data["table"]
        start_col = 21
        for j, col in enumerate(df.columns):
            ws.cell(row=1, column=start_col + j, value=str(col))
        for i, row in enumerate(df.itertuples(index=False), start=2):
            for j, val in enumerate(row):
                if isinstance(val, float) and np.isnan(val):
                    val = None
                elif isinstance(val, (list, tuple, np.ndarray)):
                    val = ",".join(map(str, val))
                ws.cell(row=i, column=start_col + j, value=val)
        ins_img_fn(ws, data.get("lollipop"), "A2")
        ins_img_fn(ws, data.get("dotplot"),  "A60")


# ==============================================================================
# 6. EXPORT EXCEL (onglets GO ajoutés au classeur principal)
# ==============================================================================

def _make_go_sheet_names(contrasts):
    """Construit des noms d'onglets GO courts, lisibles et GARANTIS uniques
    (<= 31 car., contrainte Excel).

    Problème résolu : avec de longs préfixes communs aux conditions
    (ex: 'X18_009IG01Ctrl_vs_X18_009IG01di6h'), une troncature naïve à 31 car.
    garde le préfixe commun et produit des noms IDENTIQUES pour deux contrastes
    différents -> Excel corrompt /xl/workbook.xml ('Réparations...').

    Stratégie : retirer le préfixe/suffixe commun aux conditions (partie
    discriminante), recomposer 'GO_<a>v<b>', tronquer, et dédupliquer avec un
    suffixe numérique si nécessaire.

    Retourne {contrast: sheet_name}.
    """
    # Conditions = les 2 côtés de chaque contraste
    all_conds = []
    for c in contrasts:
        all_conds.extend(c.split("_vs_"))
    cond_label = _strip_common_affix(all_conds)

    names = {}
    used = set()
    for c in contrasts:
        parts = c.split("_vs_")
        a = cond_label.get(parts[0], parts[0])[:12] if parts else c[:12]
        b = cond_label.get(parts[1], parts[1])[:12] if len(parts) > 1 else ""
        base = f"GO_{a}v{b}" if b else f"GO_{a}"
        base = base[:31]
        name = base
        i = 2
        # Dédup anti-collision (réserve 3 car. pour le suffixe #NN)
        while name in used:
            suffix = f"#{i}"
            name = (base[:31 - len(suffix)]) + suffix
            i += 1
        used.add(name)
        names[c] = name
    return names


def export_go_sheets(wb, go_results: dict, write_df_fn, ins_img_fn):
    """
    Ajoute un onglet par contraste au classeur Excel existant.
    Mise en page proche du script R : données en colonne U (21),
    plots empilés à gauche (Manhattan, lollipop, dotplot, chord).
    """
    if not go_results:
        return

    sheet_names = _make_go_sheet_names(list(go_results.keys()))

    for contrast, data in go_results.items():
        sheet_name = sheet_names[contrast]
        ws = wb.create_sheet(sheet_name)

        # Données décalées en colonne 21 (U) comme le script R
        df = data["table"]
        start_col = 21
        for j, col in enumerate(df.columns):
            ws.cell(row=1, column=start_col + j, value=str(col))
        for i, row in enumerate(df.itertuples(index=False), start=2):
            for j, val in enumerate(row):
                if isinstance(val, float) and np.isnan(val):
                    val = None
                elif isinstance(val, (list, tuple, np.ndarray)):
                    val = ",".join(map(str, val))
                ws.cell(row=i, column=start_col + j, value=val)

        # Plots empilés à gauche : lollipop (z-score) puis dotplot (gene ratio)
        # Plots empilés à gauche : lollipop (z-score) puis dotplot (gene ratio).
        # Layout vertical = images plus hautes → espacement vertical large.
        ins_img_fn(ws, data.get("lollipop"),  "A2")
        ins_img_fn(ws, data.get("dotplot"),   "A60")
