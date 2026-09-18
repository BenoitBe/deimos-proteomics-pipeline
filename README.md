<p align="center">
  <img src="deimos_icon.svg" width="90" alt="Deimos">
</p>

<h1 align="center">Deimos</h1>
<p align="center"><i>DIA Expression Integrated Multi-Omics Suite</i></p>
<p align="center">Proteogen · Université de Caen</p>

---

Pipeline d'analyse protéomique quantitative label-free en DIA (DIA-NN) :
filtrage, imputation, tests différentiels (limma/DEqMS), GO/KEGG, WGCNA, et
un dashboard HTML interactif — le tout piloté par un seul fichier de
configuration, sans dépendance à R.

## Sommaire

- [Installation](#installation)
- [Démarrage rapide](#démarrage-rapide)
- [Fonctionnalités](#fonctionnalités)
- [Configuration](#configuration)
- [Structure du projet](#structure-du-projet)
- [Écosystème](#écosystème)
- [Limites connues](#limites-connues)

## Installation

```bash
pip install -r requirements.txt
```

Le moteur statistique (`limma_ebayes.py`) est une réimplémentation Python de
`limma::eBayes`/DEqMS — **aucune dépendance à R**. Validé numériquement
contre le vrai `limma`/`DEqMS` R (précision machine sur eBayes ; voir
historique du projet pour le détail de cette validation).

Dépendances optionnelles selon les modules activés :
- `gprofiler-official` — enrichissement GO/KEGG en ligne (`go_organism`)
- `diamond-aligner` (`sudo apt install diamond-aligner` ou `conda install
  -c bioconda diamond`) — uniquement si vous utilisez **Perseverance**
  (recherche d'orthologues, voir [Écosystème](#écosystème))
- `chart.umd.js` / `xlsx.full.min.js` — à placer à côté de
  `build_dashboardv7.py` pour un dashboard 100% hors-ligne ; sinon repli
  automatique sur les CDN (cdnjs.cloudflare.com)

## Démarrage rapide

```bash
python deimos.py --config config_example.yaml
```

ou en mode interactif (pose les questions au fur et à mesure) :

```bash
python deimos.py --tsv report.pg_matrix.tsv --design ExperimentalDesign.csv
```

Entrées attendues :
- `report.pg_matrix.tsv` — export protéines DIA-NN
- `ExperimentalDesign.csv` — colonnes `label;condition;replicate[;subject]`
- `report.pr_matrix.tsv` — optionnel, requis pour DEqMS et Perseverance

Sortie : `proteomics_output/ProteomicAnalysis_Results.xlsx` (classeur complet)
+ `proteogen_dashboard.html` (dashboard interactif autonome).

## Fonctionnalités

### Cœur statistique
- Filtrage, imputation QRILC/Mixed avec **diagnostic automatique de
  missingness** (MNAR vs MAR) qui objective le choix de méthode plutôt que
  de le figer arbitrairement.
- Tests différentiels par `limma` eBayes (réimplémentation Python validée),
  avec **DEqMS optionnel** — activé seulement si un
  [diagnostic de fiabilité](#deqms) le juge exploitable sur le run en cours.
- ANOVA multi-groupes, heatmaps clusterisées, scores de robustesse
  (rééchantillonnage), Pi-score.

### Design expérimental
- **Détection automatique de la condition contrôle** (`ctl`, `ctrl`,
  `control`, `wt`, `mock`...) systématiquement placée en dénominateur des
  contrastes (`exp_vs_ctl`, jamais l'inverse) — surchargeable explicitement
  (`control_condition: "CS8"`) si le nom ne matche aucun synonyme standard.
- Design **apparié** (`subject`) pour les mesures répétées.
- **Correction de batch optionnelle** (`batch_column`) : contrôlée comme
  covariable dans le modèle statistique (pas une correction naïve appliquée
  à la matrice puis retestée — voir la mise en garde `removeBatchEffect`
  dans le code), avec diagnostic PCA avant/après (Kruskal-Wallis + ε²) pour
  objectiver le gain.
- **Détection multi-espèces** automatique (suffixe `_TAG` façon UniProt dans
  `Protein.Names`, ex. hôte + pathogène) avec analyse séparée par organisme
  sur demande.

### Enrichissement fonctionnel (GO/KEGG)
- gProfiler direct (`go_organism`), avec un **diagnostic de couverture**
  (g:Convert) qui mesure objectivement le taux de reconnaissance de vos
  accessions pour l'organisme choisi *avant* de lancer l'enrichissement.
- Si la couverture est insuffisante (espèce non-modèle mal indexée) :
  intégration avec **Perseverance**, qui traduit vos protéines
  significatives vers les orthologues d'une espèce de référence bien
  annotée par Reciprocal Best Hit (DIAMOND), interroge le **vrai gProfiler**
  sur cette espèce, puis retraduit les résultats vers vos identifiants
  d'origine — automatisable en une question posée en cours de pipeline.
- GSEA en complément de l'ORA classique.

### Dashboard interactif
- Un seul fichier HTML autonome (~850 Ko sur le jeu d'exemple), Chart.js +
  canvas fait main pour les visualisations complexes (volcano, heatmap, PCA,
  UMAP, réseaux GO/WGCNA).
- **Bouton clair/sombre** (🌙/☀️), persistant, avec détection de la
  préférence système au premier chargement.
- Export PNG par figure, tables filtrables, recherche par gène/accession.

### WGCNA
- Modules de co-expression, corrélation module-trait, réseau hub-centrique
  exporté dans le dashboard.

## Configuration

Toutes les options se pilotent via un fichier YAML (voir
`config_example.yaml` pour la version commentée complète). Paramètres les
plus utiles au-delà des seuils statistiques standards :

| Clé | Défaut | Effet |
|---|---|---|
| `use_deqms` | `false` | Active DEqMS (sous réserve du diagnostic auto) |
| `deqms_force` | `false` | Force DEqMS même si le diagnostic est défavorable |
| `control_condition` | `"__auto__"` | `"__auto__"` (détection), un nom explicite, ou `false` (désactivé) |
| `batch_column` | `null` | Colonne de `ExperimentalDesign.csv` à contrôler comme covariable |
| `go_organism` | `null` | Code espèce gProfiler (ex: `"bnapus"`) |
| `run_perseverance` | `"__auto__"` | `"__auto__"` (demande en cours de run), `true`/`false` |
| `perseverance_reference_organism` | `null` | Code gProfiler de l'espèce de référence pour l'orthologie |
| `make_dashboard` | `true` | Génère le dashboard HTML interactif |
| `make_wgcna` | `false` | Active l'analyse WGCNA |

Le chargeur de config fusionne librement toute clé YAML supplémentaire —
aucune modification de `config.py` n'est nécessaire pour utiliser un
paramètre déjà lu par `deimos.py` via `params.get(...)`.

## Structure du projet

```
deimos.py                    Orchestrateur principal
limma_ebayes.py                Moteur statistique (eBayes/DEqMS, design matrix, contrastes)
config.py                      Chargement YAML/CLI, valeurs par défaut
go_enrichment.py                gProfiler direct + ORA local + diagnostic de couverture
gsea_enrichment.py              GSEA rank-based
dashboard_integration.py        Pont Deimos -> dashboard (adaptation des données)
build_dashboardv7.py            Génération du dashboard HTML
dashboard_template.html         Squelette HTML/CSS/JS du dashboard
deimos_to_perseverance.py       Pont Deimos <-> Perseverance (FASTA, orthologues)
config_example.yaml             Configuration commentée, point de départ recommandé
example_data/                   Jeu de données minimal pour tester le pipeline
```

Outils compagnons de l'écosystème, dépôts séparés :

## Écosystème

| Outil | Rôle |
|---|---|
| **Deimos** (ce dépôt) | Pipeline DIA quantitatif complet |
| [Phobos](https://github.com/BenoitBe/phobos-peptidomics-pipeline) | Pipeline DDA/PEAKS peptidomique + PTM |
| [Pathfinder](https://github.com/BenoitBe/Pathfinder) | Synthèse cross-contrastes de plusieurs runs Deimos/Phobos, sans recalcul |
| **Perseverance** | Recherche d'orthologues (RBH/DIAMOND) pour organismes non-modèles — voir `README_PERSEVERANCE.md`, intégré directement dans Deimos |
| **Odyssey** | Génération de dépôt ProteomeXchange/PRIDE (SDRF, checksums, submission.px) |

## Limites connues

- La correction de batch et l'annotation par orthologie (Perseverance)
  partagent actuellement un seul jeu de paramètres entre organismes en mode
  multi-espèces — pas encore de configuration par organisme dans ce cas.
- Le dashboard est un fichier HTML unique auto-généré par patches successifs
  (`build_dashboardv7.py`) plutôt que par un vrai moteur de templating — voir
  `AUDIT_DASHBOARD.md` pour le détail de cette dette d'architecture et les
  recommandations associées.
- `example_data/` doit rester synchronisé manuellement avec les colonnes
  attendues par `deimos.py` (ex: `N.Sequences`) — ces colonnes sont
  désormais optionnelles (remplies à `NaN` si absentes) plutôt que
  bloquantes.

## Licence

Voir `LICENSE`.
