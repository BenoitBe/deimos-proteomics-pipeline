# Lire les résultats d'un run Deimos

Guide des trois sorties d'un run : le classeur Excel, le dashboard Deimos et le
rapport Pathfinder. Écrit pour qu'on puisse ouvrir les fichiers sans se demander
par où commencer.

```
proteomics_output/        (ou <project_dir>/output/ si vous travaillez par projet)
├── ProteomicAnalysis_Results.xlsx   les chiffres, c'est la source
├── Proteomic_dashboard.html         explorer, filtrer, exporter  (nom réglable)
├── pathfinder_report.html           la lecture d'ensemble
├── pathfinder/                      le matériel de travail de Pathfinder :
│                                    synthesis.xlsx, summary.md, log, design
├── Volcanos/                        volcanos individuels + facet, pour les figures
├── matrices/                        tables brutes pour R / Perseus
└── last_config.yaml                 ce qui a été demandé — à joindre à toute question
```

**Par où commencer** : `pathfinder_summary.md`. Il tient en une page et dit, en
phrases, ce que le run a trouvé. Ensuite le dashboard pour creuser, le classeur
pour les chiffres exacts.

---

## 1. Les trois questions à se poser avant de lire quoi que ce soit

**Quelle règle a défini « différentiel » ?** Tout en dépend. Elle est écrite dans
la feuille `Methods`, reprise en tête du dashboard et de Pathfinder. Deux runs
avec des seuils différents ne se comparent pas.

**p brute ou p ajustée ?** `volcano_use_padj: true` = p.adj (BH). C'est plus
sévère et c'est ce qu'on publie. En p brute, attends-toi à beaucoup de faux
positifs : sur 10 000 protéines, 500 passent p < 0,05 par hasard seul.

**Quelle statistique principale ?** Avec `primary_statistic: dpc`, les colonnes
sans suffixe sont le DPC sans imputation, et `*_imputed` l'analyse imputée en
second. C'est l'inverse de ce que font la plupart des pipelines : les colonnes
principales ne dépendent d'aucune valeur inventée.

---

## 2. Le classeur `ProteomicAnalysis_Results.xlsx`

C'est la source. Tout ce que montrent les deux HTML en vient.

### Les feuilles qu'on ouvre vraiment

| Feuille | Ce qu'on y trouve |
|---|---|
| **Differential_Expression** | la table centrale : une ligne par protéine, les colonnes de chaque contraste |
| **ANOVA_Results** | F-test sur toutes les conditions à la fois, avec `Cluster_ID` |
| **ANOVA_\<strate\>** | une par strate si `anova_strata` est réglé (une ANOVA par temps, par exemple) |
| **QC** | taux de manquants, CV, corrélations entre réplicats |
| **Methods** | la règle du run, à copier dans un article |
| **_manifest** | ce que le dashboard relit ; sert à diagnostiquer, pas à lire |

### Lire `Differential_Expression`

Pour chaque contraste `A_vs_B`, quatre colonnes :

| Colonne | Sens |
|---|---|
| `*_diff` | log2FC. Positif = plus abondant dans **A**, le numérateur |
| `*_p.val` | p brute |
| `*_p.adj` | p ajustée BH, dans la famille de ce contraste seul |
| `Pi_Score_*` | \|log2FC\| × −log10(p) : sert à **classer**, jamais à sélectionner |

Et par protéine : `imputed` (VRAI si au moins une valeur a été imputée),
`num_NAs` (combien), `n_precursors` (sur combien de précurseurs repose le DPC),
`Robustness_Score_*` (stabilité en laissant un échantillon de côté).

**Le piège classique** : trier sur `*_diff` seul. Un log2FC de 5 sur une protéine
vue dans deux échantillons sur huit ne vaut pas un log2FC de 1 mesuré partout.
Croise toujours avec `num_NAs` et la p.

### Les feuilles de sensibilité

- **Sensitivity_DPC** : pour chaque contraste, le log2FC du DPC, celui de
  l'analyse imputée, et leur `concordance`. Une concordance basse signale une
  protéine dont la conclusion dépend de la méthode — à regarder avant de bâtir
  une histoire dessus.
- **Differential_Detection** : présence/absence plutôt qu'abondance. `pattern`
  dit la forme (absente d'une condition entière, sporadique…), `higher_in` où
  elle est vue. Utile quand une protéine est tout ou rien.

### Les feuilles de qualité plasma

À lire **avant** `Differential_Expression` sur tout projet de plasma ou de sérum.
Elles existent si l'étape a tourné : `plasma_qc: true`, `--plasma-qc`, ou réponse « oui » à
la question que Deimos pose quand les données ressemblent à du plasma (`__auto__`, défaut).
Un échantillon `flagged` est à regarder, pas forcément à exclure : sur du plasma propre
simulé, la règle de l'article signale 5 % des échantillons sur au moins un indice.

- **Plasma_QC** : un indice par échantillon pour l'hémolyse, la contamination
  plaquettaire et la coagulation (panels de Geyer et al. 2019) ; plus haut = plus
  contaminé. `flagged` signale un échantillon hors norme, il ne le retire pas :
  l'exclure ou non est une décision à prendre, puis relancer. Le tableau *Index vs
  condition* est le plus important : un indice qui diffère entre conditions ne
  signale aucun échantillon, mais fabrique des DEP — c'est le cas le plus dangereux,
  parce qu'il ne se voit pas ailleurs.
- **Plasma_QC_DEP** : une ligne par DEP et par contraste. `high` = la protéine
  appartient à un panel ou suit un indice, **et** cet indice diffère entre les deux
  groupes dans le même sens : à considérer comme un artefact tant que rien ne prouve
  le contraire. `check` = membre d'un panel ou corrélée à un indice, sans différence
  d'indice entre les groupes : candidat possible, à valider. Une protéine hors panel
  peut être `high` : elle a été repérée par sa co-variation avec l'indice d'un
  échantillon à l'autre au sein des conditions.

---

## 3. Le dashboard `Proteomic_dashboard.html`

Un seul fichier, aucun serveur, aucune connexion : il s'ouvre par double-clic et
fonctionne hors ligne. Toutes les figures s'exportent en PNG ou SVG, toutes les
tables en Excel. Le bouton ☀️ / 🌙 bascule jour/nuit ; **le jour est le thème par
défaut**, et ton choix est retenu pour les visites suivantes.

| Page | À quoi elle sert |
|---|---|
| **Overview** | combien de DEP par contraste, la règle du run, l'état du QC |
| **Volcano** | le contraste en détail. `⇄ Inverser` permute numérateur et dénominateur, export iPath3 |
| **FC vs FC** | deux contrastes l'un contre l'autre. La pente est calculée sur les **DEP communes seulement** |
| **GO** | enrichissement : bulles, table, réseau de termes, heatmap entre contrastes |
| **GSEA** | NES par contraste, leading edge, matrice entre contrastes |
| **Heatmap** | z-scores, clusters ANOVA, curseurs de filtrage |
| **Intersections** | UpSet : quelles protéines sont partagées par quels contrastes |
| **DEP Table** | la table filtrable, avec recherche globale |
| **QC / PCA** | manquants, CV, réplicats, ACP et ses protéines motrices, UMAP |
| **WGCNA** | modules de co-expression, corrélation module-condition, hubs |

**Ce qu'il faut regarder en premier dans QC / PCA**, avant toute interprétation
biologique :

- le **raincloud observé / imputé** : si les valeurs imputées forment un second
  mode séparé, l'imputation invente une population ;
- les **motifs de manquants** (UpSet) : une protéine absente d'une condition
  entière n'est pas la même chose qu'une protéine dispersée ;
- la **corrélation entre réplicats** : un échantillon qui décroche se voit ici,
  pas dans les volcanos ;
- l'**ACP** : si les conditions ne se séparent pas sur PC1-PC2, aucune liste de
  DEP ne sauvera le run.

---

## 4. Le rapport `pathfinder_report.html`

Là où Deimos liste, Pathfinder **rapproche** : il croise les contrastes, regroupe
les termes redondants et écrit des phrases. Lui aussi s'ouvre hors ligne, en
thème jour par défaut.

| Onglet | À quoi il sert |
|---|---|
| **Summary** | le run en phrases. Le même texte que `pathfinder_summary.md` |
| **Biological reading** | programmes de réponse, groupes fonctionnels, plans factoriels |
| **Overview** | ACP, facteurs du plan, qualité, homogénéité |
| **Proteins** | signature : récurrentes, spécifiques, inversées entre contrastes |
| **Meta-GO** | termes enrichis regroupés par redondance, direction par contraste |
| **Exploration** | détail par contraste |
| **Protein UMAP** | les protéines, pas les échantillons |

### Les figures de *Biological reading*

- **Dynamic range** — les protéines rangées par abondance, DEP marquées, et la
  part de DEP par décile. Elle répond à « mes DEP sont-elles juste les plus
  abondantes ? ». Si le premier décile en contient beaucoup plus que le dernier,
  une partie du signal est de la détectabilité, pas de la biologie.
- **Enriched terms and their redundancy** — treemap (surface = −log10 q) et arbre
  de redondance. **Un gros bloc dans le treemap n'est pas une preuve forte** :
  c'est un paquet de termes qui touchent les mêmes protéines, donc une seule
  observation, pas vingt. L'arbre le montre directement.
- **Shape of each programme** — radar des programmes, de 3 à 8 conditions.
- **Order and end-to-end change** — slope graph entre les niveaux extrêmes du
  facteur ordonné, et bump chart du rang des programmes. Des lignes qui se
  croisent veulent dire que les programmes changent d'**ordre**, ce qu'une figure
  d'amplitude cache.
- **Variance shares per protein** — violons des parts de variance, DEP à droite,
  **fond de toutes les protéines quantifiées à gauche**. Les DEP sont choisies
  sur des différences entre conditions : leurs parts « facteur » sont gonflées
  par construction. **Un facteur se lit par son écart au fond, pas en absolu.**
  Un facteur qui passe sous le fond n'est pas un effet dilué, c'est une absence
  d'effet.
- **Shared and divergent responses** — matrice de toutes les paires de
  contrastes, puis nuages log2FC contre log2FC. La pente est ajustée sur les
  **DEP communes seulement** : sur l'union, les protéines significatives dans un
  seul contraste forment une croix le long des axes et tirent la pente.
- **Programmes and functional groups** — graphe biparti. Un groupe relié à
  plusieurs programmes est un groupe que le clustering ne sépare pas.

Si un bandeau rouge remplace cet onglet, la lecture biologique a échoué : le
bandeau nomme l'erreur, le reste du rapport reste valide.

### Le classeur `pathfinder_synthesis.xlsx`

Tout ce que les figures montrent, en chiffres : `Programmes`,
`Variance_partition`, `Variance_shares`, `FC_vs_FC`, `MetaGO_matrix`,
`Programme_groups`, `PCA_*`… La feuille **`Reading_columns`** décrit les autres,
colonne par colonne : c'est le point d'entrée quand on cherche où est un chiffre.

---

## 5. Erreurs de lecture les plus fréquentes

**Prendre le Pi-Score pour un test.** Il classe, il ne sélectionne pas.
L'inclusion repose sur la p ajustée.

**Lire un enrichissement GO sans regarder le fond.** Avec un fond « génome »,
tout compartiment abondant ressort, parce que la protéomique voit mieux les
protéines abondantes. Le fond « quantified » corrige ce biais. Il est écrit dans
l'en-tête de la page GO.

**Comparer des log2FC entre contrastes de tailles différentes.** Un contraste à
trois réplicats et un à sept n'ont pas la même puissance : un effet absent du
premier peut y être tout autant, sans atteindre le seuil.

**Conclure sur une protéine vue dans un seul contraste sur beaucoup.** Avec un BH
par contraste, chaque contraste est sa propre famille de tests : sur vingt
contrastes, les faux positifs s'accumulent. Une protéine isolée demande une
confirmation.

**Oublier `num_NAs`.** Un log2FC spectaculaire porté par deux valeurs mesurées et
six imputées n'est pas un résultat, c'est une question.

---

## 6. Quand quelque chose manque

| Absent | Pourquoi, le plus souvent |
|---|---|
| Page GO vide | `go_organism` non renseigné, ou aucun terme significatif |
| Page GSEA absente | `run_gsea: false`, ou pas de jeux de gènes pour l'espèce |
| Onglet *Biological reading* absent | la lecture a échoué — voir `pathfinder/pathfinder_log.txt` |
| Sections factorielles absentes | le plan n'a pas été reconnu comme factoriel ; essayer `pathfinder_options: {factor_regex: ...}` |
| Radar absent | plus de 8 conditions, c'est voulu |
| Feuilles ANOVA par strate absentes | `anova_strata` non réglé |
| WGCNA absent | `make_wgcna: false`, ou trop peu d'échantillons |

Pour toute question sur un run, joins `last_config.yaml` : il dit exactement ce
qui a été demandé.
