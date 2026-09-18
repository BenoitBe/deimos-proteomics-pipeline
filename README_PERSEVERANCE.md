# Perseverance — Transfert d'annotation fonctionnelle par orthologie

*Nommé d'après le rover Mars 2020, mandaté pour chercher des signes de vie
ancienne : déduire quelque chose d'invisible en surface à partir d'une
source bien caractérisée — exactement la logique du transfert par
orthologie.*

Cinquième outil de la suite, **indépendant** de Deimos/Phobos/Pathfinder/Odyssey.
Répond à un problème précis : pour un organisme non-modèle (ex: *Brassica
napus*), la couverture GO/KEGG de gProfiler est souvent famélique, ce qui
biaise silencieusement l'enrichissement fonctionnel (`go_enrichment.py`,
Pathfinder Meta-GO) vers le sous-ensemble de protéines déjà bien annotées.

Perseverance trouve les orthologues dans un protéome de référence bien annoté
(ex: *Arabidopsis thaliana* pour les Brassicacées) par **Reciprocal Best
Hit** (RBH, via DIAMOND) et transfère GO/KEGG/Pfam avec un **niveau de
confiance explicite** — jamais mélangé silencieusement avec une annotation
directe.

## Pourquoi RBH plutôt qu'un simple best-hit ?

Un best-hit unidirectionnel confond facilement orthologues et paralogues
proches : si votre espèce a deux gènes dupliqués récemment (fréquent chez
les polyploïdes comme *B. napus*), le meilleur hit d'un des deux peut être
le MAUVAIS paralogue de référence. RBH exige la réciprocité (A→B ET B→A)
ce qui résout correctement ce cas — validé explicitement dans les données
de test incluses (`toy_data/`, voir plus bas).

## Installation

```bash
sudo apt install diamond-aligner   # ou: conda install -c bioconda diamond
pip install pandas
```

## Usage

```bash
python perseverance.py \
    --query monorganisme.fasta \
    --reference arabidopsis.fasta \
    --ref-annotations arabidopsis_go.tsv \
    --out annotation_transferee.tsv \
    --min-identity 25 --min-coverage 50
```

`arabidopsis_go.tsv` : table `protein_id  GO  KEGG  Pfam` (GO en liste
`;`-séparée) — récupérable depuis UniProt (colonne GO) ou TAIR pour
*Arabidopsis*, ou toute base d'annotation de l'espèce de référence choisie.

Sortie (`annotation_transferee.tsv`) : `protein_id, GO, KEGG, Pfam, source
(direct|orthologue RBH), confidence (high/medium/low), ortholog_id, pident,
coverage` — directement utilisable comme annotation custom dans topGO,
clusterProfiler, ou `go_enrichment.py`.

### Combiner avec une annotation directe partielle

Si gProfiler couvre déjà une partie de votre espèce, passez
`--direct-annotations partiel.tsv` (mêmes colonnes) : ces annotations ont
priorité sur le transfert par orthologie et sont marquées `source=direct`.

## Choisir le seuil de confiance en aval

- **high** (≥50% identité, ≥70% couverture) : fiable pour de l'enrichissement
  GO standard.
- **medium** (≥30%/≥50%) : utilisable, mais à isoler dans une analyse de
  sensibilité (refaire l'enrichissement avec/sans, comparer).
- **low** : à ne garder que pour une exploration qualitative, pas pour un
  test statistique d'enrichissement — l'annotation fonctionnelle peut avoir
  divergé même quand l'homologie de séquence est détectable.

## Données de test incluses (`toy_data/`)

Séquences synthétiques construites pour valider spécifiquement le cas
piège des paralogues proches : `AT2G_KINASE1` et `AT4G_KINASE2` (85%
identiques entre eux) dans la référence, deux orthologues query distincts
correspondant à chacun, et deux protéines query "nouvelles" sans homologue
(doivent rester sans annotation). Lancer :

```bash
python perseverance.py --query toy_data/query.fasta --reference toy_data/reference.fasta \
    --ref-annotations toy_data/reference_annotations.tsv --out test_out.tsv --threads 2
```

Résultat attendu : 4/6 protéines annotées (confiance `high` partout), les
deux paralogues correctement dissociés, les 2 protéines nouvelles sans
annotation.

## Limites connues

1. **Aucun protéome de référence n'est fourni** — RBH nécessite un FASTA +
   une table d'annotation d'une espèce bien caractérisée, à récupérer
   vous-même (UniProt, TAIR, Ensembl Plants...). Le réseau de développement
   de cet outil n'a pas accès à ces bases pour un test à grande échelle ;
   seules les données synthétiques `toy_data/` ont été testées de bout en
   bout ici.
2. **RBH suppose une relation 1:1**. Pour des familles multigéniques
   complexes (WGD récent, très fréquent chez *B. napus* allotétraploïde),
   plusieurs vrais orthologues fonctionnels peuvent exister par gène de
   référence — RBH n'en retient qu'un. Une extension vers un orthogroupe
   complet (OrthoFinder, eggNOG) serait nécessaire pour ce cas, au prix
   d'une complexité et d'un temps de calcul largement supérieurs.
3. **Le choix de l'espèce de référence est manuel** — pas de sélection
   automatique de la référence la plus proche disponible.

## Architecture

```
perseverance.py           CLI (build + résumé de couverture)
rbh.py                DIAMOND makedb + blastp reciproque + calcul RBH
transfer.py            Transfert d'annotation + confiance + résumé
toy_data/               Données synthétiques de validation (avec piège paralogue)
```
