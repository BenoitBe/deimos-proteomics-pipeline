# Intégration Perseverance <-> Deimos

Perseverance devient une source d'annotation GO **branchée directement dans l'étape
d'enrichissement existante de Deimos**, en remplacement de gProfiler quand
l'espèce est mal couverte — pas un outil séparé à recoller manuellement.

## Principe

`run_local_ora()` (nouveau, dans `go_enrichment.py`) produit exactement le
même contrat de sortie que `run_gost()` (colonnes gProfiler natives). Résultat :
**`process_enrichment()`, les graphiques (lollipop, dotplot) et l'export Excel
fonctionnent sans aucune modification**, qu'ils reçoivent du gProfiler ou du
Perseverance. Testé de bout en bout et confirmé sur les données réelles.

## Workflow complet

```
1. FASTA1 (recherche DIA-NN, déjà utilisé par Deimos)
   + liste des Protein.Group réellement quantifiés (Log2_Impute)
        │
        ▼  deimos_to_perseverance.export_query_fasta()
   perseverance_query.fasta  +  mapping.tsv (accession -> Protein.Group)
        │
        ▼  perseverance.py --query perseverance_query.fasta --reference FASTA2 (ex: Arabidopsis)
   annotation_transferee.tsv (indexé par accession individuelle)
        │
        ▼  deimos_to_perseverance.remap_annotation_to_protein_groups()
   annotation_protein_group.tsv (indexé Protein.Group, convention Deimos)
        │
        ▼  config.yaml: perseverance_annotation_path: "annotation_protein_group.tsv"
   Deimos (go_enrichment.py) : run_local_ora() remplace run_gost()
        │
        ▼
   Mêmes sheets Excel GO, mêmes graphiques, mêmes seuils qu'avec gProfiler
```

## Utilisation

```bash
# 1. Extraire le FASTA query depuis Deimos
python -c "
import pandas as pd
from deimos_to_perseverance import export_query_fasta
mat = pd.read_excel('proteomics_output/deimos_results.xlsx', sheet_name='Log2_Impute')
export_query_fasta('search_fasta1.fasta', mat['Protein.Group'].tolist(),
                   'perseverance_query.fasta', 'mapping.tsv')
"

# 2. Perseverance
python perseverance.py --query perseverance_query.fasta --reference arabidopsis.fasta \
    --ref-annotations arabidopsis_go.tsv --out perseverance_annotation_raw.tsv

# 3. Remap vers la convention Deimos
python -c "
import pandas as pd
from deimos_to_perseverance import remap_annotation_to_protein_groups
annot = pd.read_csv('perseverance_annotation_raw.tsv', sep='\t')
remap_annotation_to_protein_groups(annot, 'mapping.tsv').to_csv(
    'annotation_protein_group.tsv', sep='\t', index=False)
"

# 4. Dans config.yaml (ou en CLI si tu ajoutes l'option) :
#    perseverance_annotation_path: "annotation_protein_group.tsv"
python deimos.py --config config.yaml
```

## Ce que le test de bout en bout a validé (données réelles)

- Extraction FASTA : groupes composites (`"A;B"`) correctement éclatés puis
  refusionnés après annotation.
- `run_local_ora()` détecte un enrichissement significatif (p=0.022) sur les
  deux orthologues kinases du jeu de test, avec fond non-dégénéré.
- `process_enrichment()` (fonction Deimos **non modifiée**) consomme cette
  sortie sans erreur : renommage de colonnes, `gene_ratio`, `z_score`
  (moyenne exacte des LFC), filtre `term_size≥5` — tout fonctionne comme
  avec une vraie réponse gProfiler.

## Limites connues

1. **Background statistique** : par défaut, dérivé de `df_comparaison`
   (toutes les protéines du run) si `go_params["background"]` n'est pas
   fourni explicitement — cohérent avec la pratique recommandée en ORA
   protéomique (univers = protéines quantifiables, pas génome entier).
2. **`term_name` = `term_id`** : Perseverance ne résout pas les noms lisibles des
   termes GO (pas de parsing `.obo` intégré) — les sheets Excel afficheront
   des codes `GO:0004672` plutôt que "protein kinase activity". Un module
   de résolution de noms (via `goatools` + `go-basic.obo`) serait le
   complément naturel si ce point gêne en pratique.
3. **Namespace GO (BP/CC/MF) non distingué** : tous les termes Perseverance
   partagent la valeur `source="GO:PERSEVERANCE"`, contrairement à gProfiler qui
   sépare BP/CC/MF/REAC/KEGG. Les graphiques facettés par source
   fonctionnent toujours, mais avec une seule facette au lieu de plusieurs.
4. **CLI non étendu** : le paramètre `perseverance_annotation_path` est
   utilisable via `config.yaml` mais n'a pas d'option `--perseverance-annotation`
   dédiée dans `config.py` (contrairement à `--batch-column`) — à ajouter
   si un usage en ligne de commande directe est souhaité.
