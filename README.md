<p align="center">
  <img src="deimos_icon.svg" width="90" alt="Deimos">
</p>

# Deimos

**DIA Expression Integrated Multi-Omics Suite**

Complete LFQ DIA proteomics pipeline — from a DIA-NN matrix to an interactive dashboard.

Designed for LC-MS/MS DIA workflows on TIMS-TOF (Bruker). Compatible with any experimental design.

---

## Features

- **Quality control**: detection frequency, missing values, inter-sample correlation, RLE plot
  (diagnostic only — no re-normalisation after DIA-NN)
- **Filtering and imputation** with an automatic missingness diagnostic:
  - `impute_method: auto` (default) — MIXED if ≥ 20 % of missing values are MAR, QRILC otherwise;
    the choice and its justification are written to the Methods sheet
  - MIXED imputes **every** missing value of a protein absent from a whole condition by QRILC,
    and uses kNN only for proteins with partial missingness. The earlier rule (QRILC for the
    absent condition, kNN for the partial missingness of the other one) created artificial
    differences; it is kept as `impute_method: mixed_legacy` to reproduce older analyses
  - QRILC as in Lazar *et al.* (2016): mean and SD of the left-censored distribution estimated
    by quantile regression, not from the observed values alone
  - fixed random **seed** (`seed: 42` by default): two runs on the same input give the same
    p-values; the seed is recorded in Methods
- **Primary differential statistic without imputation** (`primary_statistic: dpc`, default):
  method of the **limpa** package (Li & Smyth), independent Python implementation checked
  against limpa 1.5.2 (same significant calls, 25x faster). A detection probability curve is
  fitted on the precursors (`report.pr_matrix.tsv`), proteins are quantified from their
  precursors with standard errors — an undetected precursor contributes its probability of
  non-detection, nothing is imputed — and limma is fitted with precision weights. On null
  comparisons of real data it is the only method that kept the false discovery rate at its
  nominal level (see *Validation* below). Volcanos, UpSet, GO, GSEA, the dashboard and
  Pathfinder use it; the imputed analysis is kept in `*_imputed` columns. Without a usable
  pr_matrix, Deimos falls back to the imputed analysis and says so in the console and Methods.
  `primary_statistic: imputation` restores the previous behaviour
- **Imputed analysis**: **limma** + **DEqMS** (peptide-count weighting), with an automatic
  reliability diagnostic (variance~peptide-count relationship) that disables DEqMS on runs
  where its premise doesn't hold, rather than silently returning misleading results
- **Experimental design handling**:
  - automatic detection of the **control condition** (`ctl`, `ctrl`, `control`, `wt`, `mock`...),
    systematically placed as the denominator (`treated_vs_control`, never the reverse) —
    overridable with an explicit name
  - **contrast selection** (`contrasts`): all pairs by default, or explicit patterns such as
    `Live_{s}_vs_Heat_killed_{s}` (same value on both sides) or `Live_*_vs_Live_*`
  - **paired designs** (`subject` column): model `~0 + condition + subject` when a subject spans
    ≥ 2 conditions
  - **pseudo-replication check**: warning (console + Methods) when a subject has several samples in
    the same condition — those samples are not independent replicates. Warning only: the model is
    not modified
  - optional **batch covariate** (`batch_column`), modelled in the design (not a matrix correction
    reused for testing), with a before/after PCA diagnostic (Kruskal-Wallis + ε²)
- **QC-only mode** (`qc_only` / `--qc-only`), enabled automatically when no condition has a
  replicate: raw data, Methods and QC figures, no statistics
- Automatic **multi-organism detection** (UniProt-style `_TAG` suffix in `Protein.Names`, e.g.
  host + pathogen), with optional per-organism analysis
- **Robustness score** (`Robustness_Score_*`, 0 = disabled with `n_iter_robustness: 0`):
  - DPC statistic: percentage of **leave-one-sample-out** refits (samples of the two compared
    conditions) where the protein still passes the volcano thresholds — flags calls carried by
    a single injection or replicate;
  - imputed analysis: repeated stochastic imputation (parallelised), same method as the
    reference run, seeds derived from `seed` (`Robustness_Score_imputed_*` when the DPC is the
    primary statistic)
- **Sensitivity to missing values** — comparison of the DPC and imputed analyses, and a test for
  on/off proteins (`dpc_sensitivity: true`, `differential_detection: true`, both on by default):
  - a **detection probability curve** (DPC: logit P(detected) = β0 + β1·log2 intensity) is fitted
    by maximum likelihood under the complete-normal model (Li & Smyth 2023);
  - **imputation-free differential analysis** (method of limpa, Li et al. 2025; independent
    Python implementation reproducing limpa 1.5.2): protein abundances estimated by maximum
    posterior under the DPC, with standard errors, then limma with vooma precision weights and
    residual-df correction for wholly unobserved groups. When `report.pr_matrix.tsv` is present
    (`dpc_level: auto`), proteins are quantified **from their precursors** (`dpcQuant`), which
    is where the method is strongest; otherwise from the protein matrix (`dpcQuantByRow`).
    Each protein × contrast is labelled *both / imputation only / DPC only* against the
    volcano rule;
  - **differential detection** for on/off proteins: probability of the non-detections in one
    condition if the protein had the abundance of the other (DPC, binomial), BH-corrected.
    An absence near the detection limit is not evidence; far above it, it is — unless the
    protein is prone to sporadic dropout: a protein-specific dropout rate, with a prior
    estimated on the dataset, absorbs identification failures and individual presence/absence
  Results go to three sheets (Sensitivity_Summary, Sensitivity_DPC, Differential_Detection)
- FDR correction per contrast or globally (BH)
- Visualisations:
  - **volcano plots** with −log10(raw p) on the Y axis (limma convention), colours from the
    volcano rule (p.adj by default) and a dashed line at the raw p equivalent to the p.adj
    threshold; `volcano_y_axis: padj` plots the adjusted p instead
  - **PCA**: score plot + protein **loading plot** on the same axes; the `pca_n_drivers`
    (default 20) most extreme proteins per direction (PC1±, PC2±) are labelled outside the
    point cloud and listed with their description (sheets PCA_Loadings, PCA_Scores;
    interactive in the dashboard)
  - UMAP (skipped below 5 samples), heatmaps, UpSet, scatter plots
- Co-expression: **WGCNA** (modules, hub scores, trait correlation)
- **Functional enrichment** (GO / KEGG / Reactome), two backends:
  - **gProfiler** (`go_backend: gprofiler`, default), with an automatic **coverage diagnostic**
    (g:Convert). If the organism is poorly indexed, Deimos offers **Perseverance** (bundled):
    orthologues in a well-annotated reference species are found by Reciprocal Best Hit (DIAMOND),
    gProfiler is queried on that species, and results are mapped back to your protein IDs.
    The reference FASTA and `perseverance_reference_organism` must name the same species
    (e.g. Arabidopsis proteome + `athaliana`; rice + `osativa` for grasses); identity and
    coverage filters: `perseverance_min_identity`, `perseverance_min_coverage`
  - **STRING** (`go_backend: string`) — recommended for bacteria and non-model organisms. The
    species is identified by NCBI taxid, detected automatically from UniProt accessions
    (`go_species_taxid: auto`). If the strain itself is not indexed in STRING, Deimos switches to
    an **orthology mode** (RBH with DIAMOND against the STRING proteome of the parent species).
    The same mode is used when fewer than 60 % of the accessions are indexed in STRING for the
    species (e.g. a recent UniProt proteome): the proteins are then mapped by sequence to the
    STRING proteome of the same species, so that all of them reach the enrichment. Perseverance
    belongs to the gProfiler route and is not proposed with STRING. With RefSeq / NCBI accessions
    (XP_…), which STRING cannot place, set `string_reference_taxid` to a related species indexed
    in STRING: the orthology mode is then used directly
  - **Your own annotation** (`perseverance_annotation_path`): local enrichment (hypergeometric +
    BH, no network) on a `protein_id, GO` table. An annotation delivered by a collaborator (Excel
    Trinotate/Blast2GO, topGO background file) is converted by `annotation_import.py` and, when it
    is keyed on another gene model set, transferred to your search IDs with `perseverance.py`
    (see `README_PERSEVERANCE.md`). `perseverance_go_terms_path` (table `GO, name, namespace`,
    written by `annotation_import.py --terms-out`) names the terms and splits them into
    GO:BP / GO:CC / GO:MF
  - **statistical universe = quantified proteins** by default (`enrichment_background:
    quantified`), for gProfiler (custom background restricted to annotated genes,
    `domain_scope='custom_annotated'`), STRING and the ANOVA-cluster enrichment: only
    detectable proteins can be called, so testing against the whole genome over-calls the
    categories of abundant proteins. `enrichment_background: genome` restores the previous
    behaviour
  - ORA per ANOVA heatmap cluster (`Cluster_Enrich_N` sheets)
  - optional rank-based **GSEA** (`run_gsea: true`, `gseapy`): all tested proteins ranked by
    −log10(raw p) × sign(log2FC), matched to the gene sets by gene symbol. Gene sets: human
    Enrichr libraries (GO 2021, KEGG 2021, Reactome 2022, default), mouse MSigDB collections
    (`gsea_organism: mouse`: hallmarks, Reactome, WikiPathways, GO), or local `.gmt` files
    (`gsea_gene_sets`, pinned version, no network). The log reports the share of ranked genes
    found in the gene sets
- **Protein × sample matrices for other analyses** (`export_matrices: true`, default), folder
  `matrices/`: log2 with missing values kept (MOFA+ / MEFISTO), log2 imputed (Mfuzz, WGCNA),
  DPC abundances with standard errors and observed flags, DPC condition means, sample table and
  a README on normalisation (DIA-NN only) and multi-omics integration
- Multi-sheet Excel export + interactive HTML dashboard, with a **light/dark theme toggle**
- **Pathfinder synthesis** as an optional last step (`make_pathfinder: true` / `--pathfinder`):
  core protein signature across contrasts, Meta-GO, WGCNA cross-reading, PCA × design factors.
  Pathfinder runs in a separate process on the saved workbook (no recomputation, no question);
  a failure is logged in `pathfinder/pathfinder_log.txt` and never stops Deimos
- **InSight biological reading** as an optional very last step (`make_insight: true` /
  `--insight`): response programmes (proteins that move together), functional themes (terms
  merged by shared proteins), condition map, themed volcanos and spotlight proteins, as PNG + SVG
  and a self-contained narrative report. Any design: trajectories for time courses, labelled maps
  and bar profiles when there is no order. Separate process, no recomputation, never stops Deimos
- **Machine-readable manifest**: a hidden `_manifest` sheet (schema `deimos-results` 1.0) lists
  the contrasts, the GO sheet → contrast mapping, the organism, the p-value adjustment and the
  volcano rule (also as structured keys since schema 1.2), so downstream tools read the workbook
  without guessing from sheet names and select the same proteins
- Reusable YAML configuration

---

## Input files

| File                     | Required | Description                                          |
| ------------------------ | -------- | ----------------------------------------------------- |
| `report.pg_matrix.tsv`   | ✅        | Protein × sample matrix (DIA-NN ≥ 2.5)               |
| `ExperimentalDesign.csv` | ✅        | Columns: `label;condition;replicate` (separator `;`), optionally `subject` and a batch column |
| `report.pr_matrix.tsv`   | ⬜ (recommended) | Precursor matrix: DPC primary statistic (precursor level) and DEqMS |
| Two FASTA files          | ⬜        | Your search FASTA + a reference proteome — only for Perseverance |
| Search FASTA             | ⬜        | STRING orthology mode: query sequences (otherwise fetched from UniProt) |
| `id_mapping` source      | ⬜        | Search FASTA or table to convert Protein.Group IDs (e.g. UniParc UPI -> UniProt / locus tag) before the analysis, see `id_remap.py` |

`label` must match exactly the sample name read from each DIA-NN column: the `Sple-…` token of
the run name (from `Sple-` up to the next `_`), or the full column name if there is none.
Recognised subject columns: `subject`, `pair`, `patient`, `individual`. Other names (`animal`,
`donor`, `mouse`...) trigger a warning asking you to rename the column.

---

## Installation

```bash
pip install -r requirements.txt
```

Or manually:

```bash
pip install pandas numpy scipy scikit-learn umap-learn matplotlib seaborn \
  adjustText openpyxl pillow PyComplexHeatmap pyyaml gprofiler-official
```

Optional:

```bash
pip install gseapy                  # rank-based GSEA (run_gsea: true)
sudo apt install diamond-aligner    # or: conda install -c bioconda diamond
                                    # Windows: diamond.exe next to deimos.py or in the PATH
                                    # needed by Perseverance and by the STRING orthology mode
```

Enrichment needs network access to the gProfiler, STRING and UniProt APIs; GSEA to Enrichr or
MSigDB (data.broadinstitute.org), unless local `.gmt` files are given.

**Dashboard offline** — place `chart.umd.min.js` and `xlsx.full.min.js` alongside `build_dashboardv7.py`:

- <https://cdn.jsdelivr.net/npm/chart.js@4.4.1/dist/chart.umd.min.js>
- <https://cdn.jsdelivr.net/npm/xlsx@0.18.5/dist/xlsx.full.min.js>

(If absent, the dashboard falls back to loading both from cdnjs.cloudflare.com.)

**Rebuild the dashboard only** — after a dashboard update, regenerate it from a finished run
without re-running the pipeline: `python rebuild_dashboard.py proteomics_output`. The
dashboard applies the run's selection rule everywhere (counts, DEP table, exports) and embeds
every protein for its exports; only the volcano cloud is subsampled for display.

---

## Usage

```bash
# Interactive mode (prompted at startup)
python deimos.py

# YAML config mode
python deimos.py --config config_example.yaml

# Interactive + save config for future runs
python deimos.py --save-config

# Override paths or options without editing the config
python deimos.py --config project.yaml --tsv /data/run2/report.pg_matrix.tsv --out-dir /data/run2/out
python deimos.py --config project.yaml --batch-column batch
python deimos.py --config project.yaml --qc-only
python deimos.py --config project.yaml --pathfinder     # + Pathfinder synthesis at the end
python deimos.py --config project.yaml --insight        # + InSight figures and report at the end
```

**InSight** is looked for in `insight_path`, then `$INSIGHT_HOME`, then an `InSight` folder next
to the Deimos folder. It runs after Pathfinder with the same Python interpreter and reads the saved
workbook and `matrices/` (keep `export_matrices: true` for the abundance figures). Options go in
`insight_options` (`k`, `max_volcanos`, `volcanos`, `categorical`, `order_regex`,
`reference_fasta`, `title`); full output in `insight/insight_log.txt`. It is skipped in QC-only mode.

**Pathfinder** is looked for in `pathfinder_path`, then `$PATHFINDER_HOME`, then a `Pathfinder`
folder next to the Deimos folder (`git clone https://github.com/BenoitBe/Pathfinder` beside it).
It runs with the same Python interpreter; its dependencies are already Deimos' (plus `requests`
for the optional Publications axis). It is skipped in QC-only mode (no p-values). Pathfinder ≥ 1.2
and InSight select **the same proteins as Deimos**: they apply the volcano rule written in the
`_manifest` (statistic, threshold and ratio; structured keys since schema 1.2), so their counts
match the volcanos, UpSet and dashboard, and they read GO / pathway terms at q < 0.05 (gProfiler /
STRING). `pathfinder_options: {pthr, min_ratio}` overrides the Pathfinder protein rule (the
Methods sheet and its summary then say the counts differ). Pathfinder < 1.2 used `p.adj < 0.05`
without a fold-change criterion; the Methods sheet states whichever rule the installed version
applies.

**Unattended runs** (scheduler, cluster, no terminal): some options still ask a question at
runtime when left unset, and the run stops with `EOFError` if nobody can answer —
`run_perseverance` (gProfiler backend, asked at every run when `"__auto__"`), `contrasts` (asked
above 30 contrasts) and `species_split` (asked when several `_TAG` organism suffixes are found,
e.g. a FASTA with contaminants). Set them explicitly in the YAML (`config_example.yaml` sets
`run_perseverance: false`).

**Contrast names**: conditions starting with a digit get an `X` prefix (R convention:
`102T0Light` → `X102T0Light`); patterns such as `*T0Light_vs_TronT0Light` or
`{g}T3Dark_vs_{g}T3Light` match either way, but the Excel and dashboard names carry the `X`.

Check a paired design before running: `python check_paired.py ExperimentalDesign.csv`.

The statistical engine (`limma_ebayes.py`) is a Python reimplementation of `limma`
(`lmFit`, `contrasts.fit`, `eBayes`, BH) and of `DEqMS::spectraCounteBayes` — **no dependency
on R**. It reproduces R exactly (see *Validation*); the reference outputs of R limma 3.58.1 and
DEqMS 1.31.0 on simulated data are in `tests/fixtures/limma_r/` with the scripts that made them,
and `tests/test_limma_reference.py` checks them at every test run (plus InMoose when installed).
DEqMS uses R's loess through `scikit-misc` (same Fortran code); without it, a direct loess is
used (log-variance trend within 0.05 of R, same calls in tests).

---

## Validation

Benchmark on real DIA-NN data of the platform, without a spike-in mixture (October 2026;
scripts and full tables in the project notes). Three datasets: cell lines (P1, 5 % missing
values after filtering), a non-model marine invertebrate (P3, 23 %), a marine species with
16 replicates per condition (P4, 6 %).

- **Null comparisons**: two groups of the same biology (same cell line balanced for treatment
  and injection series, same site and time, same condition balanced for block), 3 vs 3.
  Every significant protein is a false positive.
- **In-silico spike-ins**: 10 % of the proteins lowered by 0.4-2.5 log2 in one group, values
  falling below the detection limit removed according to the dataset's own detection curve.

| BH 5 %, 3 vs 3 | P1 | P3 | P4 |
|---|---|---|---|
| *False positives per null comparison (% of comparisons with at least one)* | | | |
| MIXED, rule before 10/2026 (`mixed_legacy`) | 50 (100 %) | 116 (100 %) | 87 (100 %) |
| MIXED, corrected rule | 3.4 (79 %) | 3.4 (60 %) | 7.0 (90 %) |
| QRILC | 0.4 (17 %) | 0.1 (8 %) | 2.3 (60 %) |
| **DPC, precursor level (primary statistic)** | **0 (0 %)** | **0.1 (5 %)** | **0.05 (5 %)** |
| *Spike-ins: real FDP at nominal 5 % / power at real FDP 5 %* | | | |
| MIXED, rule before 10/2026 | 13 % / 0.1 % | 66 % / 0.1 % | 23 % / 0 % |
| MIXED, corrected rule | 4.7 % / 57 % | 22 % / 1 % | 11 % / 22 % |
| QRILC | 2.6 % / 58 % | 8 % / 1 % | 7 % / 26 % |
| **DPC, precursor level** | **0.7 % / 72 %** | **3 % / 7 %** | **4 % / 52 %** |

Fold-changes (slope of estimated vs true, 1 = exact): DPC 0.87-0.90, QRILC 1.08-1.19, MIXED
0.88-0.95. The DPC fold-changes are slightly compressed (shrinkage of weakly observed values).
The Python implementation gives the same significant calls as limpa 1.5.2 (R) on three datasets;
on the last one (B. napus embryos, 24 samples, 78,495 precursors, 16 contrasts) the precursor
DPC is identical, |Δ log2FC| ≤ 0.01, |Δ log10 p| ≤ 0.03 and the calls are identical in all 16
contrasts.

**limma and DEqMS engine** (`limma_ebayes.py`) against R limma 3.58.1, DEqMS 1.31.0 and InMoose
0.9.1, on 7 simulated designs (two groups, four groups with six contrasts, paired, exactly zero
variances, homogeneous variances with infinite prior df) and on the four real datasets above
(6 to 66 contrasts, 801 to 8,687 proteins): coefficients, moderated t, p-values and BH-adjusted
p-values identical to 1e-12, same prior df and prior variance, not a single different call.
Two corrections came out of this check: zero variances are now offset as in limma instead of
being dropped (prior df 3.4 instead of 1.25 on the test case), and when the prior df is infinite
the prior variance and total df follow limma (`mean(s²)`, pooled df). DEqMS: the previous port
approximated the loess and the prior df (275 calls out of 19,900 different on P1); it is now
exact. On two of the real datasets (P2, P3), most proteins have a single peptide and R DEqMS
itself fails (loess degenerate); Deimos's DEqMS diagnostic already switches DEqMS off there.

**Differential detection** (on/off proteins), null comparisons: without the sporadic dropout
term, 11-37 false calls per comparison on all datasets; with it, 0-0.06 at 3 vs 3 and 0.3 at
8 vs 8, with 84-97 % of the power kept on simulated "off" proteins (none on P3, where abundant
proteins often miss: the test then stays silent rather than wrong).

---

## Reproducibility

- The DPC primary statistic involves no random draw: same input, same results.
- `seed` fixes the reference imputation draw and the robustness iterations. With the MIXED
  imputation, the kNN part (deterministic) is computed once and reused by the robustness
  draws: same results, the step no longer grows with a kNN per iteration on large designs.
- The Methods sheet records the imputation method (and why `auto` chose it), the seed, the design
  matrix, the selected contrasts and the design checks.
- `last_config.yaml` in the output folder reloads the exact configuration of the last run.

---

## Repository structure

```
deimos-proteomics-pipeline/
├── deimos.py                     # Main orchestrator
├── config.py                     # YAML config management, CLI arguments
├── config_example.yaml           # Annotated YAML template
├── limma_ebayes.py               # Statistical engine (limma + DEqMS, design matrix, contrasts)
├── go_enrichment.py              # GO/KEGG enrichment: gProfiler, STRING, local ORA, coverage diagnostic
├── gsea_enrichment.py            # Rank-based GSEA
├── dashboard_integration.py      # Excel -> dashboard bridge
├── pathfinder_bridge.py          # _manifest sheet + Pathfinder launch (separate process)
├── insight_bridge.py             # InSight launch (separate process, after Pathfinder)
├── dpc_sensitivity.py            # DPC (limpa method): primary statistic, sensitivity, detection
├── build_dashboardv7.py          # HTML dashboard generator
├── dashboard_template.html       # Dashboard template (light/dark theme)
├── diagnostic_dashboard.py       # Dashboard diagnostic utility
├── rebuild_dashboard.py          # Regenerate the dashboard of a finished run
├── check_paired.py               # Pre-flight check of a paired design
├── check_gsea.py                 # GSEA integration diagnostic
├── diag_go.py                    # GO via orthologues: ID format, gProfiler coverage, backgrounds
├── perseverance.py               # Orthology search CLI (RBH via DIAMOND)
├── annotation_import.py          # Collaborator annotation (Excel/topGO) -> Perseverance/Deimos TSV
├── rbh.py / transfer.py          # Perseverance internals
├── deimos_to_perseverance.py     # Deimos <-> Perseverance bridge (FASTA, ID mapping)
├── id_remap.py                   # id_mapping: Protein.Group ID conversion in pg + pr matrices
├── README_PERSEVERANCE.md        # Perseverance documentation
├── README_INTEGRATION.md         # Perseverance <-> Deimos integration notes
├── deimos_mark.svg / deimos_icon.svg
├── CHANGELOG.md
├── requirements.txt
├── LICENSE
├── tests/                        # pytest: limma/DEqMS vs R, DPC, volcanos, PCA, GSEA, STRING plan,
│                                 #   matrices, manifest contract, Pathfinder launch
├── toy_data/                     # Toy FASTA files for Perseverance (README_PERSEVERANCE)
└── example_data/                 # Minimal synthetic dataset (smoke test, no biological signal)
    ├── report.pg_matrix.tsv
    └── ExperimentalDesign.csv
```

---

## Output

```
proteomics_output/
├── ProteomicAnalysis_Results.xlsx   # Full multi-sheet report
├── ProteomicAnalysis_QC.xlsx        # QC-only mode (instead of the above)
├── proteogen_dashboard.html         # Interactive dashboard (light/dark toggle)
├── matrices/                        # export_matrices: protein x sample TSV matrices + README
├── string_orthologs/                # STRING orthology mode: proteome cache + DIAMOND files
├── perseverance/                    # Perseverance: DIAMOND outputs, orthologue table
├── pathfinder/                      # make_pathfinder: report, synthesis workbook, summary,
│                                    #   design used, full log
├── insight/                         # make_insight: figures 01-06 (PNG + SVG),
│                                    #   insight_report.html, insight_log.txt
└── last_config.yaml                 # Last run config (reloadable)
```

**Excel sheets** (full mode): Methods_Upstream · Methods · raw_data · Log2_Impute · QC · PCA_UMAP · PCA_Loadings · PCA_Scores ·
UMAP · Scatter_Plots · Differential_Expression · Volcano_Plots · UpSet_Intersections ·
Zscore_Heatmap · ANOVA_Results · ANOVA_Clusters, plus depending on the options: WGCNA_* · GO_* ·
Cluster_Enrich_* · GSEA_* · Orthologues_Perseverance · Orthologues_STRING ·
Sensitivity_Summary · Sensitivity_DPC · Differential_Detection, and a hidden `_manifest` sheet
(last; unhide it in Excel to read it)

**QC-only mode**: Methods_Upstream · Methods · raw_data · Log2 (not imputed) · QC · PCA_Correlation
(complete cases, ≥ 3 samples)

---

## YAML configuration

Copy `config_example.yaml` and adapt; every key is documented there. The keys that most often
change between projects:

```yaml
volcano_use_padj:     true        # volcano rule on p.adj (true) or raw p (false)
volcano_p_thresh:     0.05
volcano_ratio_min:    2.0
volcano_y_axis:       pvalue      # Y axis: raw p (default) | padj
pca_n_drivers:        20          # proteins labelled on the loading plot (PC1±, PC2±)

primary_statistic:    dpc         # dpc (no imputation, precursors) | imputation
impute_method:        auto        # auto | qrilc | mixed | mixed_legacy
seed:                 42          # reproducible imputation

contrasts:            "__auto__"  # "__auto__", "all", or a list of patterns:
#  - "Live_{s}_vs_Heat_killed_{s}"   # {x} = same value on both sides
#  - "Live_*_vs_Live_*"              # *   = any value

control_condition:    "__auto__"  # "__auto__" (auto-detect), an explicit name, or false
batch_column:         null        # column of ExperimentalDesign.csv modelled as a covariate
use_deqms:            false       # DEqMS, gated by the automatic reliability diagnostic

go_organism:          null        # gProfiler species code (e.g. "bnapus"), null = disabled
go_backend:           gprofiler   # gprofiler | string
go_species_taxid:     null        # STRING: "auto" or an NCBI taxid
enrichment_background: quantified # ORA universe: quantified (tested proteins) | genome
run_perseverance:     "__auto__"  # "__auto__" (asks at runtime), true/false (gProfiler only)
perseverance_reference_organism: null  # gProfiler code of the reference FASTA (e.g. athaliana)
run_gsea:             false
gsea_organism:        human       # human (Enrichr) | mouse (MSigDB)
gsea_gene_sets:       null        # null = defaults; or .gmt files / "msigdb:<collection>"
export_matrices:      true        # matrices/ for Mfuzz, MOFA+/MEFISTO, multi-omics

qc_only:              false
make_wgcna:           false
make_dashboard:       true
dpc_sensitivity:      true        # DPC vs imputed comparison sheets
differential_detection: true      # on/off proteins, detection test (DPC)
dpc_level:            auto        # auto (precursors if pr_matrix) | precursor | protein
make_pathfinder:      false       # true = Pathfinder synthesis at the end (needs a clone)
make_insight:         false       # true = InSight biological-reading figures at the end (needs a clone)
```

Any extra YAML key is passed through to the pipeline even if not listed in `config.py`'s
defaults — no code change needed to use a parameter already read via `params.get(...)`. The
flip side: a misspelt key is ignored without error.

---

## Ecosystem

| Tool | Role |
|---|---|
| **Deimos** (this repo) | Full DIA quantitative pipeline |
| [Phobos](https://github.com/BenoitBe/phobos-peptidomics-pipeline) | DDA/PEAKS peptidomics + PTM pipeline |
| [Pathfinder](https://github.com/BenoitBe/Pathfinder) | Cross-contrast synthesis of one Deimos/Phobos workbook (core signature, meta-GO, WGCNA cross-reading), no recomputation. Can run as Deimos' last step: `make_pathfinder: true` |
| InSight | Biological reading of a finished Deimos run (response programmes, functional themes, condition map, themed volcanos, spotlight), any design, no recomputation. Can run as Deimos' very last step: `make_insight: true` |
| **Perseverance** | Orthology search (RBH/DIAMOND) for non-model organisms — bundled directly in this repository, see `README_PERSEVERANCE.md` |

---

## Name

**Deimos** — moon of Mars, small and precise.
Acronym: **D**IA **E**xpression **I**ntegrated **M**ulti-**O**mics **S**uite.

---

## License

MIT — see [LICENSE](LICENSE).
