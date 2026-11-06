# Changelog — Deimos

All notable changes to the Deimos proteomics pipeline.
Repository: https://github.com/BenoitBe/deimos-proteomics-pipeline

---

## [Unreleased] — Major update

This release bundles new statistical capabilities, several correctness fixes,
and dashboard improvements accumulated since the previous public version.

### Changed

- **Primary statistic: imputation-free DPC analysis by default** (`primary_statistic: dpc`).
  The limpa-method chain of `dpc_sensitivity.py` (precursor-level detection probability
  curve, protein quantification with standard errors, weighted limma) now feeds the main
  columns `<C>_diff`, `<C>_p.val`, `<C>_p.adj`, `Pi_Score_<C>` and `Robustness_Score_<C>`,
  so volcanos, scatter plots, UpSet, GO, GSEA, dashboard and Pathfinder follow it. Nothing
  is removed: the imputed analysis (limma or DEqMS) is kept in `<C>_diff_imputed`,
  `<C>_p.val_imputed`, `<C>_p.adj_imputed` (`_p.*_imputed_limma` with DEqMS),
  `Pi_Score_imputed_<C>`, `Robustness_Score_imputed_<C>`, and still feeds PCA, UMAP,
  correlation, ANOVA heatmaps and WGCNA. `primary_statistic: imputation` restores the
  previous behaviour exactly. Without a usable `pr_matrix` the run falls back to the imputed
  analysis (console warning + Methods line): the protein-level DPC was not calibrated on
  data with many missing values. Benchmark on real data: README, section *Validation*.
  - Proteins absent from the pr_matrix are analysed at protein level (same model, common
    BH correction) instead of being dropped.
  - Robustness score of the DPC statistic: percentage of leave-one-sample-out refits
    (samples of the two compared conditions) where the protein still passes the volcano
    thresholds (0-100). The statistic involves no random draw, so the re-imputation score
    has no direct equivalent.
  - Methods sheet, Sensitivity_Summary, console and `_manifest` (`primary_statistic`)
    state which statistic is primary; interactive mode asks for it.
- **MIXED imputation rule corrected**: a protein absent from all replicates of at least one
  condition now has *all* its missing values imputed by QRILC; kNN is used only for proteins
  with partial missingness. The former rule imputed the absent condition low (QRILC) and the
  partial missingness of the other condition near its observed values (kNN), which created
  clean, low-variance differences out of nothing: 50-120 false positives per null
  comparison on real data. Kept as `impute_method: mixed_legacy` to reproduce earlier
  analyses.
- pr_matrix runs are matched to the samples through the pg_matrix of the same search
  (identical run name, else same injection when a label was corrected in the pg_matrix
  only). A label carried by two runs of the pr_matrix is never picked at random: the
  sample stays unmatched and the run falls back, with a warning. cRAP precursors are
  removed as in the pg_matrix.

- **Differential detection calibrated on real data**: probability of non-detection
  = pi + (1 - pi) x (1 - DPC), with pi a protein-specific sporadic dropout rate
  (identification failures, proteins present in some individuals only). Its prior is
  estimated on each dataset from abundant protein x condition cells (nonparametric
  maximum likelihood on a grid, truncated binomial), then updated by the non-detections
  of the condition where the protein is present. Without this term, abundant proteins
  missing in a few runs of one group were judged impossible under the DPC: 11-37 false
  calls per null comparison on three real datasets (none in simulation, where the model
  was true). Now 0-0.06 per comparison at 3 vs 3 and 0.3 at 8 vs 8, with 84-97 % of the
  power kept on simulated "off" proteins. The prior is reported in Methods.

- **Enrichment universe = quantified proteins by default** (`enrichment_background:
  quantified`; `genome` restores the previous behaviour). In proteomics only detectable,
  mostly abundant, proteins can be called: against the whole genome, the categories of
  abundant proteins (ribosome, translation, glycolysis...) come out enriched for that
  reason alone.
  - gProfiler: custom background = all proteins tested in the differential analysis (or
    their orthologues in Perseverance mode), with `domain_scope='custom_annotated'`
    (background restricted to annotated genes). The `gprofiler-official` client forces
    `'custom'` (background genes without annotation counted in the domain, which makes
    the test optimistic); Deimos rewrites it for the duration of the call.
  - STRING follows the same key; `string_background` (now `null` by default) still forces
    the STRING universe when set, with a console note when it differs.
  - ANOVA-cluster enrichment: the universe is all tested proteins. Previously, with
    `string_background: quantified`, the STRING background was reduced to the proteins of
    the clusters themselves (ANOVA-significant only), which biased the test.
  - Methods: background row for gProfiler; interactive mode asks the universe for both
    backends (default quantified).

- **limma engine checked against R and corrected where it differed** (`limma_ebayes.py`).
  Against R limma 3.58.1 and InMoose 0.9.1 on 7 simulated designs and 4 real datasets:
  identical to 1e-12 (coefficients, moderated t, p, BH) once two edge cases were fixed:
  - exactly zero variances are offset to 1e-5 x median as in limma's `fitFDist`, instead
    of being excluded from the prior estimate (prior df 3.4 -> 1.25 on the test case, as R);
  - infinite prior df (homogeneous variances): prior variance = mean(s²) and total df =
    pooled df, as limma, instead of exp(mean log) and df + 1e6;
  - total df capped at the pooled df; residual df = n - rank(design); BH ignores NaN
    p-values as `p.adjust` does.
- **DEqMS reproduced exactly** (`spectraCounteBayes`, DEqMS 1.31.0): R's loess (via
  `scikit-misc`, same Fortran code; direct loess fallback), DEqMS's grid search for the prior
  df, prior variance with its digamma correction. The previous port used a smoothing spline,
  a different prior-df formula capped at 4 x residual df and no correction: on P1, 275 of
  19,900 calls differed from R. Now identical to R to 1e-12 on all datasets where R DEqMS runs.
- DEqMS peptide count: minimum across the **analysed** samples (runs matched through the
  pg_matrix) + 1 for every protein, as in the DEqMS LFQ vignette. Previously the minimum was
  taken over every run of the pr_matrix, including samples excluded from the design, and
  only zero counts were raised to 1. `Peptide_Count_min` keeps the raw minimum.
- Methods: "DEqMS settings" (count definition, prior df, loess engine) and "Implementation"
  rows.

- Pathfinder launch: stdin is now an empty, closed pipe instead of `DEVNULL`. On Windows
  `NUL` is a character device and `isatty()` returns True there, so Pathfinder believed it
  had a terminal (no hang: the read hits end of file at once, but a question could be
  printed). Found by the test suite on Windows.

### Added

- **RBH progress** (Perseverance and STRING orthology): DIAMOND runs quietly, so Deimos now prints
  each of the four steps (two databases, two blastp searches) with sequence counts, threads and
  duration — a long reverse search no longer looks like a freeze. DIAMOND uses all cores minus one
  by default (was 4); `perseverance_threads` sets it.

- **Collaborator annotations** (`annotation_import.py`): converts an annotation delivered as a
  table (Excel Trinotate / Blast2GO, CSV, TSV, topGO background `ID <tab> GO,GO`) into the
  Perseverance / Deimos format (`protein_id, GO, KEGG, Pfam`) and a term table
  (`GO, name, namespace`). GO found in any cell format, hits filtered on e-value (default
  1e-5), species KEGG pathways mapped to generic maps, IDs normalised (Trinity / AUGUSTUS forms,
  typos) and rewritten as in a reference FASTA; optional OBO: obsolete terms removed, alt_id
  replaced, propagation (is_a + same-namespace part_of + logical definitions), names. Warns when
  propagated and direct annotations are mixed without an OBO.
- **Named local GO** (`perseverance_go_terms_path`): the local ORA (`perseverance_annotation_path`)
  names the terms and splits them into GO:BP / GO:CC / GO:MF instead of GO IDs under a single
  `GO:PERSEVERANCE` source. Unchanged without the key.
- Local ORA: a protein group whose leader has no annotation takes the GO of its other annotated
  members (RefSeq isoforms, raw Perseverance output); nothing is removed.

- **Same proteins downstream** (manifest schema 1.2): the volcano rule is also written as
  structured keys (`volcano_p_column`, `volcano_p_thresh`, `volcano_ratio_min`) next to the
  `volcano_rule` text, so Pathfinder ≥ 1.2 and InSight select exactly the Deimos proteins
  (`p < thr & |log2FC| >= log2(ratio)`, as in `volcano_status`). New
  `pathfinder_options.min_ratio` (with `pthr`, an explicit override). The Methods row of the
  Pathfinder step now states the rule the installed Pathfinder applies (volcano rule for ≥ 1.2,
  `p.adj < pthr` without fold-change criterion for older versions).

- **InSight step** (optional, very last, after Pathfinder): `make_insight: true` / `--insight`
  runs the separate InSight module on the saved workbook and `matrices/` — response programmes,
  functional themes, condition map, themed volcanos, spotlight proteins (PNG + SVG) and a
  self-contained narrative report in `<out_dir>/insight/`. New `insight_bridge.py` (standard
  library only, same contract as `pathfinder_bridge.py`: separate process, empty stdin, timeout
  `insight_timeout`, full log, a failure never stops Deimos). Keys `make_insight`,
  `insight_path` (then `$INSIGHT_HOME`, then `../InSight`), `insight_options`
  (`k`, `max_volcanos`, `volcanos`, `categorical`, `order_regex`, `reference_fasta`, `title`,
  validated), `insight_timeout`; interactive question only when a clone is found (default no);
  Methods row "Biological reading"; skipped in QC-only mode. Nothing existing is changed.

- **Identifier mapping** (`id_mapping`, module `id_remap.py`): converts the
  Protein.Group identifiers of the pg and pr matrices (e.g. UniParc UPI -> UniProt
  accession or locus tag) from the search FASTA (`AC=`, `GN=` fields) or a table, without
  re-running DIA-NN; same mapping for both files, intensities unchanged, collisions keep
  their original ID, idempotent; converted matrices in `<out_dir>/id_mapping/`, Methods
  row. The call in `deimos.py`, lost once when the file was replaced, is now covered by a
  test.

- **Protein x sample matrices for other analyses** (`export_matrices: true`, default):
  folder `<out_dir>/matrices/` with TSV files sharing the same rows (unique `protein`
  ID, Protein.Group, Genes, Description) and the sample order of `samples.tsv`:
  `log2_filtered_NA.tsv` (log2, NA kept, `num_NA`; for MOFA+/MEFISTO, which model
  missing values), `log2_imputed.tsv` (as the Log2_Impute sheet; Mfuzz, WGCNA),
  `dpc_abundance.tsv` / `dpc_standard_error.tsv` / `dpc_observed.tsv` (precursor-level
  DPC quantification of the primary statistic: complete matrix with uncertainties) and
  `condition_means_dpc.tsv` (Mfuzz on condition means), plus `samples.tsv` and a README
  on normalisation (DIA-NN only), file choice and multi-omics integration. Methods row.

- **GSEA gene sets by organism, local GMT files.** New keys (the previous behaviour
  stays the default):
  - `gsea_organism: human | mouse`. `human` = Enrichr GO 2021 + KEGG_2021_Human +
    Reactome_2022 (unchanged); `mouse` = MSigDB mouse collections MH, M2
    CP:REACTOME, CP:WIKIPATHWAYS, M5 GO:BP/MF/CC (`gsea_msigdb_version`, default
    2025.1.Mm; MSigDB has no mouse KEGG).
  - `gsea_gene_sets`: list of local `.gmt` files (pinned version, no network),
    `msigdb:<collection>` or Enrichr library names. Gene sets are loaded once for all
    contrasts, symbols upper-cased like the ranking (mouse `Actb` = human `ACTB`); an
    unavailable source is reported and skipped.
  - `gsea_fdr`, `gsea_permutations` exposed in the configuration.
  - Log line: share of the ranked genes present in the gene sets.
  - Interactive mode: after "Also run GSEA?", question "Gene sets [1/2]" (human /
    mouse), saved in last_config.yaml as `gsea_organism`.
- **GSEA ranking on the raw p-value** (`gsea_rank_p: pvalue`, default). The ranking used
  the volcano p, i.e. the BH-adjusted p with `volcano_use_padj: true`, whose plateaus
  (many proteins with the same value) flatten the middle of the ranking, where GSEA
  picks up coordinated small shifts. `gsea_rank_p: volcano` restores the previous rule.
- `check_gsea.py`: no longer reports a missing GSEA call for the current deimos.py
  (`(go_params or {}).get("run_gsea")`).

- **PCA loading plot.** Score plot (samples) and loading plot (proteins) side by side, from
  the same PCA as `plot_pca` (same axes, same signs): a protein on the right of B is more
  abundant in the samples on the right of A. ANOVA-significant proteins are coloured by
  the condition where they are most abundant, with the score-plot colours (a single
  contrast: volcano rule, higher in A / higher in B); the coloured proteins farthest
  from the origin are labelled. Figure `plot_pca_loadings.png` (sheet PCA_UMAP) and
  sheet `PCA_Loadings` (PC1-PC3 loadings, group, ANOVA cluster, sorted by distance to
  the origin). Methods row added.
  - **Driving proteins**: the `pca_n_drivers` (default 20) most extreme loadings, split
    between PC1+, PC1-, PC2+ and PC2-, are circled and labelled, and listed under the
    figure with their description and the conditions lying on that side of the score
    plot (column `driver` in PCA_Loadings). Sample scores in sheet `PCA_Scores`.
  - **Driver labels placed outside the cloud** (no overlap): PC1+ / PC1- in a column
    right / left, PC2+ / PC2- in a row above / below (two staggered rows when space is
    short), each linked to its point, in the order of the points so leader lines do not
    cross; axis limits widened to fit. Replaces adjustText for these labels, which left
    them stacked when the extreme proteins sit against the edge of the axes (8-condition
    run: 4 PC1- labels on top of each other). Sample labels of the score plot are spread
    with adjustText (replicates often overlap), with a fixed offset as fallback.
  - **Empty gene names / descriptions** no longer crash the loading plot
    (`TypeError: object of type 'float' has no len()`): with pandas >= 3, `astype(str)`
    keeps missing values as NaN instead of writing "nan". Empty cells are now read as
    empty strings; a protein without gene name is labelled with its accession.
  - **Dashboard** (QC / PCA tab): the same score plot + loading plot, read from
    PCA_Scores / PCA_Loadings (same PCA, axes and signs as the Excel figure), with the
    driver table; hover for gene, description, group and loadings, click a protein (in
    the plot or the table) to find it in the volcano. Non-significant proteins are
    subsampled to 3,000 points (drivers always kept). Hidden for older workbooks.
- **Dashboard condition palette** had 6 colours and cycled: beyond 6 conditions two
  conditions shared a colour (8-condition run: embAviso0 = embMajor24). 12 colours now,
  the first 6 unchanged.

- `tests/test_limma_reference.py` and `tests/fixtures/limma_r/` (R reference outputs on
  simulated data, with the scripts that made them): limma and DEqMS against R, InMoose
  when installed, BH with missing values, DEqMS grid, rank-deficient design.
- `tests/test_enrichment_background.py` (8 tests, API simulated): `custom_annotated` and
  background sent to g:GOSt, universe for contrasts, orthologues and clusters, STRING
  override, config.
- **Validation section** in the README: null comparisons and in-silico spike-ins on three
  real datasets for MIXED (both rules), QRILC and the DPC statistic.
- `tests/test_primary_statistic.py` (14 tests): primary analysis with protein-level
  fallback, common BH, identity with `dpc_de`, leave-one-out score, column switch with
  nothing lost, corrected MIXED rule, run matching, fallback without pr_matrix, config;
  plus a sporadic-dropout test for the differential detection.

- **Sensitivity to missing values** (`dpc_sensitivity.py`), two analyses *added* to the
  main one — imputation + limma/DEqMS stays the reference for volcanos, robustness, GO,
  dashboard and Pathfinder; no existing sheet or column changes (checked: the 28 existing
  sheets are identical with the option on or off).
  - Detection probability curve (DPC), maximum likelihood under the complete-normal model
    (Li & Smyth 2023), fitted on the matrix *before* filtering (filtering on detection
    counts would bias it). Identical to limpa's `dpcCN` on test data.
  - Imputation-free differential analysis: maximum-posterior completion with standard
    errors, vooma weights predicted by intensity and standard error, residual-df
    correction for wholly unobserved groups, eBayes with per-protein df (method of
    limpa `dpcQuantByRow` + `dpcDE`; same calls as R limma on test data, max
    |Δlog10 p| = 0.002). Same design and contrasts as the main model.
  - Differential detection for on/off proteins: probability of the non-detections in one
    condition at the other's abundance (DPC, binomial, integrated over the uncertainty
    of that abundance), two-sided, BH over all proteins. Verdict per pair: confirmed /
    detection only / imputation only (absence compatible with the detection limit).
  - Precursor level (`dpc_level: auto`, default when `report.pr_matrix.tsv` is present):
    proteins quantified from their precursors by maximum posterior (additive precursor +
    sample model, method of limpa `dpcQuant`), with standard errors that reflect how many
    precursors were seen. Checked against limpa 1.5.2 on two real datasets: same DPC,
    abundances within 1e-2 (median 3e-5), identical standard errors, and the same 1,778
    significant protein × contrast pairs; 25x faster than limpa in R.
  - The variance trend uses exact Python ports of R `lowess` and limma `weightedLowess`
    (identical to R to 1e-14).
  - Sheets `Sensitivity_Summary` (DPC, plot, counts per contrast, robustness of
    concordant vs imputation-only calls), `Sensitivity_DPC`, `Differential_Detection`;
    Methods section; `_manifest` schema 1.1 (`sensitivity` rows); YAML keys
    `dpc_sensitivity` and `differential_detection` (default true). Non-blocking.
  - Simulation with known truth (logistic MNAR, 4 replicates): the imputation-free
    analysis recovered 73–80 % of true changes vs 52–60 % for QRILC + limma at 3–6 %
    false discoveries; differential detection made no false call on on/off proteins.

- **Pathfinder as an optional last step** (`make_pathfinder: true`, `--pathfinder`,
  or a `(y/N)` question in interactive mode when a clone is found). Pathfinder runs
  in a separate process on the workbook just saved (`pathfinder_bridge.py`): no
  question (`stdin` closed), no recomputation, UTF-8 output forced (Windows), a
  timeout (`pathfinder_timeout`), and any failure is reported as `[WARN]` with the
  log path — Deimos always completes. The GO workbook is always passed explicitly
  (`--go <workbook>` or `--go none`), so Pathfinder never picks up the GO results of
  another project lying around. The design is passed without the `replicate` index.
  Skipped in QC-only mode. One `pathfinder/` folder per organism in multi-organism
  runs. Publications axis off by default (`pathfinder_publications`); in Perseverance
  mode no species is passed (the GO is on the reference species). Methods records
  the Pathfinder version (`Downstream synthesis`).
- **`_manifest` sheet** (hidden, last): one fact per row (`section | key | value`),
  schema `deimos-results` 1.0 — Deimos commit, statistics sheet, primary statistic,
  p-value adjustment, volcano rule, sample matrix and imputation, contrasts, GO
  sheet → contrast mapping, GO backend, organism (STRING taxid or g:Profiler code),
  design factors. The shortened GO sheet names cannot be parsed back reliably
  (Pathfinder recognised 1/6 to 0/12 of them); the manifest is now the source of
  truth. `export_go_sheets()` returns the `{contrast: sheet}` mapping it wrote.
- **First tests in the repository** (`tests/`, pytest): manifest contract on five
  sheet-naming situations, round trip through Pathfinder's reader when a clone is
  present, launch step (success, failure, timeout, missing clone, arguments).

- **Paired design support (repeated measures).** When the experimental design
  contains a `subject` column (aliases: `pair`, `patient`, `individual`), Deimos
  automatically switches to a paired model `~0 + condition + subject`. Each
  subject gets its own fixed-effect intercept, absorbing inter-individual
  baseline variability (generalized paired t-test). Triggered only when a subject
  spans ≥2 conditions; otherwise the standard `~0 + condition` model is used.
  Subject IDs are normalized (trailing spaces and `.0` float suffixes removed),
  and a warning is emitted for empty `subject` values (e.g. CSV separator typos).
  `check_paired.py` lets you verify the design *before* running.

- **Rank-based GSEA (`gsea_enrichment.py`).** Complements the existing ORA. Ranks
  ALL proteins by `-log10(p) × sign(LFC)` and detects coordinated pathway shifts
  even below the DEP threshold, using `gseapy.prerank` (Broad algorithm). Gene
  sets from Enrichr (GO / KEGG / Reactome). Optional (asked at launch), fully
  non-blocking: the pipeline runs normally if `gseapy` is absent or offline.
  Output: `GSEA_<contrast>` sheets + NES bar plots.

- **GO enrichment per ANOVA heatmap cluster.** Each co-expression cluster from
  the ANOVA heatmap is functionally annotated (ORA), exported as
  `Cluster_Enrich_N` sheets (named to avoid the dashboard's GO filters).

- **RLE plot in QC.** Relative Log Expression plot added to the QC sheet, to
  document residual normalization quality after DIA-NN (diagnostic only, no
  re-normalization). Title reports the worst per-sample median deviation.

- **Utility scripts:** `check_paired.py` (pre-flight paired-design check),
  `check_gsea.py` (GSEA integration diagnostic), `regen_dashboard.py` (rebuild
  the HTML dashboard from an existing Excel workbook without rerunning analysis).

### Fixed

- **Perseverance: DIAMOND no longer hangs silently.** DIAMOND 2.2.8 (the version shipped)
  never exits when the first 10 sequences of a FASTA are nucleotide (CDS, transcripts):
  it reports "only contain DNA letters", keeps the message in its buffer when its output is
  captured, and waits forever at 0 % CPU — the run stopped at "Construction de la base
  DIAMOND (reference.fasta)..." with no message. Now (`rbh.py`):
  - both FASTA are checked before DIAMOND (`check_protein_fasta`): nucleotide content,
    UTF-16 or UTF-8 BOM encoding, CR-only line endings, empty sequences, not a FASTA,
    characters other than letters, `*` and `-` in a sequence line (DIAMOND 2.2.8 also hangs
    on a space or a `.`), and an interrupted UniProt download — the stream ends with
    "Error encountered when streaming data. Please try again later." (real case: mouse
    reference cut at 40,000 sequences) — clear message, Deimos continues without orthologues;
  - each DIAMOND call is watched: killed after 3 min without CPU time (`IDLE_KILL_S`) or
    5 s after an `Error:` line, with the last DIAMOND step and a hint in the message;
  - DIAMOND output goes to one log per step (`perseverance/diamond_<step>.log`, with the
    command line) instead of being discarded; blastp progress feeds a heartbeat line every
    30 s; Ctrl+C also stops DIAMOND (no orphan `diamond.exe`).
  Measured with DIAMOND 2.2.8: 6,703 vs 40,000 proteins, makedb < 1 s each, blastp ~5 s
  each (2 threads). Tests: `tests/test_rbh_watchdog.py`.
- **`tests/test_pca_loadings.py` restored**: the 06/11 push replaced it with an older copy,
  dropping `test_driver_labels_do_not_overlap` and `test_missing_gene_and_description_do_not_crash`.

- **GO lollipop / dotplot layout**: with a single source the general title overlapped the
  facet title (top margin as a fraction of a 2.4-inch figure); fixed 0.8-inch header now.
  Dotplot: the colour bar and the Count legend have their own columns (the legend sat on
  the colour bar), point areas scaled between the smallest and largest term (a 27-protein
  term made a 750 pt² disc), legend with 2-3 values spread from min to max.
- **STRING backend: Perseverance no longer proposed.** With `go_backend: string`
  (direct mode), the Perseverance question was asked after the statistics; its
  orthologues (UniProt accessions of another species) would then have been read as
  STRING identifiers. Perseverance belongs to the gProfiler route and is now skipped
  (`run_perseverance: true` reported as ignored).
- **STRING backend: low direct coverage -> orthology to the same species.** The direct
  mode was kept as soon as 30 % of the accessions were indexed in STRING (naked
  mole-rat, UniProt UP000694906: 16/50), so most proteins never reached the
  enrichment. Below 60 % (`STRING_DIRECT_MIN_COVERAGE`), when DIAMOND is available,
  Deimos now maps the proteins by RBH to the STRING proteome of the same species
  (identities ~100 %); without DIAMOND the direct mode is kept with a warning.

- **Methods sheet.** The GO section reported the measured species (`bnapus`) when the
  enrichment ran on the reference species through Perseverance orthologues; the row now
  names the reference organism and the translation (RBH, number of orthologues). Selected
  contrasts are no longer described as "all condition pairs"; the ANOVA correction is
  described as across proteins; the WGCNA network is described as unsigned
  (|cor|^beta, as computed).
- **ANOVA p-values** computed with the survival function (`1 - cdf` underflowed to 0 for
  very strong effects).
- **Volcanos.** Legend placed automatically again (it hid labels in the upper right
  corner); dashboard volcanos no longer floor p at 1e-10 (points piled up at
  -log10 p = 10) and their Y-axis ticks follow the range.
- **Orthologue coverage** capped at 100 % (alignment length with gaps could exceed the
  shorter sequence).

- **Perseverance with a UniProt reference FASTA.** DIAMOND reports the first word of the
  FASTA header, i.e. `sp|Q9C5R8|NAME_ARATH`, and these strings were sent as such to
  gProfiler, which does not recognise them (no enrichment). Ortholog identifiers are now
  reduced to the accession (`rbh.clean_seqid`), in the GO queries and in the
  `Orthologues_Perseverance` sheet; other header formats are unchanged.

- **Volcano Y axis back to the raw p-value.** With `volcano_use_padj: true` the individual
  volcanos plotted -log10(p.adj): without any differential protein, BH pushes every p.adj
  to ~1 and the whole cloud collapsed on the x axis (looked like a bug at early time
  points). The Y axis is now -log10(p.value), as limma's `volcanoplot()`; colours still
  follow the p.adj threshold and the dashed line sits at the raw p-value equivalent to it
  in each contrast (BH is monotone in p; no line when no protein passes). New key
  `volcano_y_axis: pvalue` (default) | `padj` (previous behaviour). The 1e-10 display floor
  is gone (p underflowing to 0 shown as triangles).
  - Facet volcano: the dotted line was drawn at p = threshold on the raw-p axis while
    colours used p.adj; now at the same equivalent raw p-value.
  - Scatter plots (Excel): status used a hard-coded raw p < 0.05; now the volcano rule
    (chosen statistic, threshold and ratio), as the UpSet and the robustness score.
  - Dashboard: the DEP status of the embedded points and the per-contrast up/down
    summaries follow the pipeline thresholds (they used raw p < 0.05 and |log2FC| > 0.26
    whatever the settings); the volcano line in p.adj mode and the comparative volcano
    lines sit at the raw p-value equivalent of each contrast; WGCNA hub rings use the
    scales of the drawn volcano (they were offset, and ignored the robustness filter and
    the LFC inversion). The interactive p.value exploration keeps all candidate points.

- **Sample-to-column matching (critical).** Replaced substring matching with
  exact matching. Previously, numeric-prefix labels collided (`Sple-2` matched
  `Sple-28`/`Sple-29`), which could silently assign a sample to the wrong run.
  **Users should re-check past analyses that had prefix-related sample labels.**

- **Robustness score under DEqMS.** The robustness loop now recomputes DEqMS
  (not limma) when DEqMS is the primary statistic, so the score reflects the same
  statistic shown in the `_p.val` / `_p.adj` columns.

- **Faster imputation robustness with the mixed imputation.** The kNN part of the
  mixed imputation does not depend on the random draw: it is now computed once and
  reused by the 100 re-imputations, only the QRILC part being redrawn. Results are
  identical (test `tests/test_knn_cache.py`); on 9,300 proteins × 180 samples a kNN
  call takes about 70 s against 0.6 s for the rest of an iteration, so the step goes
  from 100 kNN calls to one.

- **Scatter plots (dashboard).** Empty/deformed scatter fixed: intensity columns
  were not matched to contrasts for digit-prefixed conditions (missing `X`
  prefix). Axes are now square/symmetric, colors match the Excel scatter
  (Up=orange, Down=cyan), non-significant points drawn underneath, ±ratio lines
  added.

- **GO cross-contrast heatmap & network (dashboard).** (1) Empty rendering fixed:
  the term canvas exceeded the browser height limit (~32767 px) with thousands of
  terms — now capped at 150 terms (variance-ranked, per-source quota) and 150
  network nodes. (2) `REAC` (Reactome) and `KEGG` added to the source selectors
  and to the per-source panels. (3) Contrast→sheet mapping made deterministic
  (was confusing contrasts sharing a long common prefix).

- **Excel workbook corruption.** GO sheet names sharing a long common prefix
  truncated to identical 31-char names, corrupting `/xl/workbook.xml`. Names are
  now shortened (common affix stripped) and guaranteed unique.

- **UpSet plot labels.** Hard-coded `X180/X180` labels replaced by a generic
  common-prefix stripper (e.g. `Ctrl/di6h`).

- **GSEA parameter not applied.** The `run_gsea` flag chosen at launch was
  dropped in config resolution and never reached the pipeline; now transmitted
  end-to-end and persisted in `last_config.yaml`.

- **Duplicate prompt at launch.** The p-value question could appear twice due to
  a residual empty line in the input buffer (common on PowerShell). Empty inputs
  are now absorbed silently.

- **WGCNA `grey` module.** The unassigned-protein `grey` module is now shown in
  the module-size plot (hatched, labeled "unassigned", placed last) instead of
  being silently dropped — its size is a clustering-quality indicator.

### Changed

- **Full FR → EN translation** of the user-facing pipeline (prompts, progress
  messages, figure titles, Excel sheet names and headers, dashboard UI).
- **Dashboard header** homogenized with the sister suite (Phobos): icon +
  "Deimos" + expanded acronym + "Proteogen · Université de Caen".

### Dependencies

- **New:** `gseapy` (rank-based GSEA). See `requirements.txt`.
