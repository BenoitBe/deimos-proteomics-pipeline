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
