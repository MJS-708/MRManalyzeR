# MRManalyzeR 0.99.3

* **Function names are now camelCase**, following the Bioconductor style
  guide (`runMRManalyzeR()`, `processDataset()`, `plotHeatmap()`, and so on).
  This renames every exported function; there are no back-compatible aliases,
  so calling code must be updated. Acronyms stay upper case - `plotPCA()`,
  `addCVMetrics()`, `readTargetLynx()`.

* New `filterFeatures()`: drops features by detection frequency. `method =
  "within"` computes the detection fraction separately within each level of a
  grouping factor and keeps a feature that clears the threshold in any one of
  them, so a compound present in only one treatment survives; `"across"` pools
  every study sample and `"QC"` uses the QC injections. QC and blank
  injections never count towards the fraction, though they stay in the
  dataset. The fraction used is recorded per feature as `detected_frac`.

* New `filterSamples()`: drops study injections whose missing-value fraction
  exceeds a threshold, recording `na_frac` for every row. Only study samples
  are eligible - a QC over the threshold is reported but kept, since removing
  one can leave a batch with nothing for `correctBatch()` to anchor to.

* Both are off by default and are wired into `processDataset()` through new
  `filter_features:` and `filter_samples:` YAML blocks. They run once on the
  whole dataset, after blank masking and before normalisation and imputation:
  a detection fraction is a property of the study rather than of a batch, and
  imputation destroys the missingness both filters measure. Features are
  filtered before samples, so an injection is judged against features that
  are actually detectable.

* `processDataset()`'s per-batch helper was split accordingly, into a
  blank-filter pass and a normalise / concentration / impute pass, with the
  global filters in between. Per-batch imputation behaviour is unchanged.

# MRManalyzeR 0.99.2

* New `plotGroupHeatmap()`: a group-mean heatmap where the hue is the group
  mean and the opacity of each sample's sub-cell is that sample's contribution
  to it, so a mean resting on a single sample is visible rather than hidden.
  `plotHeatmap()` is unchanged and remains the right choice while the sample
  axis is still legible - it never averages.
* The colour limit defaults from the group means themselves rather than from a
  per-sample scale, since averaging replicates shrinks a scaled value by
  roughly sqrt(n) and an inherited limit leaves the panel washed out.

# MRManalyzeR 0.99.1

* Vignette and README reworked following review: repositioned around flexible
  quantitative outputs (raw areas through to calibrated concentrations) rather
  than internal-standard quantification alone, and around defined target lists
  rather than hypothesis testing alone.
* Vignette retitled, and validation promoted to its own stage in the workflow
  table and diagram.
* Corrected the description of `runStats()`: the test is named explicitly per
  comparison, not inferred from the design.
* Blanks are now excluded from the PCA examples, and the imputation fraction
  lowered to 0.2.
* "Ion ratios" renamed to "metabolite ratios" to avoid confusion with
  qualifier/quantifier transition ratios.
* Added Antonio Checa as an author, and a package logo.

# MRManalyzeR 0.99.0

* Initial Bioconductor submission.
* Post-acquisition processing of targeted LC-MS/MS lipidomics and metabolomics
  results exported from Waters TargetLynx or Skyline, or supplied as a plain
  sample-by-analyte matrix: peak-matrix assembly, S/N and LOD/LOQ filtering,
  blank filtering, missing-value imputation, normalisation, internal-standard
  and volume adjustment of vendor-reported concentrations, and batch
  correction. Chromatographic peak detection, integration and
  calibration-curve fitting are outside its scope.
* Every step is an exported function operating on a
  `struct::DatasetExperiment`; a whole workflow can additionally be driven
  from a single YAML config via `runMRManalyzeR()`.
* Self-contained HTML data-quality and statistical reports: group summaries,
  PCA, sample x feature heatmaps, volcano plots, feature-feature correlations,
  linear models, and enzyme-activity ion ratios.
* `runMRManalyzeRCombine()` merges multiple acquisition panels into a single
  analysis.
* `runExample()` runs the full workflow on the bundled example dataset
  (a subset of Kolmert et al. 2018, doi:10.1016/j.prostaglandins.2018.05.005).
* Every processing step is also an exported function operating on a
  `struct::DatasetExperiment`: `readTargetLynx()` / `readSkyline()`,
  `assembleDataset()`, `filterBlanks()`, `normaliseMatrix()`,
  `adjustConcentration()`, `imputeMissing()` and `correctBatch()`, composed
  by `processDataset()`.
* Per-compound quality metrics (`CV_QC`, `CV_sample` and their ratio) are
  written into the dataset's `variable_meta`.
* A second bundled dataset, `example_synthetic.xlsx` (plus the processed
  `example_synthetic.RDS`), is entirely simulated from a fixed seed: two
  chromatographic batches, per-sample protein amounts and reconstitution
  volumes, and a planted treatment effect on 5 of its 20 analytes, so
  normalisation, concentration adjustment, batch correction and the statistics
  can each be demonstrated against a known truth.
