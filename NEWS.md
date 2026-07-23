# MRManalyzeR 0.99.1

* Vignette and README reworked following review: repositioned around flexible
  quantitative outputs (raw areas through to calibrated concentrations) rather
  than internal-standard quantification alone, and around defined target lists
  rather than hypothesis testing alone.
* Vignette retitled, and validation promoted to its own stage in the workflow
  table and diagram.
* Corrected the description of `run_stats()`: the test is named explicitly per
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
  from a single YAML config via `run_MRManalyzeR()`.
* Self-contained HTML data-quality and statistical reports: group summaries,
  PCA, sample x feature heatmaps, volcano plots, feature-feature correlations,
  linear models, and enzyme-activity ion ratios.
* `run_MRManalyzeR_combine()` merges multiple acquisition panels into a single
  analysis.
* `run_example()` runs the full workflow on the bundled example dataset
  (a subset of Kolmert et al. 2018, doi:10.1016/j.prostaglandins.2018.05.005).
* Every processing step is also an exported function operating on a
  `struct::DatasetExperiment`: `read_targetlynx()` / `read_skyline()`,
  `assemble_dataset()`, `filter_blanks()`, `normalise_matrix()`,
  `adjust_concentration()`, `impute_missing()` and `correct_batch()`, composed
  by `process_dataset()`.
* Per-compound quality metrics (`CV_QC`, `CV_sample` and their ratio) are
  written into the dataset's `variable_meta`.
* A second bundled dataset, `example_synthetic.xlsx` (plus the processed
  `example_synthetic.RDS`), is entirely simulated from a fixed seed: two
  chromatographic batches, per-sample protein amounts and reconstitution
  volumes, and a planted treatment effect on 5 of its 20 analytes, so
  normalisation, concentration adjustment, batch correction and the statistics
  can each be demonstrated against a known truth.
