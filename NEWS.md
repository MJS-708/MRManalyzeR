# MRManalyzeR 0.99.0

* Initial Bioconductor submission.
* End-to-end processing of targeted LC-MS lipidomics and metabolomics data
  exported from Waters TargetLynx, driven from a single YAML config via
  `run_MRManalyzeR()`: peak-matrix assembly, S/N and blank filtering,
  missing-value imputation, normalisation, concentration adjustment, and
  batch correction.
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
* A second bundled workbook, `example_batch_data.xlsx`, provides a two-batch
  Skyline export with per-sample protein amounts and reconstitution volumes, so
  normalisation, concentration adjustment and batch correction can be
  demonstrated.
