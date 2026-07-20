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
* `run_example()` runs the full pipeline on the bundled example dataset
  (a subset of Kolmert et al. 2018, doi:10.1016/j.prostaglandins.2018.05.005).
