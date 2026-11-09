# MRManalyzeR 0.99.10

First Bioconductor release, in Bioconductor 3.24.

MRManalyzeR processes the quantitative results of targeted LC-MS/MS
lipidomics and metabolomics assays, as exported from Waters TargetLynx or
Skyline, or supplied as a plain sample-by-analyte matrix. Every step is an
exported function acting on a `struct::DatasetExperiment`, and a whole study
can be run from one YAML configuration.

## Reading data

* `readTargetLynx()`, `readSkyline()` and `readPeakMatrix()` read the three
  sources; `assembleDataset()` joins the matrix to the sample and feature
  metadata. A zero is read as a non-detect.

## Processing

* Signal filtering on S/N, LOD or LOQ; `filterBlanks()` against extraction
  blanks; `filterFeatures()` by detection frequency, within groups, across
  samples or in the QCs; `filterSamples()` by missing-value fraction.
* `correctBatch()` aligns batches on the median of each compound in the
  reference (QC) injections, optionally for a chosen set of compounds only.
* `normaliseMatrix()` divides by a sample amount, and `adjustConcentration()`
  applies internal-standard and volume adjustment to vendor-reported
  concentrations.
* `imputeMissing()` fills gaps with random draws below each compound's lowest
  measurement, per processing batch and from a fixed seed. The measured values
  are kept alongside, and the statistics use them unless told otherwise.
* `processDataset()` runs the steps in a fixed order and records, for every
  compound removed, which step removed it.

## Quality control

* `addCVMetrics()` adds `CV_QC`, `CV_sample` and their ratio to the feature
  metadata; `plotDrift()` and `plotQQ()` show each compound against injection
  order and against normality.
* `mrmValidateInput()`, `mrmValidateConfig()` and `mrmValidateDesign()` check
  the workbook, the configuration and the design before anything is
  processed, and name what is wrong.

## Statistics and figures

* `runStats()` runs the comparisons with the test named for each one (Welch,
  Student or Wilcoxon for two groups; ANOVA or Kruskal-Wallis with post-hoc
  tests for more), adjusted by Benjamini-Hochberg or Bonferroni, along with
  Pearson and Spearman correlations, linear models and metabolite ratios.
* `mrmRunPCA()`, `mrmPlotPCA()`, `plotLoadings()`, `mrmPlotHeatmap()`,
  `plotGroupHeatmap()`, `plotBoxplot()`, `plotVolcano()` and
  `plotCorrelation()` draw every figure in the reports.
* Box and bar plots can split into panels with their own y-axes
  (`plotBoxplot(free_y = TRUE)`; `boxplot: facet_by:` and `free_y:` in the
  configuration), for groups measured on different scales such as tissue and
  culture media. A point above the axis cap is drawn at it as a triangle,
  with its value in the tooltip, rather than left off the plot.

## YAML workflow and reports

* `newConfig()` writes a configuration to edit; `runMRManalyzeR()` runs it,
  writing the processed matrix and statistics as xlsx and RDS, and
  self-contained HTML data-quality and statistics reports.
* `runMRManalyzeRCombine()` merges separately acquired panels by sample and
  runs the statistics across them.
* `mrmRunExample()` runs the whole workflow on the bundled data: a subset of
  the BAL fluid oxylipin data of Kolmert *et al.* (2018), and a simulated
  two-batch dataset with a known treatment effect.

## Before this release

Changes made during development and review (0.99.0 to 0.99.8) are listed in
full at <https://github.com/MJS-708/MRManalyzeR/blob/1a2b9d0/NEWS.md>. Scripts
written against a pre-release version from GitHub need updating for two
renames: function names became camelCase in 0.99.3, and seven functions took
an `mrm` prefix in 0.99.6.
