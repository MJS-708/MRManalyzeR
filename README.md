# MRManalyzeR

Process and analyse targeted lipidomics / metabolomics data exported from
Waters TargetLynx. Reads the TargetLynx xlsx workbook, builds a peak-area
or concentration matrix (with optional SNR filtering, blank filtering,
missing-value imputation, normalisation and batch correction) and renders
self-contained HTML reports for data quality and statistical inference.
Driven from a single YAML config — no driver script needed.

Primarily for internal use at the Wheelock lab, Karolinska Institutet
(Sweden).

---

## Installation

### Bioconductor dependencies (do this first)

`MRManalyzeR` depends on three Bioconductor packages (`struct`,
`structToolbox`, `pmp`) that aren't on CRAN, so `install_github()` /
`devtools::install()` cannot pull them automatically. Install them once:

```r
if(!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(c("struct", "structToolbox", "pmp"))
```

### Latest release

```r
# install.packages("remotes")
remotes::install_github("MJS-708/MRManalyzeR", ref = "main")
```

Development branch:

```r
remotes::install_github("MJS-708/MRManalyzeR", ref = "dev")
```

Local install (from a cloned working copy):

```r
devtools::install("path/to/MRManalyzeR")
```

---

## Quick start — bundled example

Run the pipeline against a small bundled oxylipin dataset to verify the
install and produce both reports. No paths to configure:

```r
res <- MRManalyzeR::run_example()
attr(res, "results_dir")    # where the xlsx + html reports were written
```

`run_example()` copies the bundled YAML, points it at
`inst/extdata/example_data.xlsx` and a fresh `tempdir()`, runs
`run_MRManalyzeR()`, and (in interactive sessions) opens the rendered
HTML reports in the RStudio viewer / default browser.

### Example dataset

A subset of the BAL fluid oxylipin panel from **Kolmert *et al.* (2018),
*Prostaglandins & Other Lipid Mediators* 137, 11–18.**
DOI: <https://doi.org/10.1016/j.prostaglandins.2018.05.005>.

The bundled xlsx (`inst/extdata/example_data.xlsx`) contains a
TargetLynx-style workbook with `feature_metadata`, `sample_metadata` and
`lcms_data_*` sheets; the bundled YAML (`inst/extdata/example_config.yml`)
is annotated and ready to copy as a starting point for your own studies.

---

## Usage on your own data

### Step 1 — prepare the data xlsx

The input workbook needs at least three sheets:

* **`feature_metadata`** — one row per compound. Required columns:
  `Processing_name` (TargetLynx ID, used internally), `Compound`
  (display name), and `Report` populated with `YES` / `NO` based on
  manual peak interpretation.
* **`sample_metadata`** — one row per acquisition. The acquisition name
  must match the TargetLynx export exactly. Add any study metadata you
  want available for QC plots and statistics (e.g. `Injection_order`,
  `Chrom_Batch`, `Sex`, `Treatment`, `Sample_type`).
* **`lcms_data_1`** (and optionally `lcms_data_2`, …) — the TargetLynx
  summary table. All `lcms_data_*` sheets must share the same feature set.

### Step 2 — copy and edit the YAML template

```r
file.copy(
  system.file("extdata", "example_config.yml", package = "MRManalyzeR"),
  "config.yml"
)
```

Then edit `config.yml`. The key blocks are:

| Block | Purpose |
|---|---|
| `paths:` | input xlsx location, output directory, datatype to report (`Area` / `Response` / `ng/mL` / …) |
| `PeakMatrixProcessing:` | SNR threshold, blank/MV filters, batch correction, normalisation |
| `data_quality_report:` | per-compound intensity-vs-injection-order plots, QQ plots, QC PCA |
| `stats_report:` | comparisons (auto t-test/Wilcoxon/ANOVA), boxplots, correlations, linear models, PCA, volcano, ion ratios, heatmap |

Every option is documented inline in the bundled YAML — read it once
before editing.

### Step 3 — run

```r
res <- MRManalyzeR::run_MRManalyzeR("config.yml")
```

This writes:

* `<fn>_<datatype>_.xlsx` — processed peak matrix + the resolved YAML
  parameters
* `<fn>_<datatype>__stats.xlsx` — `stats`, `correlations`, `linear_models`
  tabs from the comparisons defined in the YAML
* `<fn>_<datatype>__data_quality_report.html`
* `<fn>_<datatype>__stats_report.html`

To re-render reports from an already-processed RDS (skipping the xlsx
ingest), set `PeakMatrixProcessing.execute: False` and re-run.

---

## Combining multiple acquisitions

To merge separately-acquired panels (e.g. GOM, cysLT, SPM, or PL-pos /
PL-neg / SL) sharing biological samples, use the combine workflow:

```r
file.copy(
  system.file("extdata", "example_combine_config.yml", package = "MRManalyzeR"),
  "combine_config.yml"
)
res <- MRManalyzeR::run_MRManalyzeR_combine("combine_config.yml")
```

The combine YAML lists per-panel paths (`.RDS` or `.xlsx`), optional QC
remappings (e.g. SL `QC1..5` ↔ PL `QC1, QC2, QC5, QC7, QC8`), feature
prefixing, and the same `stats_report:` block as a regular run.
`PeakMatrixProcessing` and the data quality report are skipped — the
inputs are already-processed matrices.

---

## Reports

Both reports are self-contained HTML with tabset navigation and (by
default) plotly hover tooltips that surface `Sample_ID` on points and
boxplot jitter. Toggle off with `interactive_plots: False` in either
report block of the YAML.

| Report | Sections |
|---|---|
| `data_quality_report` | sample summary; per-compound intensity vs injection order; per-compound Q-Q normality; QC PCA (one tab per coloring) |
| `stats_report` | global summaries; per-comparison p-value tables; per-compound boxplots (faceted across comparisons + by free-form factor); per-comparison PCA (scores + loadings); volcano; ion ratios (e.g. EpOME/DiHOME for sEH activity); samples × features heatmap; correlations; linear models |

---

## License

See `LICENSE`.
