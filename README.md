# MRManalyzeR
Tools to process and analyze targeted lipidomics and metabolomics data - based on the TargetLynx format from Waters. Primarily for internal use at Wheelock lab, Karolinska Institute, Sweden.

## Installation

Install the latest release from GitHub:

```r
# install.packages("remotes")
remotes::install_github("MJS-708/MRManalyzeR", ref = "main")
```

For the development version:

```r
remotes::install_github("MJS-708/MRManalyzeR", ref = "dev")
```

Local development install:

```r
devtools::install("path/to/MRManalyzeR")
```

## Quick start

Run the built-in example to verify your installation. The example uses a
published dataset of control and house dust mite challenged mice, with
several oxylipins (negative-ion mode):

```r
MRManalyzeR::run_example()
```

This renders a HTML report for data quality assessment
and for standard statistical analysis. These open in
directly in the RStudio viewer.
No file paths to configure.

## Usage

### Step 1 — Prepare the data .xlsx

The .xlsx file contain at least 3 sheets;

feature_metadata
- Contains all features in your analysis under `processing_name` which maps to `compound` the
name to report.
- Ensure Report header is populated with“YES” or “NO” for each compound based on manual
interpretation of the LC-MS/MS peak in TargetLynx (or similar). 

sample_metadata
- Containing all acquisition names from the TargetLynx file - they must be unchanged and match lcms_data_n
- Add any relavent study metadata for assessing data quality (e.g., injection order, extraction batch)
or statistical analysis (e.g., sex, treatment etc).

lcms_data_1
- This sheet must contain your complete summary from your Target Lynx file. 
- Additional lcms_data_n can be added, however all lcms_data sheets must have same features.


### Step 2 - Modify the YAML

Processing and analysis of the MRM data is configured via paramaeters in a YAML file
including (not exhaustive list):
- Paths and filenames to the data .xlsx curated in step 1 and where to save results
- PeakMatrixProcessing parameters including S/N thresholds, MV imputation,
batch correction, normalization
- Parameters for generating `data_quality_report` plots to assess CVs, signal drift and
correlation with study design (e.g., injection order).
- Parameters for generating `stats_report` plots to assess changes in data based upon 
study design (e.g., treatment) using PCA, boxplots etc.

A template YAML config is provided:

```r
file.copy(system.file("config_template.yaml", package = "MRManalyzeR"), "config.yaml")
```

### Step 3 — Run the workflow via YAML

The data can be processed and analyzed with a single command:

```r
res <- MRManalyzeR::run_MRManalyzeR("path/to/config.yaml")
```

This will:
1. Load and process all data defined in the YAML - saving a results.xlsx with a processed matrix
and statistics tables.
2. Generate the `data_quality_report` and `stats_report` (as .html)
