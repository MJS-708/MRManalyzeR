#' MRManalyzeR: targeted LC-MS lipidomics and metabolomics processing
#'
#' Processes and analyses targeted lipidomics and metabolomics data exported from
#' Waters TargetLynx or Skyline: builds a peak-area / concentration matrix with
#' signal-to-noise or LOD/LOQ filtering, blank filtering, normalisation,
#' calibration/internal-standard concentration adjustment, missing-value
#' imputation and batch correction, then runs the statistical analyses and
#' renders two self-contained HTML reports.
#'
#' @section Entry points:
#' [run_MRManalyzeR()] drives the whole workflow from a single YAML config;
#' [run_MRManalyzeR_combine()] merges several acquisition panels into one
#' analysis; [run_example()] runs the bundled example dataset end to end.
#'
#' @section Workflow steps:
#' [read_targetlynx()] and [read_skyline()] read a workbook into a wide
#' sample x compound matrix, which [assemble_dataset()] turns into a
#' `struct::DatasetExperiment`. Every later step operates on that dataset:
#' [filter_blanks()], [normalise_matrix()], [adjust_concentration()],
#' [impute_missing()] and [correct_batch()]. [process_dataset()] composes them
#' all, and [run_MRManalyzeR()] drives the whole workflow from a YAML config.
#'
#' @section Analysis and utilities:
#' [run_stats()], [run_pca()], [subset_dataset()], [combine_datasets()],
#' [load_config()] and [load_dataset()].
#'
#' See `vignette("MRManalyzeR")` for a step-by-step walkthrough.
#'
#' @keywords internal
#' @importFrom utils head
#' @importFrom stats median
#' @importFrom methods new
#' @importFrom struct model_apply
"_PACKAGE"

# dplyr NSE column references used across the readers / assemble step.
utils::globalVariables(c("ID", "ID2", "Name", "S/N", "Report", "Compound"))
