#' MRManalyzeR: process and analyse targeted LC-MS lipidomics and metabolomics
#'
#' Processing and analysis pipeline for targeted lipidomics and metabolomics
#' data exported from Waters TargetLynx or Skyline, driven from a single YAML
#' config via [run_MRManalyzeR()].
#'
#' @keywords internal
#' @importFrom utils head
#' @importFrom stats median
#' @importFrom methods new
#' @importFrom struct model_apply
"_PACKAGE"

# dplyr NSE column references used across the readers / assemble step.
utils::globalVariables(c("ID", "ID2", "Name", "S/N", "Report", "Compound"))
