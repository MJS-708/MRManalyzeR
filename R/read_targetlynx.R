#' Read a TargetLynx workbook into a wide sample x compound matrix
#'
#' Reads the TargetLynx sheet(s) with [extractTable()], pivots the long table
#' to a wide samples x compound matrix of the requested `datatype`, and
#' optionally masks values whose per-injection S/N is below `snr`. Columns are
#' the raw TargetLynx compound (`Processing_name`) identifiers; the
#' rename-to-`Compound` and feature/sample filtering happen later in
#' [process_dataset()].
#'
#' @param xlsx_path Path to the TargetLynx xlsx workbook.
#' @param datatype TargetLynx column to report (e.g. "Area", "Response",
#'   "ng/mL").
#' @param tl_headers TargetLynx column headers to extract.
#' @param snr S/N threshold; values below it are set to `NA`. `FALSE` to skip.
#' @param data_tab_names Sheet name(s) to read; `NULL` auto-matches sheets
#'   whose name contains "lcms_data".
#' @return A data frame, rows = samples, columns = compounds
#'   (`Processing_name`), numeric.
#' @examples
#' xlsx <- system.file("extdata", "example_data.xlsx", package = "MRManalyzeR")
#' m <- read_targetlynx(xlsx, datatype = "Area", snr = 3)
#' dim(m)
#' @family workflow steps
#' @export
read_targetlynx = function(xlsx_path, datatype = "Area",
                           tl_headers = c("ID", "Name", "Area", "ng/mL", "Response", "S/N"),
                           snr = FALSE, data_tab_names = NULL){

  lcms_table = extractTable(xlsx_path, tl_headers = tl_headers,
                            data_tab_names = data_tab_names) %>%
    subset(Name != "")

  out_table = .build_wide_matrix(lcms_table, datatype)

  if(!isFALSE(snr))
    out_table = .apply_snr_mask(out_table, lcms_table, snr)

  out_table
}
