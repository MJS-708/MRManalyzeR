#' Extract and reshape TargetLynx tables from a Waters-exported xlsx
#'
#' Reads every sheet matching `lcms_data*`, locates the compound blocks
#' (rows starting with "Compound"), and returns a long data frame with one
#' row per sample x compound. The `ID` column is overwritten with the
#' extracted compound name.
#'
#' @param lcms_wb Either a path to a TargetLynx xlsx file, or an `openxlsx`
#'   workbook object (e.g. from `openxlsx::loadWorkbook()`). A path is
#'   preferred — it avoids a C-level bug in `loadWorkbook()` that surfaces
#'   on workbooks with certain styling/drawing XML.
#' @param tl_headers Character vector of TargetLynx column headers to keep
#'   (e.g. `c("ID", "Name", "Area", "ng/mL", "Response", "S/N")`).
#' @return A data frame of stacked compound blocks with columns matching `tl_headers`.
#' @export
extractTable <- function(lcms_wb, tl_headers) {

  # Accept either a path or a Workbook object
  is_path = is.character(lcms_wb) && length(lcms_wb) == 1
  all_sheets = if(is_path) openxlsx::getSheetNames(lcms_wb) else lcms_wb[[".->sheet_names"]]

  lcms_sheetNames = all_sheets[grep("lcms_data", all_sheets)]

  out_list = vector("list", length(lcms_sheetNames))

  for(s in seq_along(lcms_sheetNames)){

    lcms_df = openxlsx::read.xlsx(lcms_wb, sheet = lcms_sheetNames[s])

    lcms = lcms_df %>%
      janitor::row_to_names(row_number = 3,
                            remove_row = FALSE, remove_rows_above = FALSE)
    # Drop any pre-existing column literally named "ID" (e.g. a metadata
    # column whose row-3 header reads "ID") so the position-based rename
    # of column 1 to "ID" can't produce duplicate column names. NA column
    # names (from blank row-3 cells) are kept; only an exact "ID" match drops.
    .nm   = as.character(colnames(lcms))
    .keep = is.na(.nm) | .nm != "ID"
    lcms  = lcms[, .keep, drop = FALSE]
    lcms = lcms %>%
      dplyr::rename("ID" = 1, "ID2" = 2) %>%
      dplyr::mutate(ID = paste(ID, ID2, sep = ", ")) %>%
      dplyr::select(dplyr::all_of(tl_headers))

    compoundsindex = grep("Compound", lcms$ID)

    if(length(compoundsindex) < 2){
      warning(sprintf("Sheet '%s': fewer than 2 compound blocks found; skipped.",
                      lcms_sheetNames[s]))
      next
    }

    # Blocks are evenly spaced; subtract 2 header rows to get sample count
    injections = compoundsindex[2] - compoundsindex[1] - 2

    # Vectorized compound-name extraction (parse "Compound: <name>" row)
    comp_names = vapply(compoundsindex, function(i){
      row_vals = unlist(lcms[i, ])
      row_vals = row_vals[!is.na(row_vals)]
      nm = paste(unlist(strsplit(row_vals, split = ":  "))[-1], collapse = ", ")
      gsub(", NA", "", nm)
    }, character(1))

    # Build a single index vector covering all sample rows across blocks
    starts  = compoundsindex + 2
    row_idx = unlist(lapply(starts, function(k) k:(k + injections - 1)))

    block = lcms[row_idx, , drop = FALSE]
    block$ID = rep(comp_names, each = injections)

    out_list[[s]] = block
  }

  dplyr::bind_rows(out_list)
}
