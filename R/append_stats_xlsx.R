#' Append/replace stats sheets in the results xlsx
#'
#' Adds (or overwrites) the following tabs in the existing results xlsx:
#' \itemize{
#'   \item `stats`         — long table of comparisons + post-hoc.
#'   \item `correlations`  — long table of feature × feature correlations.
#'   \item `linear_models` — long table of `lm()` coefficients.
#'   \item `cor_<name>_<subset>_<method>` — one wide correlation matrix tab
#'         per (correlation × subset × method). Sheet names are sanitised
#'         and truncated to 31 characters (Excel limit).
#' }
#' If the workbook does not exist (e.g. results dir was cleaned), a new one
#' is written containing only the stats tabs.
#'
#' @param out_xlsx Path to the results xlsx (from `run_MRManalyzeR()`).
#' @param stats_tables A list with elements `stats`, `correlations`,
#'   `linear_models` — as returned by [`run_stats()`].
#' @return Invisibly, the path written.
#' @export
append_stats_xlsx = function(out_xlsx, stats_tables){

  if(file.exists(out_xlsx)){
    wb = openxlsx::loadWorkbook(out_xlsx)
  } else {
    wb = openxlsx::createWorkbook()
  }

  # 1. Long-format tabs
  for(sheet_name in c("stats", "correlations", "linear_models")){
    df = stats_tables[[sheet_name]]
    if(is.null(df)) next
    .write_or_replace(wb, sheet_name, df)
  }

  # 1b. Per-comparison stats tabs
  if(!is.null(stats_tables$stats) && nrow(stats_tables$stats) > 0){
    for(cn in unique(stats_tables$stats$comparison)){
      sub = stats_tables$stats[stats_tables$stats$comparison == cn, , drop = FALSE]
      .write_or_replace(wb, .safe_sheet_name(sprintf("stats_%s", cn)), sub)
    }
  }

  # 1c. Per-model linear model tabs
  if(!is.null(stats_tables$linear_models) && nrow(stats_tables$linear_models) > 0){
    for(mn in unique(stats_tables$linear_models$model)){
      sub = stats_tables$linear_models[stats_tables$linear_models$model == mn, , drop = FALSE]
      .write_or_replace(wb, .safe_sheet_name(sprintf("lm_%s", mn)), sub)
    }
  }

  # 2. Wide-matrix tabs for correlations (one per subset × method × correlation)
  if(!is.null(stats_tables$correlations) && nrow(stats_tables$correlations) > 0){
    .write_corr_matrix_tabs(wb, stats_tables$correlations)
  }

  openxlsx::saveWorkbook(wb, file = out_xlsx, overwrite = TRUE)
  invisible(out_xlsx)
}


# ---- Internal helpers ------------------------------------------------------

#' @keywords internal
#' @noRd
.write_or_replace = function(wb, sheet_name, df){
  if(sheet_name %in% openxlsx::sheets(wb))
    openxlsx::removeWorksheet(wb, sheet_name)
  openxlsx::addWorksheet(wb, sheet_name)
  openxlsx::writeData(wb, sheet_name, df,
                      colNames = TRUE, rowNames = FALSE, keepNA = FALSE)
}

#' Write one wide correlation matrix per (correlation × subset × method)
#' @keywords internal
#' @noRd
.write_corr_matrix_tabs = function(wb, corr_long){
  combos = unique(corr_long[, c("correlation", "subset", "method")])

  for(i in seq_len(nrow(combos))){
    cn = combos$correlation[i]
    sb = combos$subset[i]
    mt = combos$method[i]

    sub = corr_long[corr_long$correlation == cn &
                    corr_long$subset      == sb &
                    corr_long$method      == mt, , drop = FALSE]
    if(nrow(sub) == 0) next

    # Preserve the import order produced by .pairwise_cor_long() (column-major
    # walk over the lower triangle => first-seen order == colnames(X) order).
    feats = unique(c(sub$feature_a, sub$feature_b))
    M = matrix(NA_real_, nrow = length(feats), ncol = length(feats),
               dimnames = list(feats, feats))
    for(k in seq_len(nrow(sub))){
      a = sub$feature_a[k]; b = sub$feature_b[k]; r = sub$estimate[k]
      M[a, b] = r; M[b, a] = r
    }
    # Self-self pairs are dropped from `correlations` (long table); leave the
    # diagonal as NA in the wide matrix so the empty cells make that obvious.
    # If the long table did include self-pairs (when include_self: True), the
    # loop above already filled the diagonal.

    sheet = .safe_sheet_name(sprintf("cor_%s_%s_%s", cn, sb, mt))
    .write_or_replace(wb, sheet, as.data.frame(M, check.names = FALSE))
    # Re-write so row names appear (writeData with rowNames=TRUE)
    openxlsx::removeWorksheet(wb, sheet)
    openxlsx::addWorksheet(wb, sheet)
    openxlsx::writeData(wb, sheet, as.data.frame(M, check.names = FALSE),
                        colNames = TRUE, rowNames = TRUE, keepNA = FALSE)
  }
}

#' Sanitise and truncate to Excel's 31-char sheet name limit.
#' @keywords internal
#' @noRd
.safe_sheet_name = function(x){
  x = gsub("[\\\\/?*\\[\\]:]", "_", x, perl = TRUE)
  x = gsub("[^A-Za-z0-9_]+", "_", x)
  x = gsub("_+", "_", x)
  x = sub("^_|_$", "", x)
  if(nchar(x) > 31) x = substr(x, 1, 31)
  if(!nzchar(x))    x = "sheet"
  x
}
