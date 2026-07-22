#' Blank-filter a dataset (mask values at or below the blank level)
#'
#' For each feature, computes `mean(blank injections) * blank_filter` and sets
#' any value below that threshold to `NA`. Blank injections are identified from
#' `blank_head` in the dataset's `sample_meta`, so nothing has to be passed
#' separately.
#'
#' @param de A `struct::DatasetExperiment`.
#' @param blank_filter Numeric fold-change multiplier applied to the blank mean
#'   (e.g. `3` keeps only peaks at least three times the blank level).
#' @param blank_head `sample_meta` column identifying blank injections.
#' @param blank_name Value in `blank_head` marking a blank injection.
#' @return `de` with sub-threshold values in `data` set to `NA`.
#' @family peak-matrix processing
#' @examples
#' de <- struct::DatasetExperiment(
#'   data = data.frame(PGE2 = c(100, 90, 5), PGD2 = c(50, 60, 1),
#'                     row.names = c("S1", "S2", "B1")),
#'   sample_meta = data.frame(Sample_type = c("Sample", "Sample", "Blank"),
#'                            row.names = c("S1", "S2", "B1")),
#'   variable_meta = data.frame(Compound = c("PGE2", "PGD2"),
#'                              row.names = c("PGE2", "PGD2")))
#' filter_blanks(de, blank_filter = 3)$data
#' @export
filter_blanks = function(de, blank_filter, blank_head = "Sample_type",
                         blank_name = "Blank"){

  df    = as.data.frame(de$data)
  smeta = as.data.frame(de$sample_meta)

  blank_samples = if(blank_head %in% colnames(smeta))
    rownames(smeta)[smeta[[blank_head]] %in% blank_name] else character(0)

  df_blank     = df[rownames(df) %in% blank_samples, , drop = FALSE]
  blank_thresh = colMeans(df_blank, na.rm = TRUE) * blank_filter

  # Features with no usable blank measurement give NaN and are left untouched.
  for(j in which(!is.na(blank_thresh))){
    col     = df[[j]]
    df[[j]] = ifelse(col < blank_thresh[j], NA, col)
  }

  de$data = df
  de
}
