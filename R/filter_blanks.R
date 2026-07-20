#' Blank-filter a matrix (mask values below the blank threshold)
#'
#' For each feature, computes `mean(blank samples) * blank_filter` and sets any
#' sample value below that threshold to `NA`.
#'
#' @param df Wide numeric matrix / data frame, rows = samples (rownames =
#'   sample names).
#' @param blank_samples Character vector of blank sample names.
#' @param blank_filter Numeric fold-change multiplier applied to the blank mean.
#' @return `df` with sub-threshold values set to `NA`.
#' @examples
#' m <- data.frame(A = c(100, 90, 5), B = c(50, 60, 1),
#'                 row.names = c("S1", "S2", "Blank1"))
#' filter_blanks(m, blank_samples = "Blank1", blank_filter = 3)
#' @export
filter_blanks = function(df, blank_samples, blank_filter){

  df_blank = df[rownames(df) %in% blank_samples, ]

  blank_df = data.frame(feature = names(df_blank),
                        blank_val = colMeans(df_blank, na.rm=TRUE))

  blank_df$blank_thresh = blank_df$blank_val * blank_filter

  for(col_idx in which(!is.na(blank_df$blank_thresh))){

    thresh = blank_df$blank_thresh[col_idx]

    temp_col = df[,col_idx]
    temp_col = ifelse(temp_col < thresh, NA, temp_col)

    df[,col_idx] = temp_col
  }

  return(df)
}
