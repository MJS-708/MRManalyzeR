#' Impute missing values with per-feature minimum * scalar
#'
#' Zeros are treated as NA. For each feature (column), NAs among non-blank
#' samples are replaced with `min(non-NA) * scalar`. Blank samples are
#' untouched and re-bound at the end.
#'
#' @param df Data frame or matrix of samples x features (rownames = sample names).
#' @param blank_samples Character vector of sample names to treat as blanks.
#' @param scalar Numeric multiplier applied to the feature minimum (e.g. 0.2).
#' @return Data frame of samples x features with NAs imputed.
#' @export
replaceNA = function(df, blank_samples, scalar){

  if(isTRUE(scalar)){
    stop("`scalar` must be numeric (e.g. 0.2), not TRUE.")
  }

  df[df == 0] = NA

  blank_df = df[rownames(df) %in% blank_samples, , drop = FALSE]
  sample_df = df[!rownames(df) %in% blank_samples, , drop = FALSE]

  # Per-feature minimum across non-blank samples
  min_vals = suppressWarnings(vapply(sample_df, min, numeric(1), na.rm = TRUE))
  fill_vals = min_vals * scalar

  # Vectorized imputation: only columns with a finite min
  for(j in which(is.finite(fill_vals))){
    col = sample_df[[j]]
    col[is.na(col)] = fill_vals[j]
    sample_df[[j]] = col
  }

  return(rbind(sample_df, blank_df))
}
