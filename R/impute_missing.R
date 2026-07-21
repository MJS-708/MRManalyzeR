#' Impute missing values with a fraction of the per-feature minimum
#'
#' Zeros are treated as missing. For each feature, `NA`s among the non-blank
#' injections are replaced with `min(non-NA) * scalar` -- the common
#' "half-minimum" style of imputation for non-detects, which places them just
#' below the lowest value actually observed. Blank injections (identified from
#' `blank_head` in `sample_meta`) are left untouched.
#'
#' This is the only imputation applied to the reported matrix, and it is
#' deliberately simple: for a targeted panel a compound that is missing in most
#' samples is usually better dropped than modelled. Where the choice matters
#' more -- ahead of PCA -- [run_pca()] offers `"none"`, `"min"`, `"half_min"`
#' and `"frac_min"`, plus per-feature and per-sample missingness filters that
#' discard rather than fill.
#'
#' @param de A `struct::DatasetExperiment`.
#' @param scalar Numeric multiplier applied to the per-feature minimum, e.g.
#'   `0.5` for half-minimum or `0.2` for a fifth of the minimum.
#' @param blank_head `sample_meta` column identifying blank injections.
#' @param blank_name Value in `blank_head` marking a blank injection.
#' @return `de` with `data` imputed, in the original row order.
#' @family workflow steps
#' @examples
#' de <- struct::DatasetExperiment(
#'   data = data.frame(PGE2 = c(10, 5, NA), PGD2 = c(2, NA, 4),
#'                     row.names = c("S1", "S2", "S3")),
#'   sample_meta = data.frame(Sample_type = rep("Sample", 3),
#'                            row.names = c("S1", "S2", "S3")),
#'   variable_meta = data.frame(Compound = c("PGE2", "PGD2"),
#'                              row.names = c("PGE2", "PGD2")))
#' impute_missing(de, scalar = 0.5)$data
#' @export
impute_missing = function(de, scalar, blank_head = "Sample_type",
                          blank_name = "Blank"){

  if(isTRUE(scalar))
    stop("`scalar` must be numeric (e.g. 0.5), not TRUE.")

  df    = as.data.frame(de$data)
  smeta = as.data.frame(de$sample_meta)

  blank_samples = if(blank_head %in% colnames(smeta))
    rownames(smeta)[smeta[[blank_head]] %in% blank_name] else character(0)

  df[df == 0] = NA
  is_blank = rownames(df) %in% blank_samples

  # Per-feature minimum across the non-blank injections only.
  sample_df = df[!is_blank, , drop = FALSE]
  fill_vals = suppressWarnings(vapply(sample_df, min, numeric(1), na.rm = TRUE)) * scalar

  # Impute in place so the dataset's row order is preserved.
  for(j in which(is.finite(fill_vals))){
    col = df[[j]]
    col[is.na(col) & !is_blank] = fill_vals[j]
    df[[j]] = col
  }

  de$data = df
  de
}
