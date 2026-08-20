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
#' more -- ahead of PCA -- [runPCA()] offers `"none"`, `"min"`, `"half_min"`
#' and `"frac_min"`, plus per-feature and per-sample missingness filters that
#' discard rather than fill.
#'
#' @param de A `struct::DatasetExperiment`.
#' @param scalar Numeric multiplier applied to the per-feature minimum, e.g.
#'   `0.5` for half-minimum or `0.2` for a fifth of the minimum.
#' @param blank_head `sample_meta` column identifying blank injections.
#' @param blank_name Value in `blank_head` marking a blank injection.
#' Imputation is the one step that invents numbers, so it records what it did:
#' `impute_fill` (the value used for each feature), `n_imputed` and
#' `frac_imputed` are added to `variable_meta`. A compound whose
#' `frac_imputed` is high is not a measured compound - it is mostly a constant
#' - and a difference found in it is an artefact of the fill value. Check that
#' column before interpreting anything, and consider excluding above
#' `warn_frac`.
#'
#' @param warn_frac Warn about features imputed in more than this fraction of
#'   the non-blank injections. `NULL` disables the warning.
#' @return `de` with `data` imputed in the original row order, and
#'   `impute_fill` / `n_imputed` / `frac_imputed` added to `variable_meta`.
#' @family peak-matrix processing
#' @examples
#' de <- struct::DatasetExperiment(
#'   data = data.frame(PGE2 = c(10, 5, NA), PGD2 = c(2, NA, 4),
#'                     row.names = c("S1", "S2", "S3")),
#'   sample_meta = data.frame(Sample_type = rep("Sample", 3),
#'                            row.names = c("S1", "S2", "S3")),
#'   variable_meta = data.frame(Compound = c("PGE2", "PGD2"),
#'                              row.names = c("PGE2", "PGD2")))
#' imputeMissing(de, scalar = 0.5)$data
#' @export
imputeMissing = function(de, scalar, blank_head = "Sample_type",
                          blank_name = "Blank", warn_frac = 0.5){

  if(isTRUE(scalar))
    stop("`scalar` must be numeric (e.g. 0.5), not TRUE.")

  df    = as.data.frame(de$data)
  smeta = as.data.frame(de$sample_meta)
  vm    = as.data.frame(de$variable_meta)

  blank_samples = if(blank_head %in% colnames(smeta))
    rownames(smeta)[smeta[[blank_head]] %in% blank_name] else character(0)

  df[df == 0] = NA
  is_blank = rownames(df) %in% blank_samples

  # Per-feature minimum across the non-blank injections only.
  sample_df = df[!is_blank, , drop = FALSE]
  fill_vals = suppressWarnings(vapply(sample_df, min, numeric(1), na.rm = TRUE)) * scalar

  # Impute in place so the dataset's row order is preserved.
  n_imp = stats::setNames(integer(ncol(df)), colnames(df))
  for(j in which(is.finite(fill_vals))){
    col      = df[[j]]
    hit      = is.na(col) & !is_blank
    n_imp[j] = sum(hit)
    col[hit] = fill_vals[j]
    df[[j]]  = col
  }

  n_sample = sum(!is_blank)
  frac_imp = if(n_sample > 0) n_imp / n_sample else n_imp * NA_real_

  key = rownames(vm)
  vm$impute_fill  = unname(fill_vals[key])
  vm$n_imputed    = unname(n_imp[key])
  vm$frac_imputed = round(unname(frac_imp[key]), 3)

  if(!is.null(warn_frac)){
    bad = names(frac_imp)[is.finite(frac_imp) & frac_imp > warn_frac]
    if(length(bad))
      warning(sprintf(
        "[imputeMissing] %d feature(s) imputed in more than %.0f%% of samples and are mostly a constant: %s",
        length(bad), warn_frac * 100,
        paste(utils::head(bad, 5), collapse = ", ")))
  }

  de$data          = df
  de$variable_meta = vm
  de
}
