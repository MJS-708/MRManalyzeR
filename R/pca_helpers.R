#' Pre-process a sample × feature matrix and run PCA
#'
#' Standard metabolomics-style preprocessing pipeline:
#' \enumerate{
#'   \item Drop all-NA features (always). Optionally apply a stricter
#'         feature missing-value filter (`feat_na_max`).
#'   \item Optionally drop samples with too many missing values (`sample_na_max`).
#'   \item Impute remaining NAs (per-feature minimum or half-minimum).
#'   \item Optional transform (log2(x+1) or sqrt(x)).
#'   \item Drop zero-variance features.
#'   \item `prcomp(center, scale.)` — mean-centre and (optionally) autoscale.
#' }
#'
#' Returns the prcomp object plus a human-readable trace of the steps that
#' were actually applied — reports embed this so the user can see at a
#' glance what transformation went into the PCA.
#'
#' @param X numeric matrix (rows = samples, cols = features). May contain NAs.
#' @param impute one of `"none"`, `"min"`, `"half_min"`, `"frac_min"` (default `"min"`).
#'   `"frac_min"` multiplies the per-feature minimum by `impute_frac` (e.g. 0.2).
#' @param impute_frac Numeric multiplier used when `impute = "frac_min"` (default 0.5).
#' @param transform one of `"none"`, `"log2"`, `"sqrt"` (default `"log2"`).
#' @param center logical, passed to [stats::prcomp()].
#' @param scale logical, passed to [stats::prcomp()] as `scale.`.
#' @param feat_na_max Maximum allowed fraction of NA values per feature (0–1).
#'   Features exceeding this threshold are dropped before imputation.
#'   Default `1` retains the original behaviour (only all-NA features are dropped).
#'   Set e.g. `0.5` to also drop features missing in more than 50 % of samples.
#' @param sample_na_max Maximum allowed fraction of NA values per sample (0–1).
#'   Samples exceeding this threshold are dropped before imputation.
#'   Default `1` applies no sample filtering (original behaviour).
#'   Set e.g. `0.8` to drop samples missing more than 80 % of features.
#' @return A list with elements
#'   \describe{
#'     \item{`pr`}{the `prcomp` object, or `NULL` if PCA could not be fit.}
#'     \item{`X_clean`}{the matrix actually fed to `prcomp`.}
#'     \item{`n_samples`, `n_features`}{dimensions after cleaning.}
#'     \item{`dropped_all_na`, `dropped_feat_filter`, `dropped_sample_filter`,
#'           `dropped_zero_var`}{counts of features/samples dropped at each step.}
#'     \item{`steps`}{character vector of human-readable steps applied.}
#'   }
#' @export
run_pca_pipeline = function(X,
                            impute        = c("min", "half_min", "frac_min", "none"),
                            impute_frac   = 0.5,
                            transform     = c("log2", "sqrt", "none"),
                            center        = TRUE,
                            scale         = TRUE,
                            feat_na_max   = 1,
                            sample_na_max = 1){

  impute    = match.arg(impute)
  transform = match.arg(transform)
  steps     = character(0)

  feat_na_max   = as.numeric(feat_na_max)
  sample_na_max = as.numeric(sample_na_max)

  X = as.matrix(X)

  # 1. Drop all-NA features (always) + optional stricter feature filter
  na_frac_feat = colMeans(is.na(X))
  all_na       = na_frac_feat >= 1
  dropped_all_na = sum(all_na)
  X = X[, !all_na, drop = FALSE]
  if(dropped_all_na > 0)
    steps = c(steps, sprintf("Dropped %d all-NA features.", dropped_all_na))

  dropped_feat_filter = 0L
  if(feat_na_max < 1){
    na_frac_feat2   = colMeans(is.na(X))
    drop_f          = na_frac_feat2 > feat_na_max
    dropped_feat_filter = sum(drop_f)
    X = X[, !drop_f, drop = FALSE]
    steps = c(steps, sprintf(
      "Feature missing filter (max %g%% NA): dropped %d feature(s).",
      feat_na_max * 100, dropped_feat_filter))
  } else {
    steps = c(steps, "No feature missing-value filter (feat_na_max = 1).")
  }

  # 2. Optional sample filter
  dropped_sample_filter = 0L
  if(sample_na_max < 1){
    na_frac_samp = rowMeans(is.na(X))
    drop_s       = na_frac_samp > sample_na_max
    dropped_sample_filter = sum(drop_s)
    X = X[!drop_s, , drop = FALSE]
    steps = c(steps, sprintf(
      "Sample missing filter (max %g%% NA): dropped %d sample(s).",
      sample_na_max * 100, dropped_sample_filter))
  } else {
    steps = c(steps, "No sample missing-value filter (sample_na_max = 1).")
  }

  # 3. Impute remaining NAs
  if(impute %in% c("min", "half_min", "frac_min")){
    fac = switch(impute,
                 min      = 1,
                 half_min = 0.5,
                 frac_min = as.numeric(impute_frac))
    for(j in seq_len(ncol(X))){
      if(any(is.na(X[, j]))){
        mn = suppressWarnings(min(X[, j], na.rm = TRUE))
        if(is.finite(mn)) X[is.na(X[, j]), j] = mn * fac
      }
    }
    fac_label = switch(impute,
                       min      = "minimum",
                       half_min = "half-minimum",
                       frac_min = sprintf("%g × minimum", as.numeric(impute_frac)))
    steps = c(steps, sprintf("Imputed NAs with per-feature %s.", fac_label))
  } else {
    steps = c(steps, "No NA imputation.")
  }

  # 3. Transform
  if(transform == "log2"){
    X[X < 0 & !is.na(X)] = 0
    X = log2(X + 1)
    steps = c(steps, "log2(x + 1) transform.")
  } else if(transform == "sqrt"){
    X[X < 0 & !is.na(X)] = 0
    X = sqrt(X)
    steps = c(steps, "sqrt(x) transform.")
  } else {
    steps = c(steps, "No transform.")
  }

  # 4. Drop zero-variance / still-NA features
  vars = apply(X, 2, stats::var, na.rm = TRUE)
  keep = is.finite(vars) & vars > 0
  dropped_zero_var = sum(!keep)
  X = X[, keep, drop = FALSE]
  if(dropped_zero_var > 0)
    steps = c(steps, sprintf("Dropped %d zero-variance features.", dropped_zero_var))

  # 5. PCA
  if(center) steps = c(steps, "Mean-centred per feature.")
  if(scale)  steps = c(steps, "Autoscaled (divide by per-feature SD).")

  pr = if(ncol(X) >= 2 && nrow(X) >= 3){
    tryCatch(stats::prcomp(X, center = center, scale. = scale),
             error = function(e){
               steps <<- c(steps, sprintf("prcomp() failed: %s.",
                                          conditionMessage(e)))
               NULL
             })
  } else {
    steps = c(steps, sprintf("Insufficient data: n_samples=%d, n_features=%d.",
                             nrow(X), ncol(X)))
    NULL
  }

  list(pr = pr, X_clean = X,
       n_samples             = nrow(X),
       n_features            = ncol(X),
       dropped_all_na        = dropped_all_na,
       dropped_feat_filter   = dropped_feat_filter,
       dropped_sample_filter = dropped_sample_filter,
       dropped_zero_var      = dropped_zero_var,
       steps                 = steps)
}
