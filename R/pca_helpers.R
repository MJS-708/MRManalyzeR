#' Pre-process a sample × feature matrix and run PCA
#'
#' Standard metabolomics-style preprocessing pipeline:
#' \enumerate{
#'   \item Drop all-NA features.
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
#' @return A list with elements
#'   \describe{
#'     \item{`pr`}{the `prcomp` object, or `NULL` if PCA could not be fit.}
#'     \item{`X_clean`}{the matrix actually fed to `prcomp`.}
#'     \item{`n_samples`, `n_features`}{dimensions after cleaning.}
#'     \item{`dropped_all_na`, `dropped_zero_var`}{counts of features dropped.}
#'     \item{`steps`}{character vector of human-readable steps applied.}
#'   }
#' @export
run_pca_pipeline = function(X,
                            impute      = c("min", "half_min", "frac_min", "none"),
                            impute_frac = 0.5,
                            transform   = c("log2", "sqrt", "none"),
                            center      = TRUE,
                            scale       = TRUE){

  impute    = match.arg(impute)
  transform = match.arg(transform)
  steps     = character(0)

  X = as.matrix(X)

  # 1. Drop all-NA features
  all_na = colSums(!is.na(X)) == 0
  dropped_all_na = sum(all_na)
  X = X[, !all_na, drop = FALSE]
  if(dropped_all_na > 0)
    steps = c(steps, sprintf("Dropped %d all-NA features.", dropped_all_na))

  # 2. Impute remaining NAs
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
       n_samples         = nrow(X),
       n_features        = ncol(X),
       dropped_all_na    = dropped_all_na,
       dropped_zero_var  = dropped_zero_var,
       steps             = steps)
}
