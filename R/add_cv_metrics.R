#' Add per-feature CV metrics to a DatasetExperiment's variable_meta
#'
#' Computes the per-feature coefficient of variation (100 * sd / mean) over the
#' QC samples (`CV_QC`) and the biological samples (`CV_sample`), plus their
#' ratio (`CV_sample_vs_QC`), and adds all three as columns of `variable_meta`,
#' matched to each feature by row (Compound). QC / sample rows are identified
#' from `sample_type_head` in `sample_meta`.
#'
#' @param de A `struct::DatasetExperiment`.
#' @param sample_type_head `sample_meta` column classifying sample type.
#' @param qc_label Value(s) in `sample_type_head` marking QC samples.
#' @param sample_labels Value(s) marking biological samples.
#' @return `de` with `CV_QC`, `CV_sample` and `CV_sample_vs_QC` columns added to
#'   `variable_meta`.
#' @keywords internal
#' @noRd
.add_cv_metrics = function(de, sample_type_head = "Sample_type",
                           qc_label = "QC", sample_labels = "Sample"){
  dm    = as.data.frame(de$data)          # samples x features
  smeta = as.data.frame(de$sample_meta)
  vm    = as.data.frame(de$variable_meta)
  if(ncol(dm) == 0 || !sample_type_head %in% colnames(smeta)) return(de)

  cv_fn  = function(x) stats::sd(x, na.rm = TRUE) / mean(x, na.rm = TRUE) * 100
  cv_col = function(idx){
    if(length(idx) < 2){
      v = rep(NA_real_, ncol(dm)); names(v) = colnames(dm); return(v)
    }
    round(vapply(dm, function(x) cv_fn(x[idx]), numeric(1)), 1)
  }
  cv_qc     = cv_col(which(smeta[[sample_type_head]] %in% qc_label))
  cv_sample = cv_col(which(smeta[[sample_type_head]] %in% sample_labels))

  key = rownames(vm)                       # Compound == colnames(dm)
  vm$CV_QC                = unname(cv_qc[key])
  vm$CV_sample            = unname(cv_sample[key])
  # Named CV_sample_vs_QC, not "CV_sample/CV_QC": assigning into a
  # DatasetExperiment passes the frame through make.names(), which would
  # silently rewrite the "/" to "." anyway.
  vm$CV_sample_vs_QC = round(unname(cv_sample[key]) / unname(cv_qc[key]), 2)

  de$variable_meta = vm
  de
}
