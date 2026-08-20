#' Batch-correct a DatasetExperiment (median-ratio)
#'
#' Convenience wrapper that builds the median-ratio batch-correction model and
#' applies it, returning the corrected `DatasetExperiment`. Each batch is scaled
#' so its per-feature median matches the grand median across the reference
#' (`qc_label`) samples.
#'
#' The reference has to be material whose true value does not differ between
#' batches, because any difference that remains is treated as technical and
#' removed. Pooled QC injections satisfy that by construction - the same
#' extract, injected repeatedly - which is why `"QC"` is the default.
#'
#' Using the biological samples as the reference is a legitimate alternative,
#' and often a more stable one because there are more of them, but it rests on
#' the study being balanced across batches: it assumes the biology is the same
#' on average in each batch. Where that holds it is fine. Where batch is
#' confounded with the study factor - all treated animals run on day 2, say -
#' the between-batch difference is partly the effect you are looking for, and
#' correcting removes it silently. Pass `check_factor` to have that assumption
#' tested rather than assumed.
#'
#' @param de A `struct::DatasetExperiment`.
#' @param qc_label Label identifying reference samples in `factor_name`.
#' @param factor_name `sample_meta` column holding `qc_label`.
#' @param batch_head `sample_meta` column identifying the batch.
#' @param check_factor Optional `sample_meta` column holding the study factor.
#'   When given, the reference samples are cross-tabulated against it and a
#'   warning is issued if any batch carries only one level - the signature of
#'   a design where correcting on biological samples would remove real
#'   differences. `NULL` skips the check.
#' @param min_ref Minimum reference samples required in each batch. Below this
#'   the batch median is not estimable and the correction would be noise.
#' @return The batch-corrected `struct::DatasetExperiment`.
#' @examples
#' m  <- data.frame(A = c(10, 12, 20, 24), B = c(5, 6, 10, 12),
#'                  row.names = c("S1", "S2", "S3", "S4"))
#' fm <- data.frame(Compound = c("A", "B"), row.names = c("A", "B"))
#' sm <- data.frame(Sample_type = "Sample",
#'                  Chrom_Batch = c("b1", "b1", "b2", "b2"),
#'                  row.names = c("S1", "S2", "S3", "S4"))
#' de <- struct::DatasetExperiment(data = m, sample_meta = sm, variable_meta = fm)
#' correctBatch(de, qc_label = "Sample", factor_name = "Sample_type",
#'               batch_head = "Chrom_Batch")
#' @family peak-matrix processing
#' @export
correctBatch = function(de, qc_label = "QC", factor_name = "Sample_type",
                         batch_head = "Chrom_Batch",
                         check_factor = NULL, min_ref = 2){

  smeta = as.data.frame(de$sample_meta)
  for(h in c(factor_name, batch_head)){
    if(!h %in% colnames(smeta))
      stop(sprintf("[correctBatch] '%s' is not a sample_meta column.", h))
  }

  is_ref = smeta[[factor_name]] %in% qc_label
  batch  = as.character(smeta[[batch_head]])

  if(!any(is_ref))
    stop(sprintf(
      "[correctBatch] no reference samples: no row has %s in {%s}.",
      factor_name, paste(qc_label, collapse = ", ")))

  # Every batch needs enough reference samples for its median to mean
  # anything; a batch with one is scaled by a single observation.
  per_batch = table(batch[is_ref])
  thin = setdiff(unique(batch), names(per_batch)[per_batch >= min_ref])
  if(length(thin))
    stop(sprintf(
      "[correctBatch] batch(es) %s have fewer than %d reference samples (%s in {%s}). Add reference samples, lower min_ref, or use a reference present in every batch.",
      paste(sprintf("'%s'", thin), collapse = ", "), min_ref,
      factor_name, paste(qc_label, collapse = ", ")))

  # Is the "biology is the same on average in each batch" assumption safe?
  if(!is.null(check_factor)){
    if(!check_factor %in% colnames(smeta)){
      warning(sprintf("[correctBatch] check_factor '%s' is not a sample_meta column; confounding check skipped.",
                      check_factor))
    } else {
      lv = table(batch[is_ref],
                 as.character(smeta[[check_factor]][is_ref]))
      n_lv = rowSums(lv > 0)
      if(ncol(lv) > 1 && any(n_lv < 2))
        warning(sprintf(
          "[correctBatch] batch(es) %s contain only one level of '%s' among the reference samples, so batch is confounded with it. Median-ratio correction will remove that difference along with the technical shift. Use pooled QCs as the reference if you have them.",
          paste(sprintf("'%s'", names(n_lv)[n_lv < 2]), collapse = ", "),
          check_factor))
    }
  }

  bc = batch_correct(qc_label = qc_label, factor_name = factor_name,
                     batch_head = batch_head)
  model_apply(bc, de)@corrected@value
}
