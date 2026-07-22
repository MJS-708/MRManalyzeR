#' Batch-correct a DatasetExperiment (median-ratio)
#'
#' Convenience wrapper that builds the median-ratio batch-correction model and
#' applies it, returning the corrected `DatasetExperiment`. Each batch is scaled
#' so its per-feature median matches the grand median across the reference
#' (`qc_label`) samples.
#'
#' @param de A `struct::DatasetExperiment`.
#' @param qc_label Label identifying reference samples in `factor_name`.
#' @param factor_name `sample_meta` column holding `qc_label`.
#' @param batch_head `sample_meta` column identifying the batch.
#' @return The batch-corrected `struct::DatasetExperiment`.
#' @examples
#' m  <- data.frame(A = c(10, 12, 20, 24), B = c(5, 6, 10, 12),
#'                  row.names = c("S1", "S2", "S3", "S4"))
#' fm <- data.frame(Compound = c("A", "B"), row.names = c("A", "B"))
#' sm <- data.frame(Sample_type = "Sample",
#'                  Chrom_Batch = c("b1", "b1", "b2", "b2"),
#'                  row.names = c("S1", "S2", "S3", "S4"))
#' de <- struct::DatasetExperiment(data = m, sample_meta = sm, variable_meta = fm)
#' correct_batch(de, qc_label = "Sample", factor_name = "Sample_type",
#'               batch_head = "Chrom_Batch")
#' @family peak-matrix processing
#' @export
correct_batch = function(de, qc_label = "QC", factor_name = "Sample_type",
                         batch_head = "Chrom_Batch"){
  bc = batch_correct(qc_label = qc_label, factor_name = factor_name,
                     batch_head = batch_head)
  model_apply(bc, de)@corrected@value
}
