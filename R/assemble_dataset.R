#' Assemble a struct::DatasetExperiment from a matrix + metadata
#'
#' Aligns `feature_meta` and `sample_meta` to the matrix (features to columns,
#' samples to rows, both by name) and builds the `DatasetExperiment` so
#' `variable_meta`/`sample_meta` rownames match the data assay exactly.
#'
#' @param x Wide numeric matrix / data frame, rows = samples, columns =
#'   compounds (canonical `Compound` names).
#' @param feature_meta Feature metadata; rows with `Report == "YES"` whose
#'   `Compound` is in `x` are kept, in matrix-column order.
#' @param sample_meta Sample metadata; rows whose `name_col` value is a row of
#'   `x` are kept.
#' @param name_col `sample_meta` column matching the matrix row names.
#' @return A `struct::DatasetExperiment`.
#' @examples
#' m  <- data.frame(A = c(10, 20, 30), B = c(1, 2, 3),
#'                  row.names = c("S1", "S2", "S3"))
#' fm <- data.frame(Compound = c("A", "B"), Report = "YES",
#'                  row.names = c("A", "B"))
#' sm <- data.frame(Name = c("S1", "S2", "S3"), Group = c("x", "y", "x"),
#'                  row.names = c("S1", "S2", "S3"))
#' assemble_dataset(m, fm, sm)
#' @family workflow steps
#' @export
assemble_dataset = function(x, feature_meta, sample_meta, name_col = "Name"){
  # Feature side: keep only feature_meta rows whose Compound is in the matrix,
  # in matrix-column order, so variable_meta rownames match assay colnames.
  fdata_output = feature_meta %>%
    dplyr::filter(Report == "YES", Compound %in% colnames(x)) %>%
    dplyr::arrange(match(Compound, colnames(x)))
  rownames(fdata_output) = fdata_output$Compound

  # Sample side: intersect first so dropped samples (blank filter, etc.)
  # don't get reintroduced as all-NA rows.
  metadata_ar = sample_meta[sample_meta[[name_col]] %in% rownames(x), , drop = FALSE]
  metadata_ar = metadata_ar[order(metadata_ar[[name_col]]), , drop = FALSE]
  rownames(metadata_ar) = metadata_ar[[name_col]]
  out_matrix_ar = x[metadata_ar[[name_col]], fdata_output$Compound, drop = FALSE]

  struct::DatasetExperiment(
    data          = out_matrix_ar,
    sample_meta   = metadata_ar,
    variable_meta = fdata_output
  )
}
