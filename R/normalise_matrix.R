#' Normalise a matrix by a per-sample metadata column
#'
#' Divides each sample (row) of `x` by that sample's value in
#' `sample_meta[[column]]`, matched by name (not row position, since the data
#' matrix and metadata row orders are not guaranteed to agree).
#'
#' @param x Wide numeric matrix / data frame, rows = samples.
#' @param sample_meta Sample metadata containing `column` and `name_col`.
#' @param column `sample_meta` column to divide by (e.g. protein amount or
#'   tissue weight).
#' @param name_col `sample_meta` column matching the row names of `x`.
#' @return `x` with each row divided by its per-sample divisor.
#' @examples
#' m  <- data.frame(A = c(10, 20), B = c(30, 40), row.names = c("S1", "S2"))
#' sm <- data.frame(Name = c("S1", "S2"), protein = c(2, 4))
#' normalise_matrix(m, sm, column = "protein")
#' @export
normalise_matrix = function(x, sample_meta, column, name_col = "Name"){
  idx = match(rownames(x), sample_meta[[name_col]])
  if(anyNA(idx))
    stop(sprintf(
      "[normalise_matrix] %d sample(s) have no matching row in sample_meta[['%s']]: %s",
      sum(is.na(idx)), name_col,
      paste(head(rownames(x)[is.na(idx)], 5), collapse = ", ")))
  x / sample_meta[[column]][idx]
}
