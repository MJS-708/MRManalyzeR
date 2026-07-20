#' Pipe operator
#'
#' Re-exports the magrittr pipe so it's available inside package code.
#'
#' @return The value returned by the right-hand side expression, applied to
#'   the left-hand side value.
#'
#' @importFrom magrittr %>%
#' @name %>%
#' @rdname pipe
#' @examples
#' c(1, 4, 9) %>% sqrt()
#' @export
NULL
