#' Internal use of the magrittr pipe
#'
#' Package code uses `%>%`, so the operator is imported into the namespace.
#' It is deliberately NOT re-exported: `%>%` is a general-purpose operator that
#' belongs to magrittr, and a domain package adding it to the user's search
#' path only creates a second place it can come from.
#'
#' @importFrom magrittr %>%
#' @name pipe-internal
#' @keywords internal
#' @noRd
NULL
