#' Load a project YAML, resolving !concat_path and !concat_str anchors
#'
#' @param path_yaml Full path to YAML file
#' @return Named list of parameters parsed from the YAML
#' @export
load_yaml = function(path_yaml){
  yaml::read_yaml(
    path_yaml,
    handlers = list(
      concat_path = function(x) paste(x, collapse = "/"),
      concat_str  = function(x) paste(x, collapse = "_")
    )
  )
}
