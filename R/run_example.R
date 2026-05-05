#' Render the bundled example MRManalyzeR reports
#'
#' Convenience wrapper that drives the full pipeline against the bundled
#' oxylipin example dataset (a subset of Kolmert et al. 2018,
#' \doi{10.1016/j.prostaglandins.2018.05.005}). It loads
#' `inst/extdata/example_data.xlsx` and `inst/extdata/example_config.yml`,
#' rewrites the YAML's `paths:` block to point at a writable results
#' directory + the bundled xlsx, and calls [`run_MRManalyzeR()`]. After the
#' run the generated HTML reports are opened in interactive sessions.
#'
#' @param results_dir Directory to write outputs (xlsx, RDS, HTML reports)
#'   to. Defaults to a fresh `tempdir()` so repeated calls do not collide.
#' @param open Logical. Open the HTML reports in the browser / RStudio
#'   viewer after rendering? Default `TRUE` in interactive sessions.
#' @return Invisibly, the list returned by [`run_MRManalyzeR()`]; the
#'   resolved results directory is attached as attribute `"results_dir"`.
#' @examples
#' \dontrun{
#'   res <- MRManalyzeR::run_example()
#'   attr(res, "results_dir")   # where the xlsx + html landed
#' }
#' @export
run_example = function(results_dir = tempfile("MRManalyzeR_example_"),
                       open        = interactive()){

  data_xlsx = system.file("extdata", "example_data.xlsx",
                          package = "MRManalyzeR")
  base_yaml = system.file("extdata", "example_config.yml",
                          package = "MRManalyzeR")
  if(!nzchar(data_xlsx) || !file.exists(data_xlsx))
    stop("[run_example] bundled example_data.xlsx not found ",
         "— reinstall MRManalyzeR.")
  if(!nzchar(base_yaml) || !file.exists(base_yaml))
    stop("[run_example] bundled example_config.yml not found ",
         "— reinstall MRManalyzeR.")

  dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

  # `load_yaml()` already resolves !concat_path so the returned `cfg` has
  # plain string paths. We override the data + result locations to point at
  # the bundled xlsx and a writable temp folder, then write the modified
  # YAML alongside the outputs so the user can inspect / edit it.
  cfg = load_yaml(base_yaml)
  cfg$project$paths$data_dir    = dirname(data_xlsx)
  cfg$project$paths$result_dir  = results_dir
  cfg$project$paths$project_dir = results_dir
  cfg$project$paths$fn          = tools::file_path_sans_ext(basename(data_xlsx))

  out_yaml = file.path(results_dir, "example_config.yml")
  yaml::write_yaml(cfg, out_yaml)

  message("[run_example] outputs -> ", results_dir)
  res = run_MRManalyzeR(out_yaml)

  if(isTRUE(open)){
    htmls = list.files(results_dir, pattern = "\\.html$", full.names = TRUE)
    for(h in htmls) utils::browseURL(h)
  }

  attr(res, "results_dir") = results_dir
  invisible(res)
}
