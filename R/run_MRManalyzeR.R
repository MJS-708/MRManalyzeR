#' Run the end-to-end MRManalyzeR workflow from a YAML config
#'
#' Single entry point that:
#' \itemize{
#'   \item loads a project YAML,
#'   \item (optionally) reads the TargetLynx xlsx workbook and builds the
#'     processed peak matrix — gated by `PeakMatrixProcessing.execute`,
#'   \item writes the processed matrix to xlsx + RDS and persists the YAML
#'     parameters alongside,
#'   \item runs the statistical analyses defined under `stats_report:`
#'     (comparisons / correlations / linear_models) and writes them as
#'     extra tabs in the results xlsx,
#'   \item (optionally) renders the bundled `data_quality_report` and
#'     `stats_report` HTML vignettes.
#' }
#'
#' Output filenames are derived from `paths.fn` + `paths.datatype` +
#' `paths.suffix`, so when `PeakMatrixProcessing.execute: False` the code
#' can still locate the existing RDS without re-reading the source xlsx.
#'
#' Backwards compatibility:
#' \itemize{
#'   \item `paths.datatype` falls back to `PeakMatrixProcessing.datatype` if absent.
#'   \item Legacy `UVA_report:` / `MVA_report:` blocks are recognised when
#'     `data_quality_report:` / `stats_report:` are missing.
#' }
#'
#' @param path_yaml Full path to a project YAML.
#' @return Invisibly, a list with the `DatasetExperiment`, `removed_features`,
#'   `stats_tables`, and the resolved output paths. When `paths.datatype` is a
#'   vector (e.g. `["Area", "Response", "Conc"]`), the function loops over each
#'   datatype and returns a list of per-datatype results.
#' @export
run_MRManalyzeR = function(path_yaml){

  stopifnot(file.exists(path_yaml))

  # --- Load & unpack YAML -------------------------------------------------
  project_params = load_yaml(path_yaml)

  project_paths = project_params$project$paths
  pmp_params    = project_params$project$PeakMatrixProcessing

  # Prefer new combined blocks; fall back to legacy UVA/MVA blocks
  dq_params = project_params$project$data_quality_report %||%
                project_params$project$UVA_report
  st_params = project_params$project$stats_report %||%
                project_params$project$MVA_report


  # --- Datatype loop ------------------------------------------------------
  # `datatype` may live in paths: (new) or PeakMatrixProcessing: (legacy).
  # Accept either a scalar ("Area") or a vector (["Area", "Response", "Conc"]).
  datatypes = project_paths$datatype %||% pmp_params$datatype %||% "Area"
  if(length(datatypes) > 1){
    results = lapply(datatypes, function(dt){
      message(sprintf("\n========== Datatype: %s ==========", dt))
      pp_one = project_paths;             pp_one$datatype = dt
      pmp_one = pmp_params;               pmp_one$datatype = dt
      pp_proj_one = project_params
      pp_proj_one$project$paths = pp_one
      pp_proj_one$project$PeakMatrixProcessing = pmp_one
      .run_MRManalyzeR_one(project_params = pp_proj_one,
                           project_paths  = pp_one,
                           pmp_params     = pmp_one,
                           dq_params      = dq_params,
                           st_params      = st_params,
                           datatype       = dt)
    })
    names(results) = datatypes
    return(invisible(results))
  }
  datatype = datatypes[[1]]

  invisible(.run_MRManalyzeR_one(
    project_params = project_params,
    project_paths  = project_paths,
    pmp_params     = pmp_params,
    dq_params      = dq_params,
    st_params      = st_params,
    datatype       = datatype
  ))
}


#' Single-datatype runner extracted from [`run_MRManalyzeR()`].
#' @keywords internal
#' @noRd
.run_MRManalyzeR_one = function(project_params, project_paths, pmp_params,
                                dq_params, st_params, datatype){

  units_tag = gsub("/", "_", datatype)
  out_stub  = sprintf("%s/%s_%s%s",
                      project_paths$result_dir,
                      project_paths$fn,
                      units_tag,
                      project_paths$suffix %||% "")

  out_xlsx       = paste0(out_stub, ".xlsx")
  out_stats_xlsx = paste0(out_stub, "_stats.xlsx")
  out_RDS        = paste0(out_stub, ".RDS")
  out_pars       = paste0(out_stub, ".txt")

  dir.create(project_paths$result_dir, recursive = TRUE, showWarnings = FALSE)

  # --- Peak-matrix processing (toggle) -----------------------------------
  if(isTRUE(pmp_params$execute)){
    out = .run_pmp(project_paths, pmp_params, datatype,
                   out_xlsx = out_xlsx, out_RDS = out_RDS, out_pars = out_pars,
                   project_params = project_params)
    combined_datamatrices = out$combined_datamatrices
    removed_features      = out$removed_features
  } else {
    if(!file.exists(out_RDS))
      stop(sprintf(
        "PeakMatrixProcessing.execute is FALSE but expected RDS not found:\n  %s\n%s",
        out_RDS,
        "Re-run with execute: True at least once, or fix paths.fn / paths.datatype."
      ))
    message("Skipping PeakMatrixProcessing (execute: False). Loading: ", out_RDS)
    combined_datamatrices = readRDS(out_RDS)
    removed_features      = data.frame()    # not available without re-processing
  }

  # --- Statistics ---------------------------------------------------------
  stats_tables = NULL
  if(isTRUE(st_params$execute) &&
     (length(st_params$comparisons) || length(st_params$correlations) ||
        length(st_params$linear_models))){
    message("Running statistics...")
    stats_tables = run_stats(combined_datamatrices, st_params)
    append_stats_xlsx(out_stats_xlsx, stats_tables)
    message("Wrote stats to: ", out_stats_xlsx)
  }

  # --- Optional reports ---------------------------------------------------
  .render_reports(project_params   = project_params,
                  project_paths    = project_paths,
                  pmp_params       = pmp_params,
                  dq_params        = dq_params,
                  st_params        = st_params,
                  out_RDS          = out_RDS,
                  out_stub         = out_stub,
                  combined_datamatrices = combined_datamatrices,
                  removed_features      = removed_features,
                  stats_tables          = stats_tables)

  invisible(list(
    datasetExperiment = combined_datamatrices,
    removed_features  = removed_features,
    stats_tables      = stats_tables,
    out_xlsx       = out_xlsx,
    out_stats_xlsx = out_stats_xlsx,
    out_RDS        = out_RDS,
    out_pars       = out_pars,
    datatype       = datatype
  ))
}


# --- Internal helpers -------------------------------------------------------

`%||%` = function(a, b) if(is.null(a)) b else a

#' @keywords internal
#' @noRd
.run_pmp = function(project_paths, pmp_params, datatype,
                    out_xlsx, out_RDS, out_pars, project_params){

  xlsx_path = sprintf("%s/%s.xlsx", project_paths$data_dir, project_paths$fn)
  message("Reading ", xlsx_path)

  # Path-based reads avoid an openxlsx::loadWorkbook() bug on workbooks
  # with certain styling/drawing XML, and are also faster.
  fdata    = openxlsx::read.xlsx(xlsx_path, sheet = "feature_metadata")
  metadata = openxlsx::read.xlsx(xlsx_path, sheet = "sample_metadata") %>%
    dplyr::filter(Include == "YES")

  message("Extracting TargetLynx tables...")
  lcms_table = extractTable(xlsx_path, tl_headers = pmp_params$tl_headers) %>%
    subset(Name != "")

  message("Generating data matrix...")
  combined_data = generateDataMatrix(
    fdata            = fdata,
    metadata         = metadata,
    lcms_table       = lcms_table,
    datatype         = datatype,
    tl_headers       = pmp_params$tl_headers,
    snr              = pmp_params$snr,
    blank_filter     = pmp_params$blank_filter,
    replace_MVs      = pmp_params$replace_MVs,
    batch_correction = pmp_params$batch_correction,
    bc_qc_label      = pmp_params$bc_qc_label,
    bc_factor_name   = pmp_params$bc_factor_name,
    bc_header        = pmp_params$bc_header,
    blank_head       = pmp_params$blank_head,
    blank_name       = pmp_params$blank_name,
    normalize        = pmp_params$normalize,
    adjust_conc      = pmp_params$adjust_conc
  )

  combined_datamatrices = combined_data[[1]]
  removed_features      = combined_data[[2]]

  scale_fac = pmp_params$scale_fac %||% 1
  combined_datamatrices$data = combined_datamatrices$data * scale_fac

  add_info = data.frame(sheet = "matrix", info = "units",
                        value = pmp_params$units %||% datatype)

  .write_output_xlsx(out_xlsx, combined_datamatrices, removed_features, add_info)
  saveRDS(combined_datamatrices, file = out_RDS)

  con = file(out_pars, open = "wt"); on.exit(close(con), add = TRUE)
  utils::capture.output(print(project_params), file = con)

  message("Wrote: ", out_xlsx)
  message("Wrote: ", out_RDS)

  list(combined_datamatrices = combined_datamatrices,
       removed_features      = removed_features)
}

#' @keywords internal
#' @noRd
.write_output_xlsx = function(out_xlsx, combined_datamatrices,
                              removed_features, add_info){
  wb = openxlsx::createWorkbook()

  openxlsx::addWorksheet(wb, "feature_metadata")
  openxlsx::writeData(wb, "feature_metadata",
                      combined_datamatrices$variable_meta,
                      colNames = TRUE, rowNames = FALSE, keepNA = FALSE)

  openxlsx::addWorksheet(wb, "sample_metadata")
  openxlsx::writeData(wb, "sample_metadata",
                      combined_datamatrices$sample_meta,
                      colNames = TRUE, rowNames = FALSE, keepNA = FALSE)

  openxlsx::addWorksheet(wb, "matrix")
  openxlsx::writeData(wb, "matrix",
                      combined_datamatrices$data %>%
                        dplyr::select(dplyr::where(function(x) any(!is.na(x)))),
                      colNames = TRUE, rowNames = TRUE, keepNA = FALSE)

  openxlsx::addWorksheet(wb, "removed_features")
  openxlsx::writeData(wb, "removed_features", removed_features,
                      colNames = TRUE, rowNames = FALSE, keepNA = FALSE)

  openxlsx::addWorksheet(wb, "key")
  openxlsx::writeData(wb, "key", add_info,
                      colNames = TRUE, rowNames = FALSE, keepNA = FALSE)

  openxlsx::saveWorkbook(wb, file = out_xlsx, overwrite = TRUE)
}


#' @keywords internal
#' @noRd
.render_reports = function(project_params, project_paths, pmp_params,
                           dq_params, st_params,
                           out_RDS, out_stub,
                           combined_datamatrices, removed_features,
                           stats_tables = NULL){

  render_one = function(rmd_file, out_html, params_block){

    rmd_path = system.file("rmd", rmd_file, package = "MRManalyzeR")
    if (rmd_path == "") {
      rmd_path = file.path("inst", "rmd", rmd_file)
    }
    if(rmd_path == ""){
      warning(sprintf("Report template '%s' not found in installed package.",
                      rmd_file))
      return(invisible(NULL))
    }

    # Parent on the package namespace so the Rmd can find package
    # functions (run_pca_pipeline, de_subset, ...) even when the caller
    # invoked us via `MRManalyzeR::run_MRManalyzeR()` without
    # `library(MRManalyzeR)`. Lookup chain: chunk env -> render env ->
    # MRManalyzeR namespace -> imports -> base.
    pkg_ns = tryCatch(asNamespace("MRManalyzeR"),
                      error = function(e) globalenv())
    env = new.env(parent = pkg_ns)
    env$project_params        = project_params
    env$project_paths         = project_paths
    env$pmp_params            = pmp_params
    env$report_pars           = params_block
    env$include_MVA           = isTRUE(params_block$include_MVA) ||
                                 isTRUE(params_block$pca$include)
    env$out_RDS               = out_RDS
    env$combined_datamatrices = combined_datamatrices
    env$removed_features      = removed_features
    env$stats_tables          = stats_tables

    message("Rendering: ", out_html)
    rmarkdown::render(
      input       = rmd_path,
      output_file = out_html,
      envir       = env,
      quiet       = TRUE
    )
  }

  if(isTRUE(dq_params$execute)){
    render_one("data_quality_report.Rmd",
               paste0(out_stub, "_data_quality_report.html"),
               dq_params)
  }
  if(isTRUE(st_params$execute)){
    render_one("stats_report.Rmd",
               paste0(out_stub, "_stats_report.html"),
               st_params)
  }
}


#' Run the combine-mode pipeline from a dedicated combine YAML
#'
#' Merges several previously-saved datasets (`.RDS` or `.xlsx` outputs of
#' [`run_MRManalyzeR()`]) into one `DatasetExperiment`, then runs the
#' `stats_report` analyses on the merged data. Skips PeakMatrixProcessing
#' and the data_quality_report.
#'
#' Defaults: samples are matched by `Sample_ID` (intersect across inputs);
#' `feature_meta_cols` and `sample_meta_cols` default to the intersection
#' of column names across inputs (specify them only to *narrow* the set).
#'
#' Outputs:
#' \itemize{
#'   \item `<output_stub>.RDS`        — merged `DatasetExperiment`
#'   \item `<output_stub>.xlsx`       — feature_metadata / sample_metadata / matrix tabs
#'   \item `<output_stub>_stats.xlsx` — stats / correlations / linear_models
#'   \item `<output_stub>_stats_report.html`
#' }
#'
#' YAML schema (see `inst/extdata/example_combine_config.yml`):
#' \preformatted{
#' combine:
#'   datasets:
#'     - path: ".../GOM_..._.RDS"      # .RDS or .xlsx, mix is fine
#'       tag:  "GOM"
#'       qc_remap: { "QC1": "QC_pooled" }
#'     - path: ".../cysLT_..._.xlsx"
#'       tag:  "cysLT"
#'   feature_meta_cols: ~              # NULL = intersect across inputs
#'   sample_meta_cols:  ~              # NULL = intersect across inputs
#'   prefix_features:   True
#'   sample_id_col:     Sample_ID
#'   output_stub:       ".../Combined/combined_HDM"
#'   units:             "combined"
#' stats_report:
#'   execute: True
#'   comparisons:   [...]
#'   correlations:  [...]
#'   linear_models: [...]
#' }
#'
#' @param path_yaml Full path to the combine YAML.
#' @return Invisibly, a list with the merged `DatasetExperiment`, the
#'   `stats_tables`, and the output paths.
#' @export
run_MRManalyzeR_combine = function(path_yaml){

  stopifnot(file.exists(path_yaml))
  cfg = load_yaml(path_yaml)

  # Tolerate either a flat top-level layout or a `project:` wrapper.
  root = cfg$project %||% cfg
  combine_params = root$combine %||%
                     stop("[combine] YAML lacks a top-level `combine:` block.")
  st_params      = root$stats_report %||% list(execute = FALSE)

  ds = combine_params$datasets %||% list()
  if(!length(ds)) stop("[combine] `combine.datasets:` is empty.")

  paths = vapply(ds, function(d) d$path, character(1))
  tags  = vapply(ds, function(d) d$tag %||% NA_character_, character(1))
  if(all(!is.na(tags) & nzchar(tags))) names(paths) = tags

  # Per-dataset config maps — keyed by tag if available, otherwise by path.
  key_for = function(i) if(!is.na(tags[i]) && nzchar(tags[i])) tags[i] else paths[i]
  build_map = function(field){
    out = list()
    for(i in seq_along(ds)){
      v = ds[[i]][[field]]
      if(!is.null(v)) out[[ key_for(i) ]] = v
    }
    if(!length(out)) NULL else out
  }
  qc_remap_list  = build_map("qc_remap")
  rename_list    = build_map("feature_meta_rename")

  output_stub = combine_params$output_stub %||%
                stop("[combine] `combine.output_stub:` is required.")
  dir.create(dirname(output_stub), recursive = TRUE, showWarnings = FALSE)

  message(sprintf("[combine] Merging %d dataset(s) ...", length(ds)))
  combined = combine_datasets(
    paths               = paths,
    feature_meta_cols   = combine_params$feature_meta_cols,
    sample_meta_cols    = combine_params$sample_meta_cols,
    feature_meta_rename = rename_list,
    qc_remap            = qc_remap_list,
    sample_id_col       = combine_params$sample_id_col   %||% "Sample_ID",
    drop_samples        = combine_params$drop_samples,
    prefix_features     = isTRUE(combine_params$prefix_features),
    combined_name       = combine_params$combined_name   %||% basename(output_stub)
  )
  message(sprintf("[combine] Result: %d samples × %d features.",
                  nrow(combined$data), ncol(combined$data)))

  out_xlsx       = paste0(output_stub, ".xlsx")
  out_stats_xlsx = paste0(output_stub, "_stats.xlsx")
  out_RDS        = paste0(output_stub, ".RDS")
  out_pars       = paste0(output_stub, ".txt")

  add_info = data.frame(
    sheet = c("matrix", "matrix"),
    info  = c("source", "n_input_datasets"),
    value = c("combine_datasets()", as.character(length(ds)))
  )
  .write_output_xlsx(out_xlsx, combined,
                     removed_features = data.frame(),
                     add_info          = add_info)
  saveRDS(combined, file = out_RDS)
  con = file(out_pars, open = "wt"); on.exit(close(con), add = TRUE)
  utils::capture.output(print(cfg), file = con)
  message("Wrote: ", out_xlsx)
  message("Wrote: ", out_RDS)

  stats_tables = NULL
  if(isTRUE(st_params$execute) &&
     (length(st_params$comparisons) || length(st_params$correlations) ||
        length(st_params$linear_models))){
    message("[combine] Running statistics on merged data ...")
    stats_tables = run_stats(combined, st_params)
    append_stats_xlsx(out_stats_xlsx, stats_tables)
    message("Wrote stats to: ", out_stats_xlsx)
  }

  if(isTRUE(st_params$execute)){
    rmd_path = system.file("rmd", "stats_report.Rmd", package = "MRManalyzeR")
    if(rmd_path == "")
      rmd_path = file.path("inst", "rmd", "stats_report.Rmd")
    pkg_ns = tryCatch(asNamespace("MRManalyzeR"),
                      error = function(e) globalenv())
    env = new.env(parent = pkg_ns)
    units_lbl = combine_params$units %||% "combined"
    env$project_params        = cfg
    env$project_paths         = list(datatype = units_lbl)
    env$pmp_params            = list(units = units_lbl, datatype = units_lbl)
    env$report_pars           = st_params
    env$out_RDS               = out_RDS
    env$combined_datamatrices = combined
    env$removed_features      = data.frame()
    env$stats_tables          = stats_tables

    out_html = paste0(output_stub, "_stats_report.html")
    message("Rendering: ", out_html)
    rmarkdown::render(input = rmd_path, output_file = out_html,
                      envir = env, quiet = TRUE)
  }

  invisible(list(
    datasetExperiment = combined,
    stats_tables      = stats_tables,
    out_xlsx          = out_xlsx,
    out_stats_xlsx    = out_stats_xlsx,
    out_RDS           = out_RDS,
    out_pars          = out_pars,
    n_input_datasets  = length(ds)
  ))
}
