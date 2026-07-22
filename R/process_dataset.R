#' Build a sample x compound matrix from TargetLynx or Skyline input
#'
#' Applies (optionally) signal filtering (SNR for TargetLynx; LOD/LOQ for
#' Skyline), blank filtering, missing-value imputation, normalisation,
#' concentration adjustment and batch correction to produce a
#' `struct::DatasetExperiment`.
#'
#' The peak-matrix workflow behind [run_MRManalyzeR()], composed of the public
#' step functions: it reads the workbook ([read_targetlynx()] /
#' [read_skyline()]) into a wide sample x compound matrix, assembles that into a
#' `struct::DatasetExperiment` ([assemble_dataset()]), then applies
#' [filter_blanks()], [normalise_matrix()], [adjust_concentration()] and
#' [impute_missing()] to the dataset per acquisition batch, and optionally
#' batch-corrects the result ([correct_batch()]).
#'
#' @param fdata Feature metadata. Must contain the columns named by
#'   `compound_col`, `processing_name_col`, `report_col` (defaults
#'   `Compound` / `Processing_name` / `Report`). For `signal_filter = "LOD"`
#'   or `"LOQ"`, must also contain the matching `LOD`/`LOQ` column.
#' @param metadata Sample metadata. Must contain the columns named by
#'   `name_col`, `include_col` and the headers referenced by `bc_header`,
#'   `bc_factor_name`, `blank_head`.
#' @param xlsx_path Path to the source workbook (read by [read_targetlynx()]
#'   or [read_skyline()]).
#' @param data_source `"targetlynx"` (default) or `"skyline"` - selects
#'   which raw-matrix reader to use.
#' @param data_tab_names Sheet name(s) holding the Skyline molecule x sample
#'   matrix; only used when `data_source = "skyline"`. Multiple names are
#'   treated as separate acquisition batches and row-bound together - see
#'   `.read_skyline_matrix()`. Sheet membership does not itself imply
#'   batch-correction grouping; that remains driven by `sample_metadata`.
#' @param datatype TargetLynx column to report (e.g. "Area", "Response",
#'   "ng/mL"). Only used when `data_source = "targetlynx"`.
#' @param tl_headers TargetLynx column headers to extract (passed to
#'   [read_targetlynx()]).
#' @param signal_filter `"SNR"`, `"LOD"`, `"LOQ"`, or `FALSE`. `"SNR"` is
#'   only valid for `data_source = "targetlynx"`; `"LOD"`/`"LOQ"` are only
#'   valid for `data_source = "skyline"`. Default `NULL` infers `"SNR"`
#'   unless `snr = FALSE`, for backwards compatibility with configs that
#'   only set `snr:`.
#' @param snr S/N threshold. Only consulted when `signal_filter = "SNR"`.
#' @param blank_filter Blank-filter factor; `FALSE` to skip.
#' @param replace_MVs Scalar passed to [`impute_missing()`]; `FALSE` to skip.
#' @param batch_correction Logical.
#' @param bc_qc_label,bc_factor_name,bc_header Batch-correction parameters.
#' @param blank_head,blank_name Metadata column and value identifying blanks.
#' @param normalize Metadata column used as divisor, or `FALSE`.
#' @param adjust_conc Logical master toggle for concentration adjustment via
#'   [adjust_concentration()], using `sample_vol_col`, `cal_vol_col`,
#'   `sample_IS_col`, `cal_IS_col` and `starting_vol_col`.
#' @param starting_vol_col `FALSE` to skip the starting-volume correction, or
#'   the `metadata` column holding each sample's pre-dry-down starting volume.
#' @param sample_vol_col,cal_vol_col,sample_IS_col,cal_IS_col `metadata`
#'   columns holding, respectively: the vial reconstitution volume, the
#'   calibration-standard vial volume, the IS volume added to the sample,
#'   and the IS volume added to the calibration standard. Only consulted
#'   when `adjust_conc = TRUE`; different studies name these columns
#'   differently, hence configurable rather than hardcoded.
#' @param name_col `sample_metadata` column whose values match the data-matrix
#'   row identifiers (LC-MS injection / sample names). Default `"Name"`.
#' @param include_col,include_value `sample_metadata` inclusion filter: rows
#'   where `metadata[[include_col]] == include_value` are processed. Defaults
#'   `"Include"` / `"YES"`.
#' @param compound_col,processing_name_col `feature_metadata` columns holding
#'   the canonical feature name (becomes the data-matrix column names) and the
#'   raw name that joins to the source data. Defaults `"Compound"` /
#'   `"Processing_name"`.
#' @param report_col,report_value `feature_metadata` inclusion filter: features
#'   where `fdata[[report_col]] == report_value` are reported. Defaults
#'   `"Report"` / `"YES"`.
#' @param comment_col `feature_metadata` column used to record why a feature
#'   was excluded. Default `"Comment"`.
#' @return A list: `[[1]]` the `DatasetExperiment`, `[[2]]` a data frame of
#'   removed features.
#' @examples
#' xlsx  <- system.file("extdata", "example_data.xlsx", package = "MRManalyzeR")
#' fdata <- openxlsx::read.xlsx(xlsx, sheet = "feature_metadata")
#' meta  <- openxlsx::read.xlsx(xlsx, sheet = "sample_metadata")
#' out <- process_dataset(fdata, meta, xlsx_path = xlsx,
#'                        data_source = "targetlynx", datatype = "Area",
#'                        data_tab_names = "lcms_data", snr = 3,
#'                        blank_filter = FALSE, bc_header = "Chrom_Batch",
#'                        blank_head = "Sample_type")
#' out[[1]]
#' @family peak-matrix processing
#' @seealso [run_MRManalyzeR()] to drive the whole workflow from a YAML config.
#' @export
process_dataset = function(fdata,
                              metadata,
                              xlsx_path = NULL,
                              data_source = "targetlynx",
                              data_tab_names = "skyline_data",
                              datatype = "Area",
                              tl_headers = c("ID", "Name", "Area", "ng/mL", "Response", "S/N"),
                              signal_filter = NULL,
                              snr = 5,
                              blank_filter = 5,
                              replace_MVs = FALSE,
                              batch_correction = FALSE,
                              bc_qc_label = "Sample",
                              bc_factor_name = "Sample_type",
                              bc_header = "extract_batch",
                              blank_head = "Sample_type1",
                              blank_name = "ExtractBlank",
                              normalize = FALSE,
                              adjust_conc = FALSE,
                              starting_vol_col = FALSE,
                              sample_vol_col = "sample_volume_uL",
                              cal_vol_col = "cal_vol_uL",
                              sample_IS_col = "sample_IS_vol_uL",
                              cal_IS_col = "cal_IS_vol_uL",
                              name_col = "Name",
                              include_col = "Include",
                              include_value = "YES",
                              compound_col = "Compound",
                              processing_name_col = "Processing_name",
                              report_col = "Report",
                              report_value = "YES",
                              comment_col = "Comment"){

  for(.c in c(name_col, include_col)){
    if(!.c %in% colnames(metadata))
      stop(sprintf("[process_dataset] sample_metadata has no '%s' column.", .c))
  }
  for(.c in c(compound_col, processing_name_col, report_col)){
    if(!.c %in% colnames(fdata))
      stop(sprintf("[process_dataset] feature_metadata has no '%s' column.", .c))
  }

  # Canonicalise feature_metadata to the internal standard columns/values
  # (Compound / Processing_name / Report=YES|NO / Comment) so the rest of the
  # package - engine, reports, subset_dataset, combine - keeps using fixed names
  # regardless of what the input sheet calls them.
  fdata$Processing_name = fdata[[processing_name_col]]
  fdata$Compound        = fdata[[compound_col]]
  fdata$Report          = ifelse(as.character(fdata[[report_col]]) ==
                                   as.character(report_value), "YES", "NO")
  if(comment_col %in% colnames(fdata)) fdata$Comment = fdata[[comment_col]]

  if(!identical(data_source, "targetlynx") && !identical(data_source, "skyline"))
    stop("[process_dataset] data_source must be 'targetlynx' or 'skyline'.")

  # Backwards compatibility: YAMLs that only set `snr:` (no `signal_filter:`)
  # keep working unchanged - infer the filter mode from `snr`.
  if(is.null(signal_filter))
    signal_filter = if(isFALSE(snr)) FALSE else "SNR"

  if(!identical(signal_filter, FALSE) && !signal_filter %in% c("SNR", "LOD", "LOQ"))
    stop("[process_dataset] signal_filter must be one of: 'SNR', 'LOD', 'LOQ', FALSE.")

  if(identical(signal_filter, "SNR") && identical(data_source, "skyline"))
    stop("[process_dataset] signal_filter='SNR' requires per-injection S/N values, which are not available for data_source='skyline'. Use signal_filter='LOD' or 'LOQ' instead.")

  if(signal_filter %in% c("LOD", "LOQ") && !identical(data_source, "skyline"))
    stop(sprintf(
      "[process_dataset] signal_filter='%s' is only supported for data_source='skyline'. For data_source='targetlynx', use signal_filter='SNR'.",
      signal_filter))

  if(signal_filter %in% c("LOD", "LOQ") && !signal_filter %in% colnames(fdata))
    stop(sprintf(
      "[process_dataset] signal_filter='%s' but feature_metadata has no '%s' column.",
      signal_filter, signal_filter))

  # --- Validation ---------------------------------------------------------
  for(h in c(bc_factor_name, bc_header, blank_head)){
    if(!h %in% colnames(metadata))
      warning(sprintf("'%s' not a metadata column.", h))
  }

  # --- Whitelists ----------------------------------------------------------
  fnames = fdata$Processing_name[fdata$Report == "YES"]
  snames = metadata[[name_col]][which(metadata[[include_col]] == include_value)]

  # --- Read raw matrix (format-specific) + apply signal filter ------------
  if(identical(data_source, "targetlynx")){

    # read_targetlynx() reads + pivots + (optionally) SNR-masks the full sheet;
    # restrict to reportable compounds (fnames) and included samples (snames)
    # here, exactly as the old pre-pivot filter did.
    out_table = read_targetlynx(xlsx_path, datatype = datatype,
                                tl_headers = tl_headers,
                                snr = if(identical(signal_filter, "SNR")) snr else FALSE,
                                data_tab_names = data_tab_names)
    out_table = out_table[rownames(out_table) %in% snames,
                          colnames(out_table) %in% fnames, drop = FALSE]

  } else {   # data_source == "skyline"

    if(is.null(xlsx_path))
      stop("[process_dataset] xlsx_path is required when data_source='skyline'.")

    out_table = read_skyline(xlsx_path, data_tab_names, fdata,
                             signal_filter = if(signal_filter %in% c("LOD", "LOQ")) signal_filter else FALSE)
    out_table = out_table[rownames(out_table) %in% snames, , drop = FALSE]
  }

  # --- Mismatch diagnostics ----------------------------------------------
  # Track features that exist in fdata$Processing_name but have no row in
  # the source data (compound never reported for any sample, regardless of
  # data_source). Tracked here so they can be surfaced in the QC report and
  # the console.
  missing_in_data = setdiff(fnames, colnames(out_table))
  if(length(missing_in_data)){
    message(sprintf(
      "[process_dataset] %d feature(s) in fdata$Processing_name have no rows in the source data:\n  %s",
      length(missing_in_data),
      paste(missing_in_data, collapse = ", ")))
  }

  # Rename Processing_name -> Compound via O(1) lookup. Any source column
  # not present in fdata$Processing_name maps to NA - those columns carry
  # data but have no feature_meta row, so we must drop them before
  # assembling the DatasetExperiment (otherwise SE complains that assay
  # colnames disagree with rowData rownames). For data_source="skyline",
  # this is also where unmatched Molecule names (e.g. internal standards
  # with no feature_metadata row) get caught.
  orig_names  = colnames(out_table)
  name_lookup = stats::setNames(fdata$Compound, fdata$Processing_name)
  new_names   = name_lookup[orig_names]
  unmapped    = is.na(new_names)
  unmapped_ids = orig_names[unmapped]   # surfaced in QC report
  if(length(unmapped_ids)){
    message(sprintf(
      "[process_dataset] %d source column(s) have no matching fdata$Processing_name and were dropped:\n  %s",
      length(unmapped_ids),
      paste(unmapped_ids, collapse = ", ")))
    out_table = out_table[, !unmapped, drop = FALSE]
    new_names = new_names[!unmapped]
  }
  colnames(out_table) = new_names

  # Drop all-NA features (vectorized via colSums)
  keep_col = colSums(!is.na(out_table)) > 0
  remove_feats = names(out_table)[!keep_col]
  out_table    = out_table[, keep_col, drop = FALSE]

  # --- Per-batch processing ----------------------------------------------
  meta_yes = metadata[which(metadata[[include_col]] == include_value), , drop = FALSE]

  # Fail loud and early if bc_header is misconfigured - silent fallback
  # would mask YAML typos and produce an empty downstream matrix.
  if(is.null(bc_header) || isFALSE(bc_header) ||
     !bc_header %in% colnames(meta_yes)){
    stop(sprintf(
      "bc_header='%s' is not a column in sample_metadata. Fix the YAML 'bc_header:' to match an existing sample_metadata column (got: %s).",
      bc_header %||% "<NULL>",
      paste(colnames(meta_yes), collapse = ", ")))
  }

  # --- Assemble the DatasetExperiment, then process it --------------------
  # The DE is built first so every processing step downstream operates on the
  # dataset object (data + sample_meta together) rather than a bare matrix.
  if(nrow(out_table) == 0 || ncol(out_table) == 0){
    stop(sprintf(
      "[process_dataset] Empty matrix after reading (%d samples x %d features). Check: bc_header column, blank/signal filters, and that the source data actually contains data for the included samples.",
      nrow(out_table), ncol(out_table)))
  }
  lcms_experiment = assemble_dataset(out_table, fdata, metadata, name_col = name_col)

  # Treat NA / "" bc_header values (e.g. blanks, QCs that weren't
  # cryosectioned) as a single "_unbatched_" group so they survive the
  # per-batch split instead of being silently dropped via `NA == "x"`.
  # Computed locally: the dataset's own sample_meta keeps the original values.
  bc_vec = as.character(as.data.frame(lcms_experiment$sample_meta)[[bc_header]])
  bc_vec[is.na(bc_vec) | !nzchar(bc_vec)] = "_unbatched_"
  batches = unique(bc_vec)

  if(length(batches) == 1L){
    # fast-path: process the whole dataset, no split/rebind
    lcms_experiment = .process_batch(
      lcms_experiment,
      blank_head       = blank_head,
      blank_name       = blank_name,
      blank_filter     = blank_filter,
      normalize        = normalize,
      adjust_conc      = adjust_conc,
      starting_vol_col = starting_vol_col,
      sample_vol_col   = sample_vol_col,
      cal_vol_col      = cal_vol_col,
      sample_IS_col    = sample_IS_col,
      cal_IS_col       = cal_IS_col,
      replace_MVs      = replace_MVs
    )
  } else {
    row_order = rownames(as.data.frame(lcms_experiment$data))
    batch_frames = lapply(batches, function(b){
      as.data.frame(.process_batch(
        .de_rows(lcms_experiment, which(bc_vec == b)),
        blank_head       = blank_head,
        blank_name       = blank_name,
        blank_filter     = blank_filter,
        normalize        = normalize,
        adjust_conc      = adjust_conc,
        starting_vol_col = starting_vol_col,
        sample_vol_col   = sample_vol_col,
        cal_vol_col      = cal_vol_col,
        sample_IS_col    = sample_IS_col,
        cal_IS_col       = cal_IS_col,
        replace_MVs      = replace_MVs
      )$data)
    })
    merged = do.call(rbind, batch_frames)
    lcms_experiment$data = merged[row_order, , drop = FALSE]
  }

  if(isTRUE(batch_correction))
    lcms_experiment = correct_batch(lcms_experiment,
                                    qc_label    = bc_qc_label,
                                    factor_name = bc_factor_name,
                                    batch_head  = bc_header)

  # --- Record removed features -------------------------------------------
  fail_reason = if(identical(signal_filter, FALSE)) "All-NA after processing"
                else sprintf("All values below %s threshold", signal_filter)
  fdata$Report [fdata$Compound %in% remove_feats] = "NO"
  fdata$Comment[fdata$Compound %in% remove_feats] = fail_reason

  # Flag features whose Processing_name appears in fdata but has no rows
  # in the source data (compound never reported for any sample).
  if(length(missing_in_data)){
    miss_row = fdata$Processing_name %in% missing_in_data
    fdata$Report [miss_row] = "NO"
    fdata$Comment[miss_row] = "No matching data in source file"
  }
  removed_features = dplyr::filter(fdata, Report == "NO")

  # Append unmapped source columns (no fdata entry at all) so the QC report
  # also surfaces those - fdata has no row for them, so we add stub rows.
  if(length(unmapped_ids)){
    extra = data.frame(
      Compound = unmapped_ids,
      Report   = "NO",
      Comment  = "No matching Processing_name in feature_metadata",
      stringsAsFactors = FALSE
    )
    # pad to the same columns as removed_features
    miss_cols = setdiff(colnames(removed_features), colnames(extra))
    for(m in miss_cols) extra[[m]] = NA
    extra = extra[, colnames(removed_features), drop = FALSE]
    removed_features = rbind(removed_features, extra)
  }

  # Console summary banner so mismatches are impossible to miss
  if(length(unmapped_ids) || length(missing_in_data)){
    message("\n--- Feature-mismatch summary ---")
    message(sprintf("  fdata Processing_name with no source data : %d", length(missing_in_data)))
    message(sprintf("  source columns with no fdata entry        : %d", length(unmapped_ids)))
    message("  (both groups listed in 'Compounds excluded' in the QC report)\n")
  }

  list(lcms_experiment, removed_features)
}


# ---- Internal helpers ------------------------------------------------------

#' Pivot the long TargetLynx table to a wide sample x compound matrix
#' @keywords internal
#' @noRd
.build_wide_matrix = function(base, datatype){
  base %>%
    dplyr::select(ID, Name, dplyr::all_of(datatype)) %>%
    dplyr::distinct() %>%
    dplyr::mutate(dplyr::across(dplyr::all_of(datatype), as.numeric)) %>%
    tidyr::pivot_wider(names_from = "ID", values_from = datatype,
                       id_cols = "Name", values_fill = 0) %>%
    tibble::column_to_rownames("Name")
}

#' Apply an S/N threshold mask to a wide matrix
#' @keywords internal
#' @noRd
.apply_snr_mask = function(out_table, base, snr){
  snr_table = base %>%
    dplyr::select(ID, Name, `S/N`) %>%
    dplyr::distinct() %>%
    dplyr::mutate(`S/N` = suppressWarnings(as.numeric(`S/N`))) %>%
    tidyr::pivot_wider(names_from = "ID", values_from = "S/N",
                       id_cols = "Name") %>%
    tibble::column_to_rownames("Name")

  snr_table = snr_table[rownames(out_table), colnames(out_table), drop = FALSE]
  out_table[snr_table < snr] = NA
  out_table
}

#' Run blank-filter / normalise / concentration / MV-impute on one batch
#'
#' Operates on a `DatasetExperiment` so each step reads what it needs
#' (blank rows, divisors, volumes) straight from `sample_meta`.
#' @keywords internal
#' @noRd
.process_batch = function(de, blank_head, blank_name,
                          blank_filter, normalize,
                          adjust_conc, starting_vol_col,
                          sample_vol_col, cal_vol_col, sample_IS_col, cal_IS_col,
                          replace_MVs){

  if(!isFALSE(blank_filter))
    de = filter_blanks(de, blank_filter = blank_filter,
                       blank_head = blank_head, blank_name = blank_name)

  if(!isFALSE(normalize))
    de = normalise_matrix(de, column = normalize)

  # Vial-level concentration correction (+ optional starting-volume
  # correction) -- see adjust_concentration.R.
  if(isTRUE(adjust_conc))
    de = adjust_concentration(de,
                              sample_vol_col   = sample_vol_col,
                              cal_vol_col      = cal_vol_col,
                              sample_IS_col    = sample_IS_col,
                              cal_IS_col       = cal_IS_col,
                              starting_vol_col = starting_vol_col)

  if(!isFALSE(replace_MVs))
    de = impute_missing(de, scalar = replace_MVs,
                        blank_head = blank_head, blank_name = blank_name)

  de
}

#' Subset a DatasetExperiment to the given sample rows, keeping all features
#' @keywords internal
#' @noRd
.de_rows = function(de, idx){
  struct::DatasetExperiment(
    data          = as.data.frame(de$data)[idx, , drop = FALSE],
    sample_meta   = as.data.frame(de$sample_meta)[idx, , drop = FALSE],
    variable_meta = as.data.frame(de$variable_meta)
  )
}
