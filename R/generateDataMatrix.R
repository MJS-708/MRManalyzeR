#' Build a sample x compound matrix from a stacked TargetLynx long table
#'
#' Applies (optionally) SNR filtering, blank filtering, missing-value
#' imputation, normalisation, concentration adjustment and batch correction
#' to produce a `struct::DatasetExperiment`.
#'
#' Top-level orchestrator. Heavy lifting is delegated to the internal
#' helpers `.build_wide_matrix()`, `.apply_snr_mask()`, `.process_batch()`
#' and `.assemble_DE()`.
#'
#' @param fdata Feature metadata. Must contain columns `Processing_name`,
#'   `Compound`, `Report`, `Comment`.
#' @param metadata Sample metadata. Must contain `Name`, `Include` and the
#'   headers referenced by `bc_header`, `bc_factor_name`, `blank_head`.
#' @param lcms_table Long table from [`extractTable()`].
#' @param datatype TargetLynx column to report (e.g. "Area", "Response", "ng/mL").
#' @param tl_headers Retained for API compatibility; not used internally.
#' @param snr S/N threshold. Peaks with S/N below this are set to NA. `FALSE` to skip.
#' @param blank_filter Blank-filter factor; `FALSE` to skip.
#' @param replace_MVs Scalar passed to [`replaceNA()`]; `FALSE` to skip.
#' @param batch_correction Logical.
#' @param bc_qc_label,bc_factor_name,bc_header Batch-correction parameters.
#' @param blank_head,blank_name Metadata column and value identifying blanks.
#' @param normalize Metadata column used as divisor, or `FALSE`.
#' @param adjust_conc Logical. If `TRUE`, volume-adjust via [`calculate_conc()`].
#' @param use_type Retained for API compatibility.
#' @return A list: `[[1]]` the `DatasetExperiment`, `[[2]]` a data frame of
#'   removed features.
#' @export
generateDataMatrix = function(fdata,
                              metadata,
                              lcms_table,
                              datatype = "Area",
                              tl_headers = c("ID", "Name", "Area", "ng/mL", "Response", "S/N"),
                              use_type = "Ratio",
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
                              adjust_conc = FALSE){

  # --- Validation ---------------------------------------------------------
  missing_tl = setdiff(tl_headers, colnames(lcms_table))
  if(length(missing_tl))
    warning("tl_headers missing from lcms_table: ", paste(missing_tl, collapse = ", "))
  for(h in c(bc_factor_name, bc_header, blank_head)){
    if(!h %in% colnames(metadata))
      warning(sprintf("'%s' not a metadata column.", h))
  }

  # --- Whitelists & pre-filter -------------------------------------------
  fnames = fdata$Processing_name[fdata$Report == "YES"]
  snames = metadata$Name[metadata$Include == "YES"]

  base = lcms_table %>%
    dplyr::filter(ID %in% fnames, Name %in% snames)

  # --- Wide datatype matrix + optional SNR mask --------------------------
  out_table = .build_wide_matrix(base, datatype)

  if(!isFALSE(snr))
    out_table = .apply_snr_mask(out_table, base, snr)

  if(length(fnames) != ncol(out_table))
    warning("Compound naming mismatch: check fdata$Processing_name vs lcms_table$ID.")

  # Rename Processing_name -> Compound via O(1) lookup
  name_lookup = stats::setNames(fdata$Compound, fdata$Processing_name)
  colnames(out_table) = name_lookup[colnames(out_table)]

  # Drop all-NA features (vectorized via colSums)
  keep_col = colSums(!is.na(out_table)) > 0
  remove_feats = names(out_table)[!keep_col]
  out_table    = out_table[, keep_col, drop = FALSE]

  # --- Per-batch processing ----------------------------------------------
  meta_yes = dplyr::filter(metadata, Include == "YES")
  batches  = unique(meta_yes[[bc_header]])

  if(length(batches) == 1L){
    # fast-path: skip the rebind
    out_matrix = .process_batch(
      out_table_batch = out_table[rownames(out_table) %in% meta_yes$Name, , drop = FALSE],
      metadata_batch  = meta_yes,
      blank_head      = blank_head,
      blank_name      = blank_name,
      blank_filter    = blank_filter,
      normalize       = normalize,
      adjust_conc     = adjust_conc,
      replace_MVs     = replace_MVs
    )
  } else {
    batch_frames = vector("list", length(batches))
    for(b_i in seq_along(batches)){
      metadata_batch = meta_yes[meta_yes[[bc_header]] == batches[b_i], , drop = FALSE]
      batch_frames[[b_i]] = .process_batch(
        out_table_batch = out_table[rownames(out_table) %in% metadata_batch$Name, , drop = FALSE],
        metadata_batch  = metadata_batch,
        blank_head      = blank_head,
        blank_name      = blank_name,
        blank_filter    = blank_filter,
        normalize       = normalize,
        adjust_conc     = adjust_conc,
        replace_MVs     = replace_MVs
      )
    }
    out_matrix = dplyr::bind_rows(batch_frames)
  }

  # --- Assemble DatasetExperiment ----------------------------------------
  lcms_experiment = .assemble_DE(out_matrix, fdata, metadata)

  if(isTRUE(batch_correction)){
    bc_wf = batch_correct(qc_label    = bc_qc_label,
                          factor_name = bc_factor_name,
                          batch_head  = bc_header)
    lcms_experiment = model_apply(bc_wf, lcms_experiment)@corrected@value
  }

  # --- Record removed features -------------------------------------------
  fdata$Report [fdata$Compound %in% remove_feats] = "NO"
  fdata$Comment[fdata$Compound %in% remove_feats] = "Fails S/NR in peak matrix processing"
  removed_features = dplyr::filter(fdata, Report == "NO")

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

#' Run blank-filter / normalise / concentration / MV-replace on one batch
#' @keywords internal
#' @noRd
.process_batch = function(out_table_batch, metadata_batch,
                          blank_head, blank_name,
                          blank_filter, normalize, adjust_conc, replace_MVs){

  blank_samples = metadata_batch$Name[metadata_batch[[blank_head]] == blank_name]

  if(!isFALSE(blank_filter))
    out_table_batch = blank_filter_apply(out_table_batch, blank_samples, blank_filter)

  if(!isFALSE(normalize))
    out_table_batch = out_table_batch / metadata_batch[[normalize]]

  if(isTRUE(adjust_conc))
    out_table_batch = calculate_conc(out_table_batch, metadata_batch)

  if(!isFALSE(replace_MVs))
    out_table_batch = replaceNA(out_table_batch, blank_samples, scalar = replace_MVs)

  out_table_batch
}

#' Assemble the final `struct::DatasetExperiment` with aligned metadata
#' @keywords internal
#' @noRd
.assemble_DE = function(out_matrix, fdata, metadata){
  fdata_output = fdata %>%
    dplyr::filter(Report == "YES", Compound %in% colnames(out_matrix)) %>%
    dplyr::arrange(match(Compound, colnames(out_matrix)))

  metadata_ar  = metadata %>% dplyr::arrange(Name)
  out_matrix_ar = out_matrix[metadata_ar$Name, , drop = FALSE]

  struct::DatasetExperiment(
    data          = out_matrix_ar,
    sample_meta   = metadata_ar,
    variable_meta = fdata_output
  )
}
