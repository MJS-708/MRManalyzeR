#' Load a `DatasetExperiment` from `.RDS` or `.xlsx`
#'
#' Used by [`combine_datasets()`] so the combine pipeline can mix RDS and
#' xlsx inputs (xlsx loading reconstructs a `DatasetExperiment` from the
#' standard sheets `feature_metadata`, `sample_metadata`, `matrix` written
#' by [`run_MRManalyzeR()`]).
#'
#' @param path Path to either an `.RDS` containing a
#'   `struct::DatasetExperiment`, or an `.xlsx` produced by
#'   [`run_MRManalyzeR()`] (with `feature_metadata`, `sample_metadata`,
#'   `matrix` tabs).
#' @param sample_id_col Column in `sample_metadata` used as the row-name
#'   key when loading from xlsx. Default `"Sample_ID"`.
#' @return A `struct::DatasetExperiment`.
#' @export
load_dataset = function(path, sample_id_col = "Sample_ID"){

  if(!file.exists(path)) stop("[load_dataset] file not found: ", path)

  ext = tolower(tools::file_ext(path))

  if(ext == "rds"){
    return(readRDS(path))
  }

  if(ext == "xlsx"){
    sheets = openxlsx::getSheetNames(path)
    need   = c("feature_metadata", "sample_metadata", "matrix")
    miss   = setdiff(need, sheets)
    if(length(miss))
      stop(sprintf("[load_dataset] %s is missing required sheet(s): %s",
                   path, paste(miss, collapse = ", ")))

    fmeta = openxlsx::read.xlsx(path, sheet = "feature_metadata")
    smeta = openxlsx::read.xlsx(path, sheet = "sample_metadata")
    M     = openxlsx::read.xlsx(path, sheet = "matrix", rowNames = TRUE)
    M     = as.matrix(M)

    if(!sample_id_col %in% colnames(smeta))
      stop(sprintf("[load_dataset] %s: sample_metadata lacks column '%s'.",
                   path, sample_id_col))

    # Align matrix rows to sample_meta in declared order, keyed by sample_id_col.
    sid = as.character(smeta[[sample_id_col]])
    if(all(sid %in% rownames(M))){
      M = M[sid, , drop = FALSE]
    } else if(nrow(M) == nrow(smeta)){
      # Fall back to row-order alignment.
      rownames(M) = sid
    } else {
      stop(sprintf(
        "[load_dataset] %s: cannot align matrix rows to sample_metadata (%d vs %d).",
        path, nrow(M), nrow(smeta)))
    }
    rownames(smeta) = sid

    # Variable_meta keyed by Compound (the convention used elsewhere in the package).
    if("Compound" %in% colnames(fmeta))
      rownames(fmeta) = as.character(fmeta$Compound)

    # Ensure variable_meta rows are in the same order as matrix columns —
    # SummarizedExperiment requires this and the two sheets aren't guaranteed
    # to be in lockstep.
    if(all(colnames(M) %in% rownames(fmeta))){
      fmeta = fmeta[colnames(M), , drop = FALSE]
    } else if(nrow(fmeta) == ncol(M)){
      rownames(fmeta) = colnames(M)
    } else {
      stop(sprintf(
        "[load_dataset] %s: cannot align feature_metadata to matrix columns (%d vs %d).",
        path, nrow(fmeta), ncol(M)))
    }

    return(struct::DatasetExperiment(
      name          = tools::file_path_sans_ext(basename(path)),
      data          = as.data.frame(M),
      sample_meta   = smeta,
      variable_meta = fmeta
    ))
  }

  stop("[load_dataset] unsupported extension '.", ext, "' for: ", path)
}
