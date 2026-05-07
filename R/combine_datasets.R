#' Combine multiple processed datasets into one
#'
#' Used to merge separately-acquired LC-MS panels (e.g. GOM, cysLT, SPM) that
#' share a common set of biological samples. Inputs can be either `.RDS`
#' files (containing a `struct::DatasetExperiment`) or `.xlsx` workbooks
#' produced by [`run_MRManalyzeR()`] — see [`load_dataset()`].
#'
#' Defaults:
#' \itemize{
#'   \item Samples matched by `Sample_ID` (intersect across inputs).
#'   \item `feature_meta_cols` / `sample_meta_cols` default to the
#'         intersection of column names across all inputs (so you only have
#'         to spell out an explicit list when you want to *narrow* the set).
#'   \item Data matrices are `cbind`'d on the common samples.
#' }
#'
#' @param paths Character vector of paths to `.RDS` or `.xlsx` files. May be
#'   named (names are used as dataset tags); otherwise the filename stem
#'   is used as the tag.
#' @param feature_meta_cols Character vector — `variable_meta` columns kept
#'   in the merged dataset. `NULL` (default) = intersect across inputs.
#'   Missing columns in any single input are filled with `NA`.
#' @param sample_meta_cols Character vector — `sample_meta` columns kept in
#'   the merged dataset. `NULL` (default) = intersect across inputs.
#' @param feature_meta_rename Optional named list. Names are dataset tags
#'   (or full paths); values are named character vectors of
#'   `c(new_name = "old_name")` used to rename `variable_meta` columns
#'   before alignment.
#' @param qc_remap Optional named list. Names are dataset tags (or paths);
#'   values are named character vectors of `c("OldQCID" = "NewQCID")`
#'   applied to `Sample_ID`. QC samples whose original ID is NOT a key in
#'   the map are dropped.
#' @param qc_regex Regex matching QC sample IDs (default `"^QC"`).
#' @param sample_id_col Column in `sample_meta` used to match samples
#'   across datasets. Default `"Sample_ID"`.
#' @param drop_samples Character vector of `Sample_ID`s to exclude from the
#'   merged result.
#' @param prefix_features Controls feature-name prefixing. `FALSE` (default):
#'   never prefix. `TRUE`: always prefix every feature with its dataset tag
#'   (e.g. `GOM__PGE2`). `"auto"`: prefix only features whose names collide
#'   across two or more inputs — note that the resulting names depend on which
#'   datasets are combined, so `"auto"` is not stable across different combine
#'   runs.
#' @param combined_name Name attribute of the resulting `DatasetExperiment`.
#' @return A single `struct::DatasetExperiment` containing the merged data.
#' @export
combine_datasets = function(paths,
                            feature_meta_cols   = NULL,
                            sample_meta_cols    = NULL,
                            feature_meta_rename = NULL,
                            qc_remap            = NULL,
                            qc_regex            = "^QC",
                            sample_id_col       = "Sample_ID",
                            drop_samples        = NULL,
                            prefix_features     = FALSE,   # FALSE | TRUE | "auto"
                            combined_name       = "combined"){

  stopifnot(length(paths) >= 1)
  stopifnot(all(file.exists(paths)))
  if(!identical(prefix_features, FALSE) &&
     !identical(prefix_features, TRUE)  &&
     !identical(prefix_features, "auto"))
    stop("[combine_datasets] prefix_features must be FALSE, TRUE, or \"auto\".")

  # Tags used for feature prefixing + lookup in the per-dataset config lists.
  tags = if(!is.null(names(paths)) && all(nzchar(names(paths))))
            names(paths)
         else tools::file_path_sans_ext(basename(paths))

  # Storage: extract to plain data.frames immediately so we can modify them
  # freely without hitting S4 slot-assignment validation on DatasetExperiment.
  dats   = vector("list", length(paths))   # samples × features
  smetas = vector("list", length(paths))   # sample metadata
  vmetas = vector("list", length(paths))   # variable (feature) metadata

  for(i in seq_along(paths)){
    p   = paths[i]
    tag = tags[i]
    de  = load_dataset(p, sample_id_col = sample_id_col)

    if(!sample_id_col %in% colnames(de$sample_meta))
      stop(sprintf("[combine_datasets] '%s' lacks column '%s' in sample_meta.",
                   p, sample_id_col))

    # Extract as plain data.frames — S4 accessors return copies anyway, and
    # we need to mutate dimnames freely before final assembly.
    dat   = as.data.frame(de$data)
    smeta = as.data.frame(de$sample_meta)
    vmeta = as.data.frame(de$variable_meta)

    # 1. Deduplication: any sample type can have a repeated Sample_ID (repeated
    # QC injections, data-entry errors, etc.). Collapse to first occurrence and
    # log which IDs were affected. The combine only needs each ID once; blanks
    # and other non-biological samples are filtered later by biological_filter.
    sid_tmp = as.character(smeta[[sample_id_col]])
    dup_any = duplicated(sid_tmp)
    if(any(dup_any)){
      dup_ids = unique(sid_tmp[dup_any])
      message(sprintf(
        "[combine] '%s': collapsed %d duplicate Sample_ID(s) to first occurrence: %s",
        tag, length(dup_ids), paste(dup_ids, collapse = ", ")))
      keep    = !dup_any
      smeta   = smeta[keep, , drop = FALSE]
      dat     = dat[keep,   , drop = FALSE]
      sid_tmp = sid_tmp[keep]
    }

    # 2. Optional QC remap: rename QC IDs and drop unmapped QCs (per-dataset).
    qm = .pick_per_dataset(qc_remap, p, tag)
    if(length(qm)){
      is_qc = grepl(qc_regex, sid_tmp)
      keep  = (!is_qc) | (sid_tmp %in% names(qm))
      smeta = smeta[keep, , drop = FALSE]
      dat   = dat[keep,   , drop = FALSE]
      sid_tmp = sid_tmp[keep]
      recoded = unname(ifelse(sid_tmp %in% names(qm), qm[sid_tmp], sid_tmp))
      if(anyDuplicated(recoded))
        stop("[combine_datasets] QC remap produced duplicate Sample_IDs: ",
             paste(unique(recoded[duplicated(recoded)]), collapse = ", "))
      smeta[[sample_id_col]] = recoded
      rownames(smeta)        = recoded
      rownames(dat)          = recoded
    }

    # 2. Optional variable_meta column renames (per-dataset)
    vmr = .pick_per_dataset(feature_meta_rename, p, tag)
    if(length(vmr)) vmeta = .rename_cols(vmeta, vmr)

    # 4. Canonicalise sample rownames to the chosen ID column
    sid = as.character(smeta[[sample_id_col]])
    if(anyNA(sid) || any(!nzchar(sid)) || anyDuplicated(sid))
      stop(sprintf("[combine_datasets] '%s' has missing/blank/duplicate %s.",
                   p, sample_id_col))
    rownames(smeta) = sid
    rownames(dat)   = sid

    dats[[i]]   = dat
    smetas[[i]] = smeta
    vmetas[[i]] = vmeta
  }

  # 3b. Apply prefixing now that all datasets are loaded.
  #     "auto": only prefix features that appear in more than one dataset.
  if(!isFALSE(prefix_features)){
    all_feats   = unlist(lapply(dats, colnames))
    dup_feats   = unique(all_feats[duplicated(all_feats)])
    for(i in seq_along(dats)){
      tag       = tags[i]
      old_feats = colnames(dats[[i]])
      to_prefix = if(isTRUE(prefix_features)) rep(TRUE, length(old_feats))
                  else old_feats %in% dup_feats
      new_feats = ifelse(to_prefix, paste0(tag, "__", old_feats), old_feats)

      # Update colnames(dat) and rownames(vmeta) every iteration so they
      # are guaranteed to be in lockstep, even when no prefixing happens.
      colnames(dats[[i]])    = new_feats
      rownames(vmetas[[i]])  = new_feats
      if("Compound" %in% colnames(vmetas[[i]])){
        # Coerce to character so factor inputs don't silently become NA when
        # we assign a new (prefixed) level.
        vmetas[[i]]$Compound = as.character(vmetas[[i]]$Compound)
        vmetas[[i]]$Compound[to_prefix] = new_feats[to_prefix]
      }
    }
  }

  # 5. Resolve column lists — intersect-by-default
  if(is.null(feature_meta_cols))
    feature_meta_cols = Reduce(intersect, lapply(vmetas, colnames))
  if(is.null(sample_meta_cols))
    sample_meta_cols  = Reduce(intersect, lapply(smetas, colnames))
  if(length(feature_meta_cols) == 0)
    stop("[combine_datasets] No common variable_meta columns across inputs ",
         "(and no `feature_meta_cols` provided).")
  if(length(sample_meta_cols) == 0)
    stop("[combine_datasets] No common sample_meta columns across inputs ",
         "(and no `sample_meta_cols` provided).")
  if(!sample_id_col %in% sample_meta_cols)
    sample_meta_cols = c(sample_id_col, sample_meta_cols)

  # 6. Align variable_meta to chosen columns (per-dataset; missing → NA)
  vmetas = lapply(vmetas, function(vm) .align_cols(vm, feature_meta_cols))

  # 7. Intersect samples by Sample_ID
  common = Reduce(intersect, lapply(dats, rownames))
  if(length(common) == 0)
    stop("[combine_datasets] No samples in common across the input datasets.")
  if(!is.null(drop_samples))
    common = setdiff(common, drop_samples)

  # 8. Assemble combined objects. All are now plain data.frames so the
  #    DatasetExperiment constructor gets consistent dimnames.
  data_combined  = do.call(cbind, lapply(dats,   function(d) d[common, , drop = FALSE]))
  vmeta_combined = do.call(rbind, vmetas)
  smeta_combined = smetas[[1]][common, sample_meta_cols, drop = FALSE]

  # Sanity check: SummarizedExperiment requires assay dimnames to line up
  # with rowData/colData rownames. Catch mismatches here with a useful
  # message rather than the cryptic SE constructor error.
  if(!identical(colnames(data_combined), rownames(vmeta_combined))){
    bad = which(colnames(data_combined) != rownames(vmeta_combined))
    stop(sprintf(
      "[combine_datasets] Internal alignment failure: %d feature(s) have ",
      length(bad)),
      "mismatching colnames(data) vs rownames(variable_meta). ",
      "First mismatch: data='", colnames(data_combined)[bad[1]],
      "' vs vmeta='", rownames(vmeta_combined)[bad[1]], "'.")
  }
  if(anyDuplicated(colnames(data_combined))){
    dups = unique(colnames(data_combined)[duplicated(colnames(data_combined))])
    stop("[combine_datasets] Combined data has duplicate feature names: ",
         paste(head(dups, 10), collapse = ", "),
         if(length(dups) > 10) sprintf(" ... (%d total)", length(dups)) else "",
         ". Use prefix_features=TRUE (or \"auto\") to disambiguate.")
  }

  struct::DatasetExperiment(
    name          = combined_name,
    data          = data_combined,
    sample_meta   = smeta_combined,
    variable_meta = vmeta_combined
  )
}


# ---- internal helpers ------------------------------------------------------

#' @keywords internal
#' @noRd
.pick_per_dataset = function(map, path, tag){
  if(is.null(map)) return(NULL)
  for(key in c(path, tag, basename(path))){
    if(!is.null(map[[key]])) return(map[[key]])
  }
  NULL
}

#' Rename data.frame columns. `rename_map` is `c(new = "old")`.
#' @keywords internal
#' @noRd
.rename_cols = function(df, rename_map){
  olds = intersect(unname(rename_map), colnames(df))
  if(length(olds)){
    news = names(rename_map)[match(olds, unname(rename_map))]
    colnames(df)[match(olds, colnames(df))] = news
  }
  df
}

#' Add missing columns as NA, keep only the requested set in the requested order.
#' @keywords internal
#' @noRd
.align_cols = function(df, cols){
  miss = setdiff(cols, colnames(df))
  for(m in miss) df[[m]] = NA
  df[, cols, drop = FALSE]
}
