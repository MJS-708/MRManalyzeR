#' Adjust a matrix to sample concentrations
#'
#' Applies the vial-level concentration correction (cal_vol/sample_vol and the
#' IS ratio), and -- when `starting_vol_col` is set --
#' chains the starting-volume correction back to the concentration in the
#' original (undried) sample.
#'
#' @param x Wide numeric matrix / data frame, rows = samples.
#' @param sample_meta Sample metadata with the volume / IS columns.
#' @param sample_vol_col,cal_vol_col,sample_IS_col,cal_IS_col `sample_meta`
#'   columns for the vial-level correction.
#' @param starting_vol_col `sample_meta` column holding the pre-dry-down
#'   (starting) volume, or `FALSE` to skip the starting-volume step.
#' @param name_col `sample_meta` column matching the row names of `x`.
#' @return `x` scaled to concentration per sample.
#' @examples
#' m  <- data.frame(A = c(10, 20), B = c(30, 40), row.names = c("S1", "S2"))
#' sm <- data.frame(Name = c("S1", "S2"),
#'                  sample_volume_uL = c(100, 100), cal_vol_uL = c(100, 100),
#'                  sample_IS_vol_uL = c(10, 10),   cal_IS_vol_uL = c(10, 10))
#' adjust_concentration(m, sm)
#' @export
adjust_concentration = function(x, sample_meta,
                                sample_vol_col   = "sample_volume_uL",
                                cal_vol_col      = "cal_vol_uL",
                                sample_IS_col    = "sample_IS_vol_uL",
                                cal_IS_col       = "cal_IS_vol_uL",
                                starting_vol_col = FALSE,
                                name_col         = "Name"){
  x = calculate_conc(x, sample_meta,
                     sample_vol_col = sample_vol_col,
                     cal_vol_col    = cal_vol_col,
                     sample_IS_col  = sample_IS_col,
                     cal_IS_col     = cal_IS_col,
                     name_col       = name_col)
  if(!isFALSE(starting_vol_col))
    x = .apply_starting_vol_correction(x, sample_meta,
                                       starting_vol_col = starting_vol_col,
                                       sample_vol_col   = sample_vol_col,
                                       name_col         = name_col)
  x
}


# Calculate the concentration of lipid in each sample.
# conc from cal curve (ng/mL) assumes the sample was prepared identically to
# the calibration standards (same total volume, same IS amount). Correct
# per-sample for any deviation from that:
#   conc = conc * (cal_vol / sample_vol) * (sample_IS / cal_IS)
# e.g. a larger sample_vol (same analyte diluted into more volume) lowers
# the corrected conc; more sample_IS (diluting the analyte/IS ratio the
# curve assumed) raises it back up to compensate.
#
# Looked up per-sample BY NAME (not row position): df's row order comes
# from the source-data reader (TargetLynx pivot / Skyline sheet order),
# metadata_batch's comes from the sample_metadata sheet -- these are not
# guaranteed to match, so positional indexing here would silently apply
# the wrong sample's volumes when df and metadata_batch are ordered
# differently.
#
# Column names are parameters, not hardcoded: different studies name these
# sample_metadata columns differently (see PeakMatrixProcessing.adjust_conc
# in the project YAML), so the defaults below are just the common case.

calculate_conc = function(df, metadata_batch,
                          sample_vol_col = "sample_volume_uL",
                          cal_vol_col    = "cal_vol_uL",
                          sample_IS_col  = "sample_IS_vol_uL",
                          cal_IS_col     = "cal_IS_vol_uL",
                          name_col       = "Name"){

  req_cols = c(sample_vol_col, cal_vol_col, sample_IS_col, cal_IS_col)
  missing_cols = setdiff(req_cols, colnames(metadata_batch))
  if(length(missing_cols))
    stop(sprintf("[calculate_conc] metadata_batch is missing required column(s): %s",
                paste(missing_cols, collapse = ", ")))

  idx = match(rownames(df), metadata_batch[[name_col]])
  if(anyNA(idx))
    stop(sprintf(
      "[calculate_conc] %d sample(s) in the data matrix have no matching row in metadata_batch[['%s']]: %s",
      sum(is.na(idx)), name_col, paste(head(rownames(df)[is.na(idx)], 5), collapse = ", ")))

  vol_fact = metadata_batch[[cal_vol_col]][idx]   / metadata_batch[[sample_vol_col]][idx]
  IS_fact  = metadata_batch[[sample_IS_col]][idx] / metadata_batch[[cal_IS_col]][idx]
  fact     = vol_fact * IS_fact

  n_missing = sum(is.na(fact))
  if(n_missing)
    message(sprintf(
      "[calculate_conc] %d sample(s) have missing volume/IS metadata and will be NA after conc adjustment: %s",
      n_missing, paste(head(rownames(df)[is.na(fact)], 5), collapse = ", ")))

  df * fact   # per-row factor recycles down each column -> each sample's own factor applies to all its compounds
}

# calculate_conc() above corrects to the concentration IN THE LC-MS VIAL
# (as if prepared like the calibration standards). This is a further,
# separate correction: the vial solution is a dried-down and reconstituted
# aliquot of `starting_vol_col` of original sample material, brought up to
# `sample_vol_col` for injection. Mass is conserved through that
# dry-down/reconstitution step, so:
#   conc_original = conc_vial * (sample_vol / starting_vol)
#                 = conc_vial / (starting_vol / sample_vol)
# e.g. more starting material dried into the same reconstitution volume
# raises conc_vial, but this correction scales it back down proportionally,
# recovering the concentration in the original (undried) sample.

.apply_starting_vol_correction = function(df, metadata_batch, starting_vol_col,
                                          sample_vol_col = "sample_volume_uL",
                                          name_col       = "Name"){

  req_cols = c(sample_vol_col, starting_vol_col)
  missing_cols = setdiff(req_cols, colnames(metadata_batch))
  if(length(missing_cols))
    stop(sprintf("[starting_vol_correction] metadata_batch is missing required column(s): %s",
                paste(missing_cols, collapse = ", ")))

  idx = match(rownames(df), metadata_batch[[name_col]])
  if(anyNA(idx))
    stop(sprintf(
      "[starting_vol_correction] %d sample(s) in the data matrix have no matching row in metadata_batch[['%s']]: %s",
      sum(is.na(idx)), name_col, paste(head(rownames(df)[is.na(idx)], 5), collapse = ", ")))

  fact = metadata_batch[[sample_vol_col]][idx] / metadata_batch[[starting_vol_col]][idx]

  n_missing = sum(is.na(fact))
  if(n_missing)
    message(sprintf(
      "[starting_vol_correction] %d sample(s) have missing %s/%s and will be NA after starting-volume correction: %s",
      n_missing, sample_vol_col, starting_vol_col,
      paste(head(rownames(df)[is.na(fact)], 5), collapse = ", ")))

  df * fact
}
