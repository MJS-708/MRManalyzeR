#' MRManalyzeR: targeted LC-MS lipidomics and metabolomics processing
#'
#' Processes and analyses targeted lipidomics and metabolomics data exported from
#' Waters TargetLynx or Skyline: builds a peak-area / concentration matrix with
#' signal-to-noise or LOD/LOQ filtering, blank filtering, normalisation,
#' calibration/internal-standard concentration adjustment, missing-value
#' imputation and batch correction, then runs the statistical analyses and
#' renders two self-contained HTML reports.
#'
#' The workflow falls into four stages. Everything from stage 1 onwards is
#' carried in a single `struct::DatasetExperiment`, so each function takes that
#' object and returns it.
#'
#' @section Entry points:
#' [run_MRManalyzeR()] drives all four stages from a single YAML config;
#' [run_MRManalyzeR_combine()] merges several acquisition panels into one
#' analysis; [run_example()] runs the bundled example dataset end to end.
#'
#' @section 1. Data parse:
#' What did the instrument report? [read_targetlynx()] and [read_skyline()]
#' read a workbook into a wide sample x compound matrix, which
#' [assemble_dataset()] turns into a `struct::DatasetExperiment`.
#' [load_config()], [load_dataset()] and [combine_datasets()] read a YAML
#' config, a stored dataset, and several datasets merged into one.
#'
#' @section 2. Peak-matrix processing:
#' What is the true concentration? [filter_blanks()], [normalise_matrix()],
#' [adjust_concentration()], [impute_missing()] and [correct_batch()], with
#' [subset_dataset()] to slice a dataset by sample metadata.
#' [process_dataset()] composes stages 1 and 2 in one call.
#'
#' @section 3. QC check:
#' Can these numbers be trusted? [run_pca()], plus the per-compound CV metrics
#' added to `variable_meta` by [run_MRManalyzeR()] and the bundled
#' data-quality report.
#'
#' @section 4. Stats:
#' What do these numbers mean? [run_stats()] runs the comparisons,
#' correlations and linear models described by the YAML `stats_report:` block
#' and returns the same tables the statistics report and the output xlsx use.
#'
#' The `plot_*()` functions are exported rather than hidden inside the report
#' templates so that a figure can be reproduced, restyled or subset in one
#' line, without editing an R Markdown file or rebuilding it from the matrix.
#' The bundled reports call the same functions.
#'
#' See `vignette("MRManalyzeR")` for a step-by-step walkthrough.
#'
#' @keywords internal
#' @importFrom rlang .data
#' @importFrom utils head
#' @importFrom stats median
#' @importFrom methods new
#' @importFrom struct model_apply
"_PACKAGE"

# dplyr NSE column references used across the readers / assemble step.
utils::globalVariables(c("ID", "ID2", "Name", "S/N", "Report", "Compound"))

#' Shared plot theme
#'
#' One theme for every plot the package draws, so the two HTML reports and
#' anything a user builds by calling the `plot_*()` functions directly look
#' like the same suite.
#'
#' @return A `ggplot2` theme object.
#' @keywords internal
#' @noRd
.mrm_theme = function(){
  ggplot2::theme_classic() +
    ggplot2::theme(
      legend.position = "bottom",
      axis.title = ggplot2::element_text(face = "bold", color = "black",
                                         size = 12),
      axis.text  = ggplot2::element_text(face = "bold", color = "black",
                                         size = 8))
}

#' Vendored ltc colour palettes
#'
#' Hex values taken from the `ltc` palettes of Loukas Theodosiou
#' (https://github.com/loukesio/ltc-color-palettes). They are inlined rather
#' than depended on because `ltc` is distributed only on GitHub, and a
#' Bioconductor package cannot require a remote install.
#'
#' `hat` is the qualitative palette used for categorical annotations - feature
#' classes, sample types, batches - and `heatmap2` the diverging scheme used
#' for z-scores and correlations.
#'
#' @keywords internal
#' @noRd
.ltc = list(
  # Reordered from the ltc original so the strong, maximally separable hues
  # come first - blue, red, green, yellow, black, purple - and the off-shades
  # only appear once a variable has more than six levels. A two- or
  # three-level factor should not be drawn in two shades of amber.
  #
  # The off-shades are themselves ordered so that each one sits as far as
  # possible from the strong hue it resembles: teal first (nothing before it
  # is close), then crimson (6 places after red), light green (6 after green)
  # and orange last (6 after yellow). At eight levels the palette therefore
  # reads blue/red/green/yellow/black/purple/teal/crimson rather than pairing
  # yellow with orange in the same legend.
  hat = c("#4e54ac", "#e8351e", "#17a769", "#efb306", "#000000", "#852f88",
          "#0f8096", "#cd023d", "#7db954", "#eb990c"),
  reading = c("#EFBC68", "#919F89", "#EDBDAE", "#57717C",
              "#5F97A4", "#CAEAC8", "#95A1AE", "#C8CFD6"),
  heatmap0 = c("#001219", "#005F73", "#0A9396", "#94D2BD", "#E9D8A6",
               "#EE9B00", "#CA6702", "#AE2012", "#9B2226"),
  heatmap2 = c("#ca0020", "#f4a582", "#f7f7f7", "#92c5de", "#0571b0")
)

#' The palette currently selected for categorical variables
#'
#' Read from `getOption("MRManalyzeR.palette")` so a report can set it once
#' from YAML rather than every plotting function taking a palette argument.
#'
#' @return Character vector of hex colours.
#' @keywords internal
#' @noRd
.mrm_palette = function(){
  nm = getOption("MRManalyzeR.palette", "hat")
  .ltc[[nm]] %||% .ltc$hat
}

#' Diverging colour ramp for z-scores and correlations
#'
#' Reversed so it runs low-blue to high-red, the direction these figures are
#' read. The scheme is `getOption("MRManalyzeR.diverging")`, default
#' `heatmap2`.
#'
#' @param n Number of steps.
#' @return Character vector of `n` hex colours.
#' @keywords internal
#' @noRd
.diverging_pal = function(n = 101){
  nm  = getOption("MRManalyzeR.diverging", "heatmap2")
  pal = .ltc[[nm]] %||% .ltc$heatmap2
  grDevices::colorRampPalette(rev(pal))(n)
}

#' Qualitative palette for n categorical levels
#'
#' The vendored `ltc` `hat` palette, which carries ten distinguishable hues -
#' more than Brewer's Set1 - and is interpolated beyond that.
#'
#' @param n Number of levels.
#' @return Character vector of `n` colours.
#' @keywords internal
#' @noRd
.qual_pal = function(n){
  if(n < 1) return(character(0))
  pal = .mrm_palette()
  if(n <= length(pal)) return(pal[seq_len(n)])
  grDevices::colorRampPalette(pal)(n)
}

#' Named colour vector for the sorted levels of a categorical variable
#'
#' Keyed on sorted levels so a given level maps to the same colour wherever it
#' appears - the heatmap feature annotation and the PCA loadings bars must
#' agree, or the two figures cannot be read against each other. `NA` maps to
#' grey; a two-level variable uses a fixed blue/red pair, which reads as an
#' ordered contrast where an arbitrary pair would not.
#'
#' @param x Vector of level labels.
#' @return Named character vector of colours.
#' @keywords internal
#' @noRd
.group_fill_colors = function(x){
  x  = as.character(x)
  lv = sort(unique(x[!is.na(x) & x != "NA"]))
  m  = if(length(lv) == 0)      character(0)
       else if(length(lv) == 2) stats::setNames(
                                  .mrm_palette()[c(1L, 2L)], lv)
       else                     stats::setNames(.qual_pal(length(lv)), lv)
  if(any(is.na(x) | x == "NA")) m = c(m, "NA" = "grey70")
  m
}

#' Annotation colour list for pheatmap
#'
#' One named colour vector per annotation column, via
#' [.group_fill_colors()] so a class keeps its colour across the heatmap
#' annotation, the PCA loadings and the correlation strips.
#'
#' @param ann Annotation data frame, or `NULL`.
#' @return A named list of colour vectors, or `NULL`.
#' @keywords internal
#' @noRd
.ann_colours = function(ann){
  if(is.null(ann)) return(NULL)
  out = lapply(ann, .group_fill_colors)
  names(out) = colnames(ann)
  out
}
