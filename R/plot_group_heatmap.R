#' Group-mean heatmap that shows which samples produced the mean
#'
#' A heatmap of groups against compounds where each cell carries two things at
#' once: the **hue** is the group's mean scaled value for that compound, and
#' the **opacity** of each sample's sub-cell is how much that sample
#' contributed to the mean. A block that is uniformly solid is an effect the
#' whole group shares; a block carrying one opaque cell among faded ones is a
#' mean resting on a single sample - which is exactly what a plain group-mean
#' heatmap hides, and what a reader needs in order to judge whether a colour
#' contrast is worth believing.
#'
#' [plot_heatmap()] draws every sample and never averages, so it does not have
#' this problem; use it while the sample axis is still legible. This function
#' is for when it is not - a large cohort, or a panel wide enough that one
#' column per sample stops fitting - where the alternative is a group mean that
#' silently discards the replicate structure.
#'
#' @section How a cell is built:
#' Values are transformed, then scaled per compound across **all** samples, so
#' a cell reads "high or low for this compound" rather than "abundant or not" -
#' without it the colour range is set by the few most abundant compounds and
#' everything else is a flat wash. Rows are therefore not comparable to each
#' other in absolute terms.
#'
#' The hue is the mean of those scores within a group. The opacity is the
#' sample's own score multiplied by the sign of its group mean, so a sample
#' that moved with its group is opaque and one that sat near zero, or pulled
#' the other way, fades. `alpha_floor` stops a non-contributing sample becoming
#' invisible and reading as missing data; genuinely missing values are drawn in
#' `na.value` grey at full opacity, so the two cannot be confused.
#'
#' Note that "contributed" is measured against the mean of **all** samples, not
#' against the group's own members, so it is not the same as "has the largest
#' value". A group of four holding one very high sample and three low ones has
#' a *negative* mean if the three sit further from the overall centre than the
#' one sits above it - and then the three are what produced that mean, and the
#' high sample is the one opposing it and fading out. That is the intended
#' reading: opacity answers "did this sample push its group where the group
#' ended up", not "is this sample large".
#'
#' @section Choosing the colour limit:
#' `cap` defaults to a quantile of the group means actually being drawn, and
#' that default matters more than it looks. Averaging *n* replicates shrinks a
#' scaled value by roughly `sqrt(n)`, so a limit chosen for per-sample scores -
#' the `+/-2` that suits [plot_heatmap()] - leaves a group-mean panel using a
#' small part of the colour range and looking washed out, whatever the opacity
#' is set to. Take the limit from the quantity being drawn rather than
#' inheriting it from the per-sample view.
#'
#' @param de A `struct::DatasetExperiment`.
#' @param group_by `sample_meta` column defining the groups to average within.
#' @param features Compounds to include; `NULL` uses all of them.
#' @param group_features_by `variable_meta` column used to block compounds into
#'   labelled bands; `NULL` leaves them unblocked.
#' @param transform `"none"`, `"log2"` or `"sqrt"`, applied before scaling.
#'   Negative values are floored at 0 first.
#' @param cap Colour limit, in standard deviations of the group means. `NULL`
#'   takes the 97th percentile of `abs()` of the group means, so a few cells
#'   saturate and the rest get the full range.
#' @param alpha_floor Opacity given to a sample that did not contribute to its
#'   group's mean. Below about 0.4 those cells read as missing data.
#' @param feature_labels `"show"`, `"small"` or `"hide"`.
#' @param sample_labels Draw the sample axis. `FALSE` for wide cohorts, where
#'   the samples stay individually visible but not individually named.
#' @param sample_id_head `sample_meta` column supplying sample labels; `NULL`
#'   uses row names.
#' @param title Plot title.
#' @return A `ggplot`, or `NULL` when there is too little data to draw (fewer
#'   than 2 samples, or fewer than 2 compounds with any variance).
#' @examples
#' de <- load_dataset(system.file("extdata", "example_synthetic.RDS",
#'                                package = "MRManalyzeR"))
#' de <- subset_dataset(de, conditions = list(Sample_type = "Sample"))
#' plot_group_heatmap(de, group_by = "Treatment",
#'                    group_features_by = "Enzymatic_pathway")
#' @family stats
#' @seealso [plot_heatmap()], which draws every sample without averaging.
#' @export
plot_group_heatmap = function(de, group_by,
                              features          = NULL,
                              group_features_by = NULL,
                              transform         = c("none", "log2", "sqrt"),
                              cap               = NULL,
                              alpha_floor       = 0.6,
                              feature_labels    = c("show", "small", "hide"),
                              sample_labels     = TRUE,
                              sample_id_head    = NULL,
                              title             = NULL){

  transform      = match.arg(transform)
  feature_labels = match.arg(feature_labels)

  smeta = as.data.frame(de$sample_meta)
  vmeta = as.data.frame(de$variable_meta)
  M     = as.matrix(as.data.frame(de$data))
  if(!is.null(features)) M = M[, intersect(features, colnames(M)),
                               drop = FALSE]
  if(nrow(M) < 2 || ncol(M) < 2) return(NULL)

  group_by = .resolve_meta_col(group_by, smeta) %||%
    stop(sprintf(
      "[plot_group_heatmap] '%s' is not a sample_meta column. Present: %s.",
      group_by, paste(colnames(smeta), collapse = ", ")))
  group_features_by = .resolve_meta_col(group_features_by, vmeta)
  sample_id_head    = .resolve_meta_col(sample_id_head, smeta)

  if(!is.null(sample_id_head))
    rownames(M) = as.character(smeta[[sample_id_head]])

  if(transform == "log2"){
    M[M < 0 & !is.na(M)] = 0; M = log2(M + 1)
  } else if(transform == "sqrt"){
    M[M < 0 & !is.na(M)] = 0; M = sqrt(M)
  }

  # A zero-variance compound has no scaled value at all - it would be a row of
  # NaN rather than a row of zeros - so drop it instead of drawing it grey.
  vars = apply(M, 2, stats::var, na.rm = TRUE)
  M    = M[, is.finite(vars) & vars > 0, drop = FALSE]
  if(ncol(M) < 2) return(NULL)

  # Scaled per compound across ALL samples, not within group: scaling within
  # group would remove the very differences the figure is about.
  Z   = scale(M, center = TRUE, scale = TRUE)
  grp = factor(as.character(smeta[[group_by]]))

  # Group means: compounds x groups.
  gm = vapply(levels(grp),
              function(g) colMeans(Z[grp == g, , drop = FALSE], na.rm = TRUE),
              numeric(ncol(Z)))
  dimnames(gm) = list(colnames(Z), levels(grp))

  # Each sample's contribution to its own group's mean: its score signed
  # towards that mean. A group mean of exactly zero gives sign() == 0, which
  # correctly leaves every sample in it with no contribution to make. Written
  # as a loop over groups - a handful of iterations, and far more legible than
  # the equivalent matrix indexing.
  contrib = Z
  for(g in levels(grp)){
    i = which(grp == g)
    contrib[i, ] = sweep(Z[i, , drop = FALSE], 2, sign(gm[, g]), `*`)
  }

  # Reference for the opacity scale: how large a contribution counts as full
  # strength. Taken from the data, so it does not need tuning per dataset.
  a_ref = stats::quantile(abs(Z), 0.95, na.rm = TRUE)
  if(!is.finite(a_ref) || a_ref <= 0) a_ref = 1
  alpha_floor = min(max(alpha_floor, 0), 1)
  alpha = pmin(pmax(contrib, 0) / as.numeric(a_ref), 1)
  alpha = alpha_floor + (1 - alpha_floor) * alpha

  if(is.null(cap)) cap = stats::quantile(abs(gm), 0.97, na.rm = TRUE)
  cap = as.numeric(cap)
  if(!is.finite(cap) || cap <= 0) cap = 1

  # Compounds keep their column order within a block, and blocks follow the
  # order the classes first appear - so the caller controls both by ordering
  # variable_meta.
  feats = colnames(Z)
  cls   = if(is.null(group_features_by)) rep(NA_character_, length(feats))
          else as.character(vmeta[[group_features_by]][
            match(feats, vmeta$Compound)])
  cls[is.na(cls)] = "unassigned"
  ord   = order(match(cls, unique(cls)))
  feats = feats[ord]
  cls   = cls[ord]

  # as.vector() on a matrix is column-major, so the sample index varies fastest
  # - which is what rep(times=) / rep(each=) below assume.
  d = data.frame(
    sample  = rep(rownames(Z),        times = length(feats)),
    feature = rep(feats,              each  = nrow(Z)),
    group   = rep(as.character(grp),  times = length(feats)),
    alpha   = as.vector(alpha[, feats, drop = FALSE]),
    stringsAsFactors = FALSE)
  d$class = cls[match(d$feature, feats)]

  # Clamped rather than left to fall outside the scale limits, which would draw
  # the most extreme cells as if they were missing.
  d$fill = pmin(pmax(gm[cbind(d$feature, d$group)], -cap), cap)

  # A missing value is grey at full opacity, so it cannot be mistaken for a
  # sample that simply did not contribute.
  na_cell = as.vector(is.na(Z[, feats, drop = FALSE]))
  d$fill[na_cell]  = NA_real_
  d$alpha[na_cell] = 1

  d$feature = factor(d$feature, levels = rev(feats))
  d$sample  = factor(d$sample,  levels = rownames(Z)[order(grp)])
  d$group   = factor(d$group,   levels = levels(grp))
  d$class   = factor(d$class,   levels = unique(cls))

  fs = switch(feature_labels, show = 9, small = 6, 6)

  p = ggplot2::ggplot(d, ggplot2::aes(.data$sample, .data$feature,
                                      fill  = .data$fill,
                                      alpha = .data$alpha)) +
    ggplot2::geom_tile(color = "white", linewidth = 0.2) +
    ggplot2::scale_fill_gradientn(colours = .diverging_pal(101),
                                  limits = c(-cap, cap), na.value = "grey88",
                                  name = "group mean") +
    ggplot2::scale_alpha_identity() +
    ggplot2::labs(x = NULL, y = NULL, title = title) +
    .mrm_theme() +
    ggplot2::theme(
      axis.text.x = if(isTRUE(sample_labels))
        ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1, size = 6)
      else ggplot2::element_blank(),
      axis.text.y = if(feature_labels == "hide") ggplot2::element_blank()
                    else ggplot2::element_text(size = fs),
      axis.ticks    = ggplot2::element_blank(),
      axis.line     = ggplot2::element_blank(),
      panel.spacing = ggplot2::unit(0.15, "lines"))

  # Facets rather than annotation strips: the group and the compound class are
  # then named in words, in blocks with gaps between them, and no second colour
  # scale is needed to compete with the one carrying the data.
  p + if(is.null(group_features_by))
    ggplot2::facet_grid(cols = ggplot2::vars(.data$group),
                        scales = "free_x", space = "free_x")
  else
    ggplot2::facet_grid(rows = ggplot2::vars(.data$class),
                        cols = ggplot2::vars(.data$group),
                        scales = "free", space = "free")
}
