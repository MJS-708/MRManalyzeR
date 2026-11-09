# The plot layer: each function returns a ggplot built from the data it was
# given, refuses input it cannot use, and returns NULL rather than erroring
# where there is nothing to draw. The bundled synthetic dataset is used so
# the plots are built from a realistic shape.

.de_plots <- function() {
  de <- loadDataset(system.file("extdata", "example_synthetic.RDS",
                                package = "MRManalyzeR"))
  subsetDataset(de, conditions = list(Sample_type = "Sample"))
}

.stats_plots <- function() {
  runStats(.de_plots(), list(
    sig_threshold = 0.05,
    comparisons = list(enabled = TRUE, entries = list(
      list(name = "PBS_vs_HDM", method = "welch",
           compare = list(factor = "Treatment", levels = c("PBS", "HDM"))))),
    correlations = list(enabled = TRUE, entries = list(
      list(name = "pairs", methods = "pearson", subsets = list(list()))))))
}

test_that("plotBoxplot draws a compound by group, as box or bar", {
  de <- .de_plots()
  # The points carry a `text` aesthetic for the plotly tooltips in the
  # reports, which ggplot2 does not know and warns about on every build.
  expect_s3_class(suppressWarnings(
    plotBoxplot(de, "Analyte_02", "Treatment")), "ggplot")
  expect_s3_class(suppressWarnings(
    plotBoxplot(de, "Analyte_02", "Treatment", type = "bar",
                error_bar = "SE")), "ggplot")
  expect_error(plotBoxplot(data.frame(a = 1), "A", "B"), "DatasetExperiment")
})

test_that("free_y gives each facet its own y-axis, and only with facets", {
  de <- .de_plots()
  # Put batch B2 on a scale ten times larger, as tissue is against media.
  b2 <- de$sample_meta$Chrom_Batch == "B2"
  de$data[b2, "Analyte_02"] <- de$data[b2, "Analyte_02"] * 10

  y_range <- function(p){
    b <- ggplot2::ggplot_build(p)
    vapply(b$layout$panel_params, function(pp) pp$y.range[2], numeric(1))
  }

  shared <- suppressWarnings(plotBoxplot(de, "Analyte_02", "Treatment",
                                         facet_by = "Chrom_Batch"))
  free   <- suppressWarnings(plotBoxplot(de, "Analyte_02", "Treatment",
                                         facet_by = "Chrom_Batch",
                                         free_y = TRUE))
  # One cap across both panels, versus one per panel from its own values.
  expect_length(unique(y_range(shared)), 1)
  cap <- tapply(de$data$Analyte_02, de$sample_meta$Chrom_Batch, function(v)
    stats::quantile(v, 0.95, na.rm = TRUE) * 1.3)
  expect_gt(y_range(free)[2], 5 * y_range(free)[1])
  # Each panel tops out at its cap, plus ggplot's 5% padding - not above.
  expect_equal(unname(y_range(free)), as.vector(cap), tolerance = 0.06)
  # Both still start from zero.
  b <- ggplot2::ggplot_build(free)
  expect_true(all(vapply(b$layout$panel_params,
                         function(pp) pp$y.range[1] <= 0, logical(1))))

  # Without panels there is nothing to free: same plot as the default.
  expect_equal(
    y_range(suppressWarnings(plotBoxplot(de, "Analyte_02", "Treatment",
                                         free_y = TRUE))),
    y_range(suppressWarnings(plotBoxplot(de, "Analyte_02", "Treatment"))))

  expect_error(plotBoxplot(de, "Analyte_02", "Treatment", free_y = "yes"),
               "free_y must be TRUE or FALSE")
})

test_that("a point above the cap is drawn at it, and the box keeps its value", {
  de <- .de_plots()
  # One injection far above the rest of its group.
  i <- which(de$sample_meta$Treatment == "HDM")[1]
  de$data[i, "Analyte_02"] <- max(de$data$Analyte_02, na.rm = TRUE) * 50

  v     <- de$data$Analyte_02
  batch <- de$sample_meta$Chrom_Batch
  cap_of <- function(x) stats::quantile(x, 0.95, na.rm = TRUE) * 1.3

  for(free in c(FALSE, TRUE)){
    p <- suppressWarnings(plotBoxplot(de, "Analyte_02", "Treatment",
                                      facet_by = "Chrom_Batch", free_y = free))
    b <- suppressWarnings(ggplot2::ggplot_build(p))
    pts <- b$data[[which(vapply(p$layers, function(l)
      inherits(l$geom, "GeomPoint"), logical(1)))]]
    top <- vapply(b$layout$panel_params, function(pp) pp$y.range[2],
                  numeric(1))

    # A triangle for every value above its cap - the planted one among
    # them - and each drawn inside its panel rather than off the top.
    cap <- if(free) ave(v, batch, FUN = cap_of) else rep(cap_of(v), length(v))
    above <- !is.na(v) & v > cap
    expect_true(above[i])
    expect_equal(sum(pts$shape == 17), sum(above))
    tri <- pts$shape == 17
    expect_true(all(pts$y[tri] <= top[as.integer(pts$PANEL[tri])]))
    expect_match(p$labels$caption, "Triangles")
  }

  # The box is computed from the value as measured: its median matches the
  # data, the extreme point included, not the value it is drawn at.
  p  <- suppressWarnings(plotBoxplot(de, "Analyte_02", "Treatment"))
  bx <- suppressWarnings(ggplot2::ggplot_build(p))$data[[1]]
  hdm <- de$data$Analyte_02[de$sample_meta$Treatment == "HDM"]
  expect_equal(bx$middle[bx$x == which(sort(unique(
    de$sample_meta$Treatment)) == "HDM")], stats::median(hdm, na.rm = TRUE))

  # Nothing above the cap: round points only, and no caption.
  p <- suppressWarnings(plotBoxplot(.de_plots(), "Analyte_02", "Treatment",
                                    yaxis_quant = 1))
  expect_null(p$labels$caption)
})

test_that("plotQQ draws one compound, transformed or not", {
  de <- .de_plots()
  expect_s3_class(plotQQ(de, "Analyte_02"), "ggplot")
  expect_s3_class(plotQQ(de, "Analyte_02", transform = "log2"), "ggplot")
  expect_error(plotQQ(de, "Analyte_02", transform = "cuberoot"), "arg")
})

test_that("plotDrift needs a numeric injection order", {
  de <- .de_plots()
  expect_s3_class(suppressWarnings(plotDrift(de, "Analyte_02")), "ggplot")

  sm <- as.data.frame(de$sample_meta)
  sm$Injection_order <- paste0("run", seq_len(nrow(sm)))
  de$sample_meta <- sm
  expect_error(plotDrift(de, "Analyte_02"), "numeric")
})

test_that("the PCA plots take the fit and the dataset", {
  de  <- .de_plots()
  pca <- mrmRunPCA(de, transform = "log2")

  expect_s3_class(mrmPlotPCA(pca, de, colour_by = "Treatment"), "ggplot")
  expect_s3_class(
    mrmPlotPCA(pca, de, colour_by = "Treatment", ellipse = "group"), "ggplot")
  expect_s3_class(plotLoadings(pca, de), "ggplot")
  expect_s3_class(plotLoadings(pca, de, colour_by = "Enzymatic_pathway"),
                  "ggplot")

  # A fit that found nothing to decompose draws nothing, and says so by
  # returning NULL rather than failing.
  expect_null(plotLoadings(list(pr = NULL), de))
})

test_that("plotVolcano takes the stats table or the whole list", {
  st <- .stats_plots()
  expect_s3_class(plotVolcano(st, comparison = "PBS_vs_HDM"), "ggplot")
  expect_s3_class(plotVolcano(st$stats, comparison = "PBS_vs_HDM"), "ggplot")
})

test_that("plotCorrelation takes the correlations table, and says when it cannot", {
  st <- .stats_plots()
  p  <- plotCorrelation(st, name = "pairs")
  expect_s3_class(p, "ggplot")

  # The dataset is not what this one takes - the message has to say so.
  expect_error(plotCorrelation(.de_plots()), "correlations table")
  expect_error(plotCorrelation("pairs"), "correlations table")

  # Nothing to draw is NULL, not an error.
  expect_null(plotCorrelation(st$correlations[0, ]))
})
