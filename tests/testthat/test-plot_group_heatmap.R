# Tests for plot_group_heatmap() - the group-mean heatmap whose opacity
# carries each sample's contribution to the mean it is drawn from.
#
# The encoding is the thing worth testing: a group whose samples all move
# together must come out uniformly opaque, and a group whose mean rests on one
# sample must come out with one opaque cell and the rest at the floor. Those
# two cases are constructed explicitly rather than asserted about real data.

.de_for_heatmap <- function(){
  de <- load_dataset(system.file("extdata", "example_synthetic.RDS",
                                 package = "MRManalyzeR"))
  subset_dataset(de, conditions = list(Sample_type = "Sample"))
}

test_that("returns a ggplot and honours the blocking arguments", {
  de <- .de_for_heatmap()

  p <- plot_group_heatmap(de, group_by = "Treatment")
  expect_s3_class(p, "ggplot")

  p2 <- plot_group_heatmap(de, group_by = "Treatment",
                           group_features_by = "Enzymatic_pathway")
  expect_s3_class(p2, "ggplot")

  # one row per sample x feature, not one per group x feature: the whole point
  # is that individual samples are still drawn
  n_feat <- ncol(as.data.frame(de$data))
  expect_equal(nrow(p$data), nrow(as.data.frame(de$data)) * n_feat)
})

test_that("alpha reaches the floor and the ceiling for the right samples", {
  de <- .de_for_heatmap()
  p  <- plot_group_heatmap(de, group_by = "Treatment", alpha_floor = 0.25)

  expect_gte(min(p$data$alpha, na.rm = TRUE), 0.25)
  expect_lte(max(p$data$alpha, na.rm = TRUE), 1)

  # the floor must actually be reachable, or it is not doing anything
  expect_true(any(abs(p$data$alpha - 0.25) < 1e-8))
})

test_that("a sample opposing its group mean fades, one carrying it does not", {
  # Two groups, three compounds. In group B every sample moves up together; in
  # group A one sample is far up and the others sit flat, so A's mean is
  # carried by that one sample.
  set.seed(42)
  M <- matrix(0, nrow = 8, ncol = 3,
              dimnames = list(paste0("S", 1:8), paste0("F", 1:3)))
  M[1, ] <- 10            # A: the single sample carrying the mean
  M[2:4, ] <- 0           # A: contributing nothing
  M[5:8, ] <- 4           # B: all moving together
  M <- M + matrix(rnorm(24, sd = 0.01), nrow = 8)   # break exact ties

  sm <- data.frame(Treatment = rep(c("A", "B"), each = 4),
                   row.names = rownames(M))
  vm <- data.frame(Compound = colnames(M), row.names = colnames(M))
  de <- struct::DatasetExperiment(data = as.data.frame(M), sample_meta = sm,
                                  variable_meta = vm, name = "toy",
                                  description = "toy")

  p <- plot_group_heatmap(de, group_by = "Treatment", alpha_floor = 0.2)
  d <- p$data

  a_carrier <- mean(d$alpha[d$sample == "S1"])
  a_flat    <- mean(d$alpha[d$sample %in% c("S2", "S3", "S4")])
  b_all     <- mean(d$alpha[d$sample %in% c("S5", "S6", "S7", "S8")])

  expect_gt(a_carrier, a_flat)     # the sample that made the mean is opaque
  expect_gt(b_all, a_flat)         # a group moving together is opaque
})

test_that("degenerate input returns NULL rather than erroring", {
  M  <- matrix(1, nrow = 4, ncol = 3,
               dimnames = list(paste0("S", 1:4), paste0("F", 1:3)))
  sm <- data.frame(Treatment = rep(c("A", "B"), each = 2),
                   row.names = rownames(M))
  vm <- data.frame(Compound = colnames(M), row.names = colnames(M))
  de <- struct::DatasetExperiment(data = as.data.frame(M), sample_meta = sm,
                                  variable_meta = vm, name = "flat",
                                  description = "flat")

  # every compound is zero-variance, so nothing is left to scale
  expect_null(plot_group_heatmap(de, group_by = "Treatment"))
})

test_that("an unknown grouping column fails with a useful message", {
  de <- .de_for_heatmap()
  expect_error(plot_group_heatmap(de, group_by = "NoSuchColumn"),
               "not a sample_meta column")
})
