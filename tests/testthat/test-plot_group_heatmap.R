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

# Two-group toy panels. The scaling is across ALL samples, so which sample
# "contributes" depends on where its group's mean sits relative to the global
# mean - not on which sample has the largest raw value. The next two tests pin
# down both sides of that, because it is the part of the encoding easiest to
# reason about wrongly.
.toy_de <- function(a_values, b_values, seed = 42L){
  set.seed(seed)
  M <- matrix(0, nrow = 8, ncol = 3,
              dimnames = list(paste0("S", 1:8), paste0("F", 1:3)))
  M[1:4, ] <- a_values
  M[5:8, ] <- b_values
  M <- M + matrix(rnorm(24, sd = 0.01), nrow = 8)   # break exact ties
  sm <- data.frame(Treatment = rep(c("A", "B"), each = 4),
                   row.names = rownames(M))
  vm <- data.frame(Compound = colnames(M), row.names = colnames(M))
  struct::DatasetExperiment(data = as.data.frame(M), sample_meta = sm,
                            variable_meta = vm, name = "toy",
                            description = "toy")
}

test_that("the sample that sets a group's mean is opaque, the rest are not", {
  # A = one sample high, three just above the global centre. S1's deviation is
  # 7x the others', so it outvotes them and A's mean comes out positive - S1 is
  # what produced it.
  de <- .toy_de(c(10, 1, 1, 1), rep(2, 4))
  d  <- plot_group_heatmap(de, group_by = "Treatment", alpha_floor = 0.2)$data

  a_carrier <- mean(d$alpha[d$sample == "S1"])
  a_rest    <- mean(d$alpha[d$sample %in% c("S2", "S3", "S4")])

  expect_gt(a_carrier, a_rest)
  expect_equal(a_rest, 0.2, tolerance = 1e-6)   # they pull the other way
  expect_gt(a_carrier, 0.9)                     # it alone carries the block
})

test_that("an extreme sample outvoted by its own group fades", {
  # The mirror case, and the one that caught me out writing these tests: A has
  # one sample far up and three far DOWN, so the three set the group mean and
  # the extreme sample is the one opposing it. Largest raw value does not mean
  # largest contribution.
  de <- .toy_de(c(10, 0, 0, 0), rep(4, 4))
  d  <- plot_group_heatmap(de, group_by = "Treatment", alpha_floor = 0.2)$data

  a_extreme  <- mean(d$alpha[d$sample == "S1"])
  a_majority <- mean(d$alpha[d$sample %in% c("S2", "S3", "S4")])

  expect_gt(a_majority, a_extreme)
  expect_equal(a_extreme, 0.2, tolerance = 1e-6)
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
