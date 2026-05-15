# Tests for run_pca_pipeline() — the PCA preprocessing pipeline used by
# both reports. Pure function on a numeric matrix, easy to unit test.

# Helper: deterministic synthetic matrix with planted structure.
.synth_matrix <- function(n_samples = 20, n_features = 12, seed = 1L){
  set.seed(seed)
  M <- matrix(rlnorm(n_samples * n_features, meanlog = 4, sdlog = 0.6),
              nrow = n_samples, ncol = n_features)
  rownames(M) <- paste0("S", seq_len(n_samples))
  colnames(M) <- paste0("F", seq_len(n_features))
  M
}

test_that("returns expected list structure", {
  M <- .synth_matrix()
  res <- run_pca_pipeline(M)

  expect_named(res, c("pr", "X_clean", "n_samples", "n_features",
                      "dropped_all_na", "dropped_feat_filter",
                      "dropped_sample_filter", "dropped_zero_var", "steps"))
  expect_s3_class(res$pr, "prcomp")
  expect_equal(res$n_samples, 20)
  expect_equal(res$n_features, 12)
  expect_type(res$steps, "character")
  expect_gt(length(res$steps), 0)
})

test_that("drops all-NA features regardless of feat_na_max", {
  M <- .synth_matrix()
  M[, 3] <- NA_real_
  res <- run_pca_pipeline(M, feat_na_max = 1)
  expect_equal(res$dropped_all_na, 1)
  expect_equal(res$n_features, 11)
})

test_that("feat_na_max drops features above threshold", {
  M <- .synth_matrix()
  # Make column 5 have 70% NAs
  M[1:14, 5] <- NA_real_
  res <- run_pca_pipeline(M, feat_na_max = 0.5, impute = "min")
  expect_equal(res$dropped_feat_filter, 1)
  expect_false("F5" %in% colnames(res$X_clean))
})

test_that("sample_na_max drops samples above threshold", {
  M <- .synth_matrix()
  # Make row 1 have 90% NAs
  M[1, 1:11] <- NA_real_
  res <- run_pca_pipeline(M, sample_na_max = 0.5, impute = "min")
  expect_equal(res$dropped_sample_filter, 1)
  expect_false("S1" %in% rownames(res$X_clean))
})

test_that("frac_min imputation uses per-feature minimum × frac", {
  M <- .synth_matrix()
  feat_min <- min(M[, 1])
  M[5, 1] <- NA_real_
  res <- run_pca_pipeline(M, impute = "frac_min", impute_frac = 0.2,
                          transform = "none", center = FALSE, scale = FALSE)
  # Imputed value should equal min × 0.2
  expect_equal(res$X_clean["S5", "F1"], feat_min * 0.2)
})

test_that("transform options apply correctly", {
  M <- matrix(c(1, 2, 4, 8), nrow = 2, ncol = 2,
              dimnames = list(c("S1", "S2"), c("F1", "F2")))

  res_none <- run_pca_pipeline(M, transform = "none", impute = "none",
                               center = FALSE, scale = FALSE)
  res_log2 <- run_pca_pipeline(M, transform = "log2", impute = "none",
                               center = FALSE, scale = FALSE)
  res_sqrt <- run_pca_pipeline(M, transform = "sqrt", impute = "none",
                               center = FALSE, scale = FALSE)

  expect_equal(unname(res_none$X_clean), unname(M))
  expect_equal(unname(res_log2$X_clean), unname(log2(M + 1)))
  expect_equal(unname(res_sqrt$X_clean), unname(sqrt(M)))
})

test_that("insufficient data yields NULL pr without error", {
  M <- matrix(1, nrow = 2, ncol = 2,
              dimnames = list(c("S1", "S2"), c("F1", "F2")))
  res <- run_pca_pipeline(M)
  expect_null(res$pr)                        # too few samples
  expect_true(any(grepl("Insufficient", res$steps)))
})

test_that("zero-variance features are dropped after transform", {
  M <- .synth_matrix()
  M[, 4] <- 5                                # constant column
  res <- run_pca_pipeline(M, transform = "none", center = FALSE, scale = FALSE)
  expect_gte(res$dropped_zero_var, 1)
  expect_false("F4" %in% colnames(res$X_clean))
})
