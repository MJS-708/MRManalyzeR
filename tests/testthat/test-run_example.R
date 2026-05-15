# End-to-end smoke test for run_example(). Renders both reports and
# checks the four expected output artefacts exist. Slow (~30-60s) and
# pulls in plotly + rmarkdown — skipped on CRAN and skipped unless the
# env var MRMANALYZER_RUN_E2E is set, so the default `devtools::test()`
# stays fast. Set the env var to force a full smoke check:
#   Sys.setenv(MRMANALYZER_RUN_E2E = "true"); devtools::test()

test_that("run_example produces both reports + xlsx + RDS", {
  skip_on_cran()
  skip_if_not(nzchar(Sys.getenv("MRMANALYZER_RUN_E2E")),
              "Set MRMANALYZER_RUN_E2E=true to run the e2e smoke test")
  skip_if_not_installed("rmarkdown")
  skip_if_not_installed("plotly")

  td <- withr::local_tempdir()
  res <- run_example(results_dir = td, open = FALSE)

  out_dir <- attr(res, "results_dir")
  expect_equal(normalizePath(out_dir), normalizePath(td))

  files <- list.files(out_dir, recursive = TRUE)
  # Two HTML reports
  htmls <- grep("\\.html$", files, value = TRUE)
  expect_length(htmls, 2)
  # Processed-matrix RDS + xlsx; stats xlsx
  expect_true(any(grepl("\\.RDS$",   files)))
  expect_true(any(grepl("\\.xlsx$",  files)))
  expect_true(any(grepl("_stats\\.xlsx$", files)))
})
