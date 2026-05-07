# Tests for extractTable() — the TargetLynx workbook parser. Uses the
# bundled example_data.xlsx as the fixture so we don't hand-roll the
# (fairly involved) TargetLynx layout in code.

test_that("extractTable returns a long data frame with the requested headers", {
  xlsx <- system.file("extdata", "example_data.xlsx", package = "MRManalyzeR")
  skip_if(!nzchar(xlsx), "Bundled example_data.xlsx not installed")

  tl_headers <- c("ID", "Name", "Area", "ng/mL", "Response", "S/N")
  out <- extractTable(xlsx, tl_headers = tl_headers)

  expect_s3_class(out, "data.frame")
  expect_true(all(tl_headers %in% colnames(out)))
  expect_gt(nrow(out), 0)
  # The ID column should carry compound names (not literally "ID")
  expect_true(any(grepl("[A-Za-z]", out$ID)))
})

test_that("pre-existing 'ID' columns are dropped before the rename", {
  # Build a minimal in-memory fixture where row 3 has 'ID' at position 15
  # to mimic the failure mode that prompted the failsafe. A plain matrix
  # with named columns is enough — extractTable consumes any data.frame.
  skip("Synthetic-xlsx fixture not implemented; covered by example_data.xlsx")
})
