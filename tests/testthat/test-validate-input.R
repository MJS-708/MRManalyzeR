fm_ok <- function()
  data.frame(Compound        = c("A", "B"),
             Processing_name = c("a", "b"),
             Report          = c("YES", "YES"),
             stringsAsFactors = FALSE)

sm_ok <- function()
  data.frame(Name    = c("S1", "S2"),
             Include = c("YES", "YES"),
             stringsAsFactors = FALSE)

test_that("a well-formed pair of sheets passes", {
  v <- mrmValidateInput(fm_ok(), sm_ok())
  expect_length(v$errors, 0)
  expect_length(v$warnings, 0)
})

test_that("report_col is optional, but its absence is announced", {
  # A feature table curated outside a TargetLynx workflow - the usual case for
  # data_source = "matrix" - has no Report column. That must not be an error,
  # because there is no exclusion information to honour.
  fm <- fm_ok()[, c("Compound", "Processing_name")]

  v <- mrmValidateInput(fm, sm_ok())
  expect_length(v$errors, 0)
  expect_match(v$warnings, "no 'Report' column", all = FALSE)

  # Silence would be wrong too: a typo in report_col looks identical from here
  # and would quietly reinstate features that were meant to be dropped.
  expect_true(any(v$summary$status == "warning"))
})

test_that("duplicate compound names are still caught without report_col", {
  fm <- data.frame(Compound        = c("A", "A"),
                  Processing_name = c("a", "b"),
                  stringsAsFactors = FALSE)
  v <- mrmValidateInput(fm, sm_ok())
  expect_match(v$errors, "duplicate", all = FALSE)
})

test_that("compound and processing-name columns are still required", {
  fm <- data.frame(Something = c("A", "B"), stringsAsFactors = FALSE)
  v <- mrmValidateInput(fm, sm_ok())
  expect_match(v$errors, "missing required column", all = FALSE)
})

test_that("an included injection's order has to be a number", {
  sm <- sm_ok()
  sm$Injection_order <- c("_48", "2")
  v <- mrmValidateInput(fm_ok(), sm, injection_order_head = "Injection_order")
  expect_match(v$errors, "non-numeric Injection_order", all = FALSE)
  expect_match(v$errors, "'_48'", all = FALSE)

  sm$Injection_order <- c(1, 2)
  v <- mrmValidateInput(fm_ok(), sm, injection_order_head = "Injection_order")
  expect_length(v$errors, 0)

  # Excluded injections are not checked.
  sm$Injection_order <- c("x", "2")
  sm$Include <- c("NO", "YES")
  v <- mrmValidateInput(fm_ok(), sm, injection_order_head = "Injection_order")
  expect_length(v$errors, 0)
})

test_that("unreported rows with blank names are not duplicate compounds", {
  # A curated workbook often fills Report only on the rows it keeps and leaves
  # the rest blank, with no Compound name to go with them. Those rows are not
  # reported, so they cannot collide with anything; comparing the blank cell
  # with == gives NA, and an NA index used to pull NA names into the
  # uniqueness check and report them as duplicates.
  fm <- data.frame(
    Compound        = c("A", "B", NA, NA),
    Processing_name = c("a", "b", NA, NA),
    Report          = c("YES", "YES", NA, ""),
    stringsAsFactors = FALSE)

  v <- mrmValidateInput(fm, sm_ok())
  expect_length(v$errors, 0)
})

test_that("report_value is honoured rather than assumed to be YES", {
  fm <- data.frame(Compound        = c("A", "A", "B"),
                  Processing_name = c("a", "b", "c"),
                  Report          = c("keep", "keep", "YES"),
                  stringsAsFactors = FALSE)

  # The duplicate is among the rows this workbook marks as reported.
  v <- mrmValidateInput(fm, sm_ok(), report_value = "keep")
  expect_match(v$errors, "duplicate", all = FALSE)

  # Under the default the same duplicate is invisible, because neither row
  # counts as reported.
  expect_length(mrmValidateInput(fm, sm_ok())$errors, 0)
})
