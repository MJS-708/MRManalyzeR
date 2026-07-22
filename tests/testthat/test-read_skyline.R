# Tests for the Skyline reader (.read_skyline_matrix, .apply_lod_loq_mask) —
# the Skyline-format counterpart to extractTable()/.apply_snr_mask() on the
# TargetLynx path. Uses small synthetic "Molecule x sample" xlsx fixtures
# rather than real Skyline exports.

.make_fdata <- function(processing_names, lod = NULL, loq = NULL){
  df <- data.frame(
    Processing_name = processing_names,
    Compound        = processing_names,
    Report          = "YES",
    Comment         = "",
    stringsAsFactors = FALSE
  )
  if(!is.null(lod)) df$LOD <- lod
  if(!is.null(loq)) df$LOQ <- loq
  df
}

# values must be a matrix with nrow = length(molecules), ncol = length(sample_names)
.write_skyline_sheet <- function(path, sheet_name, molecules, sample_names, values){
  df <- as.data.frame(values, stringsAsFactors = FALSE, check.names = FALSE)
  colnames(df) <- sample_names
  df <- cbind(data.frame(Molecule = molecules, stringsAsFactors = FALSE), df)

  wb <- if(file.exists(path)) openxlsx::loadWorkbook(path) else openxlsx::createWorkbook()
  openxlsx::addWorksheet(wb, sheet_name)
  openxlsx::writeData(wb, sheet_name, df)
  openxlsx::saveWorkbook(wb, path, overwrite = TRUE)
  path
}

test_that("reads a single sheet, transposes, and preserves non-syntactic names", {
  td   <- withr::local_tempdir()
  xlsx <- file.path(td, "skyline.xlsx")

  fdata <- .make_fdata(c("PGE2", "11(12)-EpETrE"))

  mol  <- c("PGE2", "11(12)-EpETrE", "[D4]-PGE2")   # last has no fdata match
  samp <- c("S1", "S2")
  vals <- matrix(c(10, 20,
                   30, 40,
                   50, 60), nrow = 3, byrow = TRUE)
  .write_skyline_sheet(xlsx, "skyline_data", mol, samp, vals)

  out <- .read_skyline_matrix(xlsx, "skyline_data", fdata)

  expect_equal(sort(rownames(out)), c("S1", "S2"))
  expect_true(all(c("PGE2", "11(12)-EpETrE") %in% colnames(out)))
  expect_true("[D4]-PGE2" %in% colnames(out))   # unmatched column passed through untouched
  expect_equal(out["S1", "PGE2"], 10)
  expect_equal(out["S2", "PGE2"], 20)
  expect_equal(out["S1", "11(12)-EpETrE"], 30)
  expect_equal(out["S2", "[D4]-PGE2"], 60)
})

test_that("case-insensitive match renames to the canonical Processing_name casing", {
  td   <- withr::local_tempdir()
  xlsx <- file.path(td, "skyline.xlsx")

  fdata <- .make_fdata("Mar2")
  .write_skyline_sheet(xlsx, "skyline_data",
                       molecules = "MaR2", sample_names = "S1",
                       values = matrix(5, nrow = 1))

  out <- .read_skyline_matrix(xlsx, "skyline_data", fdata)
  expect_true("Mar2" %in% colnames(out))      # canonical fdata casing, not raw "MaR2"
  expect_false("MaR2" %in% colnames(out))
})

test_that("sentinel non-detect tokens become NA", {
  td   <- withr::local_tempdir()
  xlsx <- file.path(td, "skyline.xlsx")
  fdata <- .make_fdata("PGE2")

  df <- data.frame(Molecule = "PGE2", S1 = "#N/A", S2 = "∞", S3 = "12.5",
                   stringsAsFactors = FALSE, check.names = FALSE)
  wb <- openxlsx::createWorkbook()
  openxlsx::addWorksheet(wb, "skyline_data")
  openxlsx::writeData(wb, "skyline_data", df)
  openxlsx::saveWorkbook(wb, xlsx, overwrite = TRUE)

  out <- .read_skyline_matrix(xlsx, "skyline_data", fdata)
  expect_true(is.na(out["S1", "PGE2"]))
  expect_true(is.na(out["S2", "PGE2"]))
  expect_equal(out["S3", "PGE2"], 12.5)
})

test_that(".apply_lod_loq_mask masks sub-threshold values per compound", {
  fdata <- .make_fdata(c("PGE2", "TXB2"), lod = c(5, 100), loq = c(1, 1))
  out_table <- data.frame(PGE2 = c(3, 10), TXB2 = c(50, 200),
                          row.names = c("S1", "S2"), check.names = FALSE)

  masked_lod <- .apply_lod_loq_mask(out_table, fdata, floor_col = "LOD")
  expect_true(is.na(masked_lod["S1", "PGE2"]))    # 3   < LOD 5
  expect_equal(masked_lod["S2", "PGE2"], 10)       # 10  >= LOD 5
  expect_true(is.na(masked_lod["S1", "TXB2"]))     # 50  < LOD 100
  expect_equal(masked_lod["S2", "TXB2"], 200)      # 200 >= LOD 100

  masked_loq <- .apply_lod_loq_mask(out_table, fdata, floor_col = "LOQ")
  expect_equal(masked_loq["S1", "PGE2"], 3)        # 3 >= LOQ 1, untouched
  expect_equal(masked_loq["S2", "PGE2"], 10)
})

test_that("multiple data_tab_names sheets are row-bound after validation", {
  td   <- withr::local_tempdir()
  xlsx <- file.path(td, "skyline.xlsx")
  fdata <- .make_fdata(c("PGE2", "TXB2"))

  mol <- c("PGE2", "TXB2")
  .write_skyline_sheet(xlsx, "batch1", mol, c("S1", "S2"), matrix(1:4,  nrow = 2, byrow = TRUE))
  .write_skyline_sheet(xlsx, "batch2", mol, c("S3", "S4"), matrix(5:8,  nrow = 2, byrow = TRUE))

  out <- .read_skyline_matrix(xlsx, c("batch1", "batch2"), fdata)
  expect_equal(sort(rownames(out)), c("S1", "S2", "S3", "S4"))
  expect_equal(nrow(out), 4)
})

test_that("duplicate sample names across sheets error", {
  td   <- withr::local_tempdir()
  xlsx <- file.path(td, "skyline.xlsx")
  fdata <- .make_fdata("PGE2")

  .write_skyline_sheet(xlsx, "batch1", "PGE2", c("S1", "S2"), matrix(1:2, nrow = 1))
  .write_skyline_sheet(xlsx, "batch2", "PGE2", c("S2", "S3"), matrix(3:4, nrow = 1))

  expect_error(.read_skyline_matrix(xlsx, c("batch1", "batch2"), fdata),
              "duplicate sample name")
})

test_that("mismatched Report=='YES' coverage across sheets errors", {
  td   <- withr::local_tempdir()
  xlsx <- file.path(td, "skyline.xlsx")
  fdata <- .make_fdata(c("PGE2", "TXB2"))

  .write_skyline_sheet(xlsx, "batch1", c("PGE2", "TXB2"), "S1", matrix(1:2, nrow = 2))
  .write_skyline_sheet(xlsx, "batch2", "PGE2",            "S2", matrix(3,   nrow = 1))

  expect_error(.read_skyline_matrix(xlsx, c("batch1", "batch2"), fdata),
              "disagree on which")
})
