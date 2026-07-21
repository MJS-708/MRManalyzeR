# ---------------------------------------------------------------------------
# Build inst/extdata/example_batch_data.xlsx
#
# The primary bundled dataset (inst/extdata/example_data.xlsx) is a real
# TargetLynx export: a single chromatographic batch, one reconstitution volume
# for every vial, and no per-sample amount to normalise to. It therefore cannot
# demonstrate normalise_matrix(), adjust_concentration() or correct_batch() -
# each of those steps would be a no-op on it.
#
# This script derives a second, deliberately messier dataset from it:
#   * Skyline export layout (molecule x sample) rather than TargetLynx, so the
#     read_skyline() path is exercised too;
#   * two chromatographic batches with a simulated per-compound response shift
#     between them, so batch correction has something to remove;
#   * a per-sample protein amount (protein_ug) to normalise to;
#   * a different reconstitution volume per batch, plus a starting (pre-dry-down)
#     sample volume, so the concentration adjustment is not the identity;
#   * a per-compound LOD, so signal_filter = "LOD" removes something.
#
# The peak areas are real (Kolmert et al. 2018, the same subset as
# example_data.xlsx). The batch shift, protein amounts and volumes are
# SIMULATED, with a fixed seed, purely so the documentation can show each
# processing step doing visible work. Do not read biology into them.
#
# Run from the package root, after devtools::load_all():
#   source("data-raw/make_example_batch_data.R")
# ---------------------------------------------------------------------------

set.seed(708)

n_compounds <- 15L

src   <- system.file("extdata", "example_data.xlsx", package = "MRManalyzeR")
if(!nzchar(src)) src <- "inst/extdata/example_data.xlsx"
stopifnot(file.exists(src))

fdata <- openxlsx::read.xlsx(src, sheet = "feature_metadata")
meta  <- openxlsx::read.xlsx(src, sheet = "sample_metadata")

# Real peak areas, S/N-masked exactly as the TargetLynx path would read them.
areas <- read_targetlynx(src, datatype = "Area", snr = 3)

# --- Samples: drop the calibration ladder, keep Sample / QC / Blank ---------
keep_types <- c("Sample", "QC", "Blank")
meta <- meta[meta$Sample_type %in% keep_types & meta$Name %in% rownames(areas), ]
meta <- meta[order(as.numeric(meta$Injection_order)), ]
areas <- areas[meta$Name, , drop = FALSE]

# --- Compounds: the reportable ones with the most complete coverage ---------
reportable <- fdata$Processing_name[fdata$Report == "YES"]
cand       <- intersect(reportable, colnames(areas))
is_sample  <- meta$Sample_type == "Sample"
coverage   <- colSums(!is.na(areas[is_sample, cand, drop = FALSE]))
cand       <- sort(utils::head(names(sort(coverage, decreasing = TRUE)),
                               n_compounds))

areas <- areas[, cand, drop = FALSE]
fdata <- fdata[match(cand, fdata$Processing_name), ]

# Skyline users normally name their molecules with the reporting name, so this
# workbook has Processing_name == Compound (no rename step needed). The
# TargetLynx workbook keeps the two distinct - that is what Processing_name is
# for, and process_dataset() maps one to the other.
colnames(areas) <- fdata$Compound

# --- Simulate two chromatographic batches ----------------------------------
# Split the injection sequence in half; the second half is re-analysed on a
# different day and responds slightly differently, per compound.
n_inj  <- nrow(meta)
batch  <- ifelse(seq_len(n_inj) <= ceiling(n_inj / 2), "B1", "B2")
shift  <- exp(stats::rnorm(ncol(areas), mean = 0, sd = 0.30))
names(shift) <- colnames(areas)

for(cn in colnames(areas))
  areas[batch == "B2", cn] <- areas[batch == "B2", cn] * shift[[cn]]

# --- Simulate the per-sample metadata the steps consume --------------------
# protein_ug: the amount each extract was normalised to. Blanks and QCs get
# the pooled/nominal amount so they do not become all-NA rows downstream.
protein <- round(stats::rnorm(n_inj, mean = 250, sd = 45))
protein[meta$Sample_type != "Sample"] <- 250

sample_meta <- data.frame(
  Name             = meta$Name,
  Sample_ID        = meta$Sample_ID,
  Injection_order  = seq_len(n_inj),
  Treatment        = meta$Treatment,
  Type             = meta$Type,
  Sex              = meta$Sex,
  protein_ug       = protein,
  starting_vol_uL  = 200,        # plasma volume taken into extraction
  sample_volume_uL = ifelse(batch == "B1", 50, 100),  # reconstitution volume
  cal_vol_uL       = 50,
  sample_IS_vol_uL = 10,
  cal_IS_vol_uL    = 10,
  Chrom_Batch      = batch,
  Extraction_Batch = meta$Extraction_Batch,
  Include          = "YES",
  Sample_type      = meta$Sample_type,
  stringsAsFactors = FALSE,
  check.names      = FALSE
)

# --- Feature metadata, plus a per-compound LOD ------------------------------
# LOD set to the 5th percentile of the observed sample values, so a handful of
# low measurements fall below it and signal_filter = "LOD" has an effect.
lod <- apply(areas[is_sample, , drop = FALSE], 2,
             stats::quantile, probs = 0.05, na.rm = TRUE)

feature_meta <- data.frame(
  Processing_name  = fdata$Compound,
  Compound         = fdata$Compound,
  Report           = "YES",
  Comment          = NA_character_,
  LOD              = round(unname(lod[fdata$Compound])),
  Enzymatic_pathway = fdata$Enzymatic_pathway,
  PUFA             = fdata$PUFA,
  RT               = fdata$RT,
  stringsAsFactors = FALSE,
  check.names      = FALSE
)

# --- Skyline layout: molecules down the rows, samples across the columns ----
skyline_data <- data.frame(Molecule = colnames(areas),
                           check.names = FALSE, stringsAsFactors = FALSE)
skyline_data <- cbind(skyline_data,
                      as.data.frame(t(as.matrix(areas)), check.names = FALSE))
skyline_data[, -1] <- round(skyline_data[, -1], 1)

out <- "inst/extdata/example_batch_data.xlsx"
openxlsx::write.xlsx(
  list(feature_metadata = feature_meta,
       sample_metadata  = sample_meta,
       skyline_data     = skyline_data),
  file = out, overwrite = TRUE
)

message("[make_example_batch_data] wrote ", out, " - ",
        nrow(sample_meta), " injections x ", nrow(feature_meta), " compounds")
