#!/usr/bin/env Rscript
#
# Synthetic geography lookups.

set.seed(8431)
n <- 0

out_dir <- local({
  a <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", a[grep("^--file=", a)])
  if (length(f)) dirname(normalizePath(f)) else getwd()
})

N_LSOA_E <- 1800
N_LSOA_W <- 200
lsoa11_e <- sprintf("E01%06d", seq_len(N_LSOA_E))
lsoa11_w <- sprintf("W01%06d", seq_len(N_LSOA_W))
lsoa11 <- c(lsoa11_e, lsoa11_w)
lsoa01_e <- sprintf("E01%06d", seq_len(N_LSOA_E))
lsoa01_w <- sprintf("W01%06d", seq_len(N_LSOA_W))
lsoa01 <- c(lsoa01_e, lsoa01_w)

draw_cat <- function(n, p, missing = 0, numeric = FALSE, integer = FALSE) {
  x <- sample(names(p), n, replace = TRUE, prob = as.numeric(p))
  if (missing > 0) x[runif(n) < missing] <- NA_character_
  if (integer) as.integer(x) else if (numeric) as.numeric(x) else x
}
draw_num <- function(n, mean, sd, lo, hi, missing = 0, digits = NULL) {
  x <- pmin(pmax(rnorm(n, mean, sd), lo), hi)
  if (!is.null(digits)) x <- round(x, digits)
  if (missing > 0) x[runif(n) < missing] <- NA_real_
  x
}
draw_date <- function(n, from, to, missing = 0) {
  a <- as.Date(from); b <- as.Date(to)
  x <- format(a + floor(runif(n) * as.numeric(b - a)), "%Y-%m-%d")
  if (missing > 0) x[runif(n) < missing] <- NA_character_
  x
}
write_csv_file <- function(df, path) {
  utils::write.csv(df, path, row.names = FALSE, na = "NA")
  message("wrote ", path, "  (", nrow(df), " rows, ", ncol(df), " columns)")
}

region_codes <- sprintf("E1200000%02d", 1:9)
region_of <- function(x) region_codes[1 + (seq_along(x) %% length(region_codes))]

lsoa01_to_11 <- data.frame(
  LSOA01CD = lsoa01, LSOA01NM = sprintf("SYNTH %05d", seq_along(lsoa01)),
  LSOA11CD = lsoa11, LSOA11NM = sprintf("SYNTH %05d", seq_along(lsoa11)),
  CHGIND = "U",
  LAD11CD = sprintf("E0600%04d", 1 + (seq_along(lsoa11) %% 300)),
  LAD11NM = sprintf("SYNTH LAD %04d", 1 + (seq_along(lsoa11) %% 300)),
  LAD11NMW = NA_character_, FID = seq_along(lsoa11),
  stringsAsFactors = FALSE)

lsoa11_to_region <- data.frame(
  LSOA11CD = lsoa11, LSOA11NM = sprintf("SYNTH %05d", seq_along(lsoa11)),
  BUASD11CD = NA_character_, BUASD11NM = NA_character_,
  BUA11CD = NA_character_, BUA11NM = NA_character_,
  LAD11CD = lsoa01_to_11$LAD11CD, LAD11NM = lsoa01_to_11$LAD11NM,
  LAD11NMW = NA_character_,
  RGN11CD = region_of(lsoa11), RGN11NM = NA_character_, RGN11NMW = NA_character_,
  ObjectId = seq_along(lsoa11), stringsAsFactors = FALSE)

write_csv_file(lsoa01_to_11, file.path(out_dir, "LSOA01_LSOA11_LAD11_Lookup_EW.csv"))
write_csv_file(lsoa11_to_region, file.path(out_dir, paste0(
  "Lower_Layer_Super_Output_Area_(2011)_to_Built-up_Area_Sub-division_to_",
  "Built-up_Area_to_Local_Authority_District_to_Region_(December_2011)_",
  "Lookup_in_England_and_Wales.csv")))
