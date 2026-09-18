#!/usr/bin/env Rscript
#
# Synthetic deprivation scores.

set.seed(7307)
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

if (!requireNamespace("writexl", quietly = TRUE)) {
  stop("Package 'writexl' is required to write the deprivation workbooks.",
       call. = FALSE)
}

eimd <- data.frame(.row = seq_along(lsoa01_e), check.names = FALSE,
                   stringsAsFactors = FALSE)
eimd[["LSOA CODE"]] <- lsoa01_e
eimd[["LSOA NAME"]] <- sprintf("SYNTH %04d", seq_along(lsoa01_e))
eimd[["IMD SCORE"]] <- round(pmin(pmax(rnorm(length(lsoa01_e), 22, 15), 0.5), 90), 3)
eimd$.row <- NULL

wimd <- data.frame(.row = seq_along(lsoa01_w), check.names = FALSE,
                   stringsAsFactors = FALSE)
wimd[["LSOA Code"]] <- lsoa01_w
wimd[["LSOA Name"]] <- sprintf("SYNTH %04d", seq_along(lsoa01_w))
wimd[["WIMD 2011 score"]] <- round(pmin(pmax(rnorm(length(lsoa01_w), 22, 15), 0.5), 90), 3)
wimd$.row <- NULL

write_two_sheet <- function(df, path, note) {
  writexl::write_xlsx(list(Notes = data.frame(note = note), Data = df), path)
  message("wrote ", path, "  (", nrow(df), " rows)")
}

write_two_sheet(eimd, file.path(out_dir, "EIMD2010.xlsx"),
                "Synthetic English deprivation scores. Not real.")
write_two_sheet(wimd, file.path(out_dir, "WIMD2011.xlsx"),
                "Synthetic Welsh deprivation scores. Not real.")
