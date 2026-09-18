#!/usr/bin/env Rscript
#
# Synthetic 2011 Census microdata

set.seed(5233)
N_IN <- 2789
N_OUT <- 5211
n <- N_IN + N_OUT

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

census <- data.frame(.row = seq_len(n), stringsAsFactors = FALSE)
p_country <- c("1" = 1)
census$country <- draw_cat(n, p_country, missing = 0, numeric = TRUE, integer = TRUE)
p_ageh <- c("9" = 0.195926, "10" = 0.195411, "11" = 0.17107, "13" = 0.16037,
            "12" = 0.150788, "14" = 0.126435)
census$ageh <- draw_cat(n, p_ageh, missing = 0, numeric = TRUE, integer = TRUE)
p_sex <- c("2" = 0.50703, "1" = 0.49297)
census$sex <- draw_cat(n, p_sex, missing = 0, numeric = TRUE, integer = TRUE)
p_carsnoc <- c("1" = 0.380488, "2" = 0.343331, "0" = 0.139431,
               "3" = 0.099511, "4" = 0.037238)
census$carsnoc <- draw_cat(n, p_carsnoc, missing = 0.006146, numeric = TRUE)
p_hlqupuk11 <- c("15" = 0.284136, "10" = 0.223561, "11" = 0.149373,
                 "12" = 0.141836, "14" = 0.098236, "16" = 0.057486,
                 "13" = 0.045372)
census$hlqupuk11 <- draw_cat(n, p_hlqupuk11, missing = 0.000116, numeric = TRUE)
p_health <- c("2" = 0.423797, "1" = 0.333627, "3" = 0.168762, "4" = 0.057387,
              "5" = 0.016427)
census$health <- draw_cat(n, p_health, missing = 0.000116, numeric = TRUE)
p_sizhuk11 <- c("2" = 0.383519, "3" = 0.19574, "4" = 0.170822,
                "1" = 0.157265, "5" = 0.060629, "6" = 0.021828,
                "7" = 0.009876, "0" = 0.000321)
census$sizhuk11 <- draw_cat(n, p_sizhuk11, missing = 0.006146, numeric = TRUE)
p_ethnicityew <- c("1" = 0.857196, "3" = 0.032374, "6" = 0.021589,
                   "12" = 0.015642, "2" = 0.012603, "11" = 0.012287,
                   "10" = 0.011793, "7" = 0.011348, "13" = 0.007182,
                   "9" = 0.005042, "5" = 0.004951, "4" = 0.004285,
                   "8" = 0.003707)
census$ethnicityew <- draw_cat(n, p_ethnicityew, missing = 0.000116, numeric = TRUE)
p_ecopuk11 <- c("2" = 0.383365, "10" = 0.178679, "1" = 0.155908,
                "6" = 0.062765, "13" = 0.058583, "12" = 0.038114,
                "5" = 0.032927, "7" = 0.032414, "4" = 0.025662,
                "14" = 0.019328, "3" = 0.004708, "11" = 0.003668,
                "8" = 0.003501, "9" = 0.000379)
census$ecopuk11 <- draw_cat(n, p_ecopuk11, missing = 0.000116, numeric = TRUE)
p_tenduk11 <- c("1" = 0.417191, "0" = 0.337057, "5" = 0.084725,
                "3" = 0.075648, "4" = 0.062879, "9" = 0.008956,
                "2" = 0.005397, "7" = 0.004441, "8" = 0.001896,
                "6" = 0.001811)
census$tenduk11 <- draw_cat(n, p_tenduk11, missing = 0.006462, numeric = TRUE)

p_excluded_band <- c("5" = 0.14255, "3" = 0.11224, "6" = 0.10846,
                     "8" = 0.10435, "7" = 0.10312, "1" = 0.09776,
                     "2" = 0.08825, "4" = 0.06129, "15" = 0.06054,
                     "16" = 0.04929, "17" = 0.0372, "18" = 0.023,
                     "19" = 0.01195)
excluded <- sample(n, N_OUT)
by_band <- excluded[runif(N_OUT) < 0.9163]
census$ageh[by_band] <- draw_cat(length(by_band), p_excluded_band,
                                 numeric = TRUE, integer = TRUE)
census$country[setdiff(excluded, by_band)] <- 2L

census$.row <- NULL
census <- census[sample(n), , drop = FALSE]
write_csv_file(census, file.path(out_dir, "recodev12.csv"))
