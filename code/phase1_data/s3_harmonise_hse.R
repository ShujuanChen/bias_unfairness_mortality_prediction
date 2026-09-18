#!/usr/bin/env Rscript

# Harmonise the HSE waves and the matching UKB
# extract onto a shared auxiliary-variable scheme. Columns only.

.dir <- dirname(normalizePath(sub("^--file=", "",
  grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1])))

source(file.path(.dir, "..", "common", "config.R"))
source(file.path(.dir, "_lib_shared.R"))
source(file.path(.dir, "_lib_hse.R"))

cfg <- read_framework_config(.dir)

predictors <- hse_predictors()

ukb_path <- input_path(path_data("UKB"), "UKB_extract.csv",
                       label = "UKB extract")

ukb <- prepare_ukb_for_hse(ukb_path)
hse <- prepare_hse_waves(path_data("HSE"))

absent <- setdiff(predictors, names(ukb))
if (length(absent)) {
  stop("Declared auxiliaries are not in the harmonised UKB frame: ",
       paste(absent, collapse = ", "), call. = FALSE)
}

write_table(ukb, path_temp("harmonised", "ukb_with_hse.csv"))
write_table(hse, path_temp("harmonised", "hse.csv"))
writeLines(predictors, path_temp("harmonised", "hse_auxiliaries.txt"))

message(sprintf("UKB harmonised for HSE: %d rows", nrow(ukb)))
message(sprintf("HSE reference sample:   %d rows", nrow(hse)))
message(sprintf("auxiliary variables:    %d", length(predictors)))
