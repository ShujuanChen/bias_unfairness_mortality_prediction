#!/usr/bin/env Rscript

# Harmonise the UKB outcome extract onto the shared coding scheme.
#
# Columns only. Every participant in the extract is written out, and who is in the
# study is decided in stage 5, against the rules in _lib_rows.R.

.dir <- dirname(normalizePath(sub("^--file=", "",
  grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1])))

source(file.path(.dir, "_lib_shared.R"))
source(file.path(.dir, "_lib_ukb.R"))
source(file.path(.dir, "_lib_geography.R"))
source(file.path(.dir, "_lib_rows.R"))
source(file.path(.dir, "..", "common", "config.R"))

cfg <- read_framework_config(.dir)

ukb_input_path <- input_path(path_data("UKB"), "UKB_extract.csv",
                             label = "UKB extract")
output_path <- path_temp("harmonised", "ukb_with_pmr.csv")

ukb <- read_input(ukb_input_path)
n_in <- nrow(ukb)

ukb <- ukb_prepare(ukb, verbose = TRUE)
ukb <- merge_geography(ukb, verbose = TRUE)

validate_harmonised_ethnicity(ukb$ethnicity, context = "UKB harmonised file")
assert_columns_only(n_in, ukb, "UKB harmonisation")

write_table(ukb, output_path)

message("Saved harmonised UKB dataset to ", output_path)
