#!/usr/bin/env Rscript

# Harmonise the 2011 Census five per cent microdata and the matching UKB
# extract onto a shared auxiliary-variable scheme. Columns only.

.dir <- dirname(normalizePath(sub("^--file=", "",
  grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1])))

source(file.path(.dir, "..", "common", "config.R"))
source(file.path(.dir, "_lib_rows.R"))
source(file.path(.dir, "_lib_census.R"))
source(file.path(.dir, "..", "phase3_weighting", "_lib_census.R"))

cfg <- read_framework_config(.dir)

predictors <- census_predictors()

ukb_path <- input_path(path_data("UKB"), "UKB_extract.csv",
                       label = "UKB extract")
census_path <- input_path(path_data("census"), "recodev12.csv",
                          label = "Census five per cent microdata")

ukb <- prepare_ukb_for_census(ukb_path, predictors)
census <- prepare_census_microdata(census_path)

write_table(ukb, path_temp("harmonised", "ukb_with_census.csv"))
write_table(census, path_temp("harmonised", "census.csv"))
writeLines(predictors, path_temp("harmonised", "census_auxiliaries.txt"))

message(sprintf("UKB harmonised for Census: %d rows", nrow(ukb)))
message(sprintf("Census reference sample:   %d rows", nrow(census)))
message(sprintf("auxiliary variables:       %d", length(predictors)))
