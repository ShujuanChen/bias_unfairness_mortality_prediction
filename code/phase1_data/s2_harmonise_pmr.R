#!/usr/bin/env Rscript

# Harmonise the linked Census-mortality extract, PMR, onto the shared coding
# scheme.
#
# Columns only, as for the UKB extract. The cohort rules and the
# complete-case rule are applied in stage 5, where they are counted.

suppressPackageStartupMessages({
  library(dplyr)
})

.dir <- dirname(normalizePath(sub("^--file=", "",
  grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1])))

source(file.path(.dir, "_lib_shared.R"))
source(file.path(.dir, "_lib_pmr.R"))
source(file.path(.dir, "_lib_geography.R"))
source(file.path(.dir, "_lib_rows.R"))
source(file.path(.dir, "..", "common", "config.R"))

cfg <- read_framework_config(.dir)

pmr_input_path <- input_path(path_data("PMR"), "PMR.csv",
                             label = "PMR linked census-mortality extract")
output_path <- path_temp("harmonised", "pmr_with_ukb.csv")

pmr <- read_input(pmr_input_path)
n_in <- nrow(pmr)

pmr <- pmr_prepare(pmr, verbose = TRUE)
pmr <- merge_geography(pmr, verbose = TRUE)

pmr <- pmr %>% dplyr::mutate(pmr_id = dplyr::row_number())
pmr <- pmr %>% select(-any_of("LSOA11CD_join"))

validate_harmonised_ethnicity(pmr$ethnicity, context = "PMR harmonised file")
assert_columns_only(n_in, pmr, "PMR harmonisation")

write_table(pmr, output_path)

message("Saved harmonised PMR dataset to ", output_path)
