#!/usr/bin/env Rscript

# Generate UKB participation weights against the 2011 Census 5% microdata.

options(stringsAsFactors = FALSE)

script_match <- grep("--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- normalizePath(sub("--file=", "", script_match[1]), mustWork = TRUE)
script_dir <- dirname(script_path)

source(file.path(script_dir, "_lib_census.R"))

require_pipeline("census", "The Census participation model")

read_harmonised <- function(name, stage) {
  p <- path_temp("harmonised", name, create = FALSE)
  if (!file.exists(p)) {
    stop("Missing ", basename(p), ". Run phase 1 stage ", stage, " first.", call. = FALSE)
  }
  read_table(p, stringsAsFactors = FALSE)
}

predictors <- readLines(path_temp("harmonised", "census_auxiliaries.txt", create = FALSE))
ukb_prepared <- read_harmonised("ukb_with_census.csv", 4)
missing_aux <- setdiff(predictors, names(ukb_prepared))
if (length(missing_aux)) {
  stop("The harmonised UKB frame is missing declared auxiliaries: ",
       paste(missing_aux, collapse = ", "), call. = FALSE)
}
ukb_prepared <- ukb_prepared[, c("eid", "source", "sample_weight", predictors),
                             drop = FALSE]

keep <- analytic_ids()
before <- nrow(ukb_prepared)
ukb_prepared <- ukb_prepared[ukb_prepared$eid %in% keep, , drop = FALSE]
message(sprintf("[analytic population] UKB %d of %d rows retained", nrow(ukb_prepared), before))
if (nrow(ukb_prepared) != length(keep)) {
  stop(sprintf(paste("The analytic population holds %d participants but only %d are",
                     "present in the Census-harmonised frame. Rerun phase 1."),
               length(keep), nrow(ukb_prepared)), call. = FALSE)
}

census_prepared <- read_harmonised("census.csv", 4)
before_ref <- nrow(census_prepared)
census_prepared <- census_prepared[census_prepared$ref_id %in% reference_ids("census"), , drop = FALSE]
message(sprintf("[reference cohort] Census %d of %d rows retained",
                nrow(census_prepared), before_ref))

census_prepared$eid <- max(ukb_prepared$eid, na.rm = TRUE) + seq_len(nrow(census_prepared))
census_prepared$source <- "Census"
census_prepared$sample_weight <- as.numeric(
  framework_setting("phase3", "census_design_weight"))
if ("bootstrap_multiplier" %in% names(census_prepared)) {
  census_prepared$sample_weight <- census_prepared$sample_weight *
    as.numeric(census_prepared$bootstrap_multiplier)
  message("census: bootstrap multipliers applied to the reference weights")
}
census_prepared <- census_prepared[, c("eid", "source", "sample_weight", predictors), drop = FALSE]

stacked <- rbind(ukb_prepared, census_prepared)
model_input <- stacked[complete.cases(stacked[, predictors, drop = FALSE]), , drop = FALSE]

model_input_path <- census_model_input_path()
write_table(model_input, model_input_path)

model_input <- read_table(model_input_path, stringsAsFactors = FALSE)
fit <- fit_census_superlearner(model_input, predictors)

assert_covers_analytic_population(fit$ukb_weights$eid,
                                 label = "Census Super Learner weights")

weight_path <- census_weight_path()
write_table(fit$ukb_weights, weight_path)
write_pending_winsorisation(census_weights_dir())
write_table(fit$combined[, c("eid", "source", "sample_weight", "prob_ukb",
                             "crossfit_fold", "w")],
            census_model_output_path())

message("Saved Census+UKB SuperLearner weights to ", weight_path)
