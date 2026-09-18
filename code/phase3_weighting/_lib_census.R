options(stringsAsFactors = FALSE)

# Core of the Census participation-weighting pipeline.

census_predictors <- function() {
  framework_setting("phase1", "predictors", "census_participation")
}

# ---- result-tree paths -------------------------------------------------------

census_weights_dir <- function() {
  d <- path_temp("weights", "census_superlearner", create = FALSE)
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
  d
}

census_model_input_path <- function() file.path(census_weights_dir(), "model_input.csv")
census_weight_path <- function() file.path(census_weights_dir(), "ukb_weights.csv")
census_model_output_path <- function() file.path(census_weights_dir(), "model_output.csv")

# ---- harmonisation -----------------------------------------------------------

.resolve_own_dir <- function() {
  for (i in rev(seq_len(sys.nframe()))) {
    ofile <- sys.frame(i)$ofile
    if (!is.null(ofile)) return(dirname(normalizePath(ofile)))
  }
  arg <- grep("--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(arg)) return(dirname(normalizePath(sub("--file=", "", arg[1]))))
  getwd()
}

.census_code_dir <- tryCatch(.resolve_own_dir(), error = function(e) getwd())

source(file.path(.census_code_dir, "..", "common", "config.R"))
source(file.path(.census_code_dir, "..", "phase1_data", "_lib_census.R"))
source(file.path(.census_code_dir, "_lib_superlearner.R"))
source(file.path(.census_code_dir, "_lib_weights.R"))

# ---- Super Learner fit -------------------------------------------------------

coerce_census_predictors <- function(df, predictors) {
  for (v in predictors) df[[v]] <- as.factor(df[[v]])
  df
}

census_active_predictors <- function(df, predictors) {
  keep <- vapply(predictors, function(v) {
    length(unique(df[[v]][!is.na(df[[v]])])) > 1
  }, logical(1))
  if (!all(keep)) {
    stop("The Census auxiliaries ", paste(predictors[!keep], collapse = ", "),
         " take one value on the stacked cohort and reference sample, so the ",
         "participation model cannot be fitted on the ", length(predictors),
         " that framework_config.json declares.", call. = FALSE)
  }
  predictors
}

fit_census_superlearner <- function(model_input, predictors) {
  prepared <- coerce_census_predictors(model_input, predictors)
  active <- census_active_predictors(prepared, predictors)

  x_mm <- participation_matrix(prepared, active)

  y <- as.integer(prepared$source == "UKB")
  sample_weight <- prepared$sample_weight

  x <- participation_design(as.data.frame(x_mm), character(0))

  crossfit <- crossfit_superlearner(y = y, x = x, obs_weights = sample_weight)
  prob_ukb <- crossfit$prob_ukb

  ukb_idx <- prepared$source == "UKB"
  if (any(!is.finite(prob_ukb) | prob_ukb <= 0 | prob_ukb >= 1, na.rm = TRUE)) {
    stop("SuperLearner participation probabilities produced invalid inverse-odds weights.")
  }

  inverse_odds <- (1 - prob_ukb) / prob_ukb
  prepared$prob_ukb <- prob_ukb
  prepared$inverse_odds <- inverse_odds
  prepared$crossfit_fold <- crossfit_fold_index(crossfit$fit, nrow(prepared))
  prepared$w <- NA_real_
  prepared$w[ukb_idx] <- normalise_weights(inverse_odds[ukb_idx])
  prepared$w[!ukb_idx] <- prepared$sample_weight[!ukb_idx]

  list(
    combined = prepared,
    ukb_weights = prepared[ukb_idx, c("eid", "w", "prob_ukb"), drop = FALSE]
  )
}
