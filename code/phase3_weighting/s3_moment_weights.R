#!/usr/bin/env Rscript

# Raking or entropy balancing weights for UKB against HSE.

options(stringsAsFactors = FALSE)

script_match <- grep("--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- normalizePath(sub("--file=", "", script_match[1]), mustWork = TRUE)
script_dir <- dirname(script_path)

source(file.path(script_dir, "_lib_hse.R"))
source(file.path(script_dir, "_lib_moment.R"))

require_pipeline("hse", "Moment-matched participation weights")

method <- commandArgs(trailingOnly = TRUE)
method <- if (length(method)) method[1] else "raking"
if (!method %in% MOMENT_METHODS) {
  stop("Unknown method '", method, "'. Use one of: ",
       paste(MOMENT_METHODS, collapse = ", "), call. = FALSE)
}

message("[hse_", method, "] preparing model input")
prepared <- prepare_training_data()
prep <- list(combined = prepared$combined, info = prepared$prediction_info)
is_ukb <- prep$combined$source == "UKB"
message(sprintf("[hse_%s] %d cohort rows, %d reference rows, %d predictors",
                method, sum(is_ukb), sum(!is_ukb), nrow(prep$info)))

fit <- if (identical(method, "raking")) rake_weights(prep) else entropy_balance_weights(prep)

if (!isTRUE(fit$converged)) {
  message(sprintf(paste0(
    "\n[WARNING] %s did not converge: %d iterations, largest remaining gap ",
    "%.3g against a tolerance of %.3g. The cohort has %d rows and the ",
    "reference sample %d, on %d predictors, and the two are too far apart on ",
    "some margin for any weighting to close. The weights are written and the ",
    "run continues, but they do not satisfy the constraints that define them ",
    "and no result from them should be reported.\n"),
    method, fit$iterations, fit$achieved_gap,
    if (identical(method, "raking")) RAKE_TOL else ENTROPY_TOL,
    sum(is_ukb), sum(!is_ukb), nrow(prep$info)))
}

w <- fit$w
balance <- moment_balance_table(prep, w)
ess <- effective_sample_size(w)

ukb_weights <- data.frame(eid = prep$combined$eid[is_ukb], w = w)
assert_covers_analytic_population(ukb_weights$eid,
                                  label = sprintf("HSE %s weights", method))

write_table(ukb_weights, moment_weight_path(method))
write_table(prepared$model_input, moment_model_input_path(method))

out <- data.frame(eid = prep$combined$eid,
                  source = prep$combined$source,
                  sample_weight = prep$combined$sample_weight,
                  prob_ukb = NA_real_,
                  crossfit_fold = NA_integer_,
                  w = NA_real_)
out$w[is_ukb] <- w
out$w[!is_ukb] <- prep$combined$sample_weight[!is_ukb]
write_table(out, moment_model_output_path(method))

is_con <- balance$level %in% c("mean", "sd")
worst_cat <- balance[!is_con, ][which.max(abs(balance$gap_after[!is_con])), ]

message(sprintf("[hse_%s] ESS %.2f%%, max weight %.2f, largest categorical gap %.3g pp on %s %s",
                method, 100 * ess / length(w), max(w),
                worst_cat$gap_after, worst_cat$variable, worst_cat$level))
message("Saved hse_", method, " weights to ", moment_weight_path(method))
