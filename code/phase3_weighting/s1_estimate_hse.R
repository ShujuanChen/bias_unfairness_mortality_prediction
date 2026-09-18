#!/usr/bin/env Rscript

# Generate UKB participation weights against the HSE.

options(stringsAsFactors = FALSE)

script_match <- grep("--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- normalizePath(sub("--file=", "", script_match[1]), mustWork = TRUE)
script_dir <- dirname(script_path)

source(file.path(script_dir, "_lib_hse.R"))

require_pipeline("hse", "The HSE participation model")

model <- commandArgs(trailingOnly = TRUE)
model <- if (length(model)) model[1] else "superlearner"
if (!model %in% c("superlearner", "lassologit")) {
  stop("Unknown model '", model, "'. Use 'superlearner' or 'lassologit'.")
}

sl_total_memory_gb <- function() {
  for (p in c("/sys/fs/cgroup/memory.max",
              "/sys/fs/cgroup/memory/memory.limit_in_bytes")) {
    if (file.exists(p)) {
      raw <- trimws(readLines(p, warn = FALSE)[1])
      if (!identical(raw, "max")) {
        bytes <- suppressWarnings(as.numeric(raw))
        if (!is.na(bytes) && bytes > 0 && bytes < 2^60) return(bytes / 1024^3)
      }
    }
  }
  if (file.exists("/proc/meminfo")) {
    kb <- as.numeric(gsub("[^0-9]", "",
                          grep("^MemTotal:", readLines("/proc/meminfo"),
                               value = TRUE)[1]))
    return(kb / 1024 / 1024)
  }
  bytes <- suppressWarnings(as.numeric(system2("sysctl", c("-n", "hw.memsize"),
                                               stdout = TRUE, stderr = FALSE)[1]))
  if (is.na(bytes)) NA_real_ else bytes / 1024^3
}

if (identical(model, "superlearner")) {
  need <- as.numeric(Sys.getenv("REWEIGHTING_SL_MIN_GB", "24"))
  have <- sl_total_memory_gb()
  if (!is.na(have) && have < need) {
    stop("This machine has ", format(have, digits = 3), " GB and the survey ",
         "Super Learner needs about ", need, ". A machine with less is killed ",
         "by the kernel hours into the first outer fold. Run it on a larger ",
         "machine, or set REWEIGHTING_SL_MIN_GB to override.", call. = FALSE)
  }
}

fit <- if (identical(model, "superlearner")) {
  fit_superlearner_weights()
} else {
  fit_lasso_weights()
}

model_input_path <- hse_model_input_path(model)
weight_path <- hse_weight_path(model)

write_table(fit$model_input, model_input_path)
assert_covers_analytic_population(fit$ukb_weights$eid,
                                 label = sprintf("HSE %s weights", model))
write_table(fit$ukb_weights, weight_path)
write_pending_winsorisation(hse_weights_dir(model))
write_table(fit$combined[, c("eid", "source", "sample_weight", "prob_ukb",
                             "crossfit_fold", "w")],
            hse_model_output_path(model))

message("Saved HSE ", model, " weights to ", weight_path)
