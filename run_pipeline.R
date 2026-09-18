#!/usr/bin/env Rscript
# The single entry point for the analysis. 
#
#   Phase 1  Harmonisation     who is in the study, one denominator per source
#   Phase 2  Prediction bias   the gap between predicted and observed mortality
#   Phase 3  Weighting         weights toward the target population
#   Phase 4  Bias correction   how much of the gap closes
#   Phase 5  Uncertainty       intervals on the bias and on its correction
#
# Each invocation writes into temp/<run>/ and results/<run>/, the run being a
# timestamp unless --run-id names one to continue

suppressPackageStartupMessages(library(jsonlite))

args <- commandArgs(trailingOnly = TRUE)
script_path <- sub("^--file=", "",
                   grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE))
root <- if (length(script_path)) dirname(normalizePath(script_path[1])) else getwd()
source(file.path(root, "code", "common", "paths.R"))
source(file.path(root, "code", "common", "config.R"))
cfg <- read_framework_config(root)

opt <- function(flag, default = NULL) {
  i <- match(flag, args)
  if (is.na(i) || i == length(args)) return(default)
  nxt <- args[i + 1]
  if (startsWith(nxt, "--")) default else nxt
}
has <- function(flag) flag %in% args

USAGE <- paste(
  "Usage: Rscript run_pipeline.R [options]",
  "",
  "  --phase N                   one phase",
  "  --run-id ID                 continue the run under temp/ID",
  "  --stage N                   one stage, with --phase",
  "  --from N                    that phase onward",
  "  --all                       every phase, including phase 5",
  "  (no selection)              phases 1 to 4",
  "",
  "  --jobs N                    commands to run at once within a stage",
  "  --pipelines a,b,c           which of outcome, hse and census to build a",
  "                              population for. All three is the study. Naming",
  "                              fewer estimates weights against one reference",
  "                              sample without holding the other sources.",
  "",
  "Examples",
  "  Rscript run_pipeline.R --all",
  "  Rscript run_pipeline.R --phase 3 --stage 1 --run-id <run>",
  "  Rscript run_pipeline.R --phase 1 --pipelines hse",
  "  Rscript run_pipeline.R --phase 2 --jobs 3",
  sep = "\n")

KNOWN_FLAGS <- c("--all", "--phase", "--stage", "--from",
                 "--run-id", "--jobs", "--pipelines")

unknown_flags <- setdiff(grep("^--", args, value = TRUE), KNOWN_FLAGS)
if (length(unknown_flags)) {
  cat(USAGE, "\n\n", sep = "")
  stop("Unknown option: ", paste(unknown_flags, collapse = ", "), call. = FALSE)
}

RUN_ID <- opt("--run-id")
if (is.null(RUN_ID)) {
  RUN_ID <- format(Sys.time(), "%Y%m%d_%H%M%S")
} else if (!dir.exists(file.path(root, "temp", RUN_ID))) {
  stop("No run ", RUN_ID, " under temp/. Omit --run-id to start a new one.",
       call. = FALSE)
}
Sys.setenv(REWEIGHTING_RUN_ID = RUN_ID)

jobs <- as.integer(opt("--jobs", "1"))
if (is.na(jobs) || jobs < 1L) stop("--jobs must be a positive integer.", call. = FALSE)


PIPELINES <- opt("--pipelines", "outcome,hse,census")
PIPELINE_SET <- trimws(strsplit(PIPELINES, ",", fixed = TRUE)[[1]])
unknown <- setdiff(PIPELINE_SET, c("outcome", "hse", "census"))
if (length(unknown)) {
  stop("Unknown pipeline: ", paste(unknown, collapse = ", "),
       ". Choose from outcome, hse, census.", call. = FALSE)
}
wants <- function(p) p %in% PIPELINE_SET

outcomes <- names(cfg$phase2$outcomes)
folds <- seq_len(as.integer(cfg$phase2$splits$folds)) - 1L


WEIGHT_METHOD <- as.character(cfg$phase3$weight_method)
if (!WEIGHT_METHOD %in% names(cfg$phase2$weight_sources)) {
  stop("phase3.weight_method is '", WEIGHT_METHOD, "' and phase2.weight_sources ",
       "does not declare it.", call. = FALSE)
}
WEIGHTS <- paste0("ukbw_", WEIGHT_METHOD)

PYTHON <- Sys.getenv("REWEIGHTING_PYTHON", "python3")


R  <- function(...) c("Rscript", file.path(root, "code", ...))
PY <- function(...) c(PYTHON, file.path(root, "code", ...))


# ── Phase and stage declaration ──────────────────────────────────────────────
stage <- function(n, name, script, cmds, independent = FALSE) {
  list(n = n, name = name, script = script, cmds = cmds,
       independent = isTRUE(independent))
}

train <- function(script, cause, source) {
  if (identical(source, "pmr")) {
    lapply(folds, function(f) c(PY(script), "--cause", cause,
                                "--source", source, "--fold", f))
  } else {
    list(c(PY(script), "--cause", cause, "--source", source))
  }
}

P1 <- "phase1_data"; P2 <- "phase2_prediction_bias"; P3 <- "phase3_weighting"
P4 <- "phase4_correction"; P5 <- "phase5_uncertainty"

phases <- list(
  list(
    n = 1, name = "Harmonisation",
    stages = list(
      stage(1, "Harmonise UKB", file.path(P1, "s1_harmonise_ukb.R"),
            list(R(P1, "s1_harmonise_ukb.R"))),
      stage(2, "Harmonise PMR", file.path(P1, "s2_harmonise_pmr.R"),
            if (wants("outcome")) list(R(P1, "s2_harmonise_pmr.R")) else list()),
      stage(3, "Harmonise HSE", file.path(P1, "s3_harmonise_hse.R"),
            if (wants("hse")) list(R(P1, "s3_harmonise_hse.R")) else list()),
      stage(4, "Harmonise Census", file.path(P1, "s4_harmonise_census.R"),
            if (wants("census")) list(R(P1, "s4_harmonise_census.R")) else list()),
      stage(5, "Apply the row rules and fix the population", file.path(P1, "s5_build_population.R"),
            list(c(R(P1, "s5_build_population.R"), "--pipelines", PIPELINES))),
      stage(6, "Flow table", file.path(P1, "s6_flow_table.R"),
            if (length(PIPELINE_SET) == 3L) list(R(P1, "s6_flow_table.R")) else list())
    )
  ),

  list(
    n = 2, name = "Prediction bias",
    stages = list(
      stage(1, "Assign folds", file.path(P2, "s1_assign_folds.py"),
            lapply(outcomes, function(o) c(PY(P2, "s1_assign_folds.py"), "--cause", o)),
            independent = TRUE),
      stage(2, "Train", file.path(P2, "s2_train.py"),
            local({ out <- list(); for (o in outcomes)
                      out <- c(out, train(file.path(P2, "s2_train.py"), o, "ukb"),
                                    train(file.path(P2, "s2_train.py"), o, "pmr")); out }),
            independent = TRUE),
      stage(3, "Assemble predictions", file.path(P2, "s3_assemble.py"),
            lapply(outcomes, function(o) c(PY(P2, "s3_assemble.py"), "--cause", o,
                                           "--sources", "ukb,pmr")),
            independent = TRUE),
      stage(4, "Cumulative incidence", file.path(P2, "s4_cumulative_incidence.py"),
            list(c(PY(P2, "s4_cumulative_incidence.py"), "--sources", "ukb,pmr"))),
      stage(5, "Results", file.path(P2, "s5_results.py"),
            list(PY(P2, "s5_results.py")))
    )
  ),

  list(
    n = 3, name = "Weighting",
    stages = list(
      stage(1, "HSE participation weights", file.path(P3, "s1_estimate_hse.R"),
            if (wants("hse") && WEIGHT_METHOD %in% c("hse_superlearner", "hse_lassologit")) {
              list(c(R(P3, "s1_estimate_hse.R"),
                     sub("^hse_", "", WEIGHT_METHOD)))
            } else list()),
      stage(2, "Census participation weights", file.path(P3, "s2_estimate_census.R"),
            if (wants("census") && identical(WEIGHT_METHOD, "census_superlearner")) {
              list(c(R(P3, "s2_estimate_census.R"), "superlearner"))
            } else list()),
      stage(3, "Moment-matched weights", file.path(P3, "s3_moment_weights.R"),
            if (wants("hse") && WEIGHT_METHOD %in% c("hse_raking", "hse_entropy_balancing")) {
              list(c(R(P3, "s3_moment_weights.R"),
                     sub("^hse_", "", WEIGHT_METHOD)))
            } else list()),
      stage(4, "Weight diagnostics", file.path(P3, "s4_diagnostics.R"),
            list(R(P3, "s4_diagnostics.R"))),
      stage(5, "Results", file.path(P3, "s5_results.py"),
            list(PY(P3, "s5_results.py")))
    )
  ),

  list(
    n = 4, name = "Bias correction",
    stages = list(
      stage(1, "Train weighted", file.path(P2, "s2_train.py"),
            local({ out <- list(); for (o in outcomes) for (s in WEIGHTS)
                      out <- c(out, train(file.path(P2, "s2_train.py"), o, s)); out }),
            independent = TRUE),
      stage(2, "Assemble weighted", file.path(P2, "s3_assemble.py"),
            lapply(outcomes, function(o) c(PY(P2, "s3_assemble.py"), "--cause", o,
                                           "--sources", paste(c("ukb", "pmr", WEIGHTS), collapse = ","))),
            independent = TRUE),
      stage(3, "Cumulative incidence, weighted", file.path(P2, "s4_cumulative_incidence.py"),
            list(c(PY(P2, "s4_cumulative_incidence.py"), "--sources",
                   paste(c("ukb", "pmr", WEIGHTS), collapse = ",")))),
      stage(4, "Results", file.path(P4, "s1_results.py"),
            list(PY(P4, "s1_results.py")))
    )
  ),

  list(
    n = 5, name = "Uncertainty",
    stages = list(
      stage(1, "Replicates", file.path(P5, "s1_replicates.sh"),
            list(c("bash", file.path(root, "code", P5, "s1_replicates.sh"),
                   "1", as.character(cfg$phase5$replicates)))),
      stage(2, "Intervals", file.path(P5, "s2_intervals.py"),
            list(c(PY(P5, "s2_intervals.py"),
                   "--bootstrap-dir", "bootstrap_multiplier")))
    )
  )
)

# ── Selection ────────────────────────────────────────────────────────────────
sel_phase <- opt("--phase"); sel_stage <- opt("--stage"); from <- opt("--from")
chosen <- Filter(function(p) {
  if (!is.null(sel_phase)) return(p$n == as.integer(sel_phase))
  if (!is.null(from)) return(p$n >= as.integer(from))
  if (has("--all")) return(TRUE)
  p$n <= 4
}, phases)
if (length(chosen) == 0L) stop("No phase selected.", call. = FALSE)

if (!wants("outcome") && any(vapply(chosen, function(p) p$n %in% c(2, 4, 5), logical(1)))) {
  stop("Phases 2, 4 and 5 fit the outcome model, and this run's --pipelines ",
       "names ", PIPELINES, ". Name a phase this population supports, such as ",
       "--phase 1 or --phase 3, or drop --pipelines to build the whole study.",
       call. = FALSE)
}

DEVICE <- framework_setting("phase2", "device")
if (!DEVICE %in% c("gpu", "cpu")) {
  stop("phase2.device is ", DEVICE, ", and it takes \"gpu\" or \"cpu\".", call. = FALSE)
}

if (any(vapply(chosen, function(p) p$n %in% c(2, 4, 5), logical(1)))) {
  probe <- suppressWarnings(
    system2(PYTHON, c("-c", shQuote("import torch, pandas, numpy")),
            stdout = FALSE, stderr = FALSE))
  if (!identical(probe, 0L)) {
    stop("The Python stages need an interpreter with torch, pandas and numpy, ",
         "and '", PYTHON, "' does not have them. Point REWEIGHTING_PYTHON at ",
         "the interpreter you installed requirements.txt into, for example\n",
         "  export REWEIGHTING_PYTHON=/path/to/python\n",
         "Phases 1 and 3 are R only and run without it.", call. = FALSE)
  }
  if (identical(DEVICE, "gpu")) {
    on_gpu <- suppressWarnings(system2(
      PYTHON, c("-c", shQuote("import torch, sys; sys.exit(0 if torch.cuda.is_available() else 1)")),
      stdout = FALSE, stderr = FALSE))
    if (!identical(on_gpu, 0L)) {
      stop("phase2.device is \"gpu\" and '", PYTHON, "' finds no CUDA device. ",
           "Set phase2.device to \"cpu\" in framework_config.json.", call. = FALSE)
    }
  }
}
cat("device: ", DEVICE, "\n", sep = "")

manifest <- list(); t_start <- Sys.time()

manifest_path <- file.path(root, "temp", RUN_ID, "manifest",
                           paste0(format(t_start, "%Y%m%d_%H%M%S"), ".json"))

write_manifest <- function() {
  jsonlite::write_json(
    list(started = format(t_start, "%Y-%m-%dT%H:%M:%S%z"),
         finished = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
         minutes = round(as.numeric(difftime(Sys.time(), t_start, units = "mins")), 2),
         run_id = RUN_ID,
         r_version = R.version.string, platform = R.version$platform,
         python = PYTHON,
         weight_method = WEIGHT_METHOD,
         device = DEVICE,
         commands = manifest),
    manifest_path, auto_unbox = TRUE, pretty = TRUE)
}

dir.create(dirname(manifest_path), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(root, "results", RUN_ID), recursive = TRUE, showWarnings = FALSE)
cat("run:           ", RUN_ID, "\n",
    "intermediates: temp/", RUN_ID, "/\n",
    "results:       results/", RUN_ID, "/\n", sep = "")

run_one <- function(cmd) {
  t0 <- Sys.time()
  status <- system2(cmd[1], args = cmd[-1])
  list(command = paste(cmd, collapse = " "),
       exit_status = if (is.numeric(status)) as.integer(status) else 1L,
       seconds = round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1))
}

for (p in chosen) {
  stages_to_run <- if (!is.null(sel_stage)) {
    Filter(function(s) s$n == as.integer(sel_stage), p$stages)
  } else p$stages
  for (s in stages_to_run) {
    n_par <- if (s$independent) min(jobs, length(s$cmds)) else 1L
    cat(sprintf("\n=== %d.%d %s  (%d command%s%s) ===\n", p$n, s$n, s$name,
                length(s$cmds), if (length(s$cmds) == 1L) "" else "s",
                if (n_par > 1L) sprintf(", %d at a time", n_par) else ""))
    for (cmd in s$cmds) cat("  $ ", paste(cmd, collapse = " "), "\n", sep = "")
    if (!length(s$cmds)) next

    done <- if (n_par > 1L) {
      parallel::mclapply(s$cmds, run_one, mc.cores = n_par, mc.preschedule = FALSE)
    } else lapply(s$cmds, run_one)

    for (d in done) {
      if (!is.list(d) || is.null(d$exit_status)) {
        d <- list(command = "(unknown)", exit_status = 1L, seconds = NA)
      }
      manifest[[length(manifest) + 1]] <- c(
        list(phase = p$n, stage = s$n, stage_name = s$name), d)
    }
    failed <- Filter(function(d) !identical(d$exit_status, 0L),
                     manifest[seq(length(manifest) - length(done) + 1, length(manifest))])
    if (length(failed)) {
      write_manifest()
      stop(sprintf("Phase %d stage %d: %d of %d commands failed, first was: %s",
                   p$n, s$n, length(failed), length(done), failed[[1]]$command),
           call. = FALSE)
    }
  }
}

write_manifest()
cat(sprintf("\nCompleted in %.1f minutes.\n",
            as.numeric(difftime(Sys.time(), t_start, units = "mins"))))
