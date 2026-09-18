#!/usr/bin/env Rscript
#
# Fix the analytic population. Every row rule in the pipeline is applied here,

.dir <- dirname(normalizePath(sub("^--file=", "",
  grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1])))
source(file.path(.dir, "_lib_rows.R"))
source(file.path(.dir, "..", "common", "config.R"))

cfg <- read_framework_config(.dir)

fmt <- function(x) formatC(if (is.logical(x)) sum(x) else length(x),
                           big.mark = ",", format = "d")

read_stage <- function(name, stage, select = NULL) {
  p <- path_temp("harmonised", name, create = FALSE)
  if (!file.exists(p)) {
    stop("Missing ", basename(p), ". Run phase 1 stage ", stage, " first.", call. = FALSE)
  }
  read_table(p, select = select, data.table = FALSE,
                    na.strings = c("NA", ""))
}

read_vars <- function(name, df) {
  declared <- readLines(path_temp("harmonised", name, create = FALSE))
  absent <- setdiff(declared, names(df))
  if (length(absent)) {
    stop(name, " declares auxiliaries that are not in the harmonised file: ",
         paste(absent, collapse = ", "), call. = FALSE)
  }
  declared
}

require_vars <- function(df, vars, label) {
  absent <- setdiff(vars, names(df))
  if (length(absent)) {
    stop(label, " is missing required variables: ", paste(absent, collapse = ", "),
         call. = FALSE)
  }
  invisible(TRUE)
}

outcome_vars <- c(framework_setting("phase1", "predictors", "outcome", "numeric"),
                  framework_setting("phase1", "predictors", "outcome", "categorical"))

PIPELINES <- local({
  a <- commandArgs(trailingOnly = TRUE)
  i <- match("--pipelines", a)
  chosen <- if (is.na(i) || i == length(a)) "outcome,hse,census" else a[i + 1L]
  chosen <- trimws(strsplit(chosen, ",", fixed = TRUE)[[1]])
  chosen <- chosen[nzchar(chosen)]
  unknown <- setdiff(chosen, c("outcome", "hse", "census"))
  if (length(unknown)) {
    stop("Unknown pipeline: ", paste(unknown, collapse = ", "),
         ". Choose from outcome, hse, census.", call. = FALSE)
  }
  if (!length(chosen)) stop("--pipelines named none.", call. = FALSE)
  chosen
})
has_pipeline <- function(name) name %in% PIPELINES
PIPELINES_SORTED <- sort(PIPELINES)
message("pipelines: ", paste(PIPELINES, collapse = ", "))

# ── Read ─────────────────────────────────────────────────────────────────────

outcome <- read_stage("ukb_with_pmr.csv", 1,
                      select = union(c("participant_id", "LSOA11CD", "date_of_death",
                                       "RGN11CD"), outcome_vars))
pmr <- if (has_pipeline("outcome")) {
  read_stage("pmr_with_ukb.csv", 2,
             select = union(c("pmr_id", "w", "dod_deaths", "RGN11CD"), outcome_vars))
}
hse_ukb <- if (has_pipeline("hse")) read_stage("ukb_with_hse.csv", 3)
cen_ukb <- if (has_pipeline("census")) read_stage("ukb_with_census.csv", 4)
hse_ref <- if (has_pipeline("hse")) read_stage("hse.csv", 3)
cen_ref <- if (has_pipeline("census")) read_stage("census.csv", 4)

hse_vars <- if (has_pipeline("hse")) read_vars("hse_auxiliaries.txt", hse_ukb)
cen_vars <- if (has_pipeline("census")) read_vars("census_auxiliaries.txt", cen_ukb)
require_vars(outcome, outcome_vars, "The harmonised UKB outcome frame")
if (has_pipeline("outcome")) require_vars(pmr, outcome_vars, "The harmonised PMR frame")
if (has_pipeline("hse")) require_vars(hse_ref, hse_vars, "The harmonised HSE reference sample")
if (has_pipeline("census")) require_vars(cen_ref, cen_vars, "The harmonised Census reference sample")

ukb_extracts <- c("outcome",
                  if (has_pipeline("hse")) "HSE",
                  if (has_pipeline("census")) "Census")
for (nm in ukb_extracts) {
  ids <- switch(nm, outcome = outcome$participant_id, HSE = hse_ukb$eid, cen_ukb$eid)
  if (anyDuplicated(ids)) {
    stop("The UKB ", nm, " extract has duplicated participant identifiers.", call. = FALSE)
  }
}

for (nm in setdiff(ukb_extracts, "outcome")) {
  ids <- if (nm == "HSE") hse_ukb$eid else cen_ukb$eid
  if (!setequal(as.character(outcome$participant_id), as.character(ids))) {
    stop(sprintf(paste0("The UKB outcome extract and the UKB %s extract cover ",
                        "different participants: %s against %s, %s in the first ",
                        "only and %s in the second only. Harmonisation must not ",
                        "drop rows."),
                 nm, fmt(outcome$participant_id), fmt(ids),
                 fmt(setdiff(as.character(outcome$participant_id), as.character(ids))),
                 fmt(setdiff(as.character(ids), as.character(outcome$participant_id)))),
         call. = FALSE)
  }
}

# ── 1. Cohort rules ──────────────────────────────────────────────────────────

hse_weight <- if (has_pipeline("hse")) {
  hse_window <- Reduce(`&`, lapply(cohort_rules("HSE"), function(r) {
    pass <- r$test(hse_ref); pass[is.na(pass)] <- FALSE; pass
  }))
  hse_ref$weight_individual / mean(hse_ref$weight_individual[hse_window], na.rm = TRUE)
}

ukb_r <- apply_cohort_rules(outcome, "UKB", event = !is.na(outcome$date_of_death))
pmr_r <- if (has_pipeline("outcome")) {
  apply_cohort_rules(pmr, "PMR", weight = pmr$w, event = !is.na(pmr$dod_deaths))
}
hse_r <- if (has_pipeline("hse")) apply_cohort_rules(hse_ref, "HSE", weight = hse_weight)
cen_r <- if (has_pipeline("census")) apply_cohort_rules(cen_ref, "Census")

message("UKB cohort:    ", fmt(ukb_r$mask))
if (has_pipeline("outcome")) message("PMR cohort:    ", fmt(pmr_r$mask))
if (has_pipeline("hse")) message("HSE reference: ", fmt(hse_r$mask))
if (has_pipeline("census")) message("Census reference: ", fmt(cen_r$mask))

# ── 2. Completeness ──────────────────────────────────────────────────────────

complete_on <- function(df, vars) complete.cases(df[, vars, drop = FALSE])

ukb_complete_outcome <- complete_on(outcome, outcome_vars)
ukb_complete_hse <- if (has_pipeline("hse")) {
  complete_on(hse_ukb, hse_vars)[match(outcome$participant_id, hse_ukb$eid)]
} else rep(NA, nrow(outcome))
ukb_complete_census <- if (has_pipeline("census")) {
  complete_on(cen_ukb, cen_vars)[match(outcome$participant_id, cen_ukb$eid)]
} else rep(NA, nrow(outcome))
pmr_complete_outcome <- if (has_pipeline("outcome")) complete_on(pmr, outcome_vars)

hse_ref_complete <- if (has_pipeline("hse")) complete_on(hse_ref, hse_vars)
cen_ref_complete <- if (has_pipeline("census")) complete_on(cen_ref, cen_vars)

pmr_analytic <- if (has_pipeline("outcome")) pmr_r$mask & pmr_complete_outcome
hse_ref_analytic <- if (has_pipeline("hse")) hse_r$mask & hse_ref_complete
cen_ref_analytic <- if (has_pipeline("census")) cen_r$mask & cen_ref_complete

missingness_table <- function(df, mask, vars, source) {
  d <- df[mask, vars, drop = FALSE]
  out <- data.frame(
    source = source, variable = vars,
    n_missing = vapply(vars, function(v) sum(is.na(d[[v]])), integer(1)),
    row.names = NULL, stringsAsFactors = FALSE)
  out$pct <- round(100 * out$n_missing / nrow(d), 3)
  out[order(-out$n_missing), ]
}

ukb_mask_for <- function(ids) ukb_r$mask[match(ids, outcome$participant_id)]

write_table(missingness_table(outcome, ukb_r$mask, outcome_vars, "UKB"),
          path_temp("population", "ukb_missingness.csv"))
if (has_pipeline("outcome")) {
  write_table(missingness_table(pmr, pmr_r$mask, outcome_vars, "PMR"),
            path_temp("population", "pmr_missingness.csv"))
}
if (has_pipeline("hse")) {
  write_table(rbind(
    missingness_table(hse_ukb, ukb_mask_for(hse_ukb$eid), hse_vars, "UKB"),
    missingness_table(hse_ref, hse_r$mask, hse_vars, "HSE")),
    path_temp("population", "hse_missingness.csv"))
}
if (has_pipeline("census")) {
  write_table(rbind(
    missingness_table(cen_ukb, ukb_mask_for(cen_ukb$eid), cen_vars, "UKB"),
    missingness_table(cen_ref, cen_r$mask, cen_vars, "Census")),
    path_temp("population", "census_missingness.csv"))
}

message("outcome predictors complete:  ", fmt(ukb_r$mask & ukb_complete_outcome))
if (has_pipeline("hse")) message("HSE auxiliaries complete:     ", fmt(ukb_r$mask & ukb_complete_hse))
if (has_pipeline("census")) message("Census auxiliaries complete:  ", fmt(ukb_r$mask & ukb_complete_census))

# ── Common support ───────────────────────────────────────────────────────────

is_discrete <- function(x, max_levels = 25L) {
  !is.numeric(x) || length(unique(stats::na.omit(x))) <= max_levels
}

check_common_support <- function(df, keep, ref, ref_keep, vars, label) {
  d <- df[keep, , drop = FALSE]
  r <- ref[ref_keep, , drop = FALSE]
  failures <- character(0)

  for (v in vars) {
    if (!(v %in% names(d)) || !(v %in% names(r))) next
    if (!is_discrete(d[[v]]) || !is_discrete(r[[v]])) next

    ref_levels <- unique(stats::na.omit(as.character(r[[v]])))
    bad <- !is.na(d[[v]]) & !(as.character(d[[v]]) %in% ref_levels)
    if (any(bad)) {
      failures <- c(failures, sprintf(
        "%s / %s: %s cohort members hold level(s) the %s never carries: %s",
        label, v, formatC(sum(bad), big.mark = ",", format = "d"), label,
        paste(sort(unique(as.character(d[[v]][bad]))), collapse = ", ")))
    }
  }

  failures
}

reference_only_levels <- function(df, keep, ref, ref_keep, vars, label) {
  d <- df[keep, , drop = FALSE]
  r <- ref[ref_keep, , drop = FALSE]
  rows <- list()
  for (v in intersect(vars, intersect(names(d), names(r)))) {
    if (!is_discrete(d[[v]]) || !is_discrete(r[[v]])) next
    cohort_levels <- unique(stats::na.omit(as.character(d[[v]])))
    ref_vals <- as.character(r[[v]])
    for (lev in setdiff(unique(stats::na.omit(ref_vals)), cohort_levels)) {
      n <- sum(ref_vals == lev, na.rm = TRUE)
      rows[[length(rows) + 1L]] <- data.frame(
        reference = label, variable = v, level = lev, n_reference = n,
        pct_reference = round(100 * n / sum(!is.na(ref_vals)), 4),
        stringsAsFactors = FALSE)
    }
  }
  if (!length(rows)) {
    return(data.frame(reference = character(0), variable = character(0),
                      level = character(0), n_reference = integer(0),
                      pct_reference = numeric(0), stringsAsFactors = FALSE))
  }
  do.call(rbind, rows)
}

comparisons <- list(
  if (has_pipeline("outcome"))
    list(label = "target population", df = outcome, keep = ukb_r$mask,
         ref = pmr, ref_keep = pmr_analytic, vars = outcome_vars),
  if (has_pipeline("hse"))
    list(label = "HSE reference", df = hse_ukb, keep = ukb_mask_for(hse_ukb$eid),
         ref = hse_ref, ref_keep = hse_ref_analytic, vars = hse_vars),
  if (has_pipeline("census"))
    list(label = "Census reference", df = cen_ukb, keep = ukb_mask_for(cen_ukb$eid),
         ref = cen_ref, ref_keep = cen_ref_analytic, vars = cen_vars)
)
comparisons <- Filter(Negate(is.null), comparisons)

gaps <- do.call(rbind, lapply(comparisons, function(x)
  reference_only_levels(x$df, x$keep, x$ref, x$ref_keep, x$vars, x$label)))
write_table(gaps, path_temp("population", "reference_only_levels.csv"))
if (nrow(gaps)) {
  message("levels present in a reference sample but never in the cohort, ",
          "no weight can reach these:")
  for (i in seq_len(nrow(gaps))) {
    message(sprintf("  %s / %s = %s: %s reference records (%.3f%%)",
                    gaps$reference[i], gaps$variable[i], gaps$level[i],
                    formatC(gaps$n_reference[i], big.mark = ",", format = "d"),
                    gaps$pct_reference[i]))
  }
} else {
  message("every reference level is observed in the cohort")
}

failures <- unlist(lapply(comparisons, function(x)
  check_common_support(x$df, x$keep, x$ref, x$ref_keep, x$vars, x$label)))
if (length(failures)) {
  stop("Common support fails. Harmonisation must place every cohort member on a ",
       "level the reference sample carries.\n  ",
       paste(failures, collapse = "\n  "), call. = FALSE)
}
message("every cohort level is carried by the reference samples")

# ── Membership ───────────────────────────────────────────────────────────────

membership <- data.frame(
  participant_id = outcome$participant_id,
  ukb_r$flags,
  cohort           = ukb_r$mask,
  complete_outcome = ukb_complete_outcome,
  complete_hse     = ukb_complete_hse,
  complete_census  = ukb_complete_census,
  stringsAsFactors = FALSE
)

membership$in_unified <- membership$cohort
if (has_pipeline("outcome")) membership$in_unified <- membership$in_unified & membership$complete_outcome
if (has_pipeline("hse")) membership$in_unified <- membership$in_unified & membership$complete_hse
if (has_pipeline("census")) membership$in_unified <- membership$in_unified & membership$complete_census

out <- path_temp("population", "ukb_population_membership.csv")
write_table(membership, out)

if (has_pipeline("outcome")) {
  write_table(data.frame(pmr_id = pmr$pmr_id[pmr_analytic]),
            path_temp("population", "pmr_analytic_ids.csv"))
}
if (has_pipeline("hse")) {
  write_table(data.frame(ref_id = hse_ref$ref_id[hse_ref_analytic]),
            path_temp("population", "hse_reference_ids.csv"))
}
if (has_pipeline("census")) {
  write_table(data.frame(ref_id = cen_ref$ref_id[cen_ref_analytic]),
            path_temp("population", "census_reference_ids.csv"))
}

writeLines(PIPELINES_SORTED, path_temp("population", "pipelines.txt"))

# ── Attrition ────────────────────────────────────────────────────────────────

ukb_died <- !is.na(outcome$date_of_death)

ukb_steps <- list()
running <- membership$cohort
if (has_pipeline("outcome")) {
  running <- running & membership$complete_outcome
  ukb_steps <- c(ukb_steps, list(attrition_step(
    "UKB", "complete on outcome predictors", running, event = ukb_died)))
}
if (has_pipeline("hse")) {
  running <- running & membership$complete_hse
  ukb_steps <- c(ukb_steps, list(attrition_step(
    "UKB", "and complete on HSE auxiliaries", running, event = ukb_died)))
}
if (has_pipeline("census")) {
  running <- running & membership$complete_census
  ukb_steps <- c(ukb_steps, list(attrition_step(
    "UKB", "and complete on Census auxiliaries", running, event = ukb_died)))
}
ukb_steps <- c(ukb_steps, list(attrition_step(
  "UKB", "entering analysis", membership$in_unified, event = ukb_died)))

attrition <- rbind(
  ukb_r$attrition,
  do.call(rbind, ukb_steps),
  if (has_pipeline("outcome")) rbind(
    pmr_r$attrition,
    attrition_step("PMR", "complete on outcome predictors", pmr_analytic,
                   weight = pmr$w, event = !is.na(pmr$dod_deaths))),
  if (has_pipeline("hse")) rbind(
    hse_r$attrition,
    attrition_step("HSE", "complete on auxiliaries", hse_ref_analytic,
                   weight = hse_weight)),
  if (has_pipeline("census")) rbind(
    cen_r$attrition,
    attrition_step("Census", "complete on auxiliaries", cen_ref_analytic))
)
write_table(attrition, path_temp("population", "attrition.csv"))

message("cohort:                       ", fmt(ukb_r$mask))
message("unified population:           ", fmt(membership$in_unified))
message("analytic population:          ",
        formatC(length(analytic_ids()), big.mark = ",", format = "d"))
if (has_pipeline("outcome")) message("PMR analytic population:      ", fmt(pmr_analytic))
message("wrote ", out)
