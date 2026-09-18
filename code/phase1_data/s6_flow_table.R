#!/usr/bin/env Rscript

# Every number that appears in the four population flow figures, and nothing
# else, written to flow_table.xlsx.
#
#   ukb      records, UKB
#   pmr      represented population totals, PMR
#   hse      records, plus the survey-weighted total at the end
#   census   records, 2011 Census 5% microdata

.dir <- dirname(normalizePath(sub("^--file=", "",
  grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1])))

source(file.path(.dir, "_lib_rows.R"))
source(file.path(.dir, "..", "common", "config.R"))

cfg <- read_framework_config(.dir)

pop <- function(name) read_table(path_temp("population", name, create = FALSE),
                               stringsAsFactors = FALSE)

attrition <- pop("attrition.csv")

at <- function(source, step) {
  hit <- attrition[attrition$source == source & attrition$step == step, ]
  if (nrow(hit) != 1L) {
    stop(sprintf("attrition.csv holds %d rows for %s / %s, expected exactly one.",
                 nrow(hit), source, step), call. = FALSE)
  }
  hit
}

rows <- list()
add <- function(source, item, basis, n) {
  rows[[length(rows) + 1L]] <<- data.frame(
    source = source, item = item, basis = basis, n = as.numeric(n),
    stringsAsFactors = FALSE)
  invisible(NULL)
}

# ── UKB ───────────────────────────────────────────────────────────────

memb <- read_table(path_temp("population", "ukb_population_membership.csv",
                        create = FALSE), data.table = FALSE)

n_raw_ukb <- nrow(read_input(input_path(path_data("UKB"), "UKB_extract.csv",
                                        label = "UKB extract")))
if (n_raw_ukb != nrow(memb)) {
  stop(sprintf(paste0("The UKB extract holds %s records and the membership ",
                      "table holds %s. Harmonisation changed the number of rows."),
               format(n_raw_ukb, big.mark = ","),
               format(nrow(memb), big.mark = ",")), call. = FALSE)
}

u_england <- with(memb, rule_residence_known & rule_not_scotland &
                        rule_region_known & rule_england)
u_age     <- u_england & memb$rule_age
u_alive   <- u_age & memb$rule_alive

if (!identical(u_alive, memb$cohort)) {
  stop(paste0("Reordering the UKB cohort rules changed who is in the ",
              "cohort. The figure and the population disagree."), call. = FALSE)
}

add("ukb", "raw extract", "records", nrow(memb))
add("ukb", "living in England", "records", sum(u_england))
add("ukb", "excluded, not living in England", "records", nrow(memb) - sum(u_england))
add("ukb", cohort_age_label(" at 27 March 2011"), "records", sum(u_age))
add("ukb", "excluded, outside the age window", "records", sum(u_england) - sum(u_age))
add("ukb", "alive at 27 March 2011", "records", sum(u_alive))
add("ukb", "excluded, died before 27 March 2011", "records", sum(u_age) - sum(u_alive))
add("ukb", "entering the analysis", "records", sum(memb$in_unified))
add("ukb", "excluded, incomplete on at least one of the three sets", "records",
    sum(u_alive) - sum(memb$in_unified))

branch <- function(item, complete) {
  add("ukb", paste("complete on", item), "records", sum(u_alive & complete))
  add("ukb", paste("incomplete on", item), "records", sum(u_alive & !complete))
}
branch("the outcome predictors", memb$complete_outcome)
branch("the HSE auxiliaries",    memb$complete_hse)
branch("the Census auxiliaries", memb$complete_census)

# ── PMR ──────────────────────────────────────────────────────────────────────

pmr_missing <- pop("pmr_missingness.csv")
pmr_vars <- pmr_missing$variable

pmr <- read_table(path_temp("harmonised", "pmr_with_ukb.csv", create = FALSE),
             select = union(c("pmr_id", "w", "dod_deaths", "RGN11CD"), pmr_vars),
             data.table = FALSE, na.strings = c("NA", ""))

n_raw_pmr <- nrow(read_input(input_path(path_data("PMR"), "PMR.csv",
                                        label = "PMR linked census-mortality extract")))
if (n_raw_pmr != nrow(pmr)) {
  stop(sprintf(paste0("The PMR extract holds %s records and the harmonised file ",
                      "holds %s. Harmonisation changed the number of rows."),
               format(n_raw_pmr, big.mark = ","),
               format(nrow(pmr), big.mark = ",")), call. = FALSE)
}

false_na <- function(p) { p[is.na(p)] <- FALSE; p }

p_england <- false_na(!is.na(pmr$RGN11CD) &
                      startsWith(as.character(pmr$RGN11CD), "E"))
p_age     <- p_england & false_na(.in_age_window(pmr$age))
p_alive   <- p_age & false_na(.alive_at_baseline(pmr$dod_deaths))
p_complete <- p_alive &
  false_na(complete.cases(pmr[, pmr_vars, drop = FALSE]))

n_analytic <- at("PMR", "complete on outcome predictors")$n
if (sum(p_complete) != n_analytic) {
  stop(sprintf(paste0("The reordered PMR chain ends at %s and stage 5 recorded ",
                      "%s. The figure and the population disagree."),
               format(sum(p_complete), big.mark = ","),
               format(n_analytic, big.mark = ",")), call. = FALSE)
}

W <- function(mask) round(sum(pmr$w[mask], na.rm = TRUE))

add("pmr", "raw extract", "represented", W(rep(TRUE, nrow(pmr))))
add("pmr", "living in England", "represented", W(p_england))
add("pmr", "excluded, not living in England", "represented", W(!p_england))
add("pmr", cohort_age_label(" at 27 March 2011"), "represented", W(p_age))
add("pmr", "excluded, outside the age window", "represented", W(p_england & !p_age))
add("pmr", "alive at 27 March 2011", "represented", W(p_alive))
add("pmr", "excluded, died before 27 March 2011", "represented", W(p_age & !p_alive))
add("pmr", "complete on the outcome predictors", "represented", W(p_complete))
add("pmr", "incomplete on the outcome predictors", "represented", W(p_alive & !p_complete))

# ── HSE ────────────────────────────────────────────────

hse_raw <- at("HSE", "harmonised extract")$n
hse_age <- at("HSE", cohort_age_label())
hse_in  <- at("HSE", "complete on auxiliaries")

add("hse", "raw extract", "records", hse_raw)
add("hse", cohort_age_label(" at baseline"), "records", hse_age$n)
add("hse", "excluded, outside the age window", "records", hse_raw - hse_age$n)
add("hse", "entering the participation model", "records", hse_in$n)
add("hse", "entering the participation model", "survey weighted", hse_in$represented)
add("hse", "incomplete on the auxiliaries", "records", hse_age$n - hse_in$n)
add("hse", "auxiliary variables", "variables",
    length(readLines(path_temp("harmonised", "hse_auxiliaries.txt", create = FALSE))))

# ── 2011 Census ──────────────────────────────────────────────────────────────

cen_raw <- at("Census", "harmonised extract")$n
cen_eng <- at("Census", "England only")$n
cen_age <- at("Census", cohort_age_label())$n
cen_in  <- at("Census", "complete on auxiliaries")$n

add("census", "raw extract", "records", cen_raw)
add("census", "living in England", "records", cen_eng)
add("census", "excluded, not living in England", "records", cen_raw - cen_eng)
add("census", cohort_age_label(" at 27 March 2011"), "records", cen_age)
add("census", "excluded, outside the age window", "records", cen_eng - cen_age)
add("census", "entering the participation model", "records", cen_in)
add("census", "incomplete on the auxiliaries", "records", cen_age - cen_in)
add("census", "auxiliary variables", "variables",
    length(readLines(path_temp("harmonised", "census_auxiliaries.txt", create = FALSE))))

# ── Write ────────────────────────────────────────────────────────────────────

flow <- do.call(rbind, rows)
out <- path_results("phase1_harmonisation", "flow_table.xlsx")
write_table(flow, out)

print(flow, row.names = FALSE)
message("wrote ", out)
