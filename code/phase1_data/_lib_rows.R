# The study population: every row rule, for every source, in one file.

.rows_own_dir <- local({
  for (i in rev(seq_len(sys.nframe()))) {
    ofile <- sys.frame(i)$ofile
    if (!is.null(ofile)) return(dirname(normalizePath(ofile)))
  }
  stop("_lib_rows.R must be sourced, not run on its own.", call. = FALSE)
})
if (!exists("read_framework_config", mode = "function")) {
  source(file.path(.rows_own_dir, "..", "common", "config.R"))
}

COHORT_AGE_MIN  <- as.integer(framework_setting("phase1", "age_min"))
COHORT_AGE_MAX  <- as.integer(framework_setting("phase1", "age_max"))
COHORT_BASELINE <- as.Date(framework_setting("phase1", "baseline_date"))

cohort_age_label <- function(suffix = "") {
  paste0("aged ", COHORT_AGE_MIN, " to ", COHORT_AGE_MAX, suffix)
}

.in_age_window <- function(age) {
  age <- suppressWarnings(as.integer(age))
  !is.na(age) & age >= COHORT_AGE_MIN & age <= COHORT_AGE_MAX
}

.alive_at_baseline <- function(dod) {
  d <- suppressWarnings(as.Date(dod))
  is.na(d) | d >= COHORT_BASELINE
}

.residence_known <- function(lsoa) {
  x <- as.character(lsoa)
  !is.na(x) & nzchar(trimws(x))
}

.CENSUS_AGE_BANDS_PUBLISHED <- c("9" = 40L, "10" = 45L, "11" = 50L, "12" = 55L,
                                 "13" = 60L, "14" = 65L)

stopifnot(COHORT_AGE_MIN %in% .CENSUS_AGE_BANDS_PUBLISHED,
          COHORT_AGE_MAX %in% (.CENSUS_AGE_BANDS_PUBLISHED + 4L))

.CENSUS_AGE_BAND_LOWER <- .CENSUS_AGE_BANDS_PUBLISHED[
  .CENSUS_AGE_BANDS_PUBLISHED >= COHORT_AGE_MIN &
  .CENSUS_AGE_BANDS_PUBLISHED + 4L <= COHORT_AGE_MAX]

census_age_band_from_age <- function(age) {
  lower <- 5L * (suppressWarnings(as.integer(age)) %/% 5L)
  ok <- !is.na(lower) & lower %in% .CENSUS_AGE_BANDS_PUBLISHED
  ifelse(ok, paste0(lower, "-", lower + 4L), NA_character_)
}

census_age_band <- function(ageh) {
  lower <- .CENSUS_AGE_BANDS_PUBLISHED[as.character(ageh)]
  ifelse(is.na(lower), NA_character_, paste0(lower, "-", lower + 4L))
}

census_age_band_levels <- function() {
  lower <- sort(unname(.CENSUS_AGE_BAND_LOWER))
  paste0(lower, "-", lower + 4L)
}

.census_in_age_window <- function(ageh) {
  !is.na(.CENSUS_AGE_BAND_LOWER[as.character(ageh)])
}

.RULES <- list(
  UKB = list(
    list(id = "rule_age", label = cohort_age_label(" at 27 March 2011"),
         test = function(df) .in_age_window(df$age)),
    list(id = "rule_residence_known", label = "residence known",
         test = function(df) .residence_known(df$LSOA11CD)),
    list(id = "rule_not_scotland", label = "Scotland excluded",
         test = function(df) !.residence_known(df$LSOA11CD) |
                             substr(as.character(df$LSOA11CD), 1, 1) != "S"),
    list(id = "rule_alive", label = "alive at 27 March 2011",
         test = function(df) .alive_at_baseline(df$date_of_death)),
    list(id = "rule_region_known", label = "region resolved",
         test = function(df) !is.na(df$RGN11CD)),
    list(id = "rule_england", label = "England only",
         test = function(df) is.na(df$RGN11CD) | startsWith(as.character(df$RGN11CD), "E"))
  ),
  PMR = list(
    list(id = "rule_alive", label = "alive at 27 March 2011",
         test = function(df) .alive_at_baseline(df$dod_deaths)),
    list(id = "rule_age", label = cohort_age_label(" at 27 March 2011"),
         test = function(df) .in_age_window(df$age)),
    list(id = "rule_region_known", label = "region resolved",
         test = function(df) !is.na(df$RGN11CD)),
    list(id = "rule_england", label = "England only",
         test = function(df) is.na(df$RGN11CD) | startsWith(as.character(df$RGN11CD), "E"))
  ),
  HSE = list(
    list(id = "rule_age", label = cohort_age_label(),
         test = function(df) .in_age_window(df$age))
  ),
  Census = list(
    list(id = "rule_england", label = "England only",
         test = function(df) !is.na(df$country) & df$country == 1),
    list(id = "rule_age", label = cohort_age_label(),
         test = function(df) .census_in_age_window(df$ageh))
  )
)

cohort_rules <- function(source) {
  if (is.null(.RULES[[source]])) {
    stop("No cohort rules defined for source: ", source, call. = FALSE)
  }
  .RULES[[source]]
}

attrition_step <- function(source, label, mask, weight = NULL, event = NULL) {
  data.frame(
    source = source, step = label, n = sum(mask),
    represented = if (is.null(weight)) NA_real_ else round(sum(weight[mask], na.rm = TRUE)),
    deaths = if (is.null(event)) NA_integer_ else sum(event[mask], na.rm = TRUE),
    stringsAsFactors = FALSE)
}

apply_cohort_rules <- function(df, source, weight = NULL, event = NULL) {
  rules <- cohort_rules(source)
  n <- nrow(df)
  mask <- rep(TRUE, n)
  flags <- list()
  steps <- list(attrition_step(source, "harmonised extract", mask, weight, event))

  for (r in rules) {
    pass <- r$test(df)
    if (!is.logical(pass) || length(pass) != n) {
      stop(sprintf("Cohort rule %s/%s did not return a logical vector of length %d.",
                   source, r$id, n), call. = FALSE)
    }
    pass[is.na(pass)] <- FALSE
    flags[[r$id]] <- pass
    mask <- mask & pass
    steps[[length(steps) + 1L]] <- attrition_step(source, r$label, mask, weight, event)
  }

  list(mask = mask,
       flags = as.data.frame(flags, stringsAsFactors = FALSE),
       attrition = do.call(rbind, steps))
}

assert_columns_only <- function(n_in, df, label) {
  if (nrow(df) != n_in) {
    stop(sprintf(paste0("%s changed the number of rows, %s in and %s out. ",
                        "Harmonisation maps columns onto a shared scheme and must ",
                        "not drop or duplicate a row. Row rules belong in the ",
                        "population build, phase 1 stage 5, where they are counted."),
                 label, format(n_in, big.mark = ","), format(nrow(df), big.mark = ",")),
         call. = FALSE)
  }
  message(sprintf("[%s] %s rows in, %s out, columns only", label,
                  format(n_in, big.mark = ","), format(nrow(df), big.mark = ",")))
  invisible(TRUE)
}
