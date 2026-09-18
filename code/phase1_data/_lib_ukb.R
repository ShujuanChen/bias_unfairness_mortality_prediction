# Clean and harmonise raw UKB records into PMR-aligned covariates.

suppressPackageStartupMessages({
  library(dplyr)
  library(stringr)
  library(lubridate)
})

process_ukb_sex <- function(df, col = "p31", out = "sex") {
  stopifnot(is.data.frame(df))
  if (!col %in% names(df)) stop(sprintf("Column '%s' not found.", col))

  raw <- .as_num(df[[col]])
  assert_expected_codes(raw, expected = c(0, 1), missing_codes = "", variable = col)
  mapped <- dplyr::case_when(raw == 1 ~ "1", raw == 0 ~ "2")
  df[[out]] <- factor(mapped, levels = .sex_levels)
  df
}

process_ukb_age <- function(df,
                                        year_col  = "p34",
                                        month_col = "p52",
                                        out       = "age",
                                        day       = 15L,
                                        census_date = COHORT_BASELINE) {
  stopifnot(is.data.frame(df))
  if (!year_col %in% names(df))  stop(sprintf("Column '%s' not found.", year_col))
  if (!month_col %in% names(df)) stop(sprintf("Column '%s' not found.", month_col))

  x <- df[[month_col]]
  x_trim <- str_trim(as.character(x))
  x_low  <- str_to_lower(x_trim)
  full   <- str_to_lower(month.name)
  mob <- suppressWarnings(as.integer(match(x_low, full)))

  yob <- suppressWarnings(as.integer(df[[year_col]]))
  dob <- ifelse(!is.na(yob) & !is.na(mob),
                as.Date(suppressWarnings(make_date(yob, mob, as.integer(day)))),
                as.Date(NA))

  age <- as.numeric(census_date - dob, units = "days") / 365.25
  age[!is.finite(age)] <- NA_real_

  df[[out]] <- floor(age)
  df
}

process_ukb_ethnicity <- function(df, col = "p21000_i0", out = "ethnicity") {
  stopifnot(is.data.frame(df))
  if (!col %in% names(df)) stop(sprintf("Column '%s' not found.", col))

  raw <- .as_num(df[[col]])

  assert_expected_codes(
    raw,
    expected = c(1:6, 1001:1003, 2001:2004, 3001:3004, 4001:4003),
    missing_codes = c(-1, -3, ""), variable = col)

  out_vec <- dplyr::case_when(
    raw %in% c(1, 1001:1003) ~ "1",
    raw %in% c(2, 2001:2004) ~ "2",
    raw %in% c(3, 3001:3004) ~ "3",
    raw %in% c(4, 4001:4003) ~ "4",
    raw == 5 ~ "5",
    raw == 6 ~ "6"
  )

  df[[out]] <- coerce_harmonised_ethnicity(out_vec, context = paste0(out, " mapped from ", col))
  df
}

process_ukb_tenure <- function(df, col = "p680_i0", out = "tenure") {
  stopifnot(is.data.frame(df))
  if (!col %in% names(df)) stop(sprintf("Column '%s' not found.", col))

  raw <- .as_num(df[[col]])

  assert_expected_codes(raw, expected = c(1:6, -7), missing_codes = c(-3, ""),
                        variable = col)

  mapped <- dplyr::case_when(
    raw == 1  ~ "1",
    raw == 2  ~ "2",
    raw == 3  ~ "4",
    raw == 4  ~ "5",
    raw == 5  ~ "3",
    raw == 6  ~ "6",
    raw == -7 ~ "7"
  )

  df[[out]] <- factor(mapped, levels = .tenure_levels)
  df
}

process_ukb_household_size <- function(df, col = "p709_i0", out = "household_size") {
  stopifnot(is.data.frame(df))
  if (!col %in% names(df)) stop(sprintf("Column '%s' not found.", col))

  num <- .as_num(df[[col]])
  assert_count_field(num, missing_codes = c(-1, -3, ""), variable = col)

  mapped <- dplyr::case_when(
    is.na(num)   ~ NA_character_,
    num == 1     ~ "1",
    num > 1      ~ "2",
    TRUE         ~ NA_character_
  )

  df[[out]] <- factor(mapped, levels = .household_size_levels)
  df
}

process_ukb_econstatus <- function(df,
                                   col_multi    = "p6142_i0",
                                   out          = "econstatus") {
  stopifnot(is.data.frame(df))
  if (!col_multi %in% names(df)) stop("Missing col: ", col_multi)

  col_multi_vec <- as.character(df[[col_multi]])
  col_multi_vec[is.na(col_multi_vec)] <- ""

  tokens_list <- str_split(col_multi_vec, "\\|", simplify = FALSE)
  tokens_list <- lapply(tokens_list, function(v) {
    v <- str_trim(v)
    v <- v[nzchar(v)]
    v
  })

  assert_expected_codes(unlist(tokens_list), expected = c(1:7, -7),
                        missing_codes = -3, variable = col_multi)

  is_missing <- (col_multi_vec == "") |
    vapply(tokens_list, function(tok) any(tok == "-3"), logical(1))

  has_emp <- vapply(tokens_list, function(tok) any(tok == "1"), logical(1))

  has_un  <- vapply(tokens_list, function(tok) any(tok == "5"), logical(1))

  has_ret <- vapply(tokens_list, function(tok) any(tok == "2"), logical(1))

  mapped <- character(length(col_multi_vec))
  mapped[] <- NA_character_

  mapped[!is_missing & has_emp] <- "1"
  mapped[!is_missing & !has_emp & has_un]  <- "2"
  mapped[!is_missing & !has_emp & !has_un & has_ret] <- "3"
  mapped[!is_missing & is.na(mapped)] <- "4"

  df[[out]] <- factor(mapped, levels = .econstatus_levels)
  df
}

process_ukb_education <- function(df,
                                  col = "p6138_i0",
                                  age_col = "p845_i0",
                                  out_level = "education") {
  stopifnot(is.data.frame(df))
  if (!col %in% names(df)) stop("Missing col: ", col)
  if (!age_col %in% names(df)) stop("Missing col: ", age_col)

  age_raw <- .as_num(df[[age_col]])
  assert_count_field(age_raw, missing_codes = c(-1, -2, -3, ""), variable = age_col)
  never_school <- !is.na(age_raw) & age_raw == -2

  age_clean <- age_raw
  age_clean[is.na(age_clean) | age_clean < 0] <- NA_real_

  raw <- as.character(df[[col]])
  raw[is.na(raw)] <- ""
  tokens_list <- str_split(raw, "\\|", simplify = FALSE)
  tokens_list <- lapply(tokens_list, function(v) {
    v <- str_trim(v)
    v <- v[nzchar(v)]
    v
  })

  TOK_DEGREE   <- "1"
  TOK_ALEVEL   <- "2"
  TOK_OLEVEL   <- "3"
  TOK_CSE      <- "4"
  TOK_NVQ_HND  <- "5"
  TOK_PROF     <- "6"
  TOK_NONE     <- "-7"
  TOK_PNA      <- "-3"

  assert_expected_codes(
    unlist(tokens_list),
    expected = c(TOK_DEGREE, TOK_ALEVEL, TOK_OLEVEL, TOK_CSE,
                 TOK_NVQ_HND, TOK_PROF, TOK_NONE),
    missing_codes = TOK_PNA, variable = col)

  years_list <- lapply(seq_along(tokens_list), function(i) {
    if (never_school[i]) return(0)

    toks <- tokens_list[[i]]
    if (length(toks) == 0) return(NA_real_)

    yrs <- numeric(0)

    if (any(toks == TOK_NONE))   yrs <- c(yrs, 7)
    if (any(toks == TOK_CSE))    yrs <- c(yrs, 10)
    if (any(toks == TOK_OLEVEL)) yrs <- c(yrs, 10)
    if (any(toks == TOK_ALEVEL)) yrs <- c(yrs, 13)
    if (any(toks == TOK_DEGREE)) yrs <- c(yrs, 20)

    if (any(toks == TOK_NVQ_HND)) {
      yrs <- c(yrs, ifelse(is.na(age_clean[i]), NA_real_, pmin(age_clean[i] - 5, 19)))
    }

    if (any(toks == TOK_PROF)) {
      yrs <- c(yrs, ifelse(is.na(age_clean[i]), NA_real_, pmin(age_clean[i] - 5, 15)))
    }

    if (any(toks == "" | toks == TOK_PNA)) yrs <- c(yrs, NA_real_)

    if (!length(yrs) || all(is.na(yrs))) return(NA_real_)
    max(yrs, na.rm = TRUE)
  })

  mapped_years <- unlist(years_list)
  mapped_years[is.infinite(mapped_years)] <- NA_real_

  mapped_level <- dplyr::case_when(
    is.na(mapped_years)     ~ NA_character_,
    mapped_years <= 8.5     ~ "1",
    mapped_years <= 11      ~ "2",
    mapped_years <= 17.5    ~ "3",
    mapped_years > 17.5     ~ "4",
    TRUE                    ~ NA_character_
  )

  df[[out_level]] <- factor(mapped_level, levels = .education_levels)
  df
}


process_ukb_ruralurban <- function(df, col = "p20118_i0", out = "ruralurban") {
  stopifnot(is.data.frame(df))
  if (!col %in% names(df)) stop(sprintf("Column '%s' not found.", col))

  raw <- .as_num(df[[col]])

  is_scotland <- !is.na(raw) & raw %in% 11:18

  assert_expected_codes(raw[!is_scotland], expected = 1:8,
                        missing_codes = c(9, ""), variable = col)

  mapped <- case_when(
    is_scotland ~ NA_character_,
    raw %in% c(1, 5) ~ "1",
    raw %in% c(2, 6) ~ "2",
    raw %in% c(3, 7) ~ "3",
    raw %in% c(4, 8) ~ "4",
    TRUE ~ NA_character_
  )

  df[[out]] <- factor(mapped, levels = .ruralurban_levels)
  df
}

process_ukb_health4 <- function(df, col = "p2178_i0", out = "health") {
  stopifnot(is.data.frame(df))
  if (!col %in% names(df)) stop(sprintf("Column '%s' not found.", col))

  raw <- .as_num(df[[col]])

  assert_expected_codes(raw, expected = 1:4,
                        missing_codes = c(-1, -3, ""), variable = col)
  mapped <- case_when(
    raw == 1 ~ "1",
    raw == 2 ~ "2",
    raw == 3 ~ "3",
    raw == 4 ~ "4",
    TRUE ~ NA_character_
  )

  df[[out]] <- factor(mapped, levels = .health_levels)
  df
}

process_ukb_lsoa11 <- function(df, col = "p20274_i0", out = "LSOA11CD") {
  stopifnot(is.data.frame(df))
  if (!col %in% names(df)) stop(sprintf("Column '%s' not found.", col))
  df[[col]][df[[col]] == ""] <- NA
  df[[out]] <- df[[col]]
  df
}

process_ukb_icd_code <- function(df, col, out) {
  stopifnot(is.data.frame(df))
  if (!col %in% names(df)) stop(sprintf("Column '%s' not found.", col))
  raw <- as.character(df[[col]])
  raw[is.na(raw) | trimws(raw) == ""] <- NA_character_
  code <- vapply(str_split(raw, "\\s+", n = 2L), function(x) {
    if (length(x) >= 1L && nzchar(trimws(x[1]))) trimws(x[1]) else NA_character_
  }, character(1))
  code[!nzchar(code)] <- NA_character_
  df[[out]] <- code
  df
}

.ukb_prepare_cols <- c(
  "eid", "p31", "p34", "p52", "p20274_i0",
  "p21000_i0", "p680_i0", "p709_i0", "p6142_i0",
  "p6138_i0", "p845_i0", "p20118_i0", "p2178_i0",
  "p40000_i0", "p40001_i0"
)

.ukb_prepare_out_cols <- c(
  "participant_id", "date_of_death",
  "underlying_cause_of_death_icd",
  "age", "sex", "ethnicity", "tenure", "household_size",
  "econstatus", "education", "ruralurban", "health",
  "LSOA11CD", "imd_decile"
)

ukb_prepare <- function(df, verbose = TRUE) {
  stopifnot(is.data.frame(df))

  needed <- intersect(.ukb_prepare_cols, names(df))
  missing <- setdiff(.ukb_prepare_cols, names(df))
  if (length(missing) > 0L) {
    stop("ukb_prepare: missing required columns: ", paste(missing, collapse = ", "))
  }
  out <- df %>% select(all_of(needed))

  out <- out %>%
    mutate(
      participant_id = eid,
      date_of_death = as.Date(p40000_i0)
    )

  if (isTRUE(verbose)) message("Processing ICD field: p40001_i0 -> underlying_cause_of_death_icd")
  out <- out %>% process_ukb_icd_code("p40001_i0", "underlying_cause_of_death_icd")

  out <- out %>%
    process_ukb_sex(col = "p31", out = "sex") %>%
    process_ukb_age(year_col = "p34",
                                month_col = "p52",
                                out = "age") %>%
    process_ukb_ethnicity(col = "p21000_i0", out = "ethnicity") %>%
    process_ukb_tenure(col = "p680_i0", out = "tenure") %>%
    process_ukb_household_size(col = "p709_i0", out = "household_size") %>%
    process_ukb_econstatus(col_multi = "p6142_i0", out = "econstatus") %>%
    process_ukb_education(col = "p6138_i0",
                          age_col = "p845_i0",
                          out_level = "education") %>%
    process_ukb_ruralurban(col = "p20118_i0", out = "ruralurban") %>%
    process_ukb_health4(col = "p2178_i0", out = "health") %>%
    process_ukb_lsoa11(col = "p20274_i0", out = "LSOA11CD") %>%
    process_imd_decile(lsoa11_col = "LSOA11CD", out = "imd_decile", verbose = verbose)

  out <- out %>% select(all_of(.ukb_prepare_out_cols))

  if (isTRUE(verbose)) {
    message("ukb_prepare: harmonisation summary")
    for (nm in names(out)) {
      .report(nm, out[[nm]], verbose = TRUE)
    }
  }

  tibble::as_tibble(out)
}
