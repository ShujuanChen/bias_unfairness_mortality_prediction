# Clean and harmonise raw PMR to UKB-aligned covariates. Defines pmr_prepare().

suppressPackageStartupMessages({
  library(dplyr)
  library(stringr)
})

process_pmr_icd_missing <- function(df, cols, verbose = TRUE) {
  stopifnot(is.data.frame(df))
  for (col in cols) {
    if (!col %in% names(df)) stop(sprintf("Column '%s' not found.", col))
    x <- as.character(df[[col]])
    x[is.na(x)] <- NA_character_
    x_trim <- trimws(x)
    x_trim[x_trim == "" | toupper(x_trim) == "NA"] <- NA_character_
    df[[col]] <- x_trim
    .report(col, df[[col]], verbose)
  }
  df
}

pmr_sampling_weight <- function(df, verbose = TRUE) {
  x <- .as_num(df$sampling_weight)
  .report("sampling_weight", x, verbose)
  x
}

pmr_age <- function(df, verbose = TRUE) {
  x <- .as_num(df$age_census)
  .report("age", x, verbose)
  x
}

pmr_dod_deaths <- function(df, verbose = TRUE) {
  x <- suppressWarnings(as.Date(df$dod_deaths))
  .report("dod_deaths", x, verbose)
  x
}

pmr_sex <- function(df, verbose = TRUE) {
  code <- .as_num(df$sex_census)
  assert_expected_codes(code, expected = c(1, 2), variable = "sex_census")
  lab <- dplyr::case_when(
    code == 1 ~ "1",
    code == 2 ~ "2",
    TRUE ~ NA_character_
  )
  f <- factor(lab, levels = .sex_levels)
  .report("sex", f, verbose)
  f
}

process_pmr_ethnicity <- function(df, col = "ethpuk11_census", out = "ethnicity") {
  stopifnot(is.data.frame(df))
  if (!col %in% names(df)) stop(sprintf("Column '%s' not found.", col))

  x <- toupper(trimws(as.character(df[[col]])))
  x[is.na(x)] <- ""

  assert_expected_codes(
    x, expected = sprintf("%02d", 1:18),
    missing_codes = c("", "NA", "XX"), variable = col)
  out_vec <- case_when(
    x %in% c("01", "02", "03", "04") ~ "1",
    x %in% c("05", "06", "07", "08") ~ "2",
    x %in% c("09", "10", "11", "13") ~ "3",
    x %in% c("14", "15", "16") ~ "4",
    x %in% c("12") ~ "5",
    x %in% c("17", "18") ~ "6",
    TRUE ~ NA_character_
  )

  df[[out]] <- coerce_harmonised_ethnicity(out_vec, context = paste0(out, " mapped from ", col))
  if (col != out) df <- df %>% select(-any_of(col))
  df
}

process_pmr_tenure <- function(df, col = "tenhuk11_census", out = "tenure", verbose = TRUE) {
  stopifnot(is.data.frame(df))
  if (!col %in% names(df)) stop(sprintf("Column '%s' not found.", col))

  assert_expected_codes(toupper(trimws(as.character(df[[col]]))),
                        expected = 0:9, missing_codes = c("X", "NA", ""),
                        variable = col)
  code <- .as_num(df[[col]])
  lab <- dplyr::case_when(
    code == 0 ~ "1",
    code == 1 ~ "2",
    code == 2 ~ "3",
    code %in% 3:4 ~ "4",
    code == 5 ~ "5",
    code %in% 6:8 ~ "7",
    code == 9 ~ "6",
    TRUE ~ NA_character_
  )
  f <- factor(lab, levels = .tenure_levels)
  .report(out, f, verbose)
  df[[out]] <- f
  if (col != out) df <- df %>% select(-any_of(col))
  df
}

process_pmr_econstatus <- function(df, col = "ecocatpuk11_census", out = "econstatus", verbose = TRUE) {
  stopifnot(is.data.frame(df))
  if (!col %in% names(df)) stop("Missing col: ", col)

  chr <- as.character(df[[col]])
  chr[is.na(chr)] <- NA_character_
  chr[trimws(chr) == "NA"] <- NA_character_

  assert_expected_codes(chr, expected = 1:8, missing_codes = c("X", ""),
                        variable = col)

  mapped <- dplyr::case_when(
    chr %in% c("1", "2", "3", "4") ~ "1",
    chr %in% c("5", "6")           ~ "2",
    chr == "7"                     ~ "3",
    chr == "8"                     ~ "4",
    TRUE                          ~ NA_character_
  )

  f <- factor(mapped, levels = .econstatus_levels)
  .report(out, f, verbose)
  df[[out]] <- f
  if (col != out) df <- df %>% select(-any_of(col))
  df
}

process_pmr_education <- function(df, col = "hlqpuk11_census", out_level = "education", verbose = TRUE) {
  stopifnot(is.data.frame(df))
  if (!col %in% names(df)) stop("Missing col: ", col)

  chr <- toupper(trimws(as.character(df[[col]])))
  chr[chr == "NA"] <- NA_character_

  assert_expected_codes(chr, expected = 10:16, missing_codes = c("XX", ""),
                        variable = col)

  mapped_years <- dplyr::case_when(
    chr == "10"              ~ 7,
    chr %in% c("11", "12")   ~ 10,
    chr == "13"              ~ 12,
    chr == "14"              ~ 13,
    chr == "15"              ~ 20,
    chr == "16"              ~ 15
  )

  mapped_level <- dplyr::case_when(
    is.na(mapped_years)     ~ NA_character_,
    mapped_years <= 8.5     ~ "1",
    mapped_years <= 11       ~ "2",
    mapped_years <= 17.5     ~ "3",
    mapped_years > 17.5      ~ "4",
    TRUE                     ~ NA_character_
  )

  f <- factor(mapped_level, levels = .education_levels)
  .report(out_level, f, verbose)
  df[[out_level]] <- f
  if (col != out_level) df <- df %>% select(-any_of(col))
  df
}

process_pmr_ruralurban <- function(df, col = "ruralurban_code_census", out = "ruralurban", verbose = TRUE) {
  stopifnot(is.data.frame(df))
  if (!col %in% names(df)) stop(sprintf("Column '%s' not found.", col))

  chr <- as.character(df[[col]])
  chr[is.na(chr)] <- NA_character_
  chr_trim <- trimws(chr)
  chr_trim[chr_trim == ""] <- NA_character_
  chr_up <- toupper(chr_trim)
  chr_up[chr_up == "NA"] <- NA_character_

  assert_expected_codes(
    chr_up, expected = c("A1", "A2", "C1", "C2", "D1", "D2", "E1", "E2", "F1", "F2"),
    missing_codes = "", variable = col)

  mapped <- case_when(
    chr_up %in% c("A1", "A2", "C1", "C2") ~ "1",
    chr_up %in% c("D1", "D2")             ~ "2",
    chr_up %in% c("E1", "E2")             ~ "3",
    chr_up %in% c("F1", "F2")             ~ "4",
    is.na(chr_up)                         ~ NA_character_,
    TRUE                                  ~ NA_character_
  )

  f <- factor(mapped, levels = .ruralurban_levels)
  .report(out, f, verbose)
  df[[out]] <- f
  if (col != out) df <- df %>% select(-any_of(col))
  df
}

process_pmr_household_size <- function(df, col = "hhchuk11_census", out = "household_size", verbose = TRUE) {
  stopifnot(is.data.frame(df))
  if (!col %in% names(df)) stop(sprintf("Column '%s' not found.", col))

  raw <- df[[col]]
  chr <- as.character(raw)
  chr[is.na(raw)] <- NA_character_
  chr_trim <- str_trim(chr)
  chr_trim[chr_trim == ""] <- NA_character_
  chr_trim[toupper(chr_trim) == "NA"] <- NA_character_

  norm <- ifelse(nchar(chr_trim) == 1L & chr_trim %in% as.character(1:9), paste0("0", chr_trim), chr_trim)
  norm[is.na(chr_trim)] <- NA_character_
  norm[toupper(chr_trim) == "XX" & !is.na(chr_trim)] <- "XX"

  assert_expected_codes(norm, expected = sprintf("%02d", 1:26),
                        missing_codes = c("XX", ""), variable = col)

  mapped <- dplyr::case_when(
    is.na(norm) | norm == "XX"  ~ NA_character_,
    norm %in% c("01", "02")     ~ "1",
    TRUE                        ~ "2"
  )

  f <- factor(mapped, levels = .household_size_levels)
  .report(out, f, verbose)
  df[[out]] <- f
  if (col != out) df <- df %>% select(-any_of(col))
  df
}

process_pmr_health4 <- function(df, col = "health_census", out = "health", verbose = TRUE) {
  stopifnot(is.data.frame(df))
  if (!col %in% names(df)) stop(sprintf("Column '%s' not found.", col))

  assert_expected_codes(toupper(trimws(as.character(df[[col]]))),
                        expected = 1:5, missing_codes = c("X", "NA", ""),
                        variable = col)
  code <- .as_num(df[[col]])
  mapped <- dplyr::case_when(
    code == 1 ~ "1",
    code == 2 ~ "2",
    code == 3 ~ "3",
    code %in% c(4, 5) ~ "4",
    TRUE ~ NA_character_
  )

  f <- factor(mapped, levels = .health_levels)
  .report(out, f, verbose)
  df[[out]] <- f
  if (col != out) df <- df %>% select(-any_of(col))
  df
}

NSSEC_OPERATIONAL_TO_CLASS8 <- c(
  "1" = "1",
  "2" = "1",
  "3.1" = "1",
  "3.2" = "1",
  "3.3" = "1",
  "3.4" = "1",
  "4.1" = "2",
  "4.2" = "2",
  "4.3" = "2",
  "4.4" = "2",
  "5" = "2",
  "6" = "2",
  "7.1" = "3", "7.2" = "3",
  "7.3" = "3", "7.4" = "3",
  "8.1" = "4",
  "8.2" = "4",
  "9.1" = "4",
  "9.2" = "4",
  "10" = "5",
  "11.1" = "5",
  "11.2" = "5",
  "12.1" = "6", "12.2" = "6",
  "12.3" = "6", "12.4" = "6",
  "12.5" = "6", "12.6" = "6",
  "12.7" = "6",
  "13.1" = "7", "13.2" = "7", "13.3" = "7",
  "13.4" = "7", "13.5" = "7",
  "14.1" = "8",
  "14.2" = "8",
  "15" = NA_character_, "16" = NA_character_, "17" = NA_character_,
  "XXXX" = NA_character_
)

NSSEC_CLASS8_LEVELS <- as.character(1:8)

process_pmr_nssec <- function(df, col = "nssec_census", out = "nssec_class8", verbose = TRUE) {
  stopifnot(is.data.frame(df))
  if (!col %in% names(df)) stop(sprintf("Column '%s' not found.", col))

  code <- trimws(as.character(df[[col]]))
  code[code %in% c("", "NA")] <- NA_character_
  num <- suppressWarnings(as.numeric(code))
  whole <- !is.na(num) & num == floor(num)
  code[whole] <- as.character(as.integer(num[whole]))

  unknown <- setdiff(unique(code[!is.na(code)]), names(NSSEC_OPERATIONAL_TO_CLASS8))
  if (length(unknown)) {
    stop(sprintf("Unexpected %s codes: %s", col, paste(unknown, collapse = ", ")),
         call. = FALSE)
  }

  f <- factor(unname(NSSEC_OPERATIONAL_TO_CLASS8[code]), levels = NSSEC_CLASS8_LEVELS)
  .report(out, f, verbose)
  df[[out]] <- f
  if (col != out) df <- df %>% select(-any_of(col))
  df
}

pmr_prepare <- function(df, verbose = TRUE) {
  stopifnot(is.data.frame(df))

  out <- tibble::tibble(
    w = pmr_sampling_weight(df, verbose),
    dod_deaths = pmr_dod_deaths(df, verbose),
    age = pmr_age(df, verbose),
    sex = pmr_sex(df, verbose),
    LSOA11CD = df$LSOA11CD,
    ethpuk11_census = df$ethpuk11_census,
    health_census = df$health_census,
    tenhuk11_census = df$tenhuk11_census,
    ecocatpuk11_census = df$ecocatpuk11_census,
    hlqpuk11_census = df$hlqpuk11_census,
    hhchuk11_census = df$hhchuk11_census,
    ruralurban_code_census = df$ruralurban_code_census,
    nssec_census = df$nssec_census,
    underlying_cause_of_death_icd = df$fic10und_deaths
  )

  .residence <- trimws(as.character(df$residence_type_census))
  .residence[.residence %in% c("NA", "")] <- NA_character_
  assert_expected_codes(.residence, expected = c("C", "H"),
                        variable = "residence_type_census")

  out <- out %>%
    process_pmr_icd_missing(
      cols = "underlying_cause_of_death_icd",
      verbose = verbose
    ) %>%
    process_pmr_ethnicity(col = "ethpuk11_census", out = "ethnicity") %>%
    process_pmr_health4(col = "health_census", out = "health", verbose = verbose) %>%
    process_pmr_tenure(col = "tenhuk11_census", out = "tenure", verbose = verbose) %>%
    process_pmr_econstatus(col = "ecocatpuk11_census", out = "econstatus", verbose = verbose) %>%
    process_pmr_education(col = "hlqpuk11_census", out_level = "education", verbose = verbose) %>%
    process_pmr_household_size(col = "hhchuk11_census", out = "household_size", verbose = verbose) %>%
    process_pmr_ruralurban(col = "ruralurban_code_census", out = "ruralurban", verbose = verbose) %>%
    process_pmr_nssec(col = "nssec_census", out = "nssec_class8", verbose = verbose) %>%
    process_imd_decile(lsoa11_col = "LSOA11CD", out = "imd_decile", verbose = verbose)

  out
}
