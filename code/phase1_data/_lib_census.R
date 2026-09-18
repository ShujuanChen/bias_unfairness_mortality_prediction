# Harmonise the 2011 Census 5% microdata and the UKB extract to the Census-weighting scheme.

# ---- helpers -----------------------------------------------------------------

years_of_education <- function(degree_code, age_finished) {
  years <- dplyr::case_when(
    degree_code == -7 ~ 7, degree_code == 3 ~ 10, degree_code == 4 ~ 10,
    degree_code == 2 ~ 13, degree_code == 1 ~ 20
  )
  is5 <- degree_code == 5 & !is.na(degree_code)
  years[is5] <- age_finished[is5] - 5
  years[is5 & years >= 19] <- 19
  is6 <- degree_code == 6 & !is.na(degree_code)
  years[is6] <- age_finished[is6] - 5
  years[is6 & years >= 15] <- 15
  years[years <= 7] <- 7
  years
}

census_health4 <- function(health) {
  dplyr::case_when(health %in% c(5, 4) ~ 1, health == 3 ~ 2,
                   health == 2 ~ 3, health == 1 ~ 4)
}

ukb_health4 <- function(x) {
  dplyr::case_when(x == 4 ~ 1, x == 3 ~ 2, x == 2 ~ 3, x == 1 ~ 4)
}

COLLAPSE_EMPSTAT_PRECEDENCE <- c(1, 5, 2)

collapse_empstat <- function(x) {
  raw <- as.character(x)
  raw[is.na(raw)] <- ""

  tokens <- lapply(strsplit(raw, "|", fixed = TRUE), function(v) {
    v <- trimws(v)
    suppressWarnings(as.numeric(v[nzchar(v)]))
  })
  assert_expected_codes(unlist(tokens), expected = c(1:7, -7), missing_codes = -3,
                        variable = "employment status")
  multi <- vapply(tokens, function(v) sum(!is.na(v)) >= 2, logical(1))

  out <- suppressWarnings(as.numeric(raw))
  out[raw == ""] <- NA_real_
  for (i in which(multi)) {
    v <- tokens[[i]][!is.na(tokens[[i]])]
    v <- v[v != -3]  # a refusal alongside a real answer carries no information
    hit <- COLLAPSE_EMPSTAT_PRECEDENCE[COLLAPSE_EMPSTAT_PRECEDENCE %in% v]
    out[i] <- if (!length(v)) NA_real_ else if (length(hit)) hit[1] else 8
  }
  out[out %in% c(3, 4, 6, 7, -7)] <- 8
  out[out == -3] <- NA_real_
  out
}

.CENSUS_HOUSEHOLD_SIZE_LEVELS <- c("1", "2", "3", "4", "5", "6", "7+")

census_household_size <- function(n) {
  n <- suppressWarnings(as.integer(n))
  n[!is.na(n) & n < 1] <- NA_integer_
  factor(ifelse(is.na(n), NA_character_,
                ifelse(n >= 7, "7+", as.character(n))),
         levels = .CENSUS_HOUSEHOLD_SIZE_LEVELS)
}

census_tenure <- function(code) {
  code <- suppressWarnings(as.integer(code))
  code[!is.na(code) & code == 4L] <- 3L
  code[!is.na(code) & code %in% c(7L, 8L)] <- 6L
  code
}

# ---- Census microdata --------------------------------------------------------

prepare_census_microdata <- function(path) {
  census <- read_table(
    path,
    select = c("country", "ageh", "carsnoc", "ecopuk11",
               "ethnicityew", "health", "hlqupuk11", "sizhuk11", "sex",
               "tenduk11"),
    data.table = FALSE
  )
  n_raw <- nrow(census)

  census$ref_id <- seq_len(nrow(census))

  for (spec in list(
    list("country", 1:2), list("ageh", 1:19), list("carsnoc", 0:4),
    list("ecopuk11", 1:14), list("ethnicityew", 1:13), list("health", 1:5),
    list("hlqupuk11", 10:16), list("sizhuk11", 0:7), list("sex", 1:2),
    list("tenduk11", 0:9))) {
    assert_expected_codes(census[[spec[[1]]]], expected = spec[[2]],
                          missing_codes = c("", "NA"), variable = spec[[1]])
  }

  years <- dplyr::case_when(
    census$hlqupuk11 == 10 ~ 7, census$hlqupuk11 == 11 ~ 10,
    census$hlqupuk11 == 12 ~ 10, census$hlqupuk11 == 13 ~ 12,
    census$hlqupuk11 == 14 ~ 13, census$hlqupuk11 == 15 ~ 20,
    census$hlqupuk11 == 16 ~ 15
  )

  out <- data.frame(
    ref_id = census$ref_id,
    country = census$country,
    ageh = census$ageh,
    age = census_age_band(census$ageh),
    sex = census$sex,
    cars = census$carsnoc,
    education = cut(years, breaks = c(0, 8.5, 11, 17.5, 20),
                    labels = c("1", "2", "3", "4")),
    health = census_health4(census$health),
    household_size = census_household_size(census$sizhuk11),
    ethnicity = dplyr::case_when(
      census$ethnicityew <= 3 ~ 1,
      census$ethnicityew >= 4 & census$ethnicityew <= 5 ~ 2,
      census$ethnicityew %in% c(6, 7, 8, 10) ~ 3,
      census$ethnicityew >= 11 & census$ethnicityew <= 12 ~ 4,
      census$ethnicityew == 9 ~ 5,
      census$ethnicityew == 13 ~ 6
    ),
    econstatus = dplyr::case_when(
      census$ecopuk11 <= 6 ~ 1, census$ecopuk11 == 7 ~ 5,
      census$ecopuk11 == 8 ~ 1, census$ecopuk11 == 9 ~ 8,
      census$ecopuk11 == 10 ~ 2, census$ecopuk11 == 11 ~ 8,
      census$ecopuk11 == 12 ~ 8, census$ecopuk11 == 13 ~ 8,
      census$ecopuk11 == 14 ~ 8
    ),
    tenure = census_tenure(census$tenduk11),
    stringsAsFactors = FALSE
  )

  if (nrow(out) != n_raw) {
    stop("Census harmonisation changed the row count.", call. = FALSE)
  }
  out
}

# ---- UKB extract -------------------------------------------------------------

prepare_ukb_for_census <- function(path, predictors) {
  ukb <- read_table(path, header = TRUE, sep = "auto", data.table = FALSE)
  n_raw <- nrow(ukb)

  assert_expected_codes(ukb$p31, expected = c(0, 1), missing_codes = "",
                        variable = "sex")
  assert_expected_codes(ukb$p21000_i0,
                        expected = c(1:6, 1001:1003, 2001:2004, 3001:3004, 4001:4003),
                        missing_codes = c(-1, -3, ""), variable = "ethnicity")
  assert_count_field(ukb$p709_i0, missing_codes = c(-1, -3, ""),
                     variable = "number in household")
  assert_expected_codes(ukb$p845_i0[!is.na(ukb$p845_i0) & ukb$p845_i0 < 5],
                        expected = character(0), missing_codes = c(-1, -2, -3),
                        variable = "age completed education")
  assert_expected_codes(ukb$p728_i0, expected = 1:5, missing_codes = c(-1, -3, ""),
                        variable = "vehicles in the household")
  assert_expected_codes(ukb$p2178_i0, expected = 1:4,
                        missing_codes = c(-1, -3, ""), variable = "overall health")
  assert_expected_codes(ukb$p680_i0, expected = c(1:6, -7), missing_codes = c(-3, ""),
                        variable = "tenure")

  ethnicity <- ukb$p21000_i0
  is_subgroup <- ethnicity >= 1000 & !is.na(ethnicity)
  ethnicity[is_subgroup] <- round(ethnicity[is_subgroup] / 1000)
  ethnicity[ethnicity <= 0] <- NA

  ukb$p845_i0[ukb$p845_i0 < 0] <- NA

  degree_tokens <- lapply(strsplit(as.character(ukb$p6138_i0), "|", fixed = TRUE), function(v) {
    v <- trimws(v)
    suppressWarnings(as.numeric(v[nzchar(v)]))
  })
  assert_expected_codes(unlist(degree_tokens), expected = c(1:6, -7),
                        missing_codes = -3, variable = "qualifications")
  degree_mat <- t(vapply(degree_tokens, function(v) {
    length(v) <- 6L
    v
  }, numeric(6)))
  years_cols <- vapply(seq_len(6L), function(j) {
    years_of_education(degree_mat[, j], ukb$p845_i0)
  }, numeric(nrow(ukb)))
  years <- apply(years_cols, 1, function(x) if (all(is.na(x))) NA_real_ else max(x, na.rm = TRUE))

  if (!"p52" %in% names(ukb)) stop("Column 'p52' not found.", call. = FALSE)
  raw_month <- trimws(as.character(ukb$p52))
  month <- match(tolower(raw_month), tolower(month.name))
  month[is.na(month)] <- suppressWarnings(as.integer(raw_month))[is.na(month)]
  month[!is.na(month) & (month < 1 | month > 12)] <- NA_integer_
  assert_expected_codes(raw_month[is.na(month)], expected = character(0),
                        missing_codes = "", variable = "month of birth")
  dob <- as.Date(sprintf("%04d-%02d-15", as.integer(ukb$p34), month))

  cars <- ukb$p728_i0
  cars[cars %in% c(-1, -3)] <- NA
  health <- ukb$p2178_i0
  health[health %in% c(-1, -3)] <- NA
  sex <- ukb$p31
  sex[sex == 0] <- 2

  prepared <- data.frame(
    eid = ukb$eid,
    source = "UKB",
    sample_weight = 1,
    age = census_age_band_from_age(
      floor(as.numeric(COHORT_BASELINE - dob) / 365.25)),
    sex = sex,
    cars = cars - 1,
    education = cut(years, breaks = c(0, 8.5, 11, 17.5, 20),
                    labels = c("1", "2", "3", "4")),
    health = ukb_health4(health),
    household_size = census_household_size(ukb$p709_i0),
    ethnicity = ethnicity,
    econstatus = collapse_empstat(ukb$p6142_i0),
    tenure = dplyr::case_when(
      ukb$p680_i0 == 1 ~ 0, ukb$p680_i0 == 2 ~ 1, ukb$p680_i0 == 5 ~ 2,
      ukb$p680_i0 == 3 ~ 3, ukb$p680_i0 == 4 ~ 5, ukb$p680_i0 == 6 ~ 9,
      ukb$p680_i0 == -7 ~ 6
    ),
    stringsAsFactors = FALSE
  )

  if (nrow(prepared) != n_raw) {
    stop("UKB harmonisation for the Census changed the row count.", call. = FALSE)
  }
  prepared[, c("eid", "source", "sample_weight", predictors), drop = FALSE]
}
