# Harmonise the HSE waves and the UKB extract to the survey-weighting scheme.

# ---- declarations ------------------------------------------------------------

.HSE_WAVES <- as.data.frame(rbind(
  c("2006", "UKDA-5809-spss/spss/spss12/hse06ai.sav", "ethinda", "psu"),
  c("2007", "UKDA-6112-spss/spss/spss12/hse07ai.sav", "ethinda", "area"),
  c("2008", "UKDA-6397-spss/spss/spss24/hse08ai.sav", "origin", "psu"),
  c("2009", "UKDA-6732-spss/spss/spss19/hse09ai.sav", "origin", "psu"),
  c("2010", "UKDA-6986-spss/spss/spss19/hse10ai.sav", "origin", "psu")
), stringsAsFactors = FALSE)
names(.HSE_WAVES) <- c("year", "path", "ethnicity_column", "psu_column")

# Weights are already normalised to a mean of one within each wave in raw data.
.HSE_WEIGHT_COLUMN <- "wt_int" 

.HSE_COLUMNS <- as.data.frame(rbind(
  c("sex",               "cat", "sex",     "sex",     "sex",     "sex",     "sex"),
  c("age",               "con", "age",     "age",     "age",     "age",     "age"),
  c("education_age",     "cat", "educend", "educend", "educend", "educend", "educend"),
  c("alcfrequency",      "cat", "dnoft2",  "dnoft2",  "dnoft3",  "dnoft3",  "dnoft3"),
  c("smoking_status",    "cat", "cigst1",  "cigst1",  "cigst1",  "cigst1",  "cigst1"),
  c("income",            "cat", "totinc",  "totinc",  "totinc",  "totinc",  "totinc"),
  c("household_size",    "con", "hhsize",  "hhsizeD", "hhsize",  "hhsize",  "hhsize"),
  c("econstatus", "cat", "econact", "econact", "econact", "econact", "econact"),
  c("health",     "cat", "genhelf", "genhelf", "genhelf", "genhelf", "genhelf"),
  c("height",            "con", "estht",   "estht",   "estht",   "estht",   "estht"),
  c("ruralurban",      "cat", "URINDEW", "URINDEW", "Urban",   "Urban",   "urban"),
  c("weight",            "con", "estwt",   "estwt",   "estwt",   "estwt",   "estwt"),
  c("disability",        "cat", "longill", "longill", "longill", "longill", "longill")
), stringsAsFactors = FALSE)
names(.HSE_COLUMNS) <- c("label", "type", .HSE_WAVES$year)

.UKB_FIELDS <- as.data.frame(rbind(
  c("sex",               "31",    "bin", "yes"),
  c("age",               "21022", "con", "yes"),
  c("education_age",     "845",   "con", "yes"),
  c("alcfrequency",      "1558",  "cat", "yes"),
  c("smoking_status",    "20116", "cat", "yes"),
  c("income",            "738",   "cat", "yes"),
  c("household_size",    "709",   "cat", "yes"),
  c("econstatus", "6142",  "cat", "yes"),
  c("bmi",               NA,      "con", "yes"),
  c("bmi_cat",           NA,      "cat", "yes"),
  c("health",     "2178",  "cat", "yes"),
  c("height",            "50",    "con", "yes"),
  c("ruralurban",      "20118", "cat", "yes"),
  c("weight",            "21002", "con", "yes"),
  c("ethnic_background", "21000", "cat", "no"),
  c("education_degree",  "6138",  "cat", "no"),
  c("disability",        "2188",  "cat", "yes")
), stringsAsFactors = FALSE)
names(.UKB_FIELDS) <- c("label", "field", "type", "predictor")
.UKB_FIELDS$predictor <- .UKB_FIELDS$predictor == "yes"


hse_predictor_info <- function() {
  info <- .UKB_FIELDS[.UKB_FIELDS$predictor, c("label", "type"), drop = FALSE]
  rbind(info, data.frame(label = "ethnicity", type = "cat"),
        make.row.names = FALSE)
}

hse_predictors <- function() {
  labels <- hse_predictor_info()$label
  declared <- framework_setting("phase1", "predictors", "hse_participation")
  if (!identical(labels, declared)) {
    stop("The HSE auxiliaries in _lib_hse.R and in framework_config.json differ.\n  code:   ",
         paste(labels, collapse = ", "), "\n  config: ",
         paste(declared, collapse = ", "), call. = FALSE)
  }
  labels
}

# ---- helpers -----------------------------------------------------------------
HSE_MISSING_CODES <- c(-1, -2, -6, -7, -8, -9)

hse_answer <- function(x, expected, variable) {
  code <- suppressWarnings(as.numeric(as.character(x)))
  code[code %in% HSE_MISSING_CODES] <- NA_real_
  assert_expected_codes(code, expected = expected, variable = variable)
  as.character(code)
}

labelled_to_code <- function(x) as.character(haven::zap_labels(x))

relabel_factor <- function(x, mapping, lvls = NULL, variable = "variable") {
  assert_expected_codes(x, expected = names(mapping), variable = variable)
  relabelled <- plyr::revalue(as.factor(x), mapping, warn_missing = FALSE)
  if (is.null(lvls)) relabelled else factor(relabelled, levels = lvls)
}

ukb_drop_negative <- function(x, variable) {
  assert_expected_codes(x[!is.na(x) & x < 0], expected = character(0),
                        missing_codes = c(-1, -3), variable = variable)
  ifelse(x < 0, NA, x)
}

split_codes <- function(x) {
  lapply(strsplit(as.character(x), "|", fixed = TRUE), function(tok) {
    tok <- trimws(tok)
    tok[nzchar(tok)]
  })
}

resolve_ukb_column <- function(field, label, available) {
  candidates <- c(paste0("p", field, "_i0"), paste0("p", field))
  hit <- candidates[candidates %in% available]
  if (!length(hit)) stop(sprintf("Cannot map '%s' (field %s) to a UKB column.", label, field))
  hit[1]
}

# ---- variable mappings -------------------------------------------------------
.SEX_LEVELS <- c("1", "2")

hse_recode_sex <- function(x) {
  relabel_factor(hse_answer(x, expected = 1:2, variable = "sex"),
                 c("1" = "1", "2" = "2"),
                 lvls = .SEX_LEVELS, variable = "sex")
}

ukb_recode_sex <- function(x) {
  relabel_factor(x, c("0" = "2", "1" = "1"),
                 lvls = .SEX_LEVELS, variable = "sex")
}

collapse_degree_codes <- function(x) {
  tokens <- split_codes(x)
  assert_expected_codes(unlist(tokens), expected = c(1:6, -7),
                        missing_codes = -3, variable = "education_degree")
  vapply(tokens, function(tok) {
    tok <- tok[!tok %in% c("-3")]
    if (!length(tok)) return(NA_character_)
    if ("1" %in% tok) "1" else "0"
  }, character(1))
}

hse_recode_education_age <- function(x) {
  code <- suppressWarnings(as.numeric(as.character(x)))
  code[code %in% HSE_MISSING_CODES] <- NA_real_
  assert_expected_codes(code, expected = 2:8, missing_codes = 1,
                        variable = "education_age")
  years <- rep(NA_real_, length(code))
  years[code == 2] <- 14
  years[code == 3] <- 14
  years[code == 4] <- 15
  years[code == 5] <- 16
  years[code == 6] <- 17
  years[code == 7] <- 18
  years[code == 8] <- 19
  years
}

ukb_recode_education_age <- function(x, degree) {
  degree <- suppressWarnings(as.numeric(degree))
  assert_expected_codes(degree, expected = c(0, 1), variable = "education_degree")
  raw <- suppressWarnings(as.numeric(x))
  assert_expected_codes(raw[!is.na(raw) & raw < 5], expected = character(0),
                        missing_codes = c(-1, -2, -3), variable = "education_age")

  years <- ifelse(!is.na(degree) & degree == 1, 20, raw)
  years[!is.na(raw) & raw == -2] <- 14
  years <- ifelse(!is.na(years) & years < 0, NA, years)
  years <- ifelse(years <= 14, 14, years)
  ifelse(years >= 19, 19, years)
}

.ALCOHOL_LEVELS <- c("8", "6", "5", "4", "3", "1")

hse_recode_alcohol <- function(x) {
  relabel_factor(
    hse_answer(x, expected = 1:8, variable = "alcfrequency"),
    c("8" = "8", "7" = "6", "6" = "6", "5" = "5", "4" = "4", "3" = "3",
      "2" = "1", "1" = "1"),
    lvls = .ALCOHOL_LEVELS, variable = "alcfrequency"
  )
}

ukb_recode_alcohol <- function(x) {
  relabel_factor(x, c(
    "-3" = NA, "6" = "8", "5" = "6", "4" = "5", "3" = "4", "2" = "3", "1" = "1"
  ), lvls = .ALCOHOL_LEVELS, variable = "alcfrequency")
}

.SMOKING_LEVELS <- c("1", "2", "4")

hse_recode_smoking <- function(x) {
  relabel_factor(
    hse_answer(x, expected = 1:4, variable = "smoking_status"),
    c("1" = "1", "2" = "2", "3" = "2", "4" = "4"),
    lvls = .SMOKING_LEVELS, variable = "smoking_status"
  )
}

ukb_recode_smoking <- function(x) {
  relabel_factor(x, c("-3" = NA, "0" = "1", "1" = "2", "2" = "4"),
                 lvls = .SMOKING_LEVELS, variable = "smoking_status")
}

.INCOME_LEVELS <- c("96", "1", "11", "16", "21", "26")

hse_recode_income <- function(x) {
  groups <- c(rep("1", 10), rep("11", 5), rep("16", 5), rep("21", 5),
              rep("26", 6), "96", "96")
  mapping <- stats::setNames(groups, as.character(c(1:31, 96, 97)))
  relabel_factor(hse_answer(x, expected = c(1:31, 96, 97), variable = "income"),
                 mapping, lvls = .INCOME_LEVELS, variable = "income")
}

ukb_recode_income <- function(x) {
  relabel_factor(x, c(
    "-1" = "96", "-3" = "96", "1" = "1", "2" = "11",
    "3" = "16", "4" = "21", "5" = "26"
  ), lvls = .INCOME_LEVELS, variable = "income")
}

household_size_levels <- function(n) {
  n <- suppressWarnings(as.numeric(as.character(n)))
  n[!is.na(n) & n < 1] <- NA_real_
  factor(ifelse(is.na(n), NA_character_, as.character(pmin(n, 7))),
         levels = as.character(1:7))
}

hse_recode_household_size <- function(x) {
  assert_count_field(x, variable = "household_size")
  household_size_levels(x)
}

ukb_recode_household_size <- function(x) {
  assert_count_field(x, missing_codes = c(-1, -3), variable = "household_size")
  household_size_levels(x)
}

.EMPLOYMENT_LEVELS <- c("2", "1", "4", "3")

hse_recode_employment <- function(x) {
  relabel_factor(
    hse_answer(x, expected = 1:4, variable = "econstatus"),
    c("1" = "1", "2" = "2", "3" = "3", "4" = "4"),
    lvls = .EMPLOYMENT_LEVELS, variable = "econstatus"
  )
}

collapse_employment_codes <- function(x) {
  tokens <- split_codes(x)
  assert_expected_codes(unlist(tokens), expected = c(1:7, -7),
                        missing_codes = -3, variable = "econstatus")
  vapply(tokens, function(tok) {
    tok <- tok[tok != "-3"]
    real <- tok[tok != "-7"]
    if (!length(real)) return(if ("-7" %in% tok) "-7" else NA_character_)
    if ("1" %in% real) return("1")
    if ("5" %in% real) return("5")
    if ("2" %in% real) return("2")
    if (any(real %in% c("3", "4", "6", "7"))) return("3")
    NA_character_
  }, character(1))
}

ukb_recode_employment <- function(x) {
  relabel_factor(x, c(
    "-7" = "4", "-3" = NA, "1" = "1", "2" = "3", "3" = "4", "4" = "4",
    "5" = "2", "6" = "4", "7" = "4"
  ), lvls = .EMPLOYMENT_LEVELS, variable = "econstatus")
}

height_to_whole_cm <- function(x) round(x)

bmi_from_measures <- function(height, weight) {
  h <- ifelse(!is.na(height) & height > 0, height, NA_real_)
  w <- ifelse(!is.na(weight) & weight > 0, weight, NA_real_)
  w / (h / 100)^2
}

.BMI_CATEGORY_LEVELS <- c("1", "2", "3", "4")

bmi_category <- function(bmi) {
  band <- ifelse(bmi < 18.5, "1", NA_character_)
  band <- ifelse(bmi >= 18.5 & bmi < 25, "2", band)
  band <- ifelse(bmi >= 25 & bmi < 30, "3", band)
  band <- ifelse(bmi >= 30, "4", band)
  factor(band, levels = .BMI_CATEGORY_LEVELS)
}

.HEALTH_LEVELS <- c("4", "3", "2", "1")

hse_recode_health4 <- function(x) {
  relabel_factor(
    hse_answer(x, expected = 1:5, variable = "health"),
    c("1" = "1", "2" = "2", "3" = "3", "4" = "4", "5" = "4"),
    lvls = .HEALTH_LEVELS, variable = "health"
  )
}

ukb_recode_health4 <- function(x) {
  relabel_factor(x, c("-1" = NA, "-3" = NA, "1" = "1", "2" = "2",
                      "3" = "3", "4" = "4"),
                 lvls = .HEALTH_LEVELS, variable = "health")
}

.URBANISATION_LEVELS <- c("3", "2", "1")

.URINDEW_TO_THREE <- c("1" = "1", "5" = "1", "2" = "2", "6" = "2",
                       "3" = "3", "4" = "3", "7" = "3", "8" = "3")

hse_ruralurban_to_three <- function(x, column) {
  if (identical(column, "URINDEW")) {
    code <- hse_answer(x, expected = 1:8, variable = "ruralurban")
    return(unname(.URINDEW_TO_THREE[code]))
  }
  hse_answer(x, expected = 1:3, variable = "ruralurban")
}

hse_recode_ruralurban <- function(x) {
  relabel_factor(as.character(x), c("1" = "1", "2" = "2", "3" = "3"),
                 lvls = .URBANISATION_LEVELS, variable = "ruralurban")
}

ukb_recode_ruralurban <- function(x) {
  relabel_factor(
    x,
    c(
      "1" = "1", "2" = "2", "3" = "3", "4" = "3",
      "5" = "1", "6" = "2", "7" = "3", "8" = "3",
      "9" = NA, "10" = NA, "11" = NA, "12" = NA, "13" = NA, "14" = NA,
      "15" = NA, "16" = NA, "17" = NA, "18" = NA
    ),
    lvls = .URBANISATION_LEVELS, variable = "ruralurban"
  )
}

.DISABILITY_LEVELS <- c("2", "1")

hse_recode_disability <- function(x) {
  relabel_factor(hse_answer(x, expected = 1:2, variable = "disability"),
                 c("1" = "1", "2" = "2"),
                 lvls = .DISABILITY_LEVELS, variable = "disability")
}

ukb_recode_disability <- function(x) {
  relabel_factor(x, c("-1" = NA, "-3" = NA, "0" = "2", "1" = "1"),
                 lvls = .DISABILITY_LEVELS, variable = "disability")
}

.HSE_ETHNICITY_LEVELS <- c("1", "2", "3", "4", "5")

.ORIGIN_TO_BROAD <- c("1" = "1", "2" = "1", "3" = "1",
                      "4" = "2", "5" = "2", "6" = "2", "7" = "2",
                      "8" = "3", "9" = "3", "10" = "3", "11" = "3",
                      "12" = "4", "13" = "4", "14" = "4",
                      "15" = "5", "16" = "5")

hse_ethnicity_to_broad <- function(x, column) {
  if (identical(column, "origin")) {
    code <- hse_answer(x, expected = 1:16, variable = "ethnicity_detail")
    return(unname(.ORIGIN_TO_BROAD[code]))
  }
  hse_answer(x, expected = 1:5, variable = "ethnicity_detail")
}

hse_recode_ethnicity <- function(broad) {
  assert_expected_codes(broad, expected = .HSE_ETHNICITY_LEVELS,
                        variable = "ethnicity")
  factor(broad, levels = .HSE_ETHNICITY_LEVELS)
}

.UKB_ETHNICITY <- list(
  "1" = c(1, 1001, 1002, 1003),
  "2" = c(2, 2001, 2002, 2003, 2004),
  "3" = c(3, 3001, 3002, 3003, 3004),
  "4" = c(4, 4001, 4002, 4003),
  "5" = c(5, 6)
)

ukb_recode_ethnicity <- function(codes) {
  numeric_codes <- suppressWarnings(as.numeric(as.character(codes)))
  assert_expected_codes(numeric_codes,
                        expected = unlist(.UKB_ETHNICITY, use.names = FALSE),
                        missing_codes = c(-1, -3), variable = "ethnic_background")
  group <- rep(NA_character_, length(numeric_codes))
  for (g in names(.UKB_ETHNICITY)) group[numeric_codes %in% .UKB_ETHNICITY[[g]]] <- g
  factor(group, levels = .HSE_ETHNICITY_LEVELS)
}

# ---- survey waves ------------------------------------------------------------

.read_hse_wave <- function(dir, i) {
  year <- .HSE_WAVES$year[i]
  sav <- haven::read_sav(input_path(dir, .HSE_WAVES$path[i],
                                    label = sprintf("HSE %s wave", year)))

  wave <- lapply(seq_len(nrow(.HSE_COLUMNS)), function(j) {
    column <- sav[[.HSE_COLUMNS[[year]][j]]]
    if (identical(.HSE_COLUMNS$type[j], "con")) {
      suppressWarnings(as.numeric(column))
    } else {
      labelled_to_code(column)
    }
  })

  names(wave) <- .HSE_COLUMNS$label
  wave <- as.data.frame(wave, stringsAsFactors = FALSE)

  wave$ruralurban <- hse_ruralurban_to_three(
    wave$ruralurban, .HSE_COLUMNS[[year]][match("ruralurban", .HSE_COLUMNS$label)])

  wave$weight_individual <- suppressWarnings(as.numeric(sav[[.HSE_WEIGHT_COLUMN]]))

  psu_raw <- sav[[.HSE_WAVES$psu_column[i]]]
  if (is.null(psu_raw)) {
    stop("HSE ", year, " has no column ", .HSE_WAVES$psu_column[i],
         ", which the design-respecting resample needs.", call. = FALSE)
  }
  psu_txt <- as.character(haven::zap_labels(psu_raw))
  if (all(is.na(psu_txt))) {
    stop("HSE ", year, " sampling unit ", .HSE_WAVES$psu_column[i],
         " read as all missing.", call. = FALSE)
  }
  wave$psu <- paste0(year, "_", psu_txt)
  wave$ethnicity_broad <- hse_ethnicity_to_broad(
    sav[[.HSE_WAVES$ethnicity_column[i]]], .HSE_WAVES$ethnicity_column[i])
  wave$ethnicity <- wave$ethnicity_broad
  wave$year <- year
  wave
}

prepare_hse_waves <- function(dir) {
  waves <- lapply(seq_len(nrow(.HSE_WAVES)), function(i) .read_hse_wave(dir, i))

  n_raw <- sum(vapply(waves, nrow, integer(1)))
  hse <- as.data.frame(data.table::rbindlist(waves, fill = TRUE))
  if (nrow(hse) != n_raw) {
    stop("Stacking the survey waves changed the row count.", call. = FALSE)
  }

  hse$sex <- hse_recode_sex(hse$sex)
  hse$education_age <- hse_recode_education_age(hse$education_age)
  hse$alcfrequency <- hse_recode_alcohol(hse$alcfrequency)
  hse$smoking_status <- hse_recode_smoking(hse$smoking_status)
  hse$income <- hse_recode_income(hse$income)
  hse$household_size <- hse_recode_household_size(hse$household_size)
  hse$econstatus <- hse_recode_employment(hse$econstatus)
  hse$health <- hse_recode_health4(hse$health)
  hse$ruralurban <- hse_recode_ruralurban(hse$ruralurban)
  hse$ethnicity <- hse_recode_ethnicity(hse$ethnicity)
  hse$disability <- hse_recode_disability(hse$disability)
  hse$height <- height_to_whole_cm(hse$height)
  hse$bmi <- bmi_from_measures(hse$height, hse$weight)
  hse$bmi_cat <- bmi_category(hse$bmi)
  hse$ref_id <- seq_len(nrow(hse))
  hse
}

# ---- UKB extract -------------------------------------------------------------

prepare_ukb_for_hse <- function(path) {
  raw <- read_table(path, na.strings = c("", "NA"))
  n_raw <- nrow(raw)

  with_field <- .UKB_FIELDS[!is.na(.UKB_FIELDS$field), , drop = FALSE]
  src <- data.frame(eid = raw$eid, check.names = FALSE)
  for (i in seq_len(nrow(with_field))) {
    label <- with_field$label[i]
    src[[label]] <- raw[[resolve_ukb_column(with_field$field[i], label, names(raw))]]
  }

  degree <- collapse_degree_codes(src$education_degree)
  ethnic_background <- ukb_drop_negative(src$ethnic_background, "ethnic_background")
  height <- height_to_whole_cm(ukb_drop_negative(src$height, "height"))
  bmi <- bmi_from_measures(height, src$weight)

  ukb <- data.frame(
    eid                  = src$eid,
    sex                  = ukb_recode_sex(src$sex),
    age                  = ukb_drop_negative(src$age, "age"),
    education_age        = ukb_recode_education_age(src$education_age, degree),
    alcfrequency         = ukb_recode_alcohol(src$alcfrequency),
    smoking_status       = ukb_recode_smoking(src$smoking_status),
    income               = ukb_recode_income(src$income),
    household_size       = ukb_recode_household_size(src$household_size),
    econstatus    = ukb_recode_employment(collapse_employment_codes(src$econstatus)),
    bmi                  = bmi,
    bmi_cat              = bmi_category(bmi),
    health        = ukb_recode_health4(src$health),
    height               = height,
    ruralurban         = ukb_recode_ruralurban(src$ruralurban),
    weight               = ukb_drop_negative(src$weight, "weight"),
    ethnic_background    = ethnic_background,
    disability           = ukb_recode_disability(src$disability),
    weight_individual    = 1,
    ethnicity = ukb_recode_ethnicity(ethnic_background),
    stringsAsFactors     = FALSE
  )

  if (nrow(ukb) != n_raw) {
    stop("UKB harmonisation for the survey changed the row count.", call. = FALSE)
  }
  ukb
}
