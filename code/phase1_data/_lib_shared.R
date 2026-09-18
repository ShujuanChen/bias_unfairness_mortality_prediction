# Definitions used by both the UKB and the PMR harmonisation.

suppressPackageStartupMessages({
  library(dplyr)
})

HARMONISED_ETHNICITY_LEVELS <- c("1", "2", "3", "4", "5", "6")

validate_harmonised_ethnicity <- function(x, context = "ethnicity",
                                          levels = HARMONISED_ETHNICITY_LEVELS) {
  x_chr <- as.character(x)
  x_chr[is.na(x_chr)] <- NA_character_
  x_chr[trimws(x_chr) == ""] <- NA_character_
  unexpected <- sort(unique(x_chr[!is.na(x_chr) & !x_chr %in% levels]))
  if (length(unexpected) > 0L) {
    stop(
      context, " contains unexpected ethnicity values: ",
      paste(shQuote(unexpected), collapse = ", "),
      ". Expected only: ", paste(levels, collapse = ", "),
      call. = FALSE
    )
  }
  invisible(TRUE)
}

coerce_harmonised_ethnicity <- function(x, context = "ethnicity",
                                        levels = HARMONISED_ETHNICITY_LEVELS) {
  validate_harmonised_ethnicity(x, context = context, levels = levels)
  factor(as.character(x), levels = levels)
}

.tenure_levels <- c("1", "2", "3", "4", "5", "6", "7")
.sex_levels <- c("1", "2")
.econstatus_levels <- c("1", "2", "3", "4")
.education_levels <- c("1", "2", "3", "4")
.ruralurban_levels <- c("1", "2", "3", "4")
.health_levels <- c("1", "2", "3", "4")
.household_size_levels <- c("1", "2")

.as_num <- function(x) suppressWarnings(as.numeric(as.character(x)))

.report <- function(name, x, verbose) {
  if (!isTRUE(verbose)) return(invisible(NULL))
  n_na <- sum(is.na(x))
  lev  <- if (is.factor(x)) paste(levels(x), collapse = " | ") else ""
  message(sprintf("[%-28s] N=%d, NA=%d%s",
                  name, length(x), n_na,
                  if (lev != "") paste0(", levels: ", lev) else ""))
  invisible(NULL)
}

deprivation_path <- function(stem) {
  for (ext in c("xls", "xlsx")) {
    p <- path_data("imd", paste0(stem, ".", ext))
    if (file.exists(p)) return(p)
  }
  path_data("imd", paste0(stem, ".xls"))
}

process_imd_decile <- function(df,
                               lsoa11_col = "LSOA11CD",
                               out = "imd_decile",
                               eimd_path = deprivation_path("EIMD2010"),
                               wimd_path = deprivation_path("WIMD2011"),
                               lookup_path = path_data("lookup", "LSOA01_LSOA11_LAD11_Lookup_EW.csv"),
                               verbose = TRUE) {
  stopifnot(is.data.frame(df))
  if (!lsoa11_col %in% names(df)) stop(sprintf("Column '%s' not found.", lsoa11_col))
  required_files <- c(eimd_path, wimd_path, lookup_path)
  missing_files <- required_files[!file.exists(required_files)]
  if (length(missing_files) > 0L) {
    stop("process_imd_decile: missing deprivation or lookup files: ",
         paste(missing_files, collapse = ", "), call. = FALSE)
  }

  eimd_raw <- readxl::read_excel(eimd_path, sheet = 2)
  wimd_raw <- readxl::read_excel(wimd_path, sheet = 2)
  lu <- read_table(lookup_path, stringsAsFactors = FALSE)

  need_lu <- c("LSOA01CD", "LSOA11CD")
  if (!all(need_lu %in% names(lu))) {
    stop("Lookup file must contain columns: ", paste(need_lu, collapse = ", "))
  }

  eimd_cols <- c("LSOA CODE", "IMD SCORE")
  wimd_cols <- c("LSOA Code", "WIMD 2011 score")
  if (!all(eimd_cols %in% names(eimd_raw))) {
    stop("EIMD sheet must contain columns: ", paste(eimd_cols, collapse = ", "))
  }
  if (!all(wimd_cols %in% names(wimd_raw))) {
    stop("WIMD sheet must contain columns: ", paste(wimd_cols, collapse = ", "))
  }

  e01 <- eimd_raw %>%
    transmute(
      lsoa01 = as.character(.data[["LSOA CODE"]]),
      imd_score = suppressWarnings(as.numeric(.data[["IMD SCORE"]])),
      country = "E"
    )
  w01 <- wimd_raw %>%
    transmute(
      lsoa01 = as.character(.data[["LSOA Code"]]),
      imd_score = suppressWarnings(as.numeric(.data[["WIMD 2011 score"]])),
      country = "W"
    )

  imd01 <- bind_rows(e01, w01) %>%
    filter(!is.na(lsoa01), nzchar(lsoa01), !is.na(imd_score))

  pairs <- imd01 %>%
    inner_join(unique(lu[, need_lu]), by = c("lsoa01" = "LSOA01CD")) %>%
    transmute(lsoa11 = as.character(LSOA11CD), imd_score, country)

  agg11 <- pairs %>%
    group_by(lsoa11, country) %>%
    summarise(imd_score = mean(imd_score, na.rm = TRUE), .groups = "drop") %>%
    group_by(country) %>%
    mutate(imd_decile = ntile(desc(imd_score), 10L)) %>%
    ungroup()

  out_map <- agg11 %>% select(lsoa11, imd_decile)

  if (anyDuplicated(out_map$lsoa11)) {
    stop("process_imd_decile: the IMD lookup has duplicated LSOA11 codes, so the ",
         "join would multiply rows.", call. = FALSE)
  }
  n_before <- nrow(df)
  df[[lsoa11_col]] <- as.character(df[[lsoa11_col]])
  df <- df %>% left_join(out_map, by = setNames("lsoa11", lsoa11_col))
  if (nrow(df) != n_before) {
    stop(sprintf("process_imd_decile: the join changed the row count from %d to %d.",
                 n_before, nrow(df)), call. = FALSE)
  }
  df[[out]] <- df$imd_decile
  if (out != "imd_decile") df$imd_decile <- NULL

  .report(out, df[[out]], verbose)
  df
}
