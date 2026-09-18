# Join the LSOA to Region lookup onto a harmonised dataset by LSOA11CD.

suppressPackageStartupMessages({
  library(dplyr)
})

.ukb_.read_lookup_csv <- function(path) {
  read_table(path, stringsAsFactors = FALSE, check.names = FALSE)
}

.ukb_.dedupe_lsoa <- function(df, key = "LSOA11CD", cols, name_for_errors = "lookup") {
  stopifnot(key %in% names(df))
  stopifnot(all(cols %in% names(df)))

  out <- df %>%
    transmute(
      LSOA11CD = trimws(.data[[key]]),
      across(all_of(cols), ~ trimws(.x))
    ) %>%
    filter(!is.na(LSOA11CD), LSOA11CD != "") %>%
    distinct()

  dup_keys <- out %>%
    count(LSOA11CD, name = "n") %>%
    filter(n > 1)
  if (nrow(dup_keys) > 0) {
    stop(
      sprintf(
        "Geography merge aborted: %s has duplicated LSOA11CD after column selection/distinct (showing first 10): %s",
        name_for_errors,
        paste(utils::head(dup_keys$LSOA11CD, 10), collapse = ", ")
      )
    )
  }

  out
}

merge_geography <- function(
  df,
  lsoa_to_region_path = path_data("lookup", "Lower_Layer_Super_Output_Area_(2011)_to_Built-up_Area_Sub-division_to_Built-up_Area_to_Local_Authority_District_to_Region_(December_2011)_Lookup_in_England_and_Wales.csv"),
  verbose = TRUE
) {
  stopifnot(is.data.frame(df))
  if (!"LSOA11CD" %in% names(df)) stop("merge_geography(): the data is missing `LSOA11CD`.")

  if (isTRUE(verbose)) message("Geography: reading LSOA->Region lookup…")
  reg_lookup <- .ukb_.read_lookup_csv(lsoa_to_region_path)
  needed2 <- c("LSOA11CD", "RGN11CD")
  if (!all(needed2 %in% names(reg_lookup))) {
    stop(
      "Region lookup missing required columns. Need: ",
      paste(needed2, collapse = ", "),
      ". Found: ",
      paste(names(reg_lookup), collapse = ", ")
    )
  }

  lsoa_reg <- .ukb_.dedupe_lsoa(
    reg_lookup,
    key = "LSOA11CD",
    cols = c("RGN11CD"),
    name_for_errors = "LSOA->Region lookup"
  )

  out <- df %>%
    mutate(LSOA11CD = trimws(as.character(LSOA11CD))) %>%
    left_join(lsoa_reg, by = "LSOA11CD")

  if (isTRUE(verbose)) {
    message("UKB geography merge complete: added RGN11CD.")
  }

  out
}

