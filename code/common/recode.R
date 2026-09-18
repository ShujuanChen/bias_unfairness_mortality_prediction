# Guard for value mappings.

assert_expected_codes <- function(x, expected, missing_codes = character(0),
                                  variable = "variable") {
  seen <- unique(as.character(x))
  seen <- seen[!is.na(seen)]
  unexpected <- setdiff(seen, c(as.character(expected), as.character(missing_codes)))
  if (length(unexpected)) {
    stop(sprintf(
      paste0("Unexpected codes in '%s': %s.\n",
             "Give them a branch in the harmonisation, or name them as ",
             "non-response. They must not become missing by falling through."),
      variable, paste(shQuote(sort(unexpected)), collapse = ", ")), call. = FALSE)
  }
  invisible(TRUE)
}

assert_count_field <- function(x, missing_codes = character(0),
                               variable = "variable") {
  chr <- as.character(x)
  num <- suppressWarnings(as.numeric(chr))
  assert_expected_codes(chr[is.na(num)], expected = character(0),
                        missing_codes = missing_codes, variable = variable)
  assert_expected_codes(chr[!is.na(num) & num < 1], expected = character(0),
                        missing_codes = missing_codes, variable = variable)
}
