# Reading and writing the pipeline's tables. Mirrors code/common/tabular.py.

.ID_COLUMNS <- c("eid", "participant_id", "pmr_id", "ref_id")
.ISO_DATE <- "^[0-9]{4}-[0-9]{2}-[0-9]{2}$"

.as_date <- function(days) structure(as.numeric(days), class = "Date")

.restore_dates <- function(out) {
  for (col in names(out)) {
    x <- out[[col]]
    if (inherits(x, "IDate")) {
      out[[col]] <- .as_date(unclass(x))
      next
    }
    if (!is.character(x)) next
    seen <- x[!is.na(x)]
    if (!length(seen)) next
    if (all(grepl(.ISO_DATE, seen))) out[[col]] <- .as_date(unclass(as.Date(x)))
  }
  out
}

read_table <- function(path, select = NULL, ...) {
  if (!file.exists(path)) stop("No table at ", path, call. = FALSE)
  if (grepl("\\.xlsx$", path)) {
    out <- as.data.frame(readxl::read_excel(path), stringsAsFactors = FALSE)
    if (!is.null(select)) out <- out[, select, drop = FALSE]
    return(out)
  }
  args <- list(...)
  args$file <- path
  args$data.table <- FALSE
  args$showProgress <- FALSE
  args$stringsAsFactors <- NULL
  if (is.null(args$na.strings)) args$na.strings <- "NA"
  args$keepLeadingZeros <- TRUE
  if (!is.null(select)) args$select <- select
  out <- as.data.frame(do.call(data.table::fread, args), stringsAsFactors = FALSE)
  out <- .restore_dates(out)
  for (col in intersect(.ID_COLUMNS, names(out))) {
    if (is.integer(out[[col]])) out[[col]] <- as.numeric(out[[col]])
  }
  out
}

write_table <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  x <- as.data.frame(x)

  if (grepl("\\.xlsx$", path)) {
    writexl::write_xlsx(x, path)
    return(invisible(path))
  }

  for (col in names(x)) {
    v <- x[[col]]
    if (is.factor(v)) v <- levels(v)
    if (is.character(v) && any(grepl('"', v, fixed = TRUE))) {
      stop("Column ", col, " of ", basename(path), " holds a double quote, ",
           "which this format cannot carry.", call. = FALSE)
    }
    if (is.double(x[[col]]) && !inherits(x[[col]], "Date")) {
      x[[col]] <- ifelse(is.na(x[[col]]), NA_character_,
                         sprintf("%.17g", x[[col]]))
    }
  }
  data.table::fwrite(x, path, na = "NA", quote = "auto", dateTimeAs = "ISO",
                     logical01 = FALSE)
  invisible(path)
}
