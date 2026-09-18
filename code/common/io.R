# Input resolution and reading.

`%||%` <- function(x, y) if (is.null(x)) y else x

input_path <- function(dir, name, label = NULL) {
  path <- file.path(dir, name)
  if (!file.exists(path)) {
    stop("No input file for ", label %||% "dataset", " at: ", path,
         "\nThe restricted inputs are not distributed with this repository ",
         "and are obtained under their own approvals. Place each one under ",
         "data/ with the name above.", call. = FALSE)
  }
  path
}

read_input <- function(path) {
  read_table(path)
}
