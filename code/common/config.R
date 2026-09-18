# Reader for framework_config.json

.config_own_dir <- local({
  for (i in rev(seq_len(sys.nframe()))) {
    ofile <- sys.frame(i)$ofile
    if (!is.null(ofile)) return(dirname(normalizePath(ofile)))
  }
  arg <- grep("--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(arg)) dirname(normalizePath(sub("--file=", "", arg[1]))) else getwd()
})

if (!exists("read_table", mode = "function")) source(file.path(.config_own_dir, "tabular.R"))
if (!exists(".repo_root", mode = "function")) source(file.path(.config_own_dir, "paths.R"))
if (!exists("input_path", mode = "function")) source(file.path(.config_own_dir, "io.R"))
if (!exists("analytic_ids", mode = "function")) source(file.path(.config_own_dir, "population.R"))
if (!exists("assert_expected_codes", mode = "function")) source(file.path(.config_own_dir, "recode.R"))

get_framework_root <- function(start_path = getwd()) .repo_root(start_path)

read_framework_config <- function(start_path = getwd()) {
  if (!requireNamespace("jsonlite", quietly = TRUE)) {
    stop("Package 'jsonlite' is required to read framework_config.json.", call. = FALSE)
  }

  root <- get_framework_root(start_path)
  cfg_path <- file.path(root, "framework_config.json")

  cfg <- jsonlite::fromJSON(cfg_path, simplifyVector = FALSE)
  cfg$framework_root <- root
  cfg$framework_config_path <- cfg_path
  cfg
}

.framework_settings <- new.env(parent = emptyenv())

framework_setting <- function(..., simplify = TRUE) {
  keys <- c(...)
  if (is.null(.framework_settings$cfg)) {
    .framework_settings$cfg <- read_framework_config()
  }
  value <- cfg_get(.framework_settings$cfg, keys)
  if (simplify && is.list(value) && !length(names(value))) {
    return(unlist(value, use.names = FALSE))
  }
  value
}

cfg_get <- function(cfg, keys) {
  value <- cfg
  for (key in keys) {
    if (is.null(value[[key]])) {
      stop("Missing config entry: ", paste(keys, collapse = "."), call. = FALSE)
    }
    value <- value[[key]]
  }
  value
}
