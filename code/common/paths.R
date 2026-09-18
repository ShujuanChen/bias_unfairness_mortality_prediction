# Path resolution for data/, temp/ and results/. Mirrors code/common/paths.py.

.repo_root <- function(start = getwd()) {
  p <- normalizePath(start, mustWork = TRUE)
  if (!dir.exists(p)) p <- dirname(p)
  repeat {
    if (file.exists(file.path(p, "framework_config.json"))) return(p)
    parent <- dirname(p)
    if (identical(parent, p)) stop("framework_config.json not found above ", start, call. = FALSE)
    p <- parent
  }
}

run_id <- function() {
  id <- Sys.getenv("REWEIGHTING_RUN_ID", "")
  if (!nzchar(id)) {
    stop("No run identifier. Stage scripts are launched by run_pipeline.R, ",
         "which stamps each invocation with a timestamp. To run a stage alone, ",
         "set REWEIGHTING_RUN_ID to the run you want to write into.",
         call. = FALSE)
  }
  id
}

temp_root <- function(root = .repo_root()) {
  override <- Sys.getenv("REWEIGHTING_TEMP_ROOT", "")
  if (nzchar(override)) override else file.path(root, "temp", run_id())
}

path_data <- function(..., root = .repo_root()) {
  file.path(root, "data", ...)
}

path_temp <- function(..., root = .repo_root(), create = TRUE) {
  p <- file.path(temp_root(root), ...)
  if (create) dir.create(dirname(p), recursive = TRUE, showWarnings = FALSE)
  p
}

path_results <- function(..., root = .repo_root(), create = TRUE) {
  p <- file.path(root, "results", run_id(), ...)
  if (create) dir.create(dirname(p), recursive = TRUE, showWarnings = FALSE)
  p
}

