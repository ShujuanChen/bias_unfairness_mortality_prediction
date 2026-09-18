
# Reader for the analytic population.

read_population_membership <- function(root = .repo_root()) {
  p <- path_temp("population", "ukb_population_membership.csv", root = root, create = FALSE)
  if (!file.exists(p)) {
    stop("Missing the analytic population at ", p,
         ".\nBuild it first with phase 1 stage 5.", call. = FALSE)
  }
  read_table(p, stringsAsFactors = FALSE)
}

population_pipelines <- function(root = .repo_root()) {
  p <- path_temp("population", "pipelines.txt", root = root, create = FALSE)
  if (!file.exists(p)) {
    stop("No population/pipelines.txt under temp/. Run phase 1 stage 5 first.",
         call. = FALSE)
  }
  readLines(p)
}

require_pipeline <- function(name, what, root = .repo_root()) {
  built <- population_pipelines(root)
  if (!name %in% built) {
    stop(what, " needs the ", name, " pipeline, but this run's population was ",
         "built for ", paste(built, collapse = ", "),
         ". Rebuild phase 1 with --pipelines including ", name, ".", call. = FALSE)
  }
  invisible(TRUE)
}

assert_covers_analytic_population <- function(ids, label = "weights",
                                              root = .repo_root()) {
  keep <- analytic_ids(root = root)
  ids <- unique(as.integer(ids))
  absent <- setdiff(keep, ids)
  extra <- setdiff(ids, keep)

  if (length(absent) || length(extra)) {
    stop(sprintf(
      paste0("%s cover %s participants against an analytic population of %s: ",
             "%s absent, %s unexpected.%s%s\nThe weights and the population are ",
             "out of step. Re-estimate them on the population fixed in phase 1."),
      label, format(length(ids), big.mark = ","), format(length(keep), big.mark = ","),
      format(length(absent), big.mark = ","), format(length(extra), big.mark = ","),
      if (length(absent)) paste0("\n  first absent: ",
                                 paste(utils::head(absent, 5), collapse = ", ")) else "",
      if (length(extra)) paste0("\n  first unexpected: ",
                                paste(utils::head(extra, 5), collapse = ", ")) else ""),
      call. = FALSE)
  }

  message(sprintf("[population] %s cover the analytic population exactly (%s participants)",
                  label, format(length(keep), big.mark = ",")))
  invisible(TRUE)
}

analytic_ids <- function(root = .repo_root()) {
  m <- read_population_membership(root)
  as.integer(m$participant_id[m$in_unified])
}

.read_id_file <- function(file, column, stage_hint, root = .repo_root()) {
  p <- path_temp("population", file, root = root, create = FALSE)
  if (!file.exists(p)) {
    stop("Missing ", file, " at ", p, ".\nBuild it first with ", stage_hint,
         ".", call. = FALSE)
  }
  read_table(p, stringsAsFactors = FALSE)[[column]]
}

pmr_analytic_ids <- function(root = .repo_root()) {
  as.integer(.read_id_file("pmr_analytic_ids.csv", "pmr_id",
                           "phase 1 stage 5", root))
}

reference_ids <- function(sample, root = .repo_root()) {
  spec <- switch(sample,
                 hse    = list(file = "hse_reference_ids.csv", column = "ref_id"),
                 census = list(file = "census_reference_ids.csv", column = "ref_id"),
                 stop("Unknown reference sample: ", sample, call. = FALSE))
  as.integer(.read_id_file(spec$file, spec$column, "phase 1 stage 5", root))
}
