options(stringsAsFactors = FALSE)

normalise_weights <- function(x) x / mean(x, na.rm = TRUE)

effective_sample_size <- function(w) {
  w <- w[is.finite(w) & w > 0]
  if (!length(w)) return(NA_real_)
  sum(w)^2 / sum(w^2)
}

winsorise_weights <- function(w,
                              probs = framework_setting("phase3", "winsorise_percentiles")) {
  stopifnot(length(probs) == 2L, probs[1] >= 0, probs[2] <= 1, probs[1] < probs[2])
  qs <- stats::quantile(w, probs = probs, na.rm = TRUE, names = FALSE)
  pmin(pmax(w, qs[1]), qs[2])
}

pending_winsorisation_path <- function(dir) file.path(dir, "winsorise.txt")

write_pending_winsorisation <- function(dir, apply = TRUE) {
  writeLines(if (isTRUE(apply)) "1" else "0", pending_winsorisation_path(dir))
}

winsorisation_pending <- function(dir) {
  p <- pending_winsorisation_path(dir)
  if (!file.exists(p)) return(FALSE)
  identical(trimws(readLines(p, warn = FALSE))[1], "1")
}

apply_pending_winsorisation <- function(dir, w) {
  if (!winsorisation_pending(dir)) return(w)
  normalise_weights(winsorise_weights(w))
}
