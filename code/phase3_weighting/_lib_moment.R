options(stringsAsFactors = FALSE)

# Moment-matching weights: raking and entropy balancing.

MOMENT_METHODS <- c("raking", "entropy_balancing")

.moment_cfg <- function(key) framework_setting("phase3", "moment_matching", key)

RAKE_CONTINUOUS_BINS <- as.integer(.moment_cfg("raking_continuous_bins"))
RAKE_TOL <- as.numeric(.moment_cfg("raking_tolerance"))
ENTROPY_TOL <- as.numeric(.moment_cfg("entropy_tolerance"))
MOMENT_MAX_ITER <- as.integer(.moment_cfg("max_iterations"))
ENTROPY_INTERACT_WITH <- as.character(.moment_cfg("entropy_interaction"))

moment_weights_dir <- function(method) {
  d <- path_temp("weights", paste0("hse_", method), create = FALSE)
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
  d
}
moment_weight_path <- function(method) file.path(moment_weights_dir(method), "ukb_weights.csv")
moment_model_input_path <- function(method) file.path(moment_weights_dir(method), "model_input.csv")
moment_model_output_path <- function(method) file.path(moment_weights_dir(method), "model_output.csv")

# ---- shared checks -----------------------------------------------------------

assert_common_support_levels <- function(cohort_factor, reference_factor, variable) {
  a <- sort(unique(as.character(cohort_factor[!is.na(cohort_factor)])))
  b <- sort(unique(as.character(reference_factor[!is.na(reference_factor)])))
  only_cohort <- setdiff(a, b)
  only_reference <- setdiff(b, a)
  if (length(only_cohort) || length(only_reference)) {
    stop(sprintf(paste("Positivity failure on '%s'. Cohort-only levels: %s.",
                       "Reference-only levels: %s. Moment matching cannot weight",
                       "toward a level that is empty on one side."),
                 variable,
                 if (length(only_cohort)) paste(only_cohort, collapse = ", ") else "none",
                 if (length(only_reference)) paste(only_reference, collapse = ", ") else "none"),
         call. = FALSE)
  }
  invisible(TRUE)
}

weighted_share <- function(f, w) {
  s <- tapply(w, f, sum)
  s[is.na(s)] <- 0
  s / sum(w)
}

weighted_quantile <- function(x, w, probs) {
  o <- order(x)
  x <- x[o]; w <- w[o]
  cw <- cumsum(w) / sum(w)
  stats::approx(cw, x, xout = probs, method = "constant", f = 1,
                yleft = min(x), yright = max(x))$y
}

# ---- raking ------------------------------------------------------------------

rake_design <- function(prep, bins = RAKE_CONTINUOUS_BINS) {
  combined <- prep$combined
  info <- prep$info
  is_ukb <- combined$source == "UKB"
  ref_w <- combined$sample_weight[!is_ukb]

  out <- list()
  for (i in seq_len(nrow(info))) {
    v <- info$label[i]
    if (info$type[i] == "con") {
      x <- as.numeric(combined[[v]])
      cuts <- unique(weighted_quantile(x[!is_ukb], ref_w, seq_len(bins - 1) / bins))
      breaks <- c(-Inf, cuts, Inf)
      f <- cut(x, breaks = breaks, labels = FALSE, include.lowest = TRUE)
      f <- factor(sprintf("%s_bin%d", v, f),
                  levels = sprintf("%s_bin%d", v, seq_len(length(breaks) - 1)))
      empty <- setdiff(levels(f), as.character(f[is_ukb]))
      if (length(empty)) {
        message(sprintf(paste0(
          "[WARNING] %s: no cohort row falls in %s, so raking cannot move the ",
          "weighted share onto the reference sample's target for it."),
          v, paste(empty, collapse = ", ")))
      }
    } else {
      f <- factor(as.character(combined[[v]]))
      assert_common_support_levels(f[is_ukb], f[!is_ukb], v)
    }
    f <- factor(as.character(f))
    out[[v]] <- f
  }
  out
}

rake_weights <- function(prep) {
  tol <- RAKE_TOL
  design <- rake_design(prep)
  is_ukb <- prep$combined$source == "UKB"
  ref_w <- prep$combined$sample_weight[!is_ukb]

  targets <- lapply(design, function(f) weighted_share(f[!is_ukb], ref_w))
  cohort <- lapply(design, function(f) droplevels(f[is_ukb]))

  w <- rep(1, sum(is_ukb))
  trace <- data.frame()
  for (it in seq_len(MOMENT_MAX_ITER)) {
    for (v in names(cohort)) {
      cur <- weighted_share(cohort[[v]], w)
      ratio <- targets[[v]][names(cur)] / cur
      ratio[!is.finite(ratio)] <- 1
      w <- w * as.numeric(ratio[as.character(cohort[[v]])])
    }
    gap <- max(vapply(names(cohort), function(v) {
      max(abs(weighted_share(cohort[[v]], w) - targets[[v]][levels(cohort[[v]])]))
    }, numeric(1)))
    trace <- rbind(trace, data.frame(iteration = it, max_gap = gap))
    if (gap < tol) break
  }
  list(w = normalise_weights(w), trace = trace, iterations = nrow(trace),
       converged = gap < tol, achieved_gap = gap, constraints = names(cohort),
       dropped = character(0))
}

# ---- entropy balancing -------------------------------------------------------

entropy_constraints <- function(prep) {
  interact_with <- ENTROPY_INTERACT_WITH
  combined <- prep$combined
  info <- prep$info
  is_ukb <- combined$source == "UKB"

  cols <- list()
  for (i in seq_len(nrow(info))) {
    v <- info$label[i]
    if (info$type[i] == "con") {
      x <- as.numeric(combined[[v]])
      mu <- mean(x[is_ukb]); sdv <- stats::sd(x[is_ukb])
      if (!is.finite(sdv) || sdv <= 0) next
      z <- (x - mu) / sdv
      cols[[paste0(v, "_mean")]] <- z
      cols[[paste0(v, "_var")]] <- z^2
    } else {
      f <- factor(as.character(combined[[v]]))
      assert_common_support_levels(f[is_ukb], f[!is_ukb], v)
      lev <- levels(droplevels(f[is_ukb]))
      for (l in lev[-length(lev)]) cols[[paste0(v, "_", l)]] <- as.numeric(f == l)
    }
  }
  if (!length(cols)) stop("No usable moment conditions.", call. = FALSE)

  if (!is.null(interact_with) && interact_with %in% info$label) {
    g <- factor(as.character(combined[[interact_with]]))
    ind <- as.numeric(g == levels(droplevels(g[is_ukb]))[1])
    base <- names(cols)[!startsWith(names(cols), paste0(interact_with, "_"))]
    for (nm in base) cols[[paste0(nm, ":", interact_with)]] <- cols[[nm]] * ind
  }

  C <- do.call(cbind, cols)
  colnames(C) <- names(cols)
  list(C = C[is_ukb, , drop = FALSE],
       target = as.numeric(stats::cov.wt(C[!is_ukb, , drop = FALSE],
                                         wt = combined$sample_weight[!is_ukb],
                                         method = "ML")$center),
       names = colnames(C))
}

entropy_balance_weights <- function(prep) {
  tol <- ENTROPY_TOL
  cons <- entropy_constraints(prep)
  Z <- sweep(cons$C, 2, cons$target, "-")

  infeasible <- which(cons$target < apply(cons$C, 2, min) |
                      cons$target > apply(cons$C, 2, max))
  if (length(infeasible)) {
    stop("Targets outside the cohort's range for: ",
         paste(cons$names[infeasible], collapse = ", "),
         ". No weighting can match these moments.", call. = FALSE)
  }

  qrz <- qr(Z)
  keep <- sort(qrz$pivot[seq_len(qrz$rank)])
  dropped <- setdiff(seq_along(cons$names), keep)
  Z <- Z[, keep, drop = FALSE]
  kept_names <- cons$names[keep]

  lambda <- rep(0, ncol(Z))
  trace <- data.frame()
  for (it in seq_len(MOMENT_MAX_ITER)) {
    eta <- as.numeric(Z %*% lambda)
    eta <- eta - max(eta)
    p <- exp(-eta); p <- p / sum(p)
    grad <- -as.numeric(crossprod(Z, p))
    gap <- max(abs(grad))
    trace <- rbind(trace, data.frame(iteration = it, max_gap = gap))
    if (gap < tol) break
    H <- crossprod(Z, Z * p) - tcrossprod(as.numeric(crossprod(Z, p)))
    step <- tryCatch(solve(H + diag(1e-10, ncol(Z)), grad),
                     error = function(e) grad)
    obj <- function(l) {
      e <- as.numeric(Z %*% l); log(sum(exp(-(e - max(e))))) - max(e)
    }
    f0 <- obj(lambda); alpha <- 1; stalled <- FALSE
    repeat {
      cand <- lambda - alpha * step
      if (obj(cand) <= f0) break
      alpha <- alpha / 2
      if (alpha < 1e-12) { stalled <- TRUE; break }
    }
    if (stalled) break
    lambda <- lambda - alpha * step
  }

  eta <- as.numeric(Z %*% lambda); eta <- eta - max(eta)
  list(w = normalise_weights(exp(-eta)), trace = trace, iterations = nrow(trace),
       converged = gap < tol, achieved_gap = gap, constraints = kept_names,
       dropped = if (length(dropped)) cons$names[dropped] else character(0))
}

# ---- diagnostics -------------------------------------------------------------

moment_balance_table <- function(prep, w) {
  combined <- prep$combined
  info <- prep$info
  is_ukb <- combined$source == "UKB"
  ref_w <- combined$sample_weight[!is_ukb]
  rows <- list()

  for (i in seq_len(nrow(info))) {
    v <- info$label[i]
    if (info$type[i] == "con") {
      x <- as.numeric(combined[[v]])
      xu <- x[is_ukb]; xr <- x[!is_ukb]
      tm <- stats::weighted.mean(xr, ref_w)
      ts <- sqrt(stats::weighted.mean((xr - tm)^2, ref_w))
      wm <- stats::weighted.mean(xu, w)
      ws <- sqrt(stats::weighted.mean((xu - wm)^2, w))
      rows[[length(rows) + 1]] <- data.frame(
        variable = v, level = "mean", target = tm,
        unweighted = mean(xu), weighted = wm,
        gap_before = mean(xu) - tm, gap_after = wm - tm)
      rows[[length(rows) + 1]] <- data.frame(
        variable = v, level = "sd", target = ts,
        unweighted = stats::sd(xu), weighted = ws,
        gap_before = stats::sd(xu) - ts, gap_after = ws - ts)
    } else {
      f <- factor(as.character(combined[[v]]))
      tgt <- weighted_share(f[!is_ukb], ref_w)
      unw <- weighted_share(f[is_ukb], rep(1, sum(is_ukb)))
      wtd <- weighted_share(f[is_ukb], w)
      for (l in names(tgt)) {
        t0 <- 100 * as.numeric(tgt[l])
        u0 <- 100 * as.numeric(if (l %in% names(unw)) unw[l] else 0)
        w0 <- 100 * as.numeric(if (l %in% names(wtd)) wtd[l] else 0)
        rows[[length(rows) + 1]] <- data.frame(
          variable = v, level = l, target = t0,
          unweighted = u0, weighted = w0,
          gap_before = u0 - t0, gap_after = w0 - t0)
      }
    }
  }
  do.call(rbind, rows)
}
