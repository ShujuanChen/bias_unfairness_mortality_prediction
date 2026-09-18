#!/usr/bin/env Rscript
#
# Diagnostics for the set of participation weights this run reports. 

.dir <- dirname(normalizePath(sub("^--file=", "",
  grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1])))
source(file.path(.dir, "..", "common", "config.R"))
source(file.path(.dir, "_lib_weights.R"))

cfg <- read_framework_config(.dir)
sources <- as.character(cfg$phase3$weight_method)

# ── Metrics ──────────────────────────────────────────────────────────────────

auc_weighted <- function(score, y, w) {
  ok <- is.finite(score) & is.finite(y) & is.finite(w) & w > 0
  score <- score[ok]; y <- y[ok]; w <- w[ok]
  o <- order(score)
  score <- score[o]; y <- y[o]; w <- w[o]
  grp <- cumsum(c(TRUE, diff(score) != 0))
  w1 <- as.numeric(tapply(w * y, grp, sum))
  w0 <- as.numeric(tapply(w * (1 - y), grp, sum))
  w1[is.na(w1)] <- 0; w0[is.na(w0)] <- 0
  below <- c(0, head(cumsum(w0), -1))
  if (sum(w1) <= 0 || sum(w0) <= 0) return(NA_real_)
  (sum(w1 * below) + 0.5 * sum(w1 * w0)) / (sum(w1) * sum(w0))
}

standardised_difference <- function(x, is_cohort, w) {
  levels_x <- sort(unique(as.character(x[!is.na(x)])))
  numeric_x <- suppressWarnings(!anyNA(as.numeric(x[!is.na(x)]))) && length(levels_x) > 10L

  one <- function(indicator) {
    a <- indicator[is_cohort]; wa <- w[is_cohort]
    b <- indicator[!is_cohort]; wb <- w[!is_cohort]
    ma <- sum(wa * a) / sum(wa)
    mb <- sum(wb * b) / sum(wb)
    va <- sum(wa * (a - ma)^2) / sum(wa)
    vb <- sum(wb * (b - mb)^2) / sum(wb)
    pooled <- sqrt((va + vb) / 2)
    if (!is.finite(pooled) || pooled == 0) return(c(NA_real_, ma, mb))
    c((ma - mb) / pooled, ma, mb)
  }

  if (numeric_x) {
    v <- one(as.numeric(x))
    return(data.frame(level = "", smd = v[1], cohort = v[2], reference = v[3]))
  }
  do.call(rbind, lapply(levels_x, function(lev) {
    v <- one(as.numeric(as.character(x) == lev))
    data.frame(level = lev, smd = v[1], cohort = v[2], reference = v[3])
  }))
}

support_overlap <- function(p, is_cohort, sw) {
  pc <- p[is_cohort]; pr <- p[!is_cohort]; wr <- sw[!is_cohort]
  lower <- max(min(pc), min(pr))
  upper <- min(max(pc), max(pr))
  outside_r <- pr < lower | pr > upper
  list(lower = lower, upper = upper,
       share_cohort_outside = mean(pc < lower | pc > upper),
       share_reference_outside = sum(wr[outside_r]) / sum(wr))
}

probability_density <- function(p, is_cohort, sw, breaks = seq(0, 1, by = 0.01)) {
  mid <- utils::head(breaks, -1) + diff(breaks) / 2
  one <- function(keep) {
    binned <- cut(p[keep], breaks = breaks, include.lowest = TRUE)
    mass <- tapply(sw[keep], binned, sum)
    mass[is.na(mass)] <- 0
    as.numeric(mass) / (sum(sw[keep]) * diff(breaks))
  }
  rbind(
    data.frame(sample = "cohort", probability = mid, density = one(is_cohort),
               stringsAsFactors = FALSE),
    data.frame(sample = "reference", probability = mid, density = one(!is_cohort),
               stringsAsFactors = FALSE)
  )
}

# ── Per weight source ────────────────────────────────────────────────────────

summary_rows <- list()
balance_rows <- list()
positivity_rows <- list()
density_rows <- list()

for (key in sources) {
  dir_key <- path_temp("weights", key, create = FALSE)
  out_path <- file.path(dir_key, "model_output.csv")
  in_path <- file.path(dir_key, "model_input.csv")
  if (!file.exists(out_path) || !file.exists(in_path)) {
    message("no weights for ", key, ", skipped")
    next
  }

  fitted <- read_table(out_path, stringsAsFactors = FALSE)
  design <- read_table(in_path, stringsAsFactors = FALSE)
  fitted$is_cohort <- fitted$source == "UKB"
  fitted_model <- any(is.finite(fitted$prob_ukb))
  w_raw <- fitted$w[fitted$is_cohort]
  fitted$w[fitted$is_cohort] <- apply_pending_winsorisation(
    dir_key, fitted$w[fitted$is_cohort])

  y <- as.integer(fitted$is_cohort)
  p <- fitted$prob_ukb
  sw <- fitted$sample_weight
  brier <- if (fitted_model) sum(sw * (p - y)^2) / sum(sw) else NA_real_

  w_cohort <- fitted$w[fitted$is_cohort]
  n_cohort <- length(w_cohort)
  ess <- effective_sample_size(w_cohort)
  q <- stats::quantile(w_cohort, c(0, 0.005, 0.01, 0.25, 0.5, 0.75, 0.99, 0.995, 1),
                       na.rm = TRUE)
  qr <- stats::quantile(w_raw, c(0.005, 0.01, 0.5, 0.99, 0.995, 1),
                        na.rm = TRUE)
  top <- sort(w_cohort, decreasing = TRUE)
  share_top_1pct <- sum(head(top, ceiling(0.01 * n_cohort))) / sum(w_cohort)

  summary_rows[[key]] <- data.frame(
    weight_source = key,
    n_cohort = n_cohort,
    n_reference = sum(!fitted$is_cohort),
    auc = if (fitted_model) auc_weighted(p, y, sw) else NA_real_,
    brier = brier,
    mean_prob_cohort = if (fitted_model) mean(p[fitted$is_cohort]) else NA_real_,
    mean_prob_reference = if (fitted_model) mean(p[!fitted$is_cohort]) else NA_real_,
    weight_min = q[[1]], weight_p0_5 = q[[2]], weight_p1 = q[[3]],
    weight_p25 = q[[4]], weight_median = q[[5]], weight_p75 = q[[6]],
    weight_p99 = q[[7]], weight_p99_5 = q[[8]], weight_max = q[[9]],
    weight_mean = mean(w_cohort), weight_sd = stats::sd(w_cohort),
    weight_cv = stats::sd(w_cohort) / mean(w_cohort),
    raw_p0_5 = qr[[1]], raw_p1 = qr[[2]], raw_median = qr[[3]],
    raw_p99 = qr[[4]], raw_p99_5 = qr[[5]], raw_max = qr[[6]],
    raw_cv = stats::sd(w_raw) / mean(w_raw),
    effective_sample_size = ess,
    effective_fraction = ess / n_cohort,
    share_of_weight_on_top_1pct = share_top_1pct,
    stringsAsFactors = FALSE
  )

  if (fitted_model) {
  overlap <- support_overlap(p, fitted$is_cohort, sw)
  p_cohort <- p[fitted$is_cohort]
  odds <- (1 - p_cohort) / p_cohort
  positivity_rows[[key]] <- data.frame(
    weight_source = key,
    prob_min_cohort = min(p_cohort), prob_max_cohort = max(p_cohort),
    prob_min_reference = min(p[!fitted$is_cohort]),
    prob_max_reference = max(p[!fitted$is_cohort]),
    overlap_lower = overlap$lower, overlap_upper = overlap$upper,
    share_cohort_outside_overlap = overlap$share_cohort_outside,
    share_reference_outside_overlap = overlap$share_reference_outside,
    inverse_odds_min = min(odds), inverse_odds_max = max(odds),
    stringsAsFactors = FALSE
  )

  density_rows[[key]] <- cbind(
    weight_source = key,
    probability_density(p, fitted$is_cohort, sw),
    stringsAsFactors = FALSE
  )
  }

  population <- analytic_ids()
  weighted_ids <- as.integer(fitted$eid[fitted$is_cohort])
  if (!setequal(population, weighted_ids)) {
    stop(sprintf(paste0("%s: the weight file and the analytic population differ. ",
                        "%d of %d participants have no weight and %d weighted ",
                        "records are outside the population."),
                 key, length(setdiff(population, weighted_ids)), length(population),
                 length(setdiff(weighted_ids, population))), call. = FALSE)
  }

  merged <- merge(design, fitted[, c("eid", "w", "is_cohort")], by = "eid")
  if (nrow(merged) != nrow(design)) {
    stop(sprintf("%s: %d of %d model-input rows have no fitted weight.",
                 key, nrow(design) - nrow(merged), nrow(design)), call. = FALSE)
  }
  aux <- setdiff(names(design), c("eid", "source", "sample_weight"))
  unweighted <- ifelse(merged$is_cohort, 1, merged$sample_weight)

  for (v in aux) {
    before <- standardised_difference(merged[[v]], merged$is_cohort, unweighted)
    after <- standardised_difference(merged[[v]], merged$is_cohort, merged$w)
    balance_rows[[paste(key, v)]] <- data.frame(
      weight_source = key, variable = v, level = before$level,
      cohort_unweighted = before$cohort, cohort_weighted = after$cohort,
      reference = before$reference,
      smd_before = before$smd, smd_after = after$smd,
      stringsAsFactors = FALSE
    )
  }
  message("diagnostics computed for ", key)
}

if (!length(summary_rows)) {
  stop("No weights for ", paste(sources, collapse = ", "),
       ". Run the phase 3 stage that estimates them.", call. = FALSE)
}

summary_out <- do.call(rbind, summary_rows)
balance_out <- do.call(rbind, balance_rows)

outputs <- list(
  weight_summary = summary_out,
  weight_balance = balance_out
)
if (length(positivity_rows)) {
  outputs$weight_positivity <- do.call(rbind, positivity_rows)
  outputs$weight_probability_density <- do.call(rbind, density_rows)
}

FIGURE_INPUTS <- c("weight_probability_density", "weight_balance")
for (name in names(outputs)) {
  target <- if (name %in% FIGURE_INPUTS) {
    path_temp(paste0(name, ".csv"))
  } else {
    path_results("phase3_weighting", paste0(name, ".xlsx"))
  }
  write_table(outputs[[name]], target)
  message("wrote ", target)
}

worst <- balance_out[order(-abs(balance_out$smd_after)), ][1:min(10, nrow(balance_out)), ]
message("\nlargest remaining differences after weighting:")
print(worst[, c("weight_source", "variable", "level", "smd_before", "smd_after")],
      row.names = FALSE)
