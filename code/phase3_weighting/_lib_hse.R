options(stringsAsFactors = FALSE)

# Core of the HSE participation-weighting pipeline.

.resolve_own_dir <- function() {
  for (i in rev(seq_len(sys.nframe()))) {
    ofile <- sys.frame(i)$ofile
    if (!is.null(ofile)) return(dirname(normalizePath(ofile)))
  }
  arg <- grep("--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(arg)) return(dirname(normalizePath(sub("--file=", "", arg[1]))))
  getwd()
}

.hse_lib_dir <- tryCatch(.resolve_own_dir(), error = function(e) getwd())

source(file.path(.hse_lib_dir, "..", "common", "config.R"))
source(file.path(.hse_lib_dir, "..", "phase1_data", "_lib_shared.R"))
source(file.path(.hse_lib_dir, "..", "phase1_data", "_lib_hse.R"))
source(file.path(.hse_lib_dir, "_lib_superlearner.R"))
source(file.path(.hse_lib_dir, "_lib_weights.R"))

# ---- result-tree paths -------------------------------------------------------

hse_weights_dir <- function(model) {
  d <- path_temp("weights", paste0("hse_", model), create = FALSE)
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
  d
}

hse_model_input_path <- function(model) file.path(hse_weights_dir(model), "model_input.csv")
hse_weight_path <- function(model) file.path(hse_weights_dir(model), "ukb_weights.csv")
hse_model_output_path <- function(model) file.path(hse_weights_dir(model), "model_output.csv")

# ---- stack HSE + UKB ---------------------------------------------------------

read_harmonised <- function(name, stage) {
  p <- path_temp("harmonised", name, create = FALSE)
  if (!file.exists(p)) {
    stop("Missing ", basename(p), ". Run phase 1 stage ", stage, " first.", call. = FALSE)
  }
  read_table(p, stringsAsFactors = FALSE)
}

prepare_training_data <- function() {
  hse_all <- read_harmonised("hse.csv", 3)
  ukb_cohort <- read_harmonised("ukb_with_hse.csv", 3)

  hse_cohort <- hse_all[hse_all$ref_id %in% reference_ids("hse"), , drop = FALSE]
  message(sprintf("[reference cohort] HSE %d of %d rows retained",
                  nrow(hse_cohort), nrow(hse_all)))

  usable <- is.finite(hse_cohort$weight_individual) &
    hse_cohort$weight_individual > 0
  if (any(!usable)) {
    message(sprintf("[reference cohort] HSE %d rows dropped for a missing or ",
                    sum(!usable)), "non-positive interview weight")
    hse_cohort <- hse_cohort[usable, , drop = FALSE]
  }

  aux_path <- path_temp("harmonised", "hse_auxiliaries.txt", create = FALSE)
  if (file.exists(aux_path)) {
    declared <- readLines(aux_path)
    missing_aux <- setdiff(declared, names(ukb_cohort))
    if (length(missing_aux)) {
      stop("The harmonised UKB frame is missing declared auxiliaries: ",
           paste(missing_aux, collapse = ", "), call. = FALSE)
    }
  }

  ukb_obj <- list(final_cohort = ukb_cohort)
  hse_obj <- list(final_cohort = hse_cohort)

  keep <- analytic_ids()
  before <- nrow(ukb_obj$final_cohort)
  ukb_obj$final_cohort <- ukb_obj$final_cohort[ukb_obj$final_cohort$eid %in% keep, , drop = FALSE]
  message(sprintf("[analytic population] UKB %d of %d rows retained",
                  nrow(ukb_obj$final_cohort), before))
  if (nrow(ukb_obj$final_cohort) != length(keep)) {
    stop(sprintf(paste("The analytic population holds %d participants but only %d",
                       "are present in the HSE-harmonised frame. Rerun phase 1."),
                 length(keep), nrow(ukb_obj$final_cohort)), call. = FALSE)
  }

  info <- hse_predictor_info()
  labels <- info$label

  hse_obj$final_cohort$eid <- max(ukb_obj$final_cohort$eid, na.rm = TRUE) +
    seq_len(nrow(hse_obj$final_cohort))

  hse_model <- hse_obj$final_cohort[, c("eid", "weight_individual", labels), drop = FALSE]
  hse_model$source <- "HSE"
  names(hse_model)[names(hse_model) == "weight_individual"] <- "sample_weight"
  hse_model$sample_weight <- normalise_weights(hse_model$sample_weight)

  ukb_model <- ukb_obj$final_cohort[, c("eid", labels), drop = FALSE]
  ukb_model$source <- "UKB"
  ukb_model$sample_weight <- 1

  boot_dir <- trimws(Sys.getenv("REWEIGHTING_BOOTSTRAP_DIR", ""))
  if (nzchar(boot_dir)) {
    mult_path <- file.path(boot_dir, "ukb_multiplicity.csv")
    if (!file.exists(mult_path)) {
      stop("REWEIGHTING_BOOTSTRAP_DIR is set to ", boot_dir, " but ",
           basename(mult_path), " is not there. Draw the replicate before ",
           "estimating its participation model.", call. = FALSE)
    }
    multipliers <- read_table(mult_path)
    idx <- match(ukb_model$eid, multipliers$participant_id)
    if (anyNA(idx)) {
      stop(sum(is.na(idx)), " of ", nrow(ukb_model), " cohort rows have no ",
           "multiplier in ", mult_path, ". The replicate was drawn against a ",
           "different population.", call. = FALSE)
    }
    ukb_model$sample_weight <- as.numeric(multipliers$multiplicity[idx])

    if (any(!is.finite(ukb_model$sample_weight)) ||
        any(ukb_model$sample_weight <= 0)) {
      stop(mult_path, " holds a multiplier that is not strictly positive.",
           call. = FALSE)
    }
    message(sprintf(paste0("[bootstrap] cohort multipliers applied to the ",
                           "participation model: %s rows, largest %s"),
                    format(nrow(ukb_model), big.mark = ","),
                    format(max(ukb_model$sample_weight))))
  }

  stacked <- data.table::rbindlist(list(hse_model, ukb_model), fill = TRUE)
  complete <- complete.cases(stacked[, ..labels])
  combined <- as.data.frame(stacked[complete, , drop = FALSE])

  for (label in info$label[info$type != "con"]) {
    values <- as.character(combined[[label]])
    present <- unique(values[!is.na(values)])
    numeric_codes <- suppressWarnings(as.numeric(present))
    ordered <- if (!anyNA(numeric_codes)) {
      present[order(numeric_codes)]
    } else {
      sort(present)
    }
    combined[[label]] <- factor(values, levels = ordered)
  }

  active <- vapply(info$label, function(label) {
    vals <- combined[[label]]
    length(unique(as.character(vals[!is.na(vals)]))) > 1
  }, logical(1))
  info_active <- info[active, , drop = FALSE]
  if (!nrow(info_active)) stop("No varying predictors remain after the HSE/UKB restrictions.")

  combined$sample <- as.integer(combined$source == "UKB")

  for (i in seq_len(nrow(info_active))) {
    v <- info_active$label[i]
    combined[[v]] <- if (identical(info_active$type[i], "con")) {
      as.numeric(combined[[v]])
    } else {
      factor(as.character(combined[[v]]))
    }
  }

  list(
    hse_obj = hse_obj,
    ukb_obj = ukb_obj,
    prediction_info = info_active,
    prediction_labels = info_active$label,
    combined = combined,
    model_input = combined[, c("eid", "source", "sample_weight", info_active$label), drop = FALSE]
  )
}

# ---- design matrix for the lasso (one-hot with pairwise interactions) --------

build_dummy_design <- function(df, info) {
  cat_labels <- info$label[info$type == "cat"]
  bin_labels <- info$label[info$type == "bin"]
  con_labels <- info$label[info$type == "con"]

  bin_part <- if (length(bin_labels)) {
    fastDummies::dummy_cols(subset(df, select = bin_labels), select_columns = bin_labels,
                            remove_first_dummy = TRUE, ignore_na = TRUE)
  } else data.frame()
  cat_part <- if (length(cat_labels)) {
    fastDummies::dummy_cols(subset(df, select = cat_labels), select_columns = cat_labels,
                            remove_first_dummy = FALSE, ignore_na = TRUE)
  } else data.frame()

  expanded <- c(colnames(cat_part), colnames(bin_part))
  con_part <- subset(df, select = c("eid", con_labels))
  design <- cbind(con_part, cat_part, bin_part)
  keep <- c(con_labels, unique(expanded[!expanded %in% c(cat_labels, bin_labels)]))
  list(data = design, vars = keep)
}

# ---- fits --------------------------------------------------------------------

assemble_inverse_odds_weights <- function(combined, prob_ukb) {
  if (any(!is.finite(prob_ukb) | prob_ukb <= 0 | prob_ukb >= 1, na.rm = TRUE)) {
    stop("Participation probabilities produced invalid inverse-odds weights.")
  }
  is_ukb <- combined$source == "UKB"
  combined$prob_ukb <- prob_ukb
  combined$inverse_odds <- (1 - prob_ukb) / prob_ukb
  combined$w <- NA_real_
  combined$w[is_ukb] <- normalise_weights(combined$inverse_odds[is_ukb])
  combined$w[!is_ukb] <- combined$sample_weight[!is_ukb]
  combined
}

PARTICIPATION_CROSSFIT_SEED <- as.integer(
  framework_setting("phase3", "crossfit_seed"))

participation_outer_folds <- function(y) {
  set.seed(PARTICIPATION_CROSSFIT_SEED)
  SuperLearner::CVFolds(
    N = length(y), id = NULL, Y = y,
    cvControl = list(V = CROSSFIT_V, stratifyCV = TRUE, shuffle = TRUE,
                     validRows = NULL)
  )
}

.LASSO <- framework_setting("phase3", "hse_lassologit", simplify = FALSE)

fit_lasso_weights <- function() {
  prepared <- prepare_training_data()
  combined <- prepared$combined
  design <- build_dummy_design(combined, prepared$prediction_info)
  dummy_data <- data.frame(design$data[, design$vars, drop = FALSE], check.names = TRUE)

  if (!identical(.LASSO$design, "second_order_interactions")) {
    stop("Unknown lasso design '", .LASSO$design, "'. This stage fits ",
         "second_order_interactions.", call. = FALSE)
  }
  dummy_data <- participation_design(
    dummy_data,
    prepared$prediction_info$label[prepared$prediction_info$type == "con"])
  x <- Matrix::sparse.model.matrix(stats::as.formula("~ .*."), data = dummy_data)
  y <- as.integer(combined$source == "UKB")

  folds <- participation_outer_folds(y)

  prob_ukb <- rep(NA_real_, length(y))
  for (k in seq_along(folds)) {
    held <- folds[[k]]
    train <- setdiff(seq_along(y), held)

    path <- glmnet::glmnet(
      x = x[train, , drop = FALSE], y = y[train], family = .LASSO$family,
      alpha = as.numeric(.LASSO$alpha),
      weights = combined$sample_weight[train],
      nlambda = as.integer(.LASSO$nlambda), standardize = FALSE
    )$lambda
    cvfit <- glmnet::cv.glmnet(
      x = x[train, , drop = FALSE], y = y[train], family = .LASSO$family,
      alpha = as.numeric(.LASSO$alpha), type.measure = .LASSO$type_measure,
      nfolds = as.integer(.LASSO$nfolds),
      weights = combined$sample_weight[train],
      lambda = path, parallel = FALSE, standardize = FALSE
    )
    if (length(unique(cvfit$glmnet.fit$lambda)) < 2L) {
      stop(sprintf(paste0(
        "The lasso path collapsed to %d distinct penalty value(s) on outer ",
        "fold %d. There is nothing to select over, so no probability of cohort ",
        "membership can be predicted. The design has %d columns over %d ",
        "training rows: check that it does not separate the two sources."),
        length(unique(cvfit$glmnet.fit$lambda)), k, ncol(x), length(train)),
        call. = FALSE)
    }
    prob_ukb[held] <- as.numeric(stats::predict(
      cvfit, newx = x[held, , drop = FALSE], s = .LASSO$lambda_rule,
      type = "response")[, 1])
  }
  if (anyNA(prob_ukb)) stop("The outer folds did not cover every row.")

  combined <- assemble_inverse_odds_weights(combined, prob_ukb)
  combined$crossfit_fold <- crossfit_fold_index(list(folds = folds), length(y))

  ukb_rows <- combined$source == "UKB"
  prepared$ukb_weights <- combined[ukb_rows, c("eid", "w", "prob_ukb"), drop = FALSE]
  prepared$combined <- combined
  prepared
}

fit_superlearner_weights <- function() {
  prepared <- prepare_training_data()
  combined <- prepared$combined
  vars <- prepared$prediction_labels

  stopifnot(vapply(combined[vars], function(x) is.numeric(x) || is.factor(x),
                   logical(1)))

  x_mm <- participation_matrix(combined, vars)
  if (!ncol(x_mm)) stop("The HSE/UKB design matrix has no columns.")

  x <- participation_design(
    as.data.frame(x_mm),
    prepared$prediction_info$label[prepared$prediction_info$type == "con"])

  crossfit <- crossfit_superlearner(
    y = combined$sample, x = x,
    obs_weights = combined$sample_weight,
    outer_folds = participation_outer_folds(combined$sample)
  )

  combined <- assemble_inverse_odds_weights(combined, crossfit$prob_ukb)
  combined$crossfit_fold <- crossfit_fold_index(crossfit$fit, nrow(combined))

  ukb_rows <- combined$source == "UKB"
  prepared$ukb_weights <- combined[ukb_rows, c("eid", "w", "prob_ukb"), drop = FALSE]
  prepared$combined <- combined
  prepared
}
