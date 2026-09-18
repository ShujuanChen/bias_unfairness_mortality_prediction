options(stringsAsFactors = FALSE)

suppressPackageStartupMessages(library(SuperLearner))

# Learner library and cross-fitting regime shared by both participation models.

.hp <- function(learner, key) {
  v <- framework_setting("phase3", "learners", learner, key)
  if (is.numeric(v)) as.numeric(v) else v
}

CROSSFIT_V <- as.integer(framework_setting("phase3", "crossfit_folds"))

.mtry <- function(p, rule) {
  m <- floor(eval(parse(text = rule), list(p = p)))
  if (!is.finite(m) || m < 1) {
    stop("phase3.learners.SL.ranger_participation.mtry_rule is '", rule,
         "', which gives ", m, " on ", p, " columns.", call. = FALSE)
  }
  min(p, as.integer(m))
}


# ---- resources ---------------------------------------------------------------

.env_count <- function(name, default) {
  n <- suppressWarnings(as.integer(Sys.getenv(name, as.character(default))))
  if (is.na(n) || n < 1L) default else n
}

crossfit_fold_workers <- function() {
  min(.env_count("REWEIGHTING_SL_FOLD_WORKERS", 1L), CROSSFIT_V)
}

fork_is_safe <- function(timeout_seconds = 10) {
  ok <- tryCatch({
    job <- parallel::mcparallel({
      x <- matrix(stats::rnorm(400L), 100L, 4L)
      y <- stats::rbinom(100L, 1L, 0.5)
      xgboost::xgb.train(
        params = list(objective = "binary:logistic", nthread = 1L),
        data = xgboost::xgb.DMatrix(x, label = y), nrounds = 2L, verbose = 0)
      TRUE
    })
    res <- parallel::mccollect(job, wait = FALSE, timeout = timeout_seconds)
    tools::pskill(job$pid)
    parallel::mccollect(job, wait = FALSE, timeout = 1)
    isTRUE(res[[1]])
  }, error = function(e) FALSE)
  isTRUE(ok)
}

learner_threads <- function() .env_count("REWEIGHTING_SL_THREADS", 1L)

# ---- progress ----------------------------------------------------------------

.crossfit_progress <- new.env(parent = emptyenv())
.crossfit_progress$fits <- 0L

announce_fit <- function(label, n_rows, started) {
  .crossfit_progress$fits <- .crossfit_progress$fits + 1L
  message(sprintf("[worker %d] fit %3d  %-14s n = %9s  %7.1fs",
                  Sys.getpid(), .crossfit_progress$fits, label,
                  format(n_rows, big.mark = ","),
                  as.numeric(difftime(Sys.time(), started, units = "secs"))))
}

# ---- learners ----------------------------------------------------------------

SL.xgboost_participation <- function(Y, X, newX, family, obsWeights, ...) {
  started <- Sys.time()
  dtr <- xgboost::xgb.DMatrix(as.matrix(X), label = Y, weight = obsWeights)
  fit <- xgboost::xgb.train(
    params = list(objective = .hp("SL.xgboost_participation", "objective"),
                  eval_metric = .hp("SL.xgboost_participation", "eval_metric"),
                  max_depth = as.integer(.hp("SL.xgboost_participation", "max_depth")),
                  eta = .hp("SL.xgboost_participation", "eta"),
                  min_child_weight = .hp("SL.xgboost_participation", "min_child_weight"),
                  subsample = .hp("SL.xgboost_participation", "subsample"),
                  colsample_bytree = .hp("SL.xgboost_participation", "colsample_bytree"),
                  gamma = .hp("SL.xgboost_participation", "gamma"),
                  lambda = .hp("SL.xgboost_participation", "lambda"),
                  nthread = learner_threads()),
    data = dtr, nrounds = as.integer(.hp("SL.xgboost_participation", "nrounds")), verbose = 0)
  object <- list(object = fit)
  class(object) <- "SL.xgboost_participation"
  pred <- predict(object, newdata = newX)
  announce_fit("SL.xgboost_participation", nrow(X), started)
  list(pred = pred, fit = object)
}

predict.SL.xgboost_participation <- function(object, newdata, ...) {
  as.numeric(stats::predict(object$object, xgboost::xgb.DMatrix(as.matrix(newdata))))
}

SL.ranger_participation <- function(Y, X, newX, family, obsWeights, ...) {
  started <- Sys.time()
  p <- ncol(X)
  fit <- ranger::ranger(
    x = X, y = factor(Y),
    num.trees = as.integer(.hp("SL.ranger_participation", "num_trees")),
    mtry = .mtry(p, .hp("SL.ranger_participation", "mtry_rule")),
    min.node.size = as.integer(.hp("SL.ranger_participation", "min_node_size")),
    sample.fraction = .hp("SL.ranger_participation", "sample_fraction"),
    replace = .hp("SL.ranger_participation", "replace"),
    splitrule = .hp("SL.ranger_participation", "splitrule"),
    probability = TRUE, case.weights = obsWeights,
    num.threads = learner_threads(), verbose = FALSE)
  object <- list(object = fit)
  class(object) <- "SL.ranger_participation"
  pred <- predict(object, newdata = newX)
  announce_fit("SL.ranger_participation", nrow(X), started)
  list(pred = pred, fit = object)
}

predict.SL.ranger_participation <- function(object, newdata, ...) {
  stats::predict(object$object, data = newdata,
                 num.threads = learner_threads())$predictions[, "1"]
}

SL.glmnet_participation <- function(Y, X, newX, family, obsWeights, ...) {
  started <- Sys.time()
  fit <- glmnet::cv.glmnet(as.matrix(X), Y, weights = obsWeights,
                           family = "binomial",
                           alpha = .hp("SL.glmnet_participation", "alpha"),
                           nfolds = as.integer(.hp("SL.glmnet_participation", "nfolds")),
                           nlambda = as.integer(.hp("SL.glmnet_participation", "nlambda")),
                           standardize = FALSE)
  object <- list(object = fit, lambda = .hp("SL.glmnet_participation", "lambda_rule"))
  class(object) <- "SL.glmnet_participation"
  pred <- predict(object, newdata = newX)
  announce_fit("SL.glmnet_participation", nrow(X), started)
  list(pred = pred, fit = object)
}

predict.SL.glmnet_participation <- function(object, newdata, ...) {
  as.numeric(stats::predict(object$object, newx = as.matrix(newdata),
                            s = object$lambda, type = "response"))
}

SL.bart_participation <- function(Y, X, newX, family, obsWeights, ...) {
  started <- Sys.time()
  invisible(gc(full = TRUE))
  tr <- as.data.frame(X)
  tr$.outcome <- Y
  te <- as.data.frame(newX)
  fit <- dbarts::bart2(.outcome ~ ., data = tr, test = te,
                       weights = as.numeric(obsWeights),
                       n.trees = as.integer(.hp("SL.bart_participation", "n_trees")),
                       k = .hp("SL.bart_participation", "k"),
                       power = .hp("SL.bart_participation", "power"),
                       base = .hp("SL.bart_participation", "base"),
                       n.threads = learner_threads(), verbose = FALSE,
                       keepTrainingFits = .hp("SL.bart_participation", "keep_training_fits"))
  pred <- tryCatch(as.numeric(stats::fitted(fit, type = "ev", sample = "test")),
                   error = function(e) NULL)
  if (is.null(pred) || length(pred) != nrow(newX)) {
    yt <- fit$yhat.test
    pred <- as.numeric(stats::pnorm(apply(yt, length(dim(yt)), mean)))
  }
  object <- list(object = NULL)
  class(object) <- "SL.bart_participation"
  announce_fit("SL.bart_participation", nrow(X), started)
  list(pred = pred, fit = object)
}

predict.SL.bart_participation <- function(object, newdata, ...) {
  stop("SL.bart_participation does not store its sampler, so it cannot predict new rows ",
       "after the fit. It is used only inside the cross-fit, which supplies the ",
       "rows to predict at fit time.", call. = FALSE)
}

SL.nnet_participation <- function(Y, X, newX, family, obsWeights, ...) {
  started <- Sys.time()
  fit <- nnet::nnet(x = X, y = Y, weights = obsWeights,
                    size = as.integer(.hp("SL.nnet_participation", "size")),
                    decay = .hp("SL.nnet_participation", "decay"),
                    entropy = .hp("SL.nnet_participation", "entropy"),
                    maxit = as.integer(.hp("SL.nnet_participation", "maxit")),
                    MaxNWts = as.integer(.hp("SL.nnet_participation", "max_nwts")),
                    trace = FALSE)
  object <- list(object = fit)
  class(object) <- "SL.nnet_participation"
  pred <- predict(object, newdata = newX)
  announce_fit("SL.nnet_participation", nrow(X), started)
  list(pred = pred, fit = object)
}

predict.SL.nnet_participation <- function(object, newdata, ...) {
  as.numeric(stats::predict(object$object, newdata = newdata, type = "raw"))
}

PARTICIPATION_SL_LIBRARY <- names(
  framework_setting("phase3", "learners", simplify = FALSE))
### The retry below fires on the synthetic data only. 
PARTICIPATION_SL_METHOD <- local({
  m <- do.call(framework_setting("phase3", "meta_learner"), list())
  solve_blend <- m$computeCoef
  m$computeCoef <- function(Z, Y, libraryNames, verbose, obsWeights, control, ...) {
    out <- try(solve_blend(Z = Z, Y = Y, libraryNames = libraryNames,
                           verbose = verbose, obsWeights = obsWeights,
                           control = control, ...), silent = TRUE)
    if (!inherits(out, "try-error")) return(out)
    bound <- framework_setting("phase3", "meta_learner_retry_trim")
    message(sprintf(paste0("[crossfit] the blend did not solve with the library ",
                           "predictions bounded at %g, retrying at %g"),
                    control$trimLogit, bound))
    control$trimLogit <- bound
    solve_blend(Z = Z, Y = Y, libraryNames = libraryNames, verbose = verbose,
                obsWeights = obsWeights, control = control, ...)
  }
  m
})

# ---- cross-fitted fits -------------------------------------------------------

participation_matrix <- function(df, vars) {
  d <- droplevels(df[, vars, drop = FALSE])
  is_factor <- vapply(d, is.factor, logical(1))
  contrasts_arg <- if (any(is_factor)) {
    stats::setNames(lapply(vars[is_factor], function(v) {
      stats::contr.treatment(levels(d[[v]]), contrasts = FALSE)
    }), vars[is_factor])
  } else NULL
  formula <- stats::as.formula(paste("~", paste(vars, collapse = "+")))
  m <- stats::model.matrix(formula, data = d, contrasts.arg = contrasts_arg)
  m <- m[, colnames(m) != "(Intercept)", drop = FALSE]
  colnames(m) <- make.names(colnames(m), unique = TRUE)
  m
}

participation_design <- function(x, continuous) {
  continuous <- intersect(continuous, colnames(x))
  for (v in continuous) {
    col <- as.numeric(x[[v]])
    s <- stats::sd(col)
    if (!is.finite(s) || s == 0) {
      stop("The continuous auxiliary ", v, " takes one value on the stacked ",
           "cohort and reference sample, so it cannot be standardised.",
           call. = FALSE)
    }
    x[[v]] <- (col - mean(col)) / s
  }
  message(sprintf("[design] %d columns, %d standardised: %s", ncol(x),
                  length(continuous),
                  if (length(continuous)) paste(continuous, collapse = ", ")
                  else "none, every auxiliary is categorical"))
  x
}

crossfit_superlearner <- function(y, x, obs_weights, outer_folds = NULL) {
  sl_library <- PARTICIPATION_SL_LIBRARY
  workers <- crossfit_fold_workers()
  threads <- learner_threads()
  options(mc.cores = workers)

  message(sprintf("[crossfit] %s rows, %d columns, %d learners, %d outer by %d inner folds",
                  format(length(y), big.mark = ","), ncol(x), length(sl_library),
                  CROSSFIT_V, CROSSFIT_V))
  blas <- Sys.getenv("OPENBLAS_NUM_THREADS", "unset")
  message(sprintf("[crossfit] fold workers %d, threads per fit %d, cores %d, blas threads %s, fits per outer fold %d",
                  workers, threads, workers * threads, blas,
                  length(sl_library) * (CROSSFIT_V + 1L)))

  inner <- list(V = CROSSFIT_V, stratifyCV = TRUE, shuffle = TRUE)

  if (workers > 1L && !fork_is_safe()) {
    message("[crossfit] the canary fork did not return, so this process cannot ",
            "fork safely. Falling back to one fold worker, which does not fork. ",
            "The fit will be slower and it will finish.")
    workers <- 1L
    options(mc.cores = 1L)
  }

  set.seed(as.integer(framework_setting("phase3", "learner_seed")))
  fit <- SuperLearner::CV.SuperLearner(
    Y = y, X = x, family = stats::binomial(),
    SL.library = sl_library, method = PARTICIPATION_SL_METHOD,
    cvControl = if (is.null(outer_folds)) {
      list(V = CROSSFIT_V, stratifyCV = TRUE, shuffle = TRUE)
    } else {
      list(V = length(outer_folds), validRows = outer_folds)
    },
    innerCvControl = rep(list(inner),
                         if (is.null(outer_folds)) CROSSFIT_V else length(outer_folds)),
    obsWeights = as.numeric(obs_weights),
    control = list(saveFitLibrary = FALSE), saveAll = FALSE,
    parallel = if (workers > 1L) "multicore" else "seq"
  )

  lp <- fit$library.predict
  dead <- vapply(seq_len(ncol(lp)),
                 function(k) length(unique(lp[, k])) <= 1L || all(!is.finite(lp[, k])),
                 logical(1))
  if (any(dead)) {
    stop("These learners produced no predictions and the library is not the one ",
         "this run is named after: ", paste(colnames(lp)[dead], collapse = ", "),
         ". Look for \"Error in algorithm\" in the fit log.", call. = FALSE)
  }

  prob <- as.numeric(fit$SL.predict)
  if (length(prob) != length(y)) {
    stop("Cross-fitted participation probabilities are length-mismatched.")
  }
  list(prob_ukb = prob, fit = fit)
}

crossfit_fold_index <- function(fit, n) {
  fold <- integer(n)
  for (k in seq_along(fit$folds)) fold[fit$folds[[k]]] <- k
  fold
}
