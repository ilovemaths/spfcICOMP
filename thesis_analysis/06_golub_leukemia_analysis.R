# 06_golub_leukemia_analysis.R
# Definitive Golub leukaemia binary-classification application.
# Frozen thesis analysis; compatible with spfcICOMP v0.0.0.9008, v0.0.0.9009, and the planned v0.1.0.
#
# Binary Y provides a one-dimensional inverse-response basis in the current
# validated implementation; therefore d = 1 is fixed rather than pretending
# to compare structural dimensions that are not identifiable from two classes.
#
# OAS and MEC are compared for:
#   - data-adaptive covariance shrinkage,
#   - C1F ICOMP-HD feature selection,
#   - repeated stratified five-fold classification performance.
#
# Final classifier: ridge logistic on SPFC scores.

SUPPORTED_VERSIONS <- c("0.0.0.9008", "0.0.0.9009", "0.1.0")
stopifnot(requireNamespace("spfcICOMP", quietly = TRUE))

installed_version <- as.character(utils::packageVersion("spfcICOMP"))
if (!installed_version %in% SUPPORTED_VERSIONS) {
  stop(
    "This frozen thesis script was validated for spfcICOMP versions ",
    paste(SUPPORTED_VERSIONS, collapse = ", "),
    "; found ", installed_version, "."
  )
}

output_dir <- "thesis_analysis/real_outputs_v9008"
table_dir  <- "thesis_analysis/real_tables_v9008"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)


load_golub_benchmark <- function() {
  if (!requireNamespace("multtest", quietly = TRUE)) {
    stop(
      "The Golub benchmark is not bundled with spfcICOMP. ",
      "Install Bioconductor package 'multtest' (for example with ",
      "BiocManager::install('multtest')) and rerun this script."
    )
  }
  e <- new.env(parent = emptyenv())
  utils::data(list = "golub", package = "multtest", envir = e)
  required <- c("golub", "golub.cl")
  missing_objects <- required[
    !vapply(required, exists, logical(1), envir = e, inherits = FALSE)
  ]
  if (length(missing_objects)) {
    stop(
      "Could not load the following objects from multtest::golub: ",
      paste(missing_objects, collapse = ", ")
    )
  }
  X <- t(get("golub", envir = e, inherits = FALSE))
  y <- get("golub.cl", envir = e, inherits = FALSE)
  storage.mode(X) <- "double"
  colnames(X) <- paste0("x", seq_len(ncol(X)))
  data.frame(y = y, X, check.names = FALSE)
}

extract_xy <- function(dat, dataset_name) {
  if (is.list(dat) && all(c("x", "y") %in% names(dat))) {
    X <- as.matrix(dat$x)
    y <- dat$y
  } else if (is.data.frame(dat) || is.matrix(dat)) {
    dat <- as.data.frame(dat, check.names = FALSE)
    yn <- intersect(
      c("y", "Y", "response", "Response", dataset_name),
      names(dat)
    )
    if (!length(yn)) {
      stop("Could not identify response in ", dataset_name, ".")
    }
    y <- dat[[yn[1]]]
    X <- as.matrix(dat[setdiff(names(dat), yn[1]), drop = FALSE])
  } else {
    stop("Unsupported object structure for ", dataset_name, ".")
  }

  storage.mode(X) <- "double"
  if (is.null(colnames(X))) {
    colnames(X) <- paste0("X", seq_len(ncol(X)))
  }

  list(X = X, y = y, predictor_names = colnames(X))
}


dat <- load_golub_benchmark()
xy <- extract_xy(dat, "golubleukemia")
X <- xy$X
y <- factor(xy$y)
predictor_names <- xy$predictor_names

if (!identical(dim(X), c(38L, 3051L))) {
  stop(
    "Unexpected Golub dimensions: expected 38 x 3051 predictors, found ",
    nrow(X), " x ", ncol(X), "."
  )
}
if (nrow(X) != length(y)) stop("Golub X/y size mismatch.")
if (any(!is.finite(X)) || anyNA(y)) stop("Golub data contain missing/non-finite values.")
if (nlevels(y) != 2L) stop("Golub response must have exactly two classes.")

positive_class <- levels(y)[2]
negative_class <- levels(y)[1]

cat("\n============================================================\n")
cat("GOLUB LEUKAEMIA SPFC-ICOMP ANALYSIS\n")
cat("============================================================\n")
cat("spfcICOMP version:", installed_version, "\n")
cat("n =", nrow(X), "| p =", ncol(X), "| p/n =", round(ncol(X)/nrow(X), 3), "\n")
cat("Class counts:\n")
print(table(y))
cat("Positive class:", positive_class, "\n")

cov_methods <- c("oas", "mec")
D_FIXED <- 1L
N_FOLDS <- 5L
N_REPEATS <- 20L
SEED <- 20260810L

auc_rank <- function(observed, probability, positive) {
  observed <- factor(observed)
  y01 <- as.integer(observed == positive)
  n_pos <- sum(y01 == 1L)
  n_neg <- sum(y01 == 0L)
  if (n_pos == 0L || n_neg == 0L) return(NA_real_)
  ranks <- rank(probability, ties.method = "average")
  (sum(ranks[y01 == 1L]) - n_pos * (n_pos + 1) / 2) / (n_pos * n_neg)
}

# ------------------------------------------------------------
# 1. Dataset audit
# ------------------------------------------------------------

v <- apply(X, 2, stats::var)
class_tab <- table(y)

audit <- data.frame(
  dataset = "Golub leukaemia",
  response_type = "binary",
  n = nrow(X),
  p = ncol(X),
  p_over_n = ncol(X) / nrow(X),
  missing_predictor_values = sum(!is.finite(X)),
  missing_response_values = sum(is.na(y)),
  zero_variance_predictors = sum(!is.finite(v) | v <= .Machine$double.eps),
  duplicated_rows = sum(duplicated(data.frame(y = y, X, check.names = FALSE))),
  class_1 = names(class_tab)[1],
  class_1_n = as.integer(class_tab[1]),
  class_2 = names(class_tab)[2],
  class_2_n = as.integer(class_tab[2]),
  positive_class = positive_class,
  stringsAsFactors = FALSE
)

write.csv(
  audit,
  file.path(table_dir, "golub_data_audit.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------
# 2. Full-data SPFC d=1 fit and C1F feature selection
# ------------------------------------------------------------

full_rows <- list()
variable_rows <- list()
full_objects <- list()
ii <- 0L
jj <- 0L

for (cm in cov_methods) {
  cat("Full-data d=1 fit:", toupper(cm), "\n")

  fit <- spfcICOMP::spfc_fit(
    X = X,
    y = y,
    d = D_FIXED,
    ytype = "categorical",
    cov_method = cm
  )

  # Ordinary logistic companion model is used for information-criterion scores.
  reduced_glm <- suppressWarnings(
    spfcICOMP::fit_reduced_model(
      Z = fit$Z,
      y = y,
      ytype = "categorical",
      classifier = "glm"
    )
  )

  criteria <- spfcICOMP::score_information_criteria(
    reduced_model = reduced_glm,
    n = nrow(X)
  )

  # Ridge logistic is the validated robust predictive classifier.
  reduced_ridge <- spfcICOMP::fit_reduced_model(
    Z = fit$Z,
    y = y,
    ytype = "categorical",
    classifier = "ridge_logistic"
  )

  vs <- spfcICOMP::spfc_select_variables(
    fit = fit,
    method = "adaptive_weighted_l1",
    selection_rule = "c1f",
    reduced_model = reduced_ridge,
    c1f_calibration = "icomp_hd_floor"
  )

  selected_table <- vs[vs$selected, , drop = FALSE]
  selected_idx <- selected_table$variable

  ii <- ii + 1L
  full_rows[[ii]] <- data.frame(
    dataset = "Golub leukaemia",
    cov_method = cm,
    fixed_d = D_FIXED,
    shrinkage_rho = fit$rho,
    AIC = criteria$AIC,
    BIC = criteria$BIC,
    CAIC = criteria$CAIC,
    ICOMP_IFIM = criteria$ICOMP_IFIM,
    ICOMP_MISSPEC = criteria$ICOMP_MISSPEC,
    CICOMP = criteria$CICOMP,
    c1f = attr(vs, "C1F"),
    c1f_global_penalty = attr(vs, "global_penalty"),
    n_selected = length(selected_idx),
    selected_fraction = length(selected_idx) / ncol(X),
    stringsAsFactors = FALSE
  )

  if (nrow(selected_table)) {
    for (k in seq_len(nrow(selected_table))) {
      jj <- jj + 1L
      variable_rows[[jj]] <- data.frame(
        dataset = "Golub leukaemia",
        cov_method = cm,
        fixed_d = D_FIXED,
        variable = selected_table$variable[k],
        variable_name = predictor_names[selected_table$variable[k]],
        importance = selected_table$importance[k],
        penalty = selected_table$penalty[k],
        shrunk_importance = selected_table$shrunk_importance[k],
        stringsAsFactors = FALSE
      )
    }
  }

  full_objects[[cm]] <- list(
    fit = fit,
    reduced_glm = reduced_glm,
    reduced_ridge = reduced_ridge,
    criteria = criteria,
    variable_selection = vs
  )
}

full_summary <- do.call(rbind, full_rows)
selected_variables <- if (length(variable_rows)) do.call(rbind, variable_rows) else data.frame()

write.csv(
  full_summary,
  file.path(table_dir, "golub_full_data_summary.csv"),
  row.names = FALSE
)

write.csv(
  selected_variables,
  file.path(table_dir, "golub_selected_variables.csv"),
  row.names = FALSE
)

saveRDS(
  list(
    package_version = installed_version,
    audit = audit,
    full_objects = full_objects,
    full_summary = full_summary,
    selected_variables = selected_variables,
    positive_class = positive_class
  ),
  file.path(output_dir, "golub_full_data_analysis.rds")
)

# ------------------------------------------------------------
# 3. Repeated stratified five-fold CV
# ------------------------------------------------------------

prediction_rows <- list()
kk <- 0L

for (r in seq_len(N_REPEATS)) {
  fold_id <- spfcICOMP::create_cv_folds(
    y = y,
    folds = N_FOLDS,
    seed = SEED + r - 1L,
    stratify = TRUE
  )

  for (fold in seq_len(N_FOLDS)) {
    te <- which(fold_id == fold)
    tr <- which(fold_id != fold)

    Xtr <- X[tr, , drop = FALSE]
    Xte <- X[te, , drop = FALSE]
    ytr <- droplevels(y[tr])
    yte <- factor(y[te], levels = levels(y))

    # Defensive check: every training fold must contain both classes.
    if (nlevels(ytr) != 2L) {
      stop("A Golub training fold lost one response class.")
    }

    for (cm in cov_methods) {
      cat(
        "CV repeat ", r, "/", N_REPEATS,
        " | fold ", fold, "/", N_FOLDS,
        " | ", toupper(cm), "\n",
        sep = ""
      )

      fit <- spfcICOMP::spfc_fit(
        X = Xtr,
        y = ytr,
        d = D_FIXED,
        ytype = "categorical",
        cov_method = cm
      )

      reduced <- spfcICOMP::fit_reduced_model(
        Z = fit$Z,
        y = ytr,
        ytype = "categorical",
        classifier = "ridge_logistic"
      )

      Zte <- stats::predict(fit, newdata = Xte)

      prob <- spfcICOMP::predict_reduced_model(
        object = reduced,
        Z = Zte,
        type = "response"
      )

      pred <- spfcICOMP::predict_reduced_model(
        object = reduced,
        Z = Zte,
        type = "class",
        threshold = 0.5
      )

      pred <- factor(as.character(pred), levels = levels(y))

      for (m in seq_along(te)) {
        kk <- kk + 1L
        prediction_rows[[kk]] <- data.frame(
          repeat_id = r,
          fold = fold,
          observation = te[m],
          cov_method = cm,
          fixed_d = D_FIXED,
          observed = as.character(yte[m]),
          predicted = as.character(pred[m]),
          probability_positive = prob[m],
          stringsAsFactors = FALSE
        )
      }
    }
  }
}

cv_predictions <- do.call(rbind, prediction_rows)

write.csv(
  cv_predictions,
  file.path(output_dir, "golub_cv_predictions.csv"),
  row.names = FALSE
)

# One metric vector per repeat x covariance method.
repeat_key <- interaction(
  cv_predictions$repeat_id,
  cv_predictions$cov_method,
  drop = TRUE,
  lex.order = TRUE
)

cv_by_repeat <- do.call(
  rbind,
  lapply(split(cv_predictions, repeat_key), function(z) {
    observed <- factor(z$observed, levels = levels(y))
    predicted <- factor(z$predicted, levels = levels(y))

    met <- spfcICOMP::binary_classification_metrics(
      observed = observed,
      predicted = predicted,
      positive = positive_class
    )

    data.frame(
      repeat_id = z$repeat_id[1],
      cov_method = z$cov_method[1],
      accuracy = met$accuracy,
      balanced_accuracy = met$balanced_accuracy,
      sensitivity = met$sensitivity,
      specificity = met$specificity,
      precision = met$precision,
      f1 = met$f1,
      auc = auc_rank(
        observed = observed,
        probability = z$probability_positive,
        positive = positive_class
      )
    )
  })
)
rownames(cv_by_repeat) <- NULL

write.csv(
  cv_by_repeat,
  file.path(table_dir, "golub_cv_by_repeat.csv"),
  row.names = FALSE
)

cv_summary <- do.call(
  rbind,
  lapply(split(cv_by_repeat, cv_by_repeat$cov_method), function(z) {
    data.frame(
      cov_method = z$cov_method[1],
      repeats = nrow(z),
      mean_accuracy = mean(z$accuracy),
      sd_accuracy = stats::sd(z$accuracy),
      mean_balanced_accuracy = mean(z$balanced_accuracy),
      sd_balanced_accuracy = stats::sd(z$balanced_accuracy),
      mean_sensitivity = mean(z$sensitivity, na.rm = TRUE),
      sd_sensitivity = stats::sd(z$sensitivity, na.rm = TRUE),
      mean_specificity = mean(z$specificity, na.rm = TRUE),
      sd_specificity = stats::sd(z$specificity, na.rm = TRUE),
      mean_precision = mean(z$precision, na.rm = TRUE),
      sd_precision = stats::sd(z$precision, na.rm = TRUE),
      mean_f1 = mean(z$f1, na.rm = TRUE),
      sd_f1 = stats::sd(z$f1, na.rm = TRUE),
      mean_auc = mean(z$auc, na.rm = TRUE),
      sd_auc = stats::sd(z$auc, na.rm = TRUE)
    )
  })
)
rownames(cv_summary) <- NULL

write.csv(
  cv_summary,
  file.path(table_dir, "golub_cv_summary.csv"),
  row.names = FALSE
)

saveRDS(
  list(
    predictions = cv_predictions,
    by_repeat = cv_by_repeat,
    summary = cv_summary,
    settings = list(
      fixed_d = D_FIXED,
      folds = N_FOLDS,
      repeats = N_REPEATS,
      seed = SEED,
      cov_methods = cov_methods,
      positive_class = positive_class,
      classifier = "ridge_logistic"
    )
  ),
  file.path(output_dir, "golub_cv_analysis.rds")
)

cat("\nGolub leukaemia analysis complete.\n")
print(full_summary)
print(cv_summary)
