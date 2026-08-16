# 05_riboflavin_analysis.R
# Definitive Riboflavin regression application.
# Frozen thesis analysis; compatible with spfcICOMP v0.0.0.9008, v0.0.0.9009, and the planned v0.1.0.
#
# Replaces the superseded v9006 prototype.
# Uses the actual v9008 selector interface:
#   sel$criteria, sel$selected, sel$fits.
#
# Full-data analysis:
#   OAS and MEC; d = 1,...,5;
#   AIC, BIC, CAIC, ICOMP_IFIM, CICOMP;
#   C1F ICOMP-HD adaptive weighted L1 feature selection.
#
# Prediction:
#   repeated 5-fold CV;
#   dimension selection is performed inside each training fold;
#   predictions are made only on held-out observations.

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


load_riboflavin_benchmark <- function() {
  # Prefer an installed hdi copy when available, but do not require it.
  # hdi is archived on CRAN for current R versions, so the reproducibility
  # fallback retrieves the frozen hdi 0.1-10 source archive directly.
  if (requireNamespace("hdi", quietly = TRUE)) {
    e <- new.env(parent = emptyenv())
    utils::data(list = "riboflavin", package = "hdi", envir = e)
    if (exists("riboflavin", envir = e, inherits = FALSE)) {
      return(get("riboflavin", envir = e, inherits = FALSE))
    }
  }

  archive_url <- paste0(
    "https://cran.r-project.org/src/contrib/Archive/hdi/",
    "hdi_0.1-10.tar.gz"
  )

  message(
    "Package 'hdi' is not installed; retrieving archived hdi 0.1-10 ",
    "Riboflavin data from CRAN."
  )

  td <- tempfile("hdi_archive_")
  dir.create(td, recursive = TRUE)
  tarball <- file.path(td, "hdi_0.1-10.tar.gz")

  utils::download.file(
    archive_url,
    tarball,
    mode = "wb",
    quiet = FALSE
  )
  utils::untar(tarball, exdir = td)

  rib_files <- list.files(
    td,
    pattern = "^riboflavin\\.(rda|RData|rdata)$",
    recursive = TRUE,
    full.names = TRUE,
    ignore.case = TRUE
  )

  if (length(rib_files) != 1L) {
    stop(
      "Could not uniquely identify archived Riboflavin data in hdi 0.1-10."
    )
  }

  e <- new.env(parent = emptyenv())
  loaded <- load(rib_files[[1]], envir = e)

  if (!"riboflavin" %in% loaded ||
      !exists("riboflavin", envir = e, inherits = FALSE)) {
    stop("Archived hdi 0.1-10 did not contain object 'riboflavin'.")
  }

  get("riboflavin", envir = e, inherits = FALSE)
}

extract_xy <- function(dat, dataset_name) {
  # hdi 0.1-10 stores Riboflavin as a data frame with columns y and x,
  # where x is an AsIs 71 x 4088 matrix. Handle that structure explicitly.
  if (is.data.frame(dat) &&
      all(c("x", "y") %in% names(dat)) &&
      !is.null(dim(dat[["x"]]))) {
    X <- dat[["x"]]
    if (inherits(X, "AsIs")) {
      X <- unclass(X)
    }
    X <- as.matrix(X)
    y <- dat[["y"]]
  } else if (is.list(dat) && all(c("x", "y") %in% names(dat))) {
    X <- as.matrix(dat[["x"]])
    y <- dat[["y"]]
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


dat <- load_riboflavin_benchmark()
xy <- extract_xy(dat, "riboflavin")
X <- xy$X
y <- as.numeric(xy$y)
predictor_names <- xy$predictor_names

if (!identical(dim(X), c(71L, 4088L))) {
  stop(
    "Unexpected Riboflavin dimensions: expected 71 x 4088 predictors, found ",
    nrow(X), " x ", ncol(X), "."
  )
}
if (nrow(X) != length(y)) stop("Riboflavin X/y size mismatch.")
if (any(!is.finite(X)) || any(!is.finite(y))) stop("Riboflavin contains non-finite values.")

cat("\n============================================================\n")
cat("RIBOFLAVIN SPFC-ICOMP ANALYSIS\n")
cat("============================================================\n")
cat("spfcICOMP version:", installed_version, "\n")
cat("n =", nrow(X), "| p =", ncol(X), "| p/n =", round(ncol(X)/nrow(X), 3), "\n")

criteria_keep <- c("AIC", "BIC", "CAIC", "ICOMP_IFIM", "CICOMP")
cov_methods <- c("oas", "mec")
D_GRID <- 1:5
NSLICES <- 5L
POLY_DEGREE <- 2L  # matches the frozen definitive simulation basis
N_FOLDS <- 5L
N_REPEATS <- 5L
SEED <- 20260810L

extract_selected_d <- function(sel, criterion) {
  hit <- sel$selected[sel$selected$criterion == criterion, , drop = FALSE]
  if (!nrow(hit)) return(NA_integer_)
  as.integer(hit$selected_d[1])
}

extract_candidate <- function(sel, d) {
  nm <- paste0("d", d)
  if (!is.null(sel$fits[[nm]])) return(sel$fits[[nm]])
  idx <- match(d, sel$criteria$d)
  if (is.na(idx)) stop("Candidate d=", d, " not found.")
  sel$fits[[idx]]
}

safe_r2 <- function(obs, pred) {
  den <- sum((obs - mean(obs))^2)
  if (!is.finite(den) || den <= 0) return(NA_real_)
  1 - sum((obs - pred)^2) / den
}

safe_cor <- function(obs, pred) {
  if (length(obs) < 3L || stats::sd(obs) == 0 || stats::sd(pred) == 0) return(NA_real_)
  stats::cor(obs, pred, use = "complete.obs")
}

# ------------------------------------------------------------
# 1. Dataset audit
# ------------------------------------------------------------

v <- apply(X, 2, stats::var)
audit <- data.frame(
  dataset = "Riboflavin",
  response_type = "continuous",
  n = nrow(X),
  p = ncol(X),
  p_over_n = ncol(X) / nrow(X),
  missing_predictor_values = sum(!is.finite(X)),
  missing_response_values = sum(!is.finite(y)),
  zero_variance_predictors = sum(!is.finite(v) | v <= .Machine$double.eps),
  duplicated_rows = sum(duplicated(data.frame(y = y, X, check.names = FALSE))),
  response_min = min(y),
  response_mean = mean(y),
  response_max = max(y)
)

write.csv(
  audit,
  file.path(table_dir, "riboflavin_data_audit.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------
# 2. Full-data dimension selection and feature screening
# ------------------------------------------------------------

full_rows <- list()
variable_rows <- list()
dimension_objects <- list()
ii <- 0L
jj <- 0L

for (cm in cov_methods) {
  cat("Full-data dimension selection:", toupper(cm), "\n")

  sel <- spfcICOMP::spfc_select_dimension(
    X = X,
    y = y,
    d_grid = D_GRID,
    cov_method = cm,
    ytype = "continuous",
    nslices = NSLICES,
    poly_degree = POLY_DEGREE,
    classifier = "auto",
    verbose = FALSE
  )

  dimension_objects[[cm]] <- sel

  score_table <- sel$criteria[, c(
    "cov_method", "ytype", "d", "rho", "AIC", "BIC", "CAIC",
    "ICOMP_IFIM", "ICOMP_MISSPEC", "CICOMP"
  )]

  write.csv(
    score_table,
    file.path(table_dir, paste0("riboflavin_dimension_scores_", cm, ".csv")),
    row.names = FALSE
  )

  for (cr in criteria_keep) {
    d_hat <- extract_selected_d(sel, cr)
    candidate <- extract_candidate(sel, d_hat)
    fit <- candidate$spfc_fit
    reduced_model <- candidate$reduced_model

    vs <- spfcICOMP::spfc_select_variables(
      fit = fit,
      method = "adaptive_weighted_l1",
      selection_rule = "c1f",
      reduced_model = reduced_model,
      c1f_calibration = "icomp_hd_floor"
    )

    selected_table <- vs[vs$selected, , drop = FALSE]
    selected_idx <- selected_table$variable

    ii <- ii + 1L
    full_rows[[ii]] <- data.frame(
      dataset = "Riboflavin",
      cov_method = cm,
      criterion = cr,
      selected_d = d_hat,
      criterion_value = sel$selected$minimum_value[sel$selected$criterion == cr][1],
      shrinkage_rho = fit$rho,
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
          dataset = "Riboflavin",
          cov_method = cm,
          criterion = cr,
          selected_d = d_hat,
          variable = selected_table$variable[k],
          variable_name = predictor_names[selected_table$variable[k]],
          importance = selected_table$importance[k],
          penalty = selected_table$penalty[k],
          shrunk_importance = selected_table$shrunk_importance[k],
          stringsAsFactors = FALSE
        )
      }
    }
  }
}

full_summary <- do.call(rbind, full_rows)
selected_variables <- if (length(variable_rows)) do.call(rbind, variable_rows) else data.frame()

write.csv(
  full_summary,
  file.path(table_dir, "riboflavin_full_data_summary.csv"),
  row.names = FALSE
)

write.csv(
  selected_variables,
  file.path(table_dir, "riboflavin_selected_variables.csv"),
  row.names = FALSE
)

saveRDS(
  list(
    package_version = installed_version,
    audit = audit,
    dimension_objects = dimension_objects,
    full_summary = full_summary,
    selected_variables = selected_variables
  ),
  file.path(output_dir, "riboflavin_full_data_analysis.rds")
)

# ------------------------------------------------------------
# 3. Repeated five-fold cross-validation
# ------------------------------------------------------------

prediction_rows <- list()
kk <- 0L

for (r in seq_len(N_REPEATS)) {
  fold_id <- spfcICOMP::create_cv_folds(
    y = y,
    folds = N_FOLDS,
    seed = SEED + r - 1L,
    stratify = FALSE
  )

  for (fold in seq_len(N_FOLDS)) {
    te <- which(fold_id == fold)
    tr <- which(fold_id != fold)

    Xtr <- X[tr, , drop = FALSE]
    Xte <- X[te, , drop = FALSE]
    ytr <- y[tr]
    yte <- y[te]

    for (cm in cov_methods) {
      cat(
        "CV repeat ", r, "/", N_REPEATS,
        " | fold ", fold, "/", N_FOLDS,
        " | ", toupper(cm), "\n",
        sep = ""
      )

      sel <- spfcICOMP::spfc_select_dimension(
        X = Xtr,
        y = ytr,
        d_grid = D_GRID,
        cov_method = cm,
        ytype = "continuous",
        nslices = NSLICES,
        poly_degree = POLY_DEGREE,
        classifier = "auto",
        verbose = FALSE
      )

      for (cr in criteria_keep) {
        d_hat <- extract_selected_d(sel, cr)
        candidate <- extract_candidate(sel, d_hat)
        fit <- candidate$spfc_fit
        reduced_model <- candidate$reduced_model

        Zte <- stats::predict(fit, newdata = Xte)
        pred <- spfcICOMP::predict_reduced_model(
          object = reduced_model,
          Z = Zte,
          type = "response"
        )

        for (m in seq_along(te)) {
          kk <- kk + 1L
          prediction_rows[[kk]] <- data.frame(
            repeat_id = r,
            fold = fold,
            observation = te[m],
            cov_method = cm,
            criterion = cr,
            selected_d = d_hat,
            observed = yte[m],
            predicted = pred[m],
            stringsAsFactors = FALSE
          )
        }
      }
    }
  }
}

cv_predictions <- do.call(rbind, prediction_rows)

write.csv(
  cv_predictions,
  file.path(output_dir, "riboflavin_cv_predictions.csv"),
  row.names = FALSE
)

# One performance record per repeat x covariance x criterion.
repeat_key <- interaction(
  cv_predictions$repeat_id,
  cv_predictions$cov_method,
  cv_predictions$criterion,
  drop = TRUE,
  lex.order = TRUE
)

cv_by_repeat <- do.call(
  rbind,
  lapply(split(cv_predictions, repeat_key), function(z) {
    err <- z$observed - z$predicted
    data.frame(
      repeat_id = z$repeat_id[1],
      cov_method = z$cov_method[1],
      criterion = z$criterion[1],
      mean_selected_d = mean(z$selected_d),
      rmse = sqrt(mean(err^2)),
      mae = mean(abs(err)),
      r_squared = safe_r2(z$observed, z$predicted),
      correlation = safe_cor(z$observed, z$predicted)
    )
  })
)
rownames(cv_by_repeat) <- NULL

write.csv(
  cv_by_repeat,
  file.path(table_dir, "riboflavin_cv_by_repeat.csv"),
  row.names = FALSE
)

summary_key <- interaction(
  cv_by_repeat$cov_method,
  cv_by_repeat$criterion,
  drop = TRUE,
  lex.order = TRUE
)

cv_summary <- do.call(
  rbind,
  lapply(split(cv_by_repeat, summary_key), function(z) {
    data.frame(
      cov_method = z$cov_method[1],
      criterion = z$criterion[1],
      repeats = nrow(z),
      mean_selected_d = mean(z$mean_selected_d),
      sd_selected_d = stats::sd(z$mean_selected_d),
      mean_rmse = mean(z$rmse),
      sd_rmse = stats::sd(z$rmse),
      mean_mae = mean(z$mae),
      sd_mae = stats::sd(z$mae),
      mean_r_squared = mean(z$r_squared, na.rm = TRUE),
      sd_r_squared = stats::sd(z$r_squared, na.rm = TRUE),
      mean_correlation = mean(z$correlation, na.rm = TRUE),
      sd_correlation = stats::sd(z$correlation, na.rm = TRUE)
    )
  })
)
rownames(cv_summary) <- NULL
cv_summary <- cv_summary[
  order(cv_summary$cov_method, match(cv_summary$criterion, criteria_keep)),
]

write.csv(
  cv_summary,
  file.path(table_dir, "riboflavin_cv_summary.csv"),
  row.names = FALSE
)

saveRDS(
  list(
    predictions = cv_predictions,
    by_repeat = cv_by_repeat,
    summary = cv_summary,
    settings = list(
      folds = N_FOLDS,
      repeats = N_REPEATS,
      seed = SEED,
      d_grid = D_GRID,
      nslices = NSLICES,
      poly_degree = POLY_DEGREE,
      criteria = criteria_keep,
      cov_methods = cov_methods
    )
  ),
  file.path(output_dir, "riboflavin_cv_analysis.rds")
)

cat("\nRiboflavin analysis complete.\n")
print(full_summary)
print(cv_summary)
