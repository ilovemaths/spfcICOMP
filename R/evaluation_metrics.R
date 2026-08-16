#' Compute Root Mean Squared Error
#'
#' @param observed Numeric observed values.
#' @param predicted Numeric predicted values.
#'
#' @return Numeric RMSE value.
#'
#' @noRd
rmse <- function(observed, predicted){

  observed <- as.numeric(observed)
  predicted <- as.numeric(predicted)

  if(length(observed) != length(predicted)){
    stop("observed and predicted must have the same length.")
  }

  sqrt(mean((observed - predicted)^2))
}


#' Compute Mean Absolute Error
#'
#' @param observed Numeric observed values.
#' @param predicted Numeric predicted values.
#'
#' @return Numeric MAE value.
#'
#' @noRd
mae <- function(observed, predicted){

  observed <- as.numeric(observed)
  predicted <- as.numeric(predicted)

  if(length(observed) != length(predicted)){
    stop("observed and predicted must have the same length.")
  }

  mean(abs(observed - predicted))
}


#' Compute Binary Classification Metrics
#'
#' Computes confusion matrix, accuracy, sensitivity, specificity, precision,
#' and F1 score for binary classification.
#'
#' @param observed Observed binary response.
#' @param predicted Predicted binary class labels.
#' @param positive Optional positive class label. If `NULL`, the second factor
#' level is treated as the positive class.
#'
#' @return A list containing the confusion matrix and classification metrics.
#'
#' @examples
#' obs <- factor(c('A', 'A', 'B', 'B'))
#' pred <- factor(c('A', 'B', 'B', 'B'), levels = levels(obs))
#' binary_classification_metrics(obs, pred, positive = 'B')
#' @export
binary_classification_metrics <- function(
    observed,
    predicted,
    positive = NULL
){

  observed <- as.factor(observed)
  predicted <- as.factor(predicted)

  if(length(levels(observed)) != 2){
    stop("observed must have exactly two classes.")
  }

  if(is.null(positive)){
    positive <- levels(observed)[2]
  }

  predicted <- factor(
    predicted,
    levels = levels(observed)
  )

  tab <- table(
    observed = observed,
    predicted = predicted
  )

  tn <- tab[1, 1]
  fp <- tab[1, 2]
  fn <- tab[2, 1]
  tp <- tab[2, 2]

  accuracy <- (tp + tn) / sum(tab)

  sensitivity <- ifelse((tp + fn) == 0, NA_real_, tp / (tp + fn))
  specificity <- ifelse((tn + fp) == 0, NA_real_, tn / (tn + fp))
  precision <- ifelse((tp + fp) == 0, NA_real_, tp / (tp + fp))
  f1 <- ifelse(
    is.na(precision) || is.na(sensitivity) || (precision + sensitivity) == 0,
    NA_real_,
    2 * precision * sensitivity / (precision + sensitivity)
  )

  balanced_accuracy <- mean(
    c(sensitivity, specificity),
    na.rm = TRUE
  )

  return(list(
    confusion_matrix = tab,
    accuracy = as.numeric(accuracy),
    sensitivity = as.numeric(sensitivity),
    specificity = as.numeric(specificity),
    precision = as.numeric(precision),
    f1 = as.numeric(f1),
    balanced_accuracy = as.numeric(balanced_accuracy),
    positive = positive
  ))
}


#' Predict from a Reduced Model Fit
#'
#' Generates fitted values, probabilities, or predicted classes from an object
#' returned by `fit_reduced_model()`.
#'
#' @param object Object returned by `fit_reduced_model()`.
#' @param Z Numeric reduced-score matrix.
#' @param type Character. One of `"response"`, `"class"`, or `"link"`.
#' @param threshold Classification threshold for binary probabilities.
#'
#' @return Numeric predictions or factor class labels.
#'
#' @examples
#' set.seed(1)
#' Z <- matrix(rnorm(40), nrow = 20, ncol = 2)
#' y <- Z[, 1] - 0.5 * Z[, 2] + rnorm(20, sd = 0.4)
#' mod <- fit_reduced_model(Z, y, ytype = 'continuous')
#' head(predict_reduced_model(mod, Z, type = 'response'))
#' @export
predict_reduced_model <- function(
    object,
    Z,
    type = c("response", "class", "link"),
    threshold = 0.5
){

  type <- match.arg(type)

  Z <- as.matrix(Z)

  if(is.null(colnames(Z))){
    colnames(Z) <- paste0("Z", seq_len(ncol(Z)))
  }

  if(object$ytype == "continuous"){

    newdata <- as.data.frame(Z)

    pred <- stats::predict(
      object$model,
      newdata = newdata
    )

    return(as.numeric(pred))
  }

  if(object$ytype == "categorical"){

    if(object$method %in% c("glm", "glm_single_dimension")){

      newdata <- as.data.frame(Z)

      prob <- stats::predict(
        object$model,
        newdata = newdata,
        type = "response"
      )

      prob <- as.numeric(prob)

      if(type == "response"){
        return(prob)
      }

      if(type == "class"){
        pred <- ifelse(prob >= threshold, 1, 0)
        return(factor(pred))
      }

      if(type == "link"){
        link <- stats::predict(
          object$model,
          newdata = newdata,
          type = "link"
        )
        return(as.numeric(link))
      }
    }

    if(object$method == "ridge_logistic"){

      prob <- stats::predict(
        object$model,
        newdata = Z,
        type = "response"
      )

      prob <- as.numeric(prob)

      if(type == "response"){
        return(prob)
      }

      if(type == "class"){

        pred <- ifelse(
          prob >= threshold,
          object$positive,
          setdiff(object$levels, object$positive)[1]
        )

        return(
          factor(
            pred,
            levels = object$levels
          )
        )
      }

      if(type == "link"){

        prob <- pmin(
          pmax(prob, 1e-12),
          1 - 1e-12
        )

        return(
          log(prob / (1 - prob))
        )
      }
    }
  }

  stop("Unsupported reduced model object.")
}


#' Evaluate a Reduced Model
#'
#' Evaluates a model fitted on reduced SPFC scores.
#'
#' @param object Object returned by `fit_reduced_model()`.
#' @param Z Numeric reduced-score matrix.
#' @param y Observed response.
#' @param threshold Classification threshold for binary responses.
#' @param positive Optional positive class label.
#'
#' @return A list of evaluation metrics.
#'
#' @examples
#' set.seed(1)
#' Z <- matrix(rnorm(40), nrow = 20, ncol = 2)
#' y <- Z[, 1] - 0.5 * Z[, 2] + rnorm(20, sd = 0.4)
#' mod <- fit_reduced_model(Z, y, ytype = 'continuous')
#' evaluate_reduced_model(mod, Z, y)
#' @export
evaluate_reduced_model <- function(
    object,
    Z,
    y,
    threshold = 0.5,
    positive = NULL
){

  if(object$ytype == "continuous"){

    pred <- predict_reduced_model(
      object = object,
      Z = Z,
      type = "response"
    )

    return(list(
      ytype = "continuous",
      rmse = rmse(y, pred),
      mae = mae(y, pred),
      observed = as.numeric(y),
      predicted = pred
    ))
  }

  if(object$ytype == "categorical"){

    prob <- predict_reduced_model(
      object = object,
      Z = Z,
      type = "response",
      threshold = threshold
    )

    pred_num <- ifelse(prob >= threshold, 1, 0)

    observed_factor <- as.factor(y)

    pred_factor <- factor(
      levels(observed_factor)[pred_num + 1],
      levels = levels(observed_factor)
    )

    metrics <- binary_classification_metrics(
      observed = observed_factor,
      predicted = pred_factor,
      positive = positive
    )

    metrics$ytype <- "categorical"
    metrics$probability <- prob
    metrics$predicted <- pred_factor
    metrics$observed <- observed_factor

    return(metrics)
  }

  stop("Unsupported response type.")
}
