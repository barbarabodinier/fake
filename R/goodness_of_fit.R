#' Concordance statistic
#'
#' Computes the concordance statistic given observed binary outcomes and
#' predicted probabilities of event. In logistic regression, the concordance
#' statistic is equal to the area under the Receiver Operating Characteristic
#' (ROC) curve and estimates the probability that an individual who experienced
#' the event (\eqn{Y_i=1}) had a higher probability of event than an individual
#' who did not experience the event (\eqn{Y_i=0}).
#'
#' @param observed vector of binary outcomes.
#' @param predicted vector of predicted probabilities.
#'
#' @return The concordance statistic.
#'
#' @family goodness of fit functions
#'
#' @examples
#'
#' # Data simulation
#' set.seed(1)
#' proba <- runif(n = 200)
#' ydata <- rbinom(n = length(proba), size = 1, prob = proba)
#'
#' # Observed concordance in simulated data
#' Concordance(observed = ydata, predicted = proba)
#'
#' @export
Concordance <- function(observed, predicted) {
  concordant_pairs <- 0
  number_pairs <- 0
  for (i in 1:(length(observed) - 1)) {
    for (j in (i + 1):length(observed)) {
      if (observed[i] != observed[j]) {
        number_pairs <- number_pairs + 1
        if (observed[i] > observed[j]) {
          if (predicted[i] > predicted[j]) {
            concordant_pairs <- concordant_pairs + 1
          }
        } else {
          if (predicted[j] > predicted[i]) {
            concordant_pairs <- concordant_pairs + 1
          }
        }
      }
    }
  }
  return(concordant_pairs / number_pairs)
}


#' True and False Positive Rates
#'
#' Computes the True and False Positive Rates by comparing the true (observed) and
#' predicted status. The predicted status is obtained by applying a threshold on
#' the predicted scores.
#'
#' @inheritParams ROC
#' @param thr threshold for predicted scores.
#'
#' @return True and False Positive Rates (TPR and FPR, respectively).
#'
#' @keywords internal
Rates <- function(observed, predicted, thr) {
  contingency <- table(
    factor(predicted > thr, levels = c(FALSE, TRUE)),
    factor(observed, levels = c(0, 1))
  )

  TP <- contingency[2, 2]
  P <- sum(contingency[, 2])
  TPR <- TP / P

  FP <- contingency[2, 1]
  N <- sum(contingency[, 1])
  FPR <- FP / N

  return(list(TPR = TPR, FPR = FPR))
}


#' Receiver Operating Characteristic (ROC)
#'
#' Computes the True and False Positive Rates (TPR and FPR, respectively) and
#' Area Under the Curve (AUC) by comparing the true (observed) and predicted
#' status using a range of thresholds on the predicted score.
#'
#' @param observed vector of binary outcomes.
#' @param predicted vector of predicted scores.
#' @param n_thr number of thresholds to use to construct the ROC curve. For
#'   faster computations on large data, values below \code{length(x)-1} can be
#'   used.
#'
#' @return A list with: \item{TPR}{True Positive Rate.} \item{FPR}{False
#'   Positive Rate.} \item{AUC}{Area Under the Curve.}
#'
#' @family goodness of fit functions
#'
#' @examples
#' \donttest{
#' # Data simulation
#' set.seed(1)
#' simul <- SimulateRegression(
#'   n = 500, pk = 20,
#'   family = "binomial", ev_xy = 0.8
#' )
#'
#' # Logistic regression
#' fitted <- glm(simul$ydata ~ simul$xdata, family = "binomial")$fitted.values
#'
#' # Constructing the ROC curve
#' roc <- ROC(predicted = fitted, observed = simul$ydata)
#' plot(roc)
#' }
#' @export
ROC <- function(observed, predicted, n_thr = NULL) {
  # Checking the inputs
  predicted <- as.numeric(predicted)
  if (is.factor(observed)) {
    observed <- factor(observed, levels = levels(observed), labels = c(0, 1))
  } else {
    observed <- factor(observed, levels = sort(unique(observed)), labels = c(0, 1))
  }

  # Defining the thresholds
  breaks <- sort(unique(predicted), decreasing = FALSE)
  if (length(breaks) <= 1) {
    message("The predicted value is the same for all observations.")
    FPR <- TPR <- AUC <- NA
  } else {
    breaks <- breaks[-length(breaks)]
    if (!is.null(n_thr)) {
      if (length(breaks) > n_thr) {
        breaks <- breaks[floor(seq(1, length(breaks), length.out = n_thr))]
      } else {
        breaks <- sort(c(breaks, seq(min(breaks), max(breaks), length.out = n_thr - length(breaks))))
      }
    }

    # Computing
    TPR <- FPR <- rep(NA, length(breaks) + 2)
    for (k in 1:length(breaks)) {
      out <- Rates(observed = observed, predicted = predicted, thr = breaks[k])
      TPR[k + 1] <- out$TPR
      FPR[k + 1] <- out$FPR
    }
    TPR[1] <- FPR[1] <- 1
    TPR[length(TPR)] <- FPR[length(FPR)] <- 0

    # Computing the AUC
    tmp <- apply(rbind(TPR[-1], TPR[-length(TPR)]), 2, mean)
    AUC <- abs(sum(diff(FPR) * tmp))
  }

  # Preparing output
  out <- list(FPR = rbind(FPR), TPR = rbind(TPR), AUC = AUC)

  # Defining class
  class(out) <- "roc_curve"

  return(out)
}
