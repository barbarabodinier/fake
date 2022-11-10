#' Simulation of symmetric matrix with block structure
#'
#' Simulates a symmetric matrix with block structure. If
#' \code{continuous=FALSE}, matrix entries are sampled from a discrete uniform
#' distribution taking values in \code{v_within} (for entries in the diagonal
#' block) or \code{v_between} (for entries in off-diagonal blocks). If
#' \code{continuous=TRUE}, entries are sampled from a continuous uniform
#' distribution taking values in the range given by \code{v_within} or
#' \code{v_between}.
#'
#' @param pk vector of the number of variables per group, defining the block
#'   structure.
#' @param v_within vector defining the (range of) nonzero entries in the
#'   diagonal blocks. If \code{continuous=FALSE}, \code{v_within} is the set of
#'   possible values. If \code{continuous=FALSE}, \code{v_within} is the range
#'   of possible values.
#' @param v_between vector defining the (range of) nonzero entries in the
#'   off-diagonal blocks. If \code{continuous=FALSE}, \code{v_between} is the
#'   set of possible precision values. If \code{continuous=FALSE},
#'   \code{v_between} is the range of possible precision values. This argument
#'   is only used if \code{length(pk)>1}.
#' @param v_sign vector of possible signs for matrix entries. Possible inputs
#'   are: \code{-1} for negative entries only, \code{1} for positive entries
#'   only, or \code{c(-1, 1)} for both positive and negative entries.
#' @param continuous logical indicating whether to sample precision values from
#'   a uniform distribution between the minimum and maximum values in
#'   \code{v_within} (diagonal blocks) or \code{v_between} (off-diagonal blocks)
#'   (\code{continuous=TRUE}) or from proposed values in \code{v_within}
#'   (diagonal blocks) or \code{v_between} (off-diagonal blocks)
#'   (\code{continuous=FALSE}).
#'
#' @return A symmetric matrix with uniformly distributed entries sampled from
#'   different distributions for diagonal and off-diagonal blocks.
#'
#' @keywords internal
SimulateSymmetricMatrix <- function(pk = 10,
                                    v_within = c(0.5, 1), v_between = c(0, 0.1),
                                    v_sign = c(-1, 1), continuous = FALSE) {
  # Creating matrix with block indices
  bigblocks <- BlockMatrix(pk)
  bigblocks_vect <- bigblocks[upper.tri(bigblocks)]

  # Making as factor to allow for groups with 1 variable (for clustering)
  bigblocks_vect <- factor(bigblocks_vect, levels = seq(1, max(bigblocks)))
  block_ids <- unique(as.vector(bigblocks))

  # Building absolute v matrix
  v <- bigblocks
  v_vect <- v[upper.tri(v)]
  for (k in block_ids) {
    if (k %in% v_vect) {
      if (k %in% unique(diag(bigblocks))) {
        if (continuous) {
          v_vect[bigblocks_vect == k] <- stats::runif(sum(bigblocks_vect == k), min = min(v_within), max = max(v_within))
        } else {
          v_vect[bigblocks_vect == k] <- base::sample(v_within, size = sum(bigblocks_vect == k), replace = TRUE)
        }
      } else {
        if (continuous) {
          v_vect[bigblocks_vect == k] <- stats::runif(sum(bigblocks_vect == k), min = min(v_between), max = max(v_between))
        } else {
          v_vect[bigblocks_vect == k] <- base::sample(v_between, size = sum(bigblocks_vect == k), replace = TRUE)
        }
      }
    }
  }

  # Sampling the sign of precision entries
  v_vect <- v_vect * base::sample(sort(unique(v_sign)), size = length(v_vect), replace = TRUE)

  # Building v matrix
  diag(v) <- 0
  v[upper.tri(v)] <- v_vect
  v[lower.tri(v)] <- 0
  v <- v + t(v)

  return(v)
}


#' Simulation of binary contribution status
#'
#' Simulates the binary contribution status of potential predictor variables
#' from different blocks to outcome variables. For each outcome, the set of true
#' predictors is sampled from one block of potential predictors. If the blocks
#' of variables are independent, the outcomes will be independent too.
#'
#' @inheritParams SimulateSymmetricMatrix
#' @param q number of outcome variables. By default, one block of predictor is
#'   linked to one outcome, i.e. \code{q=sum(pk)}.
#' @param nu vector of probabilities. Each entry corresponds to one block of
#'   predictors and defines the probability for each predictor within the block
#'   to be chosen as true predictor of the corresponding outcome variable.
#' @param orthogonal logical indicating if the outcomes have to be defined from
#'   independent blocks of predictors as encoded in \code{pk}.
#'
#' @return A binary matrix encoding the contribution status of each predictor
#'   variable (columns) to each outcome variable (rows).
#'
#' @keywords internal
SamplePredictors <- function(pk, q = NULL, nu = 0.1, orthogonal = TRUE) {
  # Definition of the number of outcome variables
  if (is.null(q)) {
    q <- length(pk)
  }
  if (length(nu) != q) {
    nu <- rep(nu[1], q)
  }

  # Simulation of the binary status for true predictors
  theta <- matrix(0, nrow = q, ncol = sum(pk))
  for (k in 1:q) {
    if (orthogonal) {
      if (k > 1) {
        ids <- seq(cumsum(pk)[k - 1] + 1, cumsum(pk)[k])
      } else {
        ids <- seq(1, cumsum(pk)[k])
      }
      theta[k, ids] <- stats::rbinom(pk[k], size = 1, prob = nu[k])

      # Introducing at least one true predictor
      if (sum(theta[k, ids]) == 0) {
        theta[k, sample(ids, size = 1)] <- 1
      }
    } else {
      theta[k, ] <- stats::rbinom(sum(pk), size = 1, prob = nu[k])
      theta[k, k] <- 1
    }
  }

  return(t(theta))
}


#' Maximising matrix contrast
#'
#' Computes the contrast of the correlation matrix obtained by adding u to the
#' diagonal of the precision matrix. This function is used to find the value of
#' u that maximises the contrast when constructing a diagonally dominant
#' precision matrix.
#'
#' @param u constant u added to the diagonal of the precision matrix.
#' @param omega positive semi-definite precision matrix.
#' @param digits number of digits to use in the definition of the contrast.
#'
#' @return A single number, the contrast of the generated precision matrix.
#'
#' @keywords internal
MaxContrast <- function(u, omega, digits = 3) {
  diag(omega) <- diag(omega) + u
  return(Contrast(stats::cov2cor(solve(omega)), digits = digits))
}


#' Tuning function (covariance)
#'
#' Computes the difference in absolute value between the desired and observed
#' proportion of explained variance from the first Principal Component of a
#' Principal Component Analysis applied on the covariance matrix. The precision
#' matrix is obtained by adding u to the diagonal of a positive semidefinite
#' matrix. This function is used to find the value of the constant u that
#' generates a covariance matrix with desired proportion of explained variance.
#'
#' @inheritParams MaxContrast
#' @param ev_xx desired proportion of explained variance. If \code{ev_xx=NULL},
#'   the obtained proportion of explained variance is returned.
#' @param lambda eigenvalues of the positive semidefinite precision matrix.
#'
#' @return The absolute difference in proportion of explained variance (if
#'   \code{ev_xx} is provided) or observed proportion of explained variance (if
#'   \code{ev_xx=NULL}).
#'
#' @keywords internal
TuneExplainedVarianceCov <- function(u, ev_xx = NULL, lambda) {
  lambda <- lambda + u
  lambda_inv <- 1 / lambda
  tmp_ev <- max(lambda_inv) / sum(lambda_inv)
  if (is.null(ev_xx)) {
    out <- tmp_ev
  } else {
    out <- abs(tmp_ev - ev_xx)
  }
  return(out)
}


#' Tuning function (correlation)
#'
#' Computes the difference in absolute value between the desired and observed
#' proportion of explained variance from the first Principal Component of a
#' Principal Component Analysis applied on the correlation matrix. The precision
#' matrix is obtained by adding u to the diagonal of a positive semidefinite
#' matrix. This function is used to find the value of the constant u that
#' generates a correlation matrix with desired proportion of explained variance.
#'
#' @inheritParams TuneExplainedVarianceCov
#' @param omega positive semidefinite precision matrix.
#'
#' @return The absolute difference in proportion of explained variance (if
#'   \code{ev_xx} is provided) or observed proportion of explained variance (if
#'   \code{ev_xx=NULL}).
#'
#' @keywords internal
TuneExplainedVarianceCor <- function(u, ev_xx = NULL, omega) {
  diag(omega) <- diag(omega) + u
  mycor <- stats::cov2cor(solve(omega))
  tmp_ev <- norm(mycor, type = "2") / ncol(mycor)
  if (is.null(ev_xx)) {
    out <- tmp_ev
  } else {
    out <- abs(tmp_ev - ev_xx)
  }
  return(out)
}


#' Tuning function (logistic regression)
#'
#' Computes the difference in absolute value between the desired and expected C
#' statistic given the vector of true probabilities.
#'
#' @param scaling_factor constant by which log-odds (or beta coefficients when
#'   there is no intercept) are multiplied.
#' @param crude_log_odds vector of log-odds to be multiplied by the
#'   \code{scaling_factor}.
#' @param auc desired concordance statistic (area under the ROC curve). If
#'   \code{auc=NULL}, the obtained concordance statistic is returned.
#'
#' @return The absolute difference between desired and expected concordance (if
#'   \code{auc} is provided) or expected concordance (if \code{auc=NULL}).
#'
#' @keywords internal
TuneCStatisticLogit <- function(scaling_factor, crude_log_odds, auc = NULL) {
  log_odds <- crude_log_odds * scaling_factor
  proba <- 1 / (1 + exp(-log_odds))
  if (is.null(auc)) {
    out <- ExpectedConcordance(probabilities = proba)
  } else {
    tmpauc <- ExpectedConcordance(probabilities = proba)
    out <- abs(auc - tmpauc)
  }
  return(out)
}
