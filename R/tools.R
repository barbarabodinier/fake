#' Making positive definite matrix
#'
#' Determines the diagonal entries of a symmetric matrix to make it is positive
#' definite.
#'
#' @inheritParams SimulateGraphical
#' @param omega input matrix.
#' @param scale logical indicating if the proportion of explained variance by
#'   PC1 should be computed from the correlation (\code{scale=TRUE}) or
#'   covariance (\code{scale=FALSE}) matrix.
#'
#' @return A list with: \item{omega}{positive definite matrix.} \item{u}{value
#'   of the constant u.}
#'
#' @details Two strategies are implemented to ensure positive definiteness: by
#'   diagonally dominance or using eigendecomposition.
#'
#'   A diagonally dominant symmetric matrix with positive diagonal entries is
#'   positive definite. With \code{pd_strategy="diagonally_dominant"}, the
#'   diagonal entries of the matrix are defined to be strictly higher than the
#'   sum of entries on the corresponding row in absolute value, which ensures
#'   diagonally dominance. Let \eqn{\Omega*} denote the input matrix with zeros
#'   on the diagonal and \eqn{\Omega} be the output positive definite matrix. We
#'   have:
#'
#'   \eqn{\Omega_{ii} = \sum_{j = 1}^p | \Omega_{ij}* | + u}, where \eqn{u > 0}
#'   is a parameter.
#'
#'   A matrix is positive definite if all its eigenvalues are positive. With
#'   \code{pd_strategy="diagonally_dominant"}, diagonal entries of the matrix
#'   are defined to be higher than the absolute value of the smallest eigenvalue
#'   of the same matrix with a diagonal of zeros. Let \eqn{\lambda_1} denote the
#'   smallest eigenvvalue of the input matrix \eqn{\Omega*} with a diagonal of
#'   zeros, and \eqn{v_1} be the corresponding eigenvector. Diagonal entries in
#'   the output matrix \eqn{\Omega} are defined as:
#'
#'   \eqn{\Omega_{ii} = | \lambda_1 | + u}, where \eqn{u > 0} is a parameter.
#'
#'   It can be showed that \eqn{\Omega} has stricly positive eigenvalues. Let
#'   \eqn{\lambda} and \eqn{v} denote any eigenpair of \eqn{\Omega*}:
#'
#'   \eqn{\Omega* v = \lambda v}
#'
#'   \eqn{\Omega* v + (| \lambda_1 | + u) v = \lambda v + (| \lambda_1 | + u) v}
#'
#'   \eqn{(\Omega* + (| \lambda_1 | + u) I) v = (\lambda + | \lambda_1 | + u) v}
#'
#'   \eqn{\Omega v = (\lambda + | \lambda_1 | + u) v}
#'
#'   The eigenvalues of \eqn{\Omega} are equal to the eigenvalues of
#'   \eqn{\Omega*} plus \eqn{| \lambda_1 |}. The smallest eigenvalue of
#'   \eqn{\Omega} is \eqn{(\lambda_1 + | \lambda_1 | + u) > 0}.
#'
#'
#'   Considering the matrix to make positive definite is a precision matrix, its
#'   standardised inverse matrix is the correlation matrix. In both cases, the
#'   magnitude of correlations is controlled by the constant u.
#'
#'   If \code{ev_xx=NULL}, the constant u is chosen to maximise the
#'   \code{\link{Contrast}} of the corresponding correlation matrix.
#'
#'   If \code{ev_xx} is provided, the constant u is chosen to generate a
#'   correlation matrix with required proportion of explained variance by the
#'   first Principal Component, if possible. This proportion of explained
#'   variance is equal to the largest eigenvalue of the correlation matrix
#'   divided by the sum of its eigenvalues. If \code{scale=FALSE}, the
#'   covariance matrix is used instead of the correlation matrix for faster
#'   computations.
#'
#' @references \insertRef{ourstabilityselection}{fake}
#'
#' @examples
#' # Simulation of a symmetric matrix
#' p <- 5
#' set.seed(1)
#' omega <- matrix(rnorm(p * p), ncol = p)
#' omega <- omega + t(omega)
#' diag(omega) <- 0
#'
#' # Diagonal dominance maximising contrast
#' omega_pd <- MakePositiveDefinite(omega,
#'   pd_strategy = "diagonally_dominant"
#' )
#' eigen(omega_pd$omega)$values # positive eigenvalues
#'
#' # Diagonal dominance with specific proportion of explained variance by PC1
#' omega_pd <- MakePositiveDefinite(omega,
#'   pd_strategy = "diagonally_dominant",
#'   ev_xx = 0.55
#' )
#' lambda_inv <- eigen(cov2cor(solve(omega_pd$omega)))$values
#' max(lambda_inv) / sum(lambda_inv) # expected ev
#'
#' # Version not scaled (using eigenvalues from the covariance)
#' omega_pd <- MakePositiveDefinite(omega,
#'   pd_strategy = "diagonally_dominant",
#'   ev_xx = 0.55, scale = FALSE
#' )
#' lambda_inv <- 1 / eigen(omega_pd$omega)$values
#' max(lambda_inv) / sum(lambda_inv) # expected ev
#'
#' # Non-negative eigenvalues maximising contrast
#' omega_pd <- MakePositiveDefinite(omega,
#'   pd_strategy = "min_eigenvalue"
#' )
#' eigen(omega_pd$omega)$values # positive eigenvalues
#'
#' # Non-negative eigenvalues with specific proportion of explained variance by PC1
#' omega_pd <- MakePositiveDefinite(omega,
#'   pd_strategy = "min_eigenvalue",
#'   ev_xx = 0.7
#' )
#' lambda_inv <- eigen(cov2cor(solve(omega_pd$omega)))$values
#' max(lambda_inv) / sum(lambda_inv)
#'
#' # Version not scaled (using eigenvalues from the covariance)
#' omega_pd <- MakePositiveDefinite(omega,
#'   pd_strategy = "min_eigenvalue",
#'   ev_xx = 0.7, scale = FALSE
#' )
#' lambda_inv <- 1 / eigen(omega_pd$omega)$values
#' max(lambda_inv) / sum(lambda_inv)
#' @export
MakePositiveDefinite <- function(omega, pd_strategy = "diagonally_dominant",
                                 ev_xx = NULL, scale = TRUE, u_list = c(1e-10, 1),
                                 tol = .Machine$double.eps^0.25) {

  # Making positive definite using diagonally dominance
  if (pd_strategy == "diagonally_dominant") {
    # Constructing the diagonal as the sum of entries
    diag(omega) <- apply(abs(omega), 1, sum)
    lambda <- eigen(omega, only.values = TRUE)$values
  }
  # Making positive definite using eigendecomposition
  if (pd_strategy == "min_eigenvalue") {
    # Extracting smallest eigenvalue of omega_tilde
    lambda <- eigen(omega, only.values = TRUE)$values
    lambda0 <- abs(min(lambda))

    # Making the precision matrix positive semidefinite
    lambda <- lambda + lambda0
    diag(omega) <- lambda0
  }

  if (is.null(ev_xx)) {
    # Finding u that maximises the contrast
    if (min(u_list) != max(u_list)) {
      argmax_u <- stats::optimise(MaxContrast,
        omega = omega, maximum = TRUE, tol = tol,
        lower = min(u_list), upper = max(u_list)
      )
      u <- argmax_u$maximum
    } else {
      u <- min(u_list)
    }
  } else {
    # Finding extreme values
    if (scale) {
      max_ev <- TuneExplainedVarianceCor(u = min(u_list), omega = omega)
      min_ev <- TuneExplainedVarianceCor(u = max(u_list), omega = omega)
    } else {
      max_ev <- TuneExplainedVarianceCov(u = min(u_list), lambda = lambda)
      min_ev <- TuneExplainedVarianceCov(u = max(u_list), lambda = lambda)
    }

    # Finding u corresponding to the required proportion of explained variance
    if ((ev_xx <= min_ev) | (ev_xx >= max_ev)) {
      if (ev_xx <= min_ev) {
        u <- max(u_list)
        if (ev_xx < min_ev) {
          message(paste0("The smallest proportion of explained variance by PC1 that can be obtained is ", round(min_ev, digits = 2), "."))
        }
      } else {
        u <- min(u_list)
        if (ev_xx > max_ev) {
          message(paste0("The largest proportion of explained variance by PC1 that can be obtained is ", round(max_ev, digits = 2), "."))
        }
      }
    } else {
      if (scale) {
        # Minimising the difference between requested and possible ev
        if (min(u_list) != max(u_list)) {
          argmin_u <- stats::optimise(TuneExplainedVarianceCor,
            omega = omega, ev_xx = ev_xx,
            lower = min(u_list), upper = max(u_list), tol = tol
          )
          u <- argmin_u$minimum
        } else {
          u <- min(u_list)
        }
      } else {
        # Minimising the difference between requested and possible ev
        if (min(u_list) != max(u_list)) {
          argmin_u <- stats::optimise(TuneExplainedVarianceCov,
            lambda = lambda, ev_xx = ev_xx,
            lower = min(u_list), upper = max(u_list), tol = tol
          )
          u <- argmin_u$minimum
        } else {
          u <- min(u_list)
        }
      }
    }
  }

  # Constructing the diagonal
  diag(omega) <- diag(omega) + u

  return(list(omega = omega, u = u))
}


#' Matrix contrast
#'
#' Computes matrix contrast, defined as the number of unique truncated
#' entries with a specified number of digits.
#'
#' @param mat input matrix.
#' @param digits number of digits to use.
#'
#' @return A single number, the contrast of the input matrix.
#'
#' @references \insertRef{ourstabilityselection}{fake}
#'
#' @examples
#' # Example 1
#' mat <- matrix(c(0.1, 0.2, 0.2, 0.2), ncol = 2, byrow = TRUE)
#' Contrast(mat)
#'
#' # Example 2
#' mat <- matrix(c(0.1, 0.2, 0.2, 0.3), ncol = 2, byrow = TRUE)
#' Contrast(mat)
#'
#' @export
Contrast <- function(mat, digits = 3) {
  return(length(unique(round(as.vector(abs(mat)), digits = digits))))
}
