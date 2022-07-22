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


#' Simulation of precision matrix
#'
#' Simulates a sparse precision matrix from a binary adjacency matrix
#' \code{theta} encoding conditional independence in a Gaussian Graphical Model.
#'
#' @inheritParams SimulateGraphical
#' @param theta binary and symmetric adjacency matrix encoding the conditional
#'   independence structure.
#' @param scale logical indicating if the proportion of explained variance by
#'   PC1 should be computed from the correlation (\code{scale=TRUE}) or
#'   covariance (\code{scale=FALSE}) matrix.
#'
#' @return A list with: \item{omega}{true simulated precision matrix.}
#'   \item{u}{value of the constant u used to ensure that \code{omega} is
#'   positive definite.}
#'
#' @seealso \code{\link{SimulateGraphical}}, \code{\link{MakePositiveDefinite}}
#'
#' @details Entries that are equal to zero in the adjacency matrix \code{theta}
#'   are also equal to zero in the generated precision matrix. These zero
#'   entries indicate conditional independence between the corresponding pair of
#'   variables (see \code{\link{SimulateGraphical}}).
#'
#'   Argument \code{pk} can be specified to create groups of variables and allow
#'   for nonzero precision entries to be sampled from different distributions
#'   between two variables belonging to the same group or to different groups.
#'
#'   If \code{continuous=FALSE}, nonzero off-diagonal entries of the precision
#'   matrix are sampled from a discrete uniform distribution taking values in
#'   \code{v_within} (for entries in the diagonal block) or \code{v_between}
#'   (for entries in off-diagonal blocks). If \code{continuous=TRUE}, nonzero
#'   off-diagonal entries are sampled from a continuous uniform distribution
#'   taking values in the range given by \code{v_within} or \code{v_between}.
#'
#'   Diagonal entries of the precision matrix are defined to ensure positive
#'   definiteness as described in \code{\link{MakePositiveDefinite}}.
#'
#' @references \insertRef{ourstabilityselection}{fake}
#'
#' @examples
#' # Simulation of an adjacency matrix
#' theta <- SimulateAdjacency(pk = c(5, 5), nu_within = 0.7)
#' print(theta)
#'
#' # Simulation of a precision matrix maximising the contrast
#' simul <- SimulatePrecision(theta = theta)
#' print(simul$omega)
#'
#' # Simulation of a precision matrix with specific ev by PC1
#' simul <- SimulatePrecision(
#'   theta = theta,
#'   pd_strategy = "min_eigenvalue",
#'   ev_xx = 0.3, scale = TRUE
#' )
#' print(simul$omega)
#' @export
SimulatePrecision <- function(pk = NULL, theta,
                              v_within = c(0.5, 1), v_between = c(0, 0.1),
                              v_sign = c(-1, 1), continuous = TRUE,
                              pd_strategy = "diagonally_dominant", ev_xx = NULL, scale = TRUE,
                              u_list = c(1e-10, 1), tol = .Machine$double.eps^0.25) {
  # Checking inputs and defining pk
  if (is.null(pk)) {
    pk <- ncol(theta)
  } else {
    if (sum(pk) != ncol(theta)) {
      stop("Arguments 'pk' and 'theta' are not consistent. The sum of 'pk' entries must be equal to the number of rows and columns in 'theta'.")
    }
  }

  # Checking the choice of pd_strategy
  if (!pd_strategy %in% c("diagonally_dominant", "min_eigenvalue")) {
    stop("Invalid input for argument 'pd_strategy'. Possible values are: 'diagonally_dominant' or 'min_eigenvalue'.")
  }

  # Checking other input values
  if (any((v_within < 0) | (v_within > 1))) {
    stop("Invalid input for argument 'v_within'. Values must be between 0 and 1.")
  }
  if (any((v_between < 0) | (v_between > 1))) {
    stop("Invalid input for argument 'v_between'. Values must be between 0 and 1.")
  }
  if (any(!v_sign %in% c(-1, 1))) {
    stop("Invalid input for argument 'v_sign'. Possible values are -1 and 1.")
  }

  # Ensuring that v values are lower than or equal to 1
  if (any(abs(v_within) > 1)) {
    v_within <- v_within / max(abs(v_within))
    message("The values provided in 'v_within' have been re-scaled to be lower than or equal to 1 in absolute value.")
  }

  # Ensuring that diagonal entries of theta are zero
  diag(theta) <- 0

  # Building v matrix
  v <- SimulateSymmetricMatrix(
    pk = pk, v_within = v_within, v_between = v_between,
    v_sign = v_sign, continuous = continuous
  )

  # Filling off-diagonal entries of the precision matrix
  omega_tilde <- theta * v

  # Ensuring positive definiteness
  omega_pd <- MakePositiveDefinite(
    omega = omega_tilde, pd_strategy = pd_strategy,
    ev_xx = ev_xx, scale = scale, u_list = u_list, tol = tol
  )

  # Returning the output
  return(omega_pd)
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
#' matrix. This function is used to find the value of the constant u
#' that generates a covariance matrix with desired proportion of explained
#' variance.
#'
#' @inheritParams MaxContrast
#' @param ev_xx desired proportion of explained variance. If \code{ev_xx=NULL}, the
#'   obtained proportion of explained variance is returned.
#' @param lambda eigenvalues of the positive semidefinite precision matrix.
#'
#' @return The difference in proportion of explained variance in absolute values
#'   or observed proportion of explained variance (if \code{ev_xx=NULL}).
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
#' matrix. This function is used to find the value of the constant u
#' that generates a correlation matrix with desired proportion of explained
#' variance.
#'
#' @inheritParams TuneExplainedVarianceCov
#' @param omega positive semidefinite precision matrix.
#'
#' @return The difference in proportion of explained variance in absolute values
#'   or observed proportion of explained variance (if \code{ev_xx=NULL}).
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


#' Tuning function (regression)
#'
#' Computes the absolute difference between the smallest eigenvalue and
#' requested one (parameter \code{tol}) for a precision matrix with predictors
#' and outcomes.
#'
#' @param ev_xz proportion of explained variance.
#' @param omega precision matrix.
#' @param tol requested smallest eigenvalue after transformation of the input
#'   precision matrix.
#' @param q number of outcome variables.
#' @param p number of predictor variables.
#'
#' @return The absolute difference between the smallest eigenvalue of the
#'   transformed precision matrix and requested value \code{tol}.
#'
#' @keywords internal
TuneExplainedVarianceReg <- function(ev_xz, omega, tol = 0.1, q, p) {
  ev_xz <- rep(ev_xz, q)
  for (j in 1:q) {
    pred_ids <- seq(q + 1, q + p)
    omega[j, j] <- omega[j, pred_ids, drop = FALSE] %*% solve(omega[pred_ids, pred_ids]) %*% t(omega[j, pred_ids, drop = FALSE]) * 1 / ev_xz[j]
  }
  return(abs(min(eigen(omega, only.values = TRUE)$values) - tol))
}
