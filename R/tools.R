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


#' Within-group probabilities for communities
#'
#' Computes the smallest within-group probabilities that can be used to simulate
#' a graph where communities can be expected for given probabilities of
#' between-group probabilities and group sizes.
#'
#' @inheritParams SimulateGraphical
#' @param nu_mat matrix of probabilities of having an edge between nodes
#'   belonging to a given pair of node groups defined in \code{pk}. Only
#'   off-diagonal entries are used.
#'
#' @details The vector of within-group probabilities is the smallest one that
#'   can be used to generate an expected total within degree \eqn{D^w_k}
#'   strictly higher than the expected total between degree \eqn{D^b_k} for all
#'   communities \eqn{k} (see \code{\link{ExpectedCommunities}}). Namely, using
#'   the suggested within-group probabilities would give expected \eqn{D^w_k =
#'   D^b_k + 1}.
#'
#' @return A vector of within-group probabilities.
#'
#' @seealso \code{\link{ExpectedCommunities}}, \code{\link{SimulateAdjacency}},
#'   \code{\link{SimulateGraphical}}
#'
#' @examples
#' # Simulation parameters
#' pk <- rep(20, 4)
#' nu_between <- 0.1
#'
#' # Estimating smallest nu_within
#' nu_within <- MinWithinProba(pk = pk, nu_between = nu_between)
#'
#' # Expected metrics
#' ExpectedCommunities(
#'   pk = pk,
#'   nu_within = max(nu_within),
#'   nu_between = nu_between
#' )
#'
#' @export
MinWithinProba <- function(pk, nu_between = 0, nu_mat = NULL) {
  if (is.null(nu_mat)) {
    nu_mat <- diag(length(pk))
    nu_mat[upper.tri(nu_mat)] <- nu_between
    nu_mat[lower.tri(nu_mat)] <- nu_between
    diag(nu_mat) <- NA
  } else {
    if ((ncol(nu_mat) != length(pk)) & (nrow(nu_mat) != length(pk))) {
      stop("Arguments 'pk' and 'nu_mat' are not compatible. They correspond to different numbers of communities. The number of rows and columns in 'nu_mat' must be equal to the length of the vector 'pk'.")
    }
  }

  nu_within_min <- rep(NA, length(pk))
  for (k in 1:length(pk)) {
    nu_within_min[k] <- 1 / (pk[k] - 1) * (sum(nu_mat[k, -k] * pk[-k]) + 1 / pk[k])
  }

  return(nu_within_min)
}


#' Expected community structure
#'
#' Computes expected metrics related to the community structure of a graph
#' simulated with given parameters.
#'
#' @inheritParams SimulateGraphical
#'
#' @details Given a group of nodes, the within degree \eqn{d^w_i} of node
#'   \eqn{i} is defined as the number of nodes from the same group node \eqn{i}
#'   is connected to. The between degree \eqn{d^b_i} is the number of nodes from
#'   other groups node \eqn{i} is connected to. A weak community in the network
#'   is defined as a group of nodes for which the total within degree (sum of
#'   the \eqn{d^w_i} for all nodes in the community) is stricly greater than the
#'   total between degree (sum of \eqn{d^b_i} for all nodes in the community).
#'   For more details, see
#'   \href{http://networksciencebook.com/chapter/9#basics}{Network Science} by
#'   Albert-Laszlo Barabasi.
#'
#'   The expected total within and between degrees for the groups defined in
#'   \code{pk} in a network simulated using \code{SimulateAdjacency} can be
#'   computed given the group sizes (stored in \code{pk}) and probabilities of
#'   having an edge between nodes from a given group pair (defined by
#'   \code{nu_within} and \code{nu_between} or by \code{nu_mat}). The expected
#'   presence of weak communities can be inferred from these quantities.
#'
#'   The expected modularity, measuring the difference between observed and
#'   expected number of within-community edges, is also returned. For more
#'   details on this metric, see \code{\link[igraph]{modularity}}.
#'
#' @return A list with: \item{total_within_degree_c}{total within degree by node
#'   group, i.e. sum of expected within degree over all nodes in a given group.}
#'   \item{total_between_degree}{total between degree by node group, i.e. sum of
#'   expected between degree over all nodes in a given group.}
#'   \item{weak_community}{binary indicator for a given node group to be an
#'   expected weak community.} \item{total_number_edges_c}{matrix of expected
#'   number of edges between nodes from a given node pair.}
#'   \item{modularity}{expected modularity (see
#'   \code{\link[igraph]{modularity}}).}
#'
#' @seealso \code{\link{SimulateGraphical}}, \code{\link{SimulateAdjacency}},
#'   \code{\link{MinWithinProba}}
#'
#' @examples
#' # Simulation parameters
#' pk <- rep(20, 4)
#' nu_within <- 0.8
#' nu_between <- 0.1
#'
#' # Expected metrics
#' expected <- ExpectedCommunities(
#'   pk = pk,
#'   nu_within = nu_within,
#'   nu_between = nu_between
#' )
#'
#' # Example of simulated graph
#' set.seed(1)
#' theta <- SimulateAdjacency(
#'   pk = pk,
#'   nu_within = nu_within,
#'   nu_between = nu_between
#' )
#'
#' # Comparing observed and expected numbers of edges
#' bigblocks <- BlockMatrix(pk)
#' BlockStructure(pk)
#' sum(theta[which(bigblocks == 2)]) / 2
#' expected$total_number_edges_c[1, 2]
#'
#' # Comparing observed and expected modularity
#' igraph::modularity(igraph::graph_from_adjacency_matrix(theta, mode = "undirected"),
#'   membership = rep.int(1:length(pk), times = pk)
#' )
#' expected$modularity
#'
#' @export
ExpectedCommunities <- function(pk, nu_within = 0.1, nu_between = 0, nu_mat = NULL) {
  if (is.null(nu_mat)) {
    nu_mat <- diag(length(pk))
    nu_mat <- nu_mat * nu_within
    nu_mat[upper.tri(nu_mat)] <- nu_between
    nu_mat[lower.tri(nu_mat)] <- nu_between
  } else {
    if ((ncol(nu_mat) != length(pk)) & (nrow(nu_mat) != length(pk))) {
      stop("Arguments 'pk' and 'nu_mat' are not compatible. They correspond to different numbers of communities. The number of rows and columns in 'nu_mat' must be equal to the length of the vector 'pk'.")
    }
  }

  # Calculating expected sums of within and between degrees for community k
  d_within <- d_between <- rep(NA, length(pk))
  for (k in 1:length(pk)) {
    d_within[k] <- pk[k] * nu_mat[k, k] * (pk[k] - 1) # sum of within degrees in community k
    d_between[k] <- pk[k] * sum(nu_mat[k, -k] * pk[-k]) # sum of between degrees in community k
  }
  weak_community <- ifelse(d_within > d_between, yes = 1, no = 0)

  # Calculating expected number of edges within and between communities
  mystructure <- BlockStructure(pk)
  emat <- matrix(NA, nrow = nrow(mystructure), ncol = ncol(mystructure))
  for (k in unique(mystructure[upper.tri(mystructure, diag = TRUE)])) {
    if (k %in% diag(mystructure)) {
      i <- which(diag(mystructure) == k)
      emat[which(mystructure == k)] <- nu_mat[i, i] * pk[i] * (pk[i] - 1) / 2
    } else {
      tmpids <- which(mystructure == k, arr.ind = TRUE)
      i <- tmpids[1]
      j <- tmpids[2]
      emat[which(mystructure == k)] <- nu_mat[i, j] * pk[i] * pk[j]
    }
  }
  L <- sum(emat[upper.tri(emat, diag = TRUE)]) # total number of edges

  # Calculating community modularity
  M_C <- rep(NA, length(pk))
  for (k in 1:length(pk)) {
    M_C[k] <- emat[k, k] / L - ((d_within[k] + d_between[k]) / (2 * L))^2 # modularity of community k
  }

  # Calculating overall modularity
  M <- sum(M_C)

  return(list(
    total_within_degree_c = d_within,
    total_between_degree_c = d_between,
    weak_community = weak_community,
    total_number_edges_c = emat,
    modularity = M
  ))
}


#' Expected concordance statistic
#'
#' Computes the expected concordance statistic given true probabilities of
#' event. In logistic regression, the concordance statistic is equal to the area
#' under the Receiver Operating Characteristic (ROC) curve and estimates the
#' probability that an individual who experienced the event (\eqn{Y_i=1}) had a
#' higher probability of event than an individual who did not experience the
#' event (\eqn{Y_i=0}).
#'
#' @param probabilities vector of probabilities of event.
#'
#' @return The expected concordance statistic.
#'
#' @seealso \code{\link{Concordance}}
#'
#' @examples
#'
#' # Simulation of probabilities
#' set.seed(1)
#' proba <- runif(n = 1000)
#'
#' # Expected concordance
#' ExpectedConcordance(proba)
#'
#' # Simulation of binary outcome
#' ydata <- rbinom(n = length(proba), size = 1, prob = proba)
#'
#' # Observed concordance in simulated data
#' Concordance(observed = ydata, predicted = proba)
#'
#' @export
ExpectedConcordance <- function(probabilities) {
  expected_number_pairs <- 0
  expected_concordant_pairs <- 0
  for (i in 1:(length(probabilities) - 1)) {
    for (j in (i + 1):length(probabilities)) {
      tmp1 <- probabilities[i] * (1 - probabilities[j])
      tmp2 <- probabilities[j] * (1 - probabilities[i])
      if (probabilities[i] > probabilities[j]) {
        expected_concordant_pairs <- expected_concordant_pairs + tmp1
      }
      if (probabilities[j] > probabilities[i]) {
        expected_concordant_pairs <- expected_concordant_pairs + tmp2
      }
      # Ties are not concordant pairs
      expected_number_pairs <- expected_number_pairs + tmp1 + tmp2
    }
  }
  return(expected_concordant_pairs / expected_number_pairs)
}


#' Layered Directed Acyclic Graph
#'
#' Returns the adjacency matrix of a layered Directed Acyclic Graph. In this
#' graph, arrows go from all members of a layer to all members of the following
#' layers. There are no arrows between members of the same layer.
#'
#' @param layers list of vectors. Each vector in the list corresponds to a
#'   layer. There are as many layers as items in the list. Alternatively, this
#'   argument can be a vector of the number of variables per layer.
#' @param n_manifest vector of the number of manifest (observed) variables
#'   measuring each of the latent variables. If \code{n_manifest} is provided,
#'   the variables defined in argument \code{layers} are considered latent. All
#'   entries of \code{n_manifest} must be strictly positive.
#'
#' @return The adjacency matrix of the layered Directed Acyclic Graph.
#'
#' @examples
#' # Example with 3 layers specified in a list
#' layers <- list(
#'   c("x1", "x2", "x3"),
#'   c("x4", "x5"),
#'   c("x6", "x7", "x8")
#' )
#' dag <- LayeredDAG(layers)
#' plot(dag)
#'
#' # Example with 3 layers specified in a vector
#' dag <- LayeredDAG(layers = c(3, 2, 3))
#' plot(dag)
#'
#' @export
LayeredDAG <- function(layers, n_manifest = NULL) {
  # Extracting the number of members per layer
  if (is.list(layers)) {
    pk <- sapply(layers, length)
  } else {
    pk <- layers
  }

  # Creating the adjacency matrix
  adjacency <- BlockMatrix(pk)
  adjacency[which(adjacency %in% diag(BlockStructure(pk)))] <- 0
  adjacency <- ifelse(adjacency != 0, yes = 1, no = 0)
  adjacency[lower.tri(adjacency)] <- 0

  # Addition of manifest variables for each latent variable
  if (!is.null(n_manifest)) {
    # Expanding the vector if needed
    if (length(n_manifest) != ncol(adjacency)) {
      n_manifest <- rep(n_manifest[1], ncol(adjacency))
    }

    # Adding manifest variables in adjacency matrix
    tmpfactor <- as.factor(rep.int(1:length(n_manifest), times = n_manifest))
    submatrix_manifest <- t(stats::model.matrix(~ tmpfactor - 1))
    adjacency <- cbind(submatrix_manifest, adjacency)
    adjacency <- rbind(matrix(0, ncol = ncol(adjacency), nrow = ncol(adjacency) - nrow(adjacency)), adjacency)
  }

  # Defining row and column names
  if (!is.null(n_manifest)) {
    ids_manifest <- seq(1, sum(n_manifest))
    ids_latent <- seq(sum(n_manifest) + 1, ncol(adjacency))
    rownames(adjacency)[ids_manifest] <- colnames(adjacency)[ids_manifest] <- paste0("x", seq(1, length(ids_manifest)))
    rownames(adjacency)[ids_latent] <- colnames(adjacency)[ids_latent] <- paste0("f", seq(1, length(ids_latent)))
  } else {
    if (is.list(layers)) {
      colnames(adjacency) <- rownames(adjacency) <- unlist(layers)
    } else {
      colnames(adjacency) <- rownames(adjacency) <- paste0("x", 1:sum(pk))
    }
  }

  # Defining the class of the output
  class(adjacency) <- c("matrix", "adjacency_matrix")

  return(adjacency)
}
