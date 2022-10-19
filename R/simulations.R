#' Data simulation for Gaussian Graphical Modelling
#'
#' Simulates data from a Gaussian Graphical Model (GGM).
#'
#' @param n number of observations in the simulated data.
#' @param pk vector of the number of variables per group in the simulated data.
#'   The number of nodes in the simulated graph is \code{sum(pk)}. With multiple
#'   groups, the simulated (partial) correlation matrix has a block structure,
#'   where blocks arise from the integration of the \code{length(pk)} groups.
#'   This argument is only used if \code{theta} is not provided.
#' @param theta optional binary and symmetric adjacency matrix encoding the
#'   conditional independence structure.
#' @param implementation function for simulation of the graph. By default,
#'   algorithms implemented in \code{\link[huge]{huge.generator}} are used.
#'   Alternatively, a user-defined function can be used. It must take \code{pk},
#'   \code{topology} and \code{nu} as arguments and return a
#'   \code{(sum(pk)*(sum(pk)))} binary and symmetric matrix for which diagonal
#'   entries are all equal to zero. This function is only applied if
#'   \code{theta} is not provided.
#' @param topology topology of the simulated graph. If using
#'   \code{implementation=HugeAdjacency}, possible values are listed for the
#'   argument \code{graph} of \code{\link[huge]{huge.generator}}. These are:
#'   "random", "hub", "cluster", "band" and "scale-free".
#' @param nu_within probability of having an edge between two nodes belonging to
#'   the same group, as defined in \code{pk}. If \code{length(pk)=1}, this is
#'   the expected density of the graph. If \code{implementation=HugeAdjacency},
#'   this argument is only used for \code{topology="random"} or
#'   \code{topology="cluster"} (see argument \code{prob} in
#'   \code{\link[huge]{huge.generator}}). Only used if \code{nu_mat} is not
#'   provided.
#' @param nu_between probability of having an edge between two nodes belonging
#'   to different groups, as defined in \code{pk}. By default, the same density
#'   is used for within and between blocks (\code{nu_within}=\code{nu_between}).
#'   Only used if \code{length(pk)>1}. Only used if \code{nu_mat} is not
#'   provided.
#' @param nu_mat matrix of probabilities of having an edge between nodes
#'   belonging to a given pair of node groups defined in \code{pk}.
#' @param output_matrices logical indicating if the true precision and (partial)
#'   correlation matrices should be included in the output.
#' @param v_within vector defining the (range of) nonzero entries in the
#'   diagonal blocks of the precision matrix. These values must be between -1
#'   and 1 if \code{pd_strategy="min_eigenvalue"}. If \code{continuous=FALSE},
#'   \code{v_within} is the set of possible precision values. If
#'   \code{continuous=TRUE}, \code{v_within} is the range of possible precision
#'   values.
#' @param v_between vector defining the (range of) nonzero entries in the
#'   off-diagonal blocks of the precision matrix. This argument is the same as
#'   \code{v_within} but for off-diagonal blocks. It is only used if
#'   \code{length(pk)>1}.
#' @param v_sign vector of possible signs for precision matrix entries. Possible
#'   inputs are: \code{-1} for positive partial correlations, \code{1} for
#'   negative partial correlations, or \code{c(-1, 1)} for both positive and
#'   negative partial correlations.
#' @param continuous logical indicating whether to sample precision values from
#'   a uniform distribution between the minimum and maximum values in
#'   \code{v_within} (diagonal blocks) or \code{v_between} (off-diagonal blocks)
#'   (if \code{continuous=TRUE}) or from proposed values in \code{v_within}
#'   (diagonal blocks) or \code{v_between} (off-diagonal blocks) (if
#'   \code{continuous=FALSE}).
#' @param pd_strategy method to ensure that the generated precision matrix is
#'   positive definite (and hence can be a covariance matrix). If
#'   \code{pd_strategy="diagonally_dominant"}, the precision matrix is made
#'   diagonally dominant by setting the diagonal entries to the sum of absolute
#'   values on the corresponding row and a constant u. If
#'   \code{pd_strategy="min_eigenvalue"}, diagonal entries are set to the sum of
#'   the absolute value of the smallest eigenvalue of the precision matrix with
#'   zeros on the diagonal and a constant u.
#' @param ev_xx expected proportion of explained variance by the first Principal
#'   Component (PC1) of a Principal Component Analysis. This is the largest
#'   eigenvalue of the correlation (if \code{scale_ev=TRUE}) or covariance (if
#'   \code{scale_ev=FALSE}) matrix divided by the sum of eigenvalues. If
#'   \code{ev_xx=NULL} (the default), the constant u is chosen by maximising the
#'   contrast of the correlation matrix.
#' @param scale_ev logical indicating if the proportion of explained variance by
#'   PC1 should be computed from the correlation (\code{scale_ev=TRUE}) or
#'   covariance (\code{scale_ev=FALSE}) matrix. If \code{scale_ev=TRUE}, the
#'   correlation matrix is used as parameter of the multivariate normal
#'   distribution.
#' @param u_list vector with two numeric values defining the range of values to
#'   explore for constant u.
#' @param tol accuracy for the search of parameter u as defined in
#'   \code{\link[stats]{optimise}}.
#' @param scale logical indicating if the true mean is zero and true variance is
#'   one for all simulated variables. The observed mean and variance may be
#'   slightly off by chance.
#' @param ... additional arguments passed to the graph simulation function
#'   provided in \code{implementation}.
#'
#' @seealso \code{\link{SimulatePrecision}}, \code{\link{MakePositiveDefinite}},
#'   \code{\link{Contrast}}
#'
#' @family simulation functions
#'
#' @return A list with: \item{data}{simulated data with \code{n} observation and
#'   \code{sum(pk)} variables.} \item{theta}{adjacency matrix of the simulated
#'   graph.} \item{omega}{simulated (true) precision matrix. Only returned if
#'   \code{output_matrices=TRUE}.} \item{phi}{simulated (true) partial
#'   correlation matrix. Only returned if \code{output_matrices=TRUE}.}
#'   \item{sigma}{simulated (true) covariance matrix. Only returned if
#'   \code{output_matrices=TRUE}.} \item{u}{value of the constant u used for the
#'   simulation of \code{omega}. Only returned if \code{output_matrices=TRUE}.}
#'
#' @details The simulation is done in two steps with (i) generation of a graph,
#'   and (ii) sampling from multivariate Normal distribution for which nonzero
#'   entries in the partial correlation matrix correspond to the edges of the
#'   simulated graph. This procedure ensures that the conditional independence
#'   structure between the variables corresponds to the simulated graph.
#'
#'   Step 1 is done using \code{\link{SimulateAdjacency}}.
#'
#'   In Step 2, the precision matrix (inverse of the covariance matrix) is
#'   simulated using \code{\link{SimulatePrecision}} so that (i) its nonzero
#'   entries correspond to edges in the graph simulated in Step 1, and (ii) it
#'   is positive definite (see \code{\link{MakePositiveDefinite}}). The inverse
#'   of the precision matrix is used as covariance matrix to simulate data from
#'   a multivariate Normal distribution.
#'
#'   The outputs of this function can be used to evaluate the ability of a
#'   graphical model to recover the conditional independence structure.
#'
#' @references \insertRef{ourstabilityselection}{fake}
#'
#' @examples
#' # Simulation of random graph with 50 nodes
#' set.seed(1)
#' simul <- SimulateGraphical(n = 100, pk = 50, topology = "random", nu_within = 0.05)
#' print(simul)
#' plot(simul)
#'
#' # Simulation of scale-free graph with 20 nodes
#' set.seed(1)
#' simul <- SimulateGraphical(n = 100, pk = 20, topology = "scale-free")
#' plot(simul)
#'
#' # Extracting true precision/correlation matrices
#' set.seed(1)
#' simul <- SimulateGraphical(
#'   n = 100, pk = 20,
#'   topology = "scale-free", output_matrices = TRUE
#' )
#' str(simul)
#'
#' # Simulation of multi-block data
#' set.seed(1)
#' pk <- c(20, 30)
#' simul <- SimulateGraphical(
#'   n = 100, pk = pk,
#'   pd_strategy = "min_eigenvalue"
#' )
#' mycor <- cor(simul$data)
#' Heatmap(mycor,
#'   col = c("darkblue", "white", "firebrick3"),
#'   legend_range = c(-1, 1), legend_length = 50,
#'   legend = FALSE, axes = FALSE
#' )
#' for (i in 1:2) {
#'   axis(side = i, at = c(0.5, pk[1] - 0.5), labels = NA)
#'   axis(
#'     side = i, at = mean(c(0.5, pk[1] - 0.5)),
#'     labels = ifelse(i == 1, yes = "Group 1", no = "Group 2"),
#'     tick = FALSE, cex.axis = 1.5
#'   )
#'   axis(side = i, at = c(pk[1] + 0.5, sum(pk) - 0.5), labels = NA)
#'   axis(
#'     side = i, at = mean(c(pk[1] + 0.5, sum(pk) - 0.5)),
#'     labels = ifelse(i == 1, yes = "Group 2", no = "Group 1"),
#'     tick = FALSE, cex.axis = 1.5
#'   )
#' }
#'
#' # User-defined function for graph simulation
#' CentralNode <- function(pk, hub = 1) {
#'   theta <- matrix(0, nrow = sum(pk), ncol = sum(pk))
#'   theta[hub, ] <- 1
#'   theta[, hub] <- 1
#'   diag(theta) <- 0
#'   return(theta)
#' }
#' simul <- SimulateGraphical(n = 100, pk = 10, implementation = CentralNode)
#' plot(simul) # star
#' simul <- SimulateGraphical(n = 100, pk = 10, implementation = CentralNode, hub = 2)
#' plot(simul) # variable 2 is the central node
#'
#' # User-defined adjacency matrix
#' mytheta <- matrix(c(
#'   0, 1, 1, 0,
#'   1, 0, 0, 0,
#'   1, 0, 0, 1,
#'   0, 0, 1, 0
#' ), ncol = 4, byrow = TRUE)
#' simul <- SimulateGraphical(n = 100, theta = mytheta)
#' plot(simul)
#'
#' # User-defined adjacency and block structure
#' simul <- SimulateGraphical(n = 100, theta = mytheta, pk = c(2, 2))
#' mycor <- cor(simul$data)
#' Heatmap(mycor,
#'   col = c("darkblue", "white", "firebrick3"),
#'   legend_range = c(-1, 1), legend_length = 50, legend = FALSE
#' )
#' @export
SimulateGraphical <- function(n = 100, pk = 10, theta = NULL,
                              implementation = HugeAdjacency, topology = "random",
                              nu_within = 0.1, nu_between = NULL, nu_mat = NULL,
                              v_within = c(0.5, 1), v_between = c(0.1, 0.2),
                              v_sign = c(-1, 1), continuous = TRUE,
                              pd_strategy = "diagonally_dominant", ev_xx = NULL, scale_ev = TRUE,
                              u_list = c(1e-10, 1), tol = .Machine$double.eps^0.25,
                              scale = TRUE, output_matrices = FALSE, ...) {
  # Defining number of nodes
  p <- sum(pk)
  if (!is.null(theta)) {
    if (ncol(theta) != p) {
      p <- pk <- ncol(theta)
    }
  }

  # Defining the between-block density
  if (is.null(nu_between)) {
    nu_between <- nu_within
  }

  # Building adjacency matrix
  if (is.null(theta)) {
    theta <- SimulateAdjacency(
      pk = pk,
      implementation = implementation, topology = topology,
      nu_within = nu_within, nu_between = nu_between, nu_mat = nu_mat, ...
    )
  }

  # Simulation of a precision matrix
  out <- SimulatePrecision(
    pk = pk, theta = theta,
    v_within = v_within, v_between = v_between,
    v_sign = v_sign, continuous = continuous,
    pd_strategy = pd_strategy, ev_xx = ev_xx, scale = scale_ev,
    u_list = u_list, tol = tol
  )
  omega <- out$omega

  # Computing the covariance matrix
  if (scale) {
    sigma <- stats::cov2cor(solve(omega))
  } else {
    sigma <- solve(omega)
  }

  # Computing the partial correlation matrix
  if (output_matrices) {
    phi <- -stats::cov2cor(omega) + 2 * diag(ncol(omega))
  }

  # Simulating data from multivariate normal distribution
  x <- MASS::mvrnorm(n, rep(0, p), sigma)
  colnames(x) <- paste0("var", 1:ncol(x))
  rownames(x) <- paste0("obs", 1:nrow(x))

  if (output_matrices) {
    out <- list(
      data = x, theta = theta,
      omega = omega, phi = phi, sigma = sigma,
      u = out$u
    )
  } else {
    out <- list(data = x, theta = theta)
  }

  # Defining the class
  class(out) <- "simulation_graphical_model"

  return(out)
}


#' Data simulation for sparse Principal Component Analysis
#'
#' Simulates data with with independent groups of variables.
#'
#' @inheritParams SimulateGraphical
#' @param adjacency optional binary and symmetric adjacency matrix encoding the
#'   conditional graph structure between observations. The clusters encoded in
#'   this argument must be in line with those indicated in \code{pk}. Edges in
#'   off-diagonal blocks are not allowed to ensure that the simulated orthogonal
#'   components are sparse. Corresponding entries in the precision matrix will
#'   be set to zero.
#'
#' @details The data is simulated from a centered multivariate Normal
#'   distribution with a block-diagonal covariance matrix. Independence between
#'   variables from the different blocks ensures that sparse orthogonal
#'   components can be generated.
#'
#'   The block-diagonal partial correlation matrix is obtained using a graph
#'   structure encoding the conditional independence between variables. The
#'   orthogonal latent variables are obtained from eigendecomposition of the
#'   true correlation matrix. The sparse eigenvectors contain the weights of the
#'   linear combination of variables to construct the latent variable (loadings
#'   coefficients). The proportion of explained variance by each of the latent
#'   variable is computed from eigenvalues.
#'
#'   As latent variables are defined from the true correlation matrix, the
#'   number of sparse orthogonal components is not limited by the number of
#'   observations and is equal to \code{sum(pk)}.
#'
#' @return A list with: \item{data}{simulated data with \code{n} observation and
#'   \code{sum(pk)} variables.} \item{loadings}{loadings coefficients of the
#'   orthogonal latent variables (principal components).} \item{theta}{support
#'   of the loadings coefficients.} \item{ev}{proportion of explained variance
#'   by each of the orthogonal latent variables.} \item{adjacency}{adjacency
#'   matrix of the simulated graph.} \item{omega}{simulated (true) precision
#'   matrix. Only returned if \code{output_matrices=TRUE}.} \item{phi}{simulated
#'   (true) partial correlation matrix. Only returned if
#'   \code{output_matrices=TRUE}.} \item{C}{ simulated (true) correlation
#'   matrix. Only returned if \code{output_matrices=TRUE}.}
#'
#' @seealso \code{\link{MakePositiveDefinite}}
#'
#' @family simulation functions
#'
#' @references \insertRef{ourstabilityselection}{fake}
#'
#' @examples
#' # Simulation of 3 components with high e.v.
#' set.seed(1)
#' simul <- SimulateComponents(pk = c(5, 3, 4), ev_xx = 0.4)
#' print(simul)
#' plot(simul)
#' plot(cumsum(simul$ev), ylim = c(0, 1), las = 1)
#'
#' # Simulation of 3 components with moderate e.v.
#' set.seed(1)
#' simul <- SimulateComponents(pk = c(5, 3, 4), ev_xx = 0.25)
#' print(simul)
#' plot(simul)
#' plot(cumsum(simul$ev), ylim = c(0, 1), las = 1)
#'
#' # Simulation of multiple components with low e.v.
#' pk <- sample(3:10, size = 5, replace = TRUE)
#' simul <- SimulateComponents(
#'   pk = pk,
#'   nu_within = 0.3, v_within = c(0.8, 0.5), v_sign = -1, ev_xx = 0.1
#' )
#' plot(simul)
#' plot(cumsum(simul$ev), ylim = c(0, 1), las = 1)
#' @export
SimulateComponents <- function(n = 100, pk = c(10, 10),
                               adjacency = NULL, nu_within = 1,
                               v_within = c(0.5, 1), v_sign = -1, continuous = TRUE,
                               pd_strategy = "min_eigenvalue", ev_xx = 0.1, scale_ev = TRUE,
                               u_list = c(1e-10, 1), tol = .Machine$double.eps^0.25,
                               scale = TRUE, output_matrices = FALSE) {
  # Using multi-block simulator with unconnected blocks
  out <- SimulateGraphical(
    n = n, pk = pk, theta = adjacency,
    implementation = HugeAdjacency,
    topology = "random",
    nu_within = nu_within, # fully connected components by default
    nu_between = 0, # need unconnected blocks
    v_within = v_within,
    v_between = 0,
    v_sign = v_sign,
    continuous = continuous,
    pd_strategy = pd_strategy, ev_xx = ev_xx, scale_ev = scale_ev,
    u_list = u_list, tol = tol,
    scale = scale, output_matrices = TRUE
  )

  # Eigendecomposition of the covariance
  eig <- eigen(out$sigma)

  # Definition of membership
  membership <- NULL
  for (i in 1:length(pk)) {
    membership <- c(membership, rep(i, each = pk[i]))
  }
  names(membership) <- colnames(out$data)
  out$membership <- membership

  # Re-naming the outputs
  out$adjacency <- out$theta

  # Definition of sparse principal components
  out$loadings <- round(eig$vectors, digits = 10)
  out$theta <- ifelse(out$loadings != 0, yes = 1, no = 0)
  rownames(out$theta) <- rownames(out$loadings) <- colnames(out$adjacency)
  colnames(out$theta) <- colnames(out$loadings) <- paste0("comp", 1:ncol(out$theta))

  # Definition of proportion of explained variance
  ev <- eig$values / sum(eig$values)
  names(ev) <- colnames(out$theta)
  out$ev <- ev

  # Re-arranging the output
  out <- out[c("data", "loadings", "theta", "ev", "membership", "omega", "phi", "C", "u")]
  if (!output_matrices) {
    out <- out[c("data", "loadings", "theta", "ev", "membership")]
  }

  # Defining the class
  class(out) <- "simulation_components"

  return(out)
}


#' Data simulation for multivariate regression
#'
#' Simulates data with outcome(s) and predictors, where only a subset of the
#' predictors actually contributes to the definition of the outcome(s).
#'
#' @inheritParams SimulateGraphical
#' @param pk vector with the number of predictors in each independent block of
#'   variables in \code{xdata}. The number of independent blocks, which
#'   determines the maximum number of orthogonal latent variables that can be
#'   simulated, is given by \code{length(pk)}.
#' @param family type of outcome. If \code{family="gaussian"}, normally
#'   distributed outcomes are simulated. If \code{family="binomial"} or
#'   \code{family="multinomial"}, binary outcome(s) are simulated from a
#'   multinomial distribution where the probability is defined from a linear
#'   combination of normally distributed outcomes.
#' @param N number of classes of the categorical outcome. Only used if
#'   \code{family="multinomial"}.
#' @param ev_xz vector of the expected proportions of explained variances for
#'   each of the orthogonal latent variables. It must contain values in ]0,1[,
#'   and must be a vector of length \code{length(pk)} or a single value to
#'   generate latent variables with the same expected proportion of explained
#'   variance.
#' @param adjacency_x optional matrix encoding the conditional independence
#'   structure between predictor variables in \code{xdata}. This argument must
#'   be a binary symmetric matrix of size \code{sum(pk)} with zeros on the
#'   diagonal.
#' @param nu_within expected density (number of edges over the number of node
#'   pairs) of the conditional independence graph in the within-group blocks for
#'   predictors. For independent predictors, use \code{nu_within=0}. This
#'   argument is only used if \code{adjancency_x} is not provided.
#' @param theta_xz optional binary matrix encoding the predictor variables from
#'   \code{xdata} (columns) contributing to the definition of the orthogonal
#'   latent outcomes from \code{zdata} (rows).
#' @param nu_xz expected proportion of relevant predictors over the total number
#'   of predictors to be used for the simulation of the orthogonal latent
#'   outcomes. This argument is only used if \code{theta_xz} is not provided.
#' @param theta_zy optional binary matrix encoding the latent variables from
#'   \code{zdata} (columns) contributing to the definition of the observed
#'   outcomes from \code{ydata} (rows). This argument must be a square matrix of
#'   size \code{length(pk)}. If \code{theta_zy} is a diagonal matrix, each
#'   latent variable contributes to the definition of one observed outcome so
#'   that there is a one-to-one relationship between latent and observed
#'   outcomes (i.e. they are collinear). Nonzero off-diagonal elements in
#'   \code{theta_zy} introduce some correlation between the observed outcomes by
#'   construction from linear combinations implicating common latent outcomes.
#'   This argument is only used if \code{eta} is not provided.
#' @param nu_zy probability for each of the off-diagonal elements in
#'   \code{theta_zy} to be a 1. If \code{nu_zy=0}, \code{theta_zy} is a diagonal
#'   matrix. This argument is only used if \code{theta_zy} is not provided.
#' @param eta optional matrix of coefficients used in the linear combination of
#'   latent outcomes to generate observed outcomes.
#' @param eta_set vector defining the range of values from which \code{eta} is
#'   sampled. This argument is only used if \code{eta} is not provided.
#' @param ev_xx expected proportion of explained variance by the first Principal
#'   Component (PC1) of a Principal Component Analysis applied on the
#'   predictors. This is the largest eigenvalue of the correlation (if
#'   \code{scale_ev=TRUE}) or covariance (if \code{scale_ev=FALSE}) matrix
#'   divided by the sum of eigenvalues. If \code{ev_xx=NULL} (the default), the
#'   constant u is chosen by maximising the contrast of the correlation matrix.
#'
#' @return A list with: \item{xdata}{simulated predictor data.}
#'   \item{ydata}{simulated outcome data.} \item{proba}{simulated probability of
#'   belonging to each outcome class. Only used for \code{family="binomial"} or
#'   \code{family="multinomial"}.} \item{logit_proba}{logit of the simulated
#'   probability of belonging to each outcome class. Only used for
#'   \code{family="binomial"} or \code{family="multinomial"}.}
#'   \item{zdata}{simulated data for orthogonal latent outcomes.}
#'   \item{beta}{matrix of true beta coefficients used to generate outcomes in
#'   \code{ydata} from predictors in \code{xdata}.} \item{theta}{binary matrix
#'   indicating the predictors from \code{xdata} contributing to the definition
#'   of each of the outcome variables in \code{ydata}.} \item{eta}{matrix of
#'   coefficients used in the linear combination of latent variables from
#'   \code{zdata} to define observed outcomes in \code{ydata}.}
#'   \item{theta_zy}{binary matrix indicating the latent variables from
#'   \code{zdata} used in the definition of observed outcomes in \code{ydata}.}
#'   \item{xi}{matrix of true beta coefficients used to generate orthogonal
#'   latent outcomes in \code{zdata} from predictors in \code{xdata}.}
#'   \item{theta_xz}{binary matrix indicating the predictors from \code{xdata}
#'   contributing to the definition of each of the latent outcome variables in
#'   \code{zdata}.} \item{omega_xz}{precision matrix for variables in
#'   \code{xdata} and \code{zdata}.} \item{adjacency}{binary matrix encoding the
#'   conditional independence structure between variables from \code{xdata}
#'   (\code{var}), \code{zdata} (\code{latent}) and \code{ydata}
#'   (\code{outcome}).}
#'
#' @details For a univariate outcome (\code{length(pk)=1}), the simulation is
#'   done in four steps where (i) predictors contributing to outcome definition
#'   are randomly sampled (with probability \code{nu_xz} for a given predictor
#'   to be picked), (ii) the conditional independence structure between the
#'   predictors is simulated (with probability \code{nu_within} for a given pair
#'   of predictors to be correlated, conditionally on all other variables),
#'   (iii) generation of a precision matrix (inverse covariance matrix) for all
#'   variables, where nonzero entries correspond to the predictors contributing
#'   to outcome definition or conditional correlation between the predictors,
#'   and (iv) data for both predictors and outcome is simulated from a single
#'   multivariate Normal distribution using the inverse precision matrix as
#'   covariance matrix.
#'
#'   To ensure that the generated precision matrix \eqn{\Omega} is positive
#'   definite, the diagonal entries are defined as described in
#'   \code{\link{MakePositiveDefinite}}. The conditional variance of the outcome
#'   \eqn{\Omega_{YY}} is chosen so that the proportion of variance in the
#'   outcome that is explained by the predictors is \code{ev_xz}.
#'
#'   For a multivariate outcome (\code{length(pk)>1}), we introduce independent
#'   groups of predictors and orthogonal latent variables (groups are defined in
#'   \code{pk}). Each latent variable is defined as a function of variables
#'   belonging to one group of predictors. The precision matrix is defined as
#'   described above for univariate outcomes. Subject to the re-ordering of its
#'   rows, this precision matrix is block-diagonal, encoding the independence
#'   between sets of variables made of (i) the groups of predictors, and (ii)
#'   their corresponding latent variable. The outcome variables are then
#'   constructed from a linear combination of the latent variables, allowing for
#'   contributing predictors belonging to different groups.
#'
#'   The use of latent variables in the multivariate case ensures that we can
#'   control the proportion of variance in the latent variable explained by the
#'   predictors (\code{ev_xz}).
#'
#' @family simulation functions
#'
#' @references \insertRef{ourstabilityselection}{fake}
#'
#' @examples
#' oldpar <- par(no.readonly = TRUE)
#' par(mar = c(5, 5, 5, 5))
#'
#' ## Continuous outcomes
#'
#' # Univariate outcome
#' set.seed(1)
#' simul <- SimulateRegression(pk = 15)
#' print(simul)
#' plot(simul)
#'
#' # Multivariate outcome
#' set.seed(1)
#' simul <- SimulateRegression(pk = c(5, 7, 3))
#' print(simul)
#' plot(simul)
#'
#' # Independent predictors
#' set.seed(1)
#' simul <- SimulateRegression(pk = c(5, 3), nu_within = 0)
#' print(simul)
#' plot(simul)
#'
#' # Blocks of strongly inter-connected predictors
#' set.seed(1)
#' simul <- SimulateRegression(
#'   pk = c(5, 5), nu_within = 0.5,
#'   v_within = c(0.5, 1), v_sign = -1, continuous = TRUE, pd_strategy = "min_eigenvalue"
#' )
#' print(simul)
#' Heatmap(
#'   mat = cor(simul$xdata),
#'   col = c("navy", "white", "red"),
#'   legend_range = c(-1, 1)
#' )
#' plot(simul)
#'
#'
#' ## Categorical outcomes
#'
#' # Binary outcome
#' set.seed(1)
#' simul <- SimulateRegression(pk = 20, family = "binomial")
#' print(simul)
#' table(simul$ydata[, 1])
#'
#' # Categorical outcome
#' set.seed(1)
#' simul <- SimulateRegression(pk = 20, family = "multinomial")
#' print(simul)
#' apply(simul$ydata, 2, sum)
#'
#' par(oldpar)
#' @export
SimulateRegression <- function(n = 100, pk = 10, N = 3,
                               family = "gaussian", ev_xz = 0.8,
                               adjacency_x = NULL, nu_within = 0.1,
                               theta_xz = NULL, nu_xz = 0.2,
                               theta_zy = NULL, nu_zy = 0.5,
                               eta = NULL, eta_set = c(-1, 1),
                               v_within = c(0.5, 1), v_sign = c(-1, 1), continuous = TRUE,
                               pd_strategy = "diagonally_dominant", ev_xx = NULL, scale_ev = TRUE,
                               u_list = c(1e-10, 1), tol = .Machine$double.eps^0.25) {
  # Checking the inputs
  if ((length(pk) > 1) & (family == "multinomial")) {
    stop("The simulation of multiple categorical outcomes is not possible with the current implementation.")
  }

  # Definition of the number of latent outcome variables
  q <- length(pk)
  p <- sum(pk)
  if (length(nu_xz) != q) {
    nu_xz <- rep(nu_xz[1], q)
  }
  if (length(ev_xz) != q) {
    ev_xz <- rep(ev_xz[1], q)
  }

  # Checking the values of ev_xz
  if (any(ev_xz <= 0) | any(ev_xz >= 1)) {
    stop("Invalid input for argument 'ev_xz'. Please provide values strictly between 0 and 1.")
  }

  # Simulation of the conditional independence structure with independent blocks
  if (is.null(adjacency_x)) {
    adjacency_x <- SimulateAdjacency(
      pk = pk, nu_between = 0, nu_within = nu_within,
      implementation = HugeAdjacency,
      topology = "random"
    )
  }

  # Simulation of the binary contribution status of predictors for latent outcome variables
  if (is.null(theta_xz)) {
    theta_xz <- SamplePredictors(pk = pk, q = q, nu = nu_xz, orthogonal = TRUE)
  }

  # Using the same support for all categories (multinomial only)
  if (family == "multinomial") {
    theta_xz <- matrix(rep(theta_xz, N), ncol = N)
    q <- N
    ev_xz <- rep(ev_xz, N)
  }

  # Setting row and column names
  colnames(theta_xz) <- paste0("latent", 1:ncol(theta_xz))
  rownames(theta_xz) <- paste0("var", 1:nrow(theta_xz))

  # Simulation of precision matrix for both predictors and latent outcomes
  big_theta <- cbind(
    rbind(matrix(0, nrow = q, ncol = q), theta_xz),
    rbind(t(theta_xz), adjacency_x)
  )
  rownames(big_theta) <- colnames(big_theta)
  out <- SimulatePrecision(
    theta = big_theta, v_within = v_within,
    v_sign = v_sign, continuous = continuous,
    pd_strategy = pd_strategy, ev_xx = ev_xx, scale = scale_ev, u_list = u_list, tol = tol
  )
  omega <- out$omega

  # Setting diagonal precision for latent outcomes to reach expected proportion of explained variance
  for (j in 1:q) {
    pred_ids <- seq(q + 1, q + p)
    omega[j, j] <- omega[j, pred_ids, drop = FALSE] %*% solve(omega[pred_ids, pred_ids]) %*% t(omega[j, pred_ids, drop = FALSE]) * 1 / ev_xz[j]
  }

  # Checking positive definiteness (mostly for multinomial?)
  eig <- eigen(omega, only.values = TRUE)$values
  if (any(eig <= 0)) {
    message("The requested proportion of explained variance could not be achieved.")
    ev_xz <- stats::optimise(TuneExplainedVarianceReg, interval = c(0, min(ev_xz)), omega = omega, q = q, p = p)$minimum
    ev_xz <- rep(ev_xz, q)
    for (j in 1:q) {
      pred_ids <- seq(q + 1, q + p)
      omega[j, j] <- omega[j, pred_ids, drop = FALSE] %*% solve(omega[pred_ids, pred_ids]) %*% t(omega[j, pred_ids, drop = FALSE]) * 1 / ev_xz[j]
    }
  }

  # Computing the covariance matrix
  sigma <- solve(omega)

  # Computing the regression coefficients from X to Z
  xi <- solve(sigma[grep("var", colnames(omega)), grep("var", colnames(omega))]) %*% sigma[grep("var", colnames(sigma)), grep("latent", colnames(sigma))]
  colnames(xi) <- colnames(theta_xz)

  # Simulation of data from multivariate normal distribution
  x <- MASS::mvrnorm(n, rep(0, p + q), sigma)
  rownames(x) <- paste0("obs", 1:nrow(x))

  # Separating predictors from latent outcome variables
  xdata <- x[, grep("var", colnames(x)), drop = FALSE]
  zdata <- x[, grep("latent", colnames(x)), drop = FALSE]

  # Simulation of eta coefficients to get observed outcomes from latent outcomes
  if (is.null(eta)) {
    if (family == "multinomial") {
      # Variables in Z and Y are the same for multinomial as restricted to one categorical outcome
      theta_zy <- eta <- diag(N)
    } else {
      if (is.null(theta_zy)) {
        theta_zy <- SamplePredictors(pk = q, q = q, nu = nu_zy, orthogonal = FALSE)
      }
      eta <- matrix(stats::runif(q * q, min = min(eta_set), max = max(eta_set)),
        ncol = q, nrow = q
      )
      eta <- eta * theta_zy
    }
  } else {
    theta_zy <- ifelse(eta != 0, yes = 1, no = 0)
  }
  rownames(eta) <- rownames(theta_zy) <- paste0("latent", 1:q)
  colnames(eta) <- colnames(theta_zy) <- paste0("outcome", 1:q)
  ydata <- zdata %*% eta

  # Computing the xy coefficients and binary contribution status
  beta <- xi %*% eta
  beta <- base::zapsmall(beta)
  theta <- ifelse(beta != 0, yes = 1, no = 0)

  # Compute binary outcome for logistic regression
  if (family == "binomial") {
    # Variability is coming from binomial distribution for logistic
    ydata <- xdata %*% beta

    # Sampling from (series of) binomial distributions
    proba <- matrix(NA, nrow = n, ncol = q)
    for (j in 1:q) {
      proba[, j] <- 1 / (1 + exp(-ydata[, j])) # inverse logit
    }

    ydata_cat <- matrix(0, nrow = n, ncol = q)
    for (j in 1:q) {
      for (i in 1:n) {
        ydata_cat[i, j] <- stats::rbinom(n = 1, size = 1, prob = proba[i, j])
      }
    }

    # Setting row and column names
    rownames(ydata_cat) <- rownames(proba) <- rownames(ydata)
    colnames(ydata_cat) <- colnames(proba) <- colnames(ydata)
  }

  if (family == "multinomial") {
    # Variability is coming from multinomial distribution for logistic
    ydata <- xdata %*% beta

    proba <- matrix(NA, nrow = n, ncol = q)
    for (j in 1:q) {
      proba[, j] <- 1 / (1 + exp(-ydata[, j])) # inverse logit
    }

    # Sampling from multinomial distribution
    ydata_cat <- matrix(0, nrow = n, ncol = q)
    for (i in 1:n) {
      ydata_cat[i, ] <- stats::rmultinom(n = 1, size = 1, prob = proba[i, ])[, 1]
    }

    # Setting row and column names
    rownames(ydata_cat) <- rownames(proba) <- rownames(ydata)
    colnames(ydata_cat) <- colnames(proba) <- colnames(ydata)
  }

  # Extracting the conditional independence structure between x, z and y
  adjacency <- rbind(
    cbind(matrix(0, nrow = q, ncol = q), t(theta_zy), matrix(0, nrow = q, ncol = p)),
    cbind(rbind(theta_zy, matrix(0, nrow = p, ncol = q)), big_theta)
  )
  rownames(adjacency) <- colnames(adjacency) <- c(colnames(theta_zy), rownames(big_theta))

  # Return the simulated X and Y
  if (family %in% c("binomial", "multinomial")) {
    if (family == "binomial") {
      out <- list(
        xdata = xdata, ydata = ydata_cat,
        proba = proba, logit_proba = ydata,
        zdata = zdata,
        beta = beta, theta = theta,
        eta = eta, theta_zy = theta_zy,
        xi = xi, theta_xz = theta_xz,
        ev_xz = ev_xz,
        omega_xz = omega,
        adjacency = adjacency
      )
    }
    if (family == "multinomial") {
      out <- list(
        xdata = xdata, ydata = ydata_cat,
        proba = proba, logit_proba = ydata,
        zdata = zdata,
        beta = beta, theta = theta,
        theta_zy = theta_zy,
        xi = xi, theta_xz = theta_xz,
        ev_xz = ev_xz,
        omega_xz = omega,
        adjacency = adjacency
      )
    }
  } else {
    out <- list(
      xdata = xdata, ydata = ydata, zdata = zdata,
      beta = beta, theta = theta,
      eta = eta, theta_zy = theta_zy,
      xi = xi, theta_xz = theta_xz,
      ev_xz = ev_xz,
      omega_xz = omega,
      adjacency = adjacency
    )
  }

  # Defining the class
  class(out) <- "simulation_regression"

  return(out)
}


#' Simulation of data with underlying clusters
#'
#' Simulates mixture multivariate Normal data with clusters of items (rows)
#' sharing similar profiles along (a subset of) attributes (columns).
#'
#' @inheritParams SimulateGraphical
#' @param n vector of the number of items per cluster in the simulated data. The
#'   total number of items is \code{sum(n)}.
#' @param pk vector of the number of attributes in the simulated data.
#' @param adjacency optional binary and symmetric adjacency matrix encoding the
#'   conditional independence structure between attributes.
#' @param theta_xc optional binary matrix encoding which attributes (columns)
#'   contribute to the clustering structure between which clusters (rows). If
#'   \code{theta_xc=NULL}, variables contributing to the clustering are sampled
#'   with probability \code{nu_xc}.
#' @param nu_xc expected proportion of variables contributing to the clustering
#'   over the total number of variables. Only used if \code{theta_xc} is not
#'   provided.
#' @param ev_xc vector of marginal expected proportion of explained variance by
#'   each attribute contributing to the clustering.
#' @param ev_xx expected proportion of explained variance by the first Principal
#'   Component (PC1) of a Principal Component Analysis applied on the
#'   predictors. This is the largest eigenvalue of the correlation (if
#'   \code{scale=TRUE}) or covariance (if \code{scale=FALSE}) matrix divided by
#'   the sum of eigenvalues. If \code{ev_xx=NULL} (the default), the constant u
#'   is chosen by maximising the contrast of the correlation matrix.
#'
#' @seealso \code{\link{MakePositiveDefinite}}
#' @family simulation functions
#'
#' @return A list with: \item{data}{simulated data with \code{sum(n)}
#'   observation and \code{sum(pk)} variables} \item{theta}{simulated (true)
#'   cluster membership.} \item{adjacency}{adjacency matrix of the graph
#'   encoding the conditional independence structure between variables.}
#'   \item{theta_xc}{binary vector encoding variables contributing to the
#'   clustering structure.} \item{ev}{vector of marginal expected proportions of
#'   explained variance for each variable.} \item{omega}{simulated (true)
#'   precision matrix. Only returned if \code{output_matrices=TRUE}.}
#'   \item{phi}{simulated (true) partial correlation matrix. Only returned if
#'   \code{output_matrices=TRUE}.} \item{sigma}{simulated (true) covariance
#'   matrix. Only returned if \code{output_matrices=TRUE}.} \item{u}{value of
#'   the constant u used for the simulation of \code{omega}. Only returned if
#'   \code{output_matrices=TRUE}.} \item{mu_mixture}{simulated (true)
#'   cluster-specific means. Only returned if \code{output_matrices=TRUE}.}
#'
#' @examples
#'
#' ## Example with 3 clusters
#'
#' # Data simulation
#' set.seed(1)
#' simul <- SimulateClustering(
#'   n = c(10, 30, 15),
#'   nu_xc = 1,
#'   ev_xc = 0.5
#' )
#' print(simul)
#' plot(simul)
#'
#' # Checking the proportion of explained variance
#' x <- simul$data[, 1]
#' z <- as.factor(simul$theta)
#' summary(lm(x ~ z)) # R-squared
#'
#'
#' ## Example with 2 variables contributing to clustering
#'
#' # Data simulation
#' set.seed(1)
#' simul <- SimulateClustering(
#'   n = c(20, 10, 15), pk = 10,
#'   theta_xc = c(1, 1, rep(0, 8)),
#'   ev_xc = 0.8
#' )
#' print(simul)
#' plot(simul)
#'
#' # Visualisation of the data
#' par(mar = c(5, 5, 5, 5))
#' Heatmap(
#'   mat = simul$data,
#'   colours = c("navy", "white", "red")
#' )
#' simul$ev # marginal proportions of explained variance
#'
#' # Visualisation along contributing variables
#' plot(simul$data[, 1:2], col = simul$theta)
#'
#'
#' ## Example with different levels of separation
#'
#' # Data simulation
#' set.seed(1)
#' simul <- SimulateClustering(
#'   n = c(20, 10, 15), pk = 10,
#'   theta_xc = c(1, 1, rep(0, 8)),
#'   ev_xc = c(0.99, 0.5, rep(0, 8))
#' )
#' simul$ev
#'
#' # Visualisation along contributing variables
#' plot(simul$data[, 1:2], col = simul$theta)
#'
#'
#' ## Example with correlated contributors
#'
#' # Data simulation
#' pk <- 10
#' adjacency <- matrix(0, pk, pk)
#' adjacency[1, 2] <- adjacency[2, 1] <- 1
#' set.seed(1)
#' simul <- SimulateClustering(
#'   n = c(200, 100, 150), pk = pk,
#'   theta_xc = c(1, 1, rep(0, 8)),
#'   ev_xc = c(0.9, 0.8, rep(0, 8)),
#'   adjacency = adjacency,
#'   pd_strategy = "min_eigenvalue",
#'   v_within = 0.6, v_sign = -1
#' )
#'
#' # Visualisation along contributing variables
#' plot(simul$data[, 1:2], col = simul$theta)
#'
#' # Checking marginal proportions of explained variance
#' mymodel <- lm(simul$data[, 1] ~ as.factor(simul$theta))
#' summary(mymodel)$r.squared
#' mymodel <- lm(simul$data[, 2] ~ as.factor(simul$theta))
#' summary(mymodel)$r.squared
#'
#' @export
SimulateClustering <- function(n = c(10, 10), pk = 10, adjacency = NULL,
                               theta_xc = NULL, nu_xc = 0.1, ev_xc = NULL,
                               implementation = HugeAdjacency, topology = "random",
                               nu_within = 0, nu_between = NULL, nu_mat = NULL,
                               v_within = c(0.5, 1), v_between = c(0, 0.1),
                               v_sign = c(-1, 1), continuous = TRUE,
                               pd_strategy = "diagonally_dominant", ev_xx = NULL, scale_ev = TRUE,
                               u_list = c(1e-10, 1), tol = .Machine$double.eps^0.25,
                               scale = TRUE,
                               output_matrices = FALSE) {
  # Checking the inputs
  if (!is.null(theta_xc)) {
    if (inherits(theta_xc, "numeric")) {
      theta_xc <- rbind(theta_xc)
    }
    if (nrow(theta_xc) == 1) {
      tmp_theta_xc <- theta_xc
      theta_xc <- NULL
      for (i in 1:length(n)) {
        theta_xc <- rbind(theta_xc, tmp_theta_xc)
      }
    }
    if (is.null(pk)) {
      pk <- ncol(theta_xc)
    } else {
      if (sum(pk) != ncol(theta_xc)) {
        # warning("Arguments 'pk' and 'theta_xc' are not compatible. Argument 'pk' has been set to ncol('theta_xc').")
        pk <- ncol(theta_xc)
      }
    }
  }
  if (!is.null(adjacency)) {
    if (ncol(adjacency) != nrow(adjacency)) {
      stop("Invalid input for argument 'adjacency'. It must be a square matrix (same number of rows and columns).")
    }
    if (is.null(pk)) {
      pk <- ncol(theta_xc)
    } else {
      if (sum(pk) != ncol(adjacency)) {
        warning("Arguments 'pk' and 'theta_xc' are not compatible. Argument 'pk' has been set to ncol('adjacency').")
        pk <- ncol(adjacency)
      }
    }
  }

  # Using multi-block simulator with unconnected blocks
  out <- SimulateGraphical(
    n = sum(n), pk = pk, theta = adjacency,
    implementation = implementation,
    topology = topology,
    nu_within = nu_within,
    nu_between = nu_between,
    nu_mat = nu_mat,
    v_within = v_within,
    v_between = v_between,
    v_sign = v_sign,
    continuous = continuous,
    pd_strategy = pd_strategy,
    ev_xx = ev_xx,
    scale_ev = scale_ev,
    u_list = u_list,
    tol = tol,
    scale = scale,
    output_matrices = TRUE
  )

  # Defining number of clusters
  nc <- length(n)

  # Defining variables contributing to the clustering
  if (is.null(theta_xc)) {
    theta_xc <- matrix(rep(
      SamplePredictors(pk = sum(pk), q = 1, nu = nu_xc, orthogonal = TRUE)[, 1],
      length(n)
    ),
    nrow = length(n), byrow = TRUE
    )
  }
  rownames(theta_xc) <- paste0("cluster", 1:nrow(theta_xc))

  # Simulating marginal proportions of explained variance
  if (is.null(ev_xc)) {
    ev_xc <- stats::runif(n = sum(pk))
  } else {
    if (length(ev_xc) == 1) {
      ev_xc <- rep(ev_xc, sum(pk))
    }
  }

  # Re-naming the outputs
  out$adjacency <- out$theta

  # Definition of membership
  theta <- NULL
  for (i in 1:length(n)) {
    theta <- c(theta, rep(i, each = n[i]))
  }
  names(theta) <- rownames(out$data)
  out$theta <- theta

  # Simulating the cluster-specific means
  mu_mixture <- matrix(NA, nrow = length(unique(theta)), ncol = sum(pk))
  if (length(n) > 1) {
    mu_mat <- matrix(NA, nrow = sum(n), ncol = sum(pk))
    for (k in 1:ncol(mu_mat)) {
      # Sampling initial values for cluster-specific means
      mu <- stats::rnorm(n = nc, mean = 0, sd = 1)
      mu <- mu * theta_xc[, k]

      # Attributing cluster-specific means to observations
      for (i in 1:nrow(mu_mat)) {
        mu_mat[i, k] <- mu[theta[i]]
      }

      # # Defining variance to reach expected proportion of e.v.
      # var_mu <- ev_xc[k] * 1 / (1 - ev_xc[k])
      #
      # # Scaling to ensure mean of zero and defined variance
      # mu_mat[, k] <- scale(mu_mat[, k])
      # mu_mat[, k] <- mu_mat[, k] * sqrt(var_mu)
      # # Equivalent to: sqrt(var_mu)*(mu_mat[, k]-mean(mu_mat[,k]))/sqrt(1/(nrow(mu_mat)-1)*sum((mu_mat[, k]-mean(mu_mat[, k]))^2))
      # # Equivalent to: sqrt(var_mu)*(mu_mat[, k]-1/nrow(mu_mat)*sum(table(theta)*mu))/sqrt(1/(nrow(mu_mat)-1)*sum((mu_mat[, k]-mean(mu_mat[, k]))^2))
      # # Equivalent to: sqrt(var_mu)*(mu_mat[, k]-1/nrow(mu_mat)*sum(table(theta)*mu))/sqrt(1/(nrow(mu_mat)-1)*sum(table(theta)*(mu-1/nrow(mu_mat)*sum(table(theta)*mu))^2))
      #
      # # Storing cluster-specific means
      # mu_mixture[, k] <- mu_mat[!duplicated(theta), k]
      #
      # # Adding cluster-specific means
      # out$data[, k] <- out$data[, k] + mu_mat[, k]

      # Ensuring that the grouping structure is going to represent desired proportion of variance
      if (any(theta_xc[, k] != 0)) {
        mu_mat[, k] <- scale(mu_mat[, k]) * sqrt(ev_xc[k]) * sqrt(diag(out$sigma)[k])
        out$data[, k] <- out$data[, k] * sqrt(1 - ev_xc[k])
      }

      # Adding cluster-specific means
      out$data[, k] <- out$data[, k] + mu_mat[, k]

      # Storing cluster-specific means
      mu_mixture[, k] <- mu_mat[!duplicated(theta), k]
    }
    mu_mat[is.na(mu_mat)] <- 0
  }

  # Updating the within-cluster covariance matrix
  out$sigma <- out$sigma * sqrt(cbind(1 - ev_xc) %*% rbind(1 - ev_xc))

  # Definition of contributing variables
  colnames(theta_xc) <- colnames(out$data)
  out$theta_xc <- theta_xc
  out$ev <- ev_xc * ifelse(apply(theta_xc, 2, sum) != 0, yes = 1, no = 0)

  # Returning true cluster-specific means
  if (output_matrices) {
    out$mu_mixture <- mu_mixture
  }

  # Defining the class
  class(out) <- "simulation_clustering"

  return(out)
}


#' Simulation of undirected graph with block structure
#'
#' Simulates the adjacency matrix of an unweighted, undirected graph with no
#' self-loops. If \code{topology="random"}, different densities in diagonal
#' (\code{nu_within}) compared to off-diagonal (\code{nu_between}) blocks can be
#' used.
#'
#' @inheritParams SimulateGraphical
#'
#' @return A symmetric adjacency matrix encoding an unweighted, undirected graph
#'   with no self-loops, and with different densities in diagonal compared to
#'   off-diagonal blocks.
#'
#' @details Random graphs are simulated using the Erdos-Renyi algorithm.
#'   Scale-free graphs are simulated using a preferential attachment algorithm.
#'   More details are provided in \code{\link[huge]{huge.generator}}.
#'
#' @family simulation functions
#'
#' @references \insertRef{ourstabilityselection}{fake}
#'
#' \insertRef{huge}{fake}
#'
#' @examples
#' # Simulation of a scale-free graph with 20 nodes
#' adjacency <- SimulateAdjacency(pk = 20, topology = "scale-free")
#' plot(adjacency)
#'
#' # Simulation of a random graph with three connected components
#' adjacency <- SimulateAdjacency(
#'   pk = rep(10, 3),
#'   nu_within = 0.7, nu_between = 0
#' )
#' plot(adjacency)
#'
#' # Simulation of a random graph with block structure
#' adjacency <- SimulateAdjacency(
#'   pk = rep(10, 3),
#'   nu_within = 0.7, nu_between = 0.03
#' )
#' plot(adjacency)
#'
#' # User-defined function for graph simulation
#' CentralNode <- function(pk, hub = 1) {
#'   theta <- matrix(0, nrow = sum(pk), ncol = sum(pk))
#'   theta[hub, ] <- 1
#'   theta[, hub] <- 1
#'   diag(theta) <- 0
#'   return(theta)
#' }
#' simul <- SimulateAdjacency(pk = 10, implementation = CentralNode)
#' plot(simul) # star
#' simul <- SimulateAdjacency(pk = 10, implementation = CentralNode, hub = 2)
#' plot(simul) # variable 2 is the central node
#' @export
SimulateAdjacency <- function(pk = 10,
                              implementation = HugeAdjacency,
                              topology = "random",
                              nu_within = 0.1,
                              nu_between = 0,
                              nu_mat = NULL,
                              ...) {
  # Storing all arguments
  args <- c(mget(ls()), list(...))

  # Checking the inputs
  if (topology != "random") {
    if (length(pk) > 1) {
      pk <- sum(pk)
      warning(paste0("Multi-block simulations are only allowed with topology='random'. Argument 'pk' has been set to ", pk, "."))
    }
  }

  # Creating the matrix of probabilities
  if (is.null(nu_mat)) {
    nu_mat <- diag(length(pk)) * nu_within
    nu_mat[upper.tri(nu_mat)] <- nu_between
    nu_mat[lower.tri(nu_mat)] <- nu_between
  } else {
    if ((ncol(nu_mat) != length(pk)) & (nrow(nu_mat) != length(pk))) {
      stop("Arguments 'pk' and 'nu_mat' are not compatible. They correspond to different numbers of communities. The number of rows and columns in 'nu_mat' must be equal to the length of the vector 'pk'.")
    }
  }

  # Creating matrix with block indices
  bigblocks <- BlockMatrix(pk)
  bigblocks_vect <- bigblocks[upper.tri(bigblocks)]

  # Making as factor to allow for groups with 1 variable (for clustering)
  bigblocks_vect <- factor(bigblocks_vect, levels = seq(1, max(bigblocks)))
  block_ids <- unique(as.vector(bigblocks))

  # Creating matrix with block structure
  blocks <- BlockStructure(pk)

  # Identifying relevant arguments
  if (!"..." %in% names(formals(implementation))) {
    ids <- which(names(args) %in% names(formals(implementation)))
    args <- args[ids]
  }

  # Simulation of the adjacency matrix
  if ("nu" %in% names(formals(implementation))) {
    if (length(pk) > 1) {
      # Initialising theta
      theta <- matrix(0, nrow = sum(pk), ncol = sum(pk))
      theta_vect <- theta[upper.tri(theta)]

      # # Allowing for different densities in within and between blocks
      # theta_w <- do.call(implementation, args = c(args, list(nu = nu_within)))
      # theta_w_vect <- theta_w[upper.tri(theta_w)]
      # theta_b <- do.call(implementation, args = c(args, list(nu = nu_between)))
      # theta_b_vect <- theta_b[upper.tri(theta_b)]

      # Filling within and between blocks
      for (k in block_ids) {
        tmpids <- which(blocks == k, arr.ind = TRUE)
        i <- tmpids[1]
        j <- tmpids[2]
        theta_w <- do.call(implementation, args = c(args, list(nu = nu_mat[i, j])))
        theta_w_vect <- theta_w[upper.tri(theta_w)]
        theta_vect[bigblocks_vect == k] <- theta_w_vect[bigblocks_vect == k]
        # if (k %in% unique(diag(bigblocks))) {
        #   theta_vect[bigblocks_vect == k] <- theta_w_vect[bigblocks_vect == k]
        # } else {
        #   theta_vect[bigblocks_vect == k] <- theta_b_vect[bigblocks_vect == k]
        # }
      }
      theta[upper.tri(theta)] <- theta_vect
      theta <- theta + t(theta)
    } else {
      theta <- do.call(implementation, args = c(args, list(nu = nu_within)))
    }
  } else {
    theta <- do.call(implementation, args = c(args))
  }

  # Ensuring the adjacency matrix is symmetric (undirected graph) with no self-loops
  theta <- ifelse(theta + t(theta) != 0, yes = 1, no = 0)
  diag(theta) <- 0

  # Setting variable names
  colnames(theta) <- rownames(theta) <- paste0("var", 1:ncol(theta))

  # Defining the class
  class(theta) <- c("matrix", "adjacency_matrix")

  return(theta)
}


#' Simulation of undirected graph
#'
#' Simulates the adjacency matrix encoding an unweighted, undirected graph with
#' no self-loops.
#'
#' @inheritParams SimulateGraphical
#' @param pk number of nodes.
#' @param nu expected density (number of edges over the number of node pairs) of
#'   the graph. This argument is only used for \code{topology="random"} or
#'   \code{topology="cluster"} (see argument \code{prob} in
#'   \code{\link[huge]{huge.generator}}).
#' @param ... additional arguments to be passed to
#'   \code{\link[huge]{huge.generator}}.
#'
#' @return A symmetric adjacency matrix encoding an unweighted, undirected graph
#'   with no self-loops.
#'
#' @keywords internal
HugeAdjacency <- function(pk = 10, topology = "random", nu = 0.1, ...) {
  # Storing extra arguments
  extra_args <- list(...)

  # Extracting relevant extra arguments
  tmp_extra_args <- MatchingArguments(extra_args = extra_args, FUN = huge::huge.generator)
  tmp_extra_args <- tmp_extra_args[!names(tmp_extra_args) %in% c("n", "d", "prob", "graph", "verbose")]

  # Running simulation model
  mymodel <- do.call(huge::huge.generator, args = c(
    list(
      n = 2, d = sum(pk), prob = nu,
      graph = topology, verbose = FALSE
    ),
    tmp_extra_args
  ))
  theta <- as.matrix(mymodel$theta)

  # Re-organising the variables to avoid having centrality related to variable ID (e.g. for scale-free models)
  ids <- sample(ncol(theta))
  theta <- theta[ids, ids]

  return(theta)
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
