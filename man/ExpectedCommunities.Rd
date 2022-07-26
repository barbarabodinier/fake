% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tools.R
\name{ExpectedCommunities}
\alias{ExpectedCommunities}
\title{Expected community structure}
\usage{
ExpectedCommunities(pk, nu_within = 0.1, nu_between = 0, nu_mat = NULL)
}
\arguments{
\item{pk}{vector of the number of variables per group in the simulated
dataset. The number of nodes in the simulated graph is \code{sum(pk)}. With
multiple groups, the simulated (partial) correlation matrix has a block
structure, where blocks arise from the integration of the \code{length(pk)}
groups. This argument is only used if \code{theta} is not provided.}

\item{nu_within}{probability of having an edge between two nodes belonging to
the same group, as defined in \code{pk}. If \code{length(pk)=1}, this is
the expected density of the graph. If \code{implementation=HugeAdjacency},
this argument is only used for \code{topology="random"} or
\code{topology="cluster"} (see argument \code{prob} in
\code{\link[huge]{huge.generator}}). Only used if \code{nu_mat} is not
provided.}

\item{nu_between}{probability of having an edge between two nodes belonging
to different groups, as defined in \code{pk}. By default, the same density
is used for within and between blocks (\code{nu_within}=\code{nu_between}).
Only used if \code{length(pk)>1}. Only used if \code{nu_mat} is not
provided.}

\item{nu_mat}{matrix of probabilities of having an edge between nodes
belonging to a given pair of node groups defined in \code{pk}.}
}
\value{
A list with: \item{total_within_degree_c}{total within degree by node
  group, i.e. sum of expected within degree over all nodes in a given group.}
  \item{total_between_degree}{total between degree by node group, i.e. sum of
  expected between degree over all nodes in a given group.}
  \item{weak_community}{binary indicator for a given node group to be an
  expected weak community.} \item{total_number_edges_c}{matrix of expected
  number of edges between nodes from a given node pair.}
  \item{modularity}{expected modularity (see
  \code{\link[igraph]{modularity}}).}
}
\description{
Computes expected metrics related to the community structure of a graph
simulated with given parameters.
}
\details{
Given a group of nodes, the within degree \eqn{d^w_i} of node
  \eqn{i} is defined as the number of nodes from the same group node \eqn{i}
  is connected to. The between degree \eqn{d^b_i} is the number of nodes from
  other groups node \eqn{i} is connected to. A weak community in the network
  is defined as a group of nodes for which the total within degree (sum of
  the \eqn{d^w_i} for all nodes in the community) is stricly greater than the
  total between degree (sum of \eqn{d^b_i} for all nodes in the community).
  For more details, see
  \href{http://networksciencebook.com/chapter/9#basics}{Network Science} by
  Albert-Laszlo Barabasi.

  The expected total within and between degrees for the groups defined in
  \code{pk} in a network simulated using \code{SimulateAdjacency} can be
  computed given the group sizes (stored in \code{pk}) and probabilities of
  having an edge between nodes from a given group pair (defined by
  \code{nu_within} and \code{nu_between} or by \code{nu_mat}). The expected
  presence of weak communities can be inferred from these quantities.

  The expected modularity, measuring the difference between observed and
  expected number of within-community edges, is also returned. For more
  details on this metric, see \code{\link[igraph]{modularity}}.
}
\examples{
# Simulation parameters
pk <- rep(20, 4)
nu_within <- 0.8
nu_between <- 0.1

# Expected metrics
expected <- ExpectedCommunities(
  pk = pk,
  nu_within = nu_within,
  nu_between = nu_between
)

# Example of simulated graph
set.seed(1)
theta <- SimulateAdjacency(
  pk = pk,
  nu_within = nu_within,
  nu_between = nu_between
)

# Comparing observed and expected numbers of edges
bigblocks <- BlockMatrix(pk)
BlockStructure(pk)
sum(theta[which(bigblocks == 2)]) / 2
expected$total_number_edges_c[1, 2]

# Comparing observed and expected modularity
igraph::modularity(igraph::graph_from_adjacency_matrix(theta, mode = "undirected"),
  membership = rep.int(1:length(pk), times = pk)
)
expected$modularity

}
\seealso{
\code{\link{SimulateGraphical}}, \code{\link{SimulateAdjacency}},
  \code{\link{MinWithinProba}}
}
