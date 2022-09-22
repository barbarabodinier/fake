#' @export
print.simulation_graphical_model <- function(x, ...) {
  cat(paste0("Multivariate Normal data with underlying structure of a graphical model."))
  cat("\n")
}


#' @export
print.simulation_structural_causal_model <- function(x, ...) {
  cat(paste0("Multivariate Normal data with underlying structure of a structural causal model."))
  cat("\n")
}


#' @export
print.simulation_clustering <- function(x, ...) {
  cat(paste0("Multivariate Normal data with underlying clusters of participants along (a subset of) variables."))
  cat("\n")
}


#' @export
print.simulation_components <- function(x, ...) {
  cat(paste0("Multivariate Normal data with independent groups of variables."))
  cat("\n")
}


#' @export
print.simulation_regression <- function(x, ...) {
  cat(paste0("Multivariate Normal data with predictors and outcome(s)."))
  cat("\n")
}


#' @export
summary.simulation_graphical_model <- function(object, ...) {
  cat(paste0("Number of observations: ", nrow(object$data)))
  cat("\n")
  cat(paste0("Number of variables (nodes): ", ncol(object$data)))
  cat("\n")
  cat(paste0("Number of edges: ", sum(object$theta == 1) / 2))
  cat("\n")
}


#' @export
summary.simulation_structural_causal_model <- function(object, ...) {
  cat(paste0("Number of observations: ", nrow(object$data)))
  cat("\n")
  cat(paste0("Number of variables (nodes): ", ncol(object$data)))
  cat("\n")
  cat(paste0("Number of arrows: ", sum(object$theta == 1)))
  cat("\n")
}


#' @export
summary.simulation_clustering <- function(object, ...) {
  cat(paste0("Number of observations: ", nrow(object$data)))
  cat("\n")
  cat(paste0("Number of clusters: ", max(object$theta)))
  for (k in 1:max(object$theta)) {
    cat("\n")
    cat(paste0("- Cluster ", k, " (N=", sum(object$theta == k), " observations)"))
  }
  cat("\n")
  cat("\n")
  cat(paste0("Number of variables: ", ncol(object$data)))
  cat("\n")
  cat(paste0("Number of variables contributing to the clustering: ", sum(object$theta_xc)))
  cat("\n")
}


#' @export
summary.simulation_components <- function(object, ...) {
  cat(paste0("Number of observations: ", nrow(object$data)))
  cat("\n")
  cat("\n")
  cat(paste0("Number of variables: ", ncol(object$data)))
  cat("\n")
  cat(paste0("Number of independent groups of variables: ", max(object$membership)))
  for (k in 1:max(object$membership)) {
    cat("\n")
    cat(paste0("- Group ", k, " (N=", sum(object$membership == k), " variables)"))
  }
  cat("\n")
}


#' @export
summary.simulation_regression <- function(object, ...) {
  cat(paste0("Number of observations: ", nrow(object$xdata)))
  cat("\n")
  cat(paste0("Number of outcome variable(s): ", ncol(object$ydata)))
  cat("\n")
  cat(paste0("Number of predictor variables: ", ncol(object$xdata)))
  cat("\n")
  cat(paste0(
    "Number of predictor variables contributing to the outcome(s): ",
    sum(apply(object$beta, 1, sum) != 0)
  ))
  cat("\n")
}


#' @export
plot.simulation_graphical_model <- function(x, ...) {
  mygraph <- igraph::graph_from_adjacency_matrix(x$theta, mode = "undirected")

  # Formatting vertices
  mydegrees <- igraph::degree(mygraph)
  igraph::V(mygraph)$size <- as.numeric(as.character(cut(mydegrees, breaks = 4, labels = c(3, 4, 5, 6))))
  igraph::V(mygraph)$color <- "skyblue"
  igraph::V(mygraph)$frame.color <- igraph::V(mygraph)$color
  igraph::V(mygraph)$label.family <- "sans"
  igraph::V(mygraph)$label.cex <- as.numeric(as.character(cut(mydegrees, breaks = 4, labels = c(0.4, 0.45, 0.5, 0.55))))
  igraph::V(mygraph)$label.color <- "grey20"

  # Formatting edges
  igraph::E(mygraph)$color <- "grey60"
  igraph::E(mygraph)$width <- 0.5

  igraph::plot.igraph(mygraph, ...)

  return(invisible(mygraph))
}


#' @export
plot.simulation_structural_causal_model <- function(x, ...) {
  mygraph <- igraph::graph_from_adjacency_matrix(x$theta, mode = "directed")

  # Formatting vertices
  mydegrees <- igraph::degree(mygraph)
  igraph::V(mygraph)$size <- as.numeric(as.character(cut(mydegrees, breaks = 4, labels = c(3, 4, 5, 6))))
  igraph::V(mygraph)$color <- "skyblue"
  igraph::V(mygraph)$frame.color <- igraph::V(mygraph)$color
  igraph::V(mygraph)$label.family <- "sans"
  igraph::V(mygraph)$label.cex <- as.numeric(as.character(cut(mydegrees, breaks = 4, labels = c(0.4, 0.45, 0.5, 0.55))))
  igraph::V(mygraph)$label.color <- "grey20"

  # Formatting edges
  igraph::E(mygraph)$color <- "grey60"
  igraph::E(mygraph)$width <- 0.5

  igraph::plot.igraph(mygraph, ...)

  return(invisible(mygraph))
}


#' @export
plot.adjacency_matrix <- function(x, ...) {
  mygraph <- igraph::graph_from_adjacency_matrix(x, mode = "undirected")

  # Formatting vertices
  mydegrees <- igraph::degree(mygraph)
  igraph::V(mygraph)$size <- as.numeric(as.character(cut(mydegrees, breaks = 4, labels = c(3, 4, 5, 6))))
  igraph::V(mygraph)$color <- "skyblue"
  igraph::V(mygraph)$frame.color <- igraph::V(mygraph)$color
  igraph::V(mygraph)$label.family <- "sans"
  igraph::V(mygraph)$label.cex <- as.numeric(as.character(cut(mydegrees, breaks = 4, labels = c(0.4, 0.45, 0.5, 0.55))))
  igraph::V(mygraph)$label.color <- "grey20"

  # Formatting edges
  igraph::E(mygraph)$color <- "grey60"
  igraph::E(mygraph)$width <- 0.5

  igraph::plot.igraph(mygraph, ...)
}


#' @export
plot.simulation_clustering <- function(x, ...) {
  # Visualisation of Euclidian distances along the contributing variable
  Heatmap(
    mat = as.matrix(stats::dist(x$data[, which(apply(x$theta_xc, 2, sum) != 0), drop = FALSE]))
  )
  graphics::title("Distances across variables contributing to clustering")
}


#' @export
plot.simulation_components <- function(x, ...) {
  Heatmap(
    mat = stats::cor(x$data),
    col = c("navy", "white", "red"),
    legend_range = c(-1, 1)
  )
  graphics::title("Pearson's correlations")
}


#' @export
plot.simulation_regression <- function(x, ...) {
  mygraph <- igraph::graph_from_adjacency_matrix(x$adjacency, mode = "undirected")

  # Formatting vertices
  mydegrees <- igraph::degree(mygraph)
  igraph::V(mygraph)$size <- as.numeric(as.character(cut(mydegrees, breaks = 4, labels = c(3, 4, 5, 6))))
  igraph::V(mygraph)$color <- c(
    rep("red", ncol(x$ydata)),
    rep("orange", ncol(x$zdata)),
    rep("skyblue", ncol(x$xdata))
  )
  igraph::V(mygraph)$frame.color <- igraph::V(mygraph)$color
  igraph::V(mygraph)$label.family <- "sans"
  igraph::V(mygraph)$label.cex <- as.numeric(as.character(cut(mydegrees, breaks = 4, labels = c(0.4, 0.45, 0.5, 0.55))))
  igraph::V(mygraph)$label.color <- "grey20"

  # Formatting edges
  igraph::E(mygraph)$color <- "grey60"
  igraph::E(mygraph)$width <- 0.5

  igraph::plot.igraph(mygraph, ...)
}
