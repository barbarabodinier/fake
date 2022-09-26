#' Heatmap visualisation
#'
#' Produces a heatmap for visualisation of matrix entries.
#'
#' @param mat data matrix.
#' @param col vector of colours.
#' @param resolution number of different colours to use.
#' @param bty character string indicating if the box around the plot should be
#'   drawn. Possible values include: \code{"o"} (default, the box is drawn), or
#'   \code{"n"} (no box).
#' @param axes logical indicating if the row and column names of \code{mat}
#'   should be displayed.
#' @param cex.axis font size for axes.
#' @param xlas orientation of labels on the x-axis, as \code{las} in
#'   \code{\link[graphics]{par}}.
#' @param ylas orientation of labels on the y-axis, as \code{las} in
#'   \code{\link[graphics]{par}}.
#' @param text logical indicating if numbers should be displayed.
#' @param cex font size for numbers. Only used if \code{text=TRUE}.
#' @param legend logical indicating if the colour bar should be included.
#' @param legend_length length of the colour bar.
#' @param legend_range range of the colour bar.
#' @param cex.legend font size for legend.
#' @param ... additional arguments passed to \code{\link[base]{formatC}} for
#'   number formatting. Only used if \code{text=TRUE}.
#'
#' @return A heatmap.
#'
#' @examples
#' oldpar <- par(no.readonly = TRUE)
#' par(mar = c(3, 3, 1, 5))
#'
#' # Data simulation
#' set.seed(1)
#' mat <- matrix(rnorm(100), ncol = 10)
#' rownames(mat) <- paste0("r", 1:nrow(mat))
#' colnames(mat) <- paste0("c", 1:ncol(mat))
#'
#' # Generating heatmaps
#' Heatmap(mat = mat)
#' Heatmap(mat = mat, text = TRUE, format = "f", digits = 2)
#' Heatmap(
#'   mat = mat,
#'   col = c("lightgrey", "blue", "black"),
#'   legend = FALSE
#' )
#'
#' par(oldpar)
#' @export
Heatmap <- function(mat, col = c("ivory", "navajowhite", "tomato", "darkred"),
                    resolution = 10000, bty = "o",
                    axes = TRUE, cex.axis = 1, xlas = 2, ylas = 2,
                    text = FALSE, cex = 1,
                    legend = TRUE, legend_length = NULL, legend_range = NULL, cex.legend = 1, ...) {
  oldpar <- graphics::par("xpd", "xaxs", "yaxs", no.readonly = TRUE)
  on.exit(graphics::par(oldpar))

  # Transposing the input matrix so that rows are rows
  mat <- t(mat)

  # Defining the legend length
  if (is.null(legend_length)) {
    legend_length <- ncol(mat)
  }

  # Preparing colours
  col <- grDevices::colorRampPalette(col)(resolution)
  names(col) <- 1:resolution

  # Re-formatting matrix
  mat <- mat[, ncol(mat):1, drop = FALSE]
  vect <- as.vector(mat)

  # Defining extreme values
  if (is.null(legend_range)) {
    # myrange <- c(min(vect, na.rm = TRUE), max(vect, na.rm = TRUE))
    myrange <- range(vect, na.rm = TRUE)
    myrange <- c(floor(myrange[1]), ceiling(myrange[2]))
  } else {
    myrange <- legend_range
  }

  # Getting corresponding colours
  mycol <- as.character(cut(vect, breaks = seq(myrange[1], myrange[2], length.out = resolution + 1), labels = 1:resolution, include.lowest = TRUE))
  mycol_mat <- matrix(mycol, ncol = ncol(mat))

  # Making heatmap
  withr::local_par(xaxs = "i", yaxs = "i")
  plot(NA,
    xlim = c(0, nrow(mycol_mat)), ylim = c(0, ncol(mycol_mat)),
    xlab = "", ylab = "", xaxt = "n", yaxt = "n", bty = "n"
  )
  for (i in 0:(nrow(mycol_mat) - 1)) {
    for (j in 0:(ncol(mycol_mat) - 1)) {
      graphics::polygon(
        x = c(i, i + 1, i + 1, i), y = c(j, j, j + 1, j + 1),
        col = col[mycol_mat[i + 1, j + 1]],
        border = col[mycol_mat[i + 1, j + 1]]
      )
    }
  }
  if (axes) {
    if (!is.null(rownames(mat))) {
      graphics::axis(side = 1, at = 1:nrow(mat) - 0.5, labels = rownames(mat), las = xlas, cex.axis = cex.axis)
    }
    if (!is.null(colnames(mat))) {
      graphics::axis(side = 2, at = 1:ncol(mat) - 0.5, labels = colnames(mat), las = ylas, cex.axis = cex.axis)
    }
  }
  if (bty == "o") {
    graphics::box()
  }

  # Showing text
  if (text) {
    for (i in 1:nrow(mat)) {
      for (j in 1:ncol(mat)) {
        text(i - 0.5, j - 0.5,
          cex = cex,
          labels = formatC(mat[i, j], ...)
        )
      }
    }
  }

  # Adding colour bar (legend)
  if (legend) {
    withr::local_par(list(xpd = TRUE))
    legend_width_factor <- 1.05
    mylegend_values <- grDevices::axisTicks(c(myrange[1], myrange[2]), log = FALSE)
    mylegend_ids <- as.numeric(as.character(cut(mylegend_values,
      breaks = seq(myrange[1], myrange[2], length.out = resolution + 1),
      labels = 1:resolution, include.lowest = TRUE
    )))
    ypos <- ncol(mat)
    xpos <- nrow(mat) * 1.05
    for (l in 1:length(col)) {
      graphics::polygon(
        x = c(xpos, xpos * legend_width_factor, xpos * legend_width_factor, xpos),
        y = c(
          ypos - legend_length + legend_length * l / length(col),
          ypos - legend_length + legend_length * l / length(col),
          ypos - legend_length + legend_length * (l + 1) / length(col),
          ypos - legend_length + legend_length * (l + 1) / length(col)
        ),
        col = col[l], border = col[l]
      )
      if (l %in% mylegend_ids) {
        graphics::text(
          x = xpos * legend_width_factor, y = ypos - legend_length + legend_length * (l + 0.5) / length(col),
          labels = paste0("- ", mylegend_values[which(mylegend_ids == l)]), adj = c(0, 0.5), cex = cex.legend
        )
      }
    }
    withr::local_par(list(xpd = FALSE)) # for legend
  }
}
