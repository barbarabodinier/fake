#' Block matrix
#'
#' Generates a symmetric block matrix of size (\code{sum(pk)} x \code{sum(pk)}).
#' The sizes of the submatrices is defined based on \code{pk}. For each
#' submatrix, all entries are equal to the submatrix (block) index.
#'
#' @param pk vector encoding the grouping structure.
#'
#' @return A symmetric block matrix.
#'
#' @family block matrix functions
#'
#' @examples
#' # Example 1
#' BlockMatrix(pk = c(2, 3))
#'
#' # Example 2
#' BlockMatrix(pk = c(2, 3, 2))
#' @export
BlockMatrix <- function(pk) {
  nblocks <- sum(upper.tri(matrix(NA, ncol = length(pk), nrow = length(pk)), diag = TRUE))
  blocks <- matrix(NA, nrow = length(pk), ncol = length(pk))
  blocks[upper.tri(blocks, diag = TRUE)] <- 1:nblocks

  mybreaks <- c(0, cumsum(pk))
  bigblocks <- matrix(ncol = sum(pk), nrow = sum(pk))
  row_id_start <- matrix(mybreaks[row(blocks)], ncol = length(pk)) + 1
  row_id_end <- matrix(mybreaks[row(blocks) + 1], ncol = length(pk))
  col_id_start <- matrix(mybreaks[col(blocks)], ncol = length(pk)) + 1
  col_id_end <- matrix(mybreaks[col(blocks) + 1], ncol = length(pk))

  row_id_start <- row_id_start[upper.tri(row_id_start, diag = TRUE)]
  row_id_end <- row_id_end[upper.tri(row_id_end, diag = TRUE)]
  col_id_start <- col_id_start[upper.tri(col_id_start, diag = TRUE)]
  col_id_end <- col_id_end[upper.tri(col_id_end, diag = TRUE)]

  for (block_id in blocks[upper.tri(blocks, diag = TRUE)]) {
    ids <- rbind(
      expand.grid(
        row_id_start[block_id]:row_id_end[block_id],
        col_id_start[block_id]:col_id_end[block_id]
      ),
      expand.grid(
        col_id_start[block_id]:col_id_end[block_id],
        row_id_start[block_id]:row_id_end[block_id]
      )
    )
    bigblocks[as.matrix(ids)] <- block_id
  }

  return(bigblocks)
}


#' Block diagonal matrix
#'
#' Generates a binary block diagonal matrix.
#'
#' @param pk vector encoding the grouping structure.
#'
#' @return A binary block diagonal matrix.
#'
#' @family block matrix functions
#'
#' @examples
#' # Example 1
#' BlockDiagonal(pk = c(2, 3))
#'
#' # Example 2
#' BlockDiagonal(pk = c(2, 3, 2))
#' @export
BlockDiagonal <- function(pk) {
  bigblocks <- BlockMatrix(pk)
  bigblocks[!bigblocks %in% diag(bigblocks)] <- 0
  bigblocks[bigblocks %in% diag(bigblocks)] <- 1
  return(bigblocks)
}


#' Block structure
#'
#' Generates a symmetric matrix of size (\code{length(pk)} x \code{length(pk)})
#' where entries correspond to block indices. This function can be used to
#' visualise block indices of a matrix generated with \code{\link{BlockMatrix}}.
#'
#' @inheritParams BlockMatrix
#'
#' @return A symmetric matrix of size \code{length(pk))}.
#'
#' @family block matrix functions
#'
#' @examples
#' # Example 1
#' BlockMatrix(pk = c(2, 3))
#' BlockStructure(pk = c(2, 3))
#'
#' # Example 2
#' BlockMatrix(pk = c(2, 3, 2))
#' BlockStructure(pk = c(2, 3, 2))
#' @export
BlockStructure <- function(pk) {
  blocks <- BlockMatrix(pk = rep(1, length(pk)))

  return(blocks)
}
