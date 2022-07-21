#' Matching arguments
#'
#' Returns a vector of overlapping character strings between \code{extra_args}
#' and arguments from function \code{FUN}. If \code{FUN} is taking \code{...} as
#' input, this function returns \code{extra_args}.
#'
#' @param extra_args vector of character strings.
#' @param FUN function.
#'
#' @return A vector of overlapping arguments.
#'
#' @examples
#' MatchingArguments(
#'   extra_args = list(Sigma = 1, test = FALSE),
#'   FUN = MASS::mvrnorm
#' )
#' @export
MatchingArguments <- function(extra_args, FUN) {
  if ("..." %in% names(formals(FUN))) {
    out <- extra_args
  } else {
    ids <- which(names(extra_args) %in% names(formals(FUN)))
    out <- extra_args[ids]
  }
  return(out)
}
