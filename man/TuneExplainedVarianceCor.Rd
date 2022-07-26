% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/internal.R
\name{TuneExplainedVarianceCor}
\alias{TuneExplainedVarianceCor}
\title{Tuning function (correlation)}
\usage{
TuneExplainedVarianceCor(u, ev_xx = NULL, omega)
}
\arguments{
\item{u}{constant u added to the diagonal of the precision matrix.}

\item{ev_xx}{desired proportion of explained variance. If \code{ev_xx=NULL},
the obtained proportion of explained variance is returned.}

\item{omega}{positive semidefinite precision matrix.}
}
\value{
The absolute difference in proportion of explained variance (if
  \code{ev_xx} is provided) or observed proportion of explained variance (if
  \code{ev_xx=NULL}).
}
\description{
Computes the difference in absolute value between the desired and observed
proportion of explained variance from the first Principal Component of a
Principal Component Analysis applied on the correlation matrix. The precision
matrix is obtained by adding u to the diagonal of a positive semidefinite
matrix. This function is used to find the value of the constant u that
generates a correlation matrix with desired proportion of explained variance.
}
\keyword{internal}
