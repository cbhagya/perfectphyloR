% RVcoeff from Rcpp
\name{RVcoeff}
\alias{RVcoeff}
\title{Compute RV coefficient to measure association of two distance matrices}
\usage{
RVcoeff(mDx, mDy, mC)
}
\arguments{
\item{mDx}{A numeric matrix of pairwise distances.}

\item{mDy}{A second numeric matrix of pairwise distances.}

\item{mC}{See the equation 2.4 in the reference.}
}

\description{
This function computes the RV coefficient of two distance matrices.
}

\keyword{internal}

\references{

Josse, J., & Holmes, S. (2016). Measuring multivariate association and beyond. Statistics surveys, 10, 132.
}