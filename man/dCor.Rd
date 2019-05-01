% dCorcoeff from Rcpp
\name{dCor}
\alias{dCor}
\title{Compute distance correlation (dCor) coefficient two random vectors(distance matrices).}
\usage{
dCor(mDx, mDy, mC)
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