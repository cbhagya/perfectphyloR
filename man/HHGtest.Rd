% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/HHG_test.R
\name{HHGtest}
\alias{HHGtest}
\title{HHG test for association of two distance matrices}
\usage{
HHGtest(Dx, Dy, nperm)
}
\arguments{
\item{Dx}{A numeric matrix of pairwise distances.}

\item{Dy}{A second numeric matrix of pairwise distances.}

\item{nperm}{The number of times to permute  the rows and columns of \code{Dy}.}
}
\value{
A list contains HHG coefficient and permutation P value.
}
\description{
This function performs HHG test to find the association between two distance matrices. It permutes rows and columns
of the second matrix randomly to calculate P value.
}
\examples{

x <- runif(8)
y <- runif(8)
# Distance matrices
distX = as.matrix(dist(x, upper = TRUE, diag = TRUE))
distY = as.matrix(dist(y, upper = TRUE, diag = TRUE))

HHGtest(Dx = distX, Dy = distY, nperm = 1000) 

}
\references{
Barak, B., and Shachar, K., based in part on an earlier implementation by Ruth Heller 
and Yair Heller. (2017). HHG: Heller-Heller-Gorfine Tests of Independence and Equality of Distributions. R
package version 2.2. https://CRAN.R-project.org/package=HHG
}
