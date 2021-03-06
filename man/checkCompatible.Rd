% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/find_SNV_utility_functions.R
\name{checkCompatible}
\alias{checkCompatible}
\title{Apply Four-Gamete Test to check the compatibility of a pair of SNVs.}
\usage{
checkCompatible(current, nextSNV, hapmat)
}
\arguments{
\item{current}{Current SNVs in the neighborhood. This starts from the focal SNV.}

\item{nextSNV}{A SNV whose compatability with each SNV in the current neighborhood needs to be checked.}

\item{hapmat}{A matrix of 0's and 1's, with rows representing haplotypes and columns representing SNVs.}
}
\value{
A logical vector of the same length as \code{current} with entries being \code{TRUE} if next SNV is compatible
        with the corresponding SNV in \code{current} and \code{FALSE} other wise;returned to \code{\link{findSNVs}}
}
\description{
This function applies the four-gamete test to the pair of SNVs defined next SNV with each SNV in current neighborhood,
to check the compatability.
}
\examples{

\dontshow{

data(ex_hapMatSmall_data)

# Check the compatibility of a pair of SNVs using Four-Gamete Test.

comptbility <- checkCompatible(current = c(10,10), nextSNV = 11, hapmat = ex_hapMatSmall_data$hapmat)

}

}
\seealso{
\code{\link{findSNVs}}, \code{\link{fourGamete}}
}
\keyword{internal}
