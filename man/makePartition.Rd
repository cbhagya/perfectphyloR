% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/make_dend_utility_functions.R
\name{makePartition}
\alias{makePartition}
\title{Recursively partition haplotype matrix}
\usage{
makePartition(hapmat, splitSNV)
}
\arguments{
\item{hapmat}{A hapMat object with SNVs ordered according to ancestry.}

\item{splitSNV}{The index of the SNV where the haplotype matrix from \code{\link{orderSNVs}} is partitioned.}
}
\value{
A nested partition of haplotypes, implemented as a list of nodes, each with two child nodes.
}
\description{
This function recursively partitions the SNVs in the window around the focal SNV.
}
\details{
This function makes two clades based on \code{splotSNV}. For each partition, update \code{splitSNV} and use
it to define subclades. Then, continue recursive partitioning until each partition has only one haplotype,
or there are no more SNVs to consider.
}
\examples{

\dontshow{

data(ex_hapMatSmall_data)

# First, select a window of SNVs about a focal SNV.
SNV_win <- selectWindow(hapMat = ex_hapMatSmall_data,
                        focalSNV = 10, minWindow = 1)
                        
# Then order SNVs in the window.
ordHapmat <- orderSNVs(snvWin = SNV_win)

# Recursively partition haplotype matrix.
partitions <- makePartition(hapmat = ordHapmat, splitSNV = 1)

}

}
\seealso{
\code{\link{makeDend}}, \code{\link{newNode}}, \code{\link{noVariation}}
}
\keyword{internal}
