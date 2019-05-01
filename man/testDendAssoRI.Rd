
\name{testDendAssoRI}
\alias{testDendAssoRI}
\title{Tests Rand index between a comparator dendrogram and reconstructed dendrograms}
\usage{
testDendAssoRI(rdend, cdend, hapMat, k = 2 , nperm = 0, xlab = "", ylab = "", main = "")
}
\arguments{
\item{rdend}{A multiphylo object of reconstructed dendrograms at each focal SNV.}

\item{cdend}{A phylo object of the comparator dendrogram.}

\item{hapMat}{An object of class \code{hapMat} containing SNV haplotypes.}


\item{k}{An integer that specifies the number of clusters that the dendrogram should be cut into. The              default is k=2. Clusters are defined by starting from the root of the dendrogram and moving towards the tips, cutting horizontally at any given point in the dendrogram.}

\item{nperm}{Number of permutations for the test of any association across the genomic region of interest.  The              default is \code{nperm = 0}; i.e., association will not be tested.}

\item{xlab}{An optional character string for the label on the x-axis in the plot that is returned (none by            default).} 

\item{ylab}{An optional character string for the label on the y-axis in the plot that is returned (none by            default).}

\item{main}{An optional character string for title in the plot that is returned (none by default).}

}
\value{
A list with the following components: 

\item{Stats}{A vector of observed Rand indices.}
\item{OmPval}{A permutation-based omnibus P value for the test of any association across the genomic                    region using the maximum Rand index over the genomic region as the test statistic.}
\item{mPval}{A vector of marginal P values at each SNV position}
\item{plt}{A plot of the association profile of the Rand indices over SNV locations in the region of interest.}

}
\description{
This function calculates and tests the association between a user-supplied comparator dendrogram and a list of reconstructed dendrograms across a genomic region of interest using the Rand index.}

\examples{
   See the section 'Applications' in vignette("perfectphyloR") for the detailed example.

}
\section{See Also}{

 \code{\link{RandIndexTest}}
}