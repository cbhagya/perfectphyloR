#' Example dataset
#'
#' A \code{hapMat} data object containing 200 sequences (haplotypes) with 2747 SNVs.
#'
#' @docType data
#'
#' @usage data(ex_hapMat_data)
#'
#'
#'
#' @format A list of 200 haplotypes with the physical positions of each SNV.
#' \describe{
#' \item{hapmat}{A matrix of 0's and 1's, with rows representing haplotypes and columns representing SNVs.}
#' \item{snvNames}{A vector of names of SNVs for the columns of \code{hapmat}.}
#' \item{hapNames}{A vector of names of haplotypes for the rows of \code{hapmat}.}
#' \item{posns}{ a numeric vector specifying the genomic positions (e.g. in base pairs) of SNVs in the
#'  columns of \code{hapmat}.}
#' }
#'
"ex_hapMat_data"
