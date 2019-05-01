#' Example small dataset
#'
#' A subset of \code{ex_hapMat_data}, containing 10 sequences (haplotypes) with 20 SNVs.
#'
#' @docType data
#'
#' @usage data(ex_hapMatSmall_data)
#'
#'
#'
#' @format A list of ten haplotypes with the physical positions of each SNV.
#' \describe{
#' \item{hapmat}{A matrix of 0's and 1's, with rows representing haplotypes and columns representing SNVs.}
#' \item{snvNames}{A vector of names of SNVs for the columns of \code{hapmat}.}
#' \item{hapNames}{A vector of names of haplotypes for the rows of \code{hapmat}.}
#' \item{posns}{ a numeric vector specifying the genomic positions (e.g. in base pairs) of SNVs in the
#'  columns of \code{hapmat}.}
#' }
#'
"ex_hapMatSmall_data"
