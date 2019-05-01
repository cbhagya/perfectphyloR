#' Create an object of class \code{hapMat}
#'
#' This function creates a \code{hapMat} data object, a required input for \code{\link{reconstructPP}}.
#'
#'
#' @param hapmat   A matrix of 0's and 1's, with rows representing haplotypes and columns representing 
#'                 single-nucleotide variants (SNVs).
#' @param snvNames A vector of names of SNVs for the columns of \code{hapmat}.
#' @param hapNames A vector of names of haplotypes for the rows of \code{hapmat}.
#' @param posns    A numeric vector specifying the genomic positions (e.g. in base pairs) of SNVs in the
#'                 columns of \code{hapmat}.
#'
#' @return An object of \code{hapMat}.
#' @export
#'
#' @examples
#' hapmat = matrix(c(1,1,1,0,
#'                   0,0,0,0,
#'                   1,1,1,1,
#'                   1,0,0,0,
#'                   1,1,0,0,
#'                   1,0,0,1,
#'                   1,0,0,1), byrow = TRUE, ncol = 4)
#' snvnames = c(paste("SNV", 1:4, sep = ""))
#' allhaps = c("h1", "h2", "h3", "h4", "h5", "h6", "h7")
#' # Physical positions
#' posns = c(1000, 2000, 3000, 4000)
#'
#' # Create hapMat data object
#' ex_hapMat <- createHapMat(hapmat = hapmat,
#'                                snvNames = snvnames,
#'                                hapNames = allhaps,
#'                                posns = posns)
#'
#'
createHapMat = function(hapmat, snvNames, hapNames, posns){

  if(length(snvNames)!= ncol(hapmat))
    stop("Number of SNV names not equal to number of columns of haplotype matrix")
  if(length(hapNames)!= nrow(hapmat))
    stop("Number of haplotype names not equal to number of rows of haplotype matrix")
  # Set up the data structure
  colnames(hapmat) = snvNames
  rownames(hapmat) = hapNames
  hapMat = list(hapmat = hapmat, posns = posns)
  class(hapMat) = "hapMat"
  return(hapMat)
}
