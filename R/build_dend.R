#' Build the tree for the window of SNVs.
#'
#' This function builds the perfect phylogeny for the window of SNVs which is selected from function
#' \code{\link{selectWindow}}.
#'
#' This function works in two stages. First, it orders the SNVs in the window, based on age
#' for compatible SNVs and proximity to the focal SNV for incompatible SNVs.
#' Then, it makes the perfect phylogeny for the ordered SNVs using recursive partitioning and records
#' the partitioning.
#'
#' @param snvWin A list contains these three components: (1) hapMat: the  data structure summerizing
#'               the SNV window.
#'               (2) focalSNV: the column number of the focal SNV at which to reconstruct the
#'               perfect phylogeny.
#'               (3) compat: the local vector of whether or not each SNV in the window
#'               is compatible with the focal SNV.
#'
#' @param sep A character string separator for concatenating haplotype labels in the
#'            dendrogram if they are undistingushable in the window around the focal SNV
#'            , see the arguments in \code{\link{reconstructPP}}.
#'
#' @return An object of class \code{phylo}.  
#' @keywords internal
#' @seealso \code{\link{reconstructPP}}, \code{\link{orderSNVs}}, \code{\link{makeDend}}
#'
buildDend <- function(snvWin, sep) {
  # (a) Order SNVs in the window based on (i) age for compatible
  # SNVs and (ii) proximity to the focal SNV for incompatible SNVs.
  ordHapmat <- orderSNVs(snvWin)
  # (b) Make the tree for the ordered SNVs using the function
  # makeDend().
  tree <- makeDend(ordHapmat, sep = sep)
  return(tree)
}
