#' Select a window of SNVs about the focal SNV.
#'
#' This internal function first identifies as many compatible SNVs as possible around the focal SNV. If the
#' neighborhood of compatible SNVs is smaller than a user-defined minimum number of SNVs, this function expands
#' the neighborhood by including incompatible SNVs in order of proximity to the focal SNV using the algorithm
#' of Mailund et al. (2006). Then the function subsets the columns of the \code{hapMat} data object according
#' to the resulting SNV window.
#'
#'
#'
#'
#' @param hapMat   \code{hapMat} data object.
#' @param focalSNV The column number of the focal SNV in the \code{hapMat} data object.
#' @param minWindow Minimum window size of the SNV neighborhood.
#'
#'
#'
#' @keywords internal
#'
#' @seealso \code{\link{reconstructPP}}, \code{\link{findSNVs}}, \code{\link{subsetHapMat}}
#'
#' @references  Mailund, T., Besenbacher, S., & Schierup, M. H. (2006). Whole genome association mapping
#'              by incompatibilities and local perfect phylogenies. BMC Bioinformatics, 7(1), 454.
#'              
#'
selectWindow <- function(hapMat, focalSNV, minWindow) {
  # Find the neighborhood.
  fsout = findSNVs(hapMat, focalSNV, minWindow)
  # Now subset hapMat to the selected SNV window.
  subHapMat = subsetHapMat(hapMat, fsout$SNVwindow)
  # Identify the index of the focalSNV in the subsetted data structure.
  focalSNV = which(fsout$SNVwindow == focalSNV)
  # Combine the output into a single list and return.
  out <- list(hapMat = subHapMat, focalSNV = focalSNV, compat = fsout$compat)
  return(out)
}
