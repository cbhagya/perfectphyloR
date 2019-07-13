#' Subset hapMat
#'
#' This function subsets the neighborhood of the focal SNV from the \code{hapMat} object.
#'
#' @param hapMat  A \code{hapMat} object.
#' @param SNVs    Indices of hapMat object specifying the SNVs in the neighborhood.
#'
#' @keywords internal
#'
#' @return A subset of \code{hapMat} object to be returned to \code{\link{selectWindow}}
#' @seealso \code{\link{selectWindow}}
#' 
#' @examples 
#' 
#' \dontshow{
#' 
#' data(ex_hapMatSmall_data)
#' 
#' # Find the neighborhood.
#' fsout = findSNVs(hapMat = ex_hapMatSmall_data, focalSNV = 10, minWindow = 1)
#' 
#' # Subset hapMat to the selected SNV window.
#' subHapMat = subsetHapMat(hapMat = ex_hapMatSmall_data, fsout$SNVwindow)
#' }
#'
subsetHapMat = function(hapMat,SNVs) {
  subhapmat = hapMat$hapmat[,SNVs,drop=FALSE]
  subphenos <- hapMat$phenos[SNVs]
  subposns = hapMat$posns[SNVs]
  subHapMat <- list(hapmat=subhapmat,phenos=subphenos,posns=subposns)
  class(subHapMat) <- "hapMat"
  return(subHapMat)
}
