
#' Order SNVs
#'
#' This function orders the SNVs in the window around the SNV.
#'
#' Following Gusfield (1991), the first step in building the dendogram is to order compatible SNVs by
#' their ancestry. Then, following Mailund et al. (2006) the incompatible SNVs are ordered according
#' to their proximity to the focal SNV.
#' @param snvWin See the arguments in \code{\link{buildDend}}
#'
#' @return A hapMat with ordered SNVs.
#'
#' @references  Gusfield, D. (1991). Efficient algorithms for inferring evolutionary trees.
#'              Networks, 21(1), 19-28.
#' @references  Mailund, T., Besenbacher, S., & Schierup, M. H. (2006). Whole genome association
#'              mapping by incompatibilities and local perfect phylogenies. BMC Bioinformatics, 7(1),
#'              454.
#' @keywords internal
#' @seealso \code{\link{buildDend}}, \code{\link{orderColsAncestry}}
#'
#' @examples 
#' 
#' \dontshow{
#' 
#' data(ex_hapMatSmall_data)
#' 
#' # First select a window of SNVs about a focal SNV.
#' 
#' SNV_win <- selectWindow(hapMat = ex_hapMatSmall_data,
#'                        focalSNV = 10, minWindow = 1)
#'                        
#' # Then, order the selected window of SNVs around the focal SNV.
#' 
#'  ordSNVwin <- orderSNVs(snvWin = SNV_win)
#' 
#' }
#' 
orderSNVs = function(snvWin) {
  # Input is a list with the following components.
  # - hapMat: the haplotype matrix data structure,
  # - focalSNV: the index of the focal SNV position
  #   within the hapMat$posns vector,
  # - compat: the indicator vector of whether
  #   or not each SNV is compatible with the focal SNV.
  hapMat <- snvWin$hapMat
  focalSNV <- snvWin$focalSNV
  compat <- snvWin$compat
  # Output:
  # - Just the haplotype matrix itself, *not* a hapMat
  #   data structure
  if(focalSNV > length(hapMat$posns)) stop("Index of focal SNV too large")
  # Order compatible SNVs by ancestry, with more ancestral before
  # more recent mutations.
  compatSNVs = hapMat$hapmat[, compat, drop = FALSE]
  if(ncol(compatSNVs) > 1) {
    ord = orderColsAncestry(compatSNVs)
    compatSNVs = compatSNVs[, ord]
  }
  hapmat = compatSNVs
  # Order incompatible SNVs by proximity to focal SNV(i.e. order incompatible SNVs according to the distance
  # from the focal SNV).
  distFromFocal = abs(hapMat$posns - hapMat$posns[focalSNV])
  if(any(!compat)) {
    incompatSNVs = hapMat$hapmat[, !compat, drop = FALSE]
    ord = order(distFromFocal[!compat])
    incompatSNVs = incompatSNVs[, ord, drop = FALSE]
    hapmat = cbind(hapmat, incompatSNVs)
  }
  return(hapmat)
}
