
#' Apply Four-Gamete Test to check the compatibility of a pair of SNVs.
#'
#' This function applies the four-gamete test to the pair of SNVs defined next SNV with each SNV in current neighborhood,
#' to check the compatability.
#'
#' @param current Current SNVs in the neighborhood. This starts from the focal SNV.
#' @param nextSNV A SNV whose compatability with each SNV in the current neighborhood needs to be checked.
#' @param hapmat  A matrix of 0's and 1's, with rows representing haplotypes and columns representing SNVs.
#'
#' @return A logical vector of the same length as \code{current} with entries being \code{TRUE} if next SNV is compatible
#'         with the corresponding SNV in \code{current} and \code{FALSE} other wise;returned to \code{\link{findSNVs}}
#'
#' @keywords internal
#'
#' @seealso \code{\link{findSNVs}}, \code{\link{fourGamete}}
#'
checkCompatible = function(current, nextSNV, hapmat) {
  # Apply four-gamete test to next SNV paired with each SNV in current.
  
  for(i in current[1]:current[2]) {
    if(fourGamete(hapmat[, i], hapmat[, nextSNV])) return(FALSE)
  }
  return(TRUE)
}
#' Four-Gamete Test
#'
#' This function performs the Four-Gamete Test to check the compatability of two SNVs. If all four haplotypes
#' (00,10,01,11) are observed for two SNVs, those two SNVs are incompatible.
#'
#' @param snv1 First SNV.
#' @param snv2 Second SNV.
#'
#' @return TRUE, if both SNVs are incompatible.
#' @keywords internal
#'
fourGamete = function(snv1,snv2) {
  tab = table(snv1,snv2)
  if(all(tab>0)) return(TRUE) else return(FALSE)
}

#' Number of SNVs
#'
#' This function gives the number of SNVs contained within a given pair of indices of SNVs, inclusively.
#'
#' @param cur A numeric vector of two elements indexing the first and last SNV.
#'
#' @return A numeric value specifying the number of SNVs returned to \code{\link{findSNVs}}
#' @keywords internal
#' @seealso \code{\link{findSNVs}}
#'
#' @examples
#' SNV_index = c(2,5)
#' getnSNVs(SNV_index) # return 4. i.e. number of SNVs from SNV2 to SNV5, inclusively is 4.
#'
#'
getnSNVs = function(cur) {
  # return number of SNVs from indices cur[1] to cur[2], inclusively.
  
  return(cur[2]-cur[1]+1)
}

#' Next SNV from the focal point
#'
#' This function finds the next SNV to be checked for compatability with the existing window of SNVs.
#'
#' @param current Current SNVs in the neighborhood of focal SNV.
#' @param absDist Absolute distance of all SNVs from the focal SNV in base pairs.
#'
#' @return The index of the SNV that is nearest to the focal SNV but outside the current neighborhood.
#'         To be returned to \code{\link{findSNVs}}.
#' @keywords internal
#' @seealso \code{\link{findSNVs}}
#'
getNextFromFocal = function(current, absDist) {
  # Disqualify distances within current range by setting to Inf.
  
  absDist[current[1]:current[2]] = Inf
  
  # Return index of minimum over remaining elements. In case of
  # a tie, which.min() returns the first index found.
  
  return(which.min(absDist))
}
#' Expand the neighborhood to the left
#'
#' This function tries to expand the neighborhood to the left, when the next SNV
#' from \code{\link{getNextFromFocal}} is incompatible with each SNV in
#' the neighborhood of the focal SNV, and this incompatible SNV is to the
#' right of the neighborhood.
#'
#' @param current Current SNVs in the neighborhood of focal SNV.
#' @param absDist Absolute distance of SNVs from the focal SNV in base pairs.
#'
#' @return An index of a SNV to the left of the neighborhood to be returned to \code{\link{findSNVs}}.
#' @keywords internal
#'
#' @seealso \code{\link{findSNVs}}
#'
getNextLeftFocal = function(current,absDist) {
  if(current[1]==1) {
    # There are no markers to the left, return NA.
    return(NA)
  }
  # Disqualify distances in current range and to the right of current.
  
  absDist[current[1]:length(absDist)] = Inf
  
  # Return index of minimum over remaining elements.
  return(which.min(absDist))
}
#' Expand the neighborhood to the right
#'
#' This function tries to expand the neighborhood to the right, when the next SNV
#' from \code{\link{getNextFromFocal}} function is incompatible with each SNV in
#' the neighborhood of the focal SNV, and this incompatible SNV could be to the
#' left of the neighborhood.
#'
#' @param current Current SNVs in the neighborhood of focal SNV.
#' @param absDist Absolute distance of SNVs from the focal SNV in base pairs.
#'
#' @return An index of a SNV.
#' @keywords internal
#'
getNextRightFocal = function(current,absDist) {
  if(current[2]==length(absDist)) {
    # There are no markers to the right, return NA.
    
    return(NA)
  }
  # Disqualify distances in current range and to the left of current.
  absDist[1:current[2]] = Inf
  # Return index of minimum over remaining elements.
  return(which.min(absDist))
}
