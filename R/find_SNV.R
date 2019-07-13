#' Find the window of SNVs at a focal point
#'
#' This function identifies a window of compatible SNVs about the focal SNV.
#'
#' This function starts at the focal SNV, and tries expanding the neighborhood of compatible SNVs.
#' This process is stopped when there are no more compatible SNVs to the left or right
#' of the focal SNV. If the block of compatible SNVs does not have at least a user-defined minimum number
#' of SNVs, this function expands the window to include incompatible SNVs.
#'
#' @param hapMat    \code{hapMat} data object.
#' @param focalSNV  The column number of the focal in the haplotype data matrix
#'                  at which to reconstruct the dendrogram.
#' @param minWindow Minimum window size of the SNV neighborhood.
#'
#' @return Return a list of SNVs in the window with their state of compatiblity to \code{\link{selectWindow}}.
#'
#' @keywords internal
#'
#' @seealso \code{\link{selectWindow}}, \code{\link{checkCompatible}}, \code{\link{getnSNVs}},
#'          \code{\link{getNextFromFocal}}, \code{\link{getNextLeftFocal}}, \code{\link{getNextRightFocal}}
#'          
#'
#' @examples 
#' 
#' \dontshow{
#' data(ex_hapMatSmall_data)
#' 
#' # Compatible SNVs about the focal SNV
#' 
#' compSNVs <- findSNVs(hapMat = ex_hapMatSmall_data, focalSNV = 10, minWindow = 1)
#' 
#' }                                         
#'
findSNVs = function(hapMat, focalSNV, minWindow = 1) {

  # 0. Error checking and initialization
  if(minWindow > ncol(hapMat$hapmat))
    stop("minWindow must not exceed number of SNVs")
  # Proximity to the focal SNV.
  absDistFromFocal = abs(hapMat$posns - hapMat$posns[focalSNV]) 
  # Initialize the window of compatible SNVs to just the focal SNV.
  current = c(focalSNV,focalSNV)
  nextSNV = getNextFromFocal(current,absDistFromFocal)
  ###################################
  # 1. Look for compatible SNVs in order of proximity to focal.
  while(checkCompatible(current, nextSNV, hapMat$hapmat) &&
        getnSNVs(current) < ncol(hapMat$hapmat)) {
    current = c(min(nextSNV, current[1]), max(nextSNV, current[2]))
    nextSNV = getNextFromFocal(current, absDistFromFocal)
  }
  if(getnSNVs(current) == ncol(hapMat$hapmat)) {
    # All SNVs are compatible. Return indices of all SNVs.
    SNVwindow = current[1]:current[2]
    compat = rep(TRUE, ncol(hapMat$hapmat))
    return(list(SNVwindow = SNVwindow, compat = compat))
  }
  ###################################
  # 2. At this point, nextSNV holds the index of the nearest SNV that
  # is incompatible with those in the neighborhood of the focal SNV.
  # This incompatible SNV could be to the left of the neighborhood
  # (nextSNV<current[1]) or else to the right.
  # 2a. If to the left, try to expand the neighborhood to the right.
  # 2b. If to the right, try to expand the neighborhood to the left.
  if(nextSNV < current[1]) { 
    # Incompatible on left, extend right.
    nextSNV = getNextRightFocal(current,absDistFromFocal)
    # nextRightFocal returns index to right of current range, or NA if
    # there are no SNVs to the right of current.
    while(!is.na(nextSNV) && checkCompatible(current, nextSNV, hapMat$hapmat)) {
      current[2] = nextSNV
      nextSNV = getNextRightFocal(current, absDistFromFocal)
    }
  } else { 
    # Incompatible on right, extend left.
    nextSNV = getNextLeftFocal(current, absDistFromFocal)
    # nextLeftFocal returns index to left of current range, or NA if
    # there are no SNVs to the left of current.
    while(!is.na(nextSNV) && checkCompatible(current, nextSNV, hapMat$hapmat)){
      current[1] = nextSNV
      nextSNV = getNextLeftFocal(current,absDistFromFocal)
    }
  }
  ###################################
  compatInds = (current[1]:current[2])
  # 3. The indices (current[1],current[2]) are a neighborhood of
  # compatible SNVs. If this is not a
  # big enough neighborhood, expand it by including incompatible
  # SNVs in order of proximity to the focal SNV.
  while(getnSNVs(current)<minWindow) {
    nextSNV = getNextFromFocal(current, absDistFromFocal)
    current = c(min(nextSNV, current[1]), max(nextSNV, current[2]))
  }
  # Return window of SNVs as a vector of indices of SNVs in hapMat.
  SNVwindow = current[1]:current[2]
  compat = SNVwindow %in% compatInds
  return(list(SNVwindow = SNVwindow, compat = compat))
}
