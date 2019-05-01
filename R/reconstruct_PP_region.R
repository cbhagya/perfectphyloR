
#' Reconstruct perfect phylogeny sequencce across a region
#'
#' This function reconstructs perfect phylogenies on each possible focal SNV across a genomic region.
#'
#' @param hapMat A data structure of class \code{hapMat}. See the arguments in \code{\link{reconstructPP}}
#' @param minWindow Minimum number of SNVs around the focal SNV in the window of SNVs used to
#'                  reconstruct the perfect phylogeny.
#'
#' @return An object of class \code{multiPhylo} that contains multiple \code{phylo} objects.
#' @export
#'
#' @examples
#'   
#' data(ex_hapMatSmall_data)                      
#' # Reconstruct partitions across the region of ex_hapMatSmall_data.
#' rdends <- reconsPPregion(hapMat = ex_hapMatSmall_data,
#'                       minWindow = 1)
#'
#'
reconsPPregion = function(hapMat,minWindow) {

  nSNV = ncol(hapMat$hapmat)
  treeVec = vector(length = nSNV, mode = "list")
  class(treeVec) = "multiPhylo"

  for(i in 1:nSNV) {
    treeVec[[i]] = reconstructPP(hapMat, focalSNV = i, minWindow)
  }
  return(treeVec)
}
