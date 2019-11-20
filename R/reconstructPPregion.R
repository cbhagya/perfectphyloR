
#' Reconstruct perfect phylogeny sequencce across a region
#' 
#' This function reconstructs perfect phylogenies on each possible focal SNV across a genomic region.
#'
#' 
#' 
#' @param hapMat    A data structure of class \code{hapMat}. See the arguments in \code{\link{reconstructPP}}.
#' @param minWindow Minimum number of SNVs around the focal SNV in the window of SNVs used to
#'                  reconstruct the perfect phylogeny.
#' @param posn.lb   Lower bound of the subregion of \code{hapMat} (in base pairs) within which to consider SNVs.
#' @param posn.ub   Upper bound of the subregion of \code{hapMat} (in base pairs) within which to consider SNVs.
#'
#' @return  An object of class \code{multiPhylo} that contains multiple \code{phylo} objects.
#' @export
#'
#' @examples
#' 
#' data(ex_hapMatSmall_data)   
#'                    
#' # Reconstruct partitions across the region of ex_hapMatSmall_data.
#' rdends <- reconstructPPregion(hapMat = ex_hapMatSmall_data,
#'                       minWindow = 1)
#'                       
#' # Reconstruct partitions between a given range SNV positions.
#' rdends_range <- reconstructPPregion(hapMat = ex_hapMatSmall_data, minWindow = 1,
#'                                       posn.lb = 2000, posn.ub = 7000)                      
#' 
reconstructPPregion <- function(hapMat, minWindow, posn.lb = NULL, posn.ub = NULL) {
  
  if(!is.null(posn.lb)){
    posns = subset(hapMat$posns, hapMat$posns>= posn.lb & hapMat$posns <= posn.ub )
    SNVs = which(hapMat$posns %in% posns)
    nSNV = length(SNVs)
  }else{
    nSNV = ncol(hapMat$hapmat)
    SNVs = 1:nSNV
  }
  treeVec = vector(length = nSNV, mode = "list")
  class(treeVec) = "multiPhylo"
  
  message("This may take some time. Please wait...")
  
  for(i in SNVs) {
    treeVec[[i]] = reconstructPP(hapMat, focalSNV = i, minWindow)
    print(paste0("Current focal SNV ", i))
  }
  return(treeVec)
}

