#'  Plot reconstructed dendrogram
#'
#'  This function plots reconstructed dendrogram in a genomic region.
#'
#' @param dend  An object of class \code{phylo} or of class \code{multiPhylo} returned from
#'              \code{reconstructPP} or \code{reconsPPregion}.
#' @param direction A character string specifying the direction of the dendrogram. Four values are
#'                  possible: "downwards" (the default),  "upwards", "leftwards" and "rightwards".
#'
#'
#' 
#' @export
#'
#' @examples
#'
#' data(ex_hapMat_data)
#'
#' ex_dend <- reconstructPP(hapMat = ex_hapMat_data,
#'                          focalSNV = 3,
#'                          minWindow = 1,
#'                          sep = "-")
#'
#' plotDend(dend = ex_dend, direction = "downwards")
#'
plotDend <- function(dend, direction = "downwards"){

 
  
  # Assign distance between two internal neighbouring
  # nodes as one, and the distance between an internal node and its neighbouring tip as one.
  dend$edge.length[1:dim(dend$edge)[1]] = 1
  return(plot(dend, direction = direction))
}


