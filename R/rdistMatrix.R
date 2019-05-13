#' Rank-based distances between haplotypes in a given partition
#'
#' This function computes the pairwise distances between haplotypes (tips) of the dendrogram based 
#' on the ranking of the nested partitions in the dendrogram. See the details.
#' 
#' 
#' 
#'
#' @param dend A list of nodes that represents the nested partition of haplotypes.
#' @param sep  A character string separator for concatenating haplotype labels in the
#'             dendrogram if they are undistingushable in the window around the focal SNV.
#'             See the arguments in \code{\link{reconstructPP}}.
#'
#' @return A matrix of pairwise distances between haplotypes.
#' @export
#'
#' @details We code the distance between two haplotypes of a dendrogram as the number of inner nodes that
#'          seperate the haplotypes plus one. That is, we assign the distance between two internal neighbouring 
#'          nodes as one, and the distance between an internal node and its neighbouring tip as one. To illustrate, 
#'          consider the following figure of a dendrogram. In the figure, the distance between the haplotypes 2931 and 454 is 3;
#'          the distance between other haplotypes are given in the table below.
#'
#' \if{html}{\figure{treedist.png}{options: width = "30\%" alt = "Figure: treedist.png"}}
#' \if{latex}{\figure{treedist.pdf}{options: width = 10cm}}
#'
#' \if{html}{\figure{trtable.png}{options: width = "30\%" alt = "Figure: trtable.png"}}
#' \if{latex}{\figure{trtable.pdf}{options: width = 9cm}}
#'
#'
#' @examples
#' data(ex_hapMat_data)
#' rdend <- reconstructPP(hapMat = ex_hapMat_data, focalSNV = 2, minWindow = 1, sep = "-" )
#' rdistMatrix(rdend)
#'
#'


rdistMatrix = function(dend, sep = "-"){
  
  #### 
  dend$edge.length[1:dim(dend$edge)[1]] = 1
  # Calculate rank-baesd distances between tips of the tree
  distmat = ape::cophenetic.phylo(dend)
  # A tip may be consist of multiple haplotypes. Set
  # pairwise distances between haplotypes in the same tip
  # to be zero and the distance between haplotypes in different
  # tips to be the distances between the tips.

  splitTips = strsplit(colnames(distmat), split = sep)
  # splitTips is now a list whose ith element is
  # a vector of haplotypes in tip i.
  allHaps = unlist(splitTips)
  nHaps = length(allHaps)
  
  lenTips = c()
  for(i in 1:length(splitTips)){
    lenTips[i] = length(splitTips[[i]])
  }
  
  # distMatRcpp() is inside association package. Load that package to run this function.
  newDistMat <- distMatRcpp(len = length(splitTips),
                                         len2 = lenTips, nHaplos = nHaps, mD = distmat)

  ord <- order(allHaps)
  newDistMat <- newDistMat[ord, ord]
  colnames(newDistMat) = rownames(newDistMat) = allHaps[ord]

  return(newDistMat)
}
