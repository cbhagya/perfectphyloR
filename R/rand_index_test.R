#' Rand Index Test
#' 
#'  This function performs Rand index test for association between two \code{phylo} objects.
#'
#' @param dend1 An object of type \code{phylo}.
#' @param dend2 A second object of type \code{phylo}.
#' @param k     An integer that specifies the number of clusters that the dendrogram should be cut into. 
#'              The default is \code{k = 2}. Clusters are defined by starting from the root of the dendrogram and 
#'              cutting across.
#' @param nperm The number of times to permute tips of the \code{dend2}.
#'
#' @return A numeric value between 0 and 1 and permutation P value.
#' @export
#'
#' @references Rand, W.M. (1971) Objective criteria for the evaluation of clustering methods.
#'             Journal of the American Statistical Association 66: 846-850.
#'
#' @examples
#'
#' data(ex_hapMat_data)
#' d1 <- reconstructPP(ex_hapMat_data, focalSNV = 1, minWindow = 1)
#' d2 <- reconstructPP(ex_hapMat_data, focalSNV = 5, minWindow = 1)
#' RandIndexTest(dend1 = d1,  dend2 = d2, k = 5, nperm = 100)
#'
RandIndexTest <- function(dend1, dend2, k = 2, nperm){
  
  # Ex: Let dend1 = true dend and dend2 = rconstruct dend.


  # P-value by permuting the rows and columns (haplotypes) of the original distance matrix, Dy,
  # instead of permuting the rows of hapMat matrix(rows of original haplotype matrix).
  

  g1 = dendextend::cutree(dend1, k = k)
  g2 = dendextend::cutree(dend2, k = k)
  
  # splitTips(), utility function to separate haplotypes names for haplotypes, 
  # if they can not be distinguished in the window around the focal point.   
  group1 = splitTips(group = g1)
  group2 = splitTips(group = g2)
  
  
  RIObs <- RandIndex(group1, group2)

  if(nperm > 0){
    s <-lapply(1:nperm, function(x) c(sample(group2)))

    permStats <- rep(NA, nperm)
    for(i in 1:nperm){

      permStats[i] <- RandIndex(group1, group2 = s[[i]])

  }
    pVal <- (sum(permStats > RIObs) + 1)/(nperm + 1)

    return(list(Stat = RIObs, pValue = pVal, permStats = permStats))
  }

  else{return(list(Stat = RIObs, pValue = NA))}
}

#' Rand Index
#'
#' This function calculates the Rand index, a  measure of similarity between two
#' partitions and returns a value between 0 and 1 reflecting the
#' agreement of the partitions.
#'
#' @param group1 First group.
#' @param group2 Second group.
#'
#'
#' @return A numeric value between 0 and 1.
#' @keywords internal
#'
#' @examples 
#' 
#' \dontshow{
#' 
#' g1 <- sample(1:2, size=10, replace=TRUE)
#' g2 <- sample(1:3, size=10, replace=TRUE)
#' 
#' RandIndex(g1,g2)
#'
#' 
#' }
#'

RandIndex = function(group1, group2){

    
    # Compute Rand index
    x <- abs(sapply(group1, function(x) x - group1))
    x[x > 1] <- 1
    y <- abs(sapply(group2, function(x) x - group2))
    y[y > 1] <- 1
    sg <- sum(abs(x - y))/2
    bc <- choose(dim(x)[1], 2)
    RI <- 1 - sg/bc
    
    return(RI)
}

#' Separate haplotype names for haplotypes that can not be distingushed in the window around the focal point
#' 
#' This function separates haplotypes names for haplotypes, if they can not be distinguished in the window
#' around the focal point. 
#' 
#'
#' @param group A numeric vector that represent the group indices returned by \code{cutree()}.
#'
#' @return A vector of haplotypes.
#' 
#' @keywords internal
#'
#' 
#' 
splitTips <- function(group){
  
  spltTips = strsplit(attr(group,"names"), split="-") 
  
  allHaps = unlist(spltTips) 
  
  
  spltGroup = rep(as.numeric(group[1]), length(spltTips[[1]]))
  
  for(i in 2:length(spltTips)){
    
    x2 = rep(as.numeric(group[i]), length(spltTips[[i]]))
    spltGroup = append(spltGroup, x2)
    
  }
  
  class(spltGroup)
  attr(spltGroup, "names") = allHaps
  
  return(spltGroup)
}
