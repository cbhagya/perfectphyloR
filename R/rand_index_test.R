#' Rand Index Test
#' 
#'  This function performs Rand index test for association between two \code{phylo} objects.
#'
#' @param dend1 An object of type \code{phylo}.
#' @param dend2 A second object of type \code{phylo}.
#' @param k     The number of clusters that the dendrogram should be cut into.
#' @param nperm The number of times to permute tips of the \code{dend2}.
#'
#' @return A numeric value between 0 and 1 and permutation p-value.
#' @export
#'
#' @references Rand, W.M. (1971). Objective criteria for the evaluation of clustering methods.
#'  Journal of the American Statistical Association 66: 846-850.
#'
#' @examples
#'
#' data(ex_hapMat_data)
#' d1 <- reconstructPP(ex_hapMat_data, focalSNV = 1, minWindow = 1)
#' d2 <- reconstructPP(ex_hapMat_data, focalSNV = 5, minWindow = 1)
#' RandIndexTest(dend1 = d1,  dend2 = d2, k = 5, nperm = 1000)
#'
RandIndexTest <- function(dend1, dend2, k, nperm){
  
  # Ex: Let dend1 = true dend and dend2 = rconstruct dend.


  # P-value by permuting the rows and columns (haplotypes) of the original distance matrix, Dy,
  # instead of permuting the rows of hapMat matrix(rows of original haplotype matrix).
  

  group1 = dendextend::cutree(dend1, k = k)
  group2 = dendextend::cutree(dend2, k = k)

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

RandIndex = function(group1, group2){

  if(length(group1) == length(group2)){
    
    # Compute Rand index
    x <- abs(sapply(group1, function(x) x - group1))
    x[x > 1] <- 1
    y <- abs(sapply(group2, function(x) x - group2))
    y[y > 1] <- 1
    sg <- sum(abs(x - y))/2
    bc <- choose(dim(x)[1], 2)
    RI <- 1 - sg/bc
    
    return(RI)
  }else{ 
    # Adjusted Rand Index for different sizes of groups.
    a <- length(table(group1))
    N <- length(group1)
    ctab <- matrix(NA, a, a)
    for (j in 1:a) {
      for (i in 1:a) {
        ctab[j, i] <- length(which(group2[which(group1 == 
                                                  i)] == j))
      }
    }
    sumnij <- sum(choose(ctab, 2))
    sumai <- sum(choose(colSums(ctab), 2))
    sumbj <- sum(choose(rowSums(ctab), 2))
    Ntwo <- choose(N, 2)
    ari <- abs((sumnij - (sumai * sumbj)/Ntwo)/(0.5 * (sumai + 
                                                         sumbj) - (sumai * sumbj)/Ntwo))
    return(ari)
  }
}


