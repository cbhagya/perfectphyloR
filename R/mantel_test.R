
#' Mantel test for association of two distance matirces
#' 
#' This function performs Mantel test for correlation between two distance matrices. It 
#' computes P value by randomly permuting rows and columns of the second matrix. 
#'
#' @param Dx  A numeric matrix of pairwise distances.
#' @param Dy  A second numeric matrix of pairwise distances.
#' @param nperm The number of times to permute the rows and columns of \code{Dy}.
#'
#' @return A list contains Mantel statistic and permutation P value.
#' 
#' @references Mantel, N. (1967) The detection of disease clustering and a generalized
#'             regression approach. Cancer Research, 27, 209 - 220.
#' 
#' @export
#'
#' @examples
#' 
#' x <- runif(8)
#' y <- runif(8)
#' # Distance matrices
#' distX = as.matrix(dist(x, upper = TRUE, diag = TRUE))
#' distY = as.matrix(dist(y, upper = TRUE, diag = TRUE))
#'
#' MantelTest(Dx = distX, Dy = distY, nperm = 1000)
#' 
#' 
MantelTest <- function(Dx, Dy, nperm){

  MTObs <- mantelStat(Dx, Dy)
  if(nperm != 0 ){
    permStats <- rep(NA, nperm)
    # P-value by permuting the rows and columns (haplotypes) of the original distance matrix, Dy
    # instead of permuting the rows of hapMat matrix(rows of original haplotype matrix).
    s <-lapply(1:nperm, function(x) c(sample(nrow(Dy))))
  
  
    for(i in 1:nperm){
  
      permStats[i] <- mantelStat(Dx = Dx, Dy = Dy[s[[i]], s[[i]]])
  
    }
    pVal <- (sum(permStats > MTObs) + 1)/(nperm + 1)
    return(list(Stat = MTObs, pValue = pVal, permStats = permStats))
  }
  else{
    return(list(Stat = MTObs))
  }
}

#' This function performs Mantel test for correlation between two matrices.
#'
#' @param Dx  A numeric matrix of pairwise distances.
#' @param Dy  A second numeric matrix of pairwise distances.
#'
#' @return A matrix of numeric values.
#' 
#' @keywords internal
#' 
#'
mantelStat = function(Dx, Dy){

  v1 = as.numeric(Dx[lower.tri(Dx)])
  v2 = as.numeric(Dy[lower.tri(Dy)])
  return(cor(v1,v2))

}
