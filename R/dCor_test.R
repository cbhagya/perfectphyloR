#' dCor test for similarity of two matrices
#' 
#' This function performs dCor test for association between two distance matrices and computes permutation 
#' P value. Permutation P value is computed by randomly permuting rows and columns of the 
#' second distance matrix.
#'
#' @param Dx A numeric matrix of pairwise distances.
#' @param Dy A second numeric matrix of pairwise distances.
#' @param nperm The number of times to permute the rows and columns of \code{Dy}.
#'
#' @return A list contains RV coefficient and permutation P value.
#' 
#' @references  G. J. Szekely, M. L. Rizzo, and N. K. Bakirov. (2007). Measuring and testing dependence 
#'              by correlation of distances. The Annals of Statistics, 35(6):2769 - 2794.
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
#' dCorTest(Dx = distX, Dy = distY, nperm = 1000)
#' 

dCorTest <- function(Dx, Dy, nperm){


  n <- dim(Dx)[1]
  C <- diag(n) - ((rep(1, n) %*%t (rep(1, n)))/n)

  # Function 'dCor()' in 'dCorCoeff.cpp'
  dCorObs <- as.vector(dCor(mDx = Dx, mDy = Dy, mC = C))
  
  if(nperm != 0 ){
    
    permStats <- rep(NA, nperm)
    # P-value by permuting the rows and columns (haplotypes) of the original distance matrix, Dy
    # instead of permuting the rows of hapMat matrix(rows of original haplotype matrix).
    s <-lapply(1:nperm, function(x) c(sample(nrow(Dy)))) 
  
  
    for(i in 1:nperm){
  
      permStats[i] <- as.vector(dCor(mDx = Dx, mDy = Dy[s[[i]], s[[i]]], mC = C))
  
    }
    pVal <- (sum(permStats > dCorObs) + 1)/(nperm + 1)
  
    return(list(Stat = dCorObs, pValue = pVal, permStats = permStats)) 
  }
  else{
    return(list(Stat = dCorObs))
  }
}
