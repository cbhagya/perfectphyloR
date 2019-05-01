
#' RV test for association of two distance matrices
#' 
#' This function performs RV test for similarity of two matrices. It permutes rows and columns
#' of the second matrix randomly to calculate p-value.
#'
#' @param Dx  A numeric matrix of pairwise distances.
#' @param Dy  A second numeric matrix of pairwise distances.
#' @param nperm The number of times to permute  the rows and columns of \code{Dy}.
#'
#' @return A list contains RV coefficient and permutation p-value.
#' 
#' @references Robert, P. and Escoufier, Y. (1976).A Unifying tool for linear multivariate statistical
#'             methods: the RV-coefficient. Applied Statistics, Vol.25, No.3, p. 257-265.
#' @export
#'
#' @examples
#' x <- runif(8)
#' y <- runif(8)
#' # Distance matrices
#' distX = as.matrix(dist(x, upper = TRUE, diag = TRUE))
#' distY = as.matrix(dist(y, upper = TRUE, diag = TRUE))
#'
#' RVtest(Dx = distX, Dy = distY, nperm = 1000)
#' 
RVtest <- function(Dx, Dy, nperm){

  n <- dim(Dx)[1]
  C <- diag(n) - ((rep(1, n) %*% t(rep(1, n)))/n)

  # Compute RV coefficient using RcppArmadillo and return a vector using as.vector() because RVcoeff()
  # returns a matrix(1 by 1). Note that, RcppArmadillo functions return matrices only.
  RVObs <- as.vector(RVcoeff(mDx = Dx, mDy = Dy, mC = C))
  
  if(nperm != 0 ){
    permStats <- rep(NA, nperm)
    # P-value by permuting the rows and columns (haplotypes) of the original distance matrix, Dy
    # instead of permuting the rows of hapMat matrix(rows of original haplotype matrix).
    s <-lapply(1:nperm, function(x) c(sample(nrow(Dy))))
  
  
    for(i in 1:nperm){
  
      permStats[i] <- as.vector(RVcoeff(mDx = Dx, mDy = Dy[s[[i]], s[[i]]],  mC = C))
  
    }
  
    pVal <- (sum(permStats > RVObs) + 1)/(nperm + 1)
    
   
    return(list(Stat = RVObs, pValue = pVal, permStats = permStats))
    
  }
  else{
    return(list(Stat = RVObs))
  }
}

