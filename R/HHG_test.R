#' HHG test for association of two distance matrices
#' 
#' This function performs HHG test to find the association between two distance matrices. It permutes rows and columns
#' of the second matrix randomly to calculate P value.
#'
#' @param Dx  A numeric matrix of pairwise distances.
#' @param Dy  A second numeric matrix of pairwise distances.
#' @param nperm The number of times to permute  the rows and columns of \code{Dy}.
#'
#' @return A list contains HHG coefficient and permutation P value.
#' 
#' @references Barak, B., and Shachar, K., based in part on an earlier implementation by Ruth Heller 
#' and Yair Heller. (2017). HHG: Heller-Heller-Gorfine Tests of Independence and Equality of Distributions. R
#' package version 2.2. https://CRAN.R-project.org/package=HHG
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
#' HHGtest(Dx = distX, Dy = distY, nperm = 1000) 
#'
HHGtest <- function(Dx, Dy, nperm){
 
  if (!requireNamespace("HHG", quietly = TRUE)) {
    stop("Package \"HHG\" needed for this function to work. Please install it.",
         call. = FALSE)
  }

  hhgObs <- HHG::hhg.test(Dx = Dx, Dy = Dy, nr.perm = 0)$sum.chisq
  
  if(nperm != 0 ){
    permStats <- rep(NA, nperm)
    # P-value by permuting the rows and columns (haplotypes) of the original distance matrix, Dy
    # instead of permuting the rows of hapMat matrix(rows of original haplotype matrix).
    s <-lapply(1:nperm, function(x) c(sample(nrow(Dy))))
    
  
    for(i in 1:nperm){
  
      # By HHG R package
      permStats[i] <- HHG::hhg.test(Dx = Dx, Dy = Dy[s[[i]], s[[i]]], nr.perm = 0)$sum.chisq
  
    }
    pVal <- (sum(permStats > hhgObs) + 1)/(nperm + 1)
    
    return(list(Stat = hhgObs, pValue = pVal, permStats = permStats))
  }
  else{
    return(list(Stat = hhgObs))
  }
}
