#' Compute the user-provided association statistics
#' 
#' This function computes the user-given association statistics and return the association
#' profile over a genomic region.
#' 
#'
#' @param Dx A numeric matrix of pairwise distances.
#' @param Dy A second numeric matrix of pairwise distances.
#' @param hapMat An object of \code{hapMat}.
#' @param nperm The number of times to permute the rows and columns of \code{Dy}.
#' @param method  Association methods. Use "HHG" for HHG test, "dCor" for dCor test, "Mantel" for 
#'               mantel test, "RV" for RV test and "RI" for Rand index test.
#' @param xlab An optional character string for the label on the x-axis(none by default). 
#' @param ylab An optional character string for the label on the y-axis(none by default).
#' @param main An optional character string for title(none by default).
#'
#' @return A list of permutation p-values, corresponding statistics and a plot of association profile 
#'         over SNV location.
#' @keywords internal
#'
#' 
assoTest <- function(Dx, Dy, hapMat, nperm, method, xlab = "", ylab = "", main = ""){
  # Dx: Distance matrix of the base dendrogram
  # Dy: A list of distance matrices of reconstructed dendrograms at each focal SNV
  # nperm: number of permutations
  # method: Association methods(Ex: dCor, HHG, Mantel, and RV)

  
  n = length(Dy)
  assoStats = vector(length = n, mode = "list")
  
  # All permutation statistics by each test across each SNV. Rows represent SNVs
  # and columns represent permutation number.
  
  permStatMat = matrix(NA, nrow = n, ncol = nperm)

  for(i in 1:n){

    assoStats[[i]] = performTest(Dx = Dx, Dy = Dy[[i]], testname = method, nperm = nperm)
    
    permStatMat[i, ] = assoStats[[i]]$permStats
    
  }
  
  
  # True  statistcs for association between each SNV and disease status(observed value, not the 
  # permuted one).
   trueStats = rep(NA, length(assoStats))

  for(j in 1:length(assoStats)){
    
    trueStats[j] = assoStats[[j]]$Stat
  
  }
  
  # P-value at each SNV position
  if( nperm != 0){
    
    mar_pval = rep(NA, length(assoStats))
    
    for(k in 1:length(assoStats)){
      
      mar_pval[k] = assoStats[[k]]$pValue
    
    }
  }else{
    mar_pval = NA
  }
  
  # P-value for over all association. (Omnibus p-value)
  omPvalue = (sum(apply(permStatMat, 2, max) >= max(trueStats))+1)/(nperm + 1)
  
  # Association profile over SNVs
  plt <- plot(hapMat$posns, trueStats, xlab = xlab, ylab = ylab, main = main )
  
  
  return(list(plt, Stats = trueStats, OmPval = omPvalue, mPval = mar_pval))
  
}

#' Perform the user-provided association test
#'
#' This function directs to perform the user-given method.
#' 
#' @param Dx A numeric matrix of pairwise distances.
#' @param Dy A second numeric matrix of pairwise distances.
#' @param testname Association methods. Use "HHG" for HHG test, "dCor" for dCOr test, "Mantel" for 
#'               mantel test, "RV" for RV test and "RI" for Rand index test.
#' @param nperm Numer of permutations.
#'
#' 
#' @keywords internal
#'
#'
performTest <- function(Dx, Dy, testname, nperm){

  if(testname == "HHG"){
    return( HHGtest(Dx = Dx, Dy = Dy, nperm = nperm))
  }
  if(testname == "dCor"){
    return( dCorTest(Dx = Dx, Dy = Dy, nperm = nperm))
  }
  if(testname == "Mantel"){
    return(MantelTest(Dx = Dx, Dy = Dy, nperm = nperm))
  }
  if(testname == "RV"){
    return(RVtest(Dx = Dx, Dy = Dy, nperm = nperm))
  }
}

