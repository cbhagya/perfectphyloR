#' Test the association between a comparator distance matrix, and the reconstructed dendrograms across a genomic region
#'
#' This function calculates and tests the association between a comparator distance matrix, based on
#' any pairwise distance measure, and the reconstructed dendrograms across a genomic region of
#' interest usingassociation measures such as the dCor statistic, HHG statistic, Mantel statistic, and
#' RV coefficient. See the section Applications in \code{vignette("perfectphyloR")} for the detailed example.
#' 
#' @param rdend A \code{multiPhylo} object of reconstructed dendrograms at each focal SNV.
#' @param cdmat A comparator matrix of pairwise distances (e.g. pairwise distances between haplotypes of a 
#'                comparator dendrogram).
#' @param method Association measures. Use "dCor" for dCor test, "HHG" for HHG test, "Mantel" 
#'               for mantel test, and "RV" for RV test.
#' @param hapMat An object of class \code{hapMat} containing SNV haplotypes.
#' @param nperm Number of permutations for the test of any association across the genomic region of interest.
#'              The default is \code{nperm = 0}; i.e., association will not be tested.
#' 
#' @param xlab An optional character string for the label on the x-axis  in the plot that is returned
#'             (none by default). 
#' @param ylab An optional character string for the label on the y-axis  in the plot that is returned
#'             (none by default).
#' @param main An optional character string for title  in the plot that is returned (none by default).
#'
#'
#' @return A list with the following components: 
#' @return \item{Stats}{A vector of observed statistics computed from the user-provided distance association method.}
#' @return \item{OmPval}{A permutation-based omnibus P value for the test of any association across the genomic region 
#'                    using the maximum statistic over the genomic region as the test statistic.}
#' @return \item{mPval}{A vector of marginal P values at each SNV position.}         
#' @return \item{plt}{A plot of the association profile over SNV locations in the region of interest.}
#'  
#'  
#'          
#' @export
#' @seealso \code{\link{HHGtest}}, \code{\link{dCorTest}}, \code{\link{RVtest}}, \code{\link{MantelTest}}
#'
#'               
testAssoDist <- function(rdend, cdmat, method, hapMat, nperm = 0, xlab = "", ylab = "", main = ""){
  
  if( is.null(cdmat) | class(cdmat) != "matrix" ){
      
      stop("cdmat should be a distance matrix.")
    
    
  }else{

    
        Dyy <- vector(length = length(rdend), mode = "list")
        # Compute distance matrices for all reconstructed perfect phylogenies across the region.
        for(i in 1:length(rdend)){
    
          Dyy[[i]] <-  rdistMatrix(rdend[[i]])
        
        }
    
        # Return a plot of association profile over SNVs, omnibus p-value for maximum association on profile
        # and marginal p-values at each SNV across the region.
        a <- assoTest(Dx = cdmat, Dy = Dyy, hapMat = hapMat, nperm = nperm, method = method, 
                      xlab = xlab, ylab = ylab, main = main )
        return(a)
      
  } 
}




