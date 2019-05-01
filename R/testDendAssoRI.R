
#' Tests Rand Index between a comparator dendrogram and reconstructed dendrograms
#' 
#' 
#' This function performs the Rand Index between a user-supplied comparator dendrogram and the reconstruced dendrograms
#' at each focal SNV position in a genomic region.
#'
#' @param rdend  A multiphylo object of reconstructed dendrograms at each focal SNV.
#' @param cdend  A phylo object of the comparator dendrogram.
#' @param hapMat An object of class `hapMat`containing SNV haplotypes.
#' @param k      An integer that specifies the number of clusters that the dendrogram should be cut into.
#'               The default is k=2. Clusters are defined by starting from the root of the
#'               dendrogram and moving towards the tips, cutting horizontally at any given point
#'               in the dendrogram.
#' @param nperm  Number of permutations for the test of any association across the genomic region of interest. 
#'               The default is 'nperm = 0';i.e., association will not be tested.
#' @param xlab   An optional character string for the label on the x-axis in the plot that is returned (none by
#'                default). 
#' @param ylab  An optional character string for the label on the y-axis in the plot that is returned (none by
#'                default).    
#' @param main An optional character string for title in the plot that is returned (none by default). 
#'
#' @return  A list with the following components:
#' 
#'          `Stats`  A vector of observed Rand indices.
#'          `OmPval` A permutation-based omnibus P value for the test of any association across the genomic 
#'                   region using the maximum Rand index over the genomic region as the test statistics.
#'          `mPval`  A vector of marginal P values at each SNV position.
#'          `plt`    A plot of the association profile of Rand indices over SNV locations in the
#'                   region of interest.
#'          
#' @export
#'
#' @examples  See the section 'Applications' in vignette("perfectphyloR") for the detailed example.
#' 
testDendAssoRI <- function(rdend, cdend, hapMat, k = 2, nperm = 0, xlab = "", ylab = "", main = ""){
  
   if(is.null(k) || k <= 0){
     stop("k should be > 0")
   }else{

     RI_vals = vector(length = length(rdend), mode = "list")
     
     
     # All permutation statistics by each test across each SNV. Rows represent SNVs
     # and columns represent permutation number.
     
     permStatMat = matrix(NA, nrow = length(rdend), ncol = nperm)
     
     
     for(i in 1:length(rdend)){
       
       # If k > max value that the clusters can be cut into, skip those reconstructed dendrograms. 
       if( length(rdend[[i]]$tip.label) < k){ next
         
       }else{
         RI_vals[[i]] <- RandIndexTest(dend1 = cdend, dend2 = rdend[[i]], k = k
                                    , nperm = nperm)
       
        permStatMat[i, ] = RI_vals[[i]]$permStats
       }
     }
     
     # True  statistcs for association between each SNV and disease status(observed value, not the 
     # permuted one).
     trueStats = rep(NA, length(RI_vals))
     
     for(j in 1:length(RI_vals)){
      
       if(is.null(RI_vals[[j]]$Stat)) next
       
       trueStats[j] = RI_vals[[j]]$Stat
     }
     
     # P-value at each SNV position
     
     mar_pval = rep(NA, length(RI_vals))
     
     for(k in 1:length(RI_vals)){
    
       if(is.null(RI_vals[[k]]$pValue)) next
       
       mar_pval[k] = RI_vals[[k]]$pValue
     }
     
     # P-value for over all association. (Omnibus p-value)
     omPvalue = (sum(apply(permStatMat, 2, max) >= max(trueStats))+1)/(nperm + 1)
     
     # Association profile over SNVs
     plt <- plot(hapMat$posns, trueStats, xlab = xlab, ylab = ylab, main = main )
     
    
    # Return association profile plots,
    # omnibus p value, marginal p values and true Rand indices.
     return(list(plt, Stats = trueStats, OmPval = omPvalue, mPval = mar_pval))
   }
}
