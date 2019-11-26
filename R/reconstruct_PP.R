#' Reconstruct the perfect phylogeny at a given focal SNV
#'
#' This function reconstructs the perfect phylogeny at a given focal SNV using the recursive partitioning
#' algorithm of Gusfield (1991) on compatible SNVs, and the modification of Mailund et al. (2006)
#' to include incompatible SNVs that are nearby.
#'
#' To reconstruct the perfect phylogeny from sequence data, these two steps are followed:
#' (1) Select a window of SNVs at a given focal SNV.
#' (2) Build the perfect phylogey for the window of SNVs. More details can be found in the references.
#' 
#'
#' The following figure shows the reconstructed partitions at the tenth SNV position of \code{ex_hapMatSmall_data}.
#'
#' \if{html}{\figure{tree.png}{options: width = "76\%" alt = "Figure: tree.png"}}
#' \if{latex}{\figure{tree.pdf}{options: width = 10cm}}
#'
#' @param hapMat    A data structure of class \code{hapMat}. Eg: created by the \code{\link{createHapMat}}
#'                  function.
#' @param focalSNV  The column number of the focal SNV at which to reconstruct the reconstructed partitions.
#' @param minWindow Minimum number of SNVs around the focal SNV in the window of SNVs used to reconstruct the
#'                  partitions.
#' @param sep       Character string separator to separate haplotype names for haplotypes
#'                  that can not be distingushed in the window around the focal point. For example, if a tip is comprised
#'                  of haplotypes "h1" and "h3", and sep = "-", then the tip label will be "h1-h3". The default 
#'                  value is \code{"-"}. See details.
#'
#' @return An object of class \code{phylo} with indices of the column boundaries of the \code{hapMat} object 
#'         that were used to reconstruct the partition in the window of SNVs.
#'         
#'
#' @references  Gusfield, D. (1991) Efficient algorithms for inferring evolutionary trees.
#'              Networks, 21(1), 19-28.
#' @references  Mailund, T., Besenbacher, S., and Schierup, M. H. (2006) Whole genome association
#'              mapping by incompatibilities and local perfect phylogenies. BMC Bioinformatics, 7(1),
#'              454.
#'
#' @export
#'
#' @examples
#'
#' data(ex_hapMatSmall_data)
#'
#' rdend <- reconstructPP(hapMat = ex_hapMatSmall_data,
#'                       focalSNV = 10,
#'                       minWindow = 1,
#'                       sep = "-")
#'
#' # Plot the reconstructed perfect phylogeney.
#' 
#' plotDend(rdend, direction = "down")
#' 
#' # Extract the positions of the lower and upper limits of a window of SNVs in hapMat object 
#' # to reconstruct the partition, rdend.
#' 
#' ex_hapMatSmall_data$posns[rdend$snvWinIndices]
#'   

reconstructPP = function(hapMat, focalSNV, minWindow = 1, sep = "-") {

    
  nminWin = round(length(hapMat$posns)*.02, 0)  # 2% of nSNVs for default winsize.
  
  
  # Step 1: Select a window of SNVs about a focal SNV.
  if(minWindow == 1 & nminWin > 1){
    
    snvWin <- selectWindow(hapMat, focalSNV, minWindow = nminWin)
    
  }else if(minWindow == 1 & nminWin == 1){
    
    snvWin <- selectWindow(hapMat, focalSNV, minWindow = 1)
    
  }else if(minWindow != 1 & nminWin > 1){   
    snvWin <- selectWindow(hapMat, focalSNV, minWindow = minWindow)
  }else{
    snvWin <- selectWindow(hapMat, focalSNV, minWindow = 1)
  }
  
  #  The indices of the columns of the hapMat object that were used in the window.
  snvWinIndices <- which(colnames(hapMat$hapmat) %in% colnames(snvWin$hapMat$hapmat)) 
  #-----------------------------------------------------------------#
  # Step 2: Build the tree for the window of SNVs from step 1.
  dend <- buildDend(snvWin, sep = sep)
  # Left bound and the right bound of the window of SNVs.
  dend$snvWinIndices = snvWinIndices[c(1, length(snvWinIndices))]
  
  return(dend)
  
    
}
