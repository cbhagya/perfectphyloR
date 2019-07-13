#' Recursively partition haplotype matrix
#'
#' This function recursively partitions the SNVs in the window around the focal SNV.
#'
#'
#' This function makes two clades based on \code{splotSNV}. For each partition, update \code{splitSNV} and use
#' it to define subclades. Then, continue recursive partitioning until each partition has only one haplotype,
#' or there are no more SNVs to consider.
#'
#'
#' @param hapmat A hapMat object with SNVs ordered according to ancestry.
#' @param splitSNV The index of the SNV where the haplotype matrix from \code{\link{orderSNVs}} is partitioned.
#'
#' @keywords internal
#'
#' @return A nested partition of haplotypes, implemented as a list of nodes, each with two child nodes.
#'
#' @seealso \code{\link{makeDend}}, \code{\link{newNode}}, \code{\link{noVariation}}
#' 
#' @examples 
#' 
#' \dontshow{
#' 
#' data(ex_hapMatSmall_data)
#' 
#' # First, select a window of SNVs about a focal SNV.
#' SNV_win <- selectWindow(hapMat = ex_hapMatSmall_data,
#'                         focalSNV = 10, minWindow = 1)
#'                         
#' # Then order SNVs in the window.
#' ordHapmat <- orderSNVs(snvWin = SNV_win)
#' 
#' # Recursively partition haplotype matrix.
#' partitions <- makePartition(hapmat = ordHapmat, splitSNV = 1)
#' 
#' }
#'
makePartition = function(hapmat, splitSNV) {
  if(nrow(hapmat)==1 || splitSNV>ncol(hapmat)){
    # Then we are done splitting, either because the clade is
    # a single haplotype (nrow(hapmat==1) or because we've
    # run out of SNVs (splitSNV>ncol(hapmat)). Return the
    # haplotypes in hapmat as a leaf node.
    return(newNode(hapmat))
  }
  # If not, find the next SNV to split on. To split on a SNV, there
  # must be variation. Keep searching for a SNV as long as there is
  # no variation in the current SNV.
  while(splitSNV <= ncol(hapmat) && noVariation(hapmat[, splitSNV])) {
    splitSNV = splitSNV + 1
  }
  # We may have found a SNV to split on, or we may have hit
  # the end of hapMat without finding a SNV to split on.
  if(splitSNV > ncol(hapmat)){
    # Couldn't find a SNV to split on; return hapmat as a leaf node.
    return(newNode(hapmat))
  }
  # Otherwise, we've found a SNV to split on, so split into clades.
  # The following call is to R's subset(), applied to a matrix,
  # **not** subsetHapMat() applied to a hapMat object.
  clade1 = subset(hapmat, hapmat[, splitSNV] == 1)
  child1 = makePartition(clade1, splitSNV+1)
  clade0 = subset(hapmat,hapmat[, splitSNV] == 0)
  child0 = makePartition(clade0, splitSNV+1)
  return(newNode(hapmat, child1, child0, depth = nrow(hapmat)-1))
}
#---------------------------------------------------------------#
#' Create a list of child nodes
#'
#' This function creates a pair of child nodes for a parent node.
#'
#' @param hapmat The hapMat object with columns ordered by ancestry.
#' @param child1 The child node from splitting on the mutant allele at the next SNV in the ordered
#'               neighborhood.
#' @param child0 The child node from splitting on the non-mutant allele at the next SNV in the ordered
#'               neighborhood.
#'
#' @keywords internal
#' @seealso \code{\link{makePartition}}
#' 
#' @examples 
#' 
#' \dontshow{
#'  
#' data(ex_hapMatSmall_data)  
#' 
#' # First, select a window of SNVs about a focal SNV.
#' SNV_win <- selectWindow(hapMat = ex_hapMatSmall_data,
#'                         focalSNV = 10, minWindow = 1)
#'                         
#' # Then order SNVs in the window.
#' ordHapmat <- orderSNVs(snvWin = SNV_win)
#' 
#' # Create a list of child nodes.
#' chldNodes <- newNode(hapmat = ordHapmat)
#' 
#' }
#'
newNode = function(hapmat, child1 = NULL, child0 = NULL, depth = 0) {
  return(list(haps = rownames(hapmat), child1 = child1, child0 = child0, depth = depth))
}
#---------------------------------------------------------------#
#' Check the variation in a SNV
#'
#' This function checks the variation in a specified SNV and is applied during the recursive partitioning.
#'
#' @param snv A SNV to check for variation among haplotypes.
#'
#' @return Logical:TRUE, if there is no variation in the SNV; FALSE otherwise.
#' @keywords internal
#' @seealso \code{\link{makePartition}}
#' 
#' @examples 
#' 
#' \dontshow{
#' 
#' data(ex_hapMatSmall_data)  
#' 
#' # Check the variation in a SNV.
#' noVariation(ex_hapMatSmall_data$hapmat[,1])
#' 
#' }
#'
noVariation = function(snv) {
  if(sd(snv) == 0) return(TRUE) else return(FALSE)
}

#' Convert a list data structure to Newick format
#'
#' This function traverses the dendogram to build the character string that represents the dendrogram in
#' Newick format.
#'
#'
#' @param dend A list of nodes that represents the nested partition of haplotypes.
#' @param sep A separator for haplotype names in the same clade. See the arguments in
#'           \code{\link{reconstructPP}}.
#'
#' @return A character string in Newick format.
#' @keywords internal
#' @seealso \code{\link{makeDend}}
#' 
#' @examples 
#' 
#' \dontshow{
#' 
#' data(ex_hapMatSmall_data)
#' 
#' # First, select a window of SNVs about a focal SNV.
#' SNV_win <- selectWindow(hapMat = ex_hapMatSmall_data,
#'                         focalSNV = 10, minWindow = 1)
#'                         
#' # Then order SNVs in the window.
#' ordHapmat <- orderSNVs(snvWin = SNV_win)
#' 
#' # Recursively partition haplotype matrix.
#' partitions <- makePartition(hapmat = ordHapmat, splitSNV = 1)
#' 
#' # Dendrogram in Newick format.
#' newickDend <- dendToNewick(dend = partitions, sep = "-")
#' 
#' }
#'
dendToNewick = function(dend, sep = "-"){
  # Arguments:
  #   dend is the output of makeTreeRec
  #   sep is a character string to separate haplotype names for
  #     tips comprised of multiple haplotypes (e.g, if a
  #     tip contained haplotypes C and D, the tip would
  #     appear as C-D in the Newick string).
  dendStr = makeNewickRec(dend, sep)
  # Now just append a ";" to mark the end of the dend
  return(paste(dendStr, ";", sep=""))}

#' Build the character string of nodes and haplotypes in Newick format
#'
#'
#' @param node Tree(dend) node
#' @param sep  A separator for haplotype names in the same clade. See the arguments in
#'           \code{\link{reconstructPP}}.
#'
#' @keywords internal
#' 
#' @examples 
#' 
#' \dontshow{
#' 
#' data(ex_hapMatSmall_data)
#' 
#' # First, select a window of SNVs about a focal SNV.
#' SNV_win <- selectWindow(hapMat = ex_hapMatSmall_data,
#'                         focalSNV = 10, minWindow = 1)
#'                         
#' # Then order SNVs in the window.
#' ordHapmat <- orderSNVs(snvWin = SNV_win)
#' 
#' # Recursively partition haplotype matrix.
#' partitions <- makePartition(hapmat = ordHapmat, splitSNV = 1)
#' 
#' dendStrng <- makeNewickRec(node = partitions, sep = "-")
#' 
#' }
#' 
makeNewickRec = function(node,sep) {
  # leaf nodes have two NULL children, internal nodes have
  # two non-NULL children.
  if(!is.null(node$child1)) {
    #internal -- get strings from children and parse as
    #   ( child1:len1 , child0:len2 )
    # where len1 is length of branch between node and child1
    # and len2 is length of branch between node and child2.
    len1 <- node$depth - node$child1$depth
    child1Str = makeNewickRec(node$child1,sep)
    len2 <- node$depth - node$child0$depth
    child0Str = makeNewickRec(node$child0,sep)
    return(paste("(",child1Str,":",len1,",",child0Str,":",len2,")",sep=""))
  }else{ 
    # leaf -- just return my haplotype label
    # If there are multiple labels in this leaf, separate with "sep"

    return(paste(node$haps,collapse=sep))
  }
}
