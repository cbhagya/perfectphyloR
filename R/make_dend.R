#' Recursive partitioning of sequences at a focal point.
#'
#' This function recursively partitions a hapMat object with SNVs ordered by ancestry, and records the
#' partitioning in a \code{phylo} object.
#'
#'
#' @param hapmat A hapMat object with ordered columns according to the ancestry of SNVs.
#' @param sep A separator for haplotype names in the same clade. See the arguments in
#'            \code{\link{reconstructPP}}
#'
#' @return An object of class \code{phylo}.
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
#' # Recursive partitioning of sequences. 
#' dend <- makeDend(ordHapmat, sep = "-")
#'
#' }
#'
makeDend = function(hapmat, sep = "-") {
# hapmat is the haplotype matrix, with columns ordered
#    in the order of partions we want
# sep is a separator for haplotype names in the same clade.
#   Note of caution: In Newick format, branches at an internal
#   node are separated by a comma, and multifurcating nodes
#   are allowed. Therefore using sep="," as the separator
#   will cause ape's plot function to represent the clade as
#   a multifurcating node with as many leaves as there are
#   haplotypes in the clade. This is probably not what a user wants.
ltree = makePartition(hapmat, splitSNV = 1) # our list-based data structure
ntree = dendToNewick(ltree, sep = sep) # Newick format -- see make_dend_utility_functions.R


ptree <- phytools::read.newick(text = ntree)
return(ptree)
}
