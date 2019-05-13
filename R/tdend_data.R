#' True dendrogram object
#'
#' A phylo object containing attributes of the comparator true dendrogram for the example data at
#' SNV position 975 kilo base pairs.
#'
#' @docType data
#'
#' @usage data(tdend)
#'
#'
#'
#' @format { A phylo object from the ape package containing four attributes:
#' 
#' \describe{
#' 
#' \item{edge}{A matrix containing the node labels and their child nodes.}
#' \item{Nnode}{The number of nodes.}
#' \item{tip.label}{A character vector containing the haplotype labels of the true dendrogram.}
#' \item{edge.length}{ A numeric vector giving the lengths of the branches given by \code{edge}.
#' 
#' }}}
#' 
#'
"tdend"