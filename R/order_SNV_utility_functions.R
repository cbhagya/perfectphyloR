
#' Utility function needed by \code{\link{orderSNVs}}
#'
#' This function orders the compatible SNVs in the window around the focal point by their ancestry.
#'
#' This function sorts columns of M as binary numbers in descending order.
#' The first row of M is the first digit of the binary numbers, the second row of M is the second digit, etc.
#' This guarantees that more ancestral SNVs appear before more recent SNVs.
#'
#' @param M hapMat object with compatible SNVs.
#'
#' @return A hapMat object with SNVs ordered according to their ancestry.
#' @seealso \code{\link{orderSNVs}}
#' @keywords internal
#'
#' @examples 
#' 
#' \dontshow{
#' 
#' data(ex_hapMatSmall_data)
#' 
#' # Order SNV by their ancestry.
#' 
#' ex_hapMatSmall_data$hapmat[, orderColsAncestry(ex_hapMatSmall_data$hapmat)]
#' }
#'
orderColsAncestry = function(M) {
  tem = as.data.frame(t(M))
  ord = do.call(order,tem)
  ord = rev(ord) # we want descending order
  return(ord)
}
