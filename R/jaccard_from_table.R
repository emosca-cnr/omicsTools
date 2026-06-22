#' Jaccard index from a table of intersections
#' @param x table of intersections between values of two categorical variables (intersections)
#' @export

jaccard_from_table <- function(x=NULL){
  
  
  rs <- matrix(rep(rowSums(x), ncol(x)), ncol=ncol(x))
  cs <- matrix(rep(colSums(x), nrow(x)), nrow=nrow(x), byrow = T)
  
  j <- x / (rs+cs-x)
  
  return(j)
}