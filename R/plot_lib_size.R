#' Plot library size
#' @param X genes-by-samples count matrix
#' @param col vector of colors; this must match the number of columns of X
#' @param sort whether to sort the samples by library size
#' @return  ggplot2 object
#' @importFrom ggplot2 set_theme ggplot aes geom_col theme element_text labs
#' @export 

plot_lib_size <- function(X=NULL, col=NULL, sort=FALSE){
  
  id <- N <- NULL
  
  if(is.null(col)){
    col <- "black"
  }
  lib_size <- data.frame(id=factor(colnames(X), levels=colnames(X)), N=colSums(X), col=col)
  if(sort){
    lib_size <- lib_size[order(lib_size$N), ]
    lib_size$id <- factor(lib_size$id, levels=lib_size$id)
  }
  
  set_theme(theme_science())
  
  p <- ggplot(lib_size, aes(x=id, y=N)) +
    geom_col(colour = lib_size$col) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), axis.title.y = element_text(angle = 0, vjust = 0.5)) +
    labs(x="")
  
  return(p)
  
}
