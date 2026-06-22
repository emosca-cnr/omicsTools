#' RLE plot
#' @param X Relative log expression
#' @param col vector of colors; this must match the number of columns of X
#' @return  ggplot2 object
#' @importFrom ggplot2 set_theme ggplot aes geom_boxplot theme element_text labs geom_hline stat_boxplot

#' @importFrom reshape2 melt
#' @export 

plot_RLE <- function(X=NULL, col=NULL){
  
  Var2 <- value <- NULL
  
  if(is.null(col)){
    col <- "black"
  }
  
  plot_data <- melt(X)
  
  set_theme(theme_science())
  
  p <- ggplot(plot_data, aes(x=Var2, y=value)) +
    geom_boxplot(outliers = F, fill="gray", lwd=.2) +
    stat_boxplot(geom='errorbar', linetype=1, width=0.5, lwd=.2) +
    geom_hline(yintercept=0, linetype=2, lwd=.2) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), axis.title.y = element_text(angle = 0, vjust = 0.5)) +
    labs(x="", y="RLE")
  
  return(p)
  
}
