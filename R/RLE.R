#' Relative log-expression
#' @description $y_{ij} = x_{ij} - mean(x_i)
#' @param x normalized feature-by-samples matrix in log-space.
#' @param robust if TRUE, use median instead of mean.
#' @param na.rm whether to remove NA or not.
#' @export
#' @importFrom stats median

RLE <- function(x=NULL, robust=TRUE, na.rm=TRUE){

	if(robust){
		ans <- t(apply(x, 1, function(y) y - median(y, na.rm = na.rm)))
	}else{
		ans <- t(apply(x, 1, function(y) y - mean(y, na.rm = na.rm)))
	}

	return(ans)

}
