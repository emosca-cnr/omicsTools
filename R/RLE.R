#' Relative log-expression
#' @description $y_{ij} = x_{ij} - mean(x_i)
#' @param x normalized feature-by-samples matrix in log-space.
#' @param robust if TRUE, use median instead of mean. FALSE by default.
#' @export

RLE <- function(x=NULL, robust=FALSE){

	if(robust){
		ans <- t(apply(x, 1, function(y) y - median(y)))
	}else{
		ans <- t(apply(x, 1, function(y) y - mean(y)))
	}

	return(ans)

}
