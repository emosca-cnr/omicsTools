#' RLE
#'
#'
#' @export

RLE <- function(x, robust=TRUE, use.logs=FALSE){

	if(robust){
		ans <- t(apply(x, 1, function(y) y - median(y)))
	}else{
		ans <- t(apply(x, 1, function(y) y - mean(y)))
	}

	return(ans)

}
