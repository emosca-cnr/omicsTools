#' differential_expression
#'
#' @import limma
#' @export

differential_expression <- function(dge=NULL, design=NULL, contr_mat=NULL){

	if(is.null(contr_mat)){
		cat("all-pairs contrasts...\n")
		contr_mat <- apply(t(combn(colnames(design), 2)), 1, paste, collapse="-")
	}
	cont_mat <- limma::makeContrasts(contrasts = contr_mat, levels = design)
	print(cont_mat)

	v <- limma::voom(dge, design, plot=TRUE)
	fit_l <- limma::lmFit(v, design)
	fit_l <- limma::contrasts.fit(fit_l, cont_mat)
	fit_l <- limma::eBayes(fit_l)

	tt <- lapply(1:ncol(fit_l$coefficients), function(x) limma::topTable(fit_l, coef = x, number=Inf))

	return(list(fit=fit_l, tt=tt))
}


