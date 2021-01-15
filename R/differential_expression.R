#' differential_expression
#'
#' @import limma
#' @export

differential_expression <- function(dge=NULL, design=NULL, contr_mat=NULL, out_dir="./"){

	if(is.null(contr_mat)){
		cat("all-pairs contrasts...\n")
		contr_mat <- apply(t(combn(colnames(design), 2)), 1, paste, collapse="-")
		contr_mat <- limma::makeContrasts(contrasts = contr_mat, levels = design)
	}
	print(contr_mat)

	v <- limma::voom(dge, design, plot=TRUE)
	fit_l <- limma::lmFit(v, design)
	fit_lc <- limma::contrasts.fit(fit_l, contr_mat)
	fit_lc <- limma::eBayes(fit_lc)

	tt <- lapply(1:ncol(fit_lc$coefficients), function(x) limma::topTable(fit_lc, coef = x, number=Inf))
	names(tt) <- colnames(fit_lc$coefficients)

	for(i in 1:length(tt)){

		tt[[i]] <- merge(tt[[i]], fit_l$coefficients, by=0, sort=F)
		rownames(tt[[i]]) <- tt[[i]]$Row.names
		tt[[i]]$Row.names <- NULL

		write.table(tt[[i]], file=paste0(out_dir, "/degs_", names(tt), ".txt"), sep="\t", row.names = F)

		jpeg(paste0(out_dir, "/volcano_", names(tt), ".jpg"), width = 180, height = 180, res=300, units="mm")

		par(mar=c(6, 4, 4, 1))
		plot(tt[[i]]$logFC, -log10(tt[[i]]$P.Value), pch=16, cex=0.7, xlab="logFC", ylab="-log10(p)", col="purple")
		abline(v=c(-log2(1.5), log2(1.5)), h=-log10(0.01), lty=2)

		dev.off()
	}


	return(list(fit=fit_l, tt=tt))

}


