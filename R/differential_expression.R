#' Differential_expression for count data
#' Differential expression analysis based on limma
#' @param dge DGEList
#' @param design experimental design matrix
#' @param contr_mat contrast matrix, given to makeContrasts()
#' @param out_dir output directory
#' @param top_genes number of top genes to label in the volcano plots
#' @param top_genes_column column with gene labels for volcano plots
#' @import openxlsx
#' @importFrom plotrix thigmophobe.labels
#' @importFrom limma makeContrasts lmFit voom contrasts.fit eBayes topTable
#' @importFrom utils combn read.delim write.table
#' @return A list with:
#' \enumerate{
#'   \item result of lmFit()
#'   \item tt list of differential expression tables given by TopTable()
#' }
#' @export

differential_expression <- function(dge=NULL, design=NULL, contr_mat=NULL, out_dir="./", top_genes=10, top_genes_column="symbol"){

  dir.create(out_dir, recursive = T)
  
  if(is.null(contr_mat)){
		cat("all-pairs contrasts...\n")
		contr_mat <- apply(t(combn(colnames(design), 2)), 1, paste, collapse="-")
		contr_mat <- makeContrasts(contrasts = contr_mat, levels = design)
	}
	print(contr_mat)

	v <- voom(dge, design, plot=FALSE)
	fit_l <- lmFit(v, design)
	fit_lc <- contrasts.fit(fit_l, contr_mat)
	fit_lc <- eBayes(fit_lc)

	tt <- lapply(1:ncol(fit_lc$coefficients), function(x) topTable(fit_lc, coef = x, number=Inf))
	names(tt) <- colnames(fit_lc$coefficients)

	wb <- createWorkbook()
	for(i in 1:length(tt)){

		tt[[i]] <- merge(tt[[i]], fit_l$coefficients, by=0, sort=F)
		rownames(tt[[i]]) <- tt[[i]]$Row.names
		tt[[i]]$Row.names <- NULL

		tt[[i]] <- tt[[i]][order(tt[[i]]$adj.P.Val, -abs(tt[[i]]$logFC)), ]
		
		#save results to file
		write.table(tt[[i]], file=paste0(out_dir, "/degs_", names(tt)[i], ".txt"), sep="\t", row.names = F)
		
		addWorksheet(wb, names(tt)[i])
		writeData(wb, names(tt)[i], tt[[i]])
		
		jpeg(paste0(out_dir, "/volcano_", names(tt)[i], ".jpg"), width = 180, height = 180, res=300, units="mm")

		par(mar=c(3, 3, 3, .1))
		par(mgp=c(2, 0.5, 0))
		
		plot(tt[[i]]$logFC, -log10(tt[[i]]$adj.P.Val), pch=16, cex=0.7, xlab="logFC", ylab="-log10(q)", col="black", main=names(tt)[i])
		abline(v=c(-log2(1.5), log2(1.5)), h=-log10(0.01), lty=2)
		idx_degs <- abs(tt[[i]]$logFC) > log2(1.5) & tt[[i]]$adj.P.Val < 0.05
		points(tt[[i]]$logFC[idx_degs], -log10(tt[[i]]$adj.P.Val)[idx_degs], pch=16, cex=0.7, xlab="logFC", ylab="-log10(p)", col="purple")
		
		if(!is.null(top_genes) & length(tt[[i]][, top_genes_column])>0){
		  top_genes_idx <- order(tt[[i]]$adj.P.Val, -abs(tt[[i]]$logFC))[1:top_genes]
		  thigmophobe.labels(tt[[i]]$logFC[top_genes_idx], -log10(tt[[i]]$adj.P.Val)[top_genes_idx], tt[[i]][top_genes_idx, top_genes_column], cex=0.7, font=4)
		}
		
		dev.off()
		
	}
	saveWorkbook(wb, paste0(out_dir, "/degs.xlsx"), TRUE)


	return(list(fit=fit_l, tt=tt))

}


