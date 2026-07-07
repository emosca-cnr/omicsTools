#' Differential expression analysis based on limma
#' @param dge DGEList
#' @param design experimental design matrix
#' @param contr_mat contrast matrix; limma::makeContrasts()
#' @param out_dir output directory
#' @param ... further arguments to limma::voomLmFit()
#' @importFrom openxlsx createWorkbook addWorksheet writeData saveWorkbook
#' @importFrom limma contrasts.fit eBayes topTable
#' @importFrom edgeR voomLmFit
#' @importFrom utils  write.table
#' @return A list with:
#' \enumerate{
#'   \item fit; fit object
#'   \item tt: list of differential expression tables
#' }
#' @export

differential_expression_limma <- function(dge=NULL, design=NULL, contr_mat=NULL, out_dir=NULL, ...){

  if(is.null(out_dir)){
  	out_dir <- getwd()
  }else{
  	dir.create(out_dir, recursive = T)
  }
  
  fit_l <- voomLmFit(counts = dge, design = design, ...)

	# Contrasts
	if(is.null(contr_mat)){
	  fit_lc <- fit_l
	}else{
	  fit_lc <- contrasts.fit(fit = fit_l, contrasts = contr_mat)
	}
  
  #Bayes correction
	fit_lc <- eBayes(fit = fit_lc)

	### top table and p-value adjustment
	tt <- lapply(1:ncol(fit_lc$coefficients), function(x) topTable(fit_lc, coef = x, number=Inf))
	names(tt) <- colnames(fit_lc$coefficients)

	#add the ANOVA as 
	tt.F <-  topTable(fit_lc, number=Inf)
	
	tt.F <- merge(tt.F, fit_l$coefficients, by=0, sort=F)
	rownames(tt.F) <- tt.F$Row.names
	tt.F$Row.names <- NULL
	
	tt.F <- tt.F[order(tt.F$adj.P.Val, tt.F$P.Value), ]
	
	#save results to file
	write.table(tt.F, file=file.path(out_dir, "degs_ttF.txt"), sep="\t", row.names = F)
	
	wb <- createWorkbook()
	for(i in 1:length(tt)){

		tt[[i]] <- merge(tt[[i]], fit_l$coefficients, by=0, sort=F)
		rownames(tt[[i]]) <- tt[[i]]$Row.names
		tt[[i]]$Row.names <- NULL

		tt[[i]] <- tt[[i]][order(tt[[i]]$adj.P.Val, -abs(tt[[i]]$logFC)), ]
		
		#save results to file
		write.table(tt[[i]], file=file.path(out_dir, paste0("degs_", names(tt)[i], ".txt")), sep="\t", row.names = F)
		
		addWorksheet(wb, names(tt)[i])
		writeData(wb, names(tt)[i], tt[[i]])
	
	}
	
	##Add tt F
	addWorksheet(wb, "tt.F")
	writeData(wb, "tt.F", tt.F)

	saveWorkbook(wb, file.path(out_dir, "degs.xlsx"), TRUE)

	return(list(fit=fit_lc, tt=tt, tt.F=tt.F))

}
