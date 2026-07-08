#' Differential expression analysis based on limma
#' @param counts genes-by-samples matrix of counts
#' @param sample_ann data.frame with sample annotation; row.names musst be identical to the colnames of counts
#' @param design design matrix
#' @param contrasts contrasts
#' @param type algorithm for srhinkage see the documentation of lfcShrink()
#' @importFrom DESeq2 DESeqDataSetFromMatrix DESeq resultsNames results lfcShrink
#' @return A list with:
#' \enumerate{
#'   \item dds; dds object
#'   \item res: results object
#' }
#' @export

differential_expression_DESeq <- function(counts=NULL, sample_ann=NULL, design=NULL, contrast=NULL, type="apeglm"){
	
	
	dds_int <- DESeqDataSetFromMatrix(countData = counts, colData = sample_ann, design = design) #this gives error due to the presence of NA
	dds_int <- DESeq(dds_int) #differential expression analysis
	
	cat("Coefficients: ", resultsNames(dds_int), "\n") ##have a look at the coefficients
	
	if(!is.null(contrast)){
		res_int <- results(dds_int, contrast = contrast, pAdjustMethod = "fdr")
		res_int <- as.data.frame(lfcShrink(dds_int, contrast = contrast, res=res_int, type = type))
	}else{
		res_int <- NULL
		cat("No contrast given\n")
	}	
	
	return(list(dds=dds_int, res=res_int))
	
}
