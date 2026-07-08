#' Differential expression analysis based on limma
#' @param counts genes-by-samples matrix of counts
#' @param sample_ann data.frame with sample annotation; row.names musst be identical to the colnames of counts
#' @param design design matrix
#' @param contrasts contrasts
#' @param type algorithm for srhinkage see the documentation of lfcShrink()
#' @importFrom DESeq2 DESeqDataSetFromMatrix DESeq resultsNames results lfcShrink
#' @return A list with:
#' \enumerate{
#'   \item fit; fit object
#'   \item tt: list of differential expression tables
#' }
#' @export

differential_expression_DESeq <- function(counts=NULL, sample_ann=NULL, design=NULL, contrasts=NULL, type="apeglm"){

  if(is.null(out_dir)){
  	out_dir <- getwd()
  }else{
  	dir.create(out_dir, recursive = T)
  }
  
	dds_int <- DESeqDataSetFromMatrix(countData = counts, colData = sample_ann, design = design) #this gives error due to the presence of NA
	dds_int <- DESeq(dds_int) #differential expression analysis
	
	cat("Coefficients: ", resultsNames(dds_int), "\n") ##have a look at the coefficients
	
	res_int <- results(dds_int, contrast = contrasts, pAdjustMethod = "fdr")
	res_int <- lfcShrink(dds_int, contrast = contrast, res=res_int, type = type)
	

	return(list(dds=dds_int, res=res_int))

}
