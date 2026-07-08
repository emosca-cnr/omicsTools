#' Filter and normalize count data
#' @param X genes-by-samples count matrix
#' @param gene.annotation optional data.frame with gene annotations; row names must be identical to the row names of X
#' @param sample.annotation optional data.frame with sample annotations; column names must be identical to the column names of X
#' @param design optional experimental design matrix
#' @return a list with the following elements
#'  \enumerate{
#'           \item dge DGEList
#'           \item lcpm log counts per million matrix
#'           \item dds DESeq2 object
#'           \item vsd DESeq2 transform object, resulting from vst(, blind=TRUE) for exploratory analysis
#'           \item vsdbf DESeq2 transform object, resulting from vst(..., blind=FALSE) for downstream analysis
#'        }
#' @export 
#' @importFrom edgeR calcNormFactors cpm DGEList
#' @importFrom DESeq2 DESeqDataSetFromMatrix estimateSizeFactors estimateDispersions vst
#' @importFrom SummarizedExperiment assay
#' @importFrom stats  model.matrix

normalize_counts <- function(X=NULL, gene.annotation=NULL, sample.annotation=NULL, design=NULL){
  
  if(is.null(sample.annotation)){
    sample.annotation <- data.frame(id=colnames(X), row.names = colnames(X))
  }
  if(is.null(gene.annotation)){
    gene.annotation <- data.frame(id=rownames(X), row.names = rownames(X))
  }
  
  
  do.blind.false <- TRUE
  if(is.null(design)){
    design <- model.matrix(~1, data = sample.annotation)
    cat("no design matrix provided, using ~ 1 and setting blind to TRUE.\n")
    do.blind.false<-FALSE
  }
  
  if(!is.null(gene.annotation)){
    if(!identical(rownames(X), rownames(gene.annotation))){
      stop("rownames(X) must be identical to rownames(gene.annotation).\n")
    }
  }
  
  if(!identical(colnames(X), rownames(sample.annotation))){
    stop("colnames(X) must be identical to rownames(sample.annotation).\n")
  }
  
 
  dds <- DESeqDataSetFromMatrix(countData = X, colData = sample.annotation, design=design) #no design
  
  dge <- DGEList(X, genes = gene.annotation, samples = sample.annotation)
  
  #log-cpm
  lcpm <- cpm(dge, log = T)
  
  #TMM NORMALIZATION
  cat("TMM normalization...\n")
  dge <- calcNormFactors(dge, method='TMM')

  #log-cpm_filt
  lcpm <- cpm(dge, log = T)
  names(dimnames(lcpm)) <- NULL #remove Tags and Samples
    
  ###DESEq
  dds <- estimateSizeFactors(dds)
  dds <- estimateDispersions(dds)
  
  vsd <- vst(dds, blind = TRUE)
  vsdbf <- NULL
  if(do.blind.false){
    vsdbf <- vst(dds, blind = FALSE)
  }
  
  
   return(list(dge=dge, lcpm=lcpm, dds=dds, vsd=vsd, vsdbf=vsdbf))
  
}
