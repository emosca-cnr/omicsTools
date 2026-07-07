#' Filter and normalize count data
#'
#' @param X genes-by-samples count matrix
#' @param min.cpm mininum count per million a gene must have in at least min_sample_perc percentage of samples
#' @param min.cpm mininum count per million a gene must have in at least min_sample_perc percentage of samples
#' @param min.counts minimum number of counts, used only if filter.by.cpm is FALSE
#' @param min.samples minimum number of samples
#' @param subset.remove these samples will be removed
#' @param filter.by.cpm if TRUE filter by cmp rather then raw counts
#'
#' @return  filtered count matrix
#' @export 
#' @importFrom edgeR cpm

filter_counts <- function(X=NULL, filter.by.cpm=TRUE, min.cpm=3, min.counts=10, subset.remove=NULL, min.samples=NULL){
  
  X[is.na(X)] <- 0
  
  if(!is.null(subset.remove)){
    cat("removing", sum(colnames(X) %in% subset.remove), "samples.\n")
    X <- X[, !colnames(X) %in% subset.remove]
  }
  
  if(is.null(col)){
    col <- 1:ncol(X)
  }
  
  #### FILTERING
  cat("FILTERING...\n")
  #remove rows with all elements equal to 0
  cat("removing all-zero rows...\n")
  cat("\t", sum(rowSums(X)==0), "rows...\n")
  idx_zero <- which(rowSums(X)==0)
  if(length(idx_zero) > 0){
    X <-  X[-idx_zero, ]
  }
  
  
  if(is.null(min.samples)){
    min.samples <- round(ncol(X)*0.25)
  }

  if(filter.by.cpm){
    cpm <- cpm(X)
    cat("keeping genes with at least", min.cpm, " cpm in", min.samples, "samples...\n")
    keep <- rowSums(cpm >= min.cpm) >= min.samples
    rm(cpm)
  }else{
    cat("keeping genes with at least", min.counts, "counts in", min.samples, "samples...\n")
    keep <- rowSums(X >= min.counts) >= min.samples
  }
  print(table(keep))
  X <- X[keep,]
  cat("Filtered counts", dim(X), "\n")

 
  return(X)
  
}
