#' Aggregate quant files produced by Salmon
#' @param quant_files vector of count file names
#' @param sample_names optional sample names.
#' @importFrom stats setNames
#' @export

aggregate_salmon_quant <- function(quant_files=NULL, sample_names=NULL){
  
  print(quant_files)
  print(sample_names)
  
  quant_data <- setNames(vector("list", length = length(quant_files)), quant_files)
  
  if(is.null(sample_names)){
    sample_names <- quant_files
  }
  
  for(i in 1:length(quant_files)){
    cat("Reading", quant_files[i], "...\n")
    quant_data[[i]] <- read.delim(quant_files[i], stringsAsFactors = F)
    colnames(quant_data[[i]])[3:5] <- paste(sample_names[i], colnames(quant_data[[i]])[3:5], sep="_")
  }
  
  cat("Merging files...\n")
  
  ans <- list()
  ans$l <- unique(quant_data[[1]][, c(1, 2)]) #the lenght is equal among samples
  ans$le <- merge(unique(quant_data[[1]][, c(1, 3)]), unique(quant_data[[2]][, c(1, 3)]), by=1, all=TRUE, sort=F)
  ans$tpm <- merge(unique(quant_data[[1]][, c(1, 4)]), unique(quant_data[[2]][, c(1, 4)]), by=1, all=TRUE, sort=F)
  ans$x <- merge(unique(quant_data[[1]][, c(1, 5)]), unique(quant_data[[2]][, c(1, 5)]), by=1, all=TRUE, sort=F)

  if(length(quant_data)>2){
    
    for(i in 3:length(quant_data)){
      ans$le <- merge(ans$le, unique(quant_data[[i]][, c(1, 3)]), by=1, all=TRUE, sort=F)
      ans$tpm <- merge(ans$tpm, unique(quant_data[[i]][, c(1, 4)]), by=1, all=TRUE, sort=F)
      ans$x <- merge(ans$x, unique(quant_data[[i]][, c(1, 5)]), by=1, all=TRUE, sort=F)
    }
    
  }
  
  return(ans)
}