#' Aggregate count files produced by featureCounts
#' @param count_files vector of count file names
#' @param sample_names optional sample names.
#'
aggregate_fCounts_res <- function(count_files=NULL, sample_names=NULL){
  
  print(count_files)
  print(sample_names)
  
  counts_data <- setNames(vector("list", length = length(count_files)), count_files)
  
  if(is.null(sample_names)){
    sample_names <- count_files
  }
  
  for(i in 1:length(count_files)){
    cat("Reading", count_files[i], "...\n")
    counts_data[[i]] <- read.delim(count_files[i], stringsAsFactors = F, skip = 1)
    colnames(counts_data[[i]])[7] <- sample_names[i]
    
  }
  
  cat("Merging files...\n")
  
  counts <- merge(unique(counts_data[[1]][, c(1, 7)]), unique(counts_data[[2]][, c(1, 7)]), by=1, all=TRUE, sort=F)
  if(length(counts_data)>2){
    
    for(i in 3:length(counts_data)){
      counts <- merge(counts, unique(counts_data[[i]][, c(1, 7)]), by=1, all=TRUE, sort=F)
    }
    
  }
  return(counts)
}