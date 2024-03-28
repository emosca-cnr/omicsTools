#' Read the output of STAR
#' @param file "Log.final.out" 
#' @export
#' @importFrom stats setNames
#' 
 
read_star_final_log <- function(file=NULL){
  
  preproc_log <- readLines(file)
  
  for(i in 1:length(preproc_log)){
    
    if(grepl("Number of input reads", preproc_log[i])){
      input_reads <- as.numeric(gsub(".+\t(\\d+).*", "\\1", preproc_log[i]))
    }
    if(grepl("Uniquely mapped reads number", preproc_log[i])){
      unq_reads <- as.numeric(gsub(".+\t(\\d+).*", "\\1", preproc_log[i]))
    }
    if(grepl("Number of reads mapped to multiple loci", preproc_log[i])){
      ml_reads <- as.numeric(gsub(".+\t(\\d+).*", "\\1", preproc_log[i]))
    }
    if(grepl("Number of reads mapped to too many loci", preproc_log[i])){
      tm_reads <- as.numeric(gsub(".+\t(\\d+).*", "\\1", preproc_log[i]))
    }
    if(grepl("Number of reads unmapped: too many mismatches", preproc_log[i])){
      tmm_reads <- as.numeric(gsub(".+\t(\\d+).*", "\\1", preproc_log[i]))
    }
    if(grepl("Number of reads unmapped: too short", preproc_log[i])){
      ts_reads <- as.numeric(gsub(".+\t(\\d+).*", "\\1", preproc_log[i]))
    }
    if(grepl("Number of reads unmapped: other", preproc_log[i])){
      o_reads <- as.numeric(gsub(".+\t(\\d+).*", "\\1", preproc_log[i]))
    }
    if(grepl("Number of chimeric reads", preproc_log[i])){
      c_reads <- as.numeric(gsub(".+\t(\\d+).*", "\\1", preproc_log[i]))
    }
  }
  
  out <- setNames(c(input_reads, unq_reads, ml_reads, tm_reads, tmm_reads, ts_reads, o_reads, c_reads), c("TOTAL", "uniquely_mapped", "MULTI_MAPPING_multiple_loci", "MULTI_MAPPING_too_many_loci", "UNMAPPED_too_many_mismatches", "UNMAPPED_too_short", "UNMAPPED_other", "CHIMERIC"))
  
  return(out)
  
}

