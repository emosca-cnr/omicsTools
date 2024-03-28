#' Create a table integrating the summaries of STAR and featureCounts
#' @param samples_dir directory with all samples; every sample is in a directory that contains "Log.final.out" and "counts.txt.summary"
#' @export
#' @importFrom pals brewer.set1
#' @importFrom utils write.table
#' @importFrom grDevices png
#' 
count_assigned_summary <- function(samples_dir=NULL){
  
  samples <- list.dirs(samples_dir, recursive = F)
  
  star_res <- sapply(samples, function(x) read_star_final_log(file.path(x, "Log.final.out")))
  colnames(star_res) <- gsub("[^-]+-([^_]+).+", "\\1", colnames(star_res))
  write.table(star_res, file = "mapped_reads.txt", col.names = NA)
  
  samples_names <- gsub("[^-]+-([^_]+).+", "\\1", colnames(feat_count_res))
  feat_count_res <- sapply(samples, function(x) {
    df <- read.delim(file.path(x, "counts.txt.summary"), row.names = 1)
    df <- setNames(df[, 1], rownames(df))
    return(df)
  })
  colnames(feat_count_res) <- gsub("[^-]+-([^_]+).+", "\\1", colnames(feat_count_res))
  write.table(feat_count_res, file = "assigned_reads.txt", col.names = NA)
  
  png("star_res.png", width = 200, height = 100, units="mm", res=300)
  
  par(mar=c(3, 3, .1, 8))
  par(mgp=c(2.5, .5, 0))
  
  barplot(star_res[-1, ], col = brewer.set1(nrow(star_res)-1), names.arg = samples_names, legend.text = T, las=2, cex.names = 0.7, cex.axis = 0.7, args.legend = list(x="right", inset=c(-.28), cex=0.5, xpd=T))
  
  dev.off()
  
  png("fc_res.png", width = 200, height = 100, units="mm", res=300)
  
  par(mar=c(3, 3, .1, 8))
  par(mgp=c(2.5, .5, 0))
  
  barplot(feat_count_res[-which(rownames(feat_count_res) == "Unassigned_MultiMapping"), ], col = brewer.set1(nrow(feat_count_res)-1), names.arg = samples_names, legend.text = T, las=2, cex.names = 0.7, cex.axis = 0.7, args.legend = list(x="right", inset=c(-.28), cex=0.5, xpd=T))
  
  dev.off()
  
}