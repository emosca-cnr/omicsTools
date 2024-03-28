#' Filter and normalize count data
#' @param X genes-by-samples count matrix
#' @param out_dir output directory
#' @param col vector of colors; this must match the number of columns of X
#' @param min.cpm mininum count per million a gene must have in at least min_sample_perc percentage of samples
#' @param min.counts minimum number of counts
#' @param min.samples minimum number of samples
#' @param prior.count prior.count for cpm() function
#' @param subset.remove these samples will be removed
#' @param filter.by.cpm if TRUE filter by cmp rather then raw counts
#' @param width image width
#' @param height image height
#' @param res resolution
#' @param units image units
#' @return  filtered count matrix
#' @export 
#' @import edgeR grDevices pcaMethods pals DESeq2
#' @importFrom plotrix thigmophobe.labels
#' @importFrom vsn meanSdPlot
#' @importFrom graphics abline barplot boxplot hist layout legend lines par plot.new points
#' @importFrom stats density

filter_counts <- function(X=NULL, out_dir=NULL, col=NULL, width = 180, height=180, res=300, units="mm", prior.count=3, filter.by.cpm=TRUE, min.cpm=3, min.counts=10, subset.remove=NULL, min.samples=NULL){
  
  
  if(is.null(out_dir)){
    out_dir <- getwd()
  }else{
    dir.create(out_dir, recursive = T)
  }
  
  
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
    dge <- DGEList(X)
    cpm <- cpm(dge)
    cat("keeping genes with at least", min.cpm, " cpm in", min.samples, "samples...\n")
    keep <- rowSums(cpm >= min.cpm) >= min.samples
  }else{
    cat("keeping genes with at least", min.counts, "counts in", min.samples, "samples...\n")
    keep <- rowSums(X >= min.counts) >= min.samples
  }
  print(table(keep))
  X_filt <- X[keep,]
  cat("Filtered counts", dim(X_filt), "\n")

  #log-cpm
  dge <- DGEList(X)
  lcpm <- cpm(dge, log = T, prior.count = prior.count) #for reproducibility with voom
  
  #log-cpm_filt
  dge <- DGEList(X_filt)
  lcpm_filt <- cpm(dge, log = T, prior.count = prior.count) #for reproducibility with voom
    
  #density plots of raw and filtered data
  jpeg(file.path(out_dir, "density.jpg"), width = width, height = height/2, res=res, units=units)
  par(mfrow=c(1, 2))
  par(oma=c(0, 0, 0, 4))
  
  dens <- apply(lcpm, 2, density)
  
  xlim <- summary(unlist(lapply(dens, function(x) x$x)))[c(1, 6)]
  ylim <- summary(unlist(lapply(dens, function(x) x$y)))[c(1, 6)]
  
  plot(dens[[1]]$x, dens[[1]]$y, ylim=ylim, xlim=xlim, xlab="log-cmp", ylab="d", main = "raw", type="l")
  for(i in 2:length(dens)){
    lines(dens[[i]]$x, dens[[i]]$y, col=i)
  }
  abline(v=0, lty=2)
  
  dens <- apply(lcpm_filt, 2, density)
  
  xlim <- summary(unlist(lapply(dens, function(x) x$x)))[c(1, 6)]
  ylim <- summary(unlist(lapply(dens, function(x) x$y)))[c(1, 6)]
  
  plot(dens[[1]]$x, dens[[1]]$y, ylim=ylim, xlim=xlim, xlab="log-cmp", ylab="d", main = "filtered", type="l")
  for(i in 2:length(dens)){
    lines(dens[[i]]$x, dens[[i]]$y, col=i)
  }
  abline(v=0, lty=2)
  
  par(xpd = TRUE)
  legend("right", inset=-0.4, legend = colnames(lcpm_filt), pch=16, col=col, cex=0.7, xpd=T, bty="n")
  
  
  dev.off()
  
  nS <- rowSums(sign(X_filt))
  nG <- colSums(sign(X_filt))
  
  png(file.path(out_dir, "library.size.genes.png"), width=width, height=height/2, units="mm", res=300)
  
  par(mfrow=c(1, 2))
  par(mgp=c(2, .5, 0))
  
  hist(nS, breaks = 20, xlab="n samples", ylab="n genes", main = "Gene occurrece in samples", cex.lab=0.6, cex = 0.6, cex.axis = 0.6, las=2, cex.main=0.6)
  hist(nG, breaks = 20, xlab="n genes", ylab="n samples", main = "Number of expressed genes", cex.lab=0.6, cex= 0.6, cex.axis = 0.6, las=2, cex.main=0.6)
  
  dev.off()
  
  png(file.path(out_dir, "library_sizes.ngenes.counts.png"), width=width, height=height*0.9, units="mm", res=300)
  
  par(mar=c(3, 3, .1, 5))
  par(mgp=c(1.5, .5, 0))
  
  plot(colSums(X_filt), nG, cex=2, col=col, xlab="# reads", ylab="# genes", pch=16, log="xy", cex.axis = 0.6, cex.lab=0.6)
  
  par(xpd = TRUE)
  legend("right", inset=-0.15, legend = colnames(X_filt), pch=16, col=col, cex=0.7, xpd=T, bty="n")

  dev.off()
  
  return(X_filt)
  
}
