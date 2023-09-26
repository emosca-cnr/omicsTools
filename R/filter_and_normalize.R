#' Filter and normalize count data
#' @param X genes-by-samples count matrix
#' @param annotation data.frame with gene annotations
#' @param min_sample_perc minimum percentage of samples a gene must be expressed with min_cpm counts per million
#' @param min_cpm mininum count per million a gene must have in at least min_sample_perc percentage of samples
#' @export 
#' @import edgeR grDevices pcaMethods limma pals DESeq2
#' @importFrom plotrix thigmophobe.labels
#' @importFrom vsn meanSdPlot

filter_and_normalize <- function(X=NULL, annotation=NULL, phenotypes=NULL, out_dir="./", pal=NULL, nFeatures=2000, width = 180, height=180, res=300, units="mm", voom.plot=TRUE, prior.count=3, filter.by.cpm=TRUE, min.cpm=3, min.counts=10, design=NULL, blind=TRUE, min.samples=NULL, pca_text=TRUE){
  
  
  dir.create(out_dir, recursive = T)
  
  X[is.na(X)] <- 0
  
  
  
  if(!identical(rownames(X), rownames(annotation))){
    stop("rownames(X) must be identical to rownames(annotation).\n")
  }
  
  if(!identical(colnames(X), rownames(phenotypes))){
    stop("colnames(X) must be identical to rownames(phenotypes).\n")
  }
  
  if(!all(c("condition", "batch") %in% colnames(phenotypes))){
    stop("phenotypes must contain 'condition' and 'batch' columns.\n")
  }
  class <- as.factor(pheno$condition)
  
  if(is.null(pal)){
    pal <- pals::alphabet(length(levels(class)))
  }
  
  #### FILTERING
  cat("FILTERING...\n")
  #remove rows with all elements equal to 0
  cat("removing all-zero rows...\n")
  cat("\t", sum(rowSums(X)==0), "rows...\n")
  idx_zero <- which(rowSums(X)==0)
  if(length(idx_zero) > 0){
    X <-  X[-idx_zero, ]
    annotation <- annotation[match(rownames(X), rownames(annotation)), , drop=F]
  }
  
  if(is.null(min.samples)){
    min.samples <- min(table(phenotypes$condition))
  }
  
  if(filter.by.cpm){
    dge <- edgeR::DGEList(X, genes = annotation)
    cpm <- edgeR::cpm(dge)
    cat("keeping genes with at least", min.cpm, " cpm in", min.samples, "samples...\n")
    keep <- rowSums(cpm >= min.cpm) >= min.samples
   }else{
    cat("keeping genes with at least", min.counts, "counts in", min.samples, "samples...\n")
    keep <- rowSums(X >= min.counts) >= min.samples
  }
  print(table(keep))
  X <- X[keep,]
  annotation <- annotation[keep, ]
  
  
  nS <- rowSums(sign(X))
  nG <- colSums(sign(X))
  
  png(paste0(out_dir, "/library.size.genes.png"), width=width, height=height/2, units="mm", res=300)
  
  par(mfrow=c(1, 2))
  par(mgp=c(2.5, 1.5, 0))
  
  hist(nS, breaks = 20, xlab="n samples", ylab="n genes", main = "Gene occurrece in samples", cex.lab=0.6, cex = 0.6, cex.axis = 0.6, las=2, cex.main=0.6)
  hist(nG, breaks = 20, xlab="n genes", ylab="n samples", main = "Number of expressed genes", cex.lab=0.6, cex= 0.6, cex.axis = 0.6, las=2, cex.main=0.6)
  
  dev.off()
  
  dds <- DESeq2::DESeqDataSetFromMatrix(countData = X,
                                colData = phenotypes,
                                design= ~ batch + condition)
  
  dge <- edgeR::DGEList(X, genes = annotation)
  
  #log-cpm
  cpm <- edgeR::cpm(dge)
  lcpm <- edgeR::cpm(dge, log = T, prior.count = prior.count) #for reproducibility with voom
  
  #N_min <- round(ncol(cpm) * min_sample_perc)
  #idx_keep <- rowSums(cpm > min_cpm) >= N_min
  #print(table(idx_keep))
  
  #log-cpm_filt
  lcpm_filt <- edgeR::cpm(dge, log = T, prior.count = prior.count) #for reproducibility with voom
  
  #density plots of raw and filtered data
  grDevices::jpeg(paste0(out_dir, "/density.jpg"), width = width, height = height/2, res=res, units=units)
  par(mfrow=c(1, 2))
  plot(density(lcpm[, 1]), ylim=c(0, .4), xlim=c(-5, 20), xlab="log-cmp", ylab="d", main = "raw")
  for(i in 2:ncol(lcpm)){
    lines(density(lcpm[, i]), col=i)
  }
  abline(v=0, lty=2)
  
  plot(density(lcpm_filt[, 1]), ylim=c(0, .4), xlim=c(-5, 20), xlab="log-cmp", ylab="d", main = "filtered")
  for(i in 2:ncol(lcpm_filt)){
    lines(density(lcpm_filt[, i]), col=i)
  }
  abline(v=0, lty=2)
  dev.off()
  
  
  #TMM NORMALIZATION
  cat("TMM normalization...\n")
  dge_filt_norm <- edgeR::calcNormFactors(dge, method='TMM')
  print(dge_filt_norm$samples[1:min(20, nrow(dge_filt_norm$samples)), ])
  
  #log-cpm_filt
  lcpm_filt_norm <- edgeR::cpm(dge_filt_norm, log = T, prior.count = prior.count)
  
  ###DESEq
  dds <- estimateSizeFactors(dds)
  dds <- estimateDispersions(dds)
  vsd <- vst(dds, blind = blind)
  
  ### VOOM
  if(voom.plot){
    if(!is.null(design)){
      grDevices::jpeg(paste0(out_dir, "/voom.jpg"), width = width, height = height, res=res, units=units)
      limma::voom(counts = dge_filt_norm, design = design, plot=TRUE)
      dev.off()
    }else{
      cat("VOOM not used, because design matrix is NULL.\n")
    }
  }
  
  grDevices::jpeg(paste0(out_dir, "/meansd.vst.jpg"), width = width, height = height, res=res, units=units)
  vsn::meanSdPlot(assay(vsd))
  dev.off()
  
  grDevices::jpeg(paste0(out_dir, "/meansd.logcpmt.jpg"), width = width, height = height, res=res, units=units)
  meanSdPlot(lcpm_filt_norm)
  dev.off()
  
  #library sizes
  col <- pal[as.numeric(class)]
  png(paste0(out_dir, "/library_sizes.png"), width=width, height=height, units=units, res=res)
  par(mar=c(6, 5, 4, 1))
  par(mgp=c(2, .5, 0))
  barplot(dge_filt_norm$samples$lib.size, las=2, names.arg = rownames(dge_filt_norm$samples), col=col, cex.lab=0.6, cex.names = 0.45, cex.axis = 0.6, xlab="", ylab="#", main="Library size")
  dev.off()
  
  #RLE
  jpeg(paste0(out_dir, "/RLE.jpg"), width = width, height = height/3, res=res, units=units)
  #
  par(mar=c(6, 3, 1, 1))
  par(mgp=c(2, .5, 0))
  par(mfrow=c(1, 3))
  
  boxplot(RLE(lcpm_filt, robust = F), cex=0.7, col=col, outline = F, las=2, cex.axis=0.7, ylab="relative log abundance", main="log-cpm", cex.lab=0.6)
  boxplot(RLE(lcpm_filt_norm, robust = F), cex=0.7, col=col, outline = F, las=2, cex.axis=0.7, ylab="relative log abundance", main="log-cpm + TMM", cex.lab=0.6)
  boxplot(RLE(assay(vsd), robust = F), cex=0.7, col=col, outline = F, las=2, cex.axis=0.7, ylab="relative log abundance", main="sizeFactor + vst", cex.lab=0.6)
  
  dev.off()
  
  
  #PCA
  var_feat <- NULL
  if(length(nFeatures)>0){
    
    var_feat <- apply(lcpm_filt_norm, 1, sd)
    var_feat <- order(-var_feat)[1:(min(length(var_feat), nFeatures))]
    
    var_feat_vst <- apply(assay(vsd), 1, sd)
    var_feat_vst <- order(-var_feat_vst)[1:(min(length(var_feat_vst), nFeatures))]
  }
  res_pca <- pcaMethods::pca(t(lcpm_filt_norm), center = TRUE, scale="uv", subset = var_feat)
  print(res_pca)
  
  res_pca_vst <- pcaMethods::pca(t(assay(vsd)), center = TRUE, scale="uv", subset = var_feat_vst)
  print(res_pca_vst)
  
  
  jpeg(paste0(out_dir, "/pca.jpg"), width = width, height = height/2, res=res, units=units)
  layout(matrix(c(1:3), nrow = 1, byrow = T), widths = c(0.42, 0.42, 0.16))
  
  par(mar=c(3, 3, 1, 2))
  par(mgp=c(2, .5, 0))
  
  plot(res_pca@scores, pch=21, bg=col, cex.axis=0.6, cex.lab=0.6)
  if(pca_text){
    par(xpd=TRUE)
    plotrix::thigmophobe.labels(res_pca@scores[, 1], res_pca@scores[, 2], labels = rownames(res_pca@scores), cex=0.5, xpd=T, col=col, font=2)
  }
  
  plot(res_pca_vst@scores, pch=21, bg=col, cex.axis=0.6, cex.lab=0.6)
  if(pca_text){
    par(xpd=TRUE)
    plotrix::thigmophobe.labels(res_pca_vst@scores[, 1], res_pca_vst@scores[, 2], labels = rownames(res_pca_vst@scores), cex=0.5, xpd=T, col=col, font=2)
  }
  
  par(mar=c(0.1, 0.1, 0.1, 0.1))
  plot.new()
  legend("center", legend = levels(class), col=pal, cex=0.5, pch=16)
  
  dev.off()
  
  jpeg(paste0(out_dir, "/mds.jpg"), width = width, height = height, res=res, units=units)
  
  layout(matrix(c(1, 2), nrow = 1, byrow = T), widths = c(0.85, 0.15))
  par(mar=c(3, 3, 1, 1))
  par(mgp=c(2, .5, 0))
  limma::plotMDS(lcpm_filt_norm, col=col, cex=0.6, xpd=T)
  
  par(mar=c(0.1, 0.1, 0.1, 0.1))
  plot.new()
  legend("center", legend = levels(class), col=pal, cex=0.5, pch=16)
  
  dev.off()
  
  return(list(dge=dge_filt_norm, lcpm=lcpm_filt_norm, dds=dds, vsd=vsd, res_pca=res_pca, res_pca_vst=res_pca_vst, var_feat=var_feat, var_feat_vst=var_feat_vst))
  
}
