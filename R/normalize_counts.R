#' Filter and normalize count data
#' @param X genes-by-samples count matrix
#' @param gene.annotation data.frame with gene gene.annotations
#' @param sample.annotation data.frame with sample gene.annotations
#' @param col.by a column of sample.annnotatio to color plots. it is coerced to factor
#' @param prior.count prior.count required by cpm() function
#' @param out_dir output directory
#' @param pal color palette; the numbers of colors must match factor levels indicated by col_by
#' @param nFeatures number of features that will be selected by descresing sd as input for PCA
#' @param width image width
#' @param height image height
#' @param res resolution
#' @param units image units
#' @param pca.text whether to print sample labels
#' @param design experimental design matrix
#' @return a list with the following elements
#'  \enumerate{
#'           \item dge DGEList
#'           \item lcpm log counts per million matrix
#'           \item dds DESeq2 object
#'           \item vsd DESeq2 transform object, resulting from vst(, blind=TRUE) for exploratory analysis
#'           \item vsdbf DESeq2 transform object, resulting from vst(..., blind=FALSE) for downstream analysis
#'           \item res_pca resault of PCA over lcpm
#'           \item res_pca_vst resault of PCA over values returned by vst(, blind=TRUE)
#'           \item var_feat features used to obtain res_pca
#'           \item var_feat_vst features used to obtain res_pca_vst
#'        }
#' @export 
#' @import edgeR grDevices pcaMethods pals
#' @importFrom plotrix thigmophobe.labels
#' @importFrom vsn meanSdPlot
#' @importFrom DESeq2 DESeqDataSetFromMatrix estimateSizeFactors estimateDispersions vst
#' @importFrom limma plotMDS
#' @importFrom SummarizedExperiment assay
#' @importFrom graphics abline barplot boxplot hist layout legend lines par plot.new points
#' @importFrom utils combn read.delim write.table
#' @importFrom stats sd setNames model.matrix

normalize_counts <- function(X=NULL, gene.annotation=NULL, col.by="condition", sample.annotation=NULL, out_dir="./", pal=NULL, nFeatures=2000, width = 180, height=180, res=300, units="mm", prior.count=3, pca.text=TRUE, design=NULL){
  
  dir.create(out_dir, recursive = T)
  
  if(is.null(X) | is.null(sample.annotation)){
    stop("X and sample.annotation are mandatory.\n")
  }
  
  do.blind.false <- TRUE
  if(is.null(design)){
    design <- model.matrix(~1, data = sample.annotation)
    cat("no design matrix provided, using ~ 1 and setting blind to TRUE.\n")
    blind <- TRUE
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
  
  if(!col.by %in% colnames(sample.annotation)){
    cat(col.by, "not found among colnames(sample.annotation).\n")
  }
  
  condition <- factor(setNames(sample.annotation[, col.by], rownames(sample.annotation)))
  
  if(is.null(pal)){
    pal <- alphabet(length(levels(condition)))
  }
  
  dds <- DESeqDataSetFromMatrix(countData = X, colData = sample.annotation, design=design) #no design
  
  dge <- DGEList(X, genes = gene.annotation, samples = sample.annotation)
  
  #log-cpm
  lcpm <- cpm(dge, log = T, prior.count = prior.count) #for reproducibility with voom
  
  #TMM NORMALIZATION
  cat("TMM normalization...\n")
  dge <- calcNormFactors(dge, method='TMM')

  #log-cpm_filt
  lcpm <- cpm(dge, log = T, prior.count = prior.count)
  
  ###DESEq
  dds <- estimateSizeFactors(dds)
  dds <- estimateDispersions(dds)
  
  vsd <- vst(dds, blind = TRUE)
  vsdbf <- NULL
  if(do.blind.false){
    vsdbf <- vst(dds, blind = FALSE)
  }
  
  
  jpeg(paste0(out_dir, "/meansd.vst.jpg"), width = width, height = height, res=res, units=units)
  meanSdPlot(assay(vsd))
  dev.off()
  
  jpeg(paste0(out_dir, "/meansd.logcpmt.jpg"), width = width, height = height, res=res, units=units)
  meanSdPlot(lcpm)
  dev.off()
  
  #library sizes
  col <- pal[as.numeric(condition)]
  png(paste0(out_dir, "/library_sizes.png"), width=width, height=height, units=units, res=res)
  par(mar=c(6, 5, 4, 1))
  par(mgp=c(2, .5, 0))
  barplot(dge$samples$lib.size, las=2, names.arg = rownames(dge$samples), col=col, cex.lab=0.6, cex.names = 0.45, cex.axis = 0.6, xlab="", ylab="#", main="Library size")
  dev.off()
  
  #RLE
  jpeg(paste0(out_dir, "/RLE.jpg"), width = width, height = height/2, res=res, units=units)
  #
  par(mar=c(6, 3, 1, 1))
  par(mgp=c(2, .5, 0))
  par(mfrow=c(1, 2))
  
  boxplot(RLE(lcpm, robust = F), cex=0.7, col=col, outline = F, las=2, cex.axis=0.7, ylab="RLE", main="log-cpm + TMM", cex.lab=0.6)
  boxplot(RLE(assay(vsd), robust = F), cex=0.7, col=col, outline = F, las=2, cex.axis=0.7, ylab="RLE", main="sizeFactor + vst", cex.lab=0.6)
  
  dev.off()
  
  
  #PCA
  var_feat <- NULL
  if(length(nFeatures)>0){
    
    var_feat <- apply(lcpm, 1, sd)
    var_feat <- order(-var_feat)[1:(min(length(var_feat), nFeatures))]
    
    var_feat_vst <- apply(assay(vsd), 1, sd)
    var_feat_vst <- order(-var_feat_vst)[1:(min(length(var_feat_vst), nFeatures))]
  }
  
  res_pca <- pca(t(lcpm), center = TRUE, scale="uv", subset = var_feat)
  print(res_pca)
  
  res_pca_vst <- pcaMethods::pca(t(assay(vsd)), center = TRUE, scale="uv", subset = var_feat_vst)
  print(res_pca_vst)
  
  
  jpeg(paste0(out_dir, "/pca.jpg"), width = width, height = height/2, res=res, units=units)
  layout(matrix(c(1:3), nrow = 1, byrow = T), widths = c(0.42, 0.42, 0.16))
  
  par(mar=c(3, 3, 1, 2))
  par(mgp=c(2, .5, 0))
  
  plot(res_pca@scores, pch=21, bg=col, cex.axis=0.6, cex.lab=0.6)
  if(pca.text){
    par(xpd=TRUE)
    plotrix::thigmophobe.labels(res_pca@scores[, 1], res_pca@scores[, 2], labels = rownames(res_pca@scores), cex=0.5, xpd=T, col=col, font=2)
  }
  
  plot(res_pca_vst@scores, pch=21, bg=col, cex.axis=0.6, cex.lab=0.6)
  if(pca.text){
    par(xpd=TRUE)
    thigmophobe.labels(res_pca_vst@scores[, 1], res_pca_vst@scores[, 2], labels = rownames(res_pca_vst@scores), cex=0.5, xpd=T, col=col, font=2)
  }
  
  par(mar=c(0.1, 0.1, 0.1, 0.1))
  plot.new()
  legend("center", legend = levels(condition), col=pal, cex=0.5, pch=16)
  
  dev.off()
  
  jpeg(paste0(out_dir, "/mds.jpg"), width = width, height = height, res=res, units=units)
  
  layout(matrix(c(1, 2), nrow = 1, byrow = T), widths = c(0.85, 0.15))
  par(mar=c(3, 3, 1, 1))
  par(mgp=c(2, .5, 0))
  plotMDS(lcpm, col=col, cex=0.6, xpd=T)
  
  par(mar=c(0.1, 0.1, 0.1, 0.1))
  plot.new()
  legend("center", legend = levels(condition), col=pal, cex=0.5, pch=16)
  
  dev.off()
  
  return(list(dge=dge, lcpm=lcpm, dds=dds, vsd=vsd, vsdbf=vsdbf, res_pca=res_pca, res_pca_vst=res_pca_vst, var_feat=var_feat, var_feat_vst=var_feat_vst))
  
}
