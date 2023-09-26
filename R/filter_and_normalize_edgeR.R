#' Filter and normalize count data
#' @param X genes-by-samples count matrix
#' @param annotation data.frame with gene annotations
#' @param min_sample_perc minimum percentage of samples a gene must be expressed with min_cpm counts per million
#' @param min_cpm mininum count per million a gene must have in at least min_sample_perc percentage of samples
#' @export 
#' @import edgeR grDevices pcaMethods limma pals
#' @importFrom plotrix thigmophobe.labels

filter_and_normalize_edgeR <- function(X=NULL, annotation=NULL, min_sample_perc=0.2, min_cpm=3, out_dir="./", class=NULL, pal=NULL, nFeatures=2000, width = 180, height=180, res=300, units="mm", voom.plot=TRUE, prior.count=3, design=NULL){

  
  dir.create(out_dir, recursive = T)
  
	X[is.na(X)] <- 0
	
	if(is.null(class)){
	  stop("class is mandatory.\n")
	}
	if(length(class) != ncol(X)){
	  stop("lenght of class must be equal to ncol(X)\n.")
	}
	class <- as.factor(class)
	
	if(!identical(rownames(X), rownames(annotation))){
	  stop("rownames(X) must be identical to rownames(annotation).\n")
	}
	
	if(is.null(pal)){
	  pal <- pals::alphabet(length(levels(class)))
	}
	
	#remove rows with all elements equal to 0
	cat("removing all-zero rows...\n")
	cat("\t", sum(rowSums(X)==0), "rows...\n")
	idx_zero <- which(rowSums(X)==0)
	if(length(idx_zero) > 0){
		X <-  X[-idx_zero, ]
		annotation <- annotation[match(rownames(X), rownames(annotation)), , drop=F]
	}

	cat("DGEList...\n")
	dge <- edgeR::DGEList(X, genes = annotation)

	#log-cpm
	cpm <- edgeR::cpm(dge)
	lcpm <- edgeR::cpm(dge, log = T, prior.count = prior.count) #for reproducibility with voom

	cat("filtering...\n")
	N_min <- round(ncol(cpm) * min_sample_perc)
	idx_keep <- rowSums(cpm > min_cpm) >= N_min
	print(table(idx_keep))
	counts_filt <- dge$counts[which(idx_keep), ]

	dge_filt <- edgeR::DGEList(counts_filt, genes = dge$genes[which(idx_keep), ])

	#log-cpm_filt
	cpm_filt <- edgeR::cpm(dge_filt)
	lcpm_filt <- edgeR::cpm(dge_filt, log = T, prior.count = prior.count) #for reproducibility with voom

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
	dge_filt_norm <- edgeR::calcNormFactors(dge_filt, method='TMM')
	print(dge_filt_norm$samples[1:min(20, nrow(dge_filt_norm$samples)), ])

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
	

	#log-cpm_filt
	lcpm_filt_norm <- edgeR::cpm(dge_filt_norm, log = T, prior.count = prior.count) #for reproducibility with voom

	#library sizes
	col <- pal[as.numeric(class)]
	png(paste0(out_dir, "/library_sizes.png"), width=width, height=height, units=units, res=res)
	par(mar=c(6, 5, 4, 1))
	par(mgp=c(2, .5, 0))
	barplot(dge_filt_norm$samples$lib.size, las=2, names.arg = rownames(dge_filt_norm$samples), col=col, cex.lab=0.6, cex.names = 0.45, cex.axis = 0.6, xlab="", ylab="#", main="Library size")
	dev.off()

	#boxplots
	jpeg(paste0(out_dir, "/boxplot.jpg"), width = width, height = height/2, res=res, units=units)
	par(mfrow=c(1, 2))
	par(mar=c(6, 3, 2, 1))

	boxplot(lcpm_filt, las=2, col=col, cex.axis=0.45, ylab="log-cpm", main="raw", pars=list(outcex=0.5, outpch = 8), cex.lab=0.6)
	boxplot(lcpm_filt_norm, las=2, col=col, cex.axis=0.45, ylab="log-cpm", main="normalized", pars=list(outcex=0.5, outpch = 8), cex.lab=0.6)

	dev.off()

	#RLE
	jpeg(paste0(out_dir, "/RLE.jpg"), width = width, height = height/2, res=res, units=units)
#
	par(mar=c(6, 3, 1, 1))
	par(mgp=c(2, .5, 0))
	par(mfrow=c(1, 2))

	boxplot(RLE(lcpm_filt, robust = F), cex=0.7, col=col, outline = F, las=2, cex.axis=0.7, ylab="relative log abundance", main="log-cpm filtered", cex.lab=0.6)
	boxplot(RLE(lcpm_filt_norm, robust = F), cex=0.7, col=col, outline = F, las=2, cex.axis=0.7, ylab="relative log abundance", main="log-cpm filtered, TMM", cex.lab=0.6)

	dev.off()

	#PCA
	
	var_feat <- NULL
	if(length(nFeatures)>0){
	  var_feat <- apply(lcpm_filt_norm, 1, function(x) sd(x) / mean(x))
	  var_feat <- order(-var_feat)[1:(min(length(var_feat), nFeatures))]
	}
	
	res_pca <- pcaMethods::pca(t(lcpm_filt_norm), center = TRUE, scale="uv", subset = var_feat)
	print(res_pca)

	jpeg(paste0(out_dir, "/pca.jpg"), width = width, height = height, res=res, units=units)
	layout(matrix(c(1,2), nrow = 1, byrow = T), widths = c(0.85, 0.15))
	
	par(mar=c(3, 3, 1, 2))
	par(mgp=c(2, .5, 0))
	
	plot(res_pca@scores, pch=21, bg=col, cex.axis=0.6, cex.lab=0.6)
	par(xpd=TRUE)
	plotrix::thigmophobe.labels(res_pca@scores[, 1], res_pca@scores[, 2], labels = rownames(res_pca@scores), cex=0.5, xpd=T, col=col, font=2)

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

	return(list(dge=dge_filt_norm, lcpm=lcpm_filt_norm, res_pca=res_pca, var_feat=var_feat))

}
