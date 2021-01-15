#' filter_and_normalize
#'
#'
#' @export
#' @import edgeR grDevices pcaMethods limma

filter_and_normalize <- function(genes_by_samples_matrix, annotation=NULL, samples=NULL, min_sample_perc=0.2, min_cpm=3, out_dir="./", class=NULL, pal=NULL, top_pca=500){

	genes_by_samples_matrix[is.na(genes_by_samples_matrix)] <- 0
	#remove rows with all elements equal to 0
	cat("removing all-zero rows...\n")
	cat("\t", sum(rowSums(genes_by_samples_matrix)==0), "rows...\n")
	idx_zero <- which(rowSums(genes_by_samples_matrix)==0)
	if(length(idx_zero) > 0){
		genes_by_samples_matrix <-  genes_by_samples_matrix[-idx_zero, ]
		annotation <- annotation[match(rownames(genes_by_samples_matrix), rownames(annotation)), , drop=F]
	}

	cat("DGEList...\n")
	dge <- edgeR::DGEList(genes_by_samples_matrix, genes = annotation)

	#log-cpm
	cpm <- edgeR::cpm(dge)
	lcpm <- edgeR::cpm(dge, log = T, prior.count = 0.5) #for reproducibility with voom

	cat("filtering...\n")
	N_min <- round(ncol(cpm) * min_sample_perc)
	idx_keep <- rowSums(cpm > min_cpm) >= N_min
	print(table(idx_keep))
	counts_filt <- dge$counts[which(idx_keep), ]

	dge_filt <- edgeR::DGEList(counts_filt, genes = dge$genes[which(idx_keep), ])

	#log-cpm_filt
	cpm_filt <- edgeR::cpm(dge_filt)
	lcpm_filt <- edgeR::cpm(dge_filt, log = T, prior.count = 0.5) #for reproducibility with voom

	#density plots of raw and normalized data
	grDevices::jpeg(paste0(out_dir, "/density.jpg"), width = 180, height = 90, res=300, units="mm")
	par(mfrow=c(1, 2))
	plot(density(lcpm[, 1]), ylim=c(0, 1), xlim=c(-5, 20), xlab="log-cmp", ylab="d", main = "raw")
	for(i in 2:ncol(lcpm)){
		lines(density(lcpm[, i]), col=i)
	}
	abline(v=0, lty=2)

	plot(density(lcpm_filt[, 1]), ylim=c(0, 1), xlim=c(-5, 20), xlab="log-cmp", ylab="d", main = "filtered")
	for(i in 2:ncol(lcpm_filt)){
		lines(density(lcpm_filt[, i]), col=i)
	}
	abline(v=0, lty=2)
	dev.off()

	#TMM NORMALIZATION
	cat("TMM normalization...\n")
	dge_filt_norm <- edgeR::calcNormFactors(dge_filt, method='TMM')
	print(head(dge_filt_norm$samples))

	#log-cpm_filt
	cpm_filt_norm <- edgeR::cpm(dge_filt_norm)
	lcpm_filt_norm <- edgeR::cpm(dge_filt_norm, log = T, prior.count = 0.5) #for reproducibility with voom

	#library sizes
	png(paste0(out_dir, "/library_sizes.png"), width=90, height=90, units='mm', res=300)
	par(mar=c(6, 5, 4, 1))
	barplot(dge_filt_norm$samples$lib.size, las=2, names.arg = rownames(dge_filt_norm$samples), cex.names = 0.45, cex.axis = 0.6, xlab="", ylab="#", main="library size")
	dev.off()

	#boxplots
	jpeg(paste0(out_dir, "/boxplot.jpg"), width = 180, height = 90, res=300, units="mm")
	par(mfrow=c(1, 2))
	par(mar=c(6, 4, 4, 1))

	boxplot(lcpm_filt, las=2, cex.axis=0.45, ylab="log-cpm", main="raw", pars=list(outcex=0.5))
	boxplot(lcpm_filt_norm, las=2, cex.axis=0.45, ylab="log-cpm", main="normalized", pars=list(outcex=0.5))

	dev.off()

	#RLE
	jpeg(paste0(out_dir, "/RLE.jpg"), width = 180, height = 100, res=300, units="mm")

	par(mar=c(8, 4, 1, 1))
	par(mfrow=c(1, 2))

	boxplot(RLE(lcpm_filt, robust = F), cex=0.7, outline = F, las=2, cex.axis=0.7, ylab="relative log abundance", main="log-cpm filtered")
	boxplot(RLE(lcpm_filt_norm, robust = F), cex=0.7, outline = F, las=2, cex.axis=0.7, ylab="relative log abundance", main="log-cpm filtered, TMM")

	dev.off()

	#PCA

	res_pca <- pcaMethods::pca(lcpm_filt_norm, center = TRUE, scale="uv")
	print(res_pca)

	jpeg(paste0(out_dir, "/pca.jpg"), width = 120, height = 100, res=300, units="mm")
	layout(matrix(c(1,2), nrow = 1, byrow = T), widths = c(0.85, 0.15))
	par(mar=c(4, 4, 1, 1))

	plot(res_pca@loadings, pch="")
	col <- pal[as.numeric(class)]
	text(res_pca@loadings[, 1], res_pca@loadings[, 2], labels = rownames(res_pca@loadings), cex=0.5, xpd=T, col=col)

	par(mar=c(0.1, 0.1, 0.1, 0.1))
	plot.new()
	legend("center", legend = levels(class), col=pal, cex=0.5, pch=16)

	dev.off()

	jpeg(paste0(out_dir, "/mds.jpg"), width = 120, height = 100, res=300, units="mm")

	layout(matrix(c(1,2), nrow = 1, byrow = T), widths = c(0.85, 0.15))
	par(mar=c(4, 4, 1, 1))
	limma::plotMDS(lcpm_filt_norm, col=pal[as.numeric(class)], cex=0.6, xpd=T)

	par(mar=c(0.1, 0.1, 0.1, 0.1))
	plot.new()
	legend("center", legend = levels(class), col=pal, cex=0.5, pch=16)

	dev.off()


	return(list(dge=dge_filt_norm, lcpm=lcpm_filt_norm))

}
