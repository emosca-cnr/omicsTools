#' filter_and_normalize
#'
#'
#' @export
#' @import edgeR org.Hs.eg.db grDevices

filter_and_normalize <- function(genes_by_samples_matrix, annotation=NULL, samples=NULL, min_sample_perc=0.2, min_cpm=3){

	#remove rows with all elements equal to 0
	cat("removing all-zero rows...\n")
	cat("\t", sum(rowSums(genes_by_samples_matrix)==0), "rows...\n")
	idx_zero <- which(rowSums(genes_by_samples_matrix)==0)
	genes_by_samples_matrix <-  genes_by_samples_matrix[-idx_zero, ]

	cat("DGEList...\n")
	dge <- edgeR::DGEList(genes_by_samples_matrix, genes = annotation)

	#log-cpm
	cpm <- edgeR::cpm(dge)
	lcpm <- edgeR::cpm(dge, log = T, prior.count = 0.5) #for reproducibility with voom

	cat("filtering...\n")
	N_min <- round(ncol(cpm) * min_sample_perc)
	idx_keep <- rowSums(cpm > min_cpm) >= N_min
	table(idx_keep)
	counts_filt <- dge$counts[which(idx_keep), ]

	dge_filt <- edgeR::DGEList(counts_filt, genes = dge$genes[which(idx_keep), ])

	#log-cpm_filt
	cpm_filt <- edgeR::cpm(dge_filt)
	lcpm_filt <- edgeR::cpm(dge_filt, log = T, prior.count = 0.5) #for reproducibility with voom

	#density plots of raw and normalized data
	grDevices::jpeg("density.jpg", width = 180, height = 90, res=300, units="mm")
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
	head(dge_filt_norm$samples)

	#log-cpm_filt
	cpm_filt_norm <- edgeR::cpm(dge_filt_norm)
	lcpm_filt_norm <- edgeR::cpm(dge_filt_norm, log = T, prior.count = 0.5) #for reproducibility with voom

	#library sizes
	png('library_sizes.png', width=90, height=90, units='mm', res=300)
	par(mar=c(6, 5, 4, 1))
	barplot(dge_filt_norm$samples$lib.size, las=2, names.arg = rownames(dge_filt_norm$samples), cex.names = 0.45, cex.axis = 0.6, xlab="", ylab="#", main="library size")
	dev.off()

	#boxplots
	jpeg("boxplot.jpg", width = 180, height = 90, res=300, units="mm")
	par(mfrow=c(1, 2))
	par(mar=c(6, 4, 4, 1))

	boxplot(lcpm_filt, las=2, cex.axis=0.45, ylab="log-cpm", main="raw", pars=list(outcex=0.5))
	boxplot(lcpm_filt_norm, las=2, cex.axis=0.45, ylab="log-cpm", main="normalized", pars=list(outcex=0.5))

	dev.off()

	#RLE
	jpeg("RLE.jpg", width = 180, height = 100, res=300, units="mm")

	par(mar=c(8, 4, 1, 1))
	par(mfrow=c(1, 2))

	boxplot(RLE(lcpm_filt, robust = F), cex=0.7, outline = F, las=2, cex.axis=0.7, ylab="relative log abundance")
	boxplot(RLE(lcpm_filt_norm, robust = F), cex=0.7, outline = F, las=2, cex.axis=0.7, ylab="relative log abundance")

	dev.off()

	return(list(dge=dge_filt_norm, lcpm=lcpm_filt_norm))

}
