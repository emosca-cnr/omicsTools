#' filter_and_normalize
#' 
#' 
filter_and_normalize <- function(genes_by_samples_matrix){


library(edgeR)

dir <- "/DATA/Inflammation3/ouputRNAseqSTAR/count3/"
files <- list.files(dir, pattern = "RNA", full.names = T)
files <- files[!grepl("summary", files)]

#phenotype
pheno <- read.csv("../phenotypes2.txt", sep="\t", stringsAsFactors = F)
pheno

#raw count matrix
counts_raw <- create_count_matrix(files, sep="\t", skip=1)

#adjust column names
colnames(counts_raw) <- gsub("X.DATA.Inflammation3.ouputRNAseqSTAR.alignment.", "", colnames(counts_raw))
colnames(counts_raw) <- gsub("\\..+$", "", colnames(counts_raw))
colnames(counts_raw)
samp_names <- colnames(counts_raw)[-c(1:6)]
colnames(counts_raw) <- gsub("_[ATCG]+_L00.$", "", colnames(counts_raw))

samp_names <- data.frame(name_tag=samp_names, name=colnames(counts_raw)[-c(1:6)], stringsAsFactors = F)
samp_names$tag <- gsub("^.+_([ATCG]+_L\\d+)", "\\1", samp_names$name_tag)
pheno <- merge(pheno, samp_names, by.x=1, by.y=2)
write.table(pheno, "pheno_rnaseq.txt", sep="\t", row.names = F)

#write table
write.table(counts_raw, "counts_raw.txt", sep="\t", row.names = F)

#split annotation and counts
counts_ann <- counts_raw[, 1:6]
rownames(counts_raw) <- counts_raw$Geneid
counts_raw <- counts_raw[, -c(1:6)]

#remove rows with all elements equal to 0
sum(rowSums(counts_raw)==0)
idx_zero <- which(rowSums(counts_raw)==0)
counts0 <- counts_raw[idx_zero, ]
counts <-  counts_raw[-idx_zero, ]

#DFEList
dge <- DGEList(counts, genes = counts_ann[-idx_zero, ])

#log-cpm
cpm <- cpm(dge)
lcpm <- cpm(dge, log = T, prior.count = 0.5) #for reproducibility with voom

#filtering
N_min <- round(ncol(cpm)*0.1)
N_min <- 3
idx_keep <- rowSums(cpm > 1) >= N_min
table(idx_keep)
counts_filt <- counts[which(idx_keep), ]

dge_filt <- DGEList(counts_filt, genes = dge$genes[which(idx_keep), ])

#log-cpm_filt
cpm_filt <- cpm(dge_filt)
lcpm_filt <- cpm(dge_filt, log = T, prior.count = 0.5) #for reproducibility with voom

#density plots
jpeg("density.jpg", width = 180, height = 90, res=300, units="mm")
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
dge_filt_norm <- calcNormFactors(dge_filt, method='TMM')
head(dge_filt_norm$samples)

#log-cpm_filt
cpm_filt_norm <- cpm(dge_filt_norm)
lcpm_filt_norm <- cpm(dge_filt_norm, log = T, prior.count = 0.5) #for reproducibility with voom

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


save(pheno, dge_filt_norm, lcpm_filt_norm, counts_raw, file="normalization.RData", compress = "gzip")
}