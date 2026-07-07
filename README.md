# omicsTools

An R package that contains tools for omics data analysis:

- aggregation of featureCounts output;
- aggregation of salmon output;
- filtering of counts;
- normalization of counts according to edgeR and DESeq2;
- differential expression using limma;
- Relative Log Expression;
- mapping among identifiers.

It requires the following packages:
```{r, include=TRUE, eval=FALSE}
DESeq2 (Bioconductor)
edgeR (Bioconductor)
ggplot2 (CRAN)
limma (Bioconductor)
oenxlsx (CRAN)
reshape2 (CRAN)
SummarizedExperiment (Bioconductor)
```

Once all dependencies are in place, omicsTools can be installed from github as follows:
```{r, include=TRUE, eval=FALSE}
devtools::install_github("emosca-cnr/omicsTools")
```

