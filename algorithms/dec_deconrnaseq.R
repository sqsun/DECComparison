#############################
## DeconRNASeq
## Modified Date： 2019-6-7
############################
## Gong T. et al.,
## DeconRNASeq: a statistical framework for deconvolution of heterogeneous tissue samples based on mRNA-Seq data. Bioinformatics, 29(2013):1083-1088.
## Details at https://www.bioconductor.org/packages/release/bioc/vignettes/DeconRNASeq/inst/doc/DeconRNASeq.pdf (R)

## Install package
if (!("DeconRNASeq" %in% rownames(installed.packages())) ) {
	if(!requireNamespace("BiocManager", quietly = TRUE))
	install.packages("BiocManager")
	BiocManager::install("DeconRNASeq")
}# end fi
library(DeconRNASeq)

## multi_tissue: expression profiles for 10 mixing samples from
## multiple tissues
data(multi_tissue)
datasets <- x.data[,2:11]
signatures <- x.signature.filtered.optimal[,2:6]
# proportions <- fraction
# res <- DeconRNASeq(datasets, signatures, proportions, checksig=FALSE, known.prop = TRUE, use.scale = TRUE, fig = TRUE)
res <- DeconRNASeq(datasets, signatures, checksig=FALSE, known.prop = FALSE, use.scale = TRUE, fig = TRUE)


