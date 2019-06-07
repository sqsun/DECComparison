#############################
## BSEQ-sc -- Deconvolution of Bulk Sequencing Experiments using Single Cell Data
## Modified Date： 2019-6-6
############################
## Baron M., et al. 
## A single-cell transcriptomic map of the human and mouse pancreas reveals inter- and intra-cell population structure. Cell Systems, 3(2016):346-360.
## Details at https://shenorrlab.github.io/bseqsc/vignettes/bseq-sc.html
## Two types of input data:
## 	bulk data obtained from RNA sequencing of mixed tissue samples: these are the data we want to deconvolve;
## 	data from single cell RNA sequencing: these serve to compute cell type-specific reference profiles that are used to estimate cell type proportions in the bulk data.


## Install package
if(!("bseqsc" %in% rownames(installed.packages())) ) {
	devtools::install_github('shenorrlab/bseqsc')
}# end fi
library(bseqsc)

if(!("xbioc" %in% rownames(installed.packages())) ) {
	devtools::install_github("renozao/xbioc")
}# end fi
	library(xbioc) # for pVar issue
	
	
## load bulk RNAseq data 
https://shenorrlab.github.io/bseqsc/data/GSE50244.rds

## load single-cell RNAseq data
https://shenorrlab.github.io/bseqsc/data/islet-eset.rds


############
##
eset <- readRDS('~/GSE50244.rds')
# for this analysis we only look at samples with hba1c data
eset <- droplevels(eset[, !is.na(eset$hba1c)])

eislet <- readRDS('~/islet-eset.rds')
eislet

## Building the reference basis matrix
B <- bseqsc_basis(eislet, pancreasMarkers, clusters = 'cellType', samples = 'sampleID', ct.scale = TRUE)

## Estimating proportions
## Cell type proportions are estimated using CIBERSORT, a method that was developed to investigate immune infiltrating landscape in tumour tissue (Newman et al. 2015; Gentles et al. 2015) (See section Requirements for setup instructions).

fit <- bseqsc_proportions(eset, B, verbose = TRUE)


