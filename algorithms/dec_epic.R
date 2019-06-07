#############################
## EPIC -- Estimate the Proportion of Immune and Cancer cells from bulk gene expression data
## Modified Date： 2019-6-6
############################
## Racle et al.,
## Simultaneous enumeration of cancer and immune cell types from bulk tumor gene expression data. eLife, 6(2017):e26476.
## Details at http://epic.gfellerlab.org/ (website)
##            https://github.com/scvannost/epicpy(Python)
##            https://github.com/GfellerLab/EPIC (R)

## Install package
if(!("EPIC" %in% rownames(installed.packages())) ) {
	devtools::install_github("GfellerLab/EPIC", build_vignettes=TRUE)
}# end fi
library(EPIC)
if (!("scBio" %in% rownames(installed.packages())) ) {
	devtools::install_github("amitfrish/scBio")
}# end fi
library(scBio)

#out <- EPIC(bulk = bulkSamplesMatrix, reference = referenceCellsList)
out <- EPIC(bulk = bulkSamplesMatrix, reference = referenceCellsList, mRNA_cell = mRNA_cell_vector, sigGenes = sigGenes_vector)
## ?EPIC to see details
######################################
### example code
# library(scBio)
# data(SCLabels)
# data(SCFlu)
# data(BulkFlu)
# data(SCCellSpace)



