#############################
## xCell -- cell types enrichment analysis
## Modified Date： 2019-6-6
############################
## Aran et al.,
## xCell: digitally portraying the tissue cellular heterogeneity landscape. Genome Biology, 18(2017):220.
## Details at http://xcell.ucsf.edu/ (website)
##            https://github.com/dviraran/xCell (R)

## Install package
if(!("xCell" %in% rownames(installed.packages())) ) {
	devtools::install_github("dviraran/xCell", build_vignettes=TRUE)
}# end fi
library(xCell)

exprMatrix <- read.table(expr_file,header=TRUE,row.names=1, as.is=TRUE)
xCellAnalysis(exprMatrix)


