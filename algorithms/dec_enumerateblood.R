#############################
## enumerateblood
## Modified Date： 2019-6-7
############################
## Shannon C.P., et al. 
## Enumerateblood – an R package to estimate the cellular composition of whole blood from Affymetrix Gene ST gene expression profiles. BMC Genomics, 18(2017):43.
## Details at https://cistrome.shinyapps.io/timer/ (website)
##			  https://github.com/hanfeisun/TIMER (R)
if(!("enumerateblood" %in% rownames(installed.packages())) ) {
	devtools::install_github("cashoes/enumerateblood", build_vignettes=TRUE)
}# end fi
library(enumerateblood)

# load package
library(enumerateblood)

# load package data object (the trained model)
# enumerateblood is already trained model
data('enumerateblood')

# use the trained model to predict cell proportions
# from gene expression matrix x.
#
# bulk_counts is a p (probesets) * n (samples) matrix of (RMA) normalized 
# gene expression data (assayed on the Affymetrix Gene 1.0 or 
# 1.1 ST platform), such as that returned by the `oligo::rma()` 
# function
predict(enumerateblood, bulk_counts)