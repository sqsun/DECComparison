#############################
## CPM -- Cell Population Mapping
## Date： 2019-6-6
############################
## Frishberg A, et al. 
## Cell composition analysis of bulk genomics using single-cell data. Nature Methods, 16(2019) :327–332.
## Details at https://github.com/amitfrish/scBio
## The method allows log counts or counts
## Outputs:
## 		predicted - CPM predicted cell abundance matrix (The main result of CPM). Each row represents a sample and each column a single cell.
## 		cellTypePredictions	- CPM predicted cell-type abundance matrix. Each row represnts a sample and each column a single cell-type. This is calculated if quantifyTypes = TRUE.
##		confIntervals - A matrix containing the confidence iterval for each cell and sample. Each row represnts a sample and each column a single cell. 
##			This is calculated if calculateCI = TRUE.
##		numOfRuns - The number of deconvolution repeats preformed by CPM.

DECCPM <- function(sc_counts, sc_cell_type, bulk_counts, cell_space=NULL, num_core=1){
	## Install package
	if (!("scBio" %in% rownames(installed.packages())) ) {
		devtools::install_github("amitfrish/scBio")
	}# end fi
	library(scBio)
	## using tSNE to compute the cell space based on single cell data if cell_space is not provided 
	if(is.null(cell_space)){
		if (!("Rtsne" %in% rownames(installed.packages())) ) {
			install.package("Rtsne")
		}# end fi
		library(Rtsne)
		sc_norm_counts <- log(sc_counts + 1)
		res_tsne <- Rtsne(t(sc_norm_counts))
		cell_space <- res_tsne$Y
		rm(res_tsne)
		rm(sc_norm_counts)
	}# end fi
	
	
	## absolute prediction
	res <- CPM(sc_counts, sc_cell_type, bulk_counts, cell_space, no_cores = num_core)
	
	return(res)
}# end func

######################################
### example code

# source("./dec_cpm.R")
# library(scBio)
# data(SCLabels)
# data(SCFlu)
# data(BulkFlu)
# data(SCCellSpace)

# BulkFluAbs <- exp(BulkFlu)-1
# SCFluAbs <- exp(SCFlu)-1

### run CPM func
# res <- DECCPM(sc_counts=SCFluAbs, sc_cell_type=SCLabels, bulk_counts=BulkFluAbs, cell_space=SCCellSpace, num_core = 6)




