#############################
## MuSiC -- Multi-subject Single Cell deconvolution
## Modified Date： 2019-6-6
############################
## Wang X., et al. 
## Bulk tissue cell type deconvolution with multi-subject single-cell expression reference. Nature Communications, 10(2019) :380.
## Details at https://xuranw.github.io/MuSiC/articles/MuSiC.html

CreateMusicObject <- function(sc_counts=NULL, sc_cell_type=NULL, bulk_counts=NULL){
	# ExpressionSet data format for MuSiC package
	if(!("Biobase" %in% rownames(installed.packages())) ) {
		if (!requireNamespace("BiocManager", quietly = TRUE))
		install.packages("BiocManager")
		BiocManager::install("Biobase")
	}# end fi
	require(Biobase)
	bulk_eset <- sc_eset <- NULL
	if(!is.null(bulk_counts)){
		num_bs <- ncol(bulk_counts)
		pData <- data.frame(sampleID = 1:num_bs, SubjectName = factor(colnames(bulk_counts), levels = colnames(bulk_counts)) )
		rownames(pData) = colnames(bulk_counts);
		metadata <- data.frame(labelDescription = c("Sample ID", "Subject Name"), row.names = c("sampleID", "SubjectName"))
		bulk_eset <- ExpressionSet(assayData = data.matrix(bulk_counts), phenoData = new("AnnotatedDataFrame", data = pData, varMetadata = metadata))
	}# end fi
  
	# ExpressionSet data format for MuSiC package
	if(!is.null(sc_counts)){
		num_sc <- ncol(sc_counts)
		pData <- data.frame(sampleID = 1:num_sc, SubjectName = factor(paste0("sc", 1:num_sc), levels = paste0("sc", 1:num_sc)), cellTypeID = sc_cell_type, cellType = sc_cell_type)

		rownames(pData) = colnames(sc_counts)
		metadata <- data.frame(labelDescription= c("Sample ID", "Subject Name", "Cell Type ID", "Cell Type Name"),
                           row.names=c("sampleID", "SubjectName","cellTypeID", "cellType"))
		sc_eset = ExpressionSet(assayData = data.matrix(sc_counts), phenoData = new("AnnotatedDataFrame", data = pData, varMetadata = metadata))
	}# end fi
	return(list(bulk_eset=bulk_eset, sc_eset=sc_eset) )
}# end func



DECMuSiC <- function(sc_counts, sc_cell_type, bulk_counts){
	## Install package
	if(!("MuSiC" %in% rownames(installed.packages())) ) {
		devtools::install_github("xuranw/MuSiC")
	}# end fi
	library(MuSiC)
	if(!("xbioc" %in% rownames(installed.packages())) ) {
		devtools::install_github("renozao/xbioc")
	}# end fi
	library(xbioc) # for pVar issue
	
	object <- CreateMusicObject(bulk_counts=bulk_counts, 
                                sc_counts=sc_counts, 
                                sc_cell_type=sc_cell_type)
	res_music <- music_prop(bulk.eset = object$bulk_eset, 
                         sc.eset = object$sc_eset,
                         clusters = 'cellType', 
                         samples = 'sampleID', 
                         select.ct = unique(sc_cell_type))
	#res_music$Est.prop.weighted
	
	return(res_music)
}# end func

######################################
### example code
# source("./dec_music.R")
# library(scBio)
# data(SCLabels)
# data(SCFlu)
# data(BulkFlu)
# data(SCCellSpace)

# BulkFluAbs <- exp(BulkFlu)-1
# SCFluAbs <- exp(SCFlu)-1

### run MuSiC func
# res <- DECMuSiC(sc_counts=SCFluAbs, sc_cell_type=SCLabels, bulk_counts=BulkFluAbs)




