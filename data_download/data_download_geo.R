#####################################################################################
## gene epression data: mRNA
## Modified date: 2019-6-7
###INPUT
##StudyID: "GSE1009"
##GetGPL: "TRUE", "FALSE"
##Path: the path to store the data
###OUTPUT
## an expression set object
## website to find the studies: https://www.ncbi.nlm.nih.gov/geo/

GEODownloadRNASeq <- function(StudyID, GetGPL=FALSE, Path){
  if(!("GEOquery" %in% rownames(installed.packages())) ) {
	if (!requireNamespace("BiocManager", quietly = TRUE))
	install.packages("BiocManager")
	BiocManager::install("GEOquery")
  }# end fi
  if(!("SummarizedExperiment" %in% rownames(installed.packages())) ) {
	if (!requireNamespace("BiocManager", quietly = TRUE))
	install.packages("BiocManager")
	BiocManager::install("SummarizedExperiment")
  }# end fi
  library(GEOquery)
  library(SummarizedExperiment)
  eSet <- getGEO(StudyID, destdir = Path, getGPL = GetGPL)
    
  exprSet = exprs(eSet[[1]])
  pdata = pData(eSet[[1]])
    
  write.csv(exprSet, paste0(Path, "/", StudyID, "_exprSet.csv"))
  write.csv(pdata, paste0(Path, "/", StudyID, "_metadata.csv"))
  return(eSet)
}# end func


######################
## example code
download_geo <- GEODownloadRNASeq(StudyID = "GSE1009", Path=getwd())

## return expression set object, using functions to query the gene expression data and metadata
exprSet = exprs(download_geo[[1]])
pdata = pData(download_geo[[1]])							


				
