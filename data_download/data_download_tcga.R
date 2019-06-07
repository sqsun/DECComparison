##########################################
## gene epression data: mRNA and miRNA
## Modified date: 2019-6-6
###INPUT
##CancerProject: "TCGA-LUAD"
##WorkflowType: "HTSeq - FPKM" (normalization), "HTSeq - Counts" (counts)
##DataType: "Gene Expression Quantification" (mRNA), miRNA gene quantification (miRNA)
##Path: the path to store the data
###OUTPUT
##DE expression data 
##survival information
##survival expression data
## website: https://portal.gdc.cancer.gov/

TCGADownloadRNASeq <- function(CancerProject, WorkflowType, DataType, Path){
  if(!("TCGAbiolinks" %in% rownames(installed.packages())) ) {
	if (!requireNamespace("BiocManager", quietly = TRUE))
	install.packages("BiocManager")
	BiocManager::install("TCGAbiolinks")
  }# end fi
  if(!("SummarizedExperiment" %in% rownames(installed.packages())) ) {
	if (!requireNamespace("BiocManager", quietly = TRUE))
	install.packages("BiocManager")
	BiocManager::install("SummarizedExperiment")
  }# end fi
  library(TCGAbiolinks)
  library(SummarizedExperiment)
  #create data
  DataDirectory <- paste0(Path, "/GDC/",gsub("-","_",CancerProject))
  FileNameData <- paste0(DataDirectory, "_", WorkflowType, ".rda")
  FileNameData <- gsub("-","_",FileNameData)
  FileNameData <- gsub(" ","",FileNameData)
  query <- GDCquery(project = CancerProject, data.category = "Transcriptome Profiling",
                    data.type = DataType, workflow.type = WorkflowType)
  samplesDown <- getResults(query,cols=c("cases"))
  dataSmTP <- TCGAquery_SampleTypes(barcode = samplesDown, typesample = "TP")
  dataSmNT <- TCGAquery_SampleTypes(barcode = samplesDown, typesample = "NT")
  queryDown <- GDCquery(project = CancerProject, data.category = "Transcriptome Profiling",
                        data.type = DataType, workflow.type = WorkflowType, 
                        barcode = c(dataSmTP, dataSmNT))
  GDCdownload(query = queryDown, directory = DataDirectory)
  dataPrep <- GDCprepare(query = queryDown, save = TRUE, 
                         directory =  DataDirectory, save.filename = FileNameData)
  
  ##load data
  load(FileNameData)
  exp <- assays(data)@listData[[1]]
  clin <- data.frame(colData(data))
  
  ##DE data
  norm_code <- substr(colnames(exp)[grep("-11A-", colnames(exp))], 1, 12)
  de_code <- colnames(exp)[which(substr(colnames(exp), 1, 12) %in% norm_code)]
  de_code <- de_code[c(grep("-11A-", de_code), grep("-01A-", de_code))]
  de_exp <- exp[, colnames(exp) %in% de_code]
  cat("TN vs NT: ", ncol(de_exp), "\n")
  
  ##survival data
  surv_clin <- clin[-grep("-11A-", rownames(clin)), ]
  cat("delete clin of normal: ", nrow(surv_clin), "\n")
  surv_clin <- clin[grep("-01A-", rownames(clin)), ]
  cat("delete clin of not 01A: ", nrow(surv_clin), "\n")
  samp_code <- substr(rownames(surv_clin), 1, 12)
  freq <- data.frame(table(samp_code))
  dup_idx <- which(samp_code %in% freq[which(freq[, 2] == 2), 1])
  dup <- data.frame(code = rownames(surv_clin)[dup_idx], idx = dup_idx)
  dup <- dup[order(dup$code), ]
  if (nrow(dup) != 0){
    surv_clin <- surv_clin[-dup[seq(1, nrow(dup), 2), 2], ]
    cat("delete duplicate samples: ", nrow(surv_clin), "\n")
  } 
  
  ##survival expression data
  surv_exp <- exp[, colnames(exp) %in% rownames(surv_clin)]
  
  ##
  system(paste0("rm -r ", Path, "/GDC"))
  res <- list(de_exp, surv_exp, surv_clin)
  
  return(res)
}# end func


#####
## example code
## The Cancer type can be found in https://usermanual.wiki/Pdf/Manual.1597166842.pdf, Table 1.
## "TCGA-***"
download_tcga <- TCGADownloadRNASeq(CancerProject = "TCGA-LUAD", WorkflowType="HTSeq - Counts",
									DataType="Gene Expression Quantification", Path=getwd())


		