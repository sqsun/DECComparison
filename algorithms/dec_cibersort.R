#############################
## CIBERSORT -- Cell-type Identification By Estimating Relative Subsets Of RNA Transcripts 
## Date： 2019-6-6
############################
## Frishberg A, et al. 
## Robust enumeration of cell subsets from tissue expression profiles. Nature Methods, 12(2015) : 453–457.
## Details at http://cibersort.stanford.edu/

## Data already have, 4 steps as follows
# sc_counts
# sc_cell_type
# bulk_counts

## 1. create mixture file
write.table(bulk_counts, file="~/bulk_counts.txt",col.names=T,row.names=T,quote=F,sep="\t")

## 2. create reference sample file
write.table(sc_counts, file="~/sc_counts.txt",col.names=T,row.names=T,quote=F,sep="\t")

## 3. create phenotype classes file
cellType <- unique(sc_cell_type)
class_file <- matrix(2, ncol=ncol(sc_counts),nrow=length(cellType))
for(ict in 1:length(cellType)){
	class_file[ict, which(sub_sc_celltype==cellType[ict])] <- 1
}
rownames(class_file) <- cellType
write.table(class_file, file="~/sc_classes.txt", col.names=F,row.names=T,quote=F, sep="\t")

## 4. go to website, upload these three files, and click Run
https://cibersort.stanford.edu/runcibersort.php




