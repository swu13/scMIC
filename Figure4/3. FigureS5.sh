#!/bin/bash
set -euo pipefail

############################################
# scMIC project
# Step 0. Download and organize scRNA-seq data 
#         and construct scMIC framework for 
#         MIC identification
# Dataset: GSE249057
############################################


# Project root directory
RootPath="/data/twang15/wsj/scMIC"
cd "${RootPath}"
cd Figure4/

############################################
# 1. FigureS5
#    -Input:marker DimPlot during DataProcessing.sh
############################################

Rscript -e '
    library(Seurat)
    files <- Sys.glob("GSE178318/*/*/ScData.pca.rds")
    df <- data.frame(file = files,sample = basename(dirname(dirname(files))),tissue = basename(dirname(files)),stringsAsFactors = FALSE)
    files_sel <- subset(df,sample %in% c("COL07", "COL15", "COL16"))
    genes=c("EPCAM","KRT8","KRT19","KRT18");
    each=c("null","null",0,0)
    for(i in 1:dim(files_sel)[1]){
      Data=readRDS(files_sel$file[i]);
      matrix=Data@assays$RNA@layers$data
      rownames(matrix)=rownames(Data)
      colnames(matrix)=colnames(Data)
      OutputFold=dirname(files_sel$file[i])
      cnvdata = read.table(paste0(OutputFold,"/maligant.txt"),sep="\t",header=TRUE,check.names = F)
      tumordata = subset(cnvdata,CNVtype=="malignant" | XGBtype=="malignant")
      Data$Tumor = ifelse(rownames(Data@meta.data) %in% tumordata$barcodes, "Tumor", "Others")
      index1=rownames(Data@meta.data)[which(Data@meta.data$Tumor=="Tumor")]
      index2=rownames(Data@meta.data)[which(Data@meta.data$Tumor=="Others")]
      each=rbind(each,c(files_sel$file[i],"Num",dim(tumordata)[1],ncol(Data)-dim(tumordata)[1]))
      for(gene in genes){
          each=rbind(each,c(files_sel$file[i],gene,mean(matrix[gene,index1]),mean(matrix[gene,index2])))
      }
    }
    write.table(each[-1,], "GSE178318/FigureS5.txt", quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
'