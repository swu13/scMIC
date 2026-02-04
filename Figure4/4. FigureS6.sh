#!/bin/bash
set -euo pipefail

############################################
# scMIC project
# Step 0. Download and organize scRNA-seq data 
#         and construct scMIC framework for 
#         MIC identification
# Dataset: OMIX002487
############################################


# Project root directory
RootPath="/data/twang15/wsj/scMIC"
cd "${RootPath}"
cd Figure4/

############################################
# 1. FigureS6A
#    -Input:marker DimPlot during DataProcessing.sh
############################################


Rscript -e '
    library(Seurat)
    files <- Sys.glob("OMIX002487/*/*/ScData.pca.rds")
    df <- data.frame(file = files,sample = basename(dirname(dirname(files))),tissue = basename(dirname(files)),stringsAsFactors = FALSE)
    files_sel <- subset(df,sample %in% c("P2", "P3", "P4", "P5"))
    genes=c("EPCAM","KRT8","KRT19","KRT18");
    each=c("null","null","null",0,0)
    for(file in files_sel$file){
      Data=readRDS(file);
      matrix=Data@assays$RNA@layers$data
      rownames(matrix)=rownames(Data)
      colnames(matrix)=colnames(Data)
      OutputFold=dirname(file)
      cnvdata = read.table(paste0(OutputFold,"/_copykat_prediction.txt"),sep="\t",header=TRUE,check.names = F)
      Data$Epithelial = ifelse(rownames(Data@meta.data) %in% cnvdata$cell.names, "Epithelial", "Others")
      index1=rownames(Data@meta.data)[which(Data@meta.data$Epithelial=="Epithelial")]
      index2=rownames(Data@meta.data)[which(Data@meta.data$Epithelial=="Others")]
      each=rbind(each,c(file,"Num","gene",length(index1),length(index2)))
      for(gene in genes){
          each=rbind(each,c(file,gene,mean(matrix[gene,index1]),mean(matrix[gene,index2])))
      }
    }
    write.table(each[-1,], "OMIX002487/FigureS6.txt", quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
'

############################################
# 2. FigureS6B
#    -Input: cnvscore during DataProcessing.sh
############################################


Rscript -e '
    library(Seurat)
    files <- Sys.glob("OMIX002487/*/*/ScData.pca.rds")
    df <- data.frame(file = files,sample = basename(dirname(dirname(files))),tissue = basename(dirname(files)),stringsAsFactors = FALSE)
    files_sel <- subset(df,sample %in% c("P2", "P3", "P4", "P5"))
    genes=c("EPCAM","KRT8","KRT19","KRT18");
    each=c("null","null",0,0)
    for(file in files_sel$file){
      Data=readRDS(file);
      matrix=Data@assays$RNA@layers$data
      rownames(matrix)=rownames(Data)
      colnames(matrix)=colnames(Data)
      OutputFold=dirname(file)
      cnvdata = read.table(paste0(OutputFold,"/_copykat_prediction.txt"),sep="\t",header=TRUE,check.names = F)
      Data$Epithelial = ifelse(rownames(Data@meta.data) %in% cnvdata$cell.names, "Epithelial", "Others")
      Data = subset(Data,Epithelial=="Epithelial")
      tumordata = subset(cnvdata,copykat.pred=="aneuploid")
      Data$Tumor = ifelse(rownames(Data@meta.data) %in% tumordata$cell.names, "Tumor", "Others")
      cnvexp = read.table(paste0(OutputFold,"/_copykat_CNA_results.txt"),sep="\t",header=TRUE,check.names = F)
      cnvscore=data.frame(cnvscore=colMeans((cnvexp[4:dim(cnvexp)[2]])^2))
      rownames(cnvscore) <- sub("\\\\.1$", "-1", rownames(cnvscore)) #rownames(cnvscore) <- sub("\\.1$", "-1", rownames(cnvscore))
      Data = subset(Data,cells=rownames(cnvscore))
      cnvscore <- cnvscore[rownames(Data@meta.data), , drop = FALSE]
      Data$cnvscore = cnvscore$cnvscore
      index1=rownames(Data@meta.data)[which(Data@meta.data$Tumor=="Tumor")]
      index2=rownames(Data@meta.data)[which(Data@meta.data$Tumor=="Others")]
      each=rbind(each,c(file,"Num",length(index1),length(index2)))
      each=rbind(each,c(file,"cnvscore",mean(as.numeric(Data@meta.data[index1,"cnvscore"]),na.rm=TRUE),mean(as.numeric(Data@meta.data[index2,"cnvscore"]),na.rm=TRUE)))
    }
    write.table(each[-1,], "OMIX002487/FigureS6.cnv.txt", quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
'

############################################
# 3. FigureS6C
#    -Input: "OMIX002487/FigureS6.cnv.txt"
############################################