#!/bin/bash
set -euo pipefail

############################################
# scMIC project
# Step 0. Download and organize scRNA-seq data 
#         and construct scMIC framework for 
#         MIC identification
# Dataset: GSE277783
############################################


# Project root directory
RootPath="/data/twang15/wsj/scMIC"
cd "${RootPath}"

cd Figure4/

############################################
# 1. FigureS6A
#    -Input:marker DimPlot during DataProcessing.sh
#    - Epithelial marker comparison
############################################

Rscript -e '
    library(Seurat)
    crc_dirs <- list.files(path = "GSE277783/",pattern = "[A-E]$",full.names = TRUE)
    genes=c("EPCAM","KRT8","KRT19","KRT18");
    each=c("null","null",0,0)
    for(i in 1:length(crc_dirs)){
      if (grepl("Pt-12", crc_dirs[i]) || grepl("Pt-1B", crc_dirs[i])) {next}
      Num = substr(crc_dirs[i], nchar(crc_dirs[i]), nchar(crc_dirs[i]))
      subfolder = paste0("Metastasis",Num)
      if(Num == "A"){
        subfolder = "Primary"
      }
      OutputPath=paste0("GSE277783/",sub(".*_(Pt-\\\\d+).*", "\\\\1", crc_dirs[i])) #OutputPath=paste0("GSE277783/",sub(".*_(Pt-\\d+).*", "\\1", crc_dirs[i]))
      Data = readRDS(paste0(OutputPath,"/",subfolder,"/ScData.pca.rds"))
      refscore = read.table(paste0(OutputPath,"/",subfolder,"/refcnv.txt"),sep="\t",header=TRUE,check.names = F)
      malignantscore = read.table(paste0(OutputPath,"/",subfolder,"/maligant.txt"),sep="\t",header=TRUE,check.names = F)
      Data$CNVtype = "NA"
      Data@meta.data[malignantscore$barcodes,]$CNVtype = malignantscore$CNVtype
      Data$XGBtype = "NA"
      Data@meta.data[malignantscore$barcodes,]$XGBtype = malignantscore$XGBtype

      Data$Tumor = ifelse(Data$TumorMarker1>0 & Data$cnvscore>max(refscore) & Data$TumorMarker1!="NA" & Data$cnvscore!="NA" & 
                         (Data$CNVtype == "malignant" | Data$XGBtype == "malignant"),"T","N")
      matrix=Data@assays$RNA@layers$data
      rownames(matrix)=rownames(Data)
      colnames(matrix)=colnames(Data)
      index1=rownames(Data@meta.data)[which(Data$Tumor=="T")]
      index2=rownames(Data@meta.data)[which(Data$Tumor=="N")]
      each=rbind(each,c(paste0(OutputPath,"/",subfolder),"Num",length(index1),length(index2)))
      for(gene in genes){
        each=rbind(each,c(paste0(OutputPath,"/",subfolder),gene,mean(matrix[gene,index1]),mean(matrix[gene,index2])))
      }
    }
    write.table(each[-1,], file = "GSE277783/Tumormarker.txt", sep = "\t", quote = FALSE, row.names = TRUE);
'

############################################
# 2. FigureS6B
#    -Input:marker DimPlot during DataProcessing.sh
#    - cnvscore marker comparison
############################################

Rscript -e '
    library(Seurat)
    crc_dirs <- list.files(path = "GSE277783/",pattern = "[A-E]$",full.names = TRUE)
    each=c("null","null",0,0)
    for(i in 1:length(crc_dirs)){
      if (grepl("Pt-12", crc_dirs[i]) || grepl("Pt-1B", crc_dirs[i])) {next}
      Num = substr(crc_dirs[i], nchar(crc_dirs[i]), nchar(crc_dirs[i]))
      subfolder = paste0("Metastasis",Num)
      if(Num == "A"){
        subfolder = "Primary"
      }
      OutputPath=paste0("GSE277783/",sub(".*_(Pt-\\\\d+).*", "\\\\1", crc_dirs[i])) #OutputPath=paste0("GSE277783/",sub(".*_(Pt-\\d+).*", "\\1", crc_dirs[i]))
      Data = readRDS(paste0(OutputPath,"/",subfolder,"/ScData.pca.rds"))
      refscore = read.table(paste0(OutputPath,"/",subfolder,"/refcnv.txt"),sep="\t",header=TRUE,check.names = F)
      malignantscore = read.table(paste0(OutputPath,"/",subfolder,"/maligant.txt"),sep="\t",header=TRUE,check.names = F)
      Data$CNVtype = "NA"
      Data@meta.data[malignantscore$barcodes,]$CNVtype = malignantscore$CNVtype
      Data$XGBtype = "NA"
      Data@meta.data[malignantscore$barcodes,]$XGBtype = malignantscore$XGBtype

      Data$Tumor = ifelse(Data$TumorMarker1>0 & Data$cnvscore>max(refscore) & Data$TumorMarker1!="NA" & Data$cnvscore!="NA" & 
                         (Data$CNVtype == "malignant" | Data$XGBtype == "malignant"),"T","N")
      SubData = subset(Data,TumorMarker1>0)
      matrix=SubData@assays$RNA@layers$data
      rownames(matrix)=rownames(SubData)
      colnames(matrix)=colnames(SubData)
      index1=rownames(SubData@meta.data)[which(SubData$Tumor=="T")]
      index2=rownames(SubData@meta.data)[which(SubData$Tumor=="N")]
      each=rbind(each,c(paste0(OutputPath,"/",subfolder),"Num",length(index1),length(index2)))
      each=rbind(each,c(paste0(OutputPath,"/",subfolder),"cnvscore",mean(as.numeric(SubData@meta.data[index1,"cnvscore"]),na.rm=TRUE),
                                                                    mean(as.numeric(SubData@meta.data[index2,"cnvscore"]),na.rm=TRUE)))
    }
    write.table(each[-1,], file = "GSE277783/Tumorcnv.txt", sep = "\t", quote = FALSE, row.names = TRUE);
'

############################################
# 3. FigureS6C
#    -Input:marker DimPlot during DataProcessing.sh
#    - GSE277783/Tumorcnv.txt during FigureS6B
############################################