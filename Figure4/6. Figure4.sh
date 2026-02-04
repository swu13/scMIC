#!/bin/bash
set -euo pipefail

############################################
# scMIC project
# Step 0. Download and organize scRNA-seq data 
#         and construct scMIC framework for 
#         MIC identification
# Dataset: GSE178318, OMIX002487 and GSE277783
############################################


# Project root directory
RootPath="/data/twang15/wsj/scMIC"
cd "${RootPath}"
cd Figure4/

############################################
# 1. Figure4A
#    -Input:MIC identification of GSE178318 during DataProcessing.sh
#    -Output:MIC vs. clinical phenotype
############################################

Rscript -e '
    library(Seurat);      
    samples=c("COL07","COL15","COL16");
    MICratio=c(0); #MIC ratio
    Tn=c(0); #Tumor cell number
    for(i in 1:3){
      sample=samples[i];
      Primary <- readRDS(paste0("GSE178318/",sample,"/Primary/Tumour/ScData.pca.rds"))
      Data_OT=read.table(paste0("GSE178318/",sample,"/proximity_OT.csv"));
      MICratio=c(MICratio,length(which(Data_OT$V1>0))/dim(Data_OT)[1]);
      Tn=c(Tn,dim(Primary)[2])
    }
    MICratio=MICratio[-1];
    Tn=Tn[-1];
    print(paste0("MIC Cell Ratio: ",MICratio));
    print(paste0("Tumour cell number:",Tn));

    TNM1=c(4,2,4); #Tumor size: T
    TNM2=c(1,2,2); #Lymph node: N
    Lnum=c(1,3,2); #Liver metastasis number
    Lsize=c(2,7.9,3.5); #Liver metastasis maximum size
    
    MetastasisPhenotypes = c("Lnum","Lsize");
    MICs = c("MICratio","Tn","TNM1","TNM2");
    for(MIC in MICs){
      for(pheno in MetastasisPhenotypes){
        cortest  = cor.test(get(MIC), get(pheno), method = "spearman", exact = FALSE);
        print(c(MIC,pheno));
        print(cortest);
      }
    }
'

############################################
# 2. Figure4B
#    -Input:MIC identification of OMIX002487 during DataProcessing.sh
#    -Output:MIC vs. clinical phenotype
############################################

Rscript -e '
    library(Seurat);      
    samples=c("P2","P3","P4","P5");
    MICratio=c(0); #MIC ratio
    Tn=c(0); #Tumor cell number
    for(i in 1:length(samples)){
      sample=samples[i];
      Primary <- readRDS(paste0("OMIX002487/",sample,"/Primary/Tumour/ScData.pca.rds"))
      Data_OT=read.table(paste0("OMIX002487/",sample,"/proximity_OT.csv"));
      MICratio=c(MICratio,length(which(Data_OT$V1>0))/dim(Data_OT)[1]);
      Tn=c(Tn,dim(Primary)[2])
    }
    MICratio=MICratio[-1];
    Tn=Tn[-1];
    print(paste0("MIC Cell Ratio: ",MICratio));
    print(paste0("Tumour cell number:",Tn));

    TNM=c(3,3,4,2)
    #TNM1=c(4.8,3.6,4.4,3.3); #Tumor size: max T
    #TNM2=c(4.0,3.0,3.6,2.0); #Tumor size: min T
    Lsize1=c(1.19,4.9,2.36,6.6); #Liver metastasis maximum size
    #Lsize2=c(1.07,2.44,2.09,5.2); #Liver metastasis min size
    print(TNM)
    print(Lsize1)
'

############################################
# 3. Figure4C
#    -Input:MIC identification of GSE277783 during DataProcessing.sh
#    -Output:MIC vs. clinical phenotype
############################################

Rscript -e '
    library(Seurat);     
    MICratio=c(0); #MIC ratio
    for(i in 1:13){
      if(i!=12){
        folder=paste0("GSE277783/Pt-",i,"/")
        Primary <- readRDS(paste0(folder,"/Primary/Tumour/ScData.pca.rds"))
        metastasis_dirs <- list.files(path = folder,pattern = "^Metastasis*",full.names = TRUE)
        for(j in 1:length(metastasis_dirs)){
          Data_OT=read.table(paste0(metastasis_dirs[j],"/proximity_OT.csv"));
          if(j==1){
            sum_OT=Data_OT$V1
          }else{
            sum_OT=sum_OT+Data_OT$V1
          }
        }
        ave_sum_OT = sum_OT/sum(sum_OT)
        MICratio=c(MICratio,length(which(sum_OT>20))/ncol(Primary));
      }
    }
    MICratio=MICratio[-1];
    print(MICratio);

    Stage=c(2,4,4,4,4,4,4,4,4,4,2,3)
    Treatment=c("F","F","F","F","F","T","T","T","T","T","T","T")
    print(Stage)
    print(Treatment)
'