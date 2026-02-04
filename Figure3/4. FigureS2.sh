#!/bin/bash
set -euo pipefail

############################################
# scMIC project
# Step 0. Download and organize scRNA-seq data 
#         and validate scMIC framework for 
#         MIC identification
# Dataset: GSE173958
############################################


# Project root directory
RootPath="/data/twang15/wsj/scMIC"
cd "${RootPath}"
cd Figure3/
conda activate scMIC_Figure3

############################################
# 1. FigureS2A
#    -Input:Normalized ScData.pca.rds during DataProcessing.sh
#    -Output:umap with MIC identification
############################################

Rscript -e '
    library(Seurat); 
    Data = readRDS("GSE173958/M2results/Primary/ScData.pca.rds");
    umap_embeddings <- Embeddings(Data, "umap");
    umap_embeddings = as.data.frame(umap_embeddings);    
    OPcor = read.table("GSE173958/M2results/Metastasis_Met/proximity_OT.csv",sep="\t",header=FALSE,check.names = F);
    umap_embeddings$Peritoneal = OPcor$V1
    OPcor = read.table("GSE173958/M2results/Metastasis_Liver/proximity_OT.csv",sep="\t",header=FALSE,check.names = F);
    umap_embeddings$Liver = OPcor$V1
    OPcor = read.table("GSE173958/M2results/Metastasis_Lung/proximity_OT.csv",sep="\t",header=FALSE,check.names = F);
    umap_embeddings$Lung = OPcor$V1
    umap_embeddings$MIC = ifelse(umap_embeddings$Peritoneal>5 | umap_embeddings$Liver>5 | umap_embeddings$Lung>5,"Y","N")  
    write.table(umap_embeddings, file = "GSE173958/M2results/umap.txt", sep = "\t", quote = FALSE, row.names = TRUE); 
'

############################################
# 2. FigureS2B
#    -Input:umap.txt during Figure S2A and clones.txt
#    -Output:clone distribution
############################################

Rscript -e '
    MIC = read.table("GSE173958/M2results/umap.txt",sep="\t",header=TRUE,check.names = F);
    Clones = read.table("GSE173958/GSM5283488_M2-PT-clones/clones.txt",sep="\t",header=TRUE,check.names = F);
    rownames(Clones) = Clones$Barcode;
    Clones = Clones[rownames(MIC),];
    MIC$Clones = Clones$cloneID;       
    tab = table(MIC[, c("MIC", "Clones")])
    print(c(tab[, "2"] / rowSums(tab),sum(tab[, "2"])/sum(tab)))
'

############################################
# 3. FigureS2C
#    -Input:Normalized ScData.pca.rds during DataProcessing.sh
#           umap.txt during Figure3A
#    -Output:markerexpression.txt
############################################

Rscript -e '
    library(Seurat);
    Data = readRDS("GSE173958/M2results/Primary/ScData.pca.rds");
    PrimaryMatrix=Data@assays$RNA@layers$data;
    colnames(PrimaryMatrix)=colnames(Data);
    rownames(PrimaryMatrix)=rownames(Data);
    genes=c("Epcam","Muc1","Cdh1", "Lgals4", "Ocln", "Ctse","Prrx1", "Ifitm1", "Ifitm3", "S100a4")
    MIC = read.table("GSE173958/M2results/umap.txt",sep="\t",header=TRUE,check.names = F);
    Data$MIC = ifelse(MIC$Peritoneal>5 | MIC$Liver>5 | MIC$Lung>5,"Y","N")
    MICindex = which(Data@meta.data$MIC == "Y")        
    MICcells=rownames(Data@meta.data[MICindex,]);
    nonMICcells=rownames(Data@meta.data[-MICindex,]);             
    for(gene in genes){
      ttest = t.test(PrimaryMatrix[gene,MICcells],PrimaryMatrix[gene,nonMICcells]);
      print(c(gene,ttest$p.value,as.numeric(ttest$estimate)));
    }             
    SigGeneMatrix = data.frame(gene="",nonMIC=0,MIC=0);
    SigGeneMatrix = SigGeneMatrix[-1,];
    for(gene in genes){
      SigGeneMatrix = rbind(SigGeneMatrix,data.frame(gene=rep(gene,length(nonMICcells)),nonMIC=PrimaryMatrix[gene,nonMICcells],
                                       MIC=c(PrimaryMatrix[gene,MICcells],rep("",length(nonMICcells)-length(MICcells)))));
    }              
    write.table(SigGeneMatrix, "GSE173958/M2results/markerexpression.txt", quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE);              
' 