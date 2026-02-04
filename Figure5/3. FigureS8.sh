#!/bin/bash
set -euo pipefail

############################################
# scMIC project
# Step 0. Download and organize scRNA-seq data 
#         and construct scMIC framework for 
#         MIC identification
# Dataset: GSE173958
############################################


# Project root directory
RootPath="/data/twang15/wsj/scMIC"
cd "${RootPath}"

mkdir FigureS/
cd FigureS/

############################################
# 1. FigureS8
#    -Input:MIC identification based on 
#           MIC score
############################################

Rscript -e '
    library(Seurat);
    Genes=read.table("../Figure5/MHMICgenes.txt",header=TRUE);
    colnames(Genes)=c("Mgenes","Hgenes","Hit","HumanOrNot")
    mousegenes=Genes$Mgenes
    Data = readRDS("../Figure3/GSE173958/M2results/Primary/ScData.pca.rds");
    subData=subset(Data, features = intersect(mousegenes,VariableFeatures(Data)));
    Matrix=subData@assays$RNA@layers$data;
    Meanvalue <- colMeans(Matrix);
    Data$MIC=ifelse(Meanvalue >= summary(Meanvalue)[5],"Y",ifelse(Meanvalue <= summary(Meanvalue)[2],"N","NA"));  
    write.csv(Data@meta.data, "M2score.txt", quote = FALSE, sep = ",",row.names = TRUE, col.names = TRUE)
'
awk -F"\t" 'ARGIND==1 && ($4>0 || $5>0 || $6>0){a[$1]=$1}\
            ARGIND==2{split($0,b,",");if(b[length(b)]=="Y"){if(b[1] in a){num1=num1+1}else{num2=num2+1}}}END{print num1"\t"num2}'\
            ../Figure3/GSE173958/M2results/umap.txt\
            "M2score.txt"

############################################
# 2. Clone presentation for MICs 
#    identified by unparied data
############################################

Rscript -e '
    MIC = read.csv("M2score.txt",sep=",",header=TRUE,check.names = F,row.names=1);
    Clones = read.table("../Figure3/GSE173958/GSM5283488_M2-PT-clones/clones.txt",sep="\t",header=TRUE,check.names = F);
    rownames(Clones) = Clones$Barcode;
    Clones = Clones[rownames(MIC),];
    MIC$Clones = Clones$cloneID;       
    tab = table(MIC[, c("MIC", "Clones")])
    print(c(tab[, "2"] / rowSums(tab),sum(tab[, "2"])/sum(tab)))
'
############################################
# 3. Biomarker presentation for MICs 
#    identified by unparied data
############################################
Rscript -e '
    library(Seurat);
    Data = readRDS("../Figure3/GSE173958/M2results/Primary/ScData.pca.rds");
    PrimaryMatrix=Data@assays$RNA@layers$data;
    colnames(PrimaryMatrix)=colnames(Data);
    rownames(PrimaryMatrix)=rownames(Data);
    genes=c("Epcam","Muc1","Cdh1", "Lgals4", "Ocln", "Ctse","Prrx1", "Ifitm1", "Ifitm3", "S100a4")
    MIC = read.csv("M2score.txt",sep=",",header=TRUE,check.names = F,row.names=1);
    MICindex = which(MIC$MIC == "Y")        
    MICcells=rownames(MIC[MICindex,]);
    nonMICcells=rownames(MIC[-MICindex,]);             
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
    write.table(SigGeneMatrix, "markerexpression.txt", quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE);              
' 
