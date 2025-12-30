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

############################################
# 1. FigureS1A
#    -Input:expression count during DataProcessing.sh
#    -Output:proximity_OT_test.csv: scMIC prediction
#            OTtraining.txt: training process
############################################

embeddings=("None" "None" "None" "AE" "AE" "AE" "PCA" "PCA" "PCA")
distances=("Pearson" "Spearman" "Euclidean" "Pearson" "Spearman" "Euclidean" "Pearson" "Spearman" "Euclidean")
topks=(1 1 1 1 1 1 1 1 1)
uotflags=("F" "F" "F" "F" "F" "F" "F" "F" "F")
MetFile="GSE249057/Results/Metastasis/MetastasisCount.csv"
PriFile="GSE249057/Results/Primary/PrimaryCount.csv"
outFile="GSE249057/Results/proximity_OT_test.csv"
for ((i=0; i<${#embeddings[@]}; i++))
do    
    if [[ "${uotflags[i]}" == "T" ]]; then
      python ../scMIC.py \
            --pri "$PriFile" \
            --met "$MetFile" \
            --out "$outFile" \
            --integration "${embeddings[i]}" \
            --distance "${distances[i]}" \
            --top_k "${topks[i]}" \
            --ot
    else
      python ../scMIC.py \
            --pri "$PriFile" \
            --met "$MetFile" \
            --out "$outFile" \
            --integration "${embeddings[i]}" \
            --distance "${distances[i]}" \
            --top_k "${topks[i]}"
    fi
    echo -e "${topks[i]}""\t""${embeddings[i]}""\t""${distances[i]}""\t""${uotflags[i]}" >> GSE249057/Results/OTtraining.txt
          
    Rscript -e 'library(Seurat); 
      library(dplyr);
      PrimaryData=readRDS("GSE249057/Results/Primary/ScData.pca.rds");  
      OPcor = read.table("GSE249057/Results/proximity_OT_test.csv",sep="\t",header=FALSE,check.names = F);
      OPresult=data.frame(cells=colnames(PrimaryData),value=OPcor$V1,group=PrimaryData@meta.data$seurat_clusters);
      OP_summary <- OPresult %>% group_by(group) %>% summarise(total_value = sum(value))
      OP_summary$fre = OP_summary$total_value/sum(OP_summary$total_value);
      print(as.data.frame(OP_summary))
    ' >> GSE249057/Results/OTtraining.txt
done

############################################
# 2. FigureS1B
#    -Input:expression count during DataProcessing.sh
#    -Output:proximity_OT_test.csv: scMIC prediction
#            OTtraining.txt: training process
############################################

embeddings=("None" "None" "None" "AE" "AE" "AE" "PCA" "PCA" "PCA")
distances=("Pearson" "Spearman" "Euclidean" "Pearson" "Spearman" "Euclidean" "Pearson" "Spearman" "Euclidean")
topks=(1 1 1 1 1 1 1 1 1)
uotflags=("T" "T" "T" "T" "T" "T" "T" "T" "T")
MetFile="GSE249057/Results/Metastasis/MetastasisCount.csv"
PriFile="GSE249057/Results/Primary/PrimaryCount.csv"
outFile="GSE249057/Results/proximity_OT_test.csv"
for ((i=0; i<${#embeddings[@]}; i++))
do    
    if [[ "${uotflags[i]}" == "T" ]]; then
      python ../scMIC.py \
            --pri "$PriFile" \
            --met "$MetFile" \
            --out "$outFile" \
            --integration "${embeddings[i]}" \
            --distance "${distances[i]}" \
            --top_k "${topks[i]}" \
            --ot
    else
      python ../scMIC.py \
            --pri "$PriFile" \
            --met "$MetFile" \
            --out "$outFile" \
            --integration "${embeddings[i]}" \
            --distance "${distances[i]}" \
            --top_k "${topks[i]}"
    fi
    echo -e "${topks[i]}""\t""${embeddings[i]}""\t""${distances[i]}""\t""${uotflags[i]}" >> GSE249057/Results/OTtraining.txt
          
    Rscript -e 'library(Seurat); 
      library(dplyr);
      PrimaryData=readRDS("GSE249057/Results/Primary/ScData.pca.rds");  
      OPcor = read.table("GSE249057/Results/proximity_OT_test.csv",sep="\t",header=FALSE,check.names = F);
      OPresult=data.frame(cells=colnames(PrimaryData),value=OPcor$V1,group=PrimaryData@meta.data$seurat_clusters);
      OP_summary <- OPresult %>% group_by(group) %>% summarise(total_value = sum(value))
      OP_summary$fre = OP_summary$total_value/sum(OP_summary$total_value);
      print(as.data.frame(OP_summary))
    ' >> GSE249057/Results/OTtraining.txt
done


############################################
# 3. FigureS1C
#    -Input:expression count during DataProcessing.sh
#    -Output:proximity_OT_test.csv: scMIC prediction
#            OTtraining.txt: training process
############################################

MetFile="GSE249057/Results/Metastasis/MetastasisCount.csv"
PriFile="GSE249057/Results/Primary/PrimaryCount.csv"
outFile="GSE249057/Results/proximity_OT_test.csv"
for ((i=1; i<=50; i++))
do    
    python ../scMIC.py \
        --pri "$PriFile" \
        --met "$MetFile" \
        --out "$outFile" \
        --integration "PCA" \
        --distance "Euclidean" \
        --top_k $i \
        --ot
    echo -e $i >> GSE249057/Results/OTtraining.txt
          
    Rscript -e 'library(Seurat); 
      library(dplyr);
      PrimaryData=readRDS("GSE249057/Results/Primary/ScData.pca.rds");  
      OPcor = read.table("GSE249057/Results/proximity_OT_test.csv",sep="\t",header=FALSE,check.names = F);
      OPresult=data.frame(cells=colnames(PrimaryData),value=OPcor$V1,group=PrimaryData@meta.data$seurat_clusters);
      OP_summary <- OPresult %>% group_by(group) %>% summarise(total_value = sum(value))
      OP_summary$fre = OP_summary$total_value/sum(OP_summary$total_value);
      print(as.data.frame(OP_summary))
    ' >> GSE249057/Results/OTtraining.txt
done


############################################
# 4. FigureS1D
#    -Input:expression count during DataProcessing.sh
#    -Output:proximity_OT_test.csv: scMIC prediction
#            OTtraining.txt: training process
############################################

embeddings=("None" "None" "None" "AE" "AE" "AE" "PCA" "PCA" "PCA")
distances=("Pearson" "Spearman" "Euclidean" "Pearson" "Spearman" "Euclidean" "Pearson" "Spearman" "Euclidean")
topks=(1 1 1 1 1 1 1 1 1)
uotflags=("T" "T" "T" "T" "T" "T" "T" "T" "T")
MetFile="GSE249057/Results/Metastasis/MetastasisCountHVG.csv"
PriFile="GSE249057/Results/Primary/PrimaryCountHVG.csv"
outFile="GSE249057/Results/proximity_OT_test.csv"
for ((i=0; i<${#embeddings[@]}; i++))
do    
    if [[ "${uotflags[i]}" == "T" ]]; then
      python ../scMIC.py \
            --pri "$PriFile" \
            --met "$MetFile" \
            --out "$outFile" \
            --integration "${embeddings[i]}" \
            --distance "${distances[i]}" \
            --top_k "${topks[i]}" \
            --ot
    else
      python ../scMIC.py \
            --pri "$PriFile" \
            --met "$MetFile" \
            --out "$outFile" \
            --integration "${embeddings[i]}" \
            --distance "${distances[i]}" \
            --top_k "${topks[i]}"
    fi
    echo -e "${topks[i]}""\t""${embeddings[i]}""\t""${distances[i]}""\t""${uotflags[i]}" >> GSE249057/Results/OTtraining.txt
          
    Rscript -e 'library(Seurat); 
      library(dplyr);
      PrimaryData=readRDS("GSE249057/Results/Primary/ScData.pca.rds");  
      OPcor = read.table("GSE249057/Results/proximity_OT_test.csv",sep="\t",header=FALSE,check.names = F);
      OPresult=data.frame(cells=colnames(PrimaryData),value=OPcor$V1,group=PrimaryData@meta.data$seurat_clusters);
      OP_summary <- OPresult %>% group_by(group) %>% summarise(total_value = sum(value))
      OP_summary$fre = OP_summary$total_value/sum(OP_summary$total_value);
      print(as.data.frame(OP_summary))
    ' >> GSE249057/Results/OTtraining.txt
done