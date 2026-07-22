#!/bin/bash
set -euo pipefail

############################################
# scMIC project
# Step 0. Download and organize scRNA-seq data 
#         and construct scMIC framework for 
#         MIC identification
# Dataset: OMIX002487
############################################

cd FigureS

############################################
# 1. FigureS13A
#    -Input:MIC score for OMIX002487 and LM size
#    -Output:Gene program correlation
############################################

Rscript -e '
    library(Seurat);
    Genes=read.table("../Figure5/MHMICgenes.txt",header=TRUE);
    colnames(Genes)[1:3]=c("Mgenes","Hgenes","Hit")
    Genes=subset(Genes,Hit==10)
    mousegenes=Genes$Mgenes
    humangenes=Genes$Hgenes #humangenes=subset(Genes,HumanOrNot=="Y")$Hgenes
    #OMIX002487
    pateints=c("P2","P3","P4","P5") 
    DiffScore = 0
    for(patient in pateints){
      Data=readRDS(paste0("../Figure4/OMIX002487/",patient,"/Primary/Tumour/ScData.pca.rds"));  
      OPcor = read.table(paste0("../Figure4/OMIX002487/",patient,"/proximity_OT.csv"),sep="\t",header=FALSE,check.names = F);
      Data$MIC = ifelse(OPcor$V1>0,"Y","N");
      PSelectedPrimaryData=subset(Data, features = intersect(humangenes,VariableFeatures(Data)));
      PSelectedPrimaryMatrix=PSelectedPrimaryData@assays$RNA@layers$data;
      PSelectedPrimaryavg <- colMeans(PSelectedPrimaryMatrix);
      index1 = which(Data@meta.data$MIC == "Y")  
      index2 = which(Data@meta.data$MIC == "N")
      DiffScore = c(DiffScore,mean(PSelectedPrimaryavg[index1])-mean(PSelectedPrimaryavg[index2]))
    }
    DiffScore = DiffScore[-1]
    print(DiffScore)
    Lsize1=c(1.19,4.9,2.36,6.6); #Liver metastasis maximum size
    Lsize2=c(1.07,2.44,2.09,5.2); #Liver metastasis min size
'

############################################
# 2. FigureS13B
#    -Input:tumor seurat data of OMIX002487
#    -Output:unpaired MIC prediction and gene program
############################################

# Prepare expression data with human symbol for scFoundation

Rscript -e '
  library(Seurat)
  Patients=c("P2","P3","P4","P5")
  genes=read.table("scFoundation/model/OS_scRNA_gene_index.19264.tsv",header=TRUE)
  for(i in 1:length(Patients)){
    Data=readRDS(paste0("../Figure4/OMIX002487/",Patients[i],"/Primary/Tumour/ScData.pca.rds"))
    Matrix = Data@assays$RNA@layers$data
    rownames(Matrix) = rownames(Data)  
    colnames(Matrix) = colnames(Data) 
    validgenes=intersect(rownames(Matrix),genes$gene_name)
    Matrix = Matrix[validgenes,]
    write.csv(t(Matrix), paste0(Patients[i],".data.csv"), quote = FALSE, sep = ",",row.names = TRUE, col.names = TRUE)
    OPcor = read.table(paste0("../Figure4/OMIX002487/",Patients[i],"/proximity_OT.csv"),sep="\t",header=FALSE,check.names = F);
    output=data.frame(group = ifelse(OPcor$V1>0,"Y","N"),score = OPcor$V1)
    rownames(output) = colnames(Data)
    write.csv(output, paste0(Patients[i],".pheno.csv"), quote = FALSE, sep = ",",row.names = TRUE, col.names = TRUE)
  }
'

# scFoundation embedding

cd scFoundation/model/
            
python get_embedding.py --task_name "P2" --input_type "singlecell" --output_type "cell" --pool_type "all" --data_path "../../P2.data.csv" --save_path "../../" --tgthighres "a5" --pre_normalized "F" --version "rde" --ckpt_name "models" --model_path "models/"

python get_embedding.py --task_name "P3" --input_type "singlecell" --output_type "cell" --pool_type "all" --data_path "../../P3.data.csv" --save_path "../../" --tgthighres "a5" --pre_normalized "F" --version "rde" --ckpt_name "models" --model_path "models/"

python get_embedding.py --task_name "P4" --input_type "singlecell" --output_type "cell" --pool_type "all" --data_path "../../P4.data.csv" --save_path "../../" --tgthighres "a5" --pre_normalized "F" --version "rde" --ckpt_name "models" --model_path "models/"

python get_embedding.py --task_name "P5" --input_type "singlecell" --output_type "cell" --pool_type "all" --data_path "../../P5.data.csv" --save_path "../../" --tgthighres "a5" --pre_normalized "F" --version "rde" --ckpt_name "models" --model_path "models/"
        
cd ../../

# scFoundation embedding split

python << 'EOF'
import os
import numpy as np
import pandas as pd
tasks = ["P2","P3", "P4", "P5"]
n_split = 4
for task_name in tasks:
    embed_file = f"{task_name}_models_singlecell_cell_embedding_a5_resolution.npy"
    count_file = f"{task_name}.data.csv"
    data = np.load(embed_file)                     # shape: (cells, dim)
    count = pd.read_csv(count_file, index_col=0)   # rows = cells
    splits = np.split(data, n_split, axis=1)
    for i, xi in enumerate(splits, start=1):
        df = pd.DataFrame(
            xi,
            index=count.index,
            columns=[f"X{j}" for j in range(1, xi.shape[1] + 1)]
        )
        out_file = os.path.join(f"{task_name}.x{i}.csv")
        df.to_csv(out_file)
EOF


# MIC identification based on scFoundation embedding 

python ../scMIC_unpaired.py \
    --pri P5.x4.csv \
    --Refpri P4.x4.csv \
    --Refpripheno P4.pheno.csv \
    --out proximity_OT.csv  
    
awk -F"," 'ARGIND==1{a[$1]=$0}\
           ARGIND==2 && (FNR>1){print $0","a[$1]}'\
           proximity_OT.csv\
           P5.pheno.csv > P5.mapped.predicted.csv
           
# Gene program identification based on this MIC identification result
           
Rscript -e '
    library(Seurat);
    Genes=read.table("../Figure5/MHMICgenes.txt",header=TRUE);
    colnames(Genes)[1:3]=c("Mgenes","Hgenes","Hit")
    Genes=subset(Genes,Hit==10)
    mousegenes=Genes$Mgenes
    humangenes=Genes$Hgenes #humangenes=subset(Genes,HumanOrNot=="Y")$Hgenes
    #OMIX002487
    patient="P5"
    Data=readRDS(paste0("../Figure4/OMIX002487/",patient,"/Primary/Tumour/ScData.pca.rds"));  
    OPcor = read.table("P5.mapped.predicted.csv",sep=",",header=FALSE,check.names = F);
    print(table(OPcor[,2]))
    print(table(OPcor[,c(5,7)]))
    allMICindex = which(OPcor$V5=="Y" & OPcor$V7=="Y")
    a = as.numeric(table(OPcor[allMICindex,2]))
    print(a[2]/(a[1]+a[2]))
    Data$MIC = "NA";
    Data$MIC[which(OPcor$V5=="Y" & OPcor$V7=="Y")] = "Y";
    Data$MIC[which(OPcor$V5=="N" & OPcor$V7=="N")] = "N";
    PSelectedPrimaryData=subset(Data, features = intersect(humangenes,VariableFeatures(Data)));
    PSelectedPrimaryMatrix=PSelectedPrimaryData@assays$RNA@layers$data;
    PSelectedPrimaryavg <- colMeans(PSelectedPrimaryMatrix);
    index1 = which(Data@meta.data$MIC == "Y")  
    index2 = which(Data@meta.data$MIC == "N")
    ttest=t.test(PSelectedPrimaryavg[index1],PSelectedPrimaryavg[index2]);
    wilcoxtest=wilcox.test(PSelectedPrimaryavg[index1],PSelectedPrimaryavg[index2])
    meanvalue = as.numeric(ttest$estimate)
    print(c(patient,length(index1),length(index2),ttest$p.value,wilcoxtest$p.value,meanvalue,meanvalue[1]-meanvalue[2]))
'


