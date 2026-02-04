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

conda env create -f 1.\ scMIC_Discussion.yml

cd FigureS/
mkdir scFoundation/
git clone https://github.com/biomap-research/scFoundation/ scFoundation/
#models in https://hopebio2020.sharepoint.com/:f:/s/PublicSharedfiles/IgBlEJ72TBE5Q76AmgXbgjXiAR69fzcrgzqgUYdSThPLrqk
pip install numpy pandas scanpy scipy einops local_attention torch torchvision torchaudio


############################################
# 0.1 Prepare expression data with 
#     human symbol for scFoundation
############################################

#Download mart_export.txt from ensembl biomart database
#Gene stable ID,Gene name,Mouse gene stable ID,Mouse gene name

Rscript -e '
  library(Seurat)
  Patients=c("M1results","M2results")
  Genes= read.table("mart_export.txt",sep="\t",header=TRUE,check.names = F);
  colnames(Genes)=c("ENSGHuman","GeneHuman","ENSGMouse","GeneMouse")
  Genes = subset(Genes,GeneHuman!="" & GeneMouse!="")
  rownames(Genes) = make.unique(Genes$GeneMouse)
  for(i in 1:length(Patients)){
    Data=readRDS(paste0("../Figure3/GSE173958/",Patients[i],"/Primary/ScData.pca.rds"))
    Matrix = Data@assays$RNA@layers$data
    rownames(Matrix) = rownames(Data)  
    colnames(Matrix) = colnames(Data)  
    genes=intersect(rownames(Matrix),rownames(Genes))
    Matrix = Matrix[genes,]
    rownames(Matrix) = make.unique(Genes[rownames(Matrix),]$GeneHuman)
    write.csv(t(Matrix), paste0("M",i,"data.csv"), quote = FALSE, sep = ",",row.names = TRUE, col.names = TRUE)
  }
' 

Rscript -e '
  library(Seurat)
  Patients=c("M1","M2")
  for(i in 1:length(Patients)){
    Data=readRDS(paste0("../Figure3/GSE173958/",Patients[i],"results/Primary/ScData.pca.rds"))
    OPcor1 = read.table(paste0("../Figure3/GSE173958/",Patients[i],"results/Metastasis_Met/proximity_OT.csv"),sep="\t",header=FALSE,check.names = F);
    OPcor2 = read.table(paste0("../Figure3/GSE173958/",Patients[i],"results/Metastasis_Liver/proximity_OT.csv"),sep="\t",header=FALSE,check.names = F);
    OPcor3 = read.table(paste0("../Figure3/GSE173958/",Patients[i],"results/Metastasis_Lung/proximity_OT.csv"),sep="\t",header=FALSE,check.names = F);
    output=data.frame(group = ifelse(OPcor1$V1>0 | OPcor2$V1>0 | OPcor3$V1>0,"Y","N"),score = OPcor1$V1+OPcor2$V1+OPcor3$V1)
    rownames(output) = colnames(Data)
    write.csv(output, paste0(Patients[i],"pheno.csv"), quote = FALSE, sep = ",",row.names = TRUE, col.names = TRUE)
  }
'

############################################
# 1. scFoundation embedding
############################################

cd scFoundation/model/
            
python get_embedding.py --task_name "M1" --input_type "singlecell" --output_type "cell" --pool_type "all" --data_path "../../M1data.csv" --save_path "../../" --tgthighres "a5" --pre_normalized "F" --version "rde" --ckpt_name "models" --model_path "models/"

python get_embedding.py --task_name "M2" --input_type "singlecell" --output_type "cell" --pool_type "all" --data_path "../../M2data.csv" --save_path "../../" --tgthighres "a5" --pre_normalized "F" --version "rde" --ckpt_name "models" --model_path "models/"
        
cd ../../
############################################
# 2. scFoundation embedding split
############################################

python << 'EOF'
import os
import numpy as np
import pandas as pd
tasks = ["M1", "M2"]
n_split = 4
for task_name in tasks:
    embed_file = f"{task_name}_models_singlecell_cell_embedding_a5_resolution.npy"
    count_file = f"{task_name}data.csv"
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

############################################
# 3. MIC identification based on 
#    scFoundation embedding 
############################################

rm -f scFoundation.result.txt
for((k=1;k<=200;k++))
do
  python ../scMIC_unpaired.py \
      --pri M2.x4.csv \
      --Refpri M1.x4.csv \
      --Refpripheno M1pheno.csv \
      --out proximity_OT.csv \
      --neighbork $k 
  
  awk -F"," 'ARGIND==1{a[$1]=$0}\
             ARGIND==2 && (FNR>1){print $0","a[$1]}'\
             proximity_OT.csv\
             M2pheno.csv |\
             awk -F"," -v k="$k" 'BEGIN{TP=0; FN=0; TN=0; FP=0}(NR>1){if($2=="Y"){if($(NF-2)=="Y"){TP=TP+1}else{FN=FN+1}}else{if($(NF-2)=="N"){FP=FP+1}else{TN=TN+1}}}\
             END{print "majorvote\t"k"\t"TP"\t"FN"\t"TN"\t"FP"\t"TP/(TP+FN)"\t"TN/(TN+FP)"\t"(TP+TN)/(TP+TN+FP+FN)}' >> scFoundation.result.txt
             
  awk -F"," 'ARGIND==1{a[$1]=$0}\
             ARGIND==2 && (FNR>1){print $0","a[$1]}'\
             proximity_OT.csv\
             M2pheno.csv |\
             awk -F"," -v k="$k" 'BEGIN{TP=0; FN=0; TN=0; FP=0}(NR>1){if($2=="Y"){if($(NF)=="Y"){TP=TP+1}else{FN=FN+1}}else{if($(NF)=="Y"){FP=FP+1}else{TN=TN+1}}}\
             END{print "meanvote\t"k"\t"TP"\t"FN"\t"TN"\t"FP"\t"TP/(TP+FN)"\t"TN/(TN+FP)"\t"(TP+TN)/(TP+TN+FP+FN)}' >> scFoundation.result.txt
done

############################################
# 4. Clone presentation for MICs 
#    identified by unparied data
############################################

python ../scMIC_unpaired.py \
    --pri M2.x4.csv \
    --Refpri M1.x4.csv \
    --Refpripheno M1pheno.csv \
    --out proximity_OT.csv \
    --neighbork 75

Rscript -e '
    MIC = read.csv("proximity_OT.csv",sep=",",header=TRUE,check.names = F,row.names=1);
    Clones = read.table("../Figure3/GSE173958/GSM5283488_M2-PT-clones/clones.txt",sep="\t",header=TRUE,check.names = F);
    rownames(Clones) = Clones$Barcode;
    Clones = Clones[rownames(MIC),];
    MIC$Clones = Clones$cloneID;       
    tab = table(MIC[, c("MIC_vote", "Clones")])
    print(c(tab[, "2"] / rowSums(tab),sum(tab[, "2"])/sum(tab)))
'
############################################
# 5. Biomarker presentation for MICs 
#    identified by unparied data
############################################
Rscript -e '
    library(Seurat);
    Data = readRDS("../Figure3/GSE173958/M2results/Primary/ScData.pca.rds");
    PrimaryMatrix=Data@assays$RNA@layers$data;
    colnames(PrimaryMatrix)=colnames(Data);
    rownames(PrimaryMatrix)=rownames(Data);
    genes=c("Epcam","Muc1","Cdh1", "Ocln", "Ctse","Prrx1", "Ifitm1", "S100a4")
    MIC = read.csv("proximity_OT.csv",sep=",",header=TRUE,check.names = F,row.names=1);
    MICindex = which(MIC$MIC_vote == "Y")        
    MICcells=rownames(MIC[MICindex,]);
    nonMICcells=rownames(MIC[-MICindex,]);             
    for(gene in genes){
      ttest = t.test(PrimaryMatrix[gene,MICcells],PrimaryMatrix[gene,nonMICcells]);
      print(c(gene,ttest$p.value,as.numeric(ttest$estimate)));
    }             
    SigGeneMatrix = data.frame(gene="",nonMIC=0,MIC=0);
    SigGeneMatrix = SigGeneMatrix[-1,];
    for(gene in genes){
      SigGeneMatrix = rbind(SigGeneMatrix,data.frame(gene=rep(gene,length(MICcells)),nonMIC=c(PrimaryMatrix[gene,nonMICcells],rep("",length(MICcells)-length(nonMICcells))),
                                       MIC=PrimaryMatrix[gene,MICcells]));
    }              
    write.table(SigGeneMatrix, "markerexpression.txt", quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE);              
' 
