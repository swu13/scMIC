#!/bin/bash
set -euo pipefail

############################################
# scMIC project
# Step 0. Download and organize scRNA-seq data 
#         and construct scMIC framework for 
#         MIC identification
# Dataset: GSE173958, GSE197177, GSE277783 and TCGA
############################################


# Project root directory
RootPath="/data/twang15/wsj/scMIC"
cd "${RootPath}"
mkdir Figure5
cd Figure5/

conda env create -f 1.\ scMIC_Figure5.yml
############################################
# 1. Figure5B
#    -Input:MIC latent identification
#    -Output:training process
############################################

Rscript -e '
    library(Seurat);
    PrimaryData=readRDS(paste0("../Figure3/GSE173958/M1results/Primary/ScData.pca.rds"));  
    PrimaryMatrix=PrimaryData@assays$RNA@layers$data;
    colnames(PrimaryMatrix)=colnames(PrimaryData);
    rownames(PrimaryMatrix)=rownames(PrimaryData);     
    variable_genes <- VariableFeatures(PrimaryData);         
    selecteddatamatrix = PrimaryMatrix[variable_genes,];  
    write.csv(t(selecteddatamatrix), file = paste0("scRNA.csv"), row.names = TRUE, quote = FALSE);   
    OPcor = read.table("../Figure3/GSE173958/M1results/Metastasis_Met/proximity_OT.csv",sep="\t",header=FALSE,check.names = F);
    phenotypes = data.frame(Peritoneal = ifelse(OPcor$V1>0,"MIC","non-MIC"))
    rownames(phenotypes) = rownames(PrimaryData@meta.data)
    OPcor = read.table("../Figure3/GSE173958/M1results/Metastasis_Liver/proximity_OT.csv",sep="\t",header=FALSE,check.names = F);
    phenotypes$Liver = ifelse(OPcor$V1>0,"MIC","non-MIC")
    OPcor = read.table("../Figure3/GSE173958/M1results/Metastasis_Lung/proximity_OT.csv",sep="\t",header=FALSE,check.names = F);
    phenotypes$Group = ifelse(phenotypes$Peritoneal == "MIC" | phenotypes$Liver == "MIC" | OPcor$V1>0, "MIC", "non-MIC")       
    write.csv(phenotypes, file = paste0("phenotype.csv"), row.names = TRUE, quote = FALSE);
'

#Multi-time training

for((i=1;i<=10;i++))
do
seed=$RANDOM
mkdir ./GeneProgram$i/
echo $seed > ./GeneProgram$i/seed.txt
python ../GeneProgram.py\
   --scRNA_file scRNA.csv\
   --phenotype_file phenotype.csv\
   --output_path ./GeneProgram$i/\
   --epochs 500\
   --Flag "Training"\
   --lambda_mic 1 \
   --seed $seed &
done

############################################
# 2. Figure5C
#    -Input:MIC latent identification
#    -Output:MIC latent
############################################
#get MIC latent-specific features
for((i=1;i<=10;i++))
do
python ../GeneProgram.py\
    --scRNA_file scRNA.csv\
    --phenotype_file phenotype.csv\
    --model_path ./GeneProgram$i/\
    --output_path ./GeneProgram$i/\
    --Flag "Inferencing" &
done

#DEG
Rscript -e '
    library(Seurat);
    PrimaryData=readRDS(paste0("../Figure3/GSE173958/M1results/Primary/ScData.pca.rds"));    
    OPcor1 = read.table("../Figure3/GSE173958/M1results/Metastasis_Met/proximity_OT.csv",sep="\t",header=FALSE,check.names = F);
    OPcor2 = read.table("../Figure3/GSE173958/M1results/Metastasis_Liver/proximity_OT.csv",sep="\t",header=FALSE,check.names = F);
    OPcor3 = read.table("../Figure3/GSE173958/M1results/Metastasis_Lung/proximity_OT.csv",sep="\t",header=FALSE,check.names = F);
    PrimaryData$MIC = ifelse(OPcor1$V1 > 0 | OPcor2$V1 > "MIC" | OPcor3$V1>0, "MIC", "non-MIC") 
    Idents(PrimaryData) <- PrimaryData@meta.data$MIC;
    DEG <- FindMarkers(PrimaryData, ident.1 = "MIC",ident.2 = "non-MIC", group.by="MIC",test.use = "wilcox", logfc.threshold=0);
    write.table(DEG, "DEG_MIC.txt", quote = FALSE, sep = "\t", row.names = TRUE, col.names = TRUE); 
'

#DEG and MIC-latent features overlap
Rscript -e '  
    get_top_genes_vector <- function(mat, top_k = 3) {
      top_genes_list <- apply(mat, 2, function(col) {
        ord <- order(col, decreasing = TRUE)
        idx <- ord[1:top_k]
        data.frame(
          gene = rownames(mat)[idx],
          rank = 1:top_k
        )
      })
      df <- do.call(rbind, top_genes_list)
      best_rank <- aggregate(rank ~ gene, df, min)
      colnames(best_rank) <- c("gene", "min_rank")          
      return(best_rank[order(best_rank$min_rank), ])
    } 
    OverlapGenes <- function(mergeDEG, ICfeatures) {
      genes=intersect(rownames(mergeDEG),ICfeatures$gene);
      OverlapGenes=subset(ICfeatures,gene %in% genes);
      rownames(OverlapGenes)=OverlapGenes$gene;     
      OverlapGenes=cbind(OverlapGenes,mergeDEG[rownames(OverlapGenes),])
    }
    DEG = read.table("DEG_MIC.txt")
    num=c(1,2,3,4,5,6,7,8,9,10)
    MIClatents = list(c(5,6,30),c(14,21),c(5,25,26),c(15,30),c(7,15,30),c(17,20),c(15,21),c(0,13),c(6,26),c(7,25))
    #MIClatents = list(c(11,28),c(5,9,12,18),c(1,9,20,25),c(26,27),c(2,9,24,27),c(6,17,22),c(9,10,24),c(6,14,16),c(2,31),c(9,23,26))
    Output=list()
    for(i in 1:length(num)){
      F_all=read.csv(paste0("GeneProgram",num[i],"/F_all.csv"))
      F_all = as.matrix(t(F_all))
      scRNA=read.csv("scRNA.csv",row.names=1)
      rownames(F_all)=colnames(scRNA) 
      #sd_values = t(apply(F_all, 1, function(x) {sd(x)}))    
      #index=which(sd_values>0.05)
      #F_info = F_all[index,]            
      MIClatent = MIClatents[[i]]
      MICfeatures=get_top_genes_vector(as.matrix(F_all[,(unlist(MIClatent)+1)]),500)
      a = OverlapGenes(DEG,MICfeatures)
      b = subset(a,p_val_adj<0.05 & avg_log2FC>0.1 & pct.1>0.5)
      Output[[i]] = b$gene
    }
    gene_rank_df <- do.call(rbind, lapply(seq_along(Output), function(i) {
      data.frame(
        gene = Output[[i]],
        rank = seq_along(Output[[i]]),
        iter = i,
        stringsAsFactors = FALSE
      )
    }))
    freq_df <- aggregate(iter ~ gene, gene_rank_df, length)
    colnames(freq_df)[2] <- "freq"
    best_rank_df <- aggregate(rank ~ gene, gene_rank_df, min)
    colnames(best_rank_df)[2] <- "best_rank"
    GeneSummary <- merge(freq_df, best_rank_df, by = "gene")
    GeneSummary <- GeneSummary[
      order(-GeneSummary$freq, GeneSummary$best_rank),
    ]
    write.table(GeneSummary, "MICgenes.txt", quote = FALSE, sep = "\t", row.names = TRUE, col.names = TRUE); 
'

#HMgenes.txt: mouse gene\tHuman genes
awk -F"\t" 'ARGIND==1{a[$4]=a[$4]";"$2}ARGIND==2{if(FNR==1){print "Mgenes\tHgenes\tHit"}\
                                                 else{if($2 in a){split(a[$2],b,";");for(i=2;i<=length(b);i++){print $2"\t"b[i]"\t"$3}}\
                                                      else{print $2"\tNA\t"$3}}}' mart_export.txt MICgenes.txt > MHMICgenes.txt
#check MHMICgenes.txt manually

############################################
# 3. Figure5D
#    -MIC gene enrichment
############################################
#WebGestaltR
Rscript -e '
    library(WebGestaltR);
    Genes=read.table("MHMICgenes.txt",header=TRUE,fill = TRUE);
    colnames(Genes)[1:3]=c("Mgenes","Hgenes","Hit")
    Genes=subset(Genes,Hit==10)
    score=data.frame(gene=Genes$Mgenes,score=1);
    res <- WebGestaltR(enrichMethod = "GSEA",organism = "mmusculus",interestGene = score,interestGeneType = "genesymbol",referenceSet = "genome",
           enrichDatabase = "geneontology_Biological_Process",fdrMethod = "BH",sigMethod = "fdr",fdrThr = 0.05,minNum = 3);
'
#Genes for DAVID

############################################
# 4. Figure5E
#    -Human MIC genes overlap and MIC comparison
############################################

#Rscript -e '
#    library(Seurat);
#    PrimaryData=readRDS(paste0("../Figure4/OMIX002487/P3/Primary/Tumour/ScData.pca.rds"));    
#    OPcor = read.table("../Figure4/OMIX002487/P3/proximity_OT.csv",sep="\t",header=FALSE,check.names = F);
#    PrimaryData$MIC = ifelse(OPcor$V1 > 0, "MIC", "non-MIC") 
#    Idents(PrimaryData) <- PrimaryData@meta.data$MIC;
#    DEG <- FindMarkers(PrimaryData, ident.1 = "MIC",ident.2 = "non-MIC", group.by="MIC",test.use = "wilcox", logfc.threshold=0);
#    write.table(DEG, "DEG_MIC_humanP3.txt", quote = FALSE, sep = "\t", row.names = TRUE, col.names = TRUE); 
#'
#
##MHMICgenes.txt: MICgenes.txt mapped human genes version: Mgenes, Hgenes, Hit
#
#awk -F"\t" 'ARGIND==1 && ($3>0){a[$1]=$1}ARGIND==2{if(FNR==1){print $0"\tHumanOrNot"}else if($2 in a){print $0"\tY"}else{print $0"\tN"}}' DEG_MIC_humanP3.txt MHMICgenes.txt > MHMICgenes1.txt
#mv MHMICgenes1.txt MHMICgenes.txt

############################################
# 3. Figure5E
#    -MIC comparison
############################################

Rscript -e '
    library(Seurat);
    Genes=read.table("MHMICgenes.txt",header=TRUE);
    colnames(Genes)[1:3]=c("Mgenes","Hgenes","Hit")
    Genes=subset(Genes,Hit==10)
    mousegenes=Genes$Mgenes
    humangenes=Genes$Hgenes #humangenes=subset(Genes,HumanOrNot=="Y")$Hgenes
    #GSE173958
    pateints=c("M1","M2")
    output=data.frame(sample="",MIC=0,nonMIC=0)
    for(patient in pateints){
      Data = readRDS(paste0("../Figure3/GSE173958/",patient,"results/Primary/ScData.pca.rds"));
      MIC = read.table(paste0("../Figure3/GSE173958/",patient,"results/umap.txt"),sep="\t",header=TRUE,check.names = F);
      Data$MIC = ifelse(MIC$Peritoneal>5 | MIC$Liver>5 | MIC$Lung>5,"Y","N")
      PSelectedPrimaryData=subset(Data, features = intersect(mousegenes,VariableFeatures(Data)));
      PSelectedPrimaryMatrix=PSelectedPrimaryData@assays$RNA@layers$data;
      PSelectedPrimaryavg <- colMeans(PSelectedPrimaryMatrix);  
      index = which(Data@meta.data$MIC == "Y")            
      print(t.test(PSelectedPrimaryavg[index],PSelectedPrimaryavg[-index]));
      maxnum=max(length(index),ncol(Data)-length(index))
      output=rbind(output,data.frame(sample=rep(patient,maxnum),MIC=c(PSelectedPrimaryavg[index],rep("",maxnum-length(index))),
                   nonMIC=c(PSelectedPrimaryavg[-index],rep("",maxnum-(ncol(Data)-length(index))))))
    }
    #OMIX002487
    pateints=c("P2","P3","P4"); #,"P5" 
    for(patient in pateints){
      Data=readRDS(paste0("../Figure4/OMIX002487/",patient,"/Primary/Tumour/ScData.pca.rds"));  
      OPcor = read.table(paste0("../Figure4/OMIX002487/",patient,"/proximity_OT.csv"),sep="\t",header=FALSE,check.names = F);
      Data$MIC = ifelse(OPcor$V1>0,"Y","N");
      PSelectedPrimaryData=subset(Data, features = intersect(humangenes,VariableFeatures(Data)));
      PSelectedPrimaryMatrix=PSelectedPrimaryData@assays$RNA@layers$data;
      PSelectedPrimaryavg <- colMeans(PSelectedPrimaryMatrix);
      index1 = which(Data@meta.data$MIC == "Y")  
      index2 = which(Data@meta.data$MIC == "N")
      if(length(index1)>3){
        ttest=t.test(PSelectedPrimaryavg[index1],PSelectedPrimaryavg[index2]);
        wilcoxtest=wilcox.test(PSelectedPrimaryavg[index1],PSelectedPrimaryavg[index2])
        print(c(patient,length(index1),length(index2),ttest$p.value,wilcoxtest$p.value,as.numeric(ttest$estimate)))
        maxnum=max(length(index1),length(index2))
        output=rbind(output,data.frame(sample=rep(patient,maxnum),MIC=c(PSelectedPrimaryavg[index1],rep("",maxnum-length(index1))),
                      nonMIC=c(PSelectedPrimaryavg[index2],rep("",maxnum-length(index2)))))
      }
    }
    #GSE277783
    pateints=c("Pt-1","Pt-2","Pt-3","Pt-4","Pt-5","Pt-6","Pt-7","Pt-8","Pt-9","Pt-10","Pt-11","Pt-13")
    for(patient in pateints){
      Data=readRDS(paste0("../Figure4/GSE277783/",patient,"/Primary/Tumour/ScData.pca.rds")); 
      Data$MIC = "N"
      OPfolder <- list.files(paste0("../Figure4/GSE277783/",patient,"/"),pattern = paste0("Metastasis", "[B-E]$"),full.names = TRUE)
      for(i in 1:length(OPfolder)){
        OPfile=paste0(OPfolder[i],"/proximity_OT.csv")
        OPcor = read.table(OPfile,sep="\t",header=FALSE,check.names = F);
        if(i==1){
          sum_OT=OPcor$V1
        }else{
          sum_OT=sum_OT+OPcor$V1
        }
      }
      Data$MIC  = ifelse(sum_OT>5,"Y","N");
      PSelectedPrimaryData=subset(Data, features = intersect(humangenes,VariableFeatures(Data)));
      PSelectedPrimaryMatrix=PSelectedPrimaryData@assays$RNA@layers$data;
      PSelectedPrimaryavg <- colMeans(PSelectedPrimaryMatrix);
      index = which(Data@meta.data$MIC == "Y") 
      if(length(index)>=100){ 
        ttest=t.test(PSelectedPrimaryavg[index],PSelectedPrimaryavg[-index]);
        wilcoxtest=wilcox.test(PSelectedPrimaryavg[index],PSelectedPrimaryavg[-index]);
        print(c(patient,length(index),ncol(Data)-length(index),ttest$p.value,wilcoxtest$p.value,as.numeric(ttest$estimate)))
        maxnum=max(length(index),ncol(Data)-length(index))
        output=rbind(output,data.frame(sample=rep(patient,maxnum),MIC=c(PSelectedPrimaryavg[index],rep("",maxnum-length(index))),
                   nonMIC=c(PSelectedPrimaryavg[-index],rep("",maxnum-(ncol(Data)-length(index))))))
      }
    }
    write.table(output[-1,], "MICscore.txt", sep = "\t", quote = FALSE, row.names = TRUE);
'

############################################
# 3. Figure5E
#    -TCGA MIC program comparisons
############################################
mkdir Bulk/
wget -c https://gdc-hub.s3.us-east-1.amazonaws.com/download/TCGA-PAAD.star_fpkm.tsv.gz -O Bulk/TCGA-PAAD.star_fpkm.tsv.gz
wget -c https://gdc-hub.s3.us-east-1.amazonaws.com/download/TCGA-PAAD.clinical.tsv.gz -O Bulk/TCGA-PAAD.clinical.tsv.gz
wget -c https://gdc-hub.s3.us-east-1.amazonaws.com/download/TCGA-PAAD.survival.tsv.gz -O Bulk/TCGA-PAAD.survival.tsv.gz
wget -c https://gdc-hub.s3.us-east-1.amazonaws.com/download/gencode.v36.annotation.gtf.gene.probemap -P Bulk/
gunzip Bulk/*.gz

Rscript -e '
    Genes=read.table("MHMICgenes.txt",header=TRUE);
    colnames(Genes)[1:3]=c("Mgenes","Hgenes","Hit")
    Genes=subset(Genes,Hit==10)
    humangenes=Genes$Hgenes #humangenes=subset(Genes,HumanOrNot=="Y")$Hgenes
    humangenes=Genes$Hgenes
    PAAD_fpkm_data <- read.table("Bulk/TCGA-PAAD.star_fpkm.tsv", header=T, sep="\t",check.names=FALSE);
    gene_annotation <- read.table("Bulk/gencode.v36.annotation.gtf.gene.probemap", sep = "\t", header = TRUE);
    merged_data <- merge(gene_annotation,PAAD_fpkm_data, by.x = "id", by.y = "Ensembl_ID");
    rownames(merged_data) <- make.unique(merged_data$gene);   
    TCGAnames=colnames(merged_data)[7:dim(merged_data)[2]]; 
    disease_sample_cols <- TCGAnames[grepl("-01A$", TCGAnames)]; 
    tumour_data <- merged_data[, disease_sample_cols];  
    selected_tumour_data =  tumour_data[intersect(rownames(tumour_data),humangenes),]  
    SelectedPrimaryavg <- colMeans(selected_tumour_data);
    library(survival);
    survival_data <- read.table("Bulk/TCGA-PAAD.survival.tsv", header = TRUE,sep = "\t",quote = "",comment.char = "",fill = TRUE);
    rownames(survival_data)=survival_data$sample;                                      
    data <- cbind(SelectedPrimaryavg,survival_data[colnames(selected_tumour_data),c("OS.time","OS")]);

    cox_model <- coxph(Surv(OS.time,OS)~get(colnames(data)[1]),data=data);
    beta=coef(cox_model);
    se=sqrt(diag(vcov(cox_model)));
    print(c(1 - pchisq((beta/se)^2, 1),exp(beta)));
    
    median_value <- median(data[,1], na.rm = TRUE);
    data$group <- ifelse(data[,1] > median_value, "High", "Low");                    
    surv_fit <- survfit(Surv(OS.time, OS) ~ group, data = data);
    surv_diff <- survdiff(Surv(OS.time, OS) ~ group, data = data);
    p_value <- 1 - pchisq(surv_diff$chisq, df = 1);
    print(c(p_value));
       
    data$group <- factor(data$group);
    surv_data <- data.frame(time = surv_fit$time,n.risk = surv_fit$n.risk,n.event = surv_fit$n.event,n.censor = surv_fit$n.censor,surv_prob = surv_fit$surv,
      std_err = surv_fit$std.err,upper = surv_fit$upper,lower = surv_fit$lower,group = rep(levels(data$group), surv_fit$strata));
    write.csv(surv_data, paste0("Bulk/MICscore_survival_data.txt"), row.names = FALSE);
'
