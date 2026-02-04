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
cd Figure6/

conda env create -f 1.\ scMIC_Figure6.yml

############################################
# 1. Figure6A
#    -marker comparison
############################################

Rscript -e '
    library(Seurat);
    #GSE173958
    pateints=c("M1","M2")
    output=data.frame(sample="",MIC=0,nonMIC=0)
    gene = "Ociad2"
    for(patient in pateints){
      Data = readRDS(paste0("../Figure3/GSE173958/",patient,"results/Primary/ScData.pca.rds"));
      MIC = read.table(paste0("../Figure3/GSE173958/",patient,"results/umap.txt"),sep="\t",header=TRUE,check.names = F);
      Data$MIC = ifelse(MIC$Peritoneal>5 | MIC$Liver>5 | MIC$Lung>5,"Y","N")
      PrimaryMatrix = Data@assays$RNA@layers$data;
      colnames(PrimaryMatrix) = colnames(Data);
      rownames(PrimaryMatrix) = rownames(Data); 
      index = which(Data@meta.data$MIC == "Y")            
      print(t.test(as.numeric(PrimaryMatrix[gene,index]),as.numeric(PrimaryMatrix[gene,-index])));
      maxnum=max(length(index),ncol(Data)-length(index))
      output=rbind(output,data.frame(sample=rep(patient,maxnum),MIC=c(PrimaryMatrix[gene,index],rep("",maxnum-length(index))),
                   nonMIC=c(PrimaryMatrix[gene,-index],rep("",maxnum-(ncol(Data)-length(index))))))
    }
    #OMIX002487
    pateints=c("P2","P3","P4"); #,"P5" 
    gene = "OCIAD2"
    for(patient in pateints){
      Data=readRDS(paste0("../Figure4/OMIX002487/",patient,"/Primary/Tumour/ScData.pca.rds"));  
      OPcor = read.table(paste0("../Figure4/OMIX002487/",patient,"/proximity_OT.csv"),sep="\t",header=FALSE,check.names = F);
      Data$MIC = ifelse(OPcor$V1>0,"Y","N");
      PrimaryMatrix = Data@assays$RNA@layers$data;
      colnames(PrimaryMatrix) = colnames(Data);
      rownames(PrimaryMatrix) = rownames(Data); 
      index = which(Data@meta.data$MIC == "Y") 
      if(length(index)>3){ 
        ttest=t.test(as.numeric(PrimaryMatrix[gene,index]),as.numeric(PrimaryMatrix[gene,-index]));
        print(c(patient,length(index),ncol(Data)-length(index),ttest$p.value,as.numeric(ttest$estimate)))
        maxnum=max(length(index),ncol(Data)-length(index))
        output=rbind(output,data.frame(sample=rep(patient,maxnum),MIC=c(PrimaryMatrix[gene,index],rep("",maxnum-length(index))),
                   nonMIC=c(PrimaryMatrix[gene,-index],rep("",maxnum-(ncol(Data)-length(index))))))
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
      PrimaryMatrix = Data@assays$RNA@layers$data;
      colnames(PrimaryMatrix) = colnames(Data);
      rownames(PrimaryMatrix) = rownames(Data); 
      index = which(Data@meta.data$MIC == "Y")  
      if(length(index)>=100){
        ttest=t.test(as.numeric(PrimaryMatrix[gene,index]),as.numeric(PrimaryMatrix[gene,-index]));
        wilcoxtest=wilcox.test(as.numeric(PrimaryMatrix[gene,index]),as.numeric(PrimaryMatrix[gene,-index]));
        print(c(patient,length(index),ncol(Data)-length(index),ttest$p.value,wilcoxtest$p.value,as.numeric(ttest$estimate)))
        maxnum=max(length(index),ncol(Data)-length(index))
        output=rbind(output,data.frame(sample=rep(patient,maxnum),MIC=c(PrimaryMatrix[gene,index],rep("",maxnum-length(index))),
                   nonMIC=c(PrimaryMatrix[gene,-index],rep("",maxnum-(ncol(Data)-length(index))))))
      }
    }
    write.table(output[-1,], "MICscore.txt", sep = "\t", quote = FALSE, row.names = TRUE);
'

############################################
# 2. Figure6B
#    -TCGA marker tumor normal comparison
############################################

mkdir Bulk/
wget -c https://toil-xena-hub.s3.us-east-1.amazonaws.com/download/gtex_RSEM_gene_fpkm.gz -O Bulk/gtex_RSEM_gene_fpkm.gz
wget -c https://toil-xena-hub.s3.us-east-1.amazonaws.com/download/GTEX_phenotype.gz -O Bulk/GTEX_phenotype.gz
gunzip Bulk/*.gz

Rscript -e '
    PAAD_fpkm_data <- read.table("../Figure5/Bulk/TCGA-PAAD.star_fpkm.tsv", header=T, sep="\t",check.names=FALSE);
    gene_annotation <- read.table("../Figure5/Bulk/gencode.v36.annotation.gtf.gene.probemap", sep = "\t", header = TRUE);       
    merged_data <- merge(gene_annotation,PAAD_fpkm_data, by.x = "id", by.y = "Ensembl_ID");
    rownames(merged_data) <- make.unique(merged_data$gene);   
    TCGAnames=colnames(merged_data)[7:dim(merged_data)[2]]; 
    disease_sample_cols <- TCGAnames[grepl("-01A$", TCGAnames)]; 
    tumour_data <- merged_data[, disease_sample_cols]; 
    normal_sample_cols <- TCGAnames[grepl("-11A$", TCGAnames)];
    normal_data <- merged_data[, normal_sample_cols];
    gtex_data <- read.delim("Bulk/gtex_RSEM_gene_fpkm",header = TRUE,sep = "\t",stringsAsFactors = FALSE,check.names = FALSE);
    rownames(gtex_data) <- gtex_data$sample;
    gtex_phenotype_data <- read.delim("Bulk/GTEX_phenotype", header = TRUE,sep = "\t",stringsAsFactors = FALSE,check.names = FALSE);            
    pancreas_samples <- gtex_phenotype_data[gtex_phenotype_data[["body_site_detail (SMTSD)"]] == "Pancreas", ];
    PRADnormal=gtex_data[,intersect(colnames(gtex_data),pancreas_samples$Sample)];
    gene="OCIAD2"                                      
    ttest=t.test(as.numeric(tumour_data[gene,]),as.numeric(normal_data[gene,]));
    print(c(gene,"TCGA",ttest$p.value,as.numeric(ttest$estimate),dim(tumour_data)[2],dim(normal_data)[2]));
    output=data.frame(group=c(rep("Tumor",dim(tumour_data)[2]),rep("Normal",dim(normal_data)[2])),value=c(as.numeric(tumour_data[gene,]),as.numeric(normal_data[gene,])));
    write.csv(output, paste0("Bulk/",gene,"_TCGA_TN.txt"), row.names = FALSE);
    ENSGID=gene_annotation$id[which(gene_annotation$gene==gene)];
    GTEXnormal=as.numeric(PRADnormal[which(sapply(strsplit(rownames(PRADnormal), ".", fixed = TRUE), `[`, 1)==sapply(strsplit(ENSGID, ".", fixed = TRUE), `[`, 1)),])
    GTEXnormal_changed=log2(2^GTEXnormal - 0.001+1);
    ttest=t.test(as.numeric(tumour_data[gene,]),GTEXnormal_changed);
    print(c(gene,"GTEx",ttest$p.value,as.numeric(ttest$estimate),dim(tumour_data)[2],length(GTEXnormal_changed)));
    output=data.frame(group=c(rep("Tumor",dim(tumour_data)[2]),rep("Normal",length(GTEXnormal_changed))),value=c(as.numeric(tumour_data[gene,]),GTEXnormal_changed));
    write.csv(output, paste0("Bulk/",gene,"_TCGA_GTEx.txt"), row.names = FALSE);
'

############################################
# 3. Figure6C
#    -TCGA marker survival
############################################

Rscript -e '
    library(Seurat); 
    PAAD_fpkm_data <- read.table("../Figure5/Bulk/TCGA-PAAD.star_fpkm.tsv", header=T, sep="\t",check.names=FALSE);
    gene_annotation <- read.table("../Figure5/Bulk/gencode.v36.annotation.gtf.gene.probemap", sep = "\t", header = TRUE);
    merged_data <- merge(gene_annotation,PAAD_fpkm_data, by.x = "id", by.y = "Ensembl_ID");
    rownames(merged_data) <- make.unique(merged_data$gene);   
    TCGAnames=colnames(merged_data)[7:dim(merged_data)[2]]; 
    disease_sample_cols <- TCGAnames[grepl("-01A$", TCGAnames)]; 
    tumour_data <- merged_data[, disease_sample_cols];  
    library(survival);
    survival_data <- read.table("../Figure5/Bulk/TCGA-PAAD.survival.tsv", header = TRUE,sep = "\t",quote = "",comment.char = "",fill = TRUE);
    rownames(survival_data)=survival_data$sample;   
    gene="OCIAD2" 
    data <- t(tumour_data[gene,]);
    data <- cbind(data,survival_data[rownames(data),c("OS.time","OS")]);

    cox_model <- coxph(Surv(OS.time,OS)~get(colnames(data)[1]),data=data);
    beta=coef(cox_model);
    se=sqrt(diag(vcov(cox_model)));
    print(c(1 - pchisq((beta/se)^2, 1),exp(beta)));
    
    quantiles <- quantile(data[,1],probs=c(0.25,0.75));
    data$group <- ifelse(data[,1] > as.numeric(quantiles)[2], "High", ifelse(data[,1] < as.numeric(quantiles)[1],"Low","NA"));
    data <- data[which(data$group!="NA"),]                    
    surv_fit <- survfit(Surv(OS.time, OS) ~ group, data = data);
    surv_diff <- survdiff(Surv(OS.time, OS) ~ group, data = data);
    p_value <- 1 - pchisq(surv_diff$chisq, df = 1);
    print(c(p_value));
       
    data$group <- factor(data$group);
    surv_data <- data.frame(time = surv_fit$time,n.risk = surv_fit$n.risk,n.event = surv_fit$n.event,n.censor = surv_fit$n.censor,surv_prob = surv_fit$surv,
      std_err = surv_fit$std.err,upper = surv_fit$upper,lower = surv_fit$lower,group = rep(levels(data$group), surv_fit$strata));
    write.csv(surv_data, paste0("Bulk/MICscore_survival_data.txt"), row.names = FALSE);
'
############################################
# 4. Figure6D
#    -TCGA marker stemness association
############################################

#Bulk/PAAD_mRNAsi.csv downloaded from http://bio-bigdata.hrbmu.edu.cn/CancerStemnessOnline/data.jsp

Rscript -e '
    PAAD_fpkm_data <- read.table("../Figure5/Bulk/TCGA-PAAD.star_fpkm.tsv", header=T, sep="\t",check.names=FALSE);
    gene_annotation <- read.table("../Figure5/Bulk/gencode.v36.annotation.gtf.gene.probemap", sep = "\t", header = TRUE);
    merged_data <- merge(gene_annotation,PAAD_fpkm_data, by.x = "id", by.y = "Ensembl_ID");
    rownames(merged_data) <- make.unique(merged_data$gene);   
    TCGAnames=colnames(merged_data)[7:dim(merged_data)[2]]; 
    disease_sample_cols <- TCGAnames[grepl("-01A$", TCGAnames)]; 
    tumour_data <- merged_data[, disease_sample_cols]; 
    PAAD_mRNAsi_data <- read.table("PAAD_mRNAsi.txt",row.names = 1,header=TRUE)
    samples=intersect(rownames(PAAD_mRNAsi_data),colnames(tumour_data))
    tumour_data = tumour_data[,samples]
    stemness <- PAAD_mRNAsi_data[samples, , drop = FALSE]
    stemness$gene = as.numeric(tumour_data["OCIAD2",])
    print(cor.test(stemness$gene, stemness$CSscore, method = "pearson"))
    print(cor.test(stemness$gene, stemness$CSscore, method = "spearman"))
    write.csv(stemness, paste0("stemness_data.txt"), row.names = FALSE);
'
