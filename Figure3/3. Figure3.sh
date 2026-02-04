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
# 1. Figure3A
#    -Input:Normalized ScData.pca.rds during DataProcessing.sh
#    -Output:umap with MIC identification
############################################

Rscript -e '
    library(Seurat); 
    Data = readRDS("GSE173958/M1results/Primary/ScData.pca.rds");
    umap_embeddings <- Embeddings(Data, "umap");
    umap_embeddings = as.data.frame(umap_embeddings);    
    OPcor = read.table("GSE173958/M1results/Metastasis_Met/proximity_OT.csv",sep="\t",header=FALSE,check.names = F);
    umap_embeddings$Peritoneal = OPcor$V1
    OPcor = read.table("GSE173958/M1results/Metastasis_Liver/proximity_OT.csv",sep="\t",header=FALSE,check.names = F);
    umap_embeddings$Liver = OPcor$V1
    OPcor = read.table("GSE173958/M1results/Metastasis_Lung/proximity_OT.csv",sep="\t",header=FALSE,check.names = F);
    umap_embeddings$Lung = OPcor$V1
    umap_embeddings$MIC = ifelse(umap_embeddings$Peritoneal>5 | umap_embeddings$Liver>5 | umap_embeddings$Lung>5,"Y","N")  
    write.table(umap_embeddings, file = "GSE173958/M1results/umap.txt", sep = "\t", quote = FALSE, row.names = TRUE); 
'

############################################
# 2. Figure3B
#    -Input:umap.txt during Figure 3A and clones.txt
#    -Output:clone distribution
############################################
Rscript -e '
    MIC = read.table("GSE173958/M1results/umap.txt",sep="\t",header=TRUE,check.names = F);
    Clones = read.table("GSE173958/GSM5283482_M1-PT-clones/clones.txt",sep="\t",header=TRUE,check.names = F);
    rownames(Clones) = Clones$Barcode;
    Clones = Clones[rownames(MIC),];
    MIC$Clones = Clones$cloneID;       
    tab = table(MIC[, c("MIC", "Clones")])
    print(c(tab[, "1"] / rowSums(tab),sum(tab[, "1"])/sum(tab)))
'

############################################
# 3. Figure3C
#    -Input:Normalized ScData.pca.rds during DataProcessing.sh
#           umap.txt during Figure3A
#           DEG.txt during DataProcessing.sh
#           mh.all.v2025.1.Mm.symbols.gmt download from MsigDB
#    -Output:PseudoEMT.txt with monocle and EMT scores
############################################

Rscript -e '
    source("../Processing.R");
    library(Seurat); 
    library(monocle);    
    Sys.setenv(OMP_NUM_THREADS = 1,MKL_NUM_THREADS = 1,OPENBLAS_NUM_THREADS = 1,VECLIB_MAXIMUM_THREADS = 1)
    set.seed(1234)
    Data = readRDS("GSE173958/M1results/Primary/ScData.pca.rds");    
    MIC = read.table("GSE173958/M1results/umap.txt",sep="\t",header=TRUE,check.names = F);
    #reduce the clusters with only about 20 cells£¬ different computers may have a little changes
    sccluster(Data, "GSE173958/M1results1/", resolution = 0.8, PCnums = 50, DimNum = 50, perplexity = 30,doubletflag=FALSE)
    scDEGs("GSE173958/M1results1/ScData.pca.rds",OutputFile = "GSE173958/M1results1/DEG.txt",assay = "RNA",clustergroup = "seurat_clusters",logFC = 0)
    DEGData <- read.delim("GSE173958/M1results1/DEG.txt", header = TRUE)
    SigDEGData=subset(DEGData,p_val_adj<0.05 & avg_log2FC>1);
    monoclegenes=intersect(SigDEGData$gene,rownames(Data));
    MonocleData=subset(Data,features=monoclegenes);
    counts = as.matrix(MonocleData@assays$RNA@layers$counts)
    rownames(counts) = rownames(MonocleData)
    colnames(counts) = colnames(MonocleData)
    data <- as(counts, "sparseMatrix");
    pd <- new("AnnotatedDataFrame", data = MonocleData@meta.data);
    fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data));
    fd <- new("AnnotatedDataFrame", data = fData);
    set.seed(1234);              
    mycds <- newCellDataSet(data,phenoData = pd,featureData = fd,expressionFamily = negbinomial.size());
    mycds <- estimateSizeFactors(mycds); #normalize the expression difference among cells.
    mycds <- setOrderingFilter(mycds, monoclegenes);
    mycds <- reduceDimension(mycds, max_components = 2, method = "DDRTree",auto_param_selection=FALSE);
    mycds <- orderCells(mycds);
    saveRDS(mycds,file = "GSE173958/M1results/monocle.rds")            
    MIC$Pseudotime = mycds$Pseudotime;
    MIC$Pseudotime_rev = max(mycds$Pseudotime)-mycds$Pseudotime;
    Data$Pseudotime = MIC$Pseudotime_rev;  #check the positive ones
    pdf(paste0("GSE173958/M1results/PseudoEMT.pdf"))
    figure1=FeaturePlot(Data, features = "Pseudotime", reduction = "umap",  pt.size = 1) +  theme_bw() +
            theme(
              legend.position = "right",
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.border = element_rect(fill = NA, color = "black", linewidth = 1, linetype = "solid"),
              axis.text = element_text(size = 20),
              axis.title.x = element_text(vjust = 0.5, hjust = 0.5),
              axis.title.y = element_text(vjust = 0.5, hjust = 0.5)
            )
    print(figure1)
    dev.off()    
    GeneSets=read.table("GSE173958/mh.all.v2025.1.Mm.symbols.gmt",sep="\t",header=FALSE,fill=TRUE,check.names=FALSE);#download from MSigDB
    eachgeneset = GeneSets[14, ];  #EMT
    eachgeneset = as.character(eachgeneset[3:length(eachgeneset)]);
    eachgeneset = eachgeneset[which(eachgeneset != "")];
    genes = list("genes" = eachgeneset[eachgeneset %in% rownames(Data)]); 
    set.seed(1234);
    eachgenesetscore = AddModuleScore(object = Data, features = genes, name = "EMT",slot = "counts" );
    MIC$EMT = eachgenesetscore$EMT1;
    print(cor.test(MIC$EMT,MIC$Pseudotime,method="pearson"));
    print(cor.test(MIC$EMT,MIC$Pseudotime_rev,method="pearson"));
    write.table(MIC, file = "GSE173958/M1results/PseudoEMT.txt", sep = "\t", quote = FALSE, row.names = TRUE); 
'

############################################
# 4. Figure3D
#    -Input:Normalized ScData.pca.rds during DataProcessing.sh
#           umap.txt during Figure3A
#    -Output:markerexpression.txt
############################################

Rscript -e '
    library(Seurat);
    Data = readRDS("GSE173958/M1results/Primary/ScData.pca.rds");
    PrimaryMatrix=Data@assays$RNA@layers$data;
    colnames(PrimaryMatrix)=colnames(Data);
    rownames(PrimaryMatrix)=rownames(Data);
    genes=c("Epcam","Muc1","Cdh1","Zeb2")
    MIC = read.table("GSE173958/M1results/umap.txt",sep="\t",header=TRUE,check.names = F);
    Data$MIC = ifelse(MIC$Peritoneal>5 | MIC$Liver>5 | MIC$Lung>5,"Y","N")
    Data$MIC1=ifelse(MIC$Peritoneal>0 | MIC$Liver>0 | MIC$Lung>0,"Y","N") 
    MICindex = which(Data@meta.data$MIC == "Y")        
    MICcells=rownames(Data@meta.data[MICindex,]);
    nonMICcells=rownames(Data@meta.data[-MICindex,]);
    MICindex1 = which(Data@meta.data$MIC1 == "Y")        
    MICcells1=rownames(Data@meta.data[MICindex1,]);
    nonMICcells1=rownames(Data@meta.data[-MICindex1,]);             
    for(gene in genes){
      ttest = t.test(PrimaryMatrix[gene,MICcells],PrimaryMatrix[gene,nonMICcells]);
      print(c(5,gene,ttest$p.value,as.numeric(ttest$estimate)));
      ttest = t.test(PrimaryMatrix[gene,MICcells1],PrimaryMatrix[gene,nonMICcells1]);
      print(c(0,gene,ttest$p.value,as.numeric(ttest$estimate)));
    }             
    SigGeneMatrix = data.frame(gene="",nonMIC=0,MIC=0);
    SigGeneMatrix = SigGeneMatrix[-1,];
    for(gene in genes){
      SigGeneMatrix = rbind(SigGeneMatrix,data.frame(gene=rep(gene,length(nonMICcells)),nonMIC=PrimaryMatrix[gene,nonMICcells],
                                       MIC=c(PrimaryMatrix[gene,MICcells],rep("",length(nonMICcells)-length(MICcells)))));
    }              
    write.table(SigGeneMatrix, "GSE173958/M1results/markerexpression.txt", quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE);              
' 


############################################
# 5. Figure3E
#    -Input:Normalized ScData.pca.rds during Figure2C
#           mh.all.v2025.1.Mm.symbols.gmt from MSigDB
#    -Output:cytotrace.score.txt:stemness scores
#            ttest comparison results
############################################

Rscript -e '
    MIC = read.table("GSE173958/M1results/PseudoEMT.txt",sep="\t",header=TRUE,check.names = F);
    MIC$All = MIC$Peritoneal + MIC$Liver + MIC$Lung
    PseudoEMT=data.frame("Pseudotime"=MIC$Pseudotime,"Pseudotime_rev"=MIC$Pseudotime_rev,"MICscore"=MIC$All/max(MIC$All)) #select the Pseudotime positive associated with MICscore
    write.table(PseudoEMT, "GSE173958/M1results/PseudoEMTtrajectory.txt", quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE);       
'