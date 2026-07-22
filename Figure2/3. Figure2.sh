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
cd Figure2/
############################################
# 1. Figure2A
#    -Input:10X-format files of GSE249057
#    -Output:prop_tab.txt
#            probability of cell clusters at
#            different timepoints
############################################

Rscript -e '
    source("../Processing.R"); 
    Sys.setenv(OMP_NUM_THREADS = 1,MKL_NUM_THREADS = 1,OPENBLAS_NUM_THREADS = 1,VECLIB_MAXIMUM_THREADS = 1)
    set.seed(1234)           
    scData0=scRead(Path="GSE249057/GSM7925719_0h/", Label = "0",RNAfeature=200,minRNAcount=1000,maxRNAcount=8000,mt=15,rb=100,hb=100); #Data from 0 step;
    scData6=scRead(Path="GSE249057/GSM7925720_6h/", Label = "6h",RNAfeature=200,minRNAcount=1000,maxRNAcount=8000,mt=15,rb=100,hb=100);
    scData48=scRead(Path="GSE249057/GSM7925721_48h/", Label = "48h",RNAfeature=200,minRNAcount=1000,maxRNAcount=8000,mt=15,rb=100,hb=100);
    scData2m=scRead(Path="GSE249057/GSM7925722_2mo/", Label = "2m",RNAfeature=200,minRNAcount=1000,maxRNAcount=8000,mt=15,rb=100,hb=100);
    scData4m=scRead(Path="GSE249057/GSM7925723_4mo/", Label = "4m",RNAfeature=200,minRNAcount=1000,maxRNAcount=8000,mt=15,rb=100,hb=100);   
    scDataAll = list(scData0,scData6,scData48,scData2m,scData4m);        
    NormalizeAll = scNormalize(scDataAll, NormalizeMethod = "LogNormalize", ScaleFactor = 10000); 
    anchors <- FindIntegrationAnchors(object.list = NormalizeAll, dims = 1:30);
    mergedData <- IntegrateData(anchorset = anchors, dims = 1:30)            
    sccluster(mergedData,"GSE249057/GSM7925723_4mo/");
    AData = readRDS("GSE249057/GSM7925723_4mo/ScData.pca.rds");
    tab <- table(AData@meta.data[, c("Sampleid", "seurat_clusters")]);
    prop_tab <- prop.table(tab, margin = 1); 
    write.table(prop_tab, file = "GSE249057/GSM7925723_4mo/prop_tab.txt", sep = "\t", quote = FALSE, row.names = TRUE); 
'
############################################
# 2. Figure2B
#    -Input:Normalized ScData.pca.rds during Figure2C
#    -Output:monocle.rds
#            trajectory.png
############################################

Rscript -e '
    library(Seurat); 
    library(monocle);    
    Sys.setenv(OMP_NUM_THREADS = 1,MKL_NUM_THREADS = 1,OPENBLAS_NUM_THREADS = 1,VECLIB_MAXIMUM_THREADS = 1)
    set.seed(1234)
    AData = readRDS("GSE249057/GSM7925723_4mo/ScData.pca.rds");
    PData = subset(AData,Sampleid==0)
    PData[["joined"]] <- JoinLayers(PData[["RNA"]])
    DEGs <- FindAllMarkers(PData,assay = "joined",only.pos = TRUE,min.pct = 0.1,logfc.threshold = 0.25)
    SigDEGs <- DEGs #subset(DEGs,p_val_adj<0.05)
    monoclegenes = unique(SigDEGs$gene)
    MonocleData=subset(PData,features=monoclegenes);
    counts = as.matrix(MonocleData@assays$joined@layers$counts)
    rownames(counts) = rownames(MonocleData[["joined"]])
    colnames(counts) = colnames(MonocleData[["joined"]])
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
    mycds$Pseudotime <- max(mycds$Pseudotime)-mycds$Pseudotime
    saveRDS(mycds,file = "GSE249057/GSM7925723_4mo/monocle.rds")    
    png("GSE249057/GSM7925723_4mo/trajectory.png",width=500,height=300)
    cluster_colors=c("0" = "#FBE3D6", "1" = "#C00000","2" = "#D9F2D0","3" = "#C1E5F5","4" = "#F2CFEE", "5" = "#C2F1C8", "6" = "#DCEAF7")
    print(plot_cell_trajectory(mycds,color_by="seurat_clusters") +
     scale_color_manual(values = cluster_colors))
    dev.off()
'

############################################
# 3. Figure2C
#    -Input:Normalized ScData.pca.rds during Figure2A
#           Stemness.txt from CancerSEA
#           the genes from original data are transferred to human genes
#    -Output:Stemness.txt
#            Stemness scores
############################################
wget "http://biocc.hrbmu.edu.cn/CancerSEA/download/signature/Stemness.txt"
Rscript -e '
    library(Seurat);
    Sys.setenv(OMP_NUM_THREADS = 1,MKL_NUM_THREADS = 1,OPENBLAS_NUM_THREADS = 1,VECLIB_MAXIMUM_THREADS = 1)
    set.seed(1234) 
    AData = readRDS("GSE249057/GSM7925723_4mo/ScData.pca.rds");
    PData = subset(AData,Sampleid==0)
    GeneSets=read.table("Stemness.txt",sep="\t",header=TRUE,fill=TRUE,check.names=FALSE);#download from CancerSEA
    eachgeneset = GeneSets$GeneName;
    genes = list("genes" = eachgeneset[eachgeneset %in% rownames(PData)]);
    PData[["joined"]] <- JoinLayers(PData[["RNA"]])
    PrimaryMatrix=PData@assays$RNA@layers$data;
    rownames(PrimaryMatrix) <- rownames(PData[["joined"]])
    colnames(PrimaryMatrix) <- colnames(PData[["joined"]])
    eachgenesetscore = AddModuleScore(object = PData, features = genes, name = "Stem",assay="joined",slot = "data" );
    PStem=data.frame("sample"=eachgenesetscore@meta.data$Sampleid,"group"=eachgenesetscore@meta.data$seurat_clusters,
                     "Stem"=eachgenesetscore@meta.data$Stem1,"ITGA6"=PrimaryMatrix["ITGA6",]);    
    write.table(PStem, file = "GSE249057/GSM7925723_4mo/Stemness.txt", sep = "\t", quote = FALSE, row.names = TRUE); 
    index = which(PStem$group == 1);
    for(i in 0:6){
      index1 = which(PStem$group == i);
      ttest=t.test(PStem$Stem[index],PStem$Stem[index1]);
      print(c(i,ttest$p.value,as.double(ttest$estimate)))
    }
    for(i in 0:6){
      index1 = which(PStem$group == i);
      ttest=t.test(PStem$ITGA6[index],PStem$ITGA6[index1]);
      print(c(i,ttest$p.value,as.double(ttest$estimate)))
    }
'

############################################
# 4. Figure2D
#    -Input:Normalized ScData.pca.rds during Figure2A
#           Stemness.txt from CancerSEA
#           the genes from original data are transferred to human genes
#    -Output:Metastasis.txt
#            Metastasis scores
############################################
wget "http://biocc.hrbmu.edu.cn/CancerSEA/download/signature/Metastasis.txt"

Rscript -e '
    library(Seurat);
    Sys.setenv(OMP_NUM_THREADS = 1,MKL_NUM_THREADS = 1,OPENBLAS_NUM_THREADS = 1,VECLIB_MAXIMUM_THREADS = 1)
    set.seed(1234) 
    AData = readRDS("GSE249057/GSM7925723_4mo/ScData.pca.rds");
    PData = subset(AData,Sampleid==0)
    GeneSets=read.table("Metastasis.txt",sep="\t",header=TRUE,fill=TRUE,check.names=FALSE);#download from CancerSEA
    eachgeneset = GeneSets$GeneName;
    genes = list("genes" = eachgeneset[eachgeneset %in% rownames(PData)]);
    PData[["joined"]] <- JoinLayers(PData[["RNA"]])
    PrimaryMatrix=PData@assays$RNA@layers$data;
    rownames(PrimaryMatrix) <- rownames(PData[["joined"]])
    colnames(PrimaryMatrix) <- colnames(PData[["joined"]])
    eachgenesetscore = AddModuleScore(object = PData, features = genes, name = "Meta",assay="joined",slot = "data" );
    PMeta=data.frame("sample"=eachgenesetscore@meta.data$Sampleid,"group"=eachgenesetscore@meta.data$seurat_clusters,
                     "Meta"=eachgenesetscore@meta.data$Meta1,"CD44"=PrimaryMatrix["CD44",]);    
    write.table(PMeta, file = "GSE249057/GSM7925723_4mo/Metastasis.txt", sep = "\t", quote = FALSE, row.names = TRUE); 
    index = which(PMeta$group == 1);
    for(i in 0:6){
      index1 = which(PMeta$group == i);
      ttest=t.test(PMeta$Meta[index],PMeta$Meta[index1]);
      print(c(i,ttest$p.value,as.double(ttest$estimate)))
    }
    for(i in 0:6){
      index1 = which(PMeta$group == i);
      ttest=t.test(PMeta$CD44[index],PMeta$CD44[index1]);
      print(c(i,ttest$p.value,as.double(ttest$estimate)))
    }
'

############################################
# 5. Figure2E
#    -Input:Normalized ScData.pca.rds during Figure2A
#           the genes from original data are transferred to human genes
#    -Output:EMT.png
############################################


Rscript -e '
    library(Seurat);
    library(ggplot2)
    Sys.setenv(OMP_NUM_THREADS = 1,MKL_NUM_THREADS = 1,OPENBLAS_NUM_THREADS = 1,VECLIB_MAXIMUM_THREADS = 1)
    set.seed(1234) 
    AData = readRDS("GSE249057/GSM7925723_4mo/ScData.pca.rds");
    PData = subset(AData,Sampleid==0)
    PData[["joined"]] <- JoinLayers(PData[["RNA"]])
    PrimaryMatrix=PData@assays$RNA@layers$data;
    rownames(PrimaryMatrix) <- rownames(PData[["joined"]])
    colnames(PrimaryMatrix) <- colnames(PData[["joined"]])
    PEMT=data.frame("sample"=PData@meta.data$Sampleid,"group"=PData@meta.data$seurat_clusters,
                     "CDH1"=PrimaryMatrix["CDH1",],"SNAI2"=PrimaryMatrix["SNAI2",],"TWIST1"=PrimaryMatrix["TWIST1",]); 
    png("GSE249057/GSM7925723_4mo/EMT.png",width=750,height=300)
    print(DotPlot(PData,features = c("CDH1","SNAI2","TWIST1"),group.by = "seurat_clusters",assay = "RNA")+coord_flip()+scale_color_gradient(low = "lightgrey",high = "#C00000"))
    dev.off()
'



############################################
# 6. Figure2F
#    -Input:10X-format files of timepoint 0 of GSE249057
#    -Output:prop_tab.txt
#            probability of cell clusters at
#            different timepoints
############################################

Rscript -e '
    source("../Processing.R");
    Sys.setenv(OMP_NUM_THREADS = 1,MKL_NUM_THREADS = 1,OPENBLAS_NUM_THREADS = 1,VECLIB_MAXIMUM_THREADS = 1)
    set.seed(1234)           
    scData0=scRead(Path="GSE249057/GSM7925719_0h/", Label = "0",RNAfeature=200,minRNAcount=1000,maxRNAcount=8000,mt=15,rb=100,hb=100); #Data from 0 step;          
    Normalize0 = scNormalize(scData0, NormalizeMethod = "LogNormalize", ScaleFactor = 10000);  
    sccluster(Normalize0,"GSE249057/GSM7925719_0h/"); 
    PData = readRDS("GSE249057/GSM7925719_0h/ScData.pca.rds");
    umap_embeddings <- Embeddings(PData, "umap");
    umap_embeddings = as.data.frame(umap_embeddings[rownames(PData@meta.data),]);
    umap_embeddings$seurat_clusters = as.character(PData@meta.data$seurat_clusters);
    write.table(umap_embeddings, file = "GSE249057/GSM7925719_0h/umap_embeddings.txt", sep = "\t", quote = FALSE, row.names = TRUE);  
    AData = readRDS("GSE249057/GSM7925723_4mo/ScData.pca.rds");
    PmetaAll = subset(AData@meta.data,Sampleid==0);
    Pmeta = PData@meta.data;
    rownames(Pmeta) = paste0(rownames(Pmeta),"_1");
    Pmeta = Pmeta[rownames(PmetaAll),];    
    Pmeta$Allcluster = PmetaAll$seurat_clusters;
    write.table(table(Pmeta[,c("seurat_clusters","Allcluster")]),file = "GSE249057/GSM7925719_0h/clusterfraction.txt", sep = "\t", quote = FALSE, row.names = TRUE);
' 


############################################
# 7. Figure2G
#    -Input:expression count during DataProcessing.sh
#    -Output:proximity_OT_test.csv: scMIC prediction
#            OTtraining.txt: training process
############################################

embeddings=("PCA" "None" "AE" "PCA" "PCA" "PCA" "PCA" "PCA" "PCA")
distances=("Euclidean" "Euclidean" "Euclidean" "Pearson" "Spearman" "Euclidean" "Euclidean" "Euclidean" "Euclidean")
topks=(1 1 1 1 1 1 2 10 50)
uotflags=("T" "T" "T" "T" "T" "F" "T" "T" "T")
MetFile="GSE249057/Results/Metastasis/MetastasisCount.csv"
PriFile="GSE249057/Results/Primary/PrimaryCount.csv"
outFile="GSE249057/Results/proximity_OT_test.csv"
rm -f GSE249057/Results/OTtraining.txt
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