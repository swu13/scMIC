#!/bin/bash
set -euo pipefail

############################################
# scMIC project
# Step 0. Download and organize scRNA-seq data 
#         and construct scMIC framework for 
#         MIC identification
# Dataset: GSE178318,OMIX002487 and GSE277783
############################################

# Project root directory
RootPath="/data/twang15/wsj/scMIC"
cd "${RootPath}"

mkdir Figure4/
cd Figure4/
# create environment
conda env create -f 1.\ scMIC_Figure4.yml
Rscript -e '
    devtools::install_github("kueckelj/confuns");
    devtools::install_github("theMILOlab/SPATAData");
    devtools::install_github("theMILOlab/SPATA2");
    devtools::install_github("linxihui/NNLM")
    devtools::install_github("wguo-research/scCancer", ref = "scCancer2", upgrade = "never");
    remove.packages("xgboost");
    install.packages("https://cran.r-project.org/src/contrib/Archive/xgboost/xgboost_1.5.2.1.tar.gz",repos = NULL,type = "source")
    devtools::install_github("navinlabcode/copykat")
    BiocManager::install("DropletUtils")
'

mkdir scCancer2/
wget -P scCancer2/ https://github.com/czythu/scCancer/archive/refs/heads/master.zip
unzip scCancer2/master.zip -d scCancer2/

############################################
# 0.1 GSE178318 Data download
############################################

mkdir -p GSE178318
cd GSE178318
wget -O barcodes.tsv.gz https://ftp.ncbi.nlm.nih.gov/geo/series/GSE178nnn/GSE178318/suppl/GSE178318_barcodes.tsv.gz
wget -O features.tsv.gz https://ftp.ncbi.nlm.nih.gov/geo/series/GSE178nnn/GSE178318/suppl/GSE178318_genes.tsv.gz
wget -O matrix.mtx.gz https://ftp.ncbi.nlm.nih.gov/geo/series/GSE178nnn/GSE178318/suppl/GSE178318_matrix.mtx.gz
cd ../

############################################
# 0.2 OMIX002487 Data download
############################################

mkdir OMIX002487
wget -P OMIX002487/ https://download.cncb.ac.cn/OMIX/OMIX002487/OMIX002487-01.csv
wget -P OMIX002487/ https://download.cncb.ac.cn/OMIX/OMIX002487/OMIX002487-02.txt

############################################
# 0.3 GSE277783 Data download
############################################

mkdir -p GSE277783/
cd GSE277783
wget -O GSE277783_RAW.tar  "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE277783&format=file"
tar -xvf GSE277783_RAW.tar
ls ./*tar.gz | while read file
do
    id=`echo $file | awk -F"/" '{split($NF,a,"_");print a[1]"_"a[2]}'`
    mkdir $id"/"
    tar -xzvf $file -C $id"/"
    mv $id"/"*"_spatial" $id"/spatial/"
    mv $id*"_filtered_feature_bc_matrix.h5" $id"/filtered_feature_bc_matrix.h5"
done
rm -f ./*.gz
cd ../



############################################
# 0.3 Prepare each sample 10X data of GSE178318
############################################

cd GSE178318
Rscript -e '
    library(Seurat);
    Data=Read10X("./", gene.column = 2);
    scData <- CreateSeuratObject(Data)
    sampleid=t(as.data.frame(strsplit(rownames(scData@meta.data),"_")));
    scData@meta.data$Sampleid=sampleid[,2];
    scData@meta.data$Groupid=sampleid[,3];
    samples=unique(scData@meta.data$Sampleid);
    for(sample in samples){
      EachData=subset(scData,subset=Sampleid==sample);
      organs=unique(EachData@meta.data$Groupid);
      for(organ in organs){
        OrganData=subset(EachData,subset=Groupid==organ);
        OrganCount=OrganData@assays$RNA@layers$counts;
        rownames(OrganCount)=rownames(OrganData);
        colnames(OrganCount)=colnames(OrganData);
        DropletUtils:::write10xCounts(paste0("./",sample,"_",organ,"/"),OrganCount,version="3")
      }
    }
'
cd ../

############################################
# 0.4 Prepare each sample 10X data of OMIX002487
############################################

cd OMIX002487
Rscript -e '
    library(Seurat)
    metadata=read.csv("OMIX002487-01.csv",row.names=1);
    Data=read.table("OMIX002487-02.txt",header=TRUE,check.names = F)
    Epithelial_metadata=subset(metadata,label=="Epithelial")
    for(patient in unique(Epithelial_metadata$group_patient)){
      for(tissue in unique(Epithelial_metadata$group_tissue)){
        each_metadata = subset(Epithelial_metadata,group_patient == patient & group_tissue==tissue)
        cells=rownames(each_metadata)
        each_data=Data[,cells]
        DropletUtils:::write10xCounts(paste0(patient,"_",tissue,"/"),as(each_data,"sparseMatrix"),version="3");
      }
    }
'
rm -r P6*/ #due to small number of tumor cells
cd ../

############################################
# 0.5 Prepare each sample 10X data of GSE277783
############################################

cd GSE277783
Rscript -e '
    library(Seurat)
    for(i in 1:13){
      folders <- list.files("./",pattern = paste0("Pt-", i, "[A-E]$"),full.names = TRUE)
      for(folder in folders){
        spData=Load10X_Spatial(folder);
        print(mean(spData@meta.data$nCount_Spatial))
        matrix=spData@assays$Spatial@layers$counts;
        rownames(matrix)=rownames(spData);
        colnames(matrix)=colnames(spData);
        DropletUtils:::write10xCounts(paste0(folder,"/10X/"),matrix,version="3");           
      }
    }
' > "count.txt" 
rm -r GSM8452848_Pt-1B/ #due to mean count < 200
cd ../

############################################
# 1. Primary and metastatic tumor
#    gene expression preparation of GSE178318
############################################

############################################
# 1.1 Primary tumor processing of GSE178318
#     - QC, normalization
#     - epithelial cluster identification
############################################

zcat "GSE178318/features.tsv.gz" | awk -F" " '{print $2"\t"$1}' > "GSE178318/genelist.txt"
Rscript -e '
    source("../Processing.R");
    library(Seurat)
    crc_dirs <- list.files(path = "GSE178318/",pattern = "_CRC$",full.names = TRUE)
    Sys.setenv(OMP_NUM_THREADS = 1,MKL_NUM_THREADS = 1,OPENBLAS_NUM_THREADS = 1,VECLIB_MAXIMUM_THREADS = 1)
    set.seed(1234)
    for(i in 1:length(crc_dirs)){
      PrimaryPath = crc_dirs[i]
      PrimaryData <- scRead(PrimaryPath, Label = "Primary",RNAfeature=200,minRNAcount=500,maxRNAcount=20000,mt=15,rb=100,hb=100);
      PrimaryNormalizedData <- scNormalize(PrimaryData, NormalizeMethod = "LogNormalize", ScaleFactor = 10000)
      OutputPath=sub("_CRC$", "", PrimaryPath)
      ensure_dir(OutputPath)
      sccluster(PrimaryNormalizedData, paste0(OutputPath, "/Primary/"), resolution = 0.9, PCnums = 50, DimNum = 50, perplexity = 30,doubletflag=FALSE)
      scRNAanno(paste0(OutputPath, "/Primary/ScData.pca.rds"), CellMarkerFile="GSE178318/GSE178318_markers.txt", paste0(OutputPath, "/Primary/"))
    }
'

############################################
# 1.1 Metastatic tumor processing of GSE178318
#     - QC, normalization
#     - epithelial cluster identification
############################################

Rscript -e '
    source("../Processing.R");
    library(Seurat)
    crc_dirs <- list.files(path = "GSE178318/",pattern = "_LM$",full.names = TRUE)
    Sys.setenv(OMP_NUM_THREADS = 1,MKL_NUM_THREADS = 1,OPENBLAS_NUM_THREADS = 1,VECLIB_MAXIMUM_THREADS = 1)
    set.seed(1234)
    for(i in 1:length(crc_dirs)){
      MetastasisPath = crc_dirs[i]
      MetastasisData <- scRead(MetastasisPath, Label = "Metastasis",RNAfeature=200,minRNAcount=500,maxRNAcount=20000,mt=15,rb=100,hb=100)
      MetastasisNormalizedData <- scNormalize(MetastasisData, NormalizeMethod = "LogNormalize", ScaleFactor = 10000)  
      OutputPath=sub("_LM$", "/", MetastasisPath)      
      sccluster(MetastasisNormalizedData, paste0(OutputPath, "Metastasis"), resolution = 0.9, PCnums = 50, DimNum = 50, perplexity = 30,doubletflag=FALSE)
      scRNAanno(paste0(OutputPath, "/Metastasis/ScData.pca.rds"), CellMarkerFile="GSE178318/GSE178318_markers.txt", paste0(OutputPath, "Metastasis")) 
    }
'

############################################
# 1.3 Primary tumor processing (CNV/XGboost) of GSE178318 
#     and gene expression matrix construction
############################################


Rscript -e '
    source("../Processing.R");
    library(scCancer)
    library(Seurat)
    crc_dirs <- list.files(path = "GSE178318/",pattern = "^COL[0-9]+$",full.names = TRUE)
    Sys.setenv(OMP_NUM_THREADS = 1,MKL_NUM_THREADS = 1,OPENBLAS_NUM_THREADS = 1,VECLIB_MAXIMUM_THREADS = 1)
    set.seed(1234)
    PTIndexs=c("0,3,7,16","14","7,11","2,6","15","14") #may have a little change,please see the CellmarkersDotPlot.pdf to determine the epithelial cluster
    for(i in 1:length(crc_dirs)){
      OutputPath = crc_dirs[i]
      Primary <- readRDS(paste0(OutputPath, "/Primary/ScData.pca.rds"))
      PTIndex <- as.numeric(unlist(strsplit(PTIndexs[i], ",")))
      PrimaryTumour <- subset(Primary, subset = seurat_clusters %in% PTIndex) 
      if(ncol(PrimaryTumour)>100){
        PrimaryTumour_select = PrimaryTumour
        cell.annotation <- data.frame(barcodes = colnames(PrimaryTumour_select), Cell.Type = "Epithelial")
        rownames(cell.annotation) <- colnames(PrimaryTumour_select)
        umap_coords <- Embeddings(PrimaryTumour_select, reduction = "umap")
        cell.annotation <- data.frame(cell.annotation, tSNE_1 = umap_coords[,"umap_1"], tSNE_2 = umap_coords[,"umap_2"])
        cell.annotation$Cluster=PrimaryTumour_select@meta.data$seurat_clusters
        PrimaryTumour_select@assays$RNA <- as(PrimaryTumour_select@assays$RNA, "Assay")
        cells = rownames(cell.annotation)
        genelist = read.table("GSE178318/genelist.txt",sep="\t",header=FALSE,check.names = F)
        colnames(genelist) = c("SYMBOL","EnsemblID")
        genelist <- genelist[!duplicated(genelist$SYMBOL), ]
        rownames(genelist)=genelist$SYMBOL      
        genes = intersect(rownames(PrimaryTumour_select),rownames(genelist))
        genelist = genelist[genes,]
        PrimaryTumour_select = subset(PrimaryTumour_select, features = genes)
        ensure_dir(paste0(OutputPath, "/Primary/figures/"))
        malignantscore <- runMalignancy(PrimaryTumour_select, genelist, cell.annotation = cell.annotation, savePath = paste0(OutputPath, "/Primary/"))
        cell.annotation$CNVscore = malignantscore$cell.annotation$Malign.score
        cell.annotation$CNVtype = malignantscore$cell.annotation$Malign.type
        cells_cnv = intersect(cells,rownames(subset(cell.annotation,CNVtype=="malignant")))
        malignantscore <- predMalignantCell(PrimaryTumour_select,cell.annotation = cell.annotation,ModelPath = "scCancer2/scCancer-master/inst/txt/", MALIGNANT.THRES = 0.5)
        cell.annotation$XGBscore = malignantscore$Malign.score
        cell.annotation$XGBtype = malignantscore$Malign.type
        cells_XGBoost = intersect(cells,rownames(subset(cell.annotation,XGBtype=="malignant")))
        write.table(cell.annotation, paste0(OutputPath, "/Primary/maligant.txt"), quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
        cells = union(cells_cnv,cells_XGBoost)
        PrimaryTumour = subset(PrimaryTumour,cells=cells)
        sccluster(PrimaryTumour, paste0(OutputPath, "/Primary/Tumour/"), resolution = 0.9, PCnums = 50, DimNum = 50, perplexity = 30,doubletflag=FALSE) 
        scDEGs(paste0(OutputPath, "/Primary/Tumour/ScData.pca.rds"),OutputFile = paste0(OutputPath, "/Primary/Tumour/DEG.txt"),assay = "RNA",clustergroup = "seurat_clusters",logFC = 0)
      }
    }
'
rm -r GSE178318/COL12*  #do not have enough epithelial cells for primary samples
rm -r GSE178318/COL17*  #do not have enough epithelial cells for primary samples
rm -r GSE178318/COL18*  #do not have enough epithelial cells for primary samples

############################################
# 1.4 Metastatic tumor processing (CNV/XGboost) of GSE178318 
#     and gene expression matrix construction
############################################

Rscript -e '
    source("../Processing.R");
    library(scCancer)
    library(Seurat)
    crc_dirs <- c("GSE178318//COL07","GSE178318//COL15","GSE178318//COL16")
    Sys.setenv(OMP_NUM_THREADS = 1,MKL_NUM_THREADS = 1,OPENBLAS_NUM_THREADS = 1,VECLIB_MAXIMUM_THREADS = 1)
    set.seed(1234)
    MTIndexs=c("5,16","6","10") #may have a little change,please see the CellmarkersDotPlot.pdf to determine the epithelial cluster
    for(i in 1:length(crc_dirs)){
      OutputPath = crc_dirs[i]
      Metastasis <- readRDS(paste0(OutputPath, "/Metastasis/ScData.pca.rds"))
      MTIndex <- as.numeric(unlist(strsplit(MTIndexs[i], ",")))
      MetastasisTumour <- subset(Metastasis, subset = seurat_clusters %in% MTIndex) 
      if(ncol(MetastasisTumour)>100){
        MetastasisTumour_select = MetastasisTumour
        cell.annotation <- data.frame(barcodes = colnames(MetastasisTumour_select), Cell.Type = "Epithelial")
        rownames(cell.annotation) <- colnames(MetastasisTumour_select)
        umap_coords <- Embeddings(MetastasisTumour_select, reduction = "umap")
        cell.annotation <- data.frame(cell.annotation, tSNE_1 = umap_coords[,"umap_1"], tSNE_2 = umap_coords[,"umap_2"])
        cell.annotation$Cluster=MetastasisTumour_select@meta.data$seurat_clusters
        MetastasisTumour_select@assays$RNA <- as(MetastasisTumour_select@assays$RNA, "Assay")
        cells = rownames(cell.annotation)
        genelist = read.table("GSE178318/genelist.txt",sep="\t",header=FALSE,check.names = F)
        colnames(genelist) = c("SYMBOL","EnsemblID")
        genelist <- genelist[!duplicated(genelist$SYMBOL), ]
        rownames(genelist)=genelist$SYMBOL
        genes = intersect(rownames(MetastasisTumour_select),rownames(genelist))
        genelist = genelist[genes,]
        MetastasisTumour_select = subset(MetastasisTumour_select, features = genes)
        ensure_dir(paste0(OutputPath, "/Metastasis/figures/"))
        malignantscore <- runMalignancy(MetastasisTumour_select, genelist, cell.annotation = cell.annotation, savePath = paste0(OutputPath, "/Metastasis/"))
        cell.annotation$CNVscore = malignantscore$cell.annotation$Malign.score
        cell.annotation$CNVtype = malignantscore$cell.annotation$Malign.type
        cells_cnv = intersect(cells,rownames(subset(cell.annotation,CNVtype=="malignant")))
        malignantscore <- predMalignantCell(MetastasisTumour_select,cell.annotation = cell.annotation,ModelPath = "scCancer2/scCancer-master/inst/txt/", MALIGNANT.THRES = 0.5)
        cell.annotation$XGBscore = malignantscore$Malign.score
        cell.annotation$XGBtype = malignantscore$Malign.type
        cells_XGBoost= intersect(cells,rownames(subset(cell.annotation,XGBtype=="malignant")))
        write.table(cell.annotation, paste0(OutputPath, "/Metastasis/maligant.txt"), quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)      
        cells = union(cells_cnv,cells_XGBoost)
        MetastasisTumour = subset(MetastasisTumour,cells=cells)
        PrimaryTumour <- readRDS(paste0(OutputPath, "/Primary/Tumour/ScData.pca.rds"))
        DEGData <- read.delim(paste0(OutputPath, "/Primary/Tumour/DEG.txt"), header = TRUE)
        SigDEGData <- subset(DEGData, p_val < 0.05 & avg_log2FC > 0)
        PrimaryCount <- as.matrix(PrimaryTumour@assays$RNA@layers$data)
        rownames(PrimaryCount) <- rownames(PrimaryTumour)
        colnames(PrimaryCount) <- colnames(PrimaryTumour)    
        MetastasisCount <- as.matrix(MetastasisTumour@assays$RNA@layers$data)
        rownames(MetastasisCount) <- rownames(MetastasisTumour)
        colnames(MetastasisCount) <- colnames(MetastasisTumour)
        genelist <- intersect(unique(SigDEGData$gene), intersect(rownames(PrimaryCount), rownames(MetastasisCount)))  
        PrimaryCount <- PrimaryCount[genelist, ]
        MetastasisCount <- MetastasisCount[genelist, ]
        write.csv(PrimaryCount, paste0(OutputPath, "/Primary/PrimaryCount.csv"), quote = FALSE)  
        write.csv(MetastasisCount, paste0(OutputPath, "/Metastasis/MetastasisCount.csv"), quote = FALSE)
      }
    }
'

############################################
# 2.1 Primary tumor processing of OMIX002487
#     - QC, normalization
#     - epithelial cell identification
############################################

#get genelist.txt: SYMBOL ENSGID
wget -P OMIX002487/ https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_49/gencode.v49.annotation.gtf.gz --no-check-certificate
gunzip OMIX002487/gencode.v49.annotation.gtf.gz
zcat OMIX002487/P1_Primary/features.tsv.gz | awk -F"\t" '{print $2}' > "OMIX002487/genelist.txt"
awk -F'\t' 'ARGIND==1 && $3=="gene" {match($9, /gene_id "([^"]+)"/, gid);match($9, /gene_name "([^"]+)"/, gname);\
                                    if (gid[1] && gname[1]) {split(gid[1], parts, ".");gene[gname[1]] = parts[1];}}
            ARGIND==2 && ($1 in gene){print $1"\t"gene[$1]}'\
            "OMIX002487/gencode.v49.annotation.gtf"\
            "OMIX002487/genelist.txt"  > "OMIX002487/genelist.1.txt"
mv "OMIX002487/genelist.1.txt" "OMIX002487/genelist.txt"

Rscript -e '
    source("../Processing.R");
    library(Seurat)
    Sys.setenv(OMP_NUM_THREADS = 1,MKL_NUM_THREADS = 1,OPENBLAS_NUM_THREADS = 1,VECLIB_MAXIMUM_THREADS = 1)
    set.seed(1234)
    dirs <- c("P1","P2","P3","P4","P5")
    for(i in 1:length(dirs)){
      PrimaryPath = paste0("OMIX002487/",dirs[i],"_Primary/")
      PrimaryData <- scRead(PrimaryPath, Label = "Primary",RNAfeature=0,minRNAcount=0,maxRNAcount=200000000,mt=100,rb=100,hb=100);
      PrimaryTumour <- scNormalize(PrimaryData, NormalizeMethod = "LogNormalize", ScaleFactor = 10000)
      OutputPath=paste0("OMIX002487/",dirs[i],"/")
      ensure_dir(OutputPath)
      sccluster(PrimaryTumour, paste0(OutputPath, "Primary"), resolution = 0.9, PCnums = 50, DimNum = 50, perplexity = 30,doubletflag=FALSE)
      scRNAanno(paste0(OutputPath, "/Primary/ScData.pca.rds"), CellMarkerFile="OMIX002487/OMIX002487_markers.txt", paste0(OutputPath, "/Primary/"))
    }
'

############################################
# 2.1 Primary tumor processing of OMIX002487
#     - tumor cell identification by copykat
############################################

Rscript -e '
    source("../Processing.R");
    library(Seurat)
    dirs <- c("P1","P2","P3","P4","P5")
    Sys.setenv(OMP_NUM_THREADS = 1,MKL_NUM_THREADS = 1,OPENBLAS_NUM_THREADS = 1,VECLIB_MAXIMUM_THREADS = 1)
    set.seed(1234)
    PTIndexs=c("0,1,2,3,4,5,6","0,1,2,3,4,6","0,1,2,3,4,5","0,1,2,3,4,5,6,8","0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15")
    library(copykat)
    for(i in 1:length(dirs)){
      OutputPath=paste0("OMIX002487/",dirs[i],"/")
      Primary=readRDS(paste0(OutputPath, "/Primary/ScData.pca.rds"))
      PTIndex <- as.numeric(unlist(strsplit(PTIndexs[i], ",")))
      PrimaryTumour <- subset(Primary, subset = seurat_clusters %in% PTIndex) 
      PrimaryMatrix = PrimaryTumour@assays$RNA@layers$counts
      rownames(PrimaryMatrix) =rownames(PrimaryTumour)
      colnames(PrimaryMatrix) =colnames(PrimaryTumour)
      ck <- copykat(rawmat = as.matrix(PrimaryMatrix),id.type = "S",ngene.chr = 0,n.cores=1,sam.name=paste0(OutputPath, "/Primary/"))
      rownames(ck$prediction)=ck$prediction$cell.names
      PrimaryTumour@meta.data$copykat.pred = ck$prediction[rownames(PrimaryTumour@meta.data),]$copykat.pred  
      write.table(PrimaryTumour@meta.data, paste0(OutputPath, "/Primary/maligant.txt"), quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)   
      PrimaryTumour = subset(PrimaryTumour,copykat.pred == "aneuploid")
      sccluster(PrimaryTumour, paste0(OutputPath, "/Primary/Tumour/"), resolution = 0.9, PCnums = 50, DimNum = 50, perplexity = 30,doubletflag=FALSE) 
      scDEGs(paste0(OutputPath, "/Primary/Tumour/ScData.pca.rds"),OutputFile = paste0(OutputPath, "/Primary/Tumour/DEG.txt"),assay = "RNA",clustergroup = "seurat_clusters",logFC = 0)
    }
'

############################################
# 2.3 Metastasis tumor processing of OMIX002487
#     - QC, normalization
#     - epithelial cell identification
############################################

Rscript -e '
    source("../Processing.R");
    library(Seurat)
    Sys.setenv(OMP_NUM_THREADS = 1,MKL_NUM_THREADS = 1,OPENBLAS_NUM_THREADS = 1,VECLIB_MAXIMUM_THREADS = 1)
    set.seed(1234)
    dirs <- c("P1","P2","P3","P4","P5")
    for(i in 1:length(dirs)){
        MetastasisPath = paste0("OMIX002487/",dirs[i],"_Metastasis/")
        MetastasisData <- scRead(MetastasisPath, Label = "Metastasis",RNAfeature=0,minRNAcount=0,maxRNAcount=200000000,mt=100,rb=100,hb=100)
        MetastasisNormalizedData <- scNormalize(MetastasisData, NormalizeMethod = "LogNormalize", ScaleFactor = 10000)  
        OutputPath=paste0("OMIX002487/",dirs[i],"/")
        sccluster(MetastasisNormalizedData, paste0(OutputPath, "Metastasis"), resolution = 0.9, PCnums = 50, DimNum = 50, perplexity = 30,doubletflag=FALSE)
        scRNAanno(paste0(OutputPath, "/Metastasis/ScData.pca.rds"), CellMarkerFile="OMIX002487/OMIX002487_markers.txt", paste0(OutputPath, "/Metastasis/")) 
    }
'


############################################
# 2.4 Metastasis tumor processing of OMIX002487
#     - tumor cell identification by copykat
############################################

Rscript -e '
    source("../Processing.R");
    library(Seurat)
    library(copykat)
    Sys.setenv(OMP_NUM_THREADS = 1,MKL_NUM_THREADS = 1,OPENBLAS_NUM_THREADS = 1,VECLIB_MAXIMUM_THREADS = 1)
    set.seed(1234)
    dirs <- c("P1","P2","P3","P4","P5")
    MTIndexs=c("0,1,2","0,2,3","0,1,2,3,4","0,1,2,3,4,5,6,7,10","0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15")
    for(i in 1:length(dirs)){
        OutputPath=paste0("OMIX002487/",dirs[i],"/") 
        Metastasis=readRDS(paste0(OutputPath, "/Metastasis/ScData.pca.rds"))
        MTIndex <- as.numeric(unlist(strsplit(MTIndexs[i], ",")))
        MetastasisTumour <- subset(Metastasis, subset = seurat_clusters %in% MTIndex) 
        MetastasisMatrix = MetastasisTumour@assays$RNA@layers$counts
        rownames(MetastasisMatrix) =rownames(MetastasisTumour)
        colnames(MetastasisMatrix) =colnames(MetastasisTumour)
        ck <- copykat(rawmat = as.matrix(MetastasisMatrix),id.type = "S",ngene.chr = 0,n.cores=1,sam.name=paste0(OutputPath, "/Metastasis/"))
        rownames(ck$prediction)=ck$prediction$cell.names
        MetastasisTumour@meta.data$copykat.pred = ck$prediction[rownames(MetastasisTumour@meta.data),]$copykat.pred 
        write.table(MetastasisTumour@meta.data, paste0(OutputPath, "/Metastasis/maligant.txt"), quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)   
        MetastasisTumour = subset(MetastasisTumour,copykat.pred == "aneuploid")
        PrimaryTumour <- readRDS(paste0(OutputPath, "/Primary/Tumour/ScData.pca.rds"))
        DEGData <- read.delim(paste0(OutputPath, "/Primary/Tumour/DEG.txt"), header = TRUE)
        SigDEGData <- subset(DEGData, p_val < 0.05 & avg_log2FC > 0)
        PrimaryCount <- as.matrix(PrimaryTumour@assays$RNA@layers$data)
        rownames(PrimaryCount) <- rownames(PrimaryTumour)
        colnames(PrimaryCount) <- colnames(PrimaryTumour)    
        MetastasisCount <- as.matrix(MetastasisTumour@assays$RNA@layers$data)
        rownames(MetastasisCount) <- rownames(MetastasisTumour)
        colnames(MetastasisCount) <- colnames(MetastasisTumour)
        genelist <- intersect(unique(SigDEGData$gene), intersect(rownames(PrimaryCount), rownames(MetastasisCount)))  
        PrimaryCount <- PrimaryCount[genelist, ]
        MetastasisCount <- MetastasisCount[genelist, ]
        write.csv(PrimaryCount, paste0(OutputPath, "/Primary/PrimaryCount.csv"), quote = FALSE)  
        write.csv(MetastasisCount, paste0(OutputPath, "/Metastasis/MetastasisCount.csv"), quote = FALSE)
    }
'

rm -r  OMIX002487/P1* #delelte due to the limited number of tumor cells for metastasis samples

############################################
# 3. Primary and metastatic tumor
#    gene expression preparation of GSE277783
############################################

#get genelist.txt: SYMBOL ENSGID
zcat GSE277783/GSM8452856_Pt-4A/10X/features.tsv.gz | awk -F"\t" '{print $2}' > "GSE277783/genelist.txt"
awk -F'\t' 'ARGIND==1 && $3=="gene" {match($9, /gene_id "([^"]+)"/, gid);match($9, /gene_name "([^"]+)"/, gname);\
                                    if (gid[1] && gname[1]) {split(gid[1], parts, ".");gene[gname[1]] = parts[1];}}
            ARGIND==2 && ($1 in gene){print $1"\t"gene[$1]}'\
            "OMIX002487/gencode.v49.annotation.gtf"\
            "GSE277783/genelist.txt"  > "GSE277783/genelist.1.txt"
mv "GSE277783/genelist.1.txt" "GSE277783/genelist.txt"

############################################
# 3.1 tumor processing (CNV) of GSE277783 
#     and gene expression matrix construction
############################################

Rscript -e '
    library(Seurat);
    library(ggplot2);
    library(SPATA2);
    Sys.setenv(OMP_NUM_THREADS = 1,MKL_NUM_THREADS = 1,OPENBLAS_NUM_THREADS = 1,VECLIB_MAXIMUM_THREADS = 1)
    set.seed(1234)
    source("../Processing.R");
    crc_dirs <- list.files(path = "GSE277783/",pattern = "[A-E]$",full.names = TRUE)
    genes=c("EPCAM","KRT8","KRT19","KRT18");
    for(i in 1:length(crc_dirs)){
      Num = substr(crc_dirs[i], nchar(crc_dirs[i]), nchar(crc_dirs[i]))
      subfolder = paste0("Metastasis",Num)
      if(Num == "A"){
        subfolder = "Primary"
      }
      Path = paste0(crc_dirs[i],"/10X/")
      Data <- scRead(Path, Label = subfolder,RNAfeature=0,minRNAcount=100,maxRNAcount=20000000,mt=100,rb=100,hb=100);
      NormalizedData <- scNormalize(Data, NormalizeMethod = "LogNormalize", ScaleFactor = 10000)
      NormalizedData <- AddModuleScore(object = NormalizedData,features = list(genes),name = "TumorMarker")      
      spData=Load10X_Spatial(crc_dirs[i]);                    
      spata2Data=asSPATA2(spData,sample_name = subfolder, platform = "VisiumSmall",assay_modality="gene");
      OutputPath=paste0("GSE277783/",sub(".*_(Pt-\\\\d+).*", "\\\\1", crc_dirs[i])) #OutputPath=paste0("GSE277783/",sub(".*_(Pt-\\d+).*", "\\1", crc_dirs[i]))
      ensure_dir(paste0(OutputPath,"/",subfolder,"/cnvresults/"))
      spata2Data=runCNV(object=spata2Data,directory_cnv_folder=paste0(OutputPath,"/",subfolder,"/cnvresults/"),cnv_prefix="Chr");
      cnvData = readRDS(paste0(OutputPath,"/",subfolder,"/cnvresults/infercnv-obj.RDS"))
      refnum=length(cnvData@reference_grouped_cell_indices$ref)
      cnvscore=data.frame(cnvscore=colMeans((cnvData@expr.data[,1:(dim(cnvData@expr.data)[2]-refnum)]-1)^2))
      refscore=data.frame(cnvscore=colMeans((cnvData@expr.data[,(dim(cnvData@expr.data)[2]-refnum+1):dim(cnvData@expr.data)[2]]-1)^2))
      NormalizedData$cnvscore = "NA"
      NormalizedData@meta.data[rownames(cnvscore),]$cnvscore = cnvscore$cnvscore
      saveRDS(NormalizedData,paste0(OutputPath,"/",subfolder,"/ScData.pca.rds"))
      saveRDS(spata2Data,paste0(OutputPath,"/",subfolder,"/spata2Data.rds"))
      write.table(refscore, paste0(OutputPath,"/",subfolder,"/refcnv.txt"), quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
    }
'

############################################
# 3.2 Primary tumor processing (CNV/XGBoost) of GSE277783 
#     and gene expression matrix construction
############################################

Rscript -e '
    library(Seurat);
    library(ggplot2);
    library(scCancer)
    source("../Processing.R");
    Sys.setenv(OMP_NUM_THREADS = 1,MKL_NUM_THREADS = 1,OPENBLAS_NUM_THREADS = 1,VECLIB_MAXIMUM_THREADS = 1)
    set.seed(1234)
    crc_dirs <- list.files(path = "GSE277783/",pattern = "[A]$",full.names = TRUE)
    for(i in 1:length(crc_dirs)){
      OutputPath=paste0("GSE277783/",sub(".*_(Pt-\\\\d+).*", "\\\\1", crc_dirs[i])) #OutputPath=paste0("GSE277783/",sub(".*_(Pt-\\d+).*", "\\1", crc_dirs[i]))
      subfolder = "Primary"
      Primary <- readRDS(paste0(OutputPath,"/",subfolder,"/ScData.pca.rds"))
      sccluster(Primary, paste0(OutputPath,"/",subfolder,"/"), resolution = 0.9, PCnums = 50, DimNum = 50, perplexity = 30,doubletflag=FALSE)
      Primary <- readRDS(paste0(OutputPath,"/",subfolder,"/ScData.pca.rds"))
      refscore = read.table(paste0(OutputPath,"/",subfolder,"/refcnv.txt"),sep="\t",header=TRUE,check.names = F)
      Primary$Tumor = ifelse(Primary$TumorMarker1>0 & Primary$cnvscore>max(refscore) & Primary$TumorMarker1!="NA" & Primary$cnvscore!="NA","T","N")
      PrimaryTumour = subset(Primary,Tumor=="T")
      PrimaryTumour_select = PrimaryTumour      
      cell.annotation <- data.frame(barcodes = colnames(PrimaryTumour_select), Cell.Type = "Epithelial")
      rownames(cell.annotation) <- colnames(PrimaryTumour_select)
      umap_coords <- Embeddings(PrimaryTumour_select, reduction = "umap")
      cell.annotation <- data.frame(cell.annotation, tSNE_1 = umap_coords[,"umap_1"], tSNE_2 = umap_coords[,"umap_2"])
      cell.annotation$Cluster=PrimaryTumour_select@meta.data$seurat_clusters
      PrimaryTumour_select@assays$RNA <- as(PrimaryTumour_select@assays$RNA, "Assay")
      cells = rownames(cell.annotation)
      genelist = read.table("GSE277783/genelist.txt",sep="\t",header=FALSE,check.names = F)
      colnames(genelist) = c("SYMBOL","EnsemblID")
      genelist <- genelist[!duplicated(genelist$SYMBOL), ]
      rownames(genelist)=genelist$SYMBOL      
      genes = intersect(rownames(PrimaryTumour_select),rownames(genelist))
      genelist = genelist[genes,]
      PrimaryTumour_select = subset(PrimaryTumour_select, features = genes)
      ensure_dir(paste0(OutputPath,"/",subfolder,"/figures/"))
      malignantscore <- runMalignancy(PrimaryTumour_select, genelist, cell.annotation = cell.annotation, savePath = paste0(OutputPath, "/Primary/"))
      cell.annotation$CNVscore = malignantscore$cell.annotation$Malign.score
      cell.annotation$CNVtype = malignantscore$cell.annotation$Malign.type
      cells_cnv = intersect(cells,rownames(subset(cell.annotation,CNVtype=="malignant")))
      malignantscore <- predMalignantCell(PrimaryTumour_select,cell.annotation = cell.annotation,ModelPath = "scCancer2/scCancer-master/inst/txt/", MALIGNANT.THRES = 0.5)
      cell.annotation$XGBscore = malignantscore$Malign.score
      cell.annotation$XGBtype = malignantscore$Malign.type
      cells_XGBoost = intersect(cells,rownames(subset(cell.annotation,XGBtype=="malignant")))
      write.table(cell.annotation, paste0(OutputPath,"/",subfolder,"/maligant.txt"), quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
      cells = union(cells_cnv,cells_XGBoost)
      PrimaryTumour = subset(PrimaryTumour,cells=cells)
      sccluster(PrimaryTumour, paste0(OutputPath, "/Primary/Tumour/"), resolution = 0.9, PCnums = 50, DimNum = 50, perplexity = 30,doubletflag=FALSE) 
      scDEGs(paste0(OutputPath, "/Primary/Tumour/ScData.pca.rds"),OutputFile = paste0(OutputPath, "/Primary/Tumour/DEG.txt"),assay = "RNA",clustergroup = "seurat_clusters",logFC = 0)
    }
'

rm -r GSE277783/*_Pt-12*  #removed due to the limited tumor cells in primary
rm -r GSE277783/Pt-12/ #removed due to the limited tumor cells in primary

############################################
# 3.3 Metastatic tumor processing (CNV/XGBoost) of GSE277783 
#     and gene expression matrix construction
############################################


Rscript -e '
    library(Seurat);
    library(ggplot2);
    library(scCancer)
    Sys.setenv(OMP_NUM_THREADS = 1,MKL_NUM_THREADS = 1,OPENBLAS_NUM_THREADS = 1,VECLIB_MAXIMUM_THREADS = 1)
    set.seed(1234)
    source("../Processing.R");
    crc_dirs <- list.files(path = "GSE277783/",pattern = "[B-E]$",full.names = TRUE)
    for(i in 1:length(crc_dirs)){
      Num = substr(crc_dirs[i], nchar(crc_dirs[i]), nchar(crc_dirs[i]))
      subfolder = paste0("Metastasis",Num)
      OutputPath=paste0("GSE277783/",sub(".*_(Pt-\\\\d+).*", "\\\\1", crc_dirs[i])) #OutputPath=paste0("GSE277783/",sub(".*_(Pt-\\d+).*", "\\1", crc_dirs[i]))
      Metastasis <- readRDS(paste0(OutputPath,"/",subfolder,"/ScData.pca.rds"))
      sccluster(Metastasis, paste0(OutputPath,"/",subfolder,"/"), resolution = 0.9, PCnums = 50, DimNum = 50, perplexity = 30,doubletflag=FALSE)
      Metastasis <- readRDS(paste0(OutputPath,"/",subfolder,"/ScData.pca.rds"))
      refscore = read.table(paste0(OutputPath,"/",subfolder,"/refcnv.txt"),sep="\t",header=TRUE,check.names = F)
      Metastasis$Tumor = ifelse(Metastasis$TumorMarker1>0 & Metastasis$cnvscore>max(refscore) & Metastasis$TumorMarker1!="NA" & Metastasis$cnvscore!="NA","T","N")
      MetastasisTumour = subset(Metastasis,Tumor=="T")
      MetastasisTumour_select = MetastasisTumour
      cell.annotation <- data.frame(barcodes = colnames(MetastasisTumour_select), Cell.Type = "Epithelial")
      rownames(cell.annotation) <- colnames(MetastasisTumour_select)
      umap_coords <- Embeddings(MetastasisTumour_select, reduction = "umap")
      cell.annotation <- data.frame(cell.annotation, tSNE_1 = umap_coords[,"umap_1"], tSNE_2 = umap_coords[,"umap_2"])
      cell.annotation$Cluster=MetastasisTumour_select@meta.data$seurat_clusters
      MetastasisTumour_select@assays$RNA <- as(MetastasisTumour_select@assays$RNA, "Assay")
      cells = rownames(cell.annotation)
      genelist = read.table("GSE277783/genelist.txt",sep="\t",header=FALSE,check.names = F)
      colnames(genelist) = c("SYMBOL","EnsemblID")
      genelist <- genelist[!duplicated(genelist$SYMBOL), ]
      rownames(genelist)=genelist$SYMBOL
      genes = intersect(rownames(MetastasisTumour_select),rownames(genelist))
      genelist = genelist[genes,]
      MetastasisTumour_select = subset(MetastasisTumour_select, features = genes)
      ensure_dir(paste0(OutputPath,"/",subfolder,"/figures/"))
      malignantscore <- runMalignancy(MetastasisTumour_select, genelist, cell.annotation = cell.annotation, savePath = paste0(OutputPath,"/",subfolder,"/"))
      cell.annotation$CNVscore = malignantscore$cell.annotation$Malign.score
      cell.annotation$CNVtype = malignantscore$cell.annotation$Malign.type
      cells_cnv = intersect(cells,rownames(subset(cell.annotation,CNVtype=="malignant")))
      malignantscore <- predMalignantCell(MetastasisTumour_select,cell.annotation = cell.annotation,ModelPath = "scCancer2/scCancer-master/inst/txt/", MALIGNANT.THRES = 0.5)
      cell.annotation$XGBscore = malignantscore$Malign.score
      cell.annotation$XGBtype = malignantscore$Malign.type
      cells_XGBoost= intersect(cells,rownames(subset(cell.annotation,XGBtype=="malignant")))
      write.table(cell.annotation, paste0(OutputPath,"/",subfolder,"/maligant.txt"), quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)     
      cells = union(cells_cnv,cells_XGBoost)
      MetastasisTumour = subset(MetastasisTumour,cells=cells)
      PrimaryTumour <- readRDS(paste0(OutputPath, "/Primary/Tumour/ScData.pca.rds"))
      DEGData <- read.delim(paste0(OutputPath, "/Primary/Tumour/DEG.txt"), header = TRUE)
      SigDEGData <- subset(DEGData, p_val < 0.05 & avg_log2FC > 0)
      PrimaryCount <- as.matrix(PrimaryTumour@assays$RNA@layers$data)
      rownames(PrimaryCount) <- rownames(PrimaryTumour)
      colnames(PrimaryCount) <- colnames(PrimaryTumour)    
      MetastasisCount <- as.matrix(MetastasisTumour@assays$RNA@layers$data)
      rownames(MetastasisCount) <- rownames(MetastasisTumour)
      colnames(MetastasisCount) <- colnames(MetastasisTumour)
      genelist <- intersect(unique(SigDEGData$gene), intersect(rownames(PrimaryCount), rownames(MetastasisCount)))  
      PrimaryCount <- PrimaryCount[genelist, ]
      MetastasisCount <- MetastasisCount[genelist, ]
      if(!file.exists(paste0(OutputPath, "/Primary/PrimaryCount.csv"))){
        write.csv(PrimaryCount, paste0(OutputPath, "/Primary/PrimaryCount.csv"), quote = FALSE)  
      }
      write.csv(MetastasisCount, paste0(OutputPath,"/",subfolder,"/MetastasisCount.csv"), quote = FALSE)
    }
'

############################################
# 4. MIC identification by scMIC of GSE178318 
############################################

samples=('COL07' 'COL15' 'COL16')
for i in "${!samples[@]}"
do 
  sample="${samples[i]}"  
  python ../scMIC.py \
    --pri GSE178318/$sample/Primary/PrimaryCount.csv \
    --met GSE178318/$sample/Metastasis/MetastasisCount.csv \
    --out GSE178318/$sample/proximity_OT.csv \
    --ot
done

############################################
# 5. MIC identification by scMIC of GSE178318 
############################################

samples=("P2" "P3" "P4" "P5")
for i in "${!samples[@]}"
do 
  sample="${samples[i]}"  
  python ../scMIC.py \
    --pri OMIX002487/$sample/Primary/PrimaryCount.csv \
    --met OMIX002487/$sample/Metastasis/MetastasisCount.csv \
    --out OMIX002487/$sample/proximity_OT.csv \
    --ot
done

############################################
# 6. MIC identification by scMIC of GSE277783 
############################################

ls GSE277783/Pt-*/Metastasis*/MetastasisCount.csv | while read file
do
  PrimaryPath=`echo $file | awk -F"/" '{print $1"/"$2"/"}'`
  OutputPath=`echo $file | awk -F"/" '{print $1"/"$2"/"$3}'`
  python ../scMIC.py \
    --pri $PrimaryPath/Primary/PrimaryCount.csv \
    --met $file \
    --out $OutputPath/proximity_OT.csv \
    --ot
done