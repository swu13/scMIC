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
mkdir Figure3/
cd Figure3/
# create environment
conda env create -f 1.\ scMIC_Figure3.yml
conda activate scMIC_Figure3

############################################
# 0. Data download
############################################

mkdir -p GSE173958
cd GSE173958
wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE173nnn/GSE173958/suppl/GSE173958_RAW.tar -O GSE173958_RAW.tar
tar -xvf GSE173958_RAW.tar
wget "https://prod-dcd-datasets-cache-zipfiles.s3.eu-west-1.amazonaws.com/t98pjcd7t6-1.zip" -O Clone.zip #check https://data.mendeley.com/datasets/t98pjcd7t6/1 for the download link
unzip Clone.zip
unzip M1.zip
unzip M2.zip

############################################
# 0.1 Reorganize 10X-format files
############################################

ls *"barcodes.tsv.gz" | while read file
do
  id=`echo $file | awk -F"/" '{split($NF,a,"-");print a[1]"-"a[2]}'`
  mkdir $id"/"
  mv $file $id"/barcodes.tsv.gz"
  mv ${file/"barcodes"/"features"} $id"/features.tsv.gz"
  mv ${file/"barcodes.tsv"/"matrix.mtx"} $id"/matrix.mtx.gz"
done

rm -f *txt.gz

############################################
# 0.2 Prepare clone information and 
#     select the cells with clone
############################################

Rscript -e '
  library(Seurat);
  folders=list.dirs(".",full.names=FALSE);
  samples=folders[grepl("^GSM",folders)];
  samples = samples[1:10];
  for(sample in samples){
    id=strsplit(strsplit(sample,"_")[[1]][2],"-")[[1]];
    mouseid=id[1];
    groupid=id[2];
    Data <- Read10X(sample, gene.column = 2);
    clones = read.table(paste0(mouseid,"/output-files/pdac_mouse",substr(mouseid,2,2),"_final_classifications.txt"), 
                        sep = "\t", header = TRUE, check.names = FALSE);
    split_cols <- do.call(rbind, strsplit(clones$cellID, "_"));
    clones$Sample <- split_cols[,1];
    clones$Barcode <- paste0(split_cols[,2],"-1");
    clones$Sample[which(clones$Sample=="Blood")]="CTC";
    clones$Sample[which(clones$Sample=="PTab")]="SS";
    SelectedClones=subset(clones,Sample==groupid & type=="final_singlet");
    Data=Data[,intersect(SelectedClones$Barcode,colnames(Data))];
    DropletUtils:::write10xCounts(paste0(sample,"-clones/"),Data,version="3");
    rownames(SelectedClones) = SelectedClones$Barcode;
    SelectedClones = SelectedClones[colnames(Data),];
    write.table(SelectedClones, paste0(sample,"-clones/clones.txt"), quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE);
  }
'

cd ../

############################################
# 1. Primary and metastatic tumor
#    gene expression preparation
############################################

############################################
# 1.1 Primary tumor processing
#     - QC, normalization
#     - clustering and DEG identification
############################################

Rscript -e '
  source("../Processing.R");
  Sys.setenv(OMP_NUM_THREADS = 1,MKL_NUM_THREADS = 1,OPENBLAS_NUM_THREADS = 1,VECLIB_MAXIMUM_THREADS = 1)
  set.seed(1234)
  PrimaryPaths = c("GSE173958/GSM5283482_M1-PT-clones/","GSE173958/GSM5283488_M2-PT-clones/")
  OutputPaths = c("GSE173958/M1results/","GSE173958/M2results/")
  for(i in 1:length(PrimaryPaths)){
    PrimaryPath=PrimaryPaths[i]
    PrimaryData <- scRead(PrimaryPath, Label = "Primary",RNAfeature=200,minRNAcount=0,maxRNAcount=2000000,mt=15,rb=100,hb=100);
    PrimaryNormalizedData <- scNormalize(PrimaryData, NormalizeMethod = "LogNormalize", ScaleFactor = 10000)
    OutputPath=OutputPaths[i]
    ensure_dir(OutputPath)
    sccluster(PrimaryNormalizedData, paste0(OutputPath, "/Primary/"), resolution = 0.9, PCnums = 50, DimNum = 50, perplexity = 30,doubletflag=FALSE)
    scRNAanno(paste0(OutputPath, "/Primary/ScData.pca.rds"), CellMarkerFile="", paste0(OutputPath, "/Primary/"))
    scDEGs(paste0(OutputPath, "/Primary/ScData.pca.rds"),OutputFile = paste0(OutputPath, "/Primary/DEG.txt"),assay = "RNA",clustergroup = "seurat_clusters",logFC = 0)
  }
'

############################################
# 1.2 Metastatic tumor processing and
#     gene expression matrix construction
############################################

Rscript -e '
  source("../Processing.R");
  Sys.setenv(OMP_NUM_THREADS = 1,MKL_NUM_THREADS = 1,OPENBLAS_NUM_THREADS = 1,VECLIB_MAXIMUM_THREADS = 1)
  set.seed(1234)
  MetastasisPaths=c("GSE173958/GSM5283484_M1-Met-clones/","GSE173958/GSM5283485_M1-Liver-clones/","GSE173958/GSM5283486_M1-Lung-clones/",
                    "GSE173958/GSM5283489_M2-Met-clones/","GSE173958/GSM5283490_M2-Liver-clones/","GSE173958/GSM5283491_M2-Lung-clones/")
  samples=sub(".*_(M[0-9]+)-.*", "\\\\1", MetastasisPaths)  #samples=sub(".*_(M[0-9]+)-.*", "\\1", MetastasisPaths)
  organs=sub(".*_M[0-9]+-([A-Za-z]+)-.*", "\\\\1", MetastasisPaths)  #organs=sub(".*_M[0-9]+-([A-Za-z]+)-.*", "\\1", MetastasisPaths) 
  for(i in 1:length(MetastasisPaths)){
    MetastasisPath=MetastasisPaths[i]
    sample=samples[i]
    organ=organs[i]
    MetastasisData <- scRead(MetastasisPath, Label = organ,RNAfeature=200,minRNAcount=0,maxRNAcount=2000000,mt=15,rb=100,hb=100)
    MetastasisNormalizedData <- scNormalize(MetastasisData, NormalizeMethod = "LogNormalize", ScaleFactor = 10000)  
    MetastasisTumour <- MetastasisNormalizedData
    OutputPath=paste0("GSE173958/",sample,"results/");
    PrimaryTumour <- readRDS(paste0(OutputPath, "/Primary/ScData.pca.rds"))
    DEGData <- read.delim(paste0(OutputPath, "/Primary/DEG.txt"), header = TRUE)
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
    ensure_dir(paste0(OutputPath, "/Metastasis_",organ,"/"))
    if(i %% 3 == 1){
      write.csv(PrimaryCount, paste0(OutputPath, "/Primary/PrimaryCount.csv"), quote = FALSE)  
    }
    write.csv(MetastasisCount, paste0(OutputPath, "/Metastasis_",organ,"/MetastasisCount.csv"), quote = FALSE)
  }
'

############################################
# 2. MIC identification by scMIC
############################################
OutputPaths=("M1results" "M1results" "M1results" "M2results" "M2results" "M2results")
Organs=("Met" "Liver" "Lung" "Met" "Liver" "Lung")
for ((i=0; i<${#OutputPaths[@]}; i++))
do
  OutputPath=${OutputPaths[$i]}
  Organ=${Organs[$i]}
  python ../scMIC.py \
    --pri GSE173958/$OutputPath/Primary/PrimaryCount.csv \
    --met GSE173958/$OutputPath/Metastasis"_"$Organ/MetastasisCount.csv \
    --out GSE173958/$OutputPath/Metastasis"_"$Organ/proximity_OT.csv \
    --ot
done