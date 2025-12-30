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

# create environment
conda env create -f 1.\ scMIC_Figure2.yml

############################################
# 0. Data download
############################################

mkdir Figure2/
cd Figure2/
mkdir -p GSE249057
cd GSE249057
wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE249nnn/GSE249057/suppl/GSE249057_RAW.tar
tar -xvf GSE249057_RAW.tar

############################################
# 0.1 Reorganize 10X-format files
############################################

for file in *_barcodes.tsv.gz; do
    sampleID=$(basename "${file}" | awk -F"_" '{print $1"_"$2}')
    mkdir -p "${sampleID}"
    mv "${file}" "${sampleID}/barcodes.tsv.gz"
    mv "${file/_barcodes/_features}" "${sampleID}/features.tsv.gz"
    mv "${file/_barcodes.tsv/_matrix.mtx}" "${sampleID}/matrix.mtx.gz"
done

cd ../
echo "Step 0 completed: GSE249057 data prepared."


############################################
# 1. Primary and metastatic tumor
#    gene expression preparation
############################################

############################################
# 1.1 Primary tumor processing
#     - QC, normalization
#     - clustering and DEG identification
#     or using HVGs for scMIC
############################################

#Primary tumor samples clusteringa and DEGs calculation/ or using HVGs
Rscript -e '
  source("../Processing.R");
  PrimaryPath = "GSE249057/GSM7925719_0h/"  
  PrimaryData <- scRead(PrimaryPath, Label = "Primary",RNAfeature=200,minRNAcount=1000,maxRNAcount=8000,mt=15,rb=100,hb=100);
  PrimaryNormalizedData <- scNormalize(PrimaryData, NormalizeMethod = "LogNormalize", ScaleFactor = 10000)
  OutputPath="GSE249057/Results/";
  ensure_dir(OutputPath)
  sccluster(PrimaryNormalizedData, paste0(OutputPath, "/Primary/"), resolution = 0.9, PCnums = 50, DimNum = 50, perplexity = 30,doubletflag=FALSE)
  scRNAanno(paste0(OutputPath, "/Primary/ScData.pca.rds"), CellMarkerFile="", paste0(OutputPath, "/Primary/"))
  scDEGs(paste0(OutputPath, "/Primary/ScData.pca.rds"),OutputFile = paste0(OutputPath, "/Primary/DEG.txt"),assay = "RNA",clustergroup = "seurat_clusters",logFC = 0)
'

############################################
# 1.2 Metastatic tumor processing and
#     gene expression matrix construction
############################################

Rscript -e '
  source("../Processing.R");
  MetastasisPath="GSE249057/GSM7925723_4mo/"
  MetastasisData <- scRead(MetastasisPath, Label = "Metastasis",RNAfeature=200,minRNAcount=1000,maxRNAcount=8000,mt=15,rb=100,hb=100)
  MetastasisNormalizedData <- scNormalize(MetastasisData, NormalizeMethod = "LogNormalize", ScaleFactor = 10000)  
  MetastasisTumour <- MetastasisNormalizedData
  OutputPath="GSE249057/Results/";
  PrimaryTumour <- readRDS(paste0(OutputPath, "/Primary/ScData.pca.rds"))
  DEGData <- read.delim(paste0(OutputPath, "/Primary/DEG.txt"), header = TRUE)
  SigDEGData <- subset(DEGData, p_val < 0.05 & avg_log2FC > 0)
  PrimaryCount <- as.matrix(PrimaryTumour@assays$RNA@layers$data)
  rownames(PrimaryCount) <- rownames(PrimaryTumour)
  colnames(PrimaryCount) <- colnames(PrimaryTumour)    
  MetastasisCount <- as.matrix(MetastasisTumour@assays$RNA@layers$data)
  rownames(MetastasisCount) <- rownames(MetastasisTumour)
  colnames(MetastasisCount) <- colnames(MetastasisTumour)
  #DEG
  genelist <- intersect(unique(SigDEGData$gene), intersect(rownames(PrimaryCount), rownames(MetastasisCount)))  
  PrimaryCount1 <- PrimaryCount[genelist, ]
  MetastasisCount1 <- MetastasisCount[genelist, ]
  ensure_dir(paste0(OutputPath, "/Metastasis/"))
  write.csv(PrimaryCount1, paste0(OutputPath, "/Primary/PrimaryCount.csv"), quote = FALSE)  
  write.csv(MetastasisCount1, paste0(OutputPath, "/Metastasis/MetastasisCount.csv"), quote = FALSE)
  #HVG
  genelist <- intersect(VariableFeatures(PrimaryTumour),rownames(MetastasisCount))    
  PrimaryCount2 <- PrimaryCount[genelist, ]
  MetastasisCount2 <- MetastasisCount[genelist, ]
  write.csv(PrimaryCount2, paste0(OutputPath, '/Primary/PrimaryCountHVG.csv'), quote = FALSE)
  write.csv(MetastasisCount2, paste0(OutputPath, '/Metastasis/MetastasisCountHVG.csv'), quote = FALSE)
'

############################################
# 2. MIC identification by scMIC
############################################

python ../scMIC.py \
  --pri GSE249057/Results/Primary/PrimaryCount.csv \
  --met GSE249057/Results/Metastasis/MetastasisCount.csv \
  --out GSE249057/Results/proximity_OT.csv \
  --ot