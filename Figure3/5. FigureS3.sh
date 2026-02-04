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
# 1. FigureS3
#    -Input:Normalized ScData.pca.rds during DataProcessing.sh
#    -Output:umap with EMT markers
############################################

Rscript -e '
    library(Seurat); 
    library(ggplot2);
    Data = readRDS("GSE173958/M1results/Primary/ScData.pca.rds");
    genes = c("Epcam","Muc1","Cdh1","Sparc","Zeb1","Zeb2")      
    for(gene in genes){
      pdf(paste0("GSE173958/M1results/",gene,".pdf"));
      figure1=FeaturePlot(Data, features = gene, reduction = "umap",  pt.size = 1) +  theme_bw() + ggtitle(NULL) + 
              theme(
                legend.position = "none",
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                panel.border = element_rect(fill = NA, color = "black", linewidth = 1, linetype = "solid"),
                axis.text = element_text(size = 25),
                axis.title.x = element_text(vjust = 0.5, hjust = 0.5, size = 25),
                axis.title.y = element_text(vjust = 0.5, hjust = 0.5, size = 25)
              );
      print(figure1);
      dev.off();
    }
'