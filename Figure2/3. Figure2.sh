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
#    -Input:Normalized ScData.pca.rds during Figure2A
#           h.all.v2025.1.Hs.symbols.gmt from MSigDB
#           the genes from original data are transferred to human genes
#    -Output:EMT.txt
#            EMT scores
############################################
Rscript -e '
    library(Seurat);
    Sys.setenv(OMP_NUM_THREADS = 1,MKL_NUM_THREADS = 1,OPENBLAS_NUM_THREADS = 1,VECLIB_MAXIMUM_THREADS = 1)
    set.seed(1234) 
    AData = readRDS("GSE249057/GSM7925723_4mo/ScData.pca.rds");
    GeneSets=read.table("GSE249057/h.all.v2025.1.Hs.symbols.gmt",sep="\t",header=FALSE,fill=TRUE,check.names=FALSE);
    eachgeneset = GeneSets[14, ];  #EMT
    eachgeneset = as.character(eachgeneset[3:length(eachgeneset)]);
    eachgeneset = eachgeneset[which(eachgeneset != "")];
    genes = list("genes" = eachgeneset[eachgeneset %in% rownames(AData)]);
    AData[["joined"]] <- JoinLayers(AData[["RNA"]])
    eachgenesetscore = AddModuleScore(object = AData, features = genes, name = "EMT",assay="joined",slot = "data" );
    AEMT=data.frame("sample"=eachgenesetscore@meta.data$Sampleid,"group"=eachgenesetscore@meta.data$seurat_clusters,"EMT"=eachgenesetscore@meta.data$EMT1);
    write.table(AEMT, file = "GSE249057/GSM7925723_4mo/EMT.txt", sep = "\t", quote = FALSE, row.names = TRUE); 
'

############################################
# 3. Figure2C
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
# 4. Figure2D
#    -Input:Normalized ScData.pca.rds during Figure2C
#           h.all.v2025.1.Hs.symbols.gmt from MSigDB
#    -Output:EMT.txt:EMT scores
#            ttest comparison results
############################################

Rscript -e '
    library(Seurat);
    Sys.setenv(OMP_NUM_THREADS = 1,MKL_NUM_THREADS = 1,OPENBLAS_NUM_THREADS = 1,VECLIB_MAXIMUM_THREADS = 1)
    set.seed(1234)
    PData = readRDS("GSE249057/GSM7925719_0h/ScData.pca.rds");    
    GeneSets=read.table("GSE249057/h.all.v2025.1.Hs.symbols.gmt",sep="\t",header=FALSE,fill=TRUE,check.names=FALSE);#download from MSigDB
    eachgeneset = GeneSets[14, ];  #EMT
    eachgeneset = as.character(eachgeneset[3:length(eachgeneset)]);
    eachgeneset = eachgeneset[which(eachgeneset != "")];
    genes = list("genes" = eachgeneset[eachgeneset %in% rownames(PData)]); 
    eachgenesetscore = AddModuleScore(object = PData, features = genes, name = "EMT",assay="RNA",slot = "data");
    PEMT = eachgenesetscore@meta.data[,c("seurat_clusters","EMT1")];            
    write.table(PEMT, file = "GSE249057/GSM7925719_0h/EMT.txt", sep = "\t", quote = FALSE, row.names = TRUE); 
    for(i in 0:3){
      ttresult=t.test(PEMT$EMT1[which(PEMT$seurat_clusters==4)],PEMT$EMT1[which(PEMT$seurat_clusters==i)]);
      print(c(i,ttresult$p.value,as.double(ttresult$estimate)))  
    }
'


############################################
# 5. Figure2E
#    -Input:Normalized ScData.pca.rds during Figure2C
#           mh.all.v2025.1.Mm.symbols.gmt from MSigDB
#    -Output:cytotrace.score.txt:stemness scores
#            ttest comparison results
############################################

Rscript -e '
    library(CytoTRACE2);
    Sys.setenv(OMP_NUM_THREADS = 1,MKL_NUM_THREADS = 1,OPENBLAS_NUM_THREADS = 1,VECLIB_MAXIMUM_THREADS = 1)
    set.seed(1234)
    PData = readRDS("GSE249057/GSM7925719_0h/ScData.pca.rds"); 
    PrimaryMatrix=PData@assays$RNA@layers$data;
    rownames(PrimaryMatrix)=rownames(PData);
    colnames(PrimaryMatrix)=colnames(PData);
    results <- cytotrace2(as.data.frame(as.matrix(PrimaryMatrix)),species="human",is_seurat = FALSE,slot_type = "data");
    results=results[rownames(PData@meta.data),];
    OutputData=cbind(results,PData@meta.data);
    write.table(OutputData,"GSE249057/GSM7925719_0h/cytotrace.score.txt",sep="\t",quote=FALSE,row.names=FALSE,col.names=TRUE);
    for(i in 0:3){
      ttresult=t.test(OutputData$CytoTRACE2_Score[which(OutputData$seurat_clusters==4)],OutputData$CytoTRACE2_Score[which(OutputData$seurat_clusters==i)]);
      print(c(i,ttresult$p.value,as.double(ttresult$estimate)))
    }
'

############################################
# 6. Figure2F
#    -Input:Normalized ScData.pca.rds during Figure2C
#    -Output:marker.geneexpression.txt:marker expression
#            marker.geneexpression.statistics.txt: statistics results
############################################

Rscript -e '
    library(Seurat);
    PData = readRDS("GSE249057/GSM7925719_0h/ScData.pca.rds");
    Pmatrix = PData@assays$RNA@layers$data;
    rownames(Pmatrix) = rownames(PData);
    colnames(Pmatrix) = colnames(PData);
    Genematrix = data.frame("CD44"=as.double(Pmatrix["CD44",]),"CST6"=as.double(Pmatrix["CST6",]),"C19orf33"=as.double(Pmatrix["C19orf33",]),
                            "TACSTD2"=as.double(Pmatrix["TACSTD2",]),"S100A14"=as.double(Pmatrix["S100A14",]),"RHOD"=as.double(Pmatrix["RHOD",]),
                            "TM4SF1"=as.double(Pmatrix["TM4SF1",]),"CD24"=as.double(Pmatrix["CD24",]),"ALDH1A3"=as.double(Pmatrix["ALDH1A3",]),
                            "group"=as.character(PData@meta.data$seurat_clusters));
    index = which(Genematrix$group == 4);
    Otherindex = which(Genematrix$group != 4);
    for(i in 1:(dim(Genematrix)[2]-1)){
      ttest=t.test(Genematrix[index,i],Genematrix[-index,i]);
      print(c(colnames(Genematrix)[i],ttest$p.value,as.double(ttest$estimate)))
    }
    Genematrix_new=data.frame("Gene"=c(rep("CD44",length(Otherindex)),rep("CST6",length(Otherindex)),rep("C19orf33",length(Otherindex)),rep("TACSTD2",length(Otherindex)),
                                       rep("S100A14",length(Otherindex)),rep("RHOD",length(Otherindex)),rep("TM4SF1",length(Otherindex)),rep("CD24",length(Otherindex)),
                                       rep("ALDH1A3",length(Otherindex))),
                                       "NonMIC"=c(Genematrix$CD44[Otherindex],Genematrix$CST6[Otherindex],Genematrix$C19orf33[Otherindex],
                                       Genematrix$TACSTD2[Otherindex],Genematrix$S100A14[Otherindex],Genematrix$RHOD[Otherindex],Genematrix$TM4SF1[Otherindex],
                                       Genematrix$CD24[Otherindex],Genematrix$ALDH1A3[Otherindex]),
                                       "MIC"=c(Genematrix$CD44[index],rep("",length(Otherindex)-length(index)),
                                       Genematrix$CST6[index],rep("",length(Otherindex)-length(index)),
                                       Genematrix$C19orf33[index],rep("",length(Otherindex)-length(index)),
                                       Genematrix$TACSTD2[index],rep("",length(Otherindex)-length(index)),
                                       Genematrix$S100A14[index],rep("",length(Otherindex)-length(index)),
                                       Genematrix$RHOD[index],rep("",length(Otherindex)-length(index)),
                                       Genematrix$TM4SF1[index],rep("",length(Otherindex)-length(index)),
                                       Genematrix$CD24[index],rep("",length(Otherindex)-length(index)),
                                       Genematrix$ALDH1A3[index],rep("",length(Otherindex)-length(index))));
    write.table(Genematrix_new, file = "GSE249057/GSM7925719_0h/marker.geneexpression.txt", sep = "\t", quote = FALSE, row.names = TRUE);
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