#!/usr/bin/env Rscript

ensure_dir <- function(path) {
  if (!dir.exists(path)) dir.create(path, recursive = TRUE)
}

scRead <- function(Path, Label = "Primary", sc10XGeneColumn = 2, Format = "10X", project = "Project",
                   RNAfeature = 200, minRNAcount = 1000, maxRNAcount = 8000, mt = 15, rb=100, hb = 100) {
  library(Seurat)

  if (Format == "10X") {
    Data <- Read10X(Path, gene.column = sc10XGeneColumn)
    scData <- CreateSeuratObject(Data, project = project)
  } else if (Format == "rds") {
    scData <- readRDS(Path)
  } else {
    stop("Unsupported Format. Please use '10X' or 'rds'.")
  }

  scData[["percent.mt"]] <- PercentageFeatureSet(scData, pattern = "^(?i)mt-")
  scData@meta.data$percent.mt[is.na(scData@meta.data$percent.mt)] <- 0

  scData[["percent.rb"]] <- PercentageFeatureSet(scData, pattern = "^RP[SL]")
  scData@meta.data$percent.rb[is.na(scData@meta.data$percent.rb)] <- 0

  HBgenes <- c("HBA1", "HBA2", "HBB", "HBD", "HBE1", "HBG1", "HBG2", "HBM", "HBQ1", "HBZ", "Hba1", "Hba2", "Hbb", "Hbd", "Hbe1", "Hbg1", "Hbg2", "Hbm", "Hbq1", "Hbz")
  matched_genes <- CaseMatch(HBgenes, rownames(scData))
  scData[["percent.HB"]] <- PercentageFeatureSet(scData, features = matched_genes)
  scData@meta.data$percent.HB[is.na(scData@meta.data$percent.HB)] <- 0

  scData@meta.data$Sampleid <- Label
  scData@meta.data$Groupid <- Label

  scData <- subset(scData, subset = nFeature_RNA > RNAfeature &
                                 nCount_RNA > minRNAcount & nCount_RNA < maxRNAcount &
                                 percent.rb < rb & percent.mt < mt & percent.HB < hb)
  return(scData)
}

scNormalize <- function(Data, NormalizeMethod = "LogNormalize", ScaleFactor = 10000, varfeatures = 2000, regressvar = NULL,seed = 1234) {
  library(Seurat)
  set.seed(seed)
  
  is_list <- is.list(Data) && all(sapply(Data, function(x) inherits(x, "Seurat")))
  if (!is_list) {
    Data <- list(Data)
  }
  
  NormalizeDataList <- lapply(Data, function(obj) {
    if (NormalizeMethod == "SCTtransform") {
      obj <- SCTransform(obj, assay = "RNA", variable.features.n = varfeatures,
                         vars.to.regress = regressvar, return.only.var.genes = TRUE, vst.flavor = "v2")
    } else if (NormalizeMethod %in% c("LogNormalize", "CLR", "RC")) {
      obj <- NormalizeData(obj, assay = "RNA", normalization.method = NormalizeMethod, scale.factor = ScaleFactor)
      obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = varfeatures)
    } else {
      stop("Please input one of the normalize methods: 'LogNormalize' (default), 'CLR', 'RC', or 'SCTtransform'.")
    }
    return(obj)
  })
  
  if (length(NormalizeDataList) == 1) {
    return(NormalizeDataList[[1]])
  } else {
    return(NormalizeDataList)
  }
}


sccluster <- function(Data, OutFold, varfeatures = 2000, PCnums = 50, DimNum = 50,resolution = 0.9,
                      plotreduction = "umap", otherdraw = "Sampleid", perplexity=30,doubletflag=FALSE,seed = 1234) {
  library(Seurat)
  library(ggplot2)
  
  ensure_dir(OutFold)
  set.seed(seed)
  Data <- ScaleData(object = Data)
  Data <- RunPCA(object = Data, features = VariableFeatures(Data), npcs = PCnums, approx = FALSE)
  Data <- FindNeighbors(Data, reduction = "pca", dims = 1:DimNum)
  Data <- FindClusters(Data, resolution = resolution)
  Data <- RunTSNE(Data, reduction = "pca", dims = 1:DimNum, perplexity = min(30, perplexity), seed.use = seed)
  Data <- RunUMAP(Data, reduction = "pca", dims = 1:DimNum, seed.use = seed)
  
  if(doubletflag){
    library(DoubletFinder)
    sweep.res <- paramSweep(Data, PCs = 1:20, sct = FALSE)
    sweep.stats <- summarizeSweep(sweep.res, GT = FALSE)
    bcmvn <- find.pK(sweep.stats)
    chosen.pK <- as.numeric(as.character(bcmvn$pK[which.max(bcmvn$BCmetric)]))
    annotations <- Data@meta.data$seurat_clusters
    homotypic.prop <- modelHomotypic(annotations)
    nExp_poi <- round(0.075 * nrow(Data@meta.data))
    nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
    Data <- doubletFinder(Data,PCs = 1:20,pN = 0.25,pK = chosen.pK,nExp = nExp_poi.adj)
    DFname=grep("DF.classifications", colnames(Data@meta.data), value = TRUE)
    singlet.cells <- rownames(Data@meta.data[Data@meta.data[[DFname]] == "Singlet", ])
    Data <- subset(Data, cells = singlet.cells)
  }
  
  cluster_plot_file <- file.path(OutFold, paste0("cluster.", plotreduction, ".pdf"))
  pdf(cluster_plot_file)
  tryCatch({
    print(
      DimPlot(Data, reduction = plotreduction, label = TRUE, pt.size = 1, label.size = 5) +
        theme_bw() +
        theme(
          legend.position = 'none',
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_rect(fill = NA, color = "black", linewidth = 1, linetype="solid"),
          axis.text = element_text(size = 25),
          axis.title.x = element_text(vjust = 0.5, hjust = 0.5, size = 25),
          axis.title.y = element_text(vjust = 0.5, hjust = 0.5, size = 25)
        )
    )
  }, finally = {
    dev.off()
  })
  
  saveRDS(Data,file=paste0(OutFold,"/ScData.pca.rds"))
}


scRNAanno <- function(scFile,CellMarkerFile,OutputFold){
  #install.packages("lemon")
	library(Seurat)
	library(ggplot2)
  
  ensure_dir(OutputFold)
  Data <- readRDS(scFile)
  if(CellMarkerFile!=""){
    cell.markers <- read.table(CellMarkerFile, sep = "\t", header = FALSE, check.names = FALSE)
    colnames(cell.markers)[1] <- "geneSymbol"    
    genes=unique(cell.markers$geneSymbol[cell.markers$geneSymbol %in% rownames(Data)])
    Data <- readRDS(scFile)
    
    if (length(genes) > 0) {
      pdf(paste0(OutputFold, "/CellmarkersDotPlot.pdf"),width=length(genes)*0.3,height=max(4,length(unique(Data@meta.data$seurat_clusters)) * 0.3))
      #png(paste0(OutputFold, "/CellmarkersDotPlot.png"),
      #    width = length(genes) * 20,
      #    height = length(unique(Data@meta.data$seurat_clusters)) * 30)
      plot <- DotPlot(Data, features = genes, group.by = "seurat_clusters") +
        scale_color_continuous(low = 'green', high = 'red') +
        theme(legend.position = 'right',
              legend.text = element_text(hjust = 0.5, size = 12),
              legend.title = element_text(hjust = 0.5, size = 12),
              panel.border = element_rect(fill = NA, color = "black", size = 1, linetype = "solid"),
              axis.text = element_text(hjust = 0.5, size = 12),
              axis.title.x = element_blank(),
              axis.title.y = element_blank())+ RotatedAxis()
      print(plot)
      dev.off()
    }else {
      warning("No marker genes found in dataset.")
    }
  }
  
  umap_embeddings <- Embeddings(Data, "umap")
  tsne_embeddings <- Embeddings(Data, "tsne")
  
  plotdata <- data.frame(tSNE_1 = tsne_embeddings[,1],
                         tSNE_2 = tsne_embeddings[,2],
                         UMAP_1 = umap_embeddings[,1],
                         UMAP_2 = umap_embeddings[,2],
                         Cluster = Data@meta.data$seurat_clusters,
                         Group = Data@meta.data$Groupid)

  png(paste0(OutputFold, "/Cluster.png"),width=1000,height=1000)
  umap_plot <- ggplot(plotdata, aes(x = UMAP_1, y = UMAP_2, color = Cluster)) +
    geom_point(size = 1, alpha = 1) +
    guides(colour = guide_legend(override.aes = list(size = 8))) +
    theme_bw(base_size = 30) +
    theme(axis.title.x = element_text(size = 40),
          axis.title.y = element_text(size = 40),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.text.x = element_text(size = 40, color = "black"),
          axis.text.y = element_text(size = 40, color = "black"),
          legend.position = "bottom",
          legend.title = element_text(size = 40),
          legend.text = element_text(size = 40))
  print(umap_plot)
  dev.off()
}

scDEGs <- function(InputFile, OutputFile, assay = "RNA", logFC = 1, method = "wilcox", clustergroup = "seurat_clusters") {
  library(Seurat)
  
  if (!file.exists(InputFile)) {
    stop("Input file does not exist: ", InputFile)
  }
  Data <- readRDS(InputFile)
  
  if (!(assay %in% names(Data@assays))) {
    stop("Assay ", assay, " not found in Seurat object.")
  }
  DefaultAssay(Data) <- assay
  
  if (!clustergroup %in% colnames(Data@meta.data)) {
    stop("Cluster group '", clustergroup, "' not found in meta.data.")
  }
  
  Idents(Data) <- Data@meta.data[[clustergroup]]
  DEG <- FindAllMarkers(Data, logfc.threshold = logFC, test.use = method, only.pos = FALSE)
  
  output_dir <- dirname(OutputFile)
  ensure_dir(output_dir)
  
  write.table(DEG, OutputFile, quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
}

predMalignantCell <- function(expr,cell.annotation,ModelPath,MALIGNANT.THRES = 0.5){
  #mkdir Tool/scCancer2/
  #wget -P Tool/scCancer2/ https://github.com/czythu/scCancer/archive/refs/heads/master.zip
  #unzip Tool/scCancer2/master.zip -d Tool/scCancer2/
  #install.packages("xgboost")
  library(xgboost)
  model.path <- paste0(ModelPath, "/sc_xgboost.model")
  genes.path <- paste0(ModelPath, "/genes-scRNA-tcga-sorted.txt")
  model.ref <- xgb.load(model.path)
  features <- as.list(read.table(genes.path))[[1]]
  testdata <- t(as.matrix(expr@assays$RNA@scale.data))  
  aligntestdata <- matrix(data = 0, nrow = nrow(testdata), ncol = length(features),dimnames = list(rownames(testdata), features))
  current.features <- colnames(testdata)
  for(j in 1:length(features)){
      if(features[j] %in% current.features){
          aligntestdata[,j] <- testdata[, features[j]]
      }
  }
  testdata <- xgb.DMatrix(aligntestdata)
  predict.label <- predict(model.ref, testdata)
  cell.annotation$Malign.score <- predict.label
  predict.label[which(predict.label > MALIGNANT.THRES)] <- "malignant"
  predict.label[which(predict.label <= MALIGNANT.THRES)] <- "nonMalignant"
  cell.annotation$Malign.type <- predict.label
  return(cell.annotation)
}


