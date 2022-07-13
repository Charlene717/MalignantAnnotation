##### Presetting ######
  rm(list = ls()) # Clean variable
  memory.limit(300000)

##### Load Packages #####
  #### Basic installation ####
  Package.set <- c("tidyverse","Seurat","ggplot2","ggpmisc",
                   "stringr","magrittr","dplyr")
  ## Check whether the installation of those packages is required from basic
  for (i in 1:length(Package.set)) {
    if (!requireNamespace(Package.set[i], quietly = TRUE)){
      install.packages(Package.set[i])
    }
  }
  ## Load Packages
  lapply(Package.set, library, character.only = TRUE)
  rm(Package.set,i)

  #### BiocManager installation ####
  ## Check whether the installation of those packages is required from BiocManager
  if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  Package.set <- c("AnnotationHub", "SeuratDisk","monocle",
                   "SingleR")
  for (i in 1:length(Package.set)) {
    if (!requireNamespace(Package.set[i], quietly = TRUE)){
      BiocManager::install(Package.set[i])
    }
  }
  ## Load Packages
  lapply(Package.set, library, character.only = TRUE)
  rm(Package.set,i)


##### Load Seurat dataset #####
  #### PRJCA001063 ####
  load("SeuratObject_CDS_PRJCA001063.RData")

  # Test file
  DimPlot(scRNA.SeuObj, reduction = "umap")
  DimPlot(scRNA.SeuObj, reduction = "umap",group.by = "Cell_type")
  FeaturePlot(scRNA.SeuObj, features = c("TOP2A"))
  DimPlot(scRNA.SeuObj, reduction = "umap",group.by = "RNA_snn_res.0.5")

  # Add
  scRNA.SeuObj@meta.data[["DataSetID"]] <- rep("PRJCA001063")
  scRNA.SeuObj@meta.data[["celltype"]] <- scRNA.SeuObj@meta.data[["Cell_type"]]
  DimPlot(scRNA.SeuObj, reduction = "umap",group.by = "celltype")


  # Clean
  rm(list=setdiff(ls(), "scRNA.SeuObj"))
  # Save
  save.image("scRNADataset/SeuratObject_CDS_PRJCA001063.RData")

  #### GSE131886 ####
  load("D:/Dropbox/##_GitHub/##_CAESAR/MagicDisc/2022-06-25_PDAC_GSE131886_SC/04_Perform_an_integrated_analysis.RData")

  # Test file
  DimPlot(scRNA.SeuObj, reduction = "umap")
  FeaturePlot(scRNA.SeuObj, features = c("TOP2A"))
  DimPlot(scRNA.SeuObj, reduction = "umap",group.by = "integrated_snn_res.0.5")

  # Clean
  rm(list=setdiff(ls(), "scRNA.SeuObj"))
  # Save
  save.image("scRNADataset/SeuratObject_GSE131886_PDC_SC.RData")

  #### GSE154778 ####
  load("D:/Dropbox/##_GitHub/##_CAESAR/MagicDisc/2022-06-23_PDAC_GSE154778_SC/04_Perform_an_integrated_analysis.RData")
  # Test file
  DimPlot(scRNA.SeuObj, reduction = "umap")
  FeaturePlot(scRNA.SeuObj, features = c("TOP2A"))
  DimPlot(scRNA.SeuObj, reduction = "umap",group.by = "integrated_snn_res.0.5")

  # Run the standard workflow for visualization and clustering
  scRNA.SeuObj <- NormalizeData(scRNA.SeuObj)
  scRNA.SeuObj <- FindVariableFeatures(scRNA.SeuObj, selection.method = "vst", nfeatures = 2000)
  scRNA.SeuObj <- ScaleData(scRNA.SeuObj, verbose = FALSE)
  scRNA.SeuObj <- RunPCA(scRNA.SeuObj, npcs = 200, verbose = FALSE)
  scRNA.SeuObj <- RunUMAP(scRNA.SeuObj, reduction = "pca", dims = 1:160,n.neighbors = 20,min.dist = 0.3)
  scRNA.SeuObj <- FindNeighbors(scRNA.SeuObj, reduction = "pca", dims = 1:200)
  scRNA.SeuObj <- FindClusters(scRNA.SeuObj, resolution = 0.5)

  # Test file
  DimPlot(scRNA.SeuObj, reduction = "umap")
  FeaturePlot(scRNA.SeuObj, features = c("TOP2A"))
  DimPlot(scRNA.SeuObj, reduction = "umap",group.by = "integrated_snn_res.0.5")

  # Clean
  rm(list=setdiff(ls(), "scRNA.SeuObj"))
  # Save
  save.image("scRNADataset/SeuratObject_GSE154778_PDAC_SC.RData")


