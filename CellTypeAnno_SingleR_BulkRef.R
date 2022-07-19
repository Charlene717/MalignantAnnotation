## SingleRBook Ref: http://bioconductor.org/books/release/SingleRBook/
## Example Ref: https://bioconductor.org/packages/devel/bioc/vignettes/SingleR/inst/doc/SingleR.html

##### Presetting ######
  rm(list = ls()) # Clean variable
  memory.limit(300000)


##### Load Packages #####
  if(!require("Seurat")) install.packages("Seurat")
  if(!require("tidyverse")) install.packages("tidyverse")
  if(!require("ggpubr")) install.packages("ggpubr")

  library(ggpubr)
  library(tidyverse)
  library(Seurat)

#### BiocManager installation ####
  ## Check whether the installation of those packages is required from BiocManager
  if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  Package.set <- c("SingleR","scRNAseq","celldex","scran","scater","scuttle")
  for (i in 1:length(Package.set)) {
    if (!requireNamespace(Package.set[i], quietly = TRUE)){
      BiocManager::install(Package.set[i])
    }
  }
  ## Load Packages
  lapply(Package.set, library, character.only = TRUE)
  rm(Package.set,i)


#### Load data #####
  load("D:/Dropbox/##_GitHub/##_Charlene/TrajectoryAnalysis/SeuratObject_CDS_PRJCA001063_V2.RData")
  seuratObject_1 <- scRNA.SeuObj
  seuratObject_1@meta.data[["DataSetID"]] <- rep("PRJCA001063")
  rm(list=setdiff(ls(), "seuratObject_1"))

  #load("D:/Dropbox/##_GitHub/##_CAESAR/MagicDisc/2022-06-25_PDAC_GSE131886_SC/04_Perform_an_integrated_analysis.RData")
  load("D:/Dropbox/##_GitHub/##_CAESAR/MagicDisc/2022-06-23_PDAC_GSE154778_SC/04_Perform_an_integrated_analysis.RData")
  #### Re-dimension reduction ####
  DefaultAssay(scRNA.SeuObj) <- "RNA"
  # scRNA.SeuObj <- FindVariableFeatures(scRNA.SeuObj, selection.method = "vst", nfeatures = 2000)
  scRNA.SeuObj <- FindVariableFeatures(scRNA.SeuObj)
  # Run the standard workflow for visualization and clustering
  scRNA.SeuObj <- ScaleData(scRNA.SeuObj, verbose = FALSE)
  scRNA.SeuObj <- RunPCA(scRNA.SeuObj, npcs = 200, verbose = FALSE)
  # scRNA.SeuObj <- RunUMAP(scRNA.SeuObj, reduction = "pca", dims = 1:160,n.neighbors = 20,min.dist = 0.3)
  scRNA.SeuObj <- FindNeighbors(scRNA.SeuObj, reduction = "pca", dims = 1:200)
  scRNA.SeuObj <- FindClusters(scRNA.SeuObj, resolution = 0.1)

  ##### Plot #####
  scRNA.SeuObj <- RunUMAP(scRNA.SeuObj, reduction = "pca", dims = 1:100,n.neighbors = 30, min.dist = 0.3)
  FeaturePlot(scRNA.SeuObj, features = c("TOP2A"))
  DimPlot(scRNA.SeuObj, reduction = "umap")
  DimPlot(scRNA.SeuObj, reduction = "umap",label = T)
  # DimPlot(scRNA.SeuObj, reduction = "umap",group.by = "Cell_type")

  seuratObject_2 <- scRNA.SeuObj
  DefaultAssay(seuratObject_2) <- "RNA"
  rm(list=setdiff(ls(), str_subset(objects(), pattern = "seuratObject")))

##### Function setting #####
  ## Call function
  source("FUN_Anno_SingleR.R")


##### Current path and new folder setting* #####
  ProjectName = paste0("CTAnno_singleR_PRJCA001063SSingle")
  Sampletype = "PDAC"

  Version = paste0(Sys.Date(),"_",ProjectName,"_",Sampletype)
  Save.Path = paste0(getwd(),"/",Version)
  ## Create new folder
  if (!dir.exists(Save.Path)){
    dir.create(Save.Path)
  }

##### Parameter setting* #####
  Remark1 <- "PredbyCTDB" # c("PredbyCTDB","PredbyscRNA")
  RefType <- "BuiltIn_celldex" # c("BuiltIn_celldex","BuiltIn_scRNA")
  celldexDatabase <- "BlueprintEncodeData"
  # c("BlueprintEncodeData","DatabaseImmuneCellExpressionData","HumanPrimaryCellAtlasData","ImmGenData",
  #   "MonacoImmuneData","MouseRNAseqData","NovershternHematopoieticData")
  de.method <- "classic"

  ## Parameter of classifySingleR
  quantile = 0.8
  tune.thresh = 0.05
  sd.thresh = 1

  Remark <- paste0(Remark1,"_",de.method,"_",
                   "qua",quantile,"_tun",tune.thresh,"_sd",sd.thresh)

  SmallTest = F


##### Run singleR #####
#### Presetting ####

  ## SeuObj_Tar
  scRNA.SeuObj <- seuratObject_2
  scRNA.SeuObj@meta.data[["Cell_type"]] <- ""

  ## SeuObj_Ref
  scRNA.SeuObj_Ref <- seuratObject_1

  if(SmallTest == TRUE){
    ## SeuObj_Ref for small test
    # CTFeatures.SeuObj <- scRNA.SeuObj_Ref[,scRNA.SeuObj_Ref$CELL %in% sample(scRNA.SeuObj_Ref$CELL,1000)] ## For small test
    CTFeatures.SeuObj <- scRNA.SeuObj_Ref[,scRNA.SeuObj_Ref@meta.data[[1]] %in% sample(scRNA.SeuObj_Ref@meta.data[[1]],1000)] ## For small test
    ## SeuObj_Tar for small test
    # scRNA.SeuObj <- scRNA.SeuObj[,scRNA.SeuObj$CELL %in% sample(scRNA.SeuObj$CELL,1000)] ## For small test
    scRNA.SeuObj <- scRNA.SeuObj[,scRNA.SeuObj@meta.data[["cells"]] %in% sample(scRNA.SeuObj@meta.data[["cells"]],1000)] ## For small test
  }else{
    ## SeuObj_Ref for full data
    CTFeatures.SeuObj <- scRNA.SeuObj_Ref
  }


  SingleRResult.lt <- Anno_SingleR(scRNA.SeuObj, RefType = RefType, celldexDatabase = celldexDatabase,
                                   CTFeatures.SeuObj = CTFeatures.SeuObj,
                                   quantile = quantile, tune.thresh = tune.thresh, sd.thresh = sd.thresh,
                                   de.method = de.method,
                                   Remark = Remark, Save.Path = paste0(Save.Path,"/",Remark), ProjectName = "CT")

