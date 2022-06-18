##### Presetting ######
  rm(list = ls()) # Clean variable
  memory.limit(150000)

##### Current path and new folder setting* #####
  ProjectName = "MaliAnno"
  Sampletype = "PDAC"
  #ProjSamp.Path = paste0(Sampletype,"_",ProjectName)

  Version = paste0(Sys.Date(),"_",ProjectName,"_",Sampletype)
  Save.Path = paste0(getwd(),"/",Version)
  ## Create new folder
  if (!dir.exists(Save.Path)){
    dir.create(Save.Path)
  }

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
    Package.set <- c("fgsea","AnnotationHub","ensembldb",
                     "SeuratDisk","monocle",
                     "SingleR","scRNAseq","celldex","scran")
    for (i in 1:length(Package.set)) {
      if (!requireNamespace(Package.set[i], quietly = TRUE)){
        BiocManager::install(Package.set[i])
      }
    }
    ## Load Packages
    lapply(Package.set, library, character.only = TRUE)
    rm(Package.set,i)

    options(stringsAsFactors = FALSE)

  #### GitHub installation ####
    if (!require("devtools", quietly = TRUE))
      install.packages("devtools")
    library(monocle)
    devtools::install_github("cole-trapnell-lab/garnett")
    devtools::install_github('cole-trapnell-lab/monocle3')
    devtools::install_github("LTLA/SingleR")

    library(monocle3)
    library(garnett)
    # library(SingleR)

#### Data preprocessing #####
    load("SeuratObject_PRJCA001063.RData")
    load("D:/Dropbox/##_GitHub/##_PHH_Lab/#_H5AD_PRJCA001063_PDAC/#_20220525_CleanUpS.RData")

    seuratObject_Ori <- seuratObject
    seurat_meta.df <- seuratObject@meta.data
    cds_meta.df <- as.data.frame(cds@colData@listData)
    seurat_meta.df <- left_join(seurat_meta.df, cds_meta.df)
    #seurat_meta.df[is.na(seurat_meta.df)] <- 0
    row.names(seurat_meta.df) <- seurat_meta.df[,1]
    seuratObject@meta.data <- seurat_meta.df

    seuratObjectMono_Ori <- seuratObject

    library("stringr")
    rm(list=setdiff(ls(), str_subset(objects(), pattern = "seuratObject")))
    save.image("SeuratObject_CDS_PRJCA001063.RData")

# #### Load data #####
     load("SeuratObject_CDS_PRJCA001063.RData")

#### Plot UMAP #####
    FeaturePlot(seuratObject, features = c("MS4A1", "GNLY", "CD3E", "CD14"))
    Idents(seuratObject) <- seuratObject@meta.data[["Cell_type"]]
    Idents(seuratObject) <- seuratObject@meta.data[["ReCluster"]]
    DimPlot(seuratObject, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()

    # ReDR
    seuratObject <- FindVariableFeatures(seuratObject, selection.method = "vst", nfeatures = 2000)
    seuratObject <- FindVariableFeatures(seuratObject)
    seuratObject <- RunPCA(seuratObject,npcs = 200, features = VariableFeatures(object = seuratObject))
    seuratObject <- FindNeighbors(seuratObject, dims = 1:100)
    seuratObject <- FindClusters(seuratObject, resolution = 0.5)
    seuratObject <- RunUMAP(seuratObject, dims = 1:100,n.neighbors = 20, min.dist=0.3)
    # seuratObject <- RunUMAP(seuratObject, dims = 1:100,n.neighbors = 1000, min.dist=0.1)
    # seuratObject@meta.data[["UMAP_NNei1000"]] <- seuratObject@reductions[["umap"]]@cell.embeddings
    # seuratObject@meta.data <- seuratObject@meta.data[,!colnames(seuratObject@meta.data)=="UMAP_NNei1000"]
    seuratObject@meta.data[["UMAP_NNei20_MD03"]] <- seuratObject@reductions[["umap"]]@cell.embeddings

    DimPlot(seuratObject, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()

    Idents(seuratObject) <- seuratObject@meta.data[["Cell_type"]]
    DimPlot(seuratObject, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()

    Idents(seuratObject) <- seuratObject@meta.data[["ReCluster"]]
    DimPlot(seuratObject, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
    FeaturePlot(seuratObject, features = c("TOP2A"))

#### Cell type markers #####
    ## Create cell type markers dataframe
    # DefaultAssay(scRNA.SeuObj_Small) <- "RNA"
    Idents(seuratObject) <- seuratObject@meta.data[["ReCluster"]]
    ReCellType.markers <- FindAllMarkers(seuratObject, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.25)

    write.table(ReCellType.markers, file = paste0( Save.Path,"/ReCelltypeMarker2_AllGene.txt"),
                quote = F,sep = "\t",row.names = F)
    save.image("SeuratObject_CDS_PRJCA001063_MaligAnno.RData")






