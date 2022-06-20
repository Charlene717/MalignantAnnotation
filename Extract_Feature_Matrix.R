##### Presetting ######
  rm(list = ls()) # Clean variable
  memory.limit(150000)

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



##### Load Seurat dataset #####
  load("06_Cell_type_annotation.RData")

##### Find Cell type marker #####
  Idents(scRNA.SeuObj) <- scRNA.SeuObj$celltype
  CellType.markers <- FindAllMarkers(scRNA.SeuObj, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.25)

  write.table(CellType.markers, file = paste0(PathCellType,"/CC_CelltypeMarker_AllGene.txt"),
              quote = F,sep = "\t",row.names = F)

  #### Filter marker ####
  ## By Rank
  TOPN <- 10
  CTTop.markers <-  CellType.markers %>%
                    # group_by(.,.[,Cluster]) %>% ## Bug
                    # top_n(.,n = TOPN, wt = .[,log2FC_CN]) ## Error
                    group_by(cluster) %>%
                    top_n(.,n = TOPN, wt = avg_log2FC) ## Bug

  ## By filters
  log2FC_CN <- "avg_log2FC"
  log2FC_Thr <- 1
  pValue_CN <- "p_val"
  pValue_Thr <- 0.05
  pValueAdj_CN <- "p_val_adj"
  pValueAdj_Thr <- 0.05

  CTFilter.markers <-  CellType.markers %>%
                       dplyr::filter(.,.[,log2FC_CN] > log2FC_Thr) %>%
                       dplyr::filter(.,.[pValue_CN] < pValue_Thr) %>%
                       dplyr::filter(.,.[,pValueAdj_CN] < pValueAdj_Thr)

##### Export Gene expression matrix #####
  GeneExp.df <- GetAssayData(scRNA.SeuObj, assay = "RNA", slot = "data") # normalized data matrix

##### For scSorter #####

  ##### Section 1 - Preliminaries #####
  library(scSorter)

  ##### Section 2 - Preprocessing the data #####
  ## Create anno.df


  #DoHeatmap(scRNA.SeuObj, features = CTFilter.markers$gene) + NoLegend()
  log2FC_CN <- "avg_log2FC"
  Gene_CN <- "gene"
  Cluster_CN <- "cluster"

  anno.df <- data.frame(Type = CTFilter.markers[,Cluster_CN],
                        Marker = CTFilter.markers[,Gene_CN],
                        Weight = CTFilter.markers[,log2FC_CN])

  # anno.df <- data.frame(Type = CTFilter.markers$cluster,
  #                       Marker = CTFilter.markers$gene,
  #                       Weight = CTFilter.markers$avg_log2FC)

  ##### Section 3 - Running scSorter #####
  scSorter.obj <- scSorter(GeneExp.df, anno.df)

  ## Viewing Results
  print(table(scSorter.obj$Pred_Type))

  mis_rate = 1 - mean(scSorter.obj$Pred_Type == true_type)
  round(mis_rate, 4)

  table(true_type, scSorter.obj$Pred_Type)

  ## Insert the scSOrter data to scRNA.SeuObj
  scRNA.SeuObj_Small$scSorterPred <- scSorter.obj[["Pred_Type"]]
  DimPlot(scRNA.SeuObj_Small, reduction = "umap", group.by ="scSorterPred" ,label = TRUE, pt.size = 0.5) + NoLegend()
  DimPlot(scRNA.SeuObj_Small, reduction = "umap", group.by ="celltype" ,label = TRUE, pt.size = 0.5) + NoLegend()





