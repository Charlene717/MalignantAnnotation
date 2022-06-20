## Ref: https://cran.r-project.org/web/packages/scSorter/vignettes/scSorter.html
Anno_scSorter <- function(variables) {

}



##### Section 1 - Preliminaries #####
## Instal and Load Packages
  #### Basic installation
  Package.set <- c("tidyverse","scSorter","Seurat","stringr","magrittr","dplyr")
  ## Check whether the installation of those packages is required from basic
  for (i in 1:length(Package.set)) {
    if (!requireNamespace(Package.set[i], quietly = TRUE)){
      install.packages(Package.set[i])
    }
  }
  ## Load Packages
  lapply(Package.set, library, character.only = TRUE)
  rm(Package.set,i)



##### Section 2 - Preprocessing the data #####
  ## Create anno.df
  # DoHeatmap(scRNA.SeuObj, features = CTFilter.Markers.df$gene) + NoLegend()
  log2FC_CN <- "avg_log2FC"
  Gene_CN <- "gene"
  Cluster_CN <- "cluster"

  anno.df <- data.frame(Type = CTFilter.Markers.df[,Cluster_CN],
                        Marker = CTFilter.Markers.df[,Gene_CN],
                        Weight = CTFilter.Markers.df[,log2FC_CN])

  # anno.df <- data.frame(Type = CTFilter.Markers.df$cluster,
  #                       Marker = CTFilter.Markers.df$gene,
  #                       Weight = CTFilter.Markers.df$avg_log2FC)

  ## Export Gene expression matrix
  GeneExp.df <- GetAssayData(scRNA.SeuObj, assay = "RNA", slot = "data") # normalized data matrix

##### Section 3 - Running scSorter #####
  scSorter.obj <- scSorter(GeneExp.df, anno.df)

  ## Insert the scSOrter data to scRNA.SeuObj
  # scRNA.SeuObj$scSorterPred <- scSorter.obj[["Pred_Type"]]
  # DimPlot(scRNA.SeuObj, reduction = "umap", group.by ="scSorterPred" ,label = TRUE, pt.size = 0.5) + NoLegend()
  # DimPlot(scRNA.SeuObj, reduction = "umap", group.by ="celltype" ,label = TRUE, pt.size = 0.5) + NoLegend()

  scRNA.SeuObj@meta.data[[paste0("scSorterPred","_LogFC",1,"_pV",0.05)]] <- scSorter.obj[["Pred_Type"]]
  DimPlot(scRNA.SeuObj, reduction = "umap",
          group.by = paste0("scSorterPred","_LogFC",1,"_pV",0.05),label = TRUE, pt.size = 0.5) + NoLegend()
  DimPlot(scRNA.SeuObj, reduction = "umap", group.by ="celltype" ,label = TRUE, pt.size = 0.5) + NoLegend()
