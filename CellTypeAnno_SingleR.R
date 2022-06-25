## SingleRBook Ref: http://bioconductor.org/books/release/SingleRBook/
## Example Ref: https://bioconductor.org/packages/devel/bioc/vignettes/SingleR/inst/doc/SingleR.html

##### Presetting ######
  rm(list = ls()) # Clean variable
  memory.limit(150000)


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
  Package.set <- c("SingleR","scRNAseq","celldex","scran")
  for (i in 1:length(Package.set)) {
    if (!requireNamespace(Package.set[i], quietly = TRUE)){
      BiocManager::install(Package.set[i])
    }
  }
  ## Load Packages
  lapply(Package.set, library, character.only = TRUE)
  rm(Package.set,i)

##### Current path and new folder setting* #####
  ProjectName = "MaliAnno_singleR_PRJCA001063"
  Sampletype = "PDAC"
  #ProjSamp.Path = paste0(Sampletype,"_",ProjectName)

  Version = paste0(Sys.Date(),"_",ProjectName,"_",Sampletype)
  Save.Path = paste0(getwd(),"/",Version)
  ## Create new folder
  if (!dir.exists(Save.Path)){
    dir.create(Save.Path)
  }

##### Parameter setting* #####
  singleRDatabase <- "HumanPrimaryCellAtlasData"
  # c("BlueprintEncodeData","DatabaseImmuneCellExpressionData","HumanPrimaryCellAtlasData","ImmGenData",
  #   "MonacoImmuneData","MouseRNAseqData","NovershternHematopoieticData")


##### Using built-in references #####

  #### Database setting for Cell type features ####
  library(celldex)

  if(singleRDatabase == "BlueprintEncodeData"){
    CTFeatures <- BlueprintEncodeData()
  }else if(singleRDatabase == "DatabaseImmuneCellExpressionData"){
    CTFeatures <- DatabaseImmuneCellExpressionData()
  }else if(singleRDatabase == "HumanPrimaryCellAtlasData"){
    CTFeatures <- HumanPrimaryCellAtlasData()
  }else if(singleRDatabase == "ImmGenData"){
    CTFeatures <- ImmGenData()
  }else if(singleRDatabase == "MonacoImmuneData"){
    CTFeatures <- MonacoImmuneData()
  }else if(singleRDatabase == "MouseRNAseqData"){
    CTFeatures <- MouseRNAseqData()
  }else if(singleRDatabase == "NovershternHematopoieticData"){
    CTFeatures <- NovershternHematopoieticData()
  }else{
    print("Error in database setting!")
  }

  # library(celldex)
  # hpca.se <- HumanPrimaryCellAtlasData()
  # hpca.se

  #### scRNA-seq object setting for gene expression matrix ####
  ## A numeric matrix of single-cell expression values where rows are genes and columns are cells.
  load("D:/Dropbox/##_GitHub/##_Charlene/TrajectoryAnalysis/SeuratObject_CDS_PRJCA001063_V2.RData")
  #scRNA.SeuObj <- scRNA.SeuObj[,scRNA.SeuObj$CELL %in% sample(scRNA.SeuObj$CELL,1000)] ## For small test
  scRNA <- as.SingleCellExperiment(scRNA.SeuObj)

  # library(scRNAseq)
  # hESCs <- LaMannoBrainData('human-es')
  # hESCs <- hESCs[,colSums(counts(hESCs)) > 0] # Remove libraries with no counts.
  # hESCs <- logNormCounts(hESCs)
  # hESCs <- hESCs[,1:100]

  #### Run SingleR ####
  library(SingleR)


  Pred_byCTDB <- SingleR(test = scRNA, ref = CTFeatures, assay.type.test=1,
                         labels = CTFeatures$label.main)#, de.method="wilcox") #  de.method = c("classic", "wilcox", "t")

  Pred_byCTDB

  # Summarizing the distribution:
  table(Pred_byCTDB$labels)


  ##### Annotation diagnostics #####
  p.ScoreHeatmap1 <- plotScoreHeatmap(Pred_byCTDB)
  p.ScoreHeatmap1
  p.DeltaDist1 <- plotDeltaDistribution(Pred_byCTDB, ncol = 3)
  p.DeltaDist1
  summary(is.na(Pred_byCTDB$pruned.labels))

  pdf(file = paste0(Save.Path,"/",ProjectName,"_PredbyCTDB_AnnoDiag.pdf"),
      width = 7,  height = 7
  )
    p.ScoreHeatmap1
    p.DeltaDist1
  dev.off()



  all.markers <- metadata(Pred_byCTDB)$de.genes
  scRNA$labels <- Pred_byCTDB$labels

  # Endothelial cell-related markers
  library(scater)
  plotHeatmap(scRNA, order_columns_by="labels",
              features = unique(unlist(all.markers[["Endothelial_cells"]])))



  pdf(file = paste0(Save.Path,"/",ProjectName,"_PredbyCTDB_HeatmapCTmarkers.pdf"),
      width = 12,  height = 7
  )
  for (i in 1:length(all.markers)) {
    plotHeatmap(scRNA, order_columns_by="labels",
                features=unique(unlist(all.markers[[i]]))) %>% print()
  }
  dev.off()




  ## Plot UMAP
  scRNA.SeuObj$singleRPredbyCTDB <- Pred_byCTDB$labels
  p.CTPred1 <- DimPlot(scRNA.SeuObj, reduction = "umap", group.by ="singleRPredbyCTDB" ,label = TRUE, pt.size = 0.5) + NoLegend()
  p.CTPred1
  p.CT1 <- DimPlot(scRNA.SeuObj, reduction = "umap", group.by ="Cell_type" ,label = TRUE, pt.size = 0.5) + NoLegend()
  p.CT1

  library(ggpubr)
  p.CTComp1 <- ggarrange(p.CT1, p.CTPred1, common.legend = TRUE, legend = "top")
  p.CTComp1

  pdf(file = paste0(Save.Path,"/",ProjectName,"_PredbyCTDB_CompareCTUMAP.pdf"),
      width = 12,  height = 7
  )
    p.CTComp1
  dev.off()
##### Using single-cell references   #####

  #### single-cell reference setting for Cell type features ####
  load("D:/Dropbox/##_GitHub/##_Charlene/TrajectoryAnalysis/SeuratObject_CDS_PRJCA001063_V2.RData")
  #scRNA_Ref.SeuObj <- scRNA.SeuObj[,scRNA.SeuObj$CELL %in% sample(scRNA.SeuObj$CELL,1000)] ## For small test
  scRNA_Ref.SeuObj <- scRNA.SeuObj
  scRNA_Ref <- as.SingleCellExperiment(scRNA_Ref.SeuObj)
  scRNA_Ref$label <- scRNA_Ref@colData@listData[["Cell_type"]]
  scRNA_Ref <- scRNA_Ref[,!is.na(scRNA_Ref$label)]

  # library(scRNAseq)
  # sceM <- MuraroPancreasData()
  #
  # # One should normally do cell-based quality control at this point, but for
  # # brevity's sake, we will just remove the unlabelled libraries here.
  # sceM <- sceM[,!is.na(sceM$label)]
  #
  # # SingleR() expects reference datasets to be normalized and log-transformed.
  # library(scuttle)
  # sceM <- logNormCounts(sceM)


  #### scRNA-seq object setting for gene expression matrix ####
  ## A numeric matrix of single-cell expression values where rows are genes and columns are cells.
  load("D:/Dropbox/##_GitHub/##_Charlene/TrajectoryAnalysis/SeuratObject_CDS_PRJCA001063_V2.RData")
  #scRNA.SeuObj<- scRNA.SeuObj[,scRNA.SeuObj$CELL %in% sample(scRNA.SeuObj$CELL,1000)] ## For small test
  scRNA <- as.SingleCellExperiment(scRNA.SeuObj)

  # sceG <- GrunPancreasData()
  # sceG <- sceG[,colSums(counts(sceG)) > 0] # Remove libraries with no counts.
  # sceG <- logNormCounts(sceG)
  # sceG <- sceG[,1:100]
  #
  #### sceG <- as.SingleCellExperiment(scRNA.SeuObj)

  library(SingleR)
  Pred_byscRNA <- SingleR(test = scRNA, ref=scRNA_Ref,
                          labels=scRNA_Ref$label, de.method="wilcox") #  de.method = c("classic", "wilcox", "t")
  table(Pred_byscRNA$labels)

##### Annotation diagnostics #####
  p.ScoreHeatmap2 <- plotScoreHeatmap(Pred_byscRNA)
  p.DeltaDist2 <- plotDeltaDistribution(Pred_byscRNA, ncol = 3)
  summary(is.na(Pred_byscRNA$pruned.labels))

  pdf(file = paste0(Save.Path,"/",ProjectName,"_PredbyscRNA_AnnoDiag.pdf"),
      width = 7,  height = 7
  )
    p.ScoreHeatmap2
    p.DeltaDist2
  dev.off()


  all.markers <- metadata(Pred_byscRNA)$de.genes
  scRNA$labels <- Pred_byscRNA$labels

  # B cell-related markers
  library(scater)
  plotHeatmap(scRNA, order_columns_by="labels",
              features=unique(unlist(all.markers[["B cell"]])))

  pdf(file = paste0(Save.Path,"/",ProjectName,"_PredbyscRNA_HeatmapCTmarkers.pdf"),
      width = 12,  height = 7
  )
  for (i in 1:length(all.markers)) {
    plotHeatmap(scRNA, order_columns_by="labels",
                features=unique(unlist(all.markers[[i]]))) %>% print()
  }
  dev.off()


  ## Plot UMAP
  scRNA.SeuObj$singleRPredbyscRNA <- Pred_byscRNA$labels
  p.CTPred2 <- DimPlot(scRNA.SeuObj, reduction = "umap", group.by ="singleRPredbyscRNA" ,label = TRUE, pt.size = 0.5) + NoLegend()
  p.CTPred2
  p.CT2 <- DimPlot(scRNA.SeuObj, reduction = "umap", group.by ="Cell_type" ,label = TRUE, pt.size = 0.5) + NoLegend()
  p.CT2

  library(ggpubr)
  p.CTComp2 <- ggarrange(p.CT2, p.CTPred2, common.legend = TRUE, legend = "top")
  p.CTComp2

  pdf(file = paste0(Save.Path,"/",ProjectName,"_PredbyscRNA_CompareCTUMAP.pdf"),
      width = 12,  height = 7
  )
    p.CTComp2
  dev.off()

##### Session information #####
  sessionInfo()

##### Save RData #####
  save.image(paste0(Save.Path,"/SeuratObject_",ProjectName,".RData"))




