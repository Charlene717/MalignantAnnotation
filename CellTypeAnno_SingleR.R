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
  ProjectName = "MaliAnno_singleR_PRJCA001063_Small"
  Sampletype = "PDAC"
  #ProjSamp.Path = paste0(Sampletype,"_",ProjectName)

  Version = paste0(Sys.Date(),"_",ProjectName,"_",Sampletype)
  Save.Path = paste0(getwd(),"/",Version)
  ## Create new folder
  if (!dir.exists(Save.Path)){
    dir.create(Save.Path)
  }

##### Parameter setting* #####
  celldexDatabase <- "HumanPrimaryCellAtlasData"
  # c("BlueprintEncodeData","DatabaseImmuneCellExpressionData","HumanPrimaryCellAtlasData","ImmGenData",
  #   "MonacoImmuneData","MouseRNAseqData","NovershternHematopoieticData")


  #LabelName =

##### Using built-in references #####

  #### Database: Bulk reference setting for Cell type features ####
  library(celldex)

  if(celldexDatabase == "BlueprintEncodeData"){
    CTFeatures <- BlueprintEncodeData()
  }else if(celldexDatabase == "DatabaseImmuneCellExpressionData"){
    CTFeatures <- DatabaseImmuneCellExpressionData()
  }else if(celldexDatabase == "HumanPrimaryCellAtlasData"){
    CTFeatures <- HumanPrimaryCellAtlasData()
  }else if(celldexDatabase == "ImmGenData"){
    CTFeatures <- ImmGenData()
  }else if(celldexDatabase == "MonacoImmuneData"){
    CTFeatures <- MonacoImmuneData()
  }else if(celldexDatabase == "MouseRNAseqData"){
    CTFeatures <- MouseRNAseqData()
  }else if(celldexDatabase == "NovershternHematopoieticData"){
    CTFeatures <- NovershternHematopoieticData()
  }else{
    print("Error in database setting!")
  }

    #### Demo dataset ####
    # library(celldex)
    # hpca.se <- HumanPrimaryCellAtlasData()
    # hpca.se

  #### scRNA-seq object setting for gene expression matrix ####
  ## A numeric matrix of single-cell expression values where rows are genes and columns are cells.
  load("D:/Dropbox/##_GitHub/##_Charlene/TrajectoryAnalysis/SeuratObject_CDS_PRJCA001063_V2.RData")
  rm(scRNA.Mono.SeuObj)

  ## For small test
  # scRNA.SeuObj <- scRNA.SeuObj[,scRNA.SeuObj$CELL %in% sample(scRNA.SeuObj$CELL,1000)] ## For small test
  scRNA.SeuObj <- scRNA.SeuObj[,scRNA.SeuObj@meta.data[[1]] %in% sample(scRNA.SeuObj@meta.data[[1]],1000)] ## For small test

  scRNA <- as.SingleCellExperiment(scRNA.SeuObj)

    #### Demo dataset ####
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
  CTCount_byCTDB.df <- table(Pred_byCTDB$labels) %>%
                       as.data.frame() %>%
                       dplyr::rename(Cell_Type = Var1, Count = Freq)



  ##### Annotation diagnostics #####
  p.ScoreHeatmap1 <- plotScoreHeatmap(Pred_byCTDB)
  p.ScoreHeatmap1
  p.DeltaDist1 <- plotDeltaDistribution(Pred_byCTDB, ncol = 3)
  p.DeltaDist1
  summary(is.na(Pred_byCTDB$pruned.labels))

  pdf(file = paste0(Save.Path,"/",ProjectName,"_PredbyCTDB_AnnoDiag.pdf"),
      width = 10,  height = 7
  )
    p.ScoreHeatmap1 %>% print()
    p.DeltaDist1 %>% print()
  dev.off()



  all.markers <- metadata(Pred_byCTDB)$de.genes
  scRNA$labels <- Pred_byCTDB$labels  ## scRNA@colData@listData[["labels"]] <- Pred_byCTDB$labels


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
    p.CTComp1 %>% print()
  dev.off()


##### Using single-cell references   #####
  #### Customization: single-cell reference setting for Cell type features ####

  # load("D:/Dropbox/##_GitHub/##_Charlene/TrajectoryAnalysis/SeuratObject_CDS_PRJCA001063_V2.RData")

  ## For small test
  # scRNA_Ref.SeuObj <- scRNA.SeuObj[,scRNA.SeuObj$CELL %in% sample(scRNA.SeuObj$CELL,1000)] ## For small test
  scRNA_Ref.SeuObj <- scRNA.SeuObj[,scRNA.SeuObj@meta.data[[1]] %in% sample(scRNA.SeuObj@meta.data[[1]],1000)] ## For small test

  # ## For full data
  # scRNA_Ref.SeuObj <- scRNA.SeuObj
  scRNA_Ref <- as.SingleCellExperiment(scRNA_Ref.SeuObj)
  scRNA_Ref$label <- scRNA_Ref@colData@listData[["Cell_type"]]
  scRNA_Ref <- scRNA_Ref[,!is.na(scRNA_Ref$label)]
  rm(scRNA_Ref.SeuObj)

    #### Demo dataset ####
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
  # ## A numeric matrix of single-cell expression values where rows are genes and columns are cells.
  # load("D:/Dropbox/##_GitHub/##_Charlene/TrajectoryAnalysis/SeuratObject_CDS_PRJCA001063_V2.RData")
  #
  # ## For small test
  # # scRNA.SeuObj <- scRNA.SeuObj[,scRNA.SeuObj$CELL %in% sample(scRNA.SeuObj$CELL,1000)] ## For small test
  # scRNA.SeuObj <- scRNA.SeuObj[,scRNA.SeuObj@meta.data[[1]] %in% sample(scRNA.SeuObj@meta.data[[1]],1000)] ## For small test
  #
  # scRNA <- as.SingleCellExperiment(scRNA.SeuObj)

    #### Demo dataset ####
    # sceG <- GrunPancreasData()
    # sceG <- sceG[,colSums(counts(sceG)) > 0] # Remove libraries with no counts.
    # sceG <- logNormCounts(sceG)
    # sceG <- sceG[,1:100]
    #
    #### sceG <- as.SingleCellExperiment(scRNA.SeuObj)

  #### Run SingleR ####
  library(SingleR)
  Pred_byscRNA <- SingleR(test = scRNA, ref=scRNA_Ref,
                          labels=scRNA_Ref$label, de.method="wilcox") #  de.method = c("classic", "wilcox", "t")
  Pred_byscRNA
  # Summarizing the distribution:
  CTCount_byscRNA.df <- table(Pred_byscRNA$labels) %>%
                        as.data.frame() %>%
                        dplyr::rename(Cell_Type = Var1, Count = Freq)

  ##### Annotation diagnostics #####
    p.ScoreHeatmap2 <- plotScoreHeatmap(Pred_byscRNA)
    p.DeltaDist2 <- plotDeltaDistribution(Pred_byscRNA, ncol = 3)
    summary(is.na(Pred_byscRNA$pruned.labels))

    pdf(file = paste0(Save.Path,"/",ProjectName,"_PredbyscRNA_AnnoDiag.pdf"),
        width = 10,  height = 7
    )
      p.ScoreHeatmap2 %>% print()
      p.DeltaDist2 %>% print()
    dev.off()


    all.markers2 <- metadata(Pred_byscRNA)$de.genes
    scRNA$labels2 <- Pred_byscRNA$labels

    # B cell-related markers
    library(scater)
    plotHeatmap(scRNA, order_columns_by="labels2",
                features=unique(unlist(all.markers2[["B cell"]])))

    pdf(file = paste0(Save.Path,"/",ProjectName,"_PredbyscRNA_HeatmapCTmarkers.pdf"),
        width = 12,  height = 7
    )
      for (i in 1:length(all.markers2)) {
        plotHeatmap(scRNA, order_columns_by="labels",
                    features=unique(unlist(all.markers2[[i]]))) %>% print()
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
      p.CTComp2 %>% print()
    dev.off()

##### Session information #####
  sessionInfo()
  ## Ref: https://stackoverflow.com/questions/21967254/how-to-write-a-reader-friendly-sessioninfo-to-text-file
  writeLines(capture.output(sessionInfo()), paste0(Save.Path,"/sessionInfo.txt"))

##### Verification (CellCheck) #####
  #### Install ####
  ## Check whether the installation of those packages is required
  Package.set <- c("tidyverse","caret","cvms","DescTools","devtools","ggthemes")
  for (i in 1:length(Package.set)) {
    if (!requireNamespace(Package.set[i], quietly = TRUE)){
      install.packages(Package.set[i])
    }
  }
  ## Load Packages
  # library(Seurat)
  lapply(Package.set, library, character.only = TRUE)
  rm(Package.set,i)

  ## install CellCheck
  # Install the CellCheck package
  detach("package:CellCheck", unload = TRUE)
  devtools::install_github("Charlene717/CellCheck")
  # Load CellCheck
  library(CellCheck)

  #### Run CellCheck ####
  ## Create check dataframe
  CC.df <- scRNA.SeuObj@meta.data[,c("Cell_type","singleRPredbyscRNA", "singleRPredbyCTDB")]

  CC.df <- data.frame(lapply(CC.df, as.character), stringsAsFactors=FALSE)

  colnames(CC.df) <- c("Actual","Predict1","Predict2")
  #CC.df$Actual <- as.character(CC.df$Actual)

  CC.df$Predict2 <- gsub("_", " ", CC.df$Predict2)
  CC.df$Predict2 <- gsub("cells", "cell", CC.df$Predict2)
  CC.df$Predict2 <- gsub("Macrophage", "Macrophage cell", CC.df$Predict2)
  CC.df$Predict2 <- gsub("Fibroblasts", "Fibroblast cell", CC.df$Predict2)
  CC.df$Predict2 <- gsub("Epithelial cell", "Ductal cell type 1", CC.df$Predict2)



  CC.df[!CC.df$Predict2 %in% c(CC.df$Actual %>% unique()),]$Predict2 <- "Other"
  # CC.df <- rbind(CC.df,"NotMatch")  #CC.df[nrow(CC.df)+1,1:ncol(CC.df)] <- "Other"

  CC_Anno.df <- data.frame(TestID = c("Predict1","Predict2"),
                           Tool = "singleR",
                           Type = "PDAC",
                           Set = c("singleRPredbyscRNA", "singleRPredbyCTDB"))

  ## For one prediction
  ## For one prediction
  DisCMSet.lt = list(Mode = "One", Actual = "Actual", Predict = "Predict1" , FilterSet1 = "Tool", FilterSet2 = "singleR" , Remark = "") # Mode = c("One","Multiple")
  BarChartSet.lt <- list(Mode = "One", Metrics = "Balanced.Accuracy", XValue = "Set", Group = "Tool", Remark = "")
  LinePlotSet.lt <- list(Mode = "One", Metrics = "Balanced.Accuracy", XValue = "Set", Group = "Tool", Remark = "")
  CCR_cm_DisMult.lt <- CellCheck_DisMult(CC.df, CC_Anno.df,
                                         DisCMSet.lt = DisCMSet.lt,
                                         BarChartSet.lt = BarChartSet.lt,
                                         LinePlotSet.lt = LinePlotSet.lt,
                                         Save.Path = Save.Path, ProjectName = paste0("CellCheck_",ProjectName))



##### Save RData #####
  save.image(paste0(Save.Path,"/SeuratObject_",ProjectName,".RData"))




