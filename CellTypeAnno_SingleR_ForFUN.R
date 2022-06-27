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
  Package.set <- c("SingleR","scRNAseq","celldex","scran","scater","scuttle")
  for (i in 1:length(Package.set)) {
    if (!requireNamespace(Package.set[i], quietly = TRUE)){
      BiocManager::install(Package.set[i])
    }
  }
  ## Load Packages
  lapply(Package.set, library, character.only = TRUE)
  rm(Package.set,i)

##### Function setting #####
  ## Call function
  source("FUN_Anno_SingleR.R")

##### Current path and new folder setting* #####
  ProjectName = paste0("CTAnno_singleR_PRJCA001063S")
  Sampletype = "PDAC"
  #ProjSamp.Path = paste0(Sampletype,"_",ProjectName)

  Version = paste0(Sys.Date(),"_",ProjectName,"_",Sampletype)
  Save.Path = paste0(getwd(),"/",Version)
  ## Create new folder
  if (!dir.exists(Save.Path)){
    dir.create(Save.Path)
  }

#### Load data #####
  load("SeuratObject_CDS_PRJCA001063.RData")

  ## SeuObj_Ref
  scRNA.SeuObj_Ref <- scRNA.SeuObj
  ## For small test
  # CTFeatures.SeuObj <- scRNA.SeuObj_Ref[,scRNA.SeuObj_Ref$CELL %in% sample(scRNA.SeuObj_Ref$CELL,1000)] ## For small test
  CTFeatures.SeuObj <- scRNA.SeuObj_Ref[,scRNA.SeuObj_Ref@meta.data[[1]] %in% sample(scRNA.SeuObj_Ref@meta.data[[1]],1000)] ## For small test
  # ## For full data
  # CTFeatures.SeuObj <- scRNA.SeuObj_Ref

  ## SeuObj_Tar
  ## For small test
  # scRNA.SeuObj <- scRNA.SeuObj[,scRNA.SeuObj$CELL %in% sample(scRNA.SeuObj$CELL,1000)] ## For small test
  scRNA.SeuObj <- scRNA.SeuObj[,scRNA.SeuObj@meta.data[[1]] %in% sample(scRNA.SeuObj@meta.data[[1]],1000)] ## For small test


##### Parameter setting* #####
  Remark1 <- "PredbyscRNA" # c("PredbyCTDB","PredbyscRNA")
  RefType <- "BuiltIn_scRNA" # c("BuiltIn_celldex","BuiltIn_scRNA")
  celldexDatabase <- "HumanPrimaryCellAtlasData"
  # c("BlueprintEncodeData","DatabaseImmuneCellExpressionData","HumanPrimaryCellAtlasData","ImmGenData",
  #   "MonacoImmuneData","MouseRNAseqData","NovershternHematopoieticData")
  SingleR_DE_method <- "classic"

  ## Parameter of classifySingleR
  quantile = 0.8
  tune.thresh = 0.05
  sd.thresh = 1

  Remark <- paste0(Remark1,"_",SingleR_DE_method,"_",
                   "qua",quantile,"_tun",tune.thresh,"_sd",sd.thresh)

##### Run singleR #####
  ## Presetting
  SingleRResult.lt <- Anno_SingleR(scRNA.SeuObj, RefType = RefType, celldexDatabase = celldexDatabase,
                                   CTFeatures.SeuObj = CTFeatures.SeuObj,
                                   quantile = quantile, tune.thresh = tune.thresh, sd.thresh = sd.thresh,
                                   SingleR_DE_method = SingleR_DE_method,
                                   Remark = Remark, Save.Path = paste0(Save.Path,"/",Remark), ProjectName = "CT")

  ## Try Parameter
  SingleRResult.lt <- Anno_SingleR(scRNA.SeuObj, RefType = "BuiltIn_celldex", celldexDatabase = "HumanPrimaryCellAtlasData",
                                   quantile = quantile, tune.thresh = tune.thresh, sd.thresh = sd.thresh,
                                   Remark = "PredbyCTDB",Save.Path = Save.Path, ProjectName = ProjectName)

  scRNA.SeuObj <- SingleRResult.lt[["scRNA.SeuObj"]]
  SingleRResult2.lt <- Anno_SingleR(scRNA.SeuObj, RefType = "BuiltIn_scRNA", celldexDatabase = "HumanPrimaryCellAtlasData",
                                   quantile = quantile, tune.thresh = tune.thresh, sd.thresh = sd.thresh,
                                   Remark = "PredbyscRNA",CTFeatures.SeuObj = CTFeatures.SeuObj, SingleR_DE_method = "classic",
                                   Save.Path = Save.Path, ProjectName = ProjectName)
  scRNA.SeuObj <- SingleRResult2.lt[["scRNA.SeuObj"]]
  SingleRResult2.lt <- Anno_SingleR(scRNA.SeuObj, RefType = "BuiltIn_scRNA", celldexDatabase = "HumanPrimaryCellAtlasData",
                                    quantile = quantile, tune.thresh = tune.thresh, sd.thresh = sd.thresh,
                                    Remark = "PredbyscRNABug",CTFeatures.SeuObj = CTFeatures.SeuObj, SingleR_DE_method = "classic",
                                    Save.Path = Save.Path, ProjectName = ProjectName)


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
  CC.df <- scRNA.SeuObj@meta.data[,c("Cell_type","singleR_classic_PredbyscRNA", "singleR_classic_PredbyCTDB")]

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




