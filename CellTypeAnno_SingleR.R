## SingleRBook Ref: http://bioconductor.org/books/release/SingleRBook/
## Example Ref: https://bioconductor.org/packages/devel/bioc/vignettes/SingleR/inst/doc/SingleR.html

##### Presetting ######
  rm(list = ls()) # Clean variable
  memory.limit(150000)


##### Load Packages #####
  if(!require("tidyverse")) install.packages("tidyverse")
  library(tidyverse)

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


  #### Database setting for single cell gene expression ####
  library(scRNAseq)
  hESCs <- LaMannoBrainData('human-es')
  hESCs <- hESCs[,colSums(counts(hESCs)) > 0] # Remove libraries with no counts.
  hESCs <- logNormCounts(hESCs)
  hESCs <- hESCs[,1:100]

  #### Run SingleR ####
  library(SingleR)
  pred.hesc <- SingleR(test = hESCs, ref = CTFeatures, assay.type.test=1,
                       labels = CTFeatures$label.main)

  pred.hesc

  # Summarizing the distribution:
  table(pred.hesc$labels)


  ##### Annotation diagnostics #####
  plotScoreHeatmap(pred.hesc)
  plotDeltaDistribution(pred.hesc, ncol = 3)
  summary(is.na(pred.hesc$pruned.labels))

  all.markers <- metadata(pred.hesc)$de.genes
  hESCs$labels <- pred.hesc$labels

  # Beta cell-related markers
  library(scater)
  plotHeatmap(hESCs, order_columns_by="labels",
              features = unique(unlist(all.markers[["B cells"]])))


##### Using single-cell references   #####
  library(scRNAseq)
  sceM <- MuraroPancreasData()

  # One should normally do cell-based quality control at this point, but for
  # brevity's sake, we will just remove the unlabelled libraries here.
  sceM <- sceM[,!is.na(sceM$label)]

  # SingleR() expects reference datasets to be normalized and log-transformed.
  library(scuttle)
  sceM <- logNormCounts(sceM)

  sceG <- GrunPancreasData()
  sceG <- sceG[,colSums(counts(sceG)) > 0] # Remove libraries with no counts.
  sceG <- logNormCounts(sceG)
  sceG <- sceG[,1:100]

  library(SingleR)
  pred.grun <- SingleR(test=sceG, ref=sceM, labels=sceM$label, de.method="wilcox")
  table(pred.grun$labels)

##### Annotation diagnostics #####
  plotScoreHeatmap(pred.grun)
  plotDeltaDistribution(pred.grun, ncol = 3)
  summary(is.na(pred.grun$pruned.labels))

  all.markers <- metadata(pred.grun)$de.genes
  sceG$labels <- pred.grun$labels

  # Beta cell-related markers
  library(scater)
  plotHeatmap(sceG, order_columns_by="labels",
              features=unique(unlist(all.markers$beta)))


##### Session information #####
  sessionInfo()


