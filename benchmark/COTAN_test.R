library(bench)


# COTAN

COTAN_res <- mark({
  
  
  library(COTAN)
  library(zeallot)
  library(rlang)
  library(data.table)
  library(Rtsne)
  library(CSSG.toolkit)
  library(Seurat)
  
  
  options(parallelly.fork.enable = TRUE)
  
  cond <- "benchmark test"
  
  sc_project <- create_project(sparse_matrix_path = '../benchmark/data/', sparse_name = 'matrix', rows_name = 'genes', cols_name = 'barcodes', type = 'norm')
  
  colnames(sc_project@matrices$norm) <- make.unique(colnames(sc_project@matrices$norm))
  
  obj <- COTAN(raw = as.data.frame(sc_project@matrices$norm))
  obj <- initializeMetaDataset(obj,
                               GEO = NULL,
                               sequencingMethod = "Drop_seq",
                               sampleCondition = cond)
  
  
  UMI <- CreateSeuratObject(counts = sc_project@matrices$norm)
  
  UMI@assays$RNA$data <- UMI@assays$RNA$counts
  
  thresholds <- outlires(UMI@meta.data$nFeature_RNA)
  
  rm(UMI)
  
  genesSizeHighThr <- thresholds$thresholds[length(thresholds$thresholds)]
  obj <- addElementToMetaDataset(obj, "Num genes high threshold", genesSizeHighThr)
  
  cells_to_rem <- getCells(obj)[getNumExpressedGenes(obj) > genesSizeHighThr]
  obj <- dropGenesCells(obj, cells = cells_to_rem)
  
  
  
  genesSizeLowThr <- thresholds$thresholds[1]
  obj <- addElementToMetaDataset(obj, "Num genes low threshold", genesSizeLowThr)
  
  cells_to_rem <- getCells(obj)[getNumExpressedGenes(obj) < genesSizeLowThr]
  obj <- dropGenesCells(obj, cells = cells_to_rem)
  
  
  c(mitPlot, mitSizes) %<-% mitochondrialPercentagePlot(obj, genePrefix = "^MT")
  
  
  mitPercThr <- 20
  obj <- addElementToMetaDataset(obj, "Mitoc. perc. threshold", mitPercThr)
  
  cells_to_rem <- rownames(mitSizes)[mitSizes[["mit.percentage"]] > mitPercThr]
  obj <- dropGenesCells(obj, cells = cells_to_rem)
  
  c(mitPlot, mitSizes) %<-% mitochondrialPercentagePlot(obj, genePrefix = "^MT")
  
  genes_to_rem <- getGenes(obj)[grep("^MT", getGenes(obj))]
  cells_to_rem <- getCells(obj)[which(getCellsSize(obj) == 0L)]
  
  obj <- dropGenesCells(obj, genes_to_rem, cells_to_rem)
  
  obj <- clean(obj)
  
  c(pcaCellsPlot, pcaCellsData, genesPlot,
    UDEPlot, nuPlot, zoomedNuPlot) %<-% cleanPlots(obj)
  
  
  obj <- proceedToCoex(obj, calcCoex = TRUE,
                       optimizeForSpeed = TRUE, cores = 6L, deviceStr = "cpu",
                       saveObj = FALSE, outDir = outDir)
  
  
  c(splitClusters, splitCoexDF) %<-%
    cellsUniformClustering(obj, initialResolution = 0.6, checker = NULL,
                           optimizeForSpeed = FALSE, deviceStr = "cpu",
                           cores = 6L, genesSel = "HVG_Seurat",
                           saveObj = FALSE)
  
  obj <- addClusterization(obj, clName = "split",
                           clusters = splitClusters, coexDF = splitCoexDF)
  
  
  
  c(umapPlot, cellsPCA) %<-% cellsUMAPPlot(obj,
                                           dataMethod = "LogLikelihood",
                                           genesSel = "HVG_Seurat",
                                           colors = NULL, numNeighbors = 30L,
                                           minPointsDist = 0.2)
  
  
  
  
  
}, iterations = 1)





df_to_save = data.frame(list('mem' = COTAN_res$mem_alloc, 'time' = COTAN_res$total_time, 'test' = 'COTAN'))

write.csv(df_to_save, 'test_stats/cotan.csv', quote = F, row.names = F)
