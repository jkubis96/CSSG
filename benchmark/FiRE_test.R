library(bench)



# FiRE

FiRE_res <- mark({
  
  
  library(CSSG.toolkit)
  library(Seurat)
  library(openxlsx)
  library(dplyr)
  
  
  sc_project <- create_project(sparse_matrix_path = '../benchmark/data/', sparse_name = 'matrix', rows_name = 'genes', cols_name = 'barcodes', type = 'norm')
  
  colnames(sc_project@matrices$norm) <- make.unique(colnames(sc_project@matrices$norm))
  
  UMI <- CreateSeuratObject(counts = sc_project@matrices$norm)
  
  UMI@assays$RNA$data <- UMI@assays$RNA$counts
  
  UMI[['MitoPercent']] <- PercentageFeatureSet(UMI, pattern = "^MT-") + PercentageFeatureSet(UMI, pattern = "^Mt-") + PercentageFeatureSet(UMI, pattern = "^mt-")
  UMI[['RiboPercent']] <- PercentageFeatureSet(UMI, pattern = "^Rps-") + PercentageFeatureSet(UMI, pattern = "^Rpl") + PercentageFeatureSet(UMI, pattern = "^RPS-") + PercentageFeatureSet(UMI, pattern = "^RPL")
  
  thresholds <- outlires(UMI@meta.data$nFeature_RNA)
  
  
  UMI <- subset(UMI, subset = nFeature_RNA > thresholds$thresholds[1] & nFeature_RNA <= thresholds$thresholds[length(thresholds$thresholds)]
                & MitoPercent < 20)
  
  
  UMI <- FindVariableFeatures(UMI, selection.method = "vst", nfeatures = 5000, binning.method = 'equal_frequency')
  
  
  all.genes <- rownames(UMI)
  UMI <- ScaleData(UMI, features = all.genes)
  
  UMI <- RunPCA(UMI, features = VariableFeatures(object = UMI))
  
  Elbow <- ElbowPlot(UMI, ndims = 50)
  
  dims <- as.data.frame(Elbow$data$stdev)
  
  dim <- dim_reuction_pcs(dims)
  
  UMI <- JackStraw(UMI, num.replicate = 10, dims = dim)
  UMI <- ScoreJackStraw(UMI, dims = 1:dim)
  
  jc <- as.data.frame(UMI@reductions$pca@jackstraw@overall.p.values)
  jc <- jc[jc$Score < 0.05,]
  dim <- as.vector(jc$PC)
  
  
  UMI <- RunUMAP(UMI, dims = dim)
  
  sc_project <- create_project_from_seurat(UMI)
  
  
  library('FiRE')
  
  data <- as.matrix(UMI@assays$RNA$data[VariableFeatures(object = UMI),])
  colnames(data) <- as.vector(unname(Idents(UMI)))
  data <- t(data) 
  
  genes <- colnames(data)
  
  model <- new(FiRE::FiRE, 1000, 100, 1017881, 5489, 0)
  
  model$fit(as.matrix(data))
  
  score <- model$score(as.matrix(data))
  
  
  q3 <- quantile(score, 0.75)
  iqr <- IQR(score)
  th <- q3 + (1.5 * iqr)
  
  indIqr <- which(score >= th)
  
  predictions <- integer(nrow(data))
  
  
  rownames(data)[indIqr]  <- 'Rare'
  
  
  Idents(UMI) <- rownames(data)
  
  
}, iterations = 1)



df_to_save = data.frame(list('mem' = FiRE_res$mem_alloc, 'time' = FiRE_res$total_time, 'test' = 'FiRE'))

write.csv(df_to_save, 'test_stats/fire.csv', quote = F, row.names = F)

