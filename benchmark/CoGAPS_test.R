library(bench)



# CoGAPS

CoGAPS_res <- mark({
  
  
  library(CSSG.toolkit)
  library(projectR)
  library(CoGAPS)
  library(Seurat)
  
  
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
  
  
  input <- as.data.frame(as.matrix(sc_project@matrices$norm))[VariableFeatures(object = UMI),]
  
  params <- new("CogapsParams", nPatterns = 4, nIterations = 1000)
  
  scCoGAPS_result <- CoGAPS(input, params)
  
  P_matrix <- scCoGAPS_result@featureLoadings
  
  mtx <- as.matrix(input)
  
  projected_patterns <- projectR(mtx ,loadings=P_matrix, full=TRUE,  dataNames= rownames(P_matrix))
  
  
  pval = projected_patterns$pval
  
  p = projected_patterns$projection
  
  
  p[pval > 0.05] <- NA
  
  best_patterns_per_column <- as.vector(unlist(apply(p, 2, function(x) which.max(x))))
  
  best_patterns_per_column[is.na(best_patterns_per_column)] = 'Undefined'
  
  
  Idents(UMI) <- as.vector(best_patterns_per_column)
  
  
  
}, iterations = 1)




df_to_save = data.frame(list('mem' = CoGAPS_res$mem_alloc, 'time' = CoGAPS_res$total_time, 'test' = 'CoGAPS'))

write.csv(df_to_save, 'test_stats/cogaps.csv', quote = F, row.names = F)

