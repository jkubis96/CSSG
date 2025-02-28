library(bench)

# Seurat

seurat_res <- mark({
  
  library(CSSG.toolkit)
  library(Seurat)
  library(openxlsx)
  library(dplyr)
  
  
  sc_project <- create_project(sparse_matrix_path = '../benchmark/data/', sparse_name = 'matrix', rows_name = 'genes', cols_name = 'barcodes', type = 'norm')
  
  colnames(sc_project@matrices$norm) <- make.unique(colnames(sc_project@matrices$norm))
  
  UMI <- CreateSeuratObject(counts = sc_project@matrices$norm)
  
  UMI@assays$RNA$data <- UMI@assays$RNA$counts
  
  UMI[['MitoPercent']] <- PercentageFeatureSet(UMI, pattern = "^MT-") + PercentageFeatureSet(UMI, pattern = "^Mt-") + PercentageFeatureSet(UMI, pattern = "^mt-")
  
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
  
  UMI <- FindNeighbors(UMI, dims = dim, reduction = 'pca', nn.method="rann")
  UMI <- FindClusters(UMI, resolution = 0.6, n.start = 10, n.iter = 1000)
  
  UMI <- RunUMAP(UMI, dims = dim)
  
  
}, iterations = 1)




df_to_save = data.frame(list('mem' = seurat_res$mem_alloc, 'time' = seurat_res$total_time, 'test' = 'Basic_Seurat'))

write.csv(df_to_save, 'test_stats/seurat.csv', quote = F, row.names = F)
