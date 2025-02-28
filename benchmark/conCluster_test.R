library(bench)


# conCluster

conCluster_res <- mark({
  
  
  library(CSSG.toolkit)
  library(Seurat)
  library(openxlsx)
  library(dplyr)
  library(ConsensusClusterPlus)
  library(clusterCrit) 
  
  
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
  
  Elbow <- Elbow + geom_vline(xintercept = dim, color = 'red') +   geom_text(aes(x = dim + 3, y = round(max(Elbow$data$stdev)/2,0), label = paste("Dim =", dim)), color = 'red', vjust = -1)
  
  UMI <- JackStraw(UMI, num.replicate = 10, dims = dim)
  UMI <- ScoreJackStraw(UMI, dims = 1:dim)
  
  jc <- as.data.frame(UMI@reductions$pca@jackstraw@overall.p.values)
  jc <- jc[jc$Score < 0.05,]
  dim <- as.vector(jc$PC)
  
  
  UMI <- RunUMAP(UMI, dims = dim)
  
  
  UMI <- RunTSNE(UMI, dims = dim)
  
  set.seed(123)
  
  tsne_coords <- Embeddings(UMI, "tsne")
  num_clusters <- 10
  clustering_results <- replicate(10, kmeans(tsne_coords, centers = num_clusters)$cluster)
  
  binary_matrix <- do.call(cbind, lapply(1:ncol(clustering_results), function(i) {
    table(rownames(tsne_coords), clustering_results[, i])
  }))
  
  consensus_result <- ConsensusClusterPlus(t(binary_matrix), maxK = 10, reps = 50, pItem = 1, 
                                           pFeature = 1, clusterAlg = "hc", distance = "euclidean")
  
  
  k <- c()
  inx <- c()
  for (i in 4:num_clusters) {
    
    print(i)
    cluster_assignments <- consensus_result[[i]]$consensusClass
    ch_index <- intCriteria(tsne_coords, as.integer(cluster_assignments), "Calinski_Harabasz")$calinski_harabasz
    
    k <- c(k,i)
    inx <- c(inx, ch_index)
    
  }
  
  
  inx_df <- data.frame(list('k' = k, 'inx' = inx))
  inx_df <- inx_df[order(-inx_df$inx), ]
  
  
  
  clusters <- consensus_result[[inx_df$k[1]]]$consensusClass
  
  
  Idents(UMI) <- unname(clusters)[match(names(Idents(UMI)), names(clusters))]
  
  
}, iterations = 1)




df_to_save = data.frame(list('mem' = conCluster_res$mem_alloc, 'time' = conCluster_res$total_time, 'test' = 'conCluster'))

write.csv(df_to_save, 'test_stats/concluster.csv', quote = F, row.names = F)

