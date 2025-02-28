library(bench)

#CSSG test

cssg_res <- mark({
  
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
  
  sc_project <- create_project_from_seurat(UMI)
  
  
  markers_class = read.xlsx('../markers/non_canonical.xlsx', sheet = 1)
  markers_subclass = read.xlsx('../markers/non_canonical.xlsx', sheet = 2)
  
  sc_project <- heterogeneity_select_variance(sc_project = sc_project, heterogeneity_factor = 80, max_genes = 50, min_occ = 5, min_exp = 0.3, rep_factor = 0.5, mito_content = FALSE)
  
  sc_project <- subclass_naming(sc_project = sc_project, class_markers = markers_class, subclass_markers = markers_subclass, species = 'Homo sapiens', chunk_size = 5000)
  
  sc_project <- CSSG_markers(sc_project = sc_project, max_combine = 1000, loss_val = 0.05)
  
  sc_project <- subtypes_naming(sc_project = sc_project, markers_class = markers_class, markers_subclass = markers_subclass, species = 'Homo sapiens') 
  
  thr_data <- bin_cell_test(p_val = 0.05, names = sc_project@names$repaired, min_cells = 10)
  
  Idents(UMI) <- sc_project@names$repaired$renamed_idents
  
  select_list <- thr_data$data$names[thr_data$data$test %in% c("Good marked types", "Renamed")]
  
  UMI <- subset(UMI, idents = select_list)
  
}, iterations = 1)




df_to_save = data.frame(list('mem' = cssg_res$mem_alloc, 'time' = cssg_res$total_time, 'test' = 'CSSG'))

write.csv(df_to_save, 'test_stats/cssg.csv', quote = F, row.names = F)
