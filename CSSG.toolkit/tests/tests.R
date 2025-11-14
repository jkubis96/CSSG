library(Seurat)
library(CSSG.toolkit)
library(testthat)

test_that("full CSSG pipeline runs", {

  # ====================
  # 1. LOAD PROJECT
  # ====================
  sc_project <- create_project(
    sparse_matrix_path = "../../benchmark/data/",
    sparse_name = "matrix",
    rows_name = "genes",
    cols_name = "barcodes",
    type = "norm"
  )

  expect_gt(nrow(sc_project@matrices$norm), 0)


  # ====================
  # 2. PREPARE SEURAT
  # ====================
  colnames(sc_project@matrices$norm) <- make.unique(colnames(sc_project@matrices$norm))

  UMI <- CreateSeuratObject(counts = sc_project@matrices$norm)
  UMI@assays$RNA$data <- UMI@assays$RNA$counts

  UMI[["MitoPercent"]] <-
    PercentageFeatureSet(UMI, pattern = "^MT-") +
    PercentageFeatureSet(UMI, pattern = "^Mt-") +
    PercentageFeatureSet(UMI, pattern = "^mt-")

  UMI[["RiboPercent"]] <-
    PercentageFeatureSet(UMI, pattern = "^Rps-") +
    PercentageFeatureSet(UMI, pattern = "^Rpl") +
    PercentageFeatureSet(UMI, pattern = "^RPS-") +
    PercentageFeatureSet(UMI, pattern = "^RPL")


  # ====================
  # 3. FILTER
  # ====================
  thresholds <- outlires(UMI@meta.data$nFeature_RNA)
  expect_s3_class(thresholds$plot, "ggplot")

  UMI <- subset(
    UMI,
    subset =
      nFeature_RNA > thresholds$thresholds[1] &
      nFeature_RNA <= thresholds$thresholds[length(thresholds$thresholds)] &
      MitoPercent < 20
  )


  # ====================
  # 4. PCA + DIM REDUCTION
  # ====================
  UMI <- FindVariableFeatures(
    UMI,
    selection.method = "vst",
    nfeatures = 2500,
    binning.method = "equal_frequency"
  )

  UMI <- ScaleData(UMI, features = rownames(UMI))

  UMI <- RunPCA(UMI, features = VariableFeatures(object = UMI))

  Elbow <- ElbowPlot(UMI, ndims = 50)

  # Seurat 5: stdev only in pca object
  dims <- data.frame(stdev = UMI[["pca"]]@stdev[1:50])
  colnames(dims)[1] <- "Elbow$data$stdev"

  dim <- dim_reuction_pcs(dims)
  expect_true(is.numeric(dim))


  # ====================
  # 5. JACKSTRAW + CLUSTERS + UMAP
  # ====================
  UMI <- JackStraw(UMI, num.replicate = 10, dims = dim)
  UMI <- ScoreJackStraw(UMI, dims = 1:dim)

  jc <- as.data.frame(UMI@reductions$pca@jackstraw@overall.p.values)
  jc <- jc[jc$Score < 0.05, ]
  dim <- as.vector(jc$PC)

  UMI <- FindNeighbors(UMI, dims = dim, reduction = "pca", nn.method = "rann")
  UMI <- FindClusters(UMI, resolution = 0.6, n.start = 10, n.iter = 1000)

  UMI <- RunUMAP(UMI, dims = dim)


  # ====================
  # 6. RECONSTRUCT CSSG PROJECT
  # ====================
  sc_project <- create_project_from_seurat(UMI)
  expect_gt(nrow(sc_project@matrices$norm), 0)


  # ====================
  # 7. CSSG PROCESSING
  # ====================
  sc_project <- get_cluster_stats(sc_project, type = "primary", only_pos = TRUE)
  expect_gt(nrow(sc_project@metadata$primary_markers), 0)

  sc_project <- namign_genes_selection(
    sc_project,
    type = "primary", top_n = 25,
    p_val = 0.05, select_stat = "p_val",
    mito_content = FALSE, ribo_content = FALSE
  )
  expect_gt(nrow(sc_project@metadata$naming_markers), 0)


  sc_project <- heterogeneity_select_specificity(
    sc_project = sc_project,
    type = "primary",
    heterogeneity_factor = 0.80,
    p_val = 0.05,
    max_genes = 1000,
    select_stat = "p_val",
    min_occ = 5,
    mito_content = FALSE
  )
  expect_gt(nrow(sc_project@metadata$heterogeneity_markers$heterogeneity_markers_df), 0)


  sc_project@metadata$heterogeneity_markers <- NULL


  sc_project <- heterogeneity_select_variance(
    sc_project = sc_project,
    heterogeneity_factor = 0.80,
    max_genes = 1000,
    min_occ = 5,
    min_exp = 0.1,
    rep_factor = 0.2,
    mito_content = FALSE
  )
  expect_gt(nrow(sc_project@metadata$heterogeneity_markers$heterogeneity_markers_df), 0)


  sc_project <- subclass_naming(
    sc_project = sc_project,
    class_markers = NULL,
    subclass_markers = NULL,
    species = "Homo sapiens",
    chunk_size = 5000
  )
  expect_true("subclass" %in% names(sc_project@names))


  # ====================
  # 8. SUBCLASS TESTS
  # ====================
  data <- bin_cell_test(
    p_val = 0.05,
    names = sc_project@names$subclass,
    min_cells = 20
  )
  expect_gt(nrow(data$data), 0)

  threshold <- cell_stat_graph(data$data)
  expect_s3_class(threshold, "ggplot")


  sc_project <- CSSG_markers(sc_project, max_combine = 1000, loss_val = 0.05)
  expect_gt(nrow(sc_project@metadata$cssg_markers), 0)


  # ====================
  # 9. SUBTYPES
  # ====================
  sc_project <- subtypes_naming(sc_project, species = "Homo sapiens")
  expect_true("subtypes" %in% names(sc_project@names))

  thr_data <- bin_cell_test(
    p_val = 0.05,
    names = sc_project@names$repaired,
    min_cells = 10
  )
  expect_gt(nrow(thr_data$data), 0)

  threshold <- cell_stat_graph(thr_data$data, include_ns = FALSE)
  expect_s3_class(threshold, "ggplot")

  select_list <- thr_data$data$names[
    thr_data$data$test %in% c("Good marked types", "Renamed")
  ]

  sub_sc_project <- subset_project(
    sc_project = sc_project,
    type = "subtypes",
    select_list = select_list
  )
  expect_true(length(sub_sc_project@names$subtypes) < length(sc_project@names$subtypes))


  # ====================
  # 10. SUBTYPE DATA
  # ====================
  data <- get_data(sub_sc_project, type = "subtypes", data = "norm")
  expect_gt(nrow(data), 0)

  data_avg <- get_avg_data(sub_sc_project, type = "subtypes", data = "norm")
  expect_gt(nrow(data_avg), 0)

  sub_sc_project <- get_cluster_stats(
    sc_project = sub_sc_project,
    type = "subtypes",
    only_pos = TRUE,
    min_pct = 0.05
  )
  expect_gt(nrow(sub_sc_project@metadata$subtypes_markers), 0)


  # ====================
  # 11. MARKER HEATMAPS
  # ====================
  markers <- get_names_markers(sub_sc_project, type = "subclasses")
  expect_true(is.vector(markers))

  markers <- get_names_markers(sub_sc_project, type = "subtypes")
  expect_true(is.vector(markers))

  plot1 <- marker_heatmap(
    sub_sc_project,
    type = "subtypes",
    markers = markers,
    angle_col = 270,
    fontsize_row = 7,
    fontsize_col = 7,
    font_labels = 8,
    clustering_method = "complete",
    x_axis = "Cells",
    y_axis = "Genes [log(CPM +1)]",
    scale = FALSE
  )
  expect_s3_class(plot1, "recordedplot")

  plot2 <- marker_heatmap(
    sub_sc_project,
    type = "subtypes",
    markers = markers,
    angle_col = 270,
    fontsize_row = 7,
    fontsize_col = 7,
    font_labels = 8,
    clustering_method = "complete",
    x_axis = "Cells",
    y_axis = "Scaled(Genes [log(CPM +1)])",
    scale = TRUE
  )
  expect_s3_class(plot2, "recordedplot")

})
