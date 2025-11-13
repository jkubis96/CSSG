library(CSSG.toolkit)
library(Seurat)
library(testthat)


sc_project <- create_project(sparse_matrix_path = "../../benchmark/data/", sparse_name = "matrix", rows_name = "genes", cols_name = "barcodes", type = "norm")

test_that("loading data", {
  expect_gt(nrow(sc_project@matrices$norm), 0)
})


colnames(sc_project@matrices$norm) <- make.unique(colnames(sc_project@matrices$norm))

UMI <- CreateSeuratObject(counts = sc_project@matrices$norm)

UMI@assays$RNA$data <- UMI@assays$RNA$counts

UMI[["MitoPercent"]] <- PercentageFeatureSet(UMI, pattern = "^MT-") + PercentageFeatureSet(UMI, pattern = "^Mt-") + PercentageFeatureSet(UMI, pattern = "^mt-")
UMI[["RiboPercent"]] <- PercentageFeatureSet(UMI, pattern = "^Rps-") + PercentageFeatureSet(UMI, pattern = "^Rpl") + PercentageFeatureSet(UMI, pattern = "^RPS-") + PercentageFeatureSet(UMI, pattern = "^RPL")


thresholds <- outlires(UMI@meta.data$nFeature_RNA)

test_that("threshold", {
  expect_s3_class(thresholds$plot, "ggplot")
})


UMI <- subset(UMI, subset = nFeature_RNA > thresholds$thresholds[1] & nFeature_RNA <= thresholds$thresholds[length(thresholds$thresholds)] &
  MitoPercent < 20)


UMI <- FindVariableFeatures(UMI, selection.method = "vst", nfeatures = 2500, binning.method = "equal_frequency")


UMI <- ScaleData(UMI, features = rownames(UMI))

UMI <- RunPCA(UMI, features = VariableFeatures(object = UMI))

Elbow <- ElbowPlot(UMI, ndims = 50)

dims <- as.data.frame(Elbow$data$stdev)


dim <- dim_reuction_pcs(dims)

test_that("dim_selection", {
  expect_true(is.numeric(dim))
})


UMI <- JackStraw(UMI, num.replicate = 10, dims = dim)
UMI <- ScoreJackStraw(UMI, dims = 1:dim)

jc <- as.data.frame(UMI@reductions$pca@jackstraw@overall.p.values)
jc <- jc[jc$Score < 0.05, ]
dim <- as.vector(jc$PC)

UMI <- FindNeighbors(UMI, dims = dim, reduction = "pca", nn.method = "rann")
UMI <- FindClusters(UMI, resolution = 0.6, n.start = 10, n.iter = 1000)


UMI <- RunUMAP(UMI, dims = dim)

sc_project <- create_project_from_seurat(UMI)

test_that("create project drom SEURAT", {
  expect_gt(nrow(sc_project@matrices$norm), 0)
})


sc_project <- get_cluster_stats(sc_project = sc_project, type = "primary", only_pos = TRUE)

test_that("cluster stats", {
  expect_gt(nrow(sc_project@metadata$primary_markers), 0)
})


sc_project <- namign_genes_selection(sc_project,
  type = "primary", top_n = 25,
  p_val = 0.05, select_stat = "p_val",
  mito_content = FALSE, ribo_content = FALSE
)


test_that("naming stats", {
  expect_gt(nrow(sc_project@metadata$naming_markers), 0)
})


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


test_that("naming CSSG stats", {
  expect_gt(nrow(sc_project@metadata$heterogeneity_markers$heterogeneity_markers_df), 0)
})


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


test_that("naming CSSG var", {
  expect_gt(nrow(sc_project@metadata$heterogeneity_markers$heterogeneity_markers_df), 0)
})


sc_project <- subclass_naming(sc_project = sc_project, class_markers = NULL, subclass_markers = NULL, species = "Homo sapiens", chunk_size = 5000)


test_that("subclasses naming", {
  expect_true("subclass" %in% names(sc_project@names))
})


data <- bin_cell_test(p_val = 0.05, names = sc_project@names$subclass, min_cells = 20)


test_that("bin stats subclasses", {
  expect_gt(nrow(data$data), 0)
})


threshold <- cell_stat_graph(data$data)


test_that("bin stats plot", {
  expect_s3_class(threshold, "ggplot")
})


sc_project <- CSSG_markers(sc_project = sc_project, max_combine = 1000, loss_val = 0.05)


test_that("CSSG markers", {
  expect_gt(nrow(sc_project@metadata$cssg_markers), 0)
})



sc_project <- subtypes_naming(sc_project = sc_project, markers_class = NULL, markers_subclass = NULL, species = "Homo sapiens")


test_that("subtypes naming", {
  expect_true("subtypes" %in% names(sc_project@names))
})


thr_data <- bin_cell_test(p_val = 0.05, names = sc_project@names$repaired, min_cells = 10)


test_that("bin stats subtypes", {
  expect_gt(nrow(thr_data$data), 0)
})


threshold <- cell_stat_graph(thr_data$data, include_ns = FALSE)


test_that("bin stats subtypes plot", {
  expect_s3_class(threshold, "ggplot")
})



select_list <- thr_data$data$names[thr_data$data$test %in% c("Good marked types", "Renamed")]

sub_sc_project <- subset_project(sc_project = sc_project, type = "subtypes", select_list = select_list)


test_that("project subsetting", {
  expect_true(length(sub_sc_project@names$subtypes) < length(sc_project@names$subtypes))
})


data <- get_data(sc_project = sub_sc_project, type = "subtypes", data = "norm")


test_that("get data", {
  expect_gt(nrow(data), 0)
})


data_avg <- get_avg_data(sc_project = sub_sc_project, type = "subtypes", data = "norm")

test_that("get average data", {
  expect_gt(nrow(data_avg), 0)
})


sub_sc_project <- get_cluster_stats(sc_project = sub_sc_project, type = "subtypes", only_pos = TRUE, min_pct = 0.05)

test_that("subtypes stats", {
  expect_gt(nrow(sub_sc_project@metadata$subtypes_markers), 0)
})


markers <- get_names_markers(sub_sc_project, type = "subclasses")

test_that("get merkers from names - subclases", {
  expect_true(is.vector(markers))
})


markers <- get_names_markers(sub_sc_project, type = "subtypes")

test_that("get merkers from names - subtypes", {
  expect_true(is.vector(markers))
})


plot <- marker_heatmap(sub_sc_project, type = "subtypes", markers = markers, angle_col = 270, fontsize_row = 7, fontsize_col = 7, font_labels = 8, clustering_method = "complete", x_axis = "Cells", y_axis = "Genes [log(CPM +1)]", scale = F)


test_that("heatmap norm", {
  expect_s3_class(plot, "recordedplot")
})


plot <- marker_heatmap(sub_sc_project, type = "subtypes", markers = markers, angle_col = 270, fontsize_row = 7, fontsize_col = 7, font_labels = 8, clustering_method = "complete", x_axis = "Cells", y_axis = "Scaled(Genes [log(CPM +1)])", scale = T)



test_that("heatmap scaled", {
  expect_s3_class(plot, "recordedplot")
})
