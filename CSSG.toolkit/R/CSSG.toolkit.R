#' @import grid parallel dplyr Matrix stringr ggplot2 textclean pheatmap gridExtra doParallel doSNOW foreach


firstup <- function(x) {
  set.seed(123)

  substr(x, 1, 1) <- toupper(substr(x, 1, 1))

  return(x)
}



cluster_naming <- function(matrix_a, markers) {
  set.seed(123)


  matrix_a <- as.data.frame(matrix_a)
  colname <- colnames(matrix_a)
  rownames(matrix_a) <- make.unique(toupper(rownames(matrix_a)), sep = "")

  if (length(markers) != 0) {
    list.markers <- c()
    for (l.marker in markers) {
      for (marker in l.marker) {
        TF <- (grepl("+", marker, fixed = T))
        if (TF == TRUE) {
          list.markers <- c(list.markers, textclean::mgsub(marker, c("+"), c("")))
        }
      }
    }


    index <- 0
    for (i in colnames(matrix_a)) {
      index <- index + 1
      rename_df <- as.data.frame(matrix_a[toupper(rownames(matrix_a)) %in% toupper(list.markers), index, drop = FALSE])
      rename_df <- rename_df[order(rename_df[[1]], decreasing = TRUE), , drop = FALSE]
      sum <- sum(rename_df)
      if (sum > 0) {
        colnames(matrix_a)[index] <- rownames(rename_df)[1]
      } else {
        colnames(matrix_a)[index] <- "Unknow"
      }
    }


    col <- 0
    for (c.marker in markers) {
      col <- col + 1
      for (marker in c.marker) {
        cell_n <- 0
        for (cell in colnames(matrix_a)) {
          cell_n <- cell_n + 1
          if (cell %in% textclean::mgsub(marker, c("+"), c(""))) {
            colnames(matrix_a)[cell_n] <- colnames(markers[col])
          }
        }
      }
    }
  } else {
    colnames(matrix_a) <- paste0("Cluster_", colnames(matrix_a))
  }

  return(as.vector(colnames(matrix_a)))
}





#' CSSG Cluster Naming
#'
#' Assigns names to clusters in a sc_project (sparse_matrix) object based on a CSSG data frame obtained using CSSG_markers() function.
#'
#' @param sparse_matrix A sparse matrix containing single-cell gene expression data.
#' @param CSSG_df A data frame containing CSSG information, including loss values and adjustment metrics.
#' @param species A character string specifying the species (default: 'Homo sapiens').
#' @return A character vector of cluster names.
#'
#' @examples
#' cssg_naming(sparse_matrix, CSSG_df)
#'
#' @export
cssg_naming <- function(sparse_matrix, CSSG_df, species = "Homo sapiens") {
  set.seed(123)

  cell_names <- c()
  for (cluster in sort(as.character(unique(colnames(sparse_matrix))))) {
    sub_cssg <- CSSG_df[CSSG_df$cluster %in% cluster, ]
    sub_cssg <- sub_cssg[sub_cssg$`loss_val` == min(sub_cssg$`loss_val`), ]
    sub_cssg <- sub_cssg[order(sub_cssg$`adj_hf`, decreasing = TRUE), ]
    tmp <- sparse_matrix[rownames(sparse_matrix) %in% gsub(" ", "", strsplit(rownames(sub_cssg)[1], split = " ")[[1]]), ]
    tmp <- as.matrix(tmp)

    for (i in 1:length(colnames(tmp))) {
      if (as.numeric(cluster) %in% colnames(tmp)[i] & colSums(tmp)[i] > 0) {
        tmp <- tmp[order(tmp[, i], decreasing = T), , drop = F]
        cell_names[i] <- firstup(tolower(rownames(tmp)[1]))
      } else if (as.numeric(cluster) %in% colnames(tmp)[i] & colSums(tmp)[i] == 0) {
        cell_names[i] <- "Bad!"
      }
    }
  }

  if (toupper(species) == "HOMO SAPIENS") {
    cell_names <- toupper(cell_names)
  } else {
    cell_names <- firstup(cell_names)
  }

  return(cell_names)
}



subcluster_naming <- function(average_expression, markers_subclass, cell_markers, species = "Homo sapiens") {
  set.seed(123)

  if (length(markers_subclass) != 0) {
    firstup <- function(x) {
      substr(x, 1, 1) <- toupper(substr(x, 1, 1))
      x
    }

    cell_names.1 <- c()
    for (i in 1:length(colnames(average_expression))) {
      rename_df <- average_expression[toupper(rownames(average_expression)) %in% toupper(markers_subclass), ]
      rename_df <- as.data.frame(rename_df[order(rename_df[, i], decreasing = T), , drop = F])
      if (sum(rename_df[, i]) > 0) {
        cell_names.1[i] <- toupper(rownames(rename_df[1, ]))
      } else {
        cell_names.1[i] <- ""
      }
    }


    cell_names.2 <- c()
    for (col in 1:length(colnames(average_expression))) {
      tmp_names <- cell_markers[as.character(cell_markers$cluster) %in% as.character(colnames(average_expression)[col]), ]
      tmp_names <- as.data.frame(tmp_names[order(tmp_names$pct_occurrence, decreasing = T), , drop = F])
      if (!toupper(cell_names.1[col]) %in% toupper(tmp_names$gene[1])) {
        cell_names.2[col] <- firstup(tolower(tmp_names$gene[1]))
      } else if (toupper(cell_names.1[col]) %in% toupper(tmp_names$gene[1])) {
        cell_names.2[col] <- firstup(tolower(tmp_names$gene[2]))
      }
    }
  } else if (length(markers_subclass) == 0) {
    firstup <- function(x) {
      substr(x, 1, 1) <- toupper(substr(x, 1, 1))
      x
    }


    cell_names.1 <- c()
    cell_names.2 <- c()
    for (col in 1:length(colnames(average_expression))) {
      tmp_names <- cell_markers[cell_markers$cluster %in% as.character(colnames(average_expression)[col]), ]
      tmp_names <- as.data.frame(tmp_names[order(tmp_names$pct_occurrence, decreasing = T), ])
      cell_names.1[col] <- toupper(tmp_names$gene[1])
      cell_names.2[col] <- firstup(tolower(tmp_names$gene[2]))
    }
  }


  if (toupper(species) == "HOMO SAPIENS") {
    cell_names.1 <- toupper(cell_names.1)
    cell_names.2 <- toupper(cell_names.2)
  } else {
    cell_names.1 <- firstup(tolower(cell_names.1))
    cell_names.2 <- firstup(tolower(cell_names.2))
  }

  sub_names <- paste(cell_names.1, cell_names.2)


  return(sub_names)
}





#' Subtypes Naming
#'
#' Assigns subtype names to clusters in a scRNA-seq project (sc_project) based on markers.
#'
#' @param sc_project An object containing scRNA-seq data, including expression matrices and metadata.
#' @param markers_class A list of class markers (default: NULL). If NULL the cluster will obtain native names from 'sc_project@names$primary' or 'sc_project@names$cluster'
#' @param markers_subclass A list of subclass markers (default: NULL). If NULL the markers will used from 'sc_project@metadata$heterogeneity_markers'
#' @param species A character string specifying the species (default: 'Homo sapiens'), which determines the structure of the names.
#' @param chunk_size Numeric. Size of chunks for processing large datasets (default: 5000).
#' @return The updated `sc_project` object with assigned subtypes.
#'
#' @examples
#' subtypes_naming(sc_project, markers_class, markers_subclass, species = "Homo sapiens", chunk_size = 500)
#'
#' @export
subtypes_naming <- function(sc_project, markers_class = NULL, markers_subclass = NULL, species = "Homo sapiens", chunk_size = 5000) {
  sparse_matrix <- sc_project@matrices$norm


  primary_names <- sc_project@names$subclass


  subtypes_markers <- sc_project@metadata$cssg_markers


  if (is.null(subtypes_markers) ||
    (is.atomic(subtypes_markers) && all(is.na(subtypes_markers))) ||
    (is.list(subtypes_markers) && length(subtypes_markers) == 0)) {
    stop("The value for 'subtypes_markers' is required and cannot be NA, NaN, NULL, or an empty list. Please use the 'CSSG_markers' and provide the value returned by these functions.")
  }

  if (length(primary_names) == length(colnames(sparse_matrix))) {
    cell_names <- cssg_naming(sparse_matrix = sparse_matrix, CSSG_df = subtypes_markers, species = species)

    names_to_return <- paste0(primary_names, " - ", cell_names)


    sc_project@names$subtypes <- names_to_return


    sc_project <- name_repairing(sc_project = sc_project, markers_class = markers_class, markers_subclass = markers_subclass, species = species, chunk_size = 5000)


    return(sc_project)
  }
}






#' Binary Cell Test
#'
#' Performs a binomial test to assess significance of cell type assignments.
#'
#' @param names A list or vector of renamed identifiers and marker data.
#' @param p_val Numeric. The significance threshold for p-values.
#' @param min_cells Numeric. Minimum number of cells required for significance (default: 10).
#' @return A list containing test results and categorized cell types.
#'
#' @examples
#' bin_cell_test(names, p_val = 0.05)
#'
#' @export
bin_cell_test <- function(names, p_val, min_cells = 10) {
  set.seed(123)


  if (is.atomic(names)) {
    subclass_names <- names
    bad <- subclass_names[grepl("BAD!", toupper(as.character(subclass_names)))]
    renamed.subnames <- c()
  } else if (is.list(names)) {
    subclass_names <- names$renamed_idents
    bad <- subclass_names[grepl("BAD!", toupper(as.character(subclass_names)))]

    renamed.subnames <- subclass_names[!subclass_names %in% names$original_names]
    renamed.subnames <- renamed.subnames[!as.character(renamed.subnames) %in% as.character(bad)]
  }


  new.subnames <- subclass_names[!as.character(subclass_names) %in% as.character(bad)]
  new.subnames <- new.subnames[!as.character(new.subnames) %in% as.character(renamed.subnames)]
  bad.subnames <- subclass_names[as.character(subclass_names) %in% as.character(bad)]
  subclass_names[subclass_names %in% bad.subnames] <- "Undefined"

  data <- as.data.frame(summary(as.factor(subclass_names), maxsum = length(unique(subclass_names))))
  colnames(data)[1] <- "n"
  data$names <- rownames(data)
  data$p_val <- NA


  for (n in 1:length(data$n)) {
    bin <- binom.test(data$n[n], ceiling(sum(data$n) / 2),
      p = median(data$n) / sum(data$n),
      conf.level = 0.9, alternative = "greater"
    )

    data$p_val[n] <- bin$p.value
  }


  below.names <- data$names[data$p_val > p_val]
  data$test[data$names %in% new.subnames] <- "Good marked types"
  data$test[data$names %in% renamed.subnames] <- "Renamed"
  data$test[data$names %in% below.names] <- "Non-significant"
  data$test[data$n < min_cells] <- "Non-significant"
  data$test[data$names %in% "Undefined"] <- "Undefined types"



  return(list("data" = data, "below.names" = below.names, "bad.subnames" = bad.subnames))
}





#' Cell Type Statistics Graph
#'
#' Creates a bar plot of cell type distributions and significance categories.
#'
#' @param data A data frame containing cell type counts and significance categories.
#' @param include_ns Logical. Whether to include non-significant cell types in the graph (default: TRUE).
#' @return A ggplot object representing the bar plot.
#'
#' @examples
#' cell_stat_graph(data, include_ns = TRUE)
#'
#' @export
cell_stat_graph <- function(data, include_ns = TRUE) {
  set.seed(123)

  if (include_ns == FALSE) {
    data <- data[!data$test %in% "Non-significant", , drop = FALSE]
  }

  threshold <- ggplot(data, aes(y = n, x = reorder(names, -n), fill = test, sort = test)) +
    geom_bar(stat = "identity") +
    xlab("Cells types") +
    ylab("Number of cells") +
    theme_bw() +
    theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")) +
    labs(fill = "Cells threshold") +
    coord_flip() +
    scale_fill_manual(values = c("Good marked types" = "#00BA38", "Renamed" = "#00BF7D", "Undefined types" = "#F8766D", "Non-significant" = "gray"))



  return(threshold)
}



name_repairing <- function(sc_project, markers_class, markers_subclass, species, chunk_size = 5000) {
  sparse_matrix <- sc_project@matrices$norm
  colnames(sparse_matrix) <- sc_project@names$subtypes

  cluster_heterogeneity_markers <- sc_project@metadata$naming_markers

  if (is.null(cluster_heterogeneity_markers) ||
    (is.atomic(cluster_heterogeneity_markers) && all(is.na(cluster_heterogeneity_markers))) ||
    (is.list(cluster_heterogeneity_markers) && length(cluster_heterogeneity_markers) == 0)) {
    stop("The value for 'cluster_heterogeneity_markers' is required and cannot be NA, NaN, NULL, or an empty list. Please use the 'heterogeneity_select_specificity' or 'heterogeneity_select_variance' function, and provide the '$heterogeneity_markers_df' value returned by these functions.")
  }

  #############################################################################
  # Class & Subclass naming



  if (!FALSE %in% unique(grepl("[^0-9]", colnames(sparse_matrix)))) {
    agg_subclasses <- aggregation_chr(sparse_matrix, chunk_size = chunk_size)
  } else {
    agg_subclasses <- aggregation_num(sparse_matrix, chunk_size = chunk_size)
  }




  for (n in names(sc_project@names)) {
    tmp_names <- sc_project@names[[n]]
    if (unique((unique(cluster_heterogeneity_markers$cluster) %in% tmp_names))) {
      new_names <- tmp_names
      break
    }
  }



  matching_df <- data.frame(
    old = sc_project@names$subtypes,
    new = new_names
  )


  matching_df <- distinct(matching_df)

  changed_names <- matching_df$new[match(colnames(agg_subclasses), matching_df$old)]

  colnames(agg_subclasses) <- changed_names

  clust_names <- cluster_naming(matrix_a = agg_subclasses, markers = markers_class)

  cell_names_2 <- subcluster_naming(agg_subclasses, markers_subclass, cluster_heterogeneity_markers, species)


  #############################################################################
  # Repair subclass_names

  clust_names <- paste(clust_names, cell_names_2)
  names_to_return <- new_names


  matching_df <- data.frame(
    old = as.character(colnames(agg_subclasses)),
    new = clust_names
  )

  matching_df <- distinct(matching_df)



  names_to_return <- matching_df$new[match(names_to_return, matching_df$old)]



  names_to_return <- paste(names_to_return, "-", gsub(".*- ", "", sc_project@names$subtypes))



  sc_project@names$repaired <- list("renamed_idents" = names_to_return, "original_names" = sc_project@names$subtypes)

  sc_project@names$subtypes <- names_to_return

  return(sc_project)
}




#' Aggregation of Gene Expression Data (Numeric Clustering)
#'
#' Aggregates gene expression data by grouping similar clusters based on numeric values.
#'
#' @param sparse_matrix A sparse matrix containing gene expression data.
#' @param chunk_size The maximum number of cells per processing batch (default: 5000).
#' @return A data frame containing aggregated gene expression values.
#'
#' @examples
#' aggregated_data <- aggregation_num(sparse_matrix)
#'
#' @export
aggregation_num <- function(sparse_matrix, chunk_size = 5000) {
  set.seed(123)

  subset_num <- round(length(colnames(sparse_matrix)) / chunk_size)
  cells_num <- round(length(colnames(sparse_matrix)) / subset_num)
  cols <- as.data.frame(summary(as.factor(colnames(sparse_matrix)), maxsum = length(unique(colnames(sparse_matrix)))))
  cols <- as.data.frame(cols[order(as.numeric(rownames(cols))), , drop = F])
  colnames(cols)[1] <- "n"

  if (round(length(colnames(sparse_matrix))) > (chunk_size * 1.5)) {
    for (batch in 1:subset_num) {
      if (batch == 1 & round(length(colnames(sparse_matrix))) > (chunk_size * 1.5)) {
        average_expression <- as.data.frame(as.matrix(sparse_matrix[, 1:cells_num]))
        average_expression <- t(average_expression)
        average_expression <- rowsum(average_expression, group = rownames(average_expression))
        average_expression <- t(average_expression)
      } else if (batch == subset_num & round(length(colnames(sparse_matrix))) > (chunk_size * 1.5)) {
        average_expression_tmp <- as.data.frame(as.matrix(sparse_matrix[, (((batch - 1) * cells_num) + 1):length(colnames(sparse_matrix))]))
        print((((batch - 1) * cells_num) + 1))
        print(length(colnames(sparse_matrix)))
        average_expression_tmp <- t(average_expression_tmp)
        average_expression_tmp <- rowsum(average_expression_tmp, group = rownames(average_expression_tmp))
        average_expression_tmp <- t(average_expression_tmp)
        average_expression <- cbind(average_expression, average_expression_tmp)
        rm(average_expression_tmp)
        average_expression <- as.data.frame(average_expression)
        average_expression <- t(average_expression)
        average_expression <- rowsum(average_expression, group = rownames(average_expression))
        average_expression <- t(average_expression)
        average_expression <- average_expression[, order(as.numeric(colnames(average_expression)))]
        average_expression <- as.data.frame(average_expression)

        for (col in row.names(cols)) {
          average_expression[, col] <- average_expression[, col] / cols$n[rownames(cols) %in% col]
        }
      } else if (batch > 1 & batch < subset_num & round(length(colnames(sparse_matrix))) > (chunk_size * 1.5)) {
        average_expression_tmp <- as.data.frame(as.matrix(sparse_matrix[, (((batch - 1) * cells_num) + 1):(batch * cells_num)]))
        print((((batch - 1) * cells_num) + 1))
        print(batch * cells_num)
        average_expression_tmp <- t(average_expression_tmp)
        average_expression_tmp <- rowsum(average_expression_tmp, group = rownames(average_expression_tmp))
        average_expression_tmp <- t(average_expression_tmp)
        average_expression <- cbind(average_expression, average_expression_tmp)
        rm(average_expression_tmp)
      }
    }
  }

  if (round(length(colnames(sparse_matrix))) <= (chunk_size * 1.5)) {
    average_expression <- as.data.frame(as.matrix(sparse_matrix))
    rm(sparse_matrix)
    average_expression <- t(average_expression)
    average_expression <- rowsum(average_expression, group = rownames(average_expression))
    average_expression <- t(average_expression)
    average_expression <- average_expression[, order(as.numeric(colnames(average_expression)))]
    average_expression <- as.data.frame(average_expression)


    for (col in row.names(cols)) {
      average_expression[, col] <- average_expression[, col] / cols$n[rownames(cols) %in% col]
    }
  }


  return(average_expression)
}




#' Aggregation of Gene Expression Data (Character Clustering)
#'
#' Aggregates gene expression data by grouping similar clusters based on character labels.
#'
#' @param sparse_matrix A sparse matrix containing gene expression data.
#' @param chunk_size The maximum number of cells per processing batch (default: 5000).
#' @return A data frame containing aggregated gene expression values.
#'
#' @examples
#' aggregated_data <- aggregation_chr(sparse_matrix)
#'
#' @export
aggregation_chr <- function(sparse_matrix, chunk_size = 5000) {
  set.seed(123)

  subset_num <- round(length(colnames(sparse_matrix)) / chunk_size)
  cells_num <- round(length(colnames(sparse_matrix)) / subset_num)
  cols <- as.data.frame(summary(as.factor(colnames(sparse_matrix)), maxsum = length(unique(colnames(sparse_matrix)))))
  cols <- as.data.frame(cols[order(rownames(cols)), , drop = F])
  colnames(cols)[1] <- "n"

  if (round(length(colnames(sparse_matrix))) > (chunk_size * 1.5)) {
    for (batch in 1:subset_num) {
      if (batch == 1 & round(length(colnames(sparse_matrix))) > (chunk_size * 1.5)) {
        average_expression <- as.data.frame(as.matrix(sparse_matrix[, 1:cells_num]))
        average_expression <- t(average_expression)
        average_expression <- rowsum(average_expression, group = rownames(average_expression))
        average_expression <- t(average_expression)
      } else if (batch == subset_num & round(length(colnames(sparse_matrix))) > (chunk_size * 1.5)) {
        average_expression_tmp <- as.data.frame(as.matrix(sparse_matrix[, (((batch - 1) * cells_num) + 1):length(colnames(sparse_matrix))]))
        print((((batch - 1) * cells_num) + 1))
        print(length(colnames(sparse_matrix)))
        average_expression_tmp <- t(average_expression_tmp)
        average_expression_tmp <- rowsum(average_expression_tmp, group = rownames(average_expression_tmp))
        average_expression_tmp <- t(average_expression_tmp)
        average_expression <- cbind(average_expression, average_expression_tmp)
        rm(average_expression_tmp)
        average_expression <- as.data.frame(average_expression)
        average_expression <- t(average_expression)
        average_expression <- rowsum(average_expression, group = rownames(average_expression))
        average_expression <- t(average_expression)
        average_expression <- average_expression[, order(colnames(average_expression))]
        average_expression <- as.data.frame(average_expression)

        for (col in row.names(cols)) {
          average_expression[, col] <- average_expression[, col] / cols$n[rownames(cols) %in% col]
        }
      } else if (batch > 1 & batch < subset_num & round(length(colnames(sparse_matrix))) > (chunk_size * 1.5)) {
        average_expression_tmp <- as.data.frame(as.matrix(sparse_matrix[, (((batch - 1) * cells_num) + 1):(batch * cells_num)]))
        print((((batch - 1) * cells_num) + 1))
        print(batch * cells_num)
        average_expression_tmp <- t(average_expression_tmp)
        average_expression_tmp <- rowsum(average_expression_tmp, group = rownames(average_expression_tmp))
        average_expression_tmp <- t(average_expression_tmp)
        average_expression <- cbind(average_expression, average_expression_tmp)
        rm(average_expression_tmp)
      }
    }
  }

  if (round(length(colnames(sparse_matrix))) <= (chunk_size * 1.5)) {
    average_expression <- as.data.frame(as.matrix(sparse_matrix))
    rm(sparse_matrix)
    average_expression <- t(average_expression)
    average_expression <- rowsum(average_expression, group = rownames(average_expression))
    average_expression <- t(average_expression)
    average_expression <- average_expression[, order(colnames(average_expression))]
    average_expression <- as.data.frame(average_expression)


    for (col in row.names(cols)) {
      average_expression[, col] <- average_expression[, col] / cols$n[rownames(cols) %in% col]
    }
  }


  return(average_expression)
}



#' Outlier Detection and Visualization
#'
#' Identifies outliers in a numeric dataset and visualizes their distribution using a histogram.
#'
#' @param input A numeric vector containing values to be analyzed.
#' @return A list containing a sorted vector of detected outliers and a ggplot2 histogram visualization.
#'
#' @examples
#' outlier_data <- outlires(c(1, 5, 10, 20, 30, 100, 150, 200))
#'
#' @export
outlires <- function(input) {
  set.seed(123)

  tmp_input <- input

  rang <- c()
  n <- c()

  iter <- ceiling(max(tmp_input) / 10)
  for (j in 1:iter) {
    rang <- c(rang, j * 10)
    n <- c(n, as.numeric(length(tmp_input[tmp_input <= j * 10])))
    tmp_input <- tmp_input[tmp_input > j * 10]
  }

  df <- data.frame(rang, n)

  factor <- 1.25
  bin <- c()
  for (k in 1:(nrow(df) - 1)) {
    if (df$n[k + 1] == 0) {
      bin <- c(bin, 0)
    } else if (df$n[k] > df$n[k + 1] * factor) {
      bin <- c(bin, -1)
    } else {
      bin <- c(bin, 1)
    }
  }

  bin <- c(bin[1], bin)

  df$bin <- bin

  df <- df[df$bin > 0, ]

  tf <- c()
  for (u in 1:(nrow(df) - 1)) {
    if (df$rang[u + 1] - df$rang[u] >= 100) {
      tf <- c(tf, FALSE)
    } else {
      tf <- c(tf, TRUE)
    }
  }

  df$tf <- c(tf[1], tf)

  df <- df[df$tf == TRUE, ]

  df <- df[df$n > quantile(df$n, 0.10), ]

  tmp_input <- input

  tmp_input <- tmp_input[tmp_input >= min(df$rang) & tmp_input <= max(df$rang)]

  library(ggplot2)

  min_value <- min(tmp_input)
  max_value <- max(tmp_input)
  percentile_10 <- quantile(tmp_input, 0.10)
  percentile_25 <- quantile(tmp_input, 0.25)
  percentile_33 <- quantile(tmp_input, 0.33)
  median_value <- median(tmp_input)
  percentile_66 <- quantile(tmp_input, 0.66)
  percentile_75 <- quantile(tmp_input, 0.75)
  percentile_90 <- quantile(tmp_input, 0.90)
  percentile_95 <- quantile(tmp_input, 0.95)


  plot <- ggplot(data.frame(tmp_input), aes(x = tmp_input)) +
    geom_histogram(bins = round(length(unique(tmp_input)) / 5, digits = 0), fill = "lightblue", color = "black") +
    labs(
      title = "Histogram of Features per Threshold",
      x = "Features Threshold",
      y = "Features per Threshold"
    ) +
    scale_x_continuous(
      breaks = c(min_value, percentile_10, percentile_25, percentile_33, median_value, percentile_66, percentile_75, percentile_90, percentile_95, max_value),
      labels = c(
        paste("Min:", min_value),
        paste("P10:", round(percentile_10, 2)),
        paste("P25:", round(percentile_25, 2)),
        paste("P33:", round(percentile_33, 2)),
        paste("Median:", round(median_value, 2)),
        paste("P66:", round(percentile_66, 2)),
        paste("P75:", round(percentile_75, 2)),
        paste("P90:", round(percentile_90, 2)),
        paste("P95:", round(percentile_95, 2)),
        paste("Max:", max_value)
      )
    ) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 9)) # Rotate x-axis labels by 27 degrees


  return(list("thresholds" = sort(unique(tmp_input)), "plot" = plot))
}








#' Gene Variance Calculation
#'
#' Calculates the variance, mean, and percentage occurrence of gene expressions across samples.
#'
#' @param data A matrix of gene expression values where rows represent genes and columns represent samples.
#' @param min A numeric value specifying the threshold for counting occurrences (default: 0.1).
#' @return A data frame containing gene names, variance, mean expression, and percentage occurrence.
#'
#' @examples
#' variance_data <- genes_variance_calculate(expression_matrix)
#'
#' @export
genes_variance_calculate <- function(data, min = 0.1) {
  set.seed(123)
  `%>%` <- dplyr::`%>%`
  gene_variance <- apply(data, 1, var)
  gene_mean <- apply(data, 1, mean)
  pct_occurrence <- apply(data, 1, function(row) {
    mean(row >
      min) * 100
  })
  variance_df <- data.frame(
    gene = rownames(data), variance = gene_variance,
    avg = gene_mean, pct_occurrence = pct_occurrence
  )
  variance_df <- variance_df %>% dplyr::arrange(desc(variance))
  return(variance_df)
}





#' Cluster Marker Gene Identification
#'
#' Identifies marker genes for specified cluster types in a single-cell dataset by comparing
#' gene expression across clusters.
#'
#' @param sc_project A single-cell project object containing normalized gene expression data.
#' @param type A character string specifying the cluster type. Must be one of: 'subtypes', 'subclasses', 'cluster', or 'primary'.
#' @param only_pos A logical value indicating whether to retain only positively differentially expressed genes (TRUE) or include both upregulated and downregulated genes (FALSE). Default is TRUE.
#' @param min_pct A numeric value specifying the minimum percentage of cells in a cluster that must express a gene for it to be considered. Default is 0.05.
#' @return An updated `sc_project` object with marker gene statistics stored in the corresponding
#'   metadata slot.
#'
#' @details The function performs the following steps:
#'   - Validates the `type` parameter.
#'   - Extracts the relevant cluster labels.
#'   - Iterates over each cluster to compute:
#'     - The proportion of cells expressing each gene.
#'     - The log fold change (logFC) between the target cluster and the rest.
#'     - Statistical significance (Wilcoxon test p-values).
#'   - Saves results in the appropriate metadata slot of `sc_project`.
#'
#' @examples
#' updated_project <- get_cluster_stats(sc_project, type = "subtypes", only_pos = TRUE, min_pct = 0.1)
#'
#' @export
get_cluster_stats <- function(sc_project, type = NaN, only_pos = TRUE, min_pct = 0.1) {
  set.seed(123)
  if (!type %in% c("subtypes", "subclasses", "cluster", "primary")) {
    stop("Invalid `type` provided. The `type` parameter should be either 'subtypes' or 'subclasses' or 'cluster' or 'primary'")
  }
  if (type == "subtypes") {
    slot <- "subtypes_markers"
    names <- sc_project@names$subtypes
  } else if (type == "subclasses") {
    slot <- "subclasses_markers"
    names <- sc_project@names$subclass
  } else if (type == "cluster") {
    slot <- "clusters_markers"
    names <- sc_project@names$cluster
  } else if (type == "primary") {
    slot <- "primary_markers"
    names <- sc_project@names$primary
  }
  data <- sc_project@matrices$norm
  colnames(data) <- names
  if (FALSE %in% unique(grepl("[^0-9]", colnames(data)))) {
    cluster_group <- sort(as.numeric(as.character(unique(colnames(data)))))
  } else {
    cluster_group <- as.character(unique(colnames(data)))
  }
  results <- data.frame()
  for (c in cluster_group) {
    cat(paste("\n\n Cluster ->  ", c, "- searching marker genes... \n\n"))
    tmp1 <- data[, colnames(data) %in% c]
    tmp_results <- data.frame(genes = rownames(tmp1))
    tmp2 <- data[, !colnames(data) %in% c]
    tmp_sum <- tmp1
    tmp_sum[tmp_sum > 0] <- 1
    tmp_results$pct_occurrence <- rowSums(tmp_sum) / ncol(tmp_sum)
    rm(tmp_sum)
    t1 <- rowMeans(tmp1)
    t2 <- rowMeans(tmp2)
    v1 <- apply(tmp1, 1, var, na.rm = TRUE)
    v2 <- apply(tmp2, 1, var, na.rm = TRUE)
    n1 <- ncol(tmp1)
    n2 <- ncol(tmp2)

    s_pooled <- sqrt(((n1 - 1) * v1 + (n2 - 1) * v2) / (n1 + n2 - 2))
    tmp_results$esm <- (t1 - t2) / s_pooled

    tmp_results$avg_logFC <- log2((t1 + (min(t1[t1 > 0]) / 2)) / (t2 +
      (min(t2[t2 > 0]) / 2)))
    tmp_results <- tmp_results[tmp_results$pct_occurrence >
      min_pct, , drop = FALSE]
    if (only_pos == TRUE) {
      tmp_results <- tmp_results[tmp_results$avg_logFC >
        0, ]
      tmp_data <- data[match(tmp_results$genes, rownames(data)), ]
      tmp_results$p_val <- apply(tmp_data, 1, function(gene_expression) {
        wilcox.test(gene_expression[colnames(tmp_data) %in%
          c], gene_expression[!colnames(tmp_data) %in%
          c])$p.value
      })
    } else {
      tmp_results$p_val <- apply(data, 1, function(gene_expression) {
        wilcox.test(gene_expression[colnames(data) %in%
          c], gene_expression[!colnames(data) %in% c])$p.value
      })
    }
    tmp_results$cluster <- c
    results <- rbind(results, tmp_results)
  }

  sc_project@metadata[[slot]] <- results

  return(sc_project)
}



#' Select Cluster-Specific Marker Genes Based on Heterogeneity
#'
#' Identifies cluster-specific marker genes by considering gene expression heterogeneity,
#' statistical significance, and occurrence thresholds.
#'
#' @param sc_project A single-cell project object containing normalized gene expression data.
#' @param type A character string specifying the cluster type. Must be one of:
#'   'subtypes', 'subclasses', 'cluster', or 'primary'.
#' @param heterogeneity_factor A numeric value (range 0-1) defining the maximum percentage of cells in a
#'   cluster expressing a gene for it to be considered heterogeneous. Default is 0.8.
#' @param p_val A numeric value specifying the p-value threshold for gene selection.
#'   Default is 0.05.
#' @param max_genes A numeric value defining the maximum number of genes to retain per cluster.
#'   Default is 1000.
#' @param select_stat A character string specifying the ranking statistic for gene selection.
#'   Options: 'p_val' (default), 'avg_logFC', or 'avg_log2FC'.
#' @param min_occ A numeric value specifying the minimum percentage occurrence of a gene in a cluster.
#'   Default is 5.
#' @param mito_content A logical value indicating whether to include mitochondrial genes in the selection.
#'   Default is FALSE.
#'
#' @return An updated `sc_project` object with identified heterogeneity marker genes stored
#'   in the `metadata$heterogeneity_markers` slot.
#'
#' @details The function performs the following steps:
#'   - Validates the `type` parameter.
#'   - Extracts the relevant cluster labels.
#'   - Filters marker genes based on p-value and removes mitochondrial genes if specified.
#'   - Iterates over each cluster to:
#'     - Identify genes expressed in a heterogeneous manner.
#'     - Rank genes based on `select_stat` (e.g., log fold change or p-value).
#'     - Retain a specified number of top marker genes.
#'   - Removes genes appearing in multiple clusters beyond a threshold.
#'   - Stores the final marker gene lists in `sc_project@metadata$heterogeneity_markers`.
#'
#' @examples
#' updated_project <- heterogeneity_select_specificity(sc_project,
#'   type = "subtypes",
#'   heterogeneity_factor = 0.75, p_val = 0.01, max_genes = 500, select_stat = "avg_logFC"
#' )
#'
#' @export
heterogeneity_select_specificity <- function(sc_project, type = NaN, heterogeneity_factor = 80, p_val = 0.05, max_genes = 1000, select_stat = "p_val", min_occ = 5, mito_content = FALSE) {
  set.seed(123)



  if (!type %in% c("subtypes", "subclasses", "cluster", "primary")) {
    stop("Invalid `type` provided. The `type` parameter should be either 'subtypes' or 'subclasses' or 'cluster' or 'primary'")
  }


  if (type == "subtypes") {
    slot <- "subtypes_markers"

    names <- sc_project@names$subtypes
  } else if (type == "subclasses") {
    slot <- "subclasses_markers"

    names <- sc_project@names$subclass
  } else if (type == "cluster") {
    slot <- "clusters_markers"

    names <- sc_project@names$cluster
  } else if (type == "primary") {
    slot <- "primary_markers"

    names <- sc_project@names$primary
  }



  sparse_matrix <- sc_project@matrices$norm

  colnames(sparse_matrix) <- names

  marker_df <- sc_project@metadata[[slot]]

  if (FALSE %in% unique(grepl("[^0-9]", colnames(sparse_matrix)))) {
    cluster_group <- sort(as.numeric(as.character(unique(colnames(sparse_matrix)))))
  } else {
    cluster_group <- as.character(unique(colnames(sparse_matrix)))
  }

  marker_df <- marker_df[marker_df$p_val < p_val, ]

  if (mito_content == FALSE) {
    marker_df <- marker_df[!marker_df$gene %in% marker_df$gene[grepl("MT-", toupper(marker_df$gene))], , drop = FALSE]
    marker_df <- marker_df[!marker_df$gene %in% marker_df$gene[grepl("MT.", toupper(marker_df$gene))], , drop = FALSE]
  }



  for (cluster in cluster_group) {
    cat(paste("\n\n Cluster ", cluster, "- searching marker genes... \n\n"))

    gen_cor <- unique(marker_df$gene[marker_df$cluster %in% cluster])


    tmp2 <- sparse_matrix[, colnames(sparse_matrix) %in% cluster, drop = FALSE]

    tmp2 <- tmp2[toupper(rownames(tmp2)) %in% toupper(gen_cor), , drop = FALSE]


    # expression level

    rmean <- rowMeans(tmp2)
    rmean <- rmean[rmean >= quantile(rmean, 0.1)]


    if (length(rmean) > 10) {
      tmp2 <- tmp2[toupper(rownames(tmp2)) %in% toupper(names(rmean)), , drop = FALSE]
    }


    tmp2[tmp2 > 0] <- 1L
    tmp2[tmp2 == 0] <- 0L

    perc <- (rowSums(tmp2 == 1) / ncol(tmp2)) * 100

    perc2 <- perc[perc <= heterogeneity_factor & perc >= min_occ]


    marker_df_CSSG <- marker_df[marker_df$cluster %in% cluster, ]


    if (length(perc2) > 10) {
      marker_df_CSSG <- marker_df_CSSG[marker_df_CSSG$gene %in% names(perc2), ]
    } else {
      marker_df_CSSG <- marker_df_CSSG[marker_df_CSSG$gene %in% names(perc), ]
    }


    if (exists("heterogeneity_markers_df") == FALSE) {
      heterogeneity_markers_df <- data.frame(gene = names(perc), pct_occurrence = perc)
      heterogeneity_markers_df$cluster <- cluster
    } else {
      tmp_het <- data.frame(gene = names(perc), pct_occurrence = perc)
      tmp_het$cluster <- cluster
      heterogeneity_markers_df <- rbind(heterogeneity_markers_df, tmp_het)
    }




    if (exists("marker_df_tmp") == FALSE) {
      if (select_stat %in% c("avg_logFC", "avg_log2FC")) {
        marker_df_tmp <- marker_df_CSSG[order(marker_df_CSSG$avg_logFC, decreasing = TRUE), ]
        marker_df_tmp <- marker_df_tmp[1:as.numeric(max_genes), ]
        marker_df_tmp <- marker_df_tmp[!is.na(marker_df_tmp$cluster), ]
      } else {
        marker_df_tmp <- marker_df_CSSG[order(marker_df_CSSG$p_val, decreasing = FALSE), ]
        marker_df_tmp <- marker_df_tmp[1:as.numeric(max_genes), ]
        marker_df_tmp <- marker_df_tmp[!is.na(marker_df_tmp$cluster), ]
      }
    } else {
      if (select_stat %in% c("avg_logFC", "avg_log2FC")) {
        marker_df_tmp1 <- marker_df_CSSG[order(marker_df_CSSG$avg_logFC, decreasing = TRUE), ]
        marker_df_tmp1 <- marker_df_tmp1[1:as.numeric(max_genes), ]
        marker_df_tmp1 <- marker_df_tmp1[!is.na(marker_df_tmp1$cluster), ]
      } else {
        marker_df_tmp1 <- marker_df_CSSG[order(marker_df_CSSG$p_val, decreasing = FALSE), ]
        marker_df_tmp1 <- marker_df_tmp1[1:as.numeric(max_genes), ]
        marker_df_tmp1 <- marker_df_tmp1[!is.na(marker_df_tmp1$cluster), ]
      }

      marker_df_tmp <- rbind(marker_df_tmp, marker_df_tmp1)
    }
  }


  heterogeneity_markers_df_2 <- data.frame()


  value_counts <- table(heterogeneity_markers_df$gene[heterogeneity_markers_df$pct_occurrence > 55])
  value_counts <- data.frame(
    valuename = names(value_counts),
    n = as.numeric(value_counts)
  )

  to_rm <- value_counts$valuename[value_counts$n > 1]

  for (c in unique(heterogeneity_markers_df$cluster)) {
    tmp <- heterogeneity_markers_df[heterogeneity_markers_df$cluster %in% c, , drop = FALSE]
    tmp <- tmp[!tmp$gene %in% to_rm, , drop = FALSE]
    if (nrow(tmp) < 2) {
      tmp <- heterogeneity_markers_df[heterogeneity_markers_df$cluster %in% c, , drop = FALSE]
    }

    heterogeneity_markers_df_2 <- rbind(heterogeneity_markers_df_2, tmp)
  }


  heterogeneity_markers_df_2 <- heterogeneity_markers_df_2 %>%
    dplyr::group_by(cluster) %>%
    dplyr::top_n(n = max_genes, wt = pct_occurrence)


  sc_project@metadata$heterogeneity_markers <- list("marker_df" = marker_df_tmp, "heterogeneity_markers_df" = heterogeneity_markers_df_2)


  return(sc_project)
}




#' Select Marker Genes Based on Expression Variance and Heterogeneity
#'
#' Identifies cluster-specific marker genes using expression variance, occurrence thresholds,
#' and heterogeneity factors.
#'
#' @param sc_project A single-cell project object containing normalized gene expression data.
#' @param heterogeneity_factor A numeric value (range 0-1) defining the maximum percentage of cells in a
#'   cluster expressing a gene for it to be considered heterogeneous. Default is 0.8.
#' @param max_genes A numeric value defining the maximum number of genes to retain per cluster.
#'   Default is 1000.
#' @param min_occ A numeric value specifying the minimum percentage occurrence of a gene in a cluster.
#'   Default is 5.
#' @param min_exp A numeric value specifying the minimum expression level threshold for considering
#'   a gene in the variance calculation. Default is 0.1.
#' @param rep_factor A numeric value defining the replication threshold for filtering genes appearing
#'   in multiple clusters. Default is 0.5.
#' @param mito_content A logical value indicating whether to include mitochondrial genes in the selection.
#'   Default is FALSE.
#'
#' @return An updated `sc_project` object with selected heterogeneity marker genes stored in
#'   `metadata$heterogeneity_markers`.
#'
#' @details The function performs the following steps:
#'   - Extracts and processes normalized expression data.
#'   - Calculates gene variance and occurrence using `genes_variance_calculate()`.
#'   - Removes mitochondrial genes if `mito_content = FALSE`.
#'   - Filters genes based on heterogeneity constraints (`min_occ`, `heterogeneity_factor`).
#'   - Removes genes appearing in too many clusters based on `rep_factor`.
#'   - Selects the top `max_genes` per cluster based on expression variance.
#'   - Stores the final gene list in `sc_project@metadata$heterogeneity_markers`.
#'
#' @examples
#' updated_project <- heterogeneity_select_variance(sc_project,
#'   heterogeneity_factor = 0.5,
#'   max_genes = 200, min_occ = 10, min_exp = 0.2, rep_factor = 0.4
#' )
#'
#' @export
heterogeneity_select_variance <- function(sc_project,
                                          heterogeneity_factor = 0.8,
                                          max_genes = 1000,
                                          min_occ = 5,
                                          min_exp = 0.1,
                                          rep_factor = 0.1,
                                          mito_content = FALSE) {
  set.seed(123)


  sparse_matrix <- sc_project@matrices$norm


  if (FALSE %in% unique(grepl("[^0-9]", colnames(sparse_matrix)))) {
    cluster_group <- sort(as.numeric(as.character(unique(colnames(sparse_matrix)))))
  } else {
    cluster_group <- as.character(unique(colnames(sparse_matrix)))
  }


  full_data <- data.frame()
  heterogeneity_markers_df <- data.frame()


  for (cluster in cluster_group) {
    cat(paste("\n\n Cluster ", cluster, "- searching marker genes... \n\n"))


    tmp2 <- as.data.frame(as.matrix(sparse_matrix[, colnames(sparse_matrix) %in% cluster, drop = FALSE]))

    var_data <- genes_variance_calculate(tmp2, min = min_exp)


    if (mito_content == FALSE) {
      var_data <- var_data[!var_data$gene %in% var_data$gene[grepl("MT-", toupper(var_data$gene))], , drop = FALSE]
      var_data <- var_data[!var_data$gene %in% var_data$gene[grepl("MT.", toupper(var_data$gene))], , drop = FALSE]
    }


    var_data$cluster <- cluster

    heterogeneity_markers_df <- rbind(heterogeneity_markers_df, var_data)


    var_data_tmp <- var_data[var_data$pct_occurrence > min_occ & var_data$pct_occurrence < heterogeneity_factor * 100, , drop = FALSE]

    if (nrow(var_data_tmp) > 10) {
      var_data <- var_data_tmp
    }


    full_data <- rbind(full_data, var_data)
  }




  value_counts <- table(full_data$gene)
  value_counts <- data.frame(
    valuename = names(value_counts),
    n = as.numeric(value_counts)
  )


  to_rm <- value_counts$valuename[value_counts$n <= ceiling(length(unique(full_data$cluster)) * rep_factor)]


  full_2 <- data.frame()


  for (c in unique(full_data$cluster)) {
    tmp <- full_data[full_data$cluster %in% c, , drop = FALSE]
    tmp <- tmp[!tmp$gene %in% to_rm, , drop = FALSE]
    if (nrow(tmp) < 10) {
      tmp <- full_data[full_data$cluster %in% c, , drop = FALSE]
    }

    full_2 <- rbind(full_2, tmp)
  }

  full_2 <- full_2 %>%
    dplyr::group_by(cluster) %>%
    dplyr::top_n(n = max_genes, wt = avg)


  sparse_matrix_bin <- sparse_matrix[unique(full_2$gene), ]
  sparse_matrix_bin[sparse_matrix_bin < min_exp] <- 0
  sparse_matrix_bin[sparse_matrix_bin >= min_exp] <- 1

  row_percent <- data.frame(
    gene = rownames(sparse_matrix_bin),
    percent = rowSums(sparse_matrix_bin) / ncol(sparse_matrix_bin) * 100
  )


  to_rm <- row_percent$gene[row_percent$percent <= 30]


  full_3 <- data.frame()


  for (c in unique(full_2$cluster)) {
    tmp <- full_2[full_2$cluster %in% c, , drop = FALSE]
    tmp <- tmp[!tmp$gene %in% to_rm, , drop = FALSE]
    if (nrow(tmp) < 10) {
      tmp <- full_2[full_2$cluster %in% c, , drop = FALSE]
    }

    full_3 <- rbind(full_3, tmp)
  }



  sc_project@metadata$heterogeneity_markers <- list("marker_df" = full_data, "heterogeneity_markers_df" = full_3)


  return(sc_project)
}







#' CSSG Marker Optimization
#'
#' Optimizes marker gene selection across clusters by iteratively combining genes and evaluating
#' heterogeneity and loss metrics. Utilizes parallel processing to improve computational efficiency.
#'
#' @param sc_project An object containing scRNA-seq data, including normalized expression matrices.
#' @param max_combine The maximum number of gene combinations to evaluate per cluster (default: 1000).
#' @param loss_val The maximum allowed proportion of cells in which a marker is absent (default: 0.05).
#'
#' @return The modified `sc_project` object with added CSSG marker metadata.
#'
#' @examples
#' sc_project <- CSSG_markers(sc_project, max_combine = 1000, loss_val = 0.05)
#'
#' @export
CSSG_markers <- function(sc_project, max_combine = 1000, loss_val = 0.05) {
  set.seed(123)

  sparse_matrix <- sc_project@matrices$norm
  markers <- sc_project@metadata$heterogeneity_markers$heterogeneity_markers_df

  old_warn <- options("warn") # Save the current warning level
  options(warn = -1)

  add_avg_val <- function(x, tmp2) {
    av <- mean(tmp2[x, ])
    return(av)
  }


  check_memory_par <- function(par_object, per_mem = .8) {
    os_type <- Sys.info()["sysname"]



    if (os_type == "Windows") {
      get_memory_info_windows <- function() {
        total_memory_mb <- as.numeric(system2("powershell",
          args = "-Command \"(Get-CimInstance -ClassName Win32_ComputerSystem).TotalPhysicalMemory / 1MB\"", stdout = TRUE
        ))

        free_memory_mb <- as.numeric(system2("powershell",
          args = "-Command \"(Get-CimInstance -ClassName Win32_OperatingSystem).FreePhysicalMemory / 1\"", stdout = TRUE
        ))

        used_memory_mb <- total_memory_mb - free_memory_mb

        list(
          total_memory_mb = total_memory_mb,
          free_memory_mb = free_memory_mb,
          used_memory_mb = used_memory_mb
        )
      }

      memory_info <- get_memory_info_windows()

      current_memory <- memory_info$used_memory_mb
      max_memory <- memory_info$total_memory_mb
      avaiable <- memory_info$free_memory_mb

      object_size_bytes <- 0

      for (obj in par_object) {
        object_size_bytes <- object_size_bytes + object.size(obj)
      }


      object_size_mb <- as.numeric(object_size_bytes / (1024^2)) * 2

      av_cpu <- (avaiable * per_mem) / object_size_mb

      CPU <- detectCores()

      if (av_cpu >= CPU) {
        av_cpu <- CPU - 2

        if (av_cpu < 1) {
          av_cpu <- 1
        }
      } else if (av_cpu < 1) {
        av_cpu <- 1
      }
    } else if (os_type == "Linux") {
      available_memory_bytes <- as.numeric(system("awk '/MemAvailable:/ {print $2 * 1024}' /proc/meminfo", intern = TRUE))
      avaiable <- available_memory_bytes / (1024^2)


      object_size_bytes <- 0

      for (obj in par_object) {
        object_size_bytes <- object_size_bytes + object.size(par_object)
      }


      object_size_mb <- as.numeric(object_size_bytes / (1024^2)) * 2

      av_cpu <- (avaiable * per_mem) / object_size_mb

      CPU <- detectCores()

      if (av_cpu >= CPU) {
        av_cpu <- CPU - 2
      } else if (av_cpu < 1) {
        av_cpu <- 1
      }
    } else if (os_type == "Darwin") {
      free_pages <- as.numeric(system("vm_stat | grep 'Pages free:' | awk '{print $3}' | sed 's/[^0-9]//g'", intern = TRUE))

      page_size <- as.numeric(system("sysctl -n hw.pagesize", intern = TRUE))

      available_memory_bytes <- free_pages * page_size

      avaiable <- available_memory_bytes / (1024^2)

      object_size_bytes <- 0

      for (obj in par_object) {
        object_size_bytes <- object_size_bytes + object.size(par_object)
      }


      object_size_mb <- as.numeric(object_size_bytes / (1024^2)) * 2

      av_cpu <- ((avaiable * per_mem) / object_size_mb) - 1

      CPU <- detectCores()


      if (av_cpu >= CPU) {
        av_cpu <- CPU - 2
      } else if (av_cpu <= 1) {
        av_cpu <- 1
      }
    } else {
      av_cpu <- detectCores() - 2
      print("Unknown or unsupported operating system.")
    }

    print(paste("Number of CPU cores allowed for parallel processing:", av_cpu))

    return(av_cpu)
  }



  cat("\n\n The CSSG start \n it can last several minutes depending on the number of clusters and set parameters \n\n\n")


  if (FALSE %in% unique(grepl("[^0-9]", colnames(sparse_matrix)))) {
    cluster_group <- sort(as.numeric(as.character(unique(colnames(sparse_matrix)))))
  } else {
    cluster_group <- as.character(unique(colnames(sparse_matrix)))
  }

  complete_df <- data.frame()

  for (cluster in cluster_group) {
    cat(paste("\n\n Cluster ", cluster, "- searching heterogeneity marker genes... \n\n"))


    # start data storing

    approved_df <- data.frame()


    # initisl val
    tmp_recursive <- NaN


    gen_cor <- unique(markers$gene[markers$cluster %in% cluster])


    tmp2 <- sparse_matrix[, colnames(sparse_matrix) %in% cluster, drop = FALSE]
    tmp2 <- tmp2[toupper(rownames(tmp2)) %in% toupper(gen_cor), , drop = FALSE]



    tmp3 <- tmp2
    tmp3[tmp3 > 0] <- 1L
    tmp3[tmp3 == 0] <- 0L



    tmp2 <- data.frame(row.names = rownames(tmp2), avg = rowMeans(tmp2))


    # First loop for duble combination

    res_df <- tmp3

    ncol_res_df <- ncol(res_df)
    perc0 <- rowSums(res_df == 0) / ncol_res_df
    perc1 <- rowSums(res_df == 1) / ncol_res_df
    avg_exp <- unlist(lapply(strsplit(rownames(res_df), split = " "), function(x) add_avg_val(x, tmp2)))
    multi <- 1 - (perc0 + perc1)

    # Combine perc0 and perc1 into a sparse matrix (optional)
    last_df <- as.data.frame(cbind(perc0, perc1, avg_exp, multi))

    # Filter based on perc1 quantile

    last_df <- last_df[last_df$perc1 > quantile(last_df$perc1, 0.6), , drop = FALSE]


    # first results


    approved_df <- rbind(approved_df, last_df[last_df$perc0 <= loss_val, , drop = FALSE])

    if (nrow(approved_df) > 0) {
      approved_df_tmp <- approved_df
      approved_df_tmp$het <- (1 - ((approved_df_tmp$perc1 + ((1 - (approved_df_tmp$perc1 + approved_df_tmp$perc0)) * 1.25)) / (str_count(string = rownames(approved_df_tmp), pattern = " ") + 1)))
      approved_df_tmp$het_adj <- (1 - ((approved_df_tmp$perc1 + ((1 - (approved_df_tmp$perc1 + approved_df_tmp$perc0)) * 1.25)) / (str_count(string = rownames(approved_df_tmp), pattern = " ") + 1))) - (approved_df_tmp$perc0 * 2)
      approved_df_tmp$cluster <- cluster

      complete_df <- rbind(approved_df_tmp, complete_df)

      rm(approved_df_tmp)
    }




    while (TRUE) {
      ls <- ls(all.names = TRUE)

      CPU <- check_memory_par(par_object = c("res_df", "tmp3", "tmp2", "add_avg_val"))
      cl <- makeCluster(CPU)
      registerDoParallel(cl)
      registerDoSNOW(cl)


      pb <- txtProgressBar(max = nrow(res_df), style = 3)
      progress <- function(n) setTxtProgressBar(pb, n)
      opts <- list(progress = progress)


      res_df_multi <- foreach(i = 1:nrow(res_df), .options.snow = opts, .combine = rbind, .inorder = FALSE, .packages = c("Matrix", "stringr"), .export = c("res_df", "tmp3", "tmp2", "add_avg_val"), .noexport = ls) %dopar% {
        exclude_rows <- strsplit(rownames(res_df)[i], split = " ")[[1]]
        add_tmp <- as.matrix(tmp3[!rownames(tmp3) %in% exclude_rows, , drop = FALSE])

        if (nrow(add_tmp) != 0) {
          row_vector <- paste(rownames(add_tmp), rownames(res_df)[rep(i, nrow(add_tmp))])

          row_vector <- strsplit(row_vector, split = " ")
          row_vector <- lapply(row_vector, sort)
          row_vector <- lapply(row_vector, function(x) paste(x, collapse = " "))


          # Perform the matrix addition

          results_tmp <- add_tmp + as.matrix(res_df[rep(i, nrow(add_tmp)), , drop = FALSE])
          rownames(results_tmp) <- row_vector
        } else {
          results_tmp <- as.matrix(res_df)
        }

        # top1
        #
        ncol_results_tmp <- ncol(results_tmp)
        perc0 <- rowSums(results_tmp == 0) / ncol_results_tmp
        perc1 <- rowSums(results_tmp == 1) / ncol_results_tmp
        avg_exp <- unlist(lapply(strsplit(rownames(results_tmp), split = " "), function(x) add_avg_val(x, tmp2)))
        multi <- 1 - (perc0 + perc1)

        sorting_df <- data.frame(perc1, avg_exp)

        results_tmp <- cbind(results_tmp, perc0, perc1, avg_exp, multi)


        rm(perc0, perc1, avg_exp, multi)


        results_tmp <- results_tmp[order(sorting_df$perc1, sorting_df$avg_exp, decreasing = TRUE), , drop = FALSE]


        if (nrow(add_tmp) != 0) {
          results_tmp <- results_tmp[1, , drop = FALSE]
        }

        gc()

        return(results_tmp)
      }



      res_df_multi <- res_df_multi[!duplicated(rownames(res_df_multi)), , drop = FALSE]
      final_df <- as.data.frame(res_df_multi[, !colnames(res_df_multi) %in% cluster, drop = FALSE])


      res_df <- as.matrix(res_df_multi[, colnames(res_df_multi) %in% cluster, drop = FALSE])

      if (max_combine < nrow(res_df)) {
        res_df <- res_df[1:max_combine, , drop = FALSE]
      }


      #
      colnames(res_df) <- rep(as.character(cluster), ncol(res_df))
      rm(res_df_multi)



      if (final_df$perc0[order(final_df$perc0, decreasing = FALSE)][1] == 0) {
        approved_df <- rbind(approved_df, final_df[final_df$perc0 == 0 & final_df$multi < 0.20, , drop = FALSE])
        to_exclude <- rownames(final_df)[final_df$perc0 == 0 & final_df$multi < 0.20]
        res_df <- res_df[!rownames(res_df) %in% to_exclude, , drop = FALSE]
      } else if (final_df$perc0[order(final_df$perc0, decreasing = FALSE)][1] < as.numeric(loss_val)) {
        approved_df <- rbind(approved_df, final_df[final_df$perc0 <= as.numeric(loss_val), , drop = FALSE])
      }




      if (nrow(approved_df) >= max_combine) {
        approved_df$het <- (1 - ((approved_df$perc1 + ((1 - (approved_df$perc1 + approved_df$perc0)) * 1.25)) / (str_count(string = rownames(approved_df), pattern = " ") + 1)))
        approved_df$het_adj <- (1 - ((approved_df$perc1 + ((1 - (approved_df$perc1 + approved_df$perc0)) * 1.25)) / (str_count(string = rownames(approved_df), pattern = " ") + 1))) - (approved_df$perc0 * 2)
        approved_df$cluster <- cluster

        complete_df <- rbind(approved_df, complete_df)

        break
      } else if (nrow(res_df) == 0) {
        approved_df <- rbind(approved_df, final_df[final_df$perc0 <= quantile(final_df$perc0, 0.25), , drop = FALSE])
        approved_df$het <- (1 - ((approved_df$perc1 + ((1 - (approved_df$perc1 + approved_df$perc0)) * 1.25)) / (str_count(string = rownames(approved_df), pattern = " ") + 1)))
        approved_df$het_adj <- (1 - ((approved_df$perc1 + ((1 - (approved_df$perc1 + approved_df$perc0)) * 1.25)) / (str_count(string = rownames(approved_df), pattern = " ") + 1))) - (approved_df$perc0 * 2)
        approved_df$cluster <- cluster

        complete_df <- rbind(approved_df, complete_df)

        break
      } else if (is.data.frame(tmp_recursive)) {
        if ((min(tmp_recursive$perc0) == min(final_df$perc0)) & (min(tmp_recursive & multi) <= min(final_df$multi))) {
          approved_df <- rbind(approved_df, tmp_recursive[tmp_recursive$perc0 <= quantile(tmp_recursive$perc0, 0.25), , drop = FALSE])
          approved_df$het <- (1 - ((approved_df$perc1 + ((1 - (approved_df$perc1 + approved_df$perc0)) * 1.25)) / (str_count(string = rownames(approved_df), pattern = " ") + 1)))
          approved_df$het_adj <- (1 - ((approved_df$perc1 + ((1 - (approved_df$perc1 + approved_df$perc0)) * 1.25)) / (str_count(string = rownames(approved_df), pattern = " ") + 1))) - (approved_df$perc0 * 2)
          approved_df$cluster <- cluster

          complete_df <- rbind(approved_df, complete_df)

          break
        } else {
          tmp_recursive <- final_df
        }
      } else {
        tmp_recursive <- final_df
      }

      close(pb)
      stopCluster(cl)
      gc()
    }
  }



  colnames(complete_df) <- c("loss_val", "perc_1", "avg_exp", "perc_multi", "hf", "adj_hf", "cluster")
  cat("\n\n The CSSG finish \n")

  options(old_warn)


  sc_project@metadata$cssg_markers <- complete_df

  return(sc_project)
}


###########################################################################################################################################################

#' Generate a Heatmap of Marker Genes
#'
#' This function creates a heatmap visualization of marker gene expression across cell clusters or subtypes.
#'
#' @param sc_project An object containing scRNA-seq data, including normalized expression matrices.
#' @param type A string specifying the grouping type ('subtypes' by default).
#' @param markers A character vector of gene markers to be displayed in the heatmap.
#' @param angle_col Numeric value specifying the rotation angle for column labels (default: 270).
#' @param fontsize_row Numeric value for row label font size (default: 7).
#' @param fontsize_col Numeric value for column label font size (default: 7).
#' @param font_labels Numeric value for axis label font size (default: 8).
#' @param clustering_method A string specifying the clustering method for hierarchical clustering (default: 'complete').
#' @param x_axis A string specifying the x-axis label (default: 'Cells').
#' @param y_axis A string specifying the y-axis label (default: 'Genes [log(CPM +1)]').
#' @param scale Logical. If TRUE, performs Min-Max scaling (0–1) across numeric columns (default: FALSE).
#'
#' @return A recorded plot object containing the generated heatmap.
#'
#' @export
marker_heatmap <- function(sc_project,
                           type = "subtypes",
                           markers = c(),
                           angle_col = 270,
                           fontsize_row = 7,
                           fontsize_col = 7,
                           font_labels = 8,
                           clustering_method = "complete",
                           x_axis = "Cells",
                           y_axis = "Genes [log(CPM +1)]",
                           scale = FALSE) {
  set.seed(123)

  data <- get_avg_data(sc_project = sc_project, type = type, data = "norm", chunk_size = 5000)

  data <- data[toupper(rownames(data)) %in% markers, , drop = FALSE]

  if (scale) {
    data <- data %>%
      mutate(across(where(is.numeric), ~ (. - min(.)) / (max(.) - min(.))))
  }

  pheat <- pheatmap::pheatmap(
    data,
    clustering_method = clustering_method,
    angle_col = angle_col,
    fontsize_row = fontsize_row,
    fontsize_col = fontsize_col
  )

  grid.arrange(
    pheat[[4]],
    bottom = textGrob(x_axis, gp = gpar(fontsize = font_labels)),
    left = textGrob(y_axis, rot = 90, gp = gpar(fontsize = font_labels))
  )

  pheat <- recordPlot()

  return(pheat)
}


#' Retrieve Marker Gene Names
#'
#' This function extracts marker gene names from a given `sc_project` object based on the specified cell subclass and subtype name.
#'
#' @param sc_project An object containing scRNA-seq data, including marker names.
#' @param type A string specifying the grouping type. Accepted values are 'subtypes', 'subclasses', 'cluster', or 'primary'. Default is 'subtypes'.
#'
#' @return A character vector of unique marker gene names that are present in the normalized expression matrix.
#'
#' @examples
#' markers <- get_names_markers(sc_project, type = "subtypes")
#' print(markers)
#'
#' @export
get_names_markers <- function(sc_project, type = "subtypes") {
  set.seed(123)
  if (!type %in% c("subtypes", "subclasses")) {
    stop("Invalid `type` provided. The `type` parameter should be either 'subtypes' or 'subclasses' or 'cluster' or 'primary'")
  }
  if (type == "subtypes") {
    names <- sc_project@names$subtypes
    markers_list <- unique(gsub(".* - ", "", names))
    markers_list <- unique(markers_list[toupper(markers_list) %in%
      toupper(rownames(sc_project@matrices$norm))])
  } else if (type == "subclasses") {
    slot <- "subclasses_markers"
    names <- sc_project@names$subclass
    markers_list <- unlist(strsplit(unique(names), " "))
    markers_list <- unique(markers_list[toupper(markers_list) %in%
      toupper(rownames(sc_project@matrices$norm))])
  } else if (type == "cluster") {
    slot <- "clusters_markers"
    names <- sc_project@names$cluster
    markers_list <- unlist(strsplit(unique(names), " "))
    markers_list <- unique(markers_list[toupper(markers_list) %in%
      toupper(rownames(sc_project@matrices$norm))])
  } else if (type == "primary") {
    slot <- "primary_markers"
    names <- sc_project@names$primary
    markers_list <- unlist(strsplit(unique(names), " "))
    markers_list <- unique(markers_list[toupper(markers_list) %in%
      toupper(rownames(sc_project@matrices$norm))])
  }

  return(markers_list)
}


#' Compute Averaged Expression Data
#'
#' This function retrieves and aggregates expression data from an `sc_project` object based on the specified cell type category.
#'
#' @param sc_project An object containing single-cell RNA-seq data.
#' @param type A string specifying the grouping type. Accepted values are 'subtypes' or other categorical groupings.
#' @param data A string specifying the data slot to use (e.g., 'norm' for normalized expression values). Default is 'norm'.
#' @param chunk_size An integer defining the chunk size for data aggregation. Default is 5000.
#'
#' @return A matrix or data frame with aggregated expression values.
#'
#'
#' @examples
#' avg_data <- get_avg_data(sc_project, type = "subtypes", data = "norm")
#'
#' @export
get_avg_data <- function(sc_project, type = "subtypes", data = "norm", chunk_size = 5000) {
  data <- get_data(sc_project = sc_project, type = type, data = data)

  if (!FALSE %in% unique(grepl("[^0-9]", colnames(data)))) {
    agg_subclasses <- aggregation_chr(data, chunk_size = chunk_size)
  } else {
    agg_subclasses <- aggregation_num(data, chunk_size = chunk_size)
  }


  return(agg_subclasses)
}



#' Retrieve Specific Expression Data from an sc_project Object
#'
#' This function retrieves expression data (either 'norm' or 'count') from a specified subtype or subclass grouping within a given `sc_project` object.
#'
#' @param sc_project An object containing single-cell RNA-seq data.
#' @param type A string specifying the grouping type. Accepted values are 'subtypes' or 'subclasses'.
#' @param data A string specifying the type of expression data to retrieve. Accepted values are 'norm' (normalized) or 'count' (raw counts).
#'
#' @return A sparse matrix of the specified expression data, with appropriate cell names assigned based on the chosen `type`.
#'
#' @examples
#' norm_data <- get_data(sc_project, type = "subtypes", data = "norm")
#' count_data <- get_data(sc_project, type = "subclasses", data = "count")
#'
#' @export
get_data <- function(sc_project, type = "subtypes", data = "norm") {
  set.seed(123)

  if (!type %in% c("subtypes", "subclasses")) {
    stop("Invalid `type` provided. The `type` parameter should be either 'subtypes' or 'subclasses'")
  }


  if (!data %in% c("count", "norm")) {
    stop("Invalid `data` provided. The `data` parameter should be either 'norm' or 'count'")
  }


  if (type == "subtypes") {
    names <- sc_project@names$subtypes
  } else if (type == "subclasses") {
    names <- sc_project@names$subclass
  }


  if (data == "norm") {
    sparse_matrix <- sc_project@matrices$norm
  } else if (data == "count") {
    sparse_matrix <- sc_project@matrices$count
  }


  colnames(sparse_matrix) <- names


  sparse_matrix <- sparse_matrix[, !toupper(colnames(sparse_matrix)) %in% grepl("BAD!", toupper(colnames(sparse_matrix))), drop = FALSE]

  return(sparse_matrix)
}



#' Subset a Single-Cell Project Based on Specific Groupings and Select List
#'
#' This function subsets a single-cell RNA-seq project object by selecting specific groups (subtypes, subclasses, clusters, or primary) based on the provided `select_list`.
#' Only the cells that belong to the specified `select_list` are retained in the matrices and names within the `sc_project`.
#'
#' @param sc_project An object containing single-cell RNA-seq data that has been preprocessed and categorized.
#' @param type A string indicating the groupings to subset by. Accepted values are 'subtypes', 'subclasses', 'cluster', or 'primary'.
#' @param select_list A character vector of names or identifiers to be selected. Only cells belonging to these groups will be kept in the subset.
#'
#' @return A modified `sc_project` object with only the cells that match the selected groups and identifiers.
#'
#' @export
#'
#' @examples
#' subsetted_project <- subset_project(sc_project, type = "subtypes", select_list = c("Subtype1", "Subtype3"))
#' subsetted_project <- subset_project(sc_project, type = "cluster", select_list = c("Cluster4", "Cluster7"))
subset_project <- function(sc_project, type = "subtypes", select_list = c()) {
  set.seed(123)

  if (!type %in% c("subtypes", "subclasses", "cluster", "primary")) {
    stop("Invalid `type` provided. The `type` parameter should be either 'subtypes' or 'subclasses' or 'cluster' or 'primary'")
  }


  if (type == "subtypes") {
    names <- sc_project@names$subtypes
  } else if (type == "subclasses") {
    names <- sc_project@names$subclass
  } else if (type == "cluster") {
    names <- sc_project@names$cluster
  } else if (type == "primary") {
    names <- sc_project@names$primary
  }



  logic_vec <- !toupper(names) %in% grepl("BAD!", toupper(names))
  logic_vec2 <- toupper(names) %in% toupper(select_list)

  logic_vec_sum <- logic_vec & logic_vec2


  for (i in names(sc_project@matrices)) {
    sc_project@matrices[[i]] <- sc_project@matrices[[i]][, logic_vec_sum]
  }


  for (i in names(sc_project@names)) {
    if (i != "repaired") {
      sc_project@names[[i]] <- sc_project@names[[i]][logic_vec_sum]
    }
  }


  return(sc_project)
}



#' Normalize Single-Cell RNA-seq Data
#'
#' This function normalizes the count matrix of a single-cell RNA-seq project by applying a log transformation to the data. The normalization is based on either 'counts' or 'genes' and scales the data by a specified factor (default is 1,000,000).
#' The log transformation is calculated using the formula: log2((count / gene counts) * factor + 1).
#'
#' @param sc_project An object containing single-cell RNA-seq data with pre-existing count matrices.
#' @param type A string indicating whether to normalize by 'counts' or 'genes'. Defaults to 'counts'.
#' @param factor A numeric value used to scale the normalized data. Default is 1,000,000 (for TPM-like normalization).
#'
#' @return A modified `sc_project` object with normalized data stored in the `@matrices$norm` slot.
#'
#' @examples
#' normalized_project <- normalize_data(sc_project, type = "counts", factor = 1000000)
#' normalized_project <- normalize_data(sc_project, type = "genes", factor = 10000)
#'
#' @export
normalize_data <- function(sc_project, type = "counts", factor = 1000000) {
  set.seed(123)

  data <- sc_project@matrices$count

  if (type %in% c("counts", "genes")) {
    counts <- count_genes(df = data, count = type)

    for  (r in row.names(counts)) {
      data[, r] <- log2(((data[, r] / as.integer(counts[r, ]$n)) * factor) + 1)
    }
  } else {
    stop("Invalid `type` provided. The `type` parameter should be either 'counts' or 'genes'")
  }


  data[is.na(data)] <- 0


  sc_project@matrices$norm <- data

  return(sc_project)
}



setClass(
  "scRNAproject",
  slots = list(
    metadata = "list",
    matrices = "list",
    names = "list"
  ),
  prototype = list(
    metadata = list(),
    matrices = list(),
    names = list()
  )
)



#' Create a scRNA Project from a Seurat Object
#'
#' This function converts a Seurat object into a custom `scRNAproject` object. It extracts the normalized expression data from the Seurat object, assigns it to the `norm` slot in the matrices of the custom project, and sets the active identities as the names of the cells. This facilitates further analysis using custom functions designed for `scRNAproject` objects.
#'
#' @param seurat_project A Seurat object containing single-cell RNA-seq data, including normalized gene expression data.
#'
#' @return A new `scRNAproject` object containing the matrices and names extracted from the Seurat object. The `scRNAproject` object contains:
#'   - `matrices$norm`: The normalized expression matrix from the Seurat object.
#'   - `names$primary`: The cell identities from the Seurat object.
#'
#' @examples
#' sc_project <- create_project_from_seurat(seurat_project)
#'
#' @export
create_project_from_seurat <- function(seurat_project) {
  tmp <- try(UMI@assays$RNA$data, silent = TRUE)
  if (inherits(tmp, "try-error")) {
    tmp <- try(UMI@assays$RNA@data, silent = TRUE)
  }

  colnames(tmp) <- UMI@active.ident
  matrices <- list()
  matrices[["norm"]] <- tmp
  names_list <- list()
  names_list[["primary"]] <- colnames(matrices[["norm"]])
  sc_class <- new("scRNAproject",
    metadata = list(), matrices = matrices,
    names = names_list
  )
  return(sc_class)
}


#' Create a scRNA Project from Sparse Matrix Files
#'
#' This function initializes a custom `scRNAproject` object using a sparse matrix format for single-cell RNA-seq data. The function loads the expression data from the provided matrix and row/column name files. It also assigns the appropriate names to the cells and genes, based on the provided file paths.
#'
#' @param sparse_matrix_path A character string specifying the directory containing the sparse matrix files.
#' @param sparse_name The base name of the sparse matrix file (without the file extension, e.g., '.mtx').
#' @param rows_name The base name of the file containing the gene names (without the file extension, e.g., '.tsv').
#' @param cols_name The base name of the file containing the cell names (without the file extension, e.g., '.tsv').
#' @param type A character string specifying the type of data to read. This should be either 'count' or 'norm', depending on the data format.
#'
#' @return A new `scRNAproject` object containing the loaded sparse matrix and associated cell/genes names.
#'   - `matrices$count` or `matrices$norm`: The sparse matrix loaded from the provided path.
#'   - `names$primary`: The cell names extracted from the provided column names file.
#'
#' @details
#' The function reads the sparse matrix and its associated row and column name files (typically in `.tsv` format). It creates a custom `scRNAproject` object and assigns the data to appropriate slots for easy access and analysis.
#'
#'
#' @examples
#' sc_project <- create_project("/path/to/data", "matrix_data", "genes", "cells")
#'
#' @export
create_project <- function(sparse_matrix_path, sparse_name, rows_name, cols_name, type = "count") {
  if (!type %in% c("count", "norm")) {
    stop("The parameter `type` must be either `count` or `norm`, depending on the type of data read.")
  }



  matrices <- list()
  matrices[[type]] <- readMM(file.path(sparse_matrix_path, paste0(sparse_name, ".mtx")))

  rownames(matrices[[type]]) <- readLines(file.path(sparse_matrix_path, paste0(rows_name, ".tsv")))
  colnames(matrices[[type]]) <- readLines(file.path(sparse_matrix_path, paste0(cols_name, ".tsv")))

  names_list <- list()
  names_list[["primary"]] <- colnames(matrices[[type]])

  sc_class <- new("scRNAproject", metadata = list(), matrices = matrices, names = names_list)


  return(sc_class)
}





#' Determine the Optimal Number of Principal Components (PCs) for Dimensionality Reduction
#'
#' This function estimates the optimal number of principal components (PCs) to retain in a dimensionality reduction analysis, typically used for techniques like PCA. It evaluates the standard deviation (or eigenvalue) for each principal component and identifies the point at which the explained variance starts to drop off sharply. The method applies thresholds to the decrease in standard deviation to determine the "elbow" point in the scree plot, which is often used to decide how many PCs should be retained.
#'
#' @param dim_stats A list or data frame containing the standard deviation (or eigenvalue) for each principal component. Typically, this will be from the PCA analysis output, with the relevant data stored in `dims$Elbow$data$stdev`.
#'
#' @return A numeric value indicating the optimal number of principal components (PCs) to retain based on the elbow method.
#'
#' @details
#' The function iterates through the standard deviations of each principal component and applies thresholds to detect when the variance between consecutive components begins to decrease sharply. This is often interpreted as the "elbow" in a scree plot, which represents the point where additional components contribute much less variance to the model.
#'
#'
#' @examples
#' # Assuming `dims$Elbow$data$stdev` contains standard deviation values for the PCs
#' optimal_pcs <- dim_reuction_pcs(dims$Elbow$data$stdev)
#'
#' @export
dim_reuction_pcs <- function(dim_stats) {
  set.seed(123)

  dim <- 1
  score <- c()
  element <- 0
  for (i in dims[, 1]) {
    element <- element + 1
    if (i - i * 0.01 > dims[, 1][element + 1] & element < 50 | i - i * 0.02 > dims[, 1][element + 2] & element < 49 | i - i * 0.02 > dims[, 1][element + 3] & element < 48 | i - i * 0.02 > dims[, 1][element + 4] & element < 47) {
      dim <- dim + 1
    } else {
      break
    }
  }
  dim <- as.numeric(dim)

  return(dim)
}



#' Assign Subclass and Class Names to Cells in Single-Cell RNA-Seq Data
#'
#' This function assigns meaningful names to subclasses and classes based on clustering and heterogeneity markers in a single-cell RNA sequencing project. It uses information about heterogeneity, class markers, and subclass markers to perform the annotation. It updates the subclass names in the `scRNAproject` object by combining cluster names with subclass names, based on the provided markers and species. This function assumes that clustering information and heterogeneity markers are already present in the `scRNAproject` metadata.
#' The structure of manually curated marker data, including marker class and marker subclass, is available on the GitHub page: https://github.com/jkubis96/CSSG. If set to NULL by default, the names will be based on specific markers for the clusters.
#'
#' @param sc_project An object of class `scRNAproject`, which contains the data and metadata required for subclass and class naming. It should include matrices with normalized gene expression data and metadata related to heterogeneity markers.
#' @param markers_class A list of class markers (default: NULL). If NULL the cluster will obtain native names from 'sc_project@names$primary' or 'sc_project@names$cluster'
#' @param markers_subclass A list of subclass markers (default: NULL). If NULL the markers will used from 'sc_project@metadata$heterogeneity_markers'
#' @param species A character string specifying the species (default: 'Homo sapiens'), which determines the structure of the names.
#' @param chunk_size The number of cells to process in each chunk during aggregation (default is 5000).
#'
#' @return An updated `scRNAproject` object with new subclass names assigned. The subclass names will be stored in the `names$subclass` slot of the `sc_project`.
#'
#' @details
#' The function first checks if heterogeneity markers are present in the metadata. If markers are missing or invalid, it raises an error. The subclass names are then assigned by:
#' - Using the class markers to determine class names for the clusters.
#' - Using the subclass markers and heterogeneity markers to assign meaningful subclass names to cells.
#' - Combining the class and subclass names to form the final names for each cell.
#'
#' The `aggregation_chr` and `aggregation_num` functions are used to handle different types of data (character vs numeric) during subclass aggregation.
#'
#' @examples
#' # Assuming you have a valid scRNAproject object (`sc_project`), class markers, and subclass markers
#' updated_project <- subclass_naming(sc_project, class_markers = NULL, subclass_markers = NULL, species = "Homo sapiens", chunk_size = 5000)
#'
#' @export
subclass_naming <- function(sc_project, class_markers = NULL, subclass_markers = NULL, species = "Homo sapiens", chunk_size = 5000) {
  cluster_heterogeneity_markers <- sc_project@metadata$naming_markers

  if (is.null(cluster_heterogeneity_markers) ||
    (is.atomic(cluster_heterogeneity_markers) && all(is.na(cluster_heterogeneity_markers))) ||
    (is.list(cluster_heterogeneity_markers) && length(cluster_heterogeneity_markers) == 0)) {
    stop("The value for 'cluster_heterogeneity_markers' is required and cannot be NA, NaN, NULL, or an empty list. Please use the 'heterogeneity_select_specificity' or 'heterogeneity_select_variance' function, and provide the '$heterogeneity_markers_df' value returned by these functions.")
  }

  #############################################################################
  # Class & Subclass naming


  sparse_matrix <- sc_project@matrices$norm


  for (n in names(sc_project@names)) {
    tmp_names <- sc_project@names[[n]]
    if (unique((unique(cluster_heterogeneity_markers$cluster) %in% tmp_names))) {
      new_names <- tmp_names
      break
    }
  }


  colnames(sparse_matrix) <- new_names


  if (FALSE %in% unique(grepl("[^0-9]", colnames(sparse_matrix)))) {
    agg_subclasses <- aggregation_chr(sparse_matrix, chunk_size = chunk_size)
  } else {
    agg_subclasses <- aggregation_num(sparse_matrix, chunk_size = chunk_size)
  }


  clust_names <- cluster_naming(matrix_a = agg_subclasses, markers = class_markers)

  old_clusters <- as.character(colnames(agg_subclasses))

  cell_names_2 <- subcluster_naming(agg_subclasses, subclass_markers, cluster_heterogeneity_markers, species)


  #############################################################################
  # Repair subclass_names

  new_names <- paste(clust_names, cell_names_2)
  names_to_return <- as.character(colnames(sparse_matrix))


  matching_df <- data.frame(
    old = old_clusters,
    new = new_names
  )


  sc_project@names$subclass <- matching_df$new[match(names_to_return, matching_df$old)]

  return(sc_project)
}



#' Select Marker Genes for Naming Cell Clusters
#'
#' Identifies representative marker genes for cell cluster naming by considering
#' statistical significance, gene occurrence, and removal of undesired gene
#' classes (mitochondrial/ribosomal).
#'
#' @param sc_project A single-cell project object containing marker gene metadata.
#' @param type A character string specifying the marker category to use. Must be one of:
#'   'subtypes', 'subclasses', 'cluster', or 'primary' (default).
#' @param top_n A numeric value specifying the maximum number of marker genes to retain per cluster.
#'   Default is 25.
#' @param p_val A numeric value specifying the p-value threshold for marker selection.
#'   Default is 0.05.
#' @param select_stat A character string specifying the statistic used to rank genes.
#'   Options: 'p_val' (default) or 'avg_logFC'.
#' @param mito_content A logical value indicating whether to include mitochondrial genes.
#'   Default is FALSE.
#' @param ribo_content A logical value indicating whether to include ribosomal genes.
#'   Default is FALSE.
#'
#' @return An updated `sc_project` object with selected naming marker genes stored
#'   in the `metadata$naming_markers` slot.
#'
#' @details The function performs the following steps:
#'   - Validates the `type` parameter and selects appropriate marker metadata.
#'   - Filters genes based on p-value.
#'   - Optionally removes mitochondrial and ribosomal genes.
#'   - Ranks genes by cluster based on occurrence, p-value, and log fold change.
#'   - Selects the top `top_n` genes per cluster.
#'   - Removes duplicated genes across clusters to improve naming specificity.
#'   - Ensures each cluster retains at least two marker genes.
#'
#' @examples
#' sc_project <- namign_genes_selection(sc_project,
#'   type = "primary",
#'   top_n = 25,
#'   p_val = 0.05,
#'   mito_content = FALSE,
#'   ribo_content = FALSE
#' )
#'
#' @export
namign_genes_selection <- function(sc_project, type = "primary", top_n = 25, p_val = 0.05,
                                   select_stat = "p_val", mito_content = FALSE, ribo_content = FALSE) {
  set.seed(123)

  if (!type %in% c("subtypes", "subclasses", "cluster", "primary")) {
    stop("Invalid `type` provided. The `type` parameter should be either 'subtypes' or 'subclasses' or 'cluster' or 'primary'")
  }
  if (type == "subtypes") {
    markers <- sc_project@metadata$subtypes_markers
  } else if (type == "subclasses") {
    markers <- sc_project@metadata$subclasses_markers
  } else if (type == "cluster") {
    markers <- sc_project@metadata$clusters_markers
  } else if (type == "primary") {
    markers <- sc_project@metadata$primary_markers
  }
  if (mito_content == FALSE) {
    markers <- markers %>%
      group_by(cluster) %>%
      group_modify(~ {
        mt_idx <- grepl("^(MT-|MT\\.)", toupper(.x$genes))
        cleaned <- .x[!mt_idx, , drop = FALSE]
        if (nrow(cleaned) == 0) {
          return(.x[mt_idx, , drop = FALSE])
        } else {
          return(cleaned)
        }
      }) %>%
      ungroup()
  }
  if (ribo_content == FALSE) {
    markers <- markers %>%
      group_by(cluster) %>%
      group_modify(~ {
        ribo_idx <- grepl("^(RPS|RPL|MRPL|MRPS|RS-)", toupper(.x$genes))
        cleaned <- .x[!ribo_idx, , drop = FALSE]
        if (nrow(cleaned) == 0) {
          return(.x[ribo_idx, , drop = FALSE])
        } else {
          return(cleaned)
        }
      }) %>%
      ungroup()
  }
  marker_df <- markers[markers$p_val < p_val, ]
  marker_df <- marker_df %>%
    group_by(cluster) %>%
    arrange(
      desc(pct_occurrence), desc(esm),
      desc(avg_logFC)
    ) %>%
    slice_head(n = top_n)
  dup_genes <- marker_df$genes[duplicated(marker_df$genes) |
    duplicated(marker_df$genes, fromLast = TRUE)]
  marker_clean <- marker_df %>%
    group_by(cluster) %>%
    group_modify(~ {
      cleaned <- .x %>% filter(!genes %in% dup_genes)
      if (nrow(cleaned) < 2) {
        .x %>%
          arrange(desc(pct_occurrence), p_val, desc(avg_logFC)) %>%
          slice_head(n = 2)
      } else {
        cleaned
      }
    }) %>%
    ungroup()
  sc_project@metadata$naming_markers <- marker_clean
  return(sc_project)
}
