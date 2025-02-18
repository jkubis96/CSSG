% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/CSSG.toolkit.R
\name{aggregation_chr}
\alias{aggregation_chr}
\title{Aggregation of Gene Expression Data (Character Clustering)}
\usage{
aggregation_chr(sparse_matrix, chunk_size = 5000)
}
\arguments{
\item{sparse_matrix}{A sparse matrix containing gene expression data.}

\item{chunk_size}{The maximum number of cells per processing batch (default: 5000).}
}
\value{
A data frame containing aggregated gene expression values.
}
\description{
Aggregates gene expression data by grouping similar clusters based on character labels.
}
\examples{
aggregated_data <- aggregation_chr(sparse_matrix)

}
