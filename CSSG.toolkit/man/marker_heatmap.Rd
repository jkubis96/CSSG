% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/CSSG.toolkit.R
\name{marker_heatmap}
\alias{marker_heatmap}
\title{Generate a Heatmap of Marker Genes}
\usage{
marker_heatmap(
  sc_project,
  type = "subtypes",
  markers = c(),
  angle_col = 270,
  fontsize_row = 7,
  fontsize_col = 7,
  font_labels = 8,
  clustering_method = "complete",
  x_axis = "Cells",
  y_axis = "Genes [log(CPM +1)]"
)
}
\arguments{
\item{sc_project}{An object containing scRNA-seq data, including normalized expression matrices.}

\item{type}{A string specifying the grouping type ('subtypes' by default).}

\item{markers}{A character vector of gene markers to be displayed in the heatmap.}

\item{angle_col}{Numeric value specifying the rotation angle for column labels (default: 270).}

\item{fontsize_row}{Numeric value for row label font size (default: 7).}

\item{fontsize_col}{Numeric value for column label font size (default: 7).}

\item{font_labels}{Numeric value for axis label font size (default: 8).}

\item{clustering_method}{A string specifying the clustering method for hierarchical clustering (default: 'complete').}

\item{x_axis}{A string specifying the x-axis label (default: 'Cells').}

\item{y_axis}{A string specifying the y-axis label (default: 'Genes [log(CPM +1)]').}
}
\value{
A recorded plot object containing the generated heatmap.
}
\description{
This function creates a heatmap visualization of marker gene expression across cell clusters or subtypes.
}
