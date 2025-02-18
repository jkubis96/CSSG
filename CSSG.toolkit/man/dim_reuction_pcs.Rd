% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/CSSG.toolkit.R
\name{dim_reuction_pcs}
\alias{dim_reuction_pcs}
\title{Determine the Optimal Number of Principal Components (PCs) for Dimensionality Reduction}
\usage{
dim_reuction_pcs(dim_stats)
}
\arguments{
\item{dim_stats}{A list or data frame containing the standard deviation (or eigenvalue) for each principal component. Typically, this will be from the PCA analysis output, with the relevant data stored in `dims$Elbow$data$stdev`.}
}
\value{
A numeric value indicating the optimal number of principal components (PCs) to retain based on the elbow method.
}
\description{
This function estimates the optimal number of principal components (PCs) to retain in a dimensionality reduction analysis, typically used for techniques like PCA. It evaluates the standard deviation (or eigenvalue) for each principal component and identifies the point at which the explained variance starts to drop off sharply. The method applies thresholds to the decrease in standard deviation to determine the "elbow" point in the scree plot, which is often used to decide how many PCs should be retained.
}
\details{
The function iterates through the standard deviations of each principal component and applies thresholds to detect when the variance between consecutive components begins to decrease sharply. This is often interpreted as the "elbow" in a scree plot, which represents the point where additional components contribute much less variance to the model.
}
\examples{
# Assuming `dims$Elbow$data$stdev` contains standard deviation values for the PCs
optimal_pcs <- dim_reuction_pcs(dims$Elbow$data$stdev)

}
