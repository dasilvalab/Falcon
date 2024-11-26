% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{AlcoholBrain}
\alias{AlcoholBrain}
\title{AlcoholBrain Dataset}
\format{
A list with two components:
\describe{
\item{Falcon.counts}{A matrix of gene expression counts with 8000 genes across 60 samples.}
\item{Falcon.colData}{A data frame with metadata for the 60 samples, containing 31 covariates (e.g., sample IDs, experimental conditions, etc.).}
}
}
\usage{
data(AlcoholBrain)
}
\description{
The \code{AlcoholBrain} dataset contains two key objects: \code{Falcon.counts}, which holds gene expression count data, and \code{Falcon.colData}, which contains sample metadata. This dataset has been curated for demonstration of differential expression analysis using the Falcon pipeline.
}
\details{
This dataset was initially published by XX et al. and later reanalyzed using the Falcon pipeline by MacDonald et al., 2024.
For performance reasons, the original dataset has been reduced using random sampling.
}
\examples{
# Load the AlcoholBrain dataset
data(AlcoholBrain)

# Access the gene expression counts
Falcon.counts <- AlcoholBrain$Falcon.counts
head(Falcon.counts)

# Access the sample metadata
Falcon.colData <- AlcoholBrain$Falcon.colData
head(Falcon.colData)
}
\keyword{datasets}
