#' AlcoholBrain Dataset
#'
#' @description The `AlcoholBrain` dataset contains two key objects: `Falcon.counts`, which holds gene expression count data, and `Falcon.colData`, which contains sample metadata. This dataset has been curated for demonstration of differential expression analysis using the Falcon pipeline.
#'
#' @usage data(AlcoholBrain)
#'
#' @format A list with two components:
#' \describe{
#'   \item{Falcon.counts}{A matrix of gene expression counts with 8000 genes across 60 samples.}
#'   \item{Falcon.colData}{A data frame with metadata for the 60 samples, containing 31 covariates (e.g., sample IDs, experimental conditions, etc.).}
#' }
#'
#' @details This dataset was initially published by XX et al. and later reanalyzed using the Falcon pipeline by MacDonald et al., 2024.
#' For performance reasons, the original dataset has been reduced using random sampling.
#'
#' @examples
#' # Load the AlcoholBrain dataset
#' data(AlcoholBrain)
#'
#' # Access the gene expression counts
#' Falcon.counts <- AlcoholBrain$Falcon.counts
#' head(Falcon.counts)
#'
#' # Access the sample metadata
#' Falcon.colData <- AlcoholBrain$Falcon.colData
#' head(Falcon.colData)
"AlcoholBrain"
