#' @title MatrixOrDataFrame Class Union
#' @description Defines a class union for objects that are either
#' matrices or data frames. This allows functions
#' and methods to handle both types of objects interchangeably.
#' @importFrom methods new
#' @export
setClassUnion("MatrixOrDataFrame", c("matrix", "data.frame"))

#' @title FalconNest Class Definition
#' @description An S4 class to represent the `FalconNest` object, which encapsulates
#' counts data, sample metadata, and results from various analyses.
#' The class also includes internal slots for design formulas and
#' normalization data.
#' @slot counts A matrix or data frame containing counts data.
#' @slot colData A data frame with sample metadata.
#' @slot NormalizedCounts A matrix or data frame of normalized counts.
#' @slot AdjustedCounts A matrix or data frame of adjusted counts.
#' @slot Results A matrix or data frame containing results from analysis.
#' @slot .design A formula representing the experimental design (internal slot).
#' @slot .GroupVariable A character string for the group variable used in analysis (internal slot).
#' @slot .InternalNormList A list containing internal normalization data (internal slot).
#' @slot .FalconHatch A list for internal use (internal slot).
#' @export
setClass(
  "FalconNest",
  slots = list(
    counts = "MatrixOrDataFrame",          # Public slot for counts data
    colData = "data.frame",                # Public slot for sample metadata
    .design = "formula",                   # Internal slot for experimental design
    NormalizedCounts = "MatrixOrDataFrame", # Public slot for normalized counts
    .GroupVariable = "character",          # Internal slot for group variable
    AdjustedCounts = "MatrixOrDataFrame",  # Public slot for adjusted counts
    Results = "MatrixOrDataFrame",         # Public slot for results
    .InternalNormList = "list",            # Internal slot for internal normalization data
    .FalconHatch = "list"                  # Internal slot for additional use
  )
)

#' Constructor for FalconNest Class
#'
#' Initializes a new `FalconNest` object with counts data, sample metadata, and an experimental design formula.
#' Internal slots are set to default values, and the object is ready for further analysis and processing.
#'
#' @param counts A matrix or data frame of counts data.
#' @param colData A data frame containing sample metadata.
#' @param design A formula representing the experimental design.
#'
#' @return An instance of the `FalconNest` class.
#' @export
#' @examples
#' # Create example data
#' colData <- data.frame(sample = factor(rep(c("A", "B"), each = 5)))
#' counts <- matrix(data = abs(rpois(50 * 10, lambda = 5)), nrow = 50, ncol = 10)
#' colnames(counts) <- rownames(colData)
#' rownames(counts) <- 1:50
#' falcon_obj <- FalconNest(counts = counts, colData = colData, design = ~sample)
FalconNest <- function(counts, colData, design) {
  # Initialize the FalconNest object
  object <- new("FalconNest",
                counts = counts,                    # Initialize counts slot
                colData = colData,                  # Initialize colData slot
                .design = design,                   # Initialize design slot (internal)
                .GroupVariable = character(),       # Initialize .GroupVariable slot (internal, empty)
                .InternalNormList = list(),         # Initialize .InternalNormList slot (internal, empty)
                NormalizedCounts = matrix(),        # Initialize NormalizedCounts slot (default empty)
                AdjustedCounts = matrix(),          # Initialize AdjustedCounts slot (default empty)
                Results = matrix(),                  # Initialize Results slot (default empty)
                .FalconHatch = list()
                )

  # Extract variable names from the design formula
  design_vars <- all.vars(design)

  # Check if all variables are present in colData
  missing_vars <- setdiff(design_vars, colnames(colData))
  if (length(missing_vars) > 0) {
    stop("The following variables in the design are not present in colData: ", paste(missing_vars, collapse = ", "))
  }

  # Identify rows with NA in the specified columns
  na_rows <- rowSums(is.na(counts)) > 0
  if (any(na_rows)) {
    # Print the row names that will be removed due to NA values
    message("Warning: Rows removed due to NA values. Count matrix has been updated accordingly:")
    cat(paste(row.names(counts)[na_rows], collapse = ", "), "\n")

    # Update the object with cleaned data
    counts <- counts[!na_rows, ]
  }

  # Check for negative values in the counts data
  if (any(counts < 0)) {
    stop("Some values in the count data are negative. Provide raw count data.")
  }

  # Check for non-integer values in the counts data
  if (any(counts %% 1 != 0)) {
    stop("Some values in the count data are non-integer. Provide raw count data.")
  }

  same_order <- identical(colnames(object@counts), rownames(object@colData))
  
  if (!same_order) {
    print("Warning: The columns of counts are NOT in the same order as the rows of colData. Reordering them, but it is highly recommended that the user checks the data order.")
    object@counts<-object@counts[,dimnames(object@colData)[[1]]]
  }
  
  return(object)
}

#' Show Method for FalconNest
#'
#' Displays a summary of the `FalconNest` object. This method provides an overview
#' of the object's main components, including the dimensions of the counts matrix,
#' the design formula, and the number of samples in the metadata.
#'
#' @param object An object of class `FalconNest`. The object should contain the following slots:
#' \describe{
#'   \item{\code{counts}}{A matrix of counts data.}
#'   \item{\code{colData}}{A data frame containing metadata about the samples.}
#'   \item{\code{design}}{A formula specifying the design of the experiment.}
#' }
#' @export
setMethod("show", "FalconNest", function(object) {
  cat("An object of class 'FalconNest'\n")
  cat("Counts matrix dimensions: ", dim(object@counts), "\n")
  cat("Design formula: ", deparse(object@.design), "\n")
  cat("Number of samples in colData: ", nrow(object@colData), "\n")
})

#' Normalization of count data
#'
#' Normalizes the counts data in the `FalconNest` object based on specified parameters.
#' The normalization is applied according to the provided group variable and other thresholds.
#'
#' @param object An instance of the `FalconNest` class.
#' @param group A character string specifying the group variable for normalization.
#' @param mean.pick An optional numeric value for selecting the mean threshold. If NULL (default),
#' mean.pick will be automatically calculated based on normalized count distribution.
#' @param cv.pick A numeric value specifying the coefficient of variation threshold (default is 20).
#'
#' @return An updated `FalconNest` object with normalized counts.
#' @export
#' @examples
#' # Create example data
#' colData <- data.frame(sample = factor(rep(c("A", "B"), each = 5)))
#' counts <- matrix(data = abs(rpois(50 * 10, lambda = 5)), nrow = 50, ncol = 10)
#' colnames(counts) <- rownames(colData)
#' rownames(counts) <- 1:50
#' falcon_obj <- FalconNest(counts = counts, colData = colData, design = ~sample)
#' # Normalize counts
#' normalizedObject <- FalconPreen(object = falcon_obj, group = "sample", mean.pick = 1, cv.pick = 100)

setGeneric("FalconPreen", function(object, group, ...) {
  standardGeneric("FalconPreen")
})

setMethod("FalconPreen", "FalconNest", function(object, group, mean.pick = NULL, cv.pick = NULL, flooring.option = "soft", flooring.ratio = 0.6, lowess.f = 1/64) {
  # Ensure the input is of the correct class
  if (!inherits(object, "FalconNest")) {
    stop("Input must be of class 'FalconNest'")
  }

  object <- normData(object, group, mean.pick, cv.pick, flooring.option, flooring.ratio, lowess.f)
  return(object)
})

#' Falcon differential analysis with covariate selection
#'
#' Applies ANCOVA analysis on the `FalconNest` object and updates the results.
#' The analysis uses the normalized counts and the specified contrasts.
#'
#' @param object An instance of the `FalconNest` class.
#' @param contrasts A character string or matrix specifying the contrasts for post-hoc testing.
#' Default is "all".
#'
#' @return The updated `FalconNest` object with results from the FalconHatch analysis.
#' @export
#' @examples
#' # Create example data
#' colData <- data.frame(sample = factor(rep(c("A", "B"), each = 5)))
#' counts <- matrix(data = abs(rpois(50 * 10, lambda = 5)), nrow = 50, ncol = 10)
#' colnames(counts) <- rownames(colData)
#' rownames(counts) <- 1:50
#' falcon_obj <- FalconNest(counts = counts, colData = colData, design = ~sample)
#' normalizedObject <- FalconPreen(object = falcon_obj, group = "sample", mean.pick = 1, cv.pick = 100)
#' # Perform analysis
#' falconResults <- FalconHatch(object = normalizedObject, random.covar = FALSE, shuffle.group = FALSE, contrasts = "all", FalconOutput= FALSE)

setGeneric("FalconHatch", function(object, ...) {
  standardGeneric("FalconHatch")
})
setMethod("FalconHatch", "FalconNest", function(object, contrasts = "all", p.adjust="BH", random.covar = TRUE, shuffle.group = TRUE, FalconOutput=TRUE) {
  # Ensure the input is of the correct class
  if (!inherits(object, "FalconNest")) {
    stop("Input must be of class 'FalconNest'")
  }

  # Ensure the NormalizedCounts slot is not NULL and has length greater than 0
  # Check if normalization has been completed
  if (is.null(object@.InternalNormList) || length(object@.InternalNormList) == 0) {
    stop("The normalization process has not been completed. Run the 'FalconPreen' method before proceeding.")
  }
  object <- runANCOVA(object, contrasts, p.adjust, random.covar, shuffle.group, FalconOutput)
  return(object)
})

#' Extract posthoc group comparisons
#'
#' Retrieves the results from the Falcon analysis based on specified p-value adjustment threshold and data type.
#' See \code{Falcon}, \code{FalconPreen}, \code{FalconNest}.
#'
#' @param object An instance of the `FalconNest` class.
#' @param p.adjusted A numeric value specifying the threshold for adjusted p-values (default is 0.05).
#' @param extendedData A logical value indicating whether to return a detailed results table,
#' including degrees of freedom, confidence interval, and unadjusted expression (default is FALSE).
#'
#' @return A data frame of results that meet the specified criteria.
#' @export
#'
#' @examples
#' # Example data
#' colData <- data.frame(sample = factor(rep(c("A", "B"), each = 5)))
#' counts <- matrix(data = abs(rpois(50 * 10, lambda = 5)), nrow = 50, ncol = 10)
#' colnames(counts) <- rownames(colData)
#' rownames(counts) <- 1:50
#'
#' # Initialize FalconNest object
#' falcon_obj <- FalconNest(counts = counts, colData = colData, design = ~sample)
#'
#' # Normalize counts
#' normalizedObject <- FalconPreen(object = falcon_obj, group = "sample", mean.pick = 1, cv.pick = 100)
#'
#' # Perform analysis
#' falconResults <- FalconHatch(object = normalizedObject, random.covar = FALSE, shuffle.group = FALSE, contrasts = "all", FalconOutput= FALSE)
#' falconResults <- getResults(object = falconResults, p.adjusted = 0.05, extendedData = FALSE)
#'
#' # Print results
#' print(falconResults)

setGeneric("getResults", function(object, p.adjusted = 0.05, extendedData = FALSE)
  standardGeneric("getResults")
)

setMethod("getResults", "FalconNest", function(object, p.adjusted = 0.05, extendedData = FALSE) {

  if (!inherits(object, "FalconNest")) {
    stop("Input must be of class 'FalconNest'")
  }

  # Check if normalization has been completed
  if (is.null(object@Results) || length(object@Results) == 0) {
    stop("Differential expression analysis has not been completed. Run the 'Falcon' method before proceeding.")
  }

  # Filter results based on p.adjusted value
  output.results <- object@Results[object@Results$FDR < p.adjusted, ]

  if (extendedData) {
    if (nrow(output.results) > 0) {
      return(output.results)  # Return results if extended data is requested
    } else {
      cat("No DEG found with p.adjusted < ", p.adjusted, "\n")
    }
  } else {
    # Define patterns for matching column names
    patterns <- c("Fstat.*_Pvalue$", "^UnadjMean_", "SumSq")

    # Function to check if a column name matches any pattern
    matches_pattern <- function(column_name, patterns) {
      any(sapply(patterns, function(pattern) grepl(pattern, column_name)))
    }

    # Find columns matching any of the patterns
    matching_columns <- names(object@Results)[sapply(names(object@Results), matches_pattern, patterns = patterns)]
    matching_indices <- which(names(output.results) %in% matching_columns)

    if (nrow(output.results) > 0) {
      return(output.results[, -matching_indices])  # Return filtered results with specific columns
    } else {
      cat("No DEG found with p.adjusted < ", p.adjusted, "\n")
    }
  }

  return(NULL)  # Return NULL if no results meet the criteria
})

#' Retrieve expression values
#'
#' Retrieves the counts data from the `FalconNest` object based on the specified type.
#' See \code{Falcon}, \code{FalconPreen}, \code{FalconNest}.
#'
#' @param object An instance of the `FalconNest` class.
#' @param type A character string specifying the type of data to retrieve. Options are "raw", "normalized", or "adjusted".
#'
#' @return The requested type of counts data.
#'
#' @export
#'
#' @examples
#' # Example data
#' colData <- data.frame(sample = factor(rep(c("A", "B"), each = 5)))
#' counts <- matrix(data = abs(rpois(50 * 10, lambda = 5)), nrow = 50, ncol = 10)
#' colnames(counts) <- rownames(colData)
#' rownames(counts) <- 1:50
#'
#' # Initialize FalconNest object
#' falcon_obj <- FalconNest(counts = counts, colData = colData, design = ~sample)
#'
#' # Normalize counts
#' normalizedObject <- FalconPreen(object = falcon_obj, group = "sample", mean.pick = 1, cv.pick = 100)
#'
#' # Perform analysis
#' falconResults <- FalconHatch(object = normalizedObject, random.covar = FALSE, shuffle.group = FALSE, contrasts = "all", FalconOutput= FALSE)
#' adjExpression <- getExpression(object = falconResults, type = "adjusted")

setGeneric("getExpression", function(object, type)
  standardGeneric("getExpression")
)

setMethod("getExpression", "FalconNest", function(object, type = "adjusted") {

  if (!inherits(object, "FalconNest")) {
    stop("Input must be of class 'FalconNest'")
  }

  if (type == "adjusted") {
    # Check if normalization has been completed
    if (is.null(object@Results) || length(object@Results) == 0) {
      stop("Differential expression analysis has not been completed. Run the 'Falcon' method before proceeding.")
    }
    return(object@AdjustedCounts)

  } else if (type == "normalized") {
    if (is.null(object@.InternalNormList) || length(object@.InternalNormList) == 0) {
      stop("The normalization process has not been completed. Run the 'FalconPreen' method before proceeding.")
    }
    return(object@NormalizedCounts)

  } else if (type == "raw") {
    return(object@counts)

  } else {
    stop("Unknown type specified. Use 'raw', 'normalized', or 'adjusted'.")
  }
})


#' Drop Covariate
#'
#' Retrieves the counts data from the `FalconNest` object based on the specified type.
#' See \code{Falcon}, \code{FalconPreen}, \code{FalconNest}.
#'
#' @param object An instance of the `FalconNest` class.
#' @param covariate A vector containing one or more covariates to be dropped.
#'
#'
#' @export
#'
#' @examples
#' # Example data
#' colData <- data.frame(sample = factor(rep(c("A", "B"), each = 5)),
#' cov1 = factor(rep(c("C", "D"), each = 5)))
#' counts <- matrix(data = abs(rpois(50 * 10, lambda = 5)), nrow = 50, ncol = 10)
#' colnames(counts) <- rownames(colData)
#' rownames(counts) <- 1:50
#'
#' # Initialize FalconNest object
#' falcon_obj <- FalconNest(counts = counts, colData = colData, design = ~sample+cov1)
#'
#' # Drop covariate
#' falcon_obj <- dropCovariate(object = falcon_obj, covariate = "cov1")
#'
#' # Perform analysis
#' falconResults <- FalconHatch(object = normalizedObject, random.covar = FALSE, shuffle.group = FALSE, contrasts = "all", FalconOutput= FALSE)
#' adjExpression <- getExpression(object = falconResults, type = "adjusted")
#' 

setGeneric("dropCovariate", function(object, covariate)
  standardGeneric("dropCovariate")
)

setMethod("dropCovariate", "FalconNest", function(object, covariate = NULL) {
  
  # Check if covariate is "NULL"
  if (is.null(covariate)) {
    message("No covariate removed.")
    return(object)
  }
  
  # Ensure covariate is provided as a character string or vector
  if (!is.character(covariate)) {
    stop("Covariate must be provided as a character string or vector.")
  }
  
  # Extract all variables from the design formula
  design.nest <- all.vars(object@.design)
  
  # Check if there's only one covariate remaining
  if (length(design.nest) == 1) {
    stop("Cannot remove covariate from a single covariate formula.")
  }
  
  # Ensure all specified covariates exist in the design formula
  if (all(covariate %in% design.nest)) {
    
    # Remove the specified covariates
    updated_covariates <- setdiff(design.nest, covariate)  
    
    # Create the updated formula manually
    updated_formula <- as.formula(paste("~", paste(updated_covariates, collapse = " + ")))
    
    # Update the design in the object
    object@.design <- updated_formula
    
    return(object)
    
  } else {
    missing_covariates <- covariate[!(covariate %in% design.nest)]
    stop(paste("The following covariates are not in the design formula:", paste(missing_covariates, collapse = ", ")))
  }
})