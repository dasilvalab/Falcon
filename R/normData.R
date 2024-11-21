#' R/normData.R
#' Normalize Counts and Perform Quality Checks
#'
#' This function normalizes count data from a `FalconNest` object using a specified normalization method. It performs various quality checks
#' and filtering based on mean and coefficient of variation (CV) thresholds. The function also updates the `FalconNest` object with the
#' normalized data and stores internal results for further analysis.
#'
#' @importFrom edgeR cpm
#' @importFrom limma normalizeBetweenArrays
#' @importFrom stats sd lowess setNames
#' @importFrom stats reorder princomp loadings
#' @param object A `FalconNest` object containing count data, metadata, and design information.
#' @param group The name of the grouping variable in the `colData` slot of the `FalconNest` object, which will be used for normalization and
#' filtering.
#' @param mean.pick An optional threshold for mean values. If `NULL`, the function will compute this threshold based on combined data from all groups.
#' @param cv.pick A numeric value representing the maximum acceptable coefficient of variation (CV) for filtering. Default is 20.
#' @param flooring.option A character string specifying the flooring method: "soft" or "hard".
#' @param flooring.ratio A numeric value between 0 and 1 used in soft flooring to determine the ratio of samples above `mean.pick`.
#' @param lowess.f A numeric value representing the span for LOWESS smoothing.
#'
#' @return The updated `FalconNest` object with normalized counts and additional internal results stored.
#'
#' @details
#' - The function first extracts count data and metadata from the `FalconNest` object and checks the design formula for required variables.
#' - It normalizes the count data using counts per million (CPM) and log transformation, followed by CyclicLoess normalization between samples.
#' - Quality checks are performed on the design formula, including checking the standard deviation of numeric variables.
#' - The function applies binary transformation to factor covariates and calculates mean, standard deviation, CV, and performs LOWESS smoothing for each group.
#' - If `mean.pick` is `NULL`, the function computes this threshold based on combined data from all groups. Otherwise, it uses the provided `mean.pick`.
#' - Data is filtered based on the mean and CV thresholds, and the results are stored in the `FalconNest` object.
#'
#' @keywords internal
#'
normData <- function(object, group, mean.pick = NULL, cv.pick = NULL, flooring.option = "soft", flooring.ratio = 0.6, lowess.f = 1/64) {
  
  if (all(!is.na(object@NormalizedCounts))) {
    
    if (any(object@.InternalNormList$Outlier == TRUE)) {
      outlier.name<-as.character(rownames(object@colData[object@.InternalNormList$Outlier,]))
      object@counts<-object@counts[,!object@.InternalNormList$Outlier]
      object@colData<-object@colData[!object@.InternalNormList$Outlier,]
      message("Warning: Removing the following outliers:")
      cat(paste(outlier.name, collapse = ", "), "\n")
    
      }else{
      stop("Normalization step already performed with no outliers found")
    }
  }
  
  # Extract counts, metadata, and design from the FalconNest object
  FalconCount <- object@counts
  FalconInfo <- object@colData
  FalconDesign <- object@.design
  object@.GroupVariable <- group  # Store the group variable into the FalconNest object
  
  # Validate count data
  if (any(FalconCount < 0)) stop("Some values in the count data are negative. Provide raw count data.")
  if (any(FalconCount %% 1 != 0)) stop("Some values in the count data are non-integer. Provide raw count data.")

  # Validate design variables
  variables <- all.vars(FalconDesign)
  if (!all(variables %in% colnames(FalconInfo))) stop("Some specified variables are not present in the data frame")

  # Ensure input class
  if (!inherits(object, "FalconNest")) stop("Input must be of class 'FalconNest'")
  if (!(group %in% variables)) stop("Variable ", group, " not found in design matrix")

  # Validate flooring ratio and option
  if (!(is.numeric(flooring.ratio) && length(flooring.ratio) == 1 && (flooring.ratio >= 0 && flooring.ratio <= 1))) {
    stop("Invalid 'flooring.ratio' selection. Choose a number between 0 and 1.")
  }
  if (!(flooring.option %in% c("soft", "hard"))) stop("No valid option selected for 'flooring.option'")

  # Check if group variable is suitable
  if (is.numeric(FalconInfo[[group]])) {
    if (all(FalconInfo[[group]] %% 1 == 0)) {
      value_counts <- table(FalconInfo[[group]])
      if (all(value_counts > 2) && length(value_counts) < 5) {
        message("Warning: It looks like you are using a numeric variable to define group levels. The data has been transformed into a factor.")
      } else {
        stop(paste(group, "variable is numeric and is not suitable for ANCOVA. Use a factor variable for group levels."))
      }
    } else {
      stop(paste(group, " is a continuous variable. ANCOVA is not suitable for continuous independent variables."))
    }
  }


  # Check for NA values in FalconInfo
  na_rows <- !complete.cases(FalconInfo[, variables])  # Identify rows with NA in the specified columns
  if (any(na_rows)) {
    message("Warning: Rows removed due to NA values. Count matrix and colData have been updated accordingly:")
    cat(paste(row.names(FalconInfo)[na_rows], collapse = ", "), "\n")

    # Remove NA rows
    FalconInfo <- FalconInfo[!na_rows, ]
    FalconCount <- FalconCount[, !na_rows]
  }

  FalconInfo <- check_categorical_factors(FalconInfo, variables)  # Check categorical variables

  if (!all(variables %in% colnames(FalconInfo))) {
    missing.data <- variables[!(variables %in% colnames(FalconInfo))]
    stop(missing.data, " not present in colData")
  }

  # Perform quality checks on the design formula
  numeric_vars <- variables[sapply(FalconInfo[variables], is.numeric)]
  if (length(numeric_vars) > 0) {
    numeric_data <- FalconInfo[, numeric_vars, drop = FALSE]
    check_sd(numeric_data)  # Perform standard deviation check
  }

  # Normalize counts data
  if (!(group %in% colnames(FalconInfo))) stop("Variable ", group, " not found in colData")
  factor.level <- as.factor(FalconInfo[[group]])
  unique_level <- unique(factor.level)  # Create a unique vector of levels in the column

  temp1 <- apply(FalconCount, 1, max) > 0  # Filter out rows with all zeros
  data.ancova <- FalconCount[temp1, ]  # Subset counts data
  within.sample.normalized.data <- cpm(data.ancova)  # Calculate counts per million (CPM)
  within.sample.normalized.data <- log(within.sample.normalized.data + 2, 2)  # Pseudo log2 normalization
  running.normalized <- normalizeBetweenArrays(within.sample.normalized.data, "cyclicloess")  # Normalize between samples using CyclicLoess

  #Calculate estimated variance per gene
  combined_means.total <- apply(running.normalized, 1, mean)  # Calculate mean over filtered data (not group-specific)
  estimated_sd.total <- apply(running.normalized, 1, sd)  # Estimate SD across filtered data
  
  # Check mean.pick handling
  if (is.null(mean.pick)) {
    mean.pick <- find_linear_point(combined_means.total, (estimated_sd.total / combined_means.total) * 100, lowess.f)
  } else if (!is.numeric(mean.pick)) {
    stop("Invalid entry for 'mean.pick'. Use either 'NULL' to auto-compute or a numeric value.")
  }

  # Initialize lists to store results
  means <- list()
  sds <- list()
  cvs <- list()
  lowess_results <- list()
  percentage_data <- list()  # To store percentage of samples above mean.pick
 

  # Loop over each unique group to calculate statistics
  for (index in unique_level) {
    group_data <- running.normalized[, factor.level == index]  # Subset data for the current group
    group_mean <- apply(group_data, 1, mean)  # Calculate mean
    group_sd <- apply(group_data, 1, sd)  # Calculate standard deviation
    group_cv <- (group_sd / group_mean) * 100  # Calculate coefficient of variation
    lowess_result <- lowess(group_cv ~ group_mean, f = lowess.f)  # Perform LOWESS smoothing

    # Store the results
    means[[as.character(index)]] <- group_mean
    sds[[as.character(index)]] <- group_sd
    cvs[[as.character(index)]] <- group_cv
    lowess_results[[as.character(index)]] <- lowess_result

    # Calculate the percentage of samples above mean.pick
    percentage_above <- apply(group_data, 1, function(row) mean(row >= mean.pick))
    percentage_data[[index]] <- percentage_above
  }

  # Combine percentage data for all groups
  percentage_df <- do.call(cbind, percentage_data)
  colnames(percentage_df) <- unique_level

  # Filter genes based on the percentage cutoff
  keepers.mean <- apply(percentage_df, 1, function(row) any(row >= flooring.ratio))
  keepers.mean <- names(keepers.mean[keepers.mean == TRUE])

  if (is.null(cv.pick)) {
    # If cv.pick is NULL, calculate it using the custom function
    cv.pick <- calculate_cv(mean.pick, lowess_results)
  } else if (is.numeric(cv.pick)) {
    # If cv.pick is numeric, use it as is
    cv.pick <- cv.pick
  } else {
    # If cv.pick is neither NULL nor numeric, stop with an error message
    stop("Select a valid 'cv.pick' option. Either 'NULL' or a number")
  }

  # Filter CV values for each group and apply CV threshold
  cv_filtered <- lapply(unique_level, function(factor.level) {
    cv <- cvs[[as.character(factor.level)]]
    cv_f <- cv[rownames(running.normalized)]  # Subset CV based on filtered data
    cv_f <= cv.pick  # Apply CV threshold
  })

  cv_x <- do.call(cbind, cv_filtered)  # Combine filtered CVs
  cv_x <- apply(cv_x, 1, sum)  # Sum CVs across groups
  keepers.cv <- names(cv_x[cv_x == length(unique_level)])  # Identify rows to keep

  # Check if no genes passed the filter parameters
  if (length(keepers.cv) == 0) {
    stop("No genes passed the filter parameters. Rerun the analysis with different filter settings.")
  }

  final.keepers <- intersect(keepers.mean,keepers.cv)
  
  # Filter based on Mean and CV
  Final.working.data <- running.normalized[final.keepers, ]
  # Flooring
  Final.working.data <- ifelse(Final.working.data < mean.pick, mean.pick, Final.working.data)

  # Apply jitter to floored values if flooring option is "soft"
  if (flooring.option == "soft") {
    # Call the add_jitter_to_mean_pick function
    Final.working.data <- add_jitter_to_mean_pick(Final.working.data, mean.pick, means, sds)

    # Re-calculate means after jittering
    means.post <- list()
    for (index in unique_level) {
      group_data.post <- Final.working.data[, factor.level == index]  # Subset data for the current group
      group_mean.post <- apply(group_data.post, 1, mean)  # Calculate mean
      means.post[[as.character(index)]] <- group_mean.post  # Store the results
    }
  } else if (flooring.option == "hard") {
    # Just return Final.working.data as is
    means.post <- list()
    for (index in unique_level) {
      group_data.post <- Final.working.data[, factor.level == index]  # Subset data for the current group
      group_mean.post <- apply(group_data.post, 1, mean)  # Calculate mean
      means.post[[as.character(index)]] <- group_mean.post  # Store the results
    }
  }

  # Perform PCA to calculate outlier
  pca_results <- suppressWarnings(princomp(Final.working.data, cor = TRUE))
  pca_loadings <- loadings(pca_results)
  group_assignment <- as.factor(FalconInfo[[group]])
  
  # Calculate explained variance for the first two principal components
  total_variance <- (pca_results$sdev^2)
  variance_percentages <- round(total_variance / sum(total_variance) * 100, 1)[1:2]
  
  # Prepare PCA data frame
  pca_df <- data.frame(
    PC1 = pca_loadings[, 1],
    PC2 = pca_loadings[, 2],
    GroupAssignment = group_assignment
  )
  
  # Identify outliers
  pca_df <- mahalanobis_by_group(pca_df, group_var = "GroupAssignment", cols = c("PC1", "PC2"))
  
  #Calculate estimated variance per gene
  combined_means <- apply(Final.working.data, 1, mean)  # Calculate mean over filtered data (not group-specific)
  estimated_sd <- apply(Final.working.data, 1, sd)  # Estimate SD across filtered data

estimated_variance <- lowess(estimated_sd ~ combined_means, f = lowess.f)  # Perform LOWESS smoothing for estimated variance
smoothed_variance <- suppressWarnings(approx(estimated_variance$x, estimated_variance$y, xout = combined_means)$y)
smoothed_variance <- smoothed_variance^2

  # Update the object with normalized counts and processed data

 InternalList <- list(
    means = means,
    sds = sds,
    cvs = cvs,
    lowess = lowess_results,
    keepers = final.keepers,
    cv.pick = cv.pick,
    mean.pick = mean.pick,
    PostNoiseSoftMeans = means.post,
    smoothed_variance= smoothed_variance,
    Outlier=pca_df$Outlier
  )

  FalconHatch<-list(
    flooring.option = flooring.option,
    flooring.ratio = flooring.ratio,
    lowess.f = lowess.f,
    cv.pick = cv.pick,
    mean.pick = mean.pick,
    FinalSampleList = colnames(FalconCount),
    Outlier = pca_df$Outlier
  )
  
  object@.InternalNormList <- InternalList
  object@.FalconHatch <- FalconHatch
  object@NormalizedCounts <- Final.working.data
  object@colData <- as.data.frame(FalconInfo)
  object@counts <- FalconCount
  

  return(object)
}
