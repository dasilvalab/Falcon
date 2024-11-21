#' Check and Convert Categorical Variables to Factors
#'
#' This function checks specified variables in a data frame to ensure they are factors.
#' If any specified categorical variables are not factors, they are converted to factors.
#'
#' @param data_frame A data frame containing the variables to be checked and converted.
#' @param variables A character vector of variable names in `data_frame` that are intended to be categorical.
#'
#' @return The original data frame with specified categorical variables converted to factors if they were not already.
#'
#' @details
#' The function identifies columns in the data frame that are categorical (non-numeric). It then checks
#' whether these categorical columns are factors. If any of the specified categorical variables are not factors,
#' they are converted to factors, and a message is displayed listing the converted variables. If all specified
#' categorical variables are already factors, a message indicating this is printed.
#'
#' @keywords internal
#'
#'
check_categorical_factors <- function(data_frame, variables) {
  # Identify non-numeric columns (potential categorical variables)
  categorical_cols <- !sapply(data_frame, is.numeric)
  categorical_data <- data_frame[, categorical_cols, drop = FALSE]

  # Identify non-factor columns among categorical columns
  non_factors <- !sapply(categorical_data, is.factor)

  # Filter variables that are categorical but not factors
  non_factorsIndesign <- non_factors[variables]
  non_factorsIndesign <- non_factorsIndesign[!is.na(non_factorsIndesign)] # Remove NA values

  # Check if there are any non-factor categorical columns to convert
  if (length(non_factorsIndesign) > 0 && any(non_factorsIndesign)) {
    non_factor_cols <- names(non_factorsIndesign)[non_factorsIndesign]
    message("Warning: The following categorical variables have been converted to factors:")
    cat(paste(non_factor_cols, collapse = ", "), "\n")

    # Convert non-factor categorical columns to factors
    data_frame[non_factor_cols] <- lapply(data_frame[non_factor_cols], as.factor)
  } 

  return(data_frame)  # Return the modified data frame
}
