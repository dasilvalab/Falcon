#' Check Standard Deviation and Display Warnings
#'
#' This function checks the standard deviation (SD) of numeric columns in a data frame and displays warnings
#' if any column has an SD greater than a specified threshold. High SD values can affect the fit of statistical
#' models like ANCOVA, and the function suggests scaling these variables if needed.
#'
#' @param data_frame A data frame containing the variables to be checked.
#' @param threshold A numeric value specifying the threshold for SD. Columns with SD greater than this value
#' will trigger a warning. Default is 5.
#'
#' @return The original data frame. The function returns the data frame invisibly and does not modify it.
#'
#' @details
#' The function first ensures that the input is a data frame. It then identifies numeric columns and
#' calculates their standard deviations. If any numeric column has an SD greater than the specified threshold,
#' a warning message is printed, listing these variables. If no columns exceed the threshold, the function does
#' nothing.
#'
#' @keywords internal
#'
#' @importFrom stats sd
#'
check_sd <- function(data_frame, threshold = 5) {
  # Ensure input is a data frame
  if (!is.data.frame(data_frame)) {
    stop("Input must be a data frame")
  }

  # Identify numeric columns
  numeric_cols <- sapply(data_frame, is.numeric)

  # Check if there are numeric columns to process
  if (all(!numeric_cols)) {
    return(invisible(NULL))  # Return invisibly if no numeric columns are present
  }

  # Extract numeric columns
  numeric_data <- data_frame[, numeric_cols, drop = FALSE]

  # Calculate standard deviations
  sds <- sapply(numeric_data, sd, na.rm = TRUE)

  # Identify columns with SD greater than the threshold
  problematic_cols <- names(sds)[sds > threshold]

  # Display warning if any columns exceed the threshold
  if (length(problematic_cols) > 0) {
    message("Warning: The following numeric variables have not been scaled. It is recommended to scale these variables to enhance the model's performance and accuracy.")
    cat(paste(problematic_cols, collapse = ", "), "\n")
  }

  return(data_frame)  # Return the original data frame
}
