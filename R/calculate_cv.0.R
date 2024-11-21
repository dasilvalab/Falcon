#' @title Calculate Coefficient of Variation
#' This function calculates the Coefficient of Variation (CV) based on a provided mean and LOESS smoothing results.
#' It identifies the standard deviations corresponding to a specific mean from the LOESS results
#' and computes the CV as twice the mean of these standard deviations.
#'
#' @param mean.pick A numeric value representing the target mean for which to calculate the CV.
#' @param lowess A list containing LOESS results, where each entry includes 'x' (means) and 'y' (standard deviations).
#'
#' @return A numeric value representing the Coefficient of Variation (CV) calculated as twice the mean of the selected standard deviations.
#' If no valid indices are found, it returns a message indicating so.
#'
#' @keywords internal
calculate_cv <- function(mean.pick, lowess) {

  # Check if mean.pick is numeric
  if (!is.numeric(mean.pick)) {
    stop("mean.pick must be numeric")
  }

  # Initialize a list to store selected standard deviations
  lowess.pick <- numeric()

  # Loop through each group in the lowess results
  for (index in names(lowess)) {
    # Find the index of the closest mean to mean.pick
    closest_index <- which.min(abs(lowess[[index]]$x - mean.pick))
    s.pick <- as.numeric(lowess[[index]]$y[closest_index])
    lowess.pick[[index]] <- s.pick
  }

  # Calculate the CV as twice the mean of the selected standard deviations
  cv_2sd_values <- mean(unlist(lowess.pick), na.rm = TRUE) * 2

  # Return the CV value
  return(cv_2sd_values)
}
