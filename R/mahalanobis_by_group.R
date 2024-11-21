#' @name mahalanobis_by_group
#' @title Plot Outlier Check for Principal Components and Normalized Counts
#' This function generates quality check plots for an object of class `FalconNest`.
#' It produces two plots: one for Principal Component Analysis (PCA) and one for a boxplot of normalized counts.
#' These plots help in assessing the distribution and variation in the data after noise correction.
#'
#' @param object An object of class \code{FalconNest} containing normalized counts and other relevant data.
#' @param group A character string specifying the column in \code{colData} to be used for grouping in the plots.
#'
#' @return A grid of plots is displayed, including:
#' \itemize{
#'   \item A PCA plot showing the distribution of samples along the first two principal components, colored by group.
#'   \item A boxplot of normalized counts for each variable, colored by group.
#' }
#' The function does not return a value but displays the plots directly.
#'
#' @details
#' The function first ensures that the input object is of class \code{FalconNest} and contains the required slots.
#' It then extracts the relevant data and checks if the specified group is present in the \code{colData} slot.
#' PCA is performed on the normalized counts to visualize the variation in the data, and a boxplot is created to show
#' the distribution of normalized counts across variables. Both plots use color palettes to differentiate groups and
#' are arranged side-by-side with a title indicating that they are quality check plots.
#'
#' @importFrom ggplot2 ggplot aes geom_point geom_boxplot labs theme_minimal theme ggtitle
#' @importFrom grid textGrob gpar
#' @importFrom gridExtra arrangeGrob grid.arrange
#' @importFrom reshape2 melt
#' @import RColorBrewer
#' @importFrom stats reorder setNames princomp loadings qchisq cov mahalanobis
#' @keywords internal

mahalanobis_by_group <- function(data, group_var, cols) {
  # Initialize Mahalanobis distance column
  data$Mahalanobis_Dist <- NA

  # Calculate threshold for Mahalanobis distance (chi-squared distribution)
  threshold <- qchisq(0.99, df = length(cols))

  # Loop over each unique group to compute Mahalanobis distance
  for (group in unique(data[[group_var]])) {
    # Filter data for the current group
    group_data <- data[data[[group_var]] == group, cols, drop = FALSE]

    # Compute the group mean and covariance matrix
    group_mean <- colMeans(group_data, na.rm = TRUE)
    group_cov <- cov(group_data, use = "complete.obs")  # Use complete cases only

    # Calculate Mahalanobis distance for each point in the group
    data$Mahalanobis_Dist[data[[group_var]] == group] <- mahalanobis(group_data, group_mean, group_cov)
  }

  # Identify outliers based on the threshold
  data$Outlier <- data$Mahalanobis_Dist > threshold

  return(data)
}
