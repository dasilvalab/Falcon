#' Add Jitter to Mean Pick Values
#'
#' This function adds jitter to the values in the provided data matrix where the values
#' equal a specified mean pick. The jitter is generated based on the standard deviation
#' of the means corresponding to the groups. This is useful for softening values that
#' are exactly equal to the mean pick, helping to mitigate potential issues with
#' downstream analysis.
#'
#' @param Final.working.data A numeric matrix containing the normalized counts data
#' where the mean pick values will be modified.
#' @param mean.pick A numeric value indicating the mean threshold below which values
#' will be modified.
#' @param means A list of numeric vectors where each vector contains the means for each
#' group.
#' @param sds A list of numeric vectors where each vector contains the standard deviations
#' corresponding to each group.
#'
#' @return A numeric matrix with modified values where the mean pick has been jittered.
#' The output matrix maintains the same dimensions as the input matrix.
#'
#' @details
#' The function identifies the locations in the input matrix where the values are equal to
#' the mean pick, calculates the standard deviations associated with these mean values,
#' generates jitter based on these standard deviations, and replaces the mean pick values
#' with the generated jitter. If no corresponding standard deviations are found, a warning
#' message is displayed, and the function proceeds without jittering the mean pick values.
#'
#' @importFrom stats rnorm
#' @keywords internal
#'
add_jitter_to_mean_pick <- function(Final.working.data, mean.pick, means, sds, verbose=TRUE) {
  # Identify rows where any value equals mean.pick
  mean_pick_rows <- Final.working.data == mean.pick

  # Create a vector to store SD values corresponding to mean.pick
  sd_values <- numeric(0)

  # Loop through each group in means and sds
  for (group in names(means)) {
    # Round means and sds for comparison
    group_means <- round(means[[group]], 1)
    group_sds <- round(sds[[group]], 1)

    # Find indices where the group_means equal mean.pick
    match_indices <- which(group_means == round(mean.pick), 1)

    # If there are matching mean.pick values, extract the corresponding SDs
    if (length(match_indices) > 0) {
      sd_values <- c(sd_values, group_sds[match_indices])
    }
  }

  # Check if any SD values were found
  if (length(sd_values) == 0) {
    if(verbose) {
    message("Warning: No corresponding SD value for mean.pick was found. Hard flooring will be used instead.")
  }
  } else {
    # Calculate the average SD for the mean.pick
    avg_sd <- mean(sd_values)

    # Generate jitter values using rnorm with mean.pick as mean and avg_sd as sd
    n_jitter <- sum(mean_pick_rows)  # Count how many values need to be jittered
    jitter_values <- rnorm(n_jitter, mean = mean.pick, sd = avg_sd)

    # Replace the mean.pick values in Final.working.data with jitter values
    Final.working.data[mean_pick_rows] <- jitter_values

    # Calculate and display the percentage of values that were softened
    total_values <- prod(dim(Final.working.data))  # Total number of values in the data
    percentage_softened <- (n_jitter / total_values) * 100
    if(verbose) {
    message(sprintf("%.2f%% of values have been softened.", percentage_softened))
    }
  }

  # Return the modified data
  return(Final.working.data)
}
