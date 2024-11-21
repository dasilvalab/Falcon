#' Create Contrast Matrix for Falcon Package Analysis
#'
#' Generates a contrast matrix for use in differential expression analysis within the Falcon package.
#' The contrast matrix is designed to facilitate hypothesis testing between specified levels
#' of a categorical variable in a linear model fit.
#'
#' @param fit An object of class `lm` from a model fit. This object is used to extract
#' the terms from the model matrix to construct the contrast matrix.
#' @param contrasts A character vector specifying the contrasts to be created. Each contrast should be
#' defined as a comparison between two levels of a factor, separated by a colon (e.g., "level1:level2").
#' The levels should correspond to the terms in the model matrix.
#' @param reference.level A character string indicating the reference level for contrasts. Levels specified
#' as the reference level will be compared against other levels in the contrasts.
#' @param n.levels A character vector containing the levels of the categorical variable in the model.
#' This vector is used to locate the relevant columns for contrast vectors.
#'
#' @return A matrix where each row represents a specified contrast. The columns correspond to the model
#' terms, and each row contains values of 1, -1, and 0 that define the contrast vector for comparing
#' the specified levels.
#'
#' @details
#' This function is used within the Falcon package to construct contrast matrices for differential expression
#' analysis. The process includes:
#' \itemize{
#'   \item Initializing an empty list to store contrast vectors.
#'   \item Iterating through the provided contrast specifications.
#'   \item Splitting each contrast into two levels.
#'   \item Creating and populating a contrast vector with values corresponding to the levels in the model matrix.
#'   \item Combining all individual contrast vectors into a single matrix.
#'   \item Setting row names of the matrix to the contrast specifications.
#' }
#'
#' @keywords internal
#' @importFrom stats model.matrix
#'
create_contrast_matrix <- function(fit, contrasts, reference.level, n.levels) {
  # Initialize list to store contrast vectors
  contrast_list <- list()

  # Loop through each contrast specification
  for (contrast in contrasts) {
    # Split the contrast specification into two levels
    contrast_split <- unlist(strsplit(contrast, ":"))

    # Initialize contrast vector with zeros
    contrast_vector <- rep(0, length(n.levels))

    # Identify the levels for comparison
    level1 <- contrast_split[1]
    level2 <- contrast_split[2]

    # Find the positions of the levels in the model matrix
    pos_level1 <- match(level1, n.levels)
    pos_level2 <- match(level2, n.levels)

    # Check and assign values to the contrast vector
    if (!is.na(pos_level1)) {
      contrast_vector[pos_level1] <- 1  # Assign +1 to level1
    }
    if (!is.na(pos_level2)) {
      contrast_vector[pos_level2] <- -1  # Assign -1 to level2
    }

    # Store the contrast vector in the list with the contrast name as the key
    contrast_list[[contrast]] <- contrast_vector
  }

  # Combine all contrast vectors into a matrix
  contrast_matrix <- do.call(rbind, contrast_list)

  # Set row names for the contrast matrix based on contrast specifications
  rownames(contrast_matrix) <- contrasts

  return(contrast_matrix)  # Return the constructed contrast matrix
}
