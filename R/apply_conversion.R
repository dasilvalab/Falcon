#' Apply Binary Conversion to Categorical Covariates
#'
#' This internal function applies binary conversion to specified categorical columns in a DataFrame.
#' For factors with two levels, it converts them to binary (0/1). For factors with more than two levels,
#' it creates new binary columns for each level and removes the original column.
#'
#' @param df A data frame containing the covariates to be transformed. This data frame should consist of
#' categorical columns specified in the `cols` parameter.
#' @param cols A character vector of column names to be converted. These should be categorical columns
#' in the data frame `df`.
#' @param verbose A logical value indicating whether to print messages about the transformation process
#' to the console. Defaults to TRUE.
#'
#' @return A data frame with the specified columns converted to binary format. A message is printed to
#' the console listing all the columns that were transformed or notifying if a specified column was not found.
#'
#' @keywords internal
#'

apply_conversion <- function(df, cols, verbose=T) {
  # Ensure the input data frame is valid

  if(!(is.data.frame(df))){
    df<-as.data.frame(df)
  }

  # Iterate over the specified columns
  for (colname in cols) {
    if (colname %in% colnames(df)) {
      col <- df[[colname]]

      if (is.factor(col)) {
        levels <- levels(col)

        if (length(levels) == 2) {
          # Two-level factor: Convert to 0 and 1
          df[[colname]] <- ifelse(col == levels[1], 1, 0)

        } else if (length(levels) > 2) {
          # Multiple-level factor: Create new binary columns
          existing_cols <- colnames(df)

          for (level in levels) {
            new_colname <- paste(colname, level, sep = "_")

            # Only create new column if it does not already exist
            if (!(new_colname %in% existing_cols)) {
              df[[new_colname]] <- ifelse(col == level, 1, 0)
            }
          }

          # Remove the original column
          df[[colname]] <- NULL
        }
      }
    } else {
      if(verbose){
      cat("Column not found:", colname, "\n")
      }
    }
  }

  return(df)
}
