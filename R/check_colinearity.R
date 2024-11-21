#' Check for Collinearity among Covariates
#'
#' This function checks for collinearity among specified covariates using Variance Inflation Factor (VIF).
#' It identifies high collinearity terms based on a given threshold.
#'
#' @param variables A character vector of variable names to consider for collinearity checks.
#' @param ccc_report_final A character string specifying the covariate names, separated by '+'.
#' @param Final_working_data A data frame containing the working data for analysis.
#' @param i An integer index representing the row in `Final_working_data` to evaluate.
#' @param threshold A numeric value defining the threshold above which VIF indicates high collinearity (default is 10).
#'
#' @return A character string containing the names of covariates exhibiting high collinearity.
#'         If no high collinearity is detected, it returns "FALSE".
#'         If there are insufficient terms for a VIF calculation, it returns NA.
#'@import stats
#'@importFrom stringr str_split
#'
#' @keywords internal

check_colinearity <- function(variables, ccc_report_final, Final_working_data, i, threshold = 10) {

  utils::globalVariables(c("FalconInfo", "group"))
                           
  # Split the ccc_report_final string into individual covariate names
  ccc_report_final <- unlist(str_split(ccc_report_final, "\\+"))

  # Initialize the output vector for covariates
  output_covariate <- character()

  # Identify dummy variable prefixes (those in ccc_report_final with underscores)
  dummy_covariates <- ccc_report_final[grepl("_", ccc_report_final)]
  dummy_prefixes <- unique(sub("_.*$", "", dummy_covariates))

  # Add corresponding variables for dummy prefixes
  for (prefix in dummy_prefixes) {
    matching_vars <- variables[grepl(paste0("^", prefix, "$"), variables)]
    output_covariate <- unique(c(output_covariate, matching_vars))
  }

  # Add all non-dummy variables that match the names in ccc_report_final
  output_covariate <- unique(c(output_covariate, ccc_report_final[ccc_report_final %in% variables]))

  # Format the output correctly and substitute 'GROUP' if applicable
  output_covariate <- unique(sub("GROUP", group, output_covariate))

  # Prepare data for VIF calculation
  perm_VIF <- Final_working_data[i, , drop = FALSE] %>%
    as.data.frame() %>%
    setNames("DATA")  # Rename the column to 'DATA'

  perm_VIF <- cbind(perm_VIF, FalconInfo[, output_covariate, drop = FALSE])

  # Check if there is more than one term in the formula
  if (length(output_covariate) > 1) {
    # Construct the formula for VIF
    formula_str <- paste("DATA ~", paste(output_covariate, collapse = "+"))

    # Create the model and compute VIF
    VIF_model <- lm(as.formula(formula_str), data = perm_VIF)
    vif_results <- vif(VIF_model)

    # Identify collinear covariates
    high_vif_terms <- names(vif_results)[vif_results > threshold]

    if (length(high_vif_terms) > 0) {
      colinear_covariate <- paste(high_vif_terms, collapse = "+")
    } else {
      colinear_covariate <- "FALSE"
    }
  } else {
    colinear_covariate <- NA  # Not enough terms for VIF
  }

  return(colinear_covariate)
}
