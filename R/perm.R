# Load necessary packages
library(dplyr)
library(car)
library(purrr)

# Function to compute the ANOVA F-statistic for a single row

permutation_test_all <- function(object, num_permutations = 1000) {
  # Extract necessary components from the object
  Model_selection <- object@Results$Model
  expression.data <- object@NormalizedCounts
  perm.temp <- object@.InternalNormList$fixed_data
  
  # Reorder terms in Model_selection vector
  Model_selection <- sapply(Model_selection, function(term) {
    # Split the term into components based on "+"
    components <- unlist(strsplit(term, "\\+"))
    
    # Check if "GROUP" is present
    if ("GROUP" %in% components) {
      # Remove "GROUP" from components
      components <- components[components != "GROUP"]
      # Reconstruct the term with "GROUP" at the start
      return(paste(c("GROUP", components), collapse = "+"))
    }
    return(term)  # Return original term if no "GROUP"
  })
  
  # Initialize vectors to store observed statistics and p-values
  observed_statistics <- numeric(nrow(expression.data)) 
  p_values <- numeric(nrow(expression.data))
  
  # Initialize progress bar for rows
  pb <- progress_bar$new(
    format = "  Processing [:bar] :percent ETA: :eta",
    total = nrow(expression.data), clear = FALSE, width = 60
  )
  
  # Run permutation tests for each row
  for (i in 1:nrow(expression.data)) {
    # Extract the covariates for the current row
    row_covariates <- Model_selection[i]
    perm.temp$DATA <- as.numeric(expression.data[i, ])
    
    # Create the formula dynamically using GROUP and the covariates for this row
    formula <- as.formula(paste("DATA ~", row_covariates))
    
    # Compute the observed F-statistic for this row
    aov_model <- aov(formula, data = perm.temp)
    anova_result <- car::Anova(aov_model, type = "III", singular.ok = TRUE)
    observed_statistics[i] <- anova_result$`F value`[1]  # Save observed F-statistic
    
    # Create a vector to store permuted statistics for this row
    permuted_statistics <- numeric(num_permutations)
    
    # Initialize progress bar for permutations
    pb$tick() 
    
    for (j in 1:num_permutations) {
      # Shuffle the response variable
      permuted_temp <- perm.temp %>%
        mutate(DATA = sample(DATA))
      
      # Compute the F-statistic for the permuted data
      aov_model_perm <- aov(formula, data = permuted_temp)
      anova_result_perm <- car::Anova(aov_model_perm, type = "III", singular.ok = TRUE)
      permuted_statistics[j] <- anova_result_perm$`F value`[1]  # Save permuted F-statistic
      
      
    }
    
    # Calculate the p-value
    p_values[i] <- mean(permuted_statistics >= observed_statistics[i])
    
  }
  
 
  
  # Return a data frame with observed statistics and p-values
  return(data.frame(observed_statistic = observed_statistics, p_value = p_values))
}
