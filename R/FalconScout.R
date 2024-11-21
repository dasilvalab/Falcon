#' FalconScout: Visualize PCA and Density Plots for Informed Covariate Selection
#'
#' This function generates PCA plots for categorical variables and density plots for continuous variables
#' from a FalconNest object. It normalizes the counts and organizes plots
#' based on the specified group variable passed to FalconNest().
#'
#' @param object A FalconNest class object.
#' @param group A character string specifying the grouping variable in colData to be used for plotting.
#'
#' @return A grid of PCA and density plots organized by the variables in the design.
#' @importFrom RColorBrewer brewer.pal
#' @importFrom ggplot2 ggplot aes_string geom_density geom_point labs theme_minimal
#' @importFrom stats princomp
#' @importFrom edgeR cpm
#' @importFrom limma normalizeBetweenArrays
#' @import gridExtra
#' @export

FalconScout <- function(object, group) {
  
  plot_list <- list()  # List to store the plots
  FalconCount <- object@counts
  FalconInfo <- object@colData
  FalconDesign <- object@.design
  variables <- all.vars(FalconDesign)
  
  # Normalize counts data
  if (!(group %in% variables)) {
    stop("Variable ", group, " not found in the design formula of FalconNest")
  }
  
  # Validate count data
  if (any(FalconCount < 0)) stop("Some values in the count data are negative. Provide raw count data.")
  if (any(FalconCount %% 1 != 0)) stop("Some values in the count data are non-integer. Provide raw count data.")
  
  # Ensure input class
  if (!inherits(object, "FalconNest")) stop("Input must be of class 'FalconNest'")
  if (!(group %in% variables)) stop("Variable ", group, " not found in design matrix")
  
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
  
  
  # Define custom colors for the group levels (used in density plot)
  group_levels <- levels(FalconInfo[[group]])
  if (length(group_levels) < 3) {
    group_colors <- c("darkorange", "deepskyblue")[1:length(group_levels)]  # Manually define for 1 or 2 levels
  } else {
    group_colors <- RColorBrewer::brewer.pal(min(8, length(group_levels)), "Set1")
  }
  names(group_colors) <- group_levels

  
  factor.level <- as.factor(FalconInfo[[group]])
  unique_level <- unique(factor.level)  # Create a unique vector of levels in the column
  
  temp1 <- apply(FalconCount, 1, max) > 0  # Filter out rows with all zeros
  data.ancova <- FalconCount[temp1, ]  # Subset counts data
  within.sample.normalized.data <- cpm(data.ancova)  # Calculate counts per million (CPM)
  within.sample.normalized.data <- log(within.sample.normalized.data + 2, 2)  # Pseudo log2 normalization
  running.normalized <- normalizeBetweenArrays(within.sample.normalized.data, "cyclicloess")  # Normalize between samples using CyclicLoess
  
  # Iterate through each variable in the design vector
  for (variable in variables) {
    # Skip if the variable is the group variable
    if (variable == group) {
      next
    }
    
    # Check if the variable is factorial or continuous
    if (is.factor(FalconInfo[[variable]]) || is.character(FalconInfo[[variable]])) {
      
      # Categorical variable: Perform PCA and color by group
      pca.results <- princomp(running.normalized, cor = FALSE)  # Perform PCA on sample_info (numeric)
      pca.loadings <- pca.results$loadings
      totalvar <- (pca.results$sdev^2)
      variancePer <- round(totalvar / sum(totalvar) * 100, 1)
      
      # Ensure colData[[variable]] is a factor
      FalconInfo[[variable]] <- as.factor(FalconInfo[[variable]])
      
      # Define custom colors for the variable levels (used in PCA plot)
      var_levels <- levels(FalconInfo[[variable]])
      
      # Handle cases where there are fewer than 3 levels
      if (length(var_levels) < 3) {
        var_colors <- c("darkorange", "deepskyblue")[1:length(var_levels)]  # Manually define for 1 or 2 levels
      } else {
        var_colors <- RColorBrewer::brewer.pal(min(8, length(var_levels)), "Set1")
      }
      names(var_colors) <- var_levels
      
      # Create PCA plot
      pca_df <- data.frame(PC1 = pca.loadings[,1], PC2 = pca.loadings[,2], Group = FalconInfo[[variable]])
      p <- ggplot(pca_df, aes(x = PC1, y = PC2, color = Group)) +
        geom_point(size = 3) +
        scale_color_manual(values = var_colors) +  # Apply the color palette for colData[[variable]]
        labs(title = paste("PCA Plot for", variable), 
             x = paste0("PC1 (", variancePer[1], "%)"), 
             y = paste0("PC2 (", variancePer[2], "%)")) +
        theme_minimal() +
        theme(plot.title = element_text(hjust = 0.5, face = "bold"))
      
      plot_list[[variable]] <- p  # Add the plot to the list
      
    } else if (is.numeric(FalconInfo[[variable]])) {
      # Continuous variable: Plot density histogram of colData[[variable]]
      temp1 <- FalconInfo[[variable]]
      FalconInfo[[group]] <- as.factor(FalconInfo[[group]])  # Ensure the group variable is a factor
      
      # Calculate density for each group level
      p <- ggplot(FalconInfo, aes(x = !!sym(variable), fill = !!sym(group))) +
        geom_density(alpha = 0.5) +
        labs(title = paste("Density Plot for", variable), x = variable, y = "Density") +
        theme_minimal() +
        theme(plot.title = element_text(hjust = 0.5, face = "bold")) +
        theme(legend.title = element_blank())
      
      plot_list[[variable]] <- p  # Add the plot to the list
    } else {
      stop("Variable type not recognized.")
    }
  }
  
  # Arrange the plots in a grid
  grid.arrange(grobs = plot_list, ncol = 2)  # Adjust ncol to your preference
}

  