#' Plot Results of Gene Expression Analysis
#'
#' This function generates PCA plots and heatmaps from gene expression data stored in a
#' FalconNest object. It visualizes the results based on the adjusted counts and group information.
#'
#' @param object An object of class 'FalconNest' containing the gene expression data, normalization information,
#' and results for analysis.
#' @param group A character string specifying the name of the grouping variable to be used for coloring the points
#' in the PCA plot and annotating the heatmap.
#' @param p.adjusted A numeric value indicating the adjusted p-value threshold for selecting significant genes.
#'
#' @return This function does not return a value. It generates and displays a PCA plot and a heatmap.
#' If the number of significant genes is less than the number of samples, it will display a warning message
#' and only show the heatmap.
#'
#' @examples
#' # Assuming 'falcon_nest_object' is a properly constructed FalconNest object
#' plotResults(falcon_nest_object, group = "Group", p.adjusted = 0.05)
#'
#' @importFrom RColorBrewer brewer.pal
#' @importFrom ggplot2 ggplot aes geom_point labs theme_minimal theme scale_fill_manual
#' @importFrom pheatmap pheatmap
#' @importFrom gridExtra grid.arrange
#' @importFrom grDevices colorRampPalette
#'
#' @export
plotResults <- function(object, group, p.adjusted) {

  # Ensure the input is of the correct class
  if (!inherits(object, "FalconNest")) {
    stop("Input must be of class 'FalconNest'")
  }

  if (is.null(object@.InternalNormList) || length(object@.InternalNormList) == 0) {
    stop("The normalization process has not been completed. Please run the 'normCounts' method before proceeding.")
  }

  FalconInfo <- object@colData
  GeneList <- object@Results[which(object@Results$FDR < p.adjusted), 1]
  adjusted.data <- object@AdjustedCounts[which(object@AdjustedCounts[, 1] %in% GeneList), ]

  # Check for duplicate row names
  if (any(duplicated(object@AdjustedCounts[, 1]))) {
    message("Warning: Row names of count data are not unique and will not be plotted")

    # Remove the first column and convert data
    adjusted.data <- adjusted.data[, -1, drop = FALSE]  # Remove the first column
    adjusted.data <- data.frame(lapply(adjusted.data, function(x) {
      if (is.character(x)) {
        return(as.numeric(x))  # Convert character to numeric
      }
      return(x)  # Keep numeric columns as they are
    }))

  } else {
    # Extract the first column for row names
    row.names <- adjusted.data[, 1]

    # Remove the first column and convert data
    adjusted.data <- adjusted.data[, -1, drop = FALSE]  # Remove the first column
    adjusted.data <- data.frame(lapply(adjusted.data, function(x) {
      if (is.character(x)) {
        return(as.numeric(x))  # Convert character to numeric
      }
      return(x)  # Keep numeric columns as they are
    }))

    # Set the row names
    rownames(adjusted.data) <- row.names
  }

  adjusted.data <- data.frame(lapply(adjusted.data, function(x) {
    if (is.character(x)) {
      return(as.numeric(x))  # Convert character to numeric
    } else {
      return(x)  # Keep numeric columns as they are
    }
  }))

  if (!(group %in% colnames(FalconInfo))) {
    stop(group, " not present in colData")
  }

  factor.level <- as.factor(FalconInfo[[group]])
  unique_level <- unique(FalconInfo[[group]])  # Create a unique vector of levels in the column
  numeric_mapping <- seq_along(unique_level)  # Generate a numeric sequence for the levels
  level_to_numeric <- setNames(numeric_mapping, unique_level)  # Create a named vector for mapping levels to numeric values

  # Define the total figure dimensions
  total_width <- 11  # Total width in inches
  total_height <- 6  # Total height in inches

  # Custom color palette
  group_colors <- suppressWarnings(RColorBrewer::brewer.pal(n = length(unique_level), name = "Set1"))[1:length(unique_level)]
  names(group_colors) <- unique_level

  # Dynamically calculate the number of rows in the legend
  nrow_legend <- ceiling(length(unique_level) / 3)

  if (length(GeneList) > dim(adjusted.data)[2]) {

    # PCA Plot
    pca.results <- suppressWarnings(princomp(adjusted.data, cor = TRUE))  # Perform PCA
    pca.loadings <- suppressWarnings(loadings(pca.results))  # Extract PCA loadings
    totalvar <- (pca.results$sdev^2)  # Calculate total variance
    variancePer <- round(totalvar / sum(totalvar) * 100, 1)  # Calculate percentage variance
    variancePer <- variancePer[1:3]  # Select top 3 components
    pca.df <- data.frame(PC1 = pca.loadings[, 1], PC2 = pca.loadings[, 2], Color = as.factor(factor.level))

    pca_width <- 0.4 * total_width  # 40% of total width for PCA plot

    p1 <- ggplot(data = pca.df, aes(x = PC1, y = PC2)) +
      geom_point(aes(fill = Color), size = 3, shape = 21) +
      labs(
        x = paste("Principal Component 1 (", variancePer[1], "%)", sep = ""),
        y = paste("Principal Component 2 (", variancePer[2], "%)", sep = ""),
        fill = paste(group)
      ) +
      theme_minimal() +
      theme(
        legend.position = "top",
        legend.title = element_text(hjust = 0.5),
        legend.title.position = "top",
        plot.margin = unit(c(1, 1, 1, 1), "cm"),
        plot.title = element_text(hjust = 0.5, margin = margin(b = 10)),
        legend.margin = margin(t = 10, b = 10),
        legend.box.margin = margin(t = 10)
      ) +
      scale_fill_manual(values = group_colors)

  }

  # Pheatmap plot
  # Custom color palette for the heatmap
  heatmap_colors <- colorRampPalette(c("navy", "white", "firebrick3"))(50)
  
  # Define relative dimensions for the heatmap
  heatmap_width <- 0.6 * total_width  # 60% of total width for heatmap
  cellwidth <- heatmap_width / dim(adjusted.data)[2]  # Cell width based on number of samples
  cellheight <- total_height / dim(adjusted.data)[1]  # Cell height based on number of genes

  output <- c()
  covariates <- all.vars(object@.design)
  for (index in covariates) {
    res <- is.factor(object@colData[[index]])
    output <- c(output, res)
  }
  covariates <- covariates[output]
  annotation_data <- FalconInfo[covariates]
  
  #Choosing colors to match with PCA
  annotation_colors <- list()
  annotation_colors[[group]] <- setNames(group_colors, levels(factor.level))

  # Select the reference annotation for clustering
  reference_annotation <- annotation_data[[group]]  # Change to your annotation of interest
  # Create a clustering order based on the reference annotation
  clustering_order <- order(reference_annotation)

  # Generate pheatmap with reduced resolution
  p2 <- pheatmap(
    adjusted.data,
    color = heatmap_colors,
    scale = "row",
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    fontsize_row = 8,
    fontsize_col = 8,
    fontsize = 10,
    border_color = "grey60",
    cellwidth = cellwidth * 40,
    cellheight = cellheight * 35,
    annotation_col = annotation_data,
    annotation_colors = annotation_colors,
    column_order = clustering_order,
    annotation_legend = FALSE,
    show_rownames = FALSE,
    show_colnames = FALSE,
    treeheight_row = 0,
    main = "Heatmap of Gene Expression",
    silent = TRUE  # Suppresses automatic plot printing
  )


  if (length(GeneList) > dim(adjusted.data)[2]) {
    # Combine plots using gridExtra
    grid.arrange(p1, p2[[4]], ncol = 2)
  } else {
    plot(p2[[4]])
    message("The number of selected genes is less than the number of samples, so PCA cannot be plotted.")
  }
}
