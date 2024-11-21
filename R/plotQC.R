#' @name plotQC
#' @title Plot Quality Check for Principal Components and Normalized Counts
#' This function generates quality check plots for an object of class `FalconNest`.
#' It produces two plots: one for Principal Component Analysis (PCA) and one for boxplot of normalized counts.
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
#' @importFrom ggplot2 ggplot aes geom_point geom_boxplot labs theme_minimal theme ggtitle margin
#' @importFrom grid textGrob gpar
#' @importFrom gridExtra arrangeGrob grid.arrange
#' @importFrom reshape2 melt
#' @import RColorBrewer
#' @importFrom stats reorder princomp loadings
#'
#' @examples
#' colData <- data.frame(sample = factor(rep(c("A", "B"), each = 5)))
#' counts <- matrix(data = abs(rpois(50 * 10, lambda = 5)), nrow = 50, ncol = 10)
#' colnames(counts) <- rownames(colData)
#' rownames(counts) <- 1:50
#'
#' # Initialize FalconNest object
#' falcon_obj <- FalconNest(counts = counts, colData = colData, design = ~sample)
#'
#' # Normalize counts
#' normalizedObject <- normCounts(object = falcon_obj, group = "sample", mean.pick = 1, cv.pick = 100)
#' plotQC(normalizedObject, "sample")
#'
#' @export

utils::globalVariables(c("Group", "PC1", "PC2", "Color", "variable", "value", "Sample", "Expression"))

plotQC <- function(object, group) {

  set.seed(07142016)
  
  # Ensure the input is of the correct class
  if (!inherits(object, "FalconNest")) {
    stop("Input must be of class 'FalconNest'")
  }

  if (is.null(object@.InternalNormList) || length(object@.InternalNormList) == 0) {
    stop("The normalization process has not been completed. Please run the 'normCounts' method before proceeding.")
  }

  # Extract necessary information from the object
  FalconInfo <- object@colData
  normalized_counts <- object@NormalizedCounts

  if (!(group %in% colnames(FalconInfo))) {
    stop(group, " not present in colData")
  }

  # Prepare data for PCA
  group_factor <- as.factor(FalconInfo[[group]])
  group_assignment <- as.factor(FalconInfo[[group]])
  unique_groups <- unique(FalconInfo[[group]])

  
  
  # Generate colors for the groups
  if (length(unique_groups) < 3) {
    group_colors <- RColorBrewer::brewer.pal(3, "Set1")[1:length(unique_groups)]
  } else {
    group_colors <- RColorBrewer::brewer.pal(length(unique_groups), "Set1")
  }
  sample_colors <- setNames(group_colors, unique_groups)

  # Perform PCA
  pca_results <- suppressWarnings(princomp(normalized_counts, cor = TRUE))
  pca_loadings <- loadings(pca_results)

  # Calculate explained variance for the first two principal components
  total_variance <- (pca_results$sdev^2)
  variance_percentages <- round(total_variance / sum(total_variance) * 100, 1)[1:2]

  # Prepare PCA data frame
  pca_df <- data.frame(
    PC1 = pca_loadings[, 1],
    PC2 = pca_loadings[, 2],
    Color = group_factor,
    GroupAssignment = group_assignment
  )

  # Identify outliers
  pca_df <- mahalanobis_by_group(pca_df, group_var = "GroupAssignment", cols = c("PC1", "PC2"))
  outlier_coords <- pca_df[pca_df$Outlier, c("PC1", "PC2")]

  # Calculate offsets for outlier annotations
  x_offset <- diff(range(pca_df$PC1)) * 0.01
  y_offset <- diff(range(pca_df$PC2)) * 0.01

  # Create PCA plot
  p1 <- ggplot(data = pca_df, aes(x = PC1, y = PC2)) +
    geom_point(aes(fill = Color), size = 3, shape = 21) +
    labs(
      x = paste("Principal Component 1 (", variance_percentages[1], "%)", sep = ""),
      y = paste("Principal Component 2 (", variance_percentages[2], "%)", sep = ""),
      fill = group
    ) +
    theme_minimal() +
    ggtitle("PCA after noise correction") +
    theme(
      legend.position = "top",
      legend.title = element_text(hjust = 0.5),
      plot.margin = unit(c(1, 1, 1, 1), "cm"),
      plot.title = element_text(hjust = 0.5, margin = margin(b = 10)),
      legend.margin = margin(t = 10, b = 10),
      legend.box.margin = margin(t = 10)
    ) +
    scale_fill_manual(values = sample_colors) +
    # Annotate outliers
    annotate("segment",
             x = outlier_coords$PC1 - x_offset, y = outlier_coords$PC2 - y_offset,
             xend = outlier_coords$PC1, yend = outlier_coords$PC2,
             arrow = arrow(length = unit(0.3, "cm")), color = "blue", linetype = "solid", linewidth = 1.5) +
    annotate("text", x = outlier_coords$PC1, y = outlier_coords$PC2,
             label = "Outlier", hjust = -0.2, vjust = -0.5, color = "blue") +
    geom_text_repel(data = outlier_coords, aes(label = rownames(outlier_coords)),
                    nudge_y = -y_offset * 2,  # Adjust the nudge to move the label below the point
                    color = "red")

  # Prepare data for Boxplot
  normalized_counts_df <- as.data.frame(normalized_counts)
  normalized_counts_df$GeneID <- rownames(normalized_counts_df)
  boxplot_data <- melt(normalized_counts_df, id.vars = "GeneID",
                       variable.name = "Sample", value.name = "Expression")
  
  boxplot_data$group_assignment <- rep(group_assignment, each = table(boxplot_data$Sample)[1])

  # Create Boxplot
  p2 <- ggplot(boxplot_data, aes(x = Sample, y = Expression)) +
    geom_boxplot(aes(fill = group_assignment), outlier.size = 0.5, alpha = 1) +
    labs(x = group, y = "Expression (CyclicLoess(Log2(CPM+2)))", fill = group) +
    scale_fill_manual(values = sample_colors) +
    theme_minimal() +
    ggtitle("Normalized Counts After Noise Correction") +
     theme(
      legend.position = "top",
      legend.title = element_text(hjust = 0.5),
      plot.margin = unit(c(1, 1, 1, 1), "cm"),
      plot.title = element_text(hjust = 0.5, margin = margin(b = 10)),
      legend.margin = margin(t = 10, b = 10),
      legend.box.margin = margin(t = 10),
      axis.text.x = element_blank()
    )

  # Create a grob for the title
  title_grob <- textGrob("Quality Control Plots",
                         gp = gpar(fontsize = 16, fontface = "bold", col = "black"))

  # Arrange the plots side-by-side with a title
  combined_grob <- arrangeGrob(p2, p1, ncol = 2, top = title_grob)

  # Display the arranged plots with the title
  grid.arrange(combined_grob)

}
