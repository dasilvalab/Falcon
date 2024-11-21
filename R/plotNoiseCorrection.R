#' @name plotNoiseCorrection
#' @title Plot Pre and Post Noise Correction Data
#' This function generates plots to visually inspect the quality of noise correction in the data.
#' It creates two plots: one showing the data post-noise filtering and one pre-noise filtering,
#' with LOESS fits for the different groups. The plots help in assessing the effectiveness of the noise
#' correction and identifying any remaining issues.
#'
#' @param object An object of class \code{FalconNest} which contains slots with normalized counts and other relevant data.
#'
#' @return A grid of plots is displayed, including:
#' \itemize{
#'   \item A plot showing post noise filtering LOESS with points colored by group.
#'   \item A plot showing pre noise filtering LOESS with points and LOESS smooth lines.
#' }
#' The function does not return a value but displays the plots directly.
#'
#' @details
#' The function checks that the input object is of class \code{FalconNest} and contains the required slots.
#' It then extracts relevant data such as normalized counts, internal normalization lists, and group variables.
#' It creates two types of plots: one showing post-noise filtering data and another showing pre-noise filtering data.
#' Both plots include dashed lines for the mean and coefficient of variation thresholds. The plots are arranged
#' side-by-side with a title indicating that they are noise correction plots.
#'
#' @importFrom ggplot2 ggplot aes geom_point geom_abline geom_vline geom_hline geom_smooth labs theme_minimal theme ggtitle margin
#' @importFrom grid textGrob gpar
#' @importFrom gridExtra arrangeGrob grid.arrange
#' @import RColorBrewer
#' @importFrom stats setNames
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
#' normalizedObject <- normCounts(object = falcon_obj,
#' group = "sample", mean.pick = 1, cv.pick = 100)
#' plotNoiseCorrection(normalizedObject)
#'
#' @export

utils::globalVariables(c("Mean", "CV", "Group", "Color", "variable", "value", "Lowess.X", "Lowess.Y"))

# Plotting for quality check
plotNoiseCorrection <- function(object) {

  # Ensure the input is of the correct class
  if (!inherits(object, "FalconNest")) {
    stop("Input must be of class 'FalconNest'")
  }

  if (is.null(object@.InternalNormList) || length(object@.InternalNormList) == 0) {
    stop("The normalization process has not been completed. Please run the 'normCounts' method before proceeding.")
  }

  # Extract common elements from object
  FalconInfo <- object@colData
  Final.working.data <- object@NormalizedCounts
  group <- object@.GroupVariable
  unique_levels <- unique(FalconInfo[[group]])

    cvs <- object@.InternalNormList$cvs
    means <- object@.InternalNormList$means
    lowess<- object@.InternalNormList$lowess
    keepers <- object@.InternalNormList$keepers
    mean.pick <- object@.InternalNormList$mean.pick
    cv.pick <- object@.InternalNormList$cv.pick
    PostNoiseSoftMeans <- object@.InternalNormList$PostNoiseSoftMeans


  # Prepare data for post noise filtering plot
  p1_data <- do.call(rbind, lapply(unique_levels, function(factor.level) {
    data.frame(
      Mean = PostNoiseSoftMeans[[as.character(factor.level)]],
      CV = cvs[[as.character(factor.level)]][keepers],
      Group = as.factor(factor.level)
    )
  }))

  # Create a numeric mapping for group colors
  factor.level <- as.factor(FalconInfo[[group]])
    numeric_mapping <- seq_along(unique_levels)
    level_to_numeric <- setNames(numeric_mapping, unique_levels)
    as.numeric(level_to_numeric[factor.level])

  
    # Generate colors for the groups
    if (length(unique_levels) < 3) {
      # RColorBrewer's palettes require at least 3 colors, so duplicate or adjust for fewer groups
      group_colors <- RColorBrewer::brewer.pal(3, "Set1")[1:length(unique_levels)]
    } else {
      group_colors <- RColorBrewer::brewer.pal(length(unique_levels), "Set1")
    }

    names(group_colors) <- unique_levels

  # Plot 1: Post noise filtering LOESS
  post_filter_count <- nrow(p1_data)

  p1 <- suppressWarnings(ggplot(p1_data, aes(x = Mean, y = CV, fill = Group)) +
    geom_point(shape = 21, alpha = 0.4) +
      geom_abline(intercept = 0, slope = 0, linetype = "dashed") +
      geom_vline(xintercept = mean.pick, linetype = "dashed", linewidth = 1) +
      geom_hline(yintercept = cv.pick, linetype = "dashed", linewidth = 1)  +
    scale_fill_manual(values = group_colors) +
    labs(x = "Mean (CyclicLoess(Log2(CPM+2)))", y = "Coefficient of Variation", fill = group) +
    theme_minimal() +
    ggtitle("Post Noise Filtering") +
    theme(
      legend.position = "top",
      legend.title = element_text(hjust = 0.5),
      legend.title.position = "top",
      plot.margin = unit(c(1, 1, 1, 1), "cm"),
      plot.title = element_text(hjust = 0.5, margin = margin(b = 10)),
      legend.margin = margin(t = 10, b = 10),
      legend.box.margin = margin(t = 10)
    ) +
    annotate("text", x = Inf, y = Inf,
             label = paste("Post-filtering # genes:", post_filter_count),
             hjust = 1.1, vjust = 1.1, size = 3, color = "black", fontface = "bold"))

  # Plot 2: Pre noise filtering LOESS
  p2_data <- do.call(rbind, lapply(unique_levels, function(factor.level) {
    data.frame(
      Mean = means[[as.character(factor.level)]],
      CV = cvs[[as.character(factor.level)]],
      Lowess.X = lowess[[as.character(factor.level)]]$x,
      Lowess.Y = lowess[[as.character(factor.level)]]$y,
      Group = as.factor(factor.level)
    )
  }))

  pre_filter_count <- nrow(p2_data)

  p2 <-suppressWarnings(ggplot(p2_data, aes(x = Mean, y = CV, fill = Group)) +
    geom_point(shape = 21, alpha = 0.4) +
    #geom_smooth(aes(color = Group), method = "loess", se = FALSE, formula = y ~ x) +
      geom_line(aes(x = Lowess.X, y = Lowess.Y, color = Group), size = 1) +
    geom_abline(intercept = 0, slope = 0, linetype = "dashed") +
    geom_vline(xintercept = mean.pick, linetype = "dashed", linewidth = 1) +
    geom_hline(yintercept = cv.pick, linetype = "dashed", linewidth = 1) +
    labs(x = "Mean (CyclicLoess(Log2(CPM+2)))", y = "Coefficient of Variation", fill = group) +
    scale_fill_manual(values = group_colors) +
    scale_color_manual(values = group_colors) +
    theme_minimal() +
    ggtitle("Pre Noise Filtering") +
    theme(
      legend.position = "none",
      plot.margin = unit(c(1, 1, 1, 1), "cm"),
      plot.title = element_text(hjust = 0.5, margin = margin(b = 10))
    ) +
    annotate("text", x = Inf, y = Inf,
             label = paste("Pre-filtering # genes:", pre_filter_count),
             hjust = 1.1, vjust = 1.1, size = 3, color = "black", fontface = "bold"))

  # Create a grob for the title
  title_grob <- textGrob("Noise Correction Plots",
                         gp = gpar(fontsize = 16, fontface = "bold", col = "black"))

  # Arrange the plots and add the title
  combined_grob <- arrangeGrob(p2, p1, ncol = 2, top = title_grob)

  # Display the arranged plots with the title
  grid.arrange(combined_grob)
}
