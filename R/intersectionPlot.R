#' @import dplyr
#' @import dplyr
#' @importFrom patchwork plot_layout

utils::globalVariables(c("Combination", "colname", "color", "Intersection"))

plotIntersection <- function(input_data) {
  # Count unique combinations in the input data
  unique_rows <- unique(input_data)
  counts <- as.data.frame(table(apply(input_data, 1, paste, collapse = "_")))
  colnames(counts) <- c("Combination", "Count")  # Changed Pattern to Combination
  counts <- counts %>%
    arrange(desc(Count)) %>%
    mutate(Combination = fct_reorder(Combination, Count, .desc = TRUE))

  # Combine unique_rows with counts
  output_df <- cbind(unique_rows, Count = counts$Count)
  output_df$Combination <- apply(unique_rows, 1, paste, collapse = "_")  # Ensure 'Combination' column is added correctly

  # Create the binary matrix directly from the unique input data
  matrix_data <- as.matrix(output_df)

  # Create layout function
  Create_layout <- function(setup, unique_counts) {
    # Generate base layout
    Matrix_layout <- expand.grid(y = seq(ncol(setup) - 2), x = seq(nrow(setup)))
    Matrix_layout <- data.frame(Matrix_layout, value = as.vector(setup[, 1:(ncol(setup) - 2)]))

    # Initialize count and combination columns
    Matrix_layout$count <- NA
    Matrix_layout$Combination <- NA  # New column for combination reference
    Matrix_layout$colname <- rep(colnames(setup)[1:(ncol(setup) - 2)], nrow(unique_counts))

    # Assign colors and alpha values
    for (i in 1:nrow(Matrix_layout)) {
      if (!is.null(Matrix_layout$value[i]) && Matrix_layout$value[i] > 0) {
        Matrix_layout$color[i] <- "black"
        Matrix_layout$alpha[i] <- 1
        Matrix_layout$Intersection[i] <- paste(Matrix_layout$x[i], "yes", sep = "")
      } else {
        Matrix_layout$color[i] <- "black"
        Matrix_layout$alpha[i] <- 0.8
        Matrix_layout$Intersection[i] <- paste(i, "No", sep = "")
      }
    }

    # Match counts with corresponding combinations based on x values
    for (i in seq_len(nrow(unique_counts))) {
      combination <- unique_counts$Combination[i]  # Get the combination string
      count <- unique_counts$Count[i]  # Get the corresponding count

      # Assign combination and count based on x
      Matrix_layout$count[Matrix_layout$x == i] <- count
      Matrix_layout$Combination[Matrix_layout$x == i] <- combination
    }

    return(Matrix_layout)
  }

  # Generate layout
  unique_counts <- counts  # Ensure unique_counts is defined from counts
  mat_data <- Create_layout(matrix_data, unique_counts)

  # Create the main dot plot
  p <- ggplot(data = mat_data, aes(x = Combination, y = colname)) +
    geom_point(aes(fill = color, alpha=alpha), size = 6, shape = 21, stroke = 0) +  # Main points without border
    geom_line(aes(group = Intersection), color = "black", size = 1) +  # Connect points
    scale_x_discrete(expand = expansion(add = 0.5)) +  # Add space to the x-axis
    scale_y_discrete(expand = expansion(add = 0.5)) +  # Add space to the y-axis
    scale_fill_identity() +  # Keep fill colors as defined
    labs(x = "Combinations", y = "Covariates") +
    theme_minimal() +  # Minimal theme for a clean look
    theme(
      legend.position = "none",  # Remove legend
      axis.title.y = element_text(margin = margin(r = 10)),  # Adjust margin for y-axis title
      axis.text.y = element_text(margin = margin(r = 5)),  # Adjust margin for y-axis text
      axis.ticks = element_blank(),  # Remove both x and y axis ticks
      panel.grid.major = element_blank(),  # Remove major grid lines
      panel.grid.minor = element_blank()  # Remove minor grid lines
    )

  # Create a bar plot for counts aligned with the main plot
  count_plot <- ggplot(data = unique_counts, aes(x = fct_reorder(Combination, -Count), y = Count)) +
    geom_bar(stat = "identity", fill = "lightblue", alpha = 0.5) +  # Bar plot for counts
    labs(y = "Count") +  # Label for count axis
    scale_x_discrete(expand = expansion(add = 0.5)) +  # Ensure x-axis alignment
    theme_minimal() +  # Minimal theme
    scale_y_continuous(limits = c(0, max(unique_counts$Count) * 1.1)) +  # Show y=0
    theme(
      axis.title.x = element_blank(),  # Remove x-axis title
      axis.ticks.x = element_blank(),  # Remove x-axis ticks
      axis.text.x = element_blank(),    # Remove x-axis text
      panel.grid.major = element_blank(),  # Remove major grid lines
      panel.grid.minor = element_blank()   # Remove minor grid lines
    )

  # Combine the plots, ensuring the x-axes align
  combined_plot <- count_plot / p +
    plot_layout(heights = c(1, 2))  # Adjust proportions as needed

  return(combined_plot)
}
