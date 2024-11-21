#' @name plotCovariates
#' @title Process and Visualize Covariate Counts
#' @description This function processes model covariates, calculates their occurrence percentage,
#'              and visualizes the results using bar and line plots.
#'
#' # Initialize FalconNest object
#' falcon_obj <- FalconNest(counts = counts, colData = colData, design = ~sample)
#'
#' # Normalize counts
#' normalizedObject <- normCounts(object = falcon_obj,
#'                                 group = "sample", mean.pick = 1, cv.pick = 100)
#' falconResults <- Falcon(object = normalizedObject, contrasts = "all")
#' plotCovariates(falconResults)
#'
#' @param object A FalconNest object containing the `Results` data frame with the `Model` column.
#' @return A grid plot with two visualizations of covariate selection frequency.
#' @title Process and Visualize Covariate Counts
#' @description This function processes model covariates, calculates their occurrence percentage,
#'              and visualizes the results using bar and line plots.
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
#'                                 group = "sample", mean.pick = 1, cv.pick = 100)
#' falconResults <- Falcon(object = normalizedObject, contrasts = "all")
#' plotCovariates(falconResults)
#'
#' @param object A FalconNest object containing the `Results` data frame with the `Model` column.
#' @return A grid plot with two visualizations of covariate selection frequency.
#' @import dplyr
#' @import forcats
#' @import ggplot2
#' @import ggrepel
#' @importFrom cowplot plot_grid draw_label ggdraw
#' @importFrom tidyr unnest
#' @importFrom stringr str_trim
#' @importFrom stats reorder setNames
#' @importFrom gridExtra arrangeGrob grid.arrange
#' @importFrom tidyr unnest
#' @importFrom stringr str_trim
#' @importFrom stats reorder setNames
#' @export

utils::globalVariables(c("Model", "Models_List", "Covariates", "Percentage", "Count", "value", 
                         "freq", "set", "Root", "P.value", "Model_with_count", "bin", "percentage", "Partial_R2"))

# Process covariate counts
plotCovariates <- function(object, expanded=TRUE) {

  # Ensure the input is of the correct class
  if (!inherits(object, "FalconNest")) {
    stop("Input must be of class 'FalconNest'")
  }

  # Check if the differential expression analysis has been completed
  if (is.null(object@.InternalNormList) || length(object@.InternalNormList) == 0) {
    stop("The differential expression analysis has not been completed. Please run the 'Falcon' method before proceeding.")
  }
  
  Partial_R2<-object@.InternalNormList$Partial_R2
  
  # Process covariate counts
  covariate_counts <- object@Results %>%
    # Step 1: Split the 'Model' column by the "+" sign to separate the covariates within each model
    mutate(Models_List = strsplit(as.character(Model), "\\+")) %>%

    # Step 2: Expand the list of covariates into individual rows
    unnest(Models_List) %>%

    # Step 3: Trim any leading or trailing spaces from the covariate names
    mutate(Models_List = str_trim(Models_List)) %>%

    # Step 4: Count the occurrence of each unique covariate
    count(Models_List) %>%

    # Step 5: Calculate the percentage of each covariate's occurrence relative to the total number of models
    mutate(Percentage = n / nrow(object@Results) * 100) %>%

    # Step 6: Rename the columns for clarity
    dplyr::rename(Covariates = Models_List, Count = n) %>%

    # Replace the placeholder "GROUP" with the actual group variable name
    mutate(Covariates = sub("GROUP", object@.GroupVariable, Covariates))

  # Create a line plot to visualize the frequency of covariates
  plot1 <- ggplot(covariate_counts, aes(x = reorder(Covariates, -Count), y = Percentage, fill = Covariates)) +
    # Add bars to represent each covariate's percentage
    geom_bar(stat = "identity", position = "dodge", alpha = 0.7) +
    # Add a line connecting the points to show frequency trend
    geom_line(aes(group = 1), color = "black", linewidth = 1, linetype = "solid") +
    # Label the x-axis and y-axis
    labs(x = element_blank(), y = "Percentage selection (%)") +
    # Apply a minimal theme for simplicity and focus on the data
    theme_minimal() +
    # Rotate x-axis labels for readability
    theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none")

  # Create the Upsetplot for shared covariates
  Model <- as.vector(object@Results$Model)
  df <- data.frame(strings = Model, stringsAsFactors = FALSE)

  # Get unique elements (entities) from the data
  unique_elements <- unique(unlist(strsplit(Model, "\\+")))

  # Initialize the binary matrix
  binary_matrix <- as.data.frame(matrix(0, nrow = length(Model), ncol = length(unique_elements)))
  colnames(binary_matrix) <- unique_elements
  rownames(binary_matrix) <- 1:length(Model)

  # Populate the binary matrix
  for (i in 1:length(Model)) {
    elements <- strsplit(Model[i], "\\+")[[1]]
    binary_matrix[i, elements] <- 1
  }

  # Convert to data frame for UpSetR
  binary_matrix_df <- as.data.frame(binary_matrix)

  # Calculate the frequency of each intersection
  intersection_freq <- colSums(binary_matrix_df)
  intersection_freq_df <- data.frame(
    set = names(intersection_freq),
    freq = intersection_freq,
    stringsAsFactors = FALSE
  )

  # Get the top 5 most frequent intersections
  top_intersections <- intersection_freq_df %>%
    arrange(desc(freq)) %>%
    head(5) %>%
    pull(set)

  # Filter the binary matrix to keep only the top 5 intersections
  # Using base R approach to avoid `all_of()` issue
  binary_matrix_df_filtered <- binary_matrix_df[, top_intersections, drop = FALSE]

  p2 <- plotIntersection(binary_matrix_df_filtered)

  # Partial R2 plot
  # Combine the data frames into one
  combined_df <- do.call(rbind, lapply(Partial_R2, function(df) {
    df %>% mutate(Gene = deparse(substitute(df)))  # Add a column to identify the source
  }))

  # Remove any rows where Partial_R2 is NA or not finite
  combined_df <- combined_df %>%
    dplyr::filter(is.finite(Partial_R2))

  # Calculate the average Partial R² for each covariate
  average_partial_r2 <- combined_df %>%
    group_by(Root) %>%
    summarize(Average_Partial_R2 = mean(Partial_R2, na.rm = TRUE))

  # Order the Covariate factor levels by the average Partial R²
  combined_df$Root <- factor(combined_df$Root,
                              levels = average_partial_r2$Root[order(average_partial_r2$Average_Partial_R2)])

  # Create the violin plot with individual dots
  p3 = ggplot(combined_df, aes(x = Root, y = Partial_R2 * 100, fill = Root)) +
    geom_violin(trim = TRUE, alpha = 0.5, show.legend = FALSE) +  # No legend for violin plot
    geom_jitter(color = "black", size = 0.5, width = 0.1, alpha = 0.4, shape = 21, show.legend = FALSE) +  # No legend for jitter points
    labs(
      x = element_blank(),
      y = "% of Variance Explained") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "none")  # Ensure legend position is set to none

  # Adjust margins for individual plots
  plot1 <- plot1 + theme(plot.margin = margin(10, 10, 10, 10))
  p2 <- p2 + theme(plot.margin = margin(10, 10, 10, 10))
  p3 <- p3 + theme(plot.margin = margin(10, 10, 10, 10))

  # Combine plot1 and p3 into the top row
  top_row <- plot_grid(plot1, p3, align = 'h', axis = 'tb', ncol = 2)

  # Bottom row is just p2
  bottom_row <- p2

  # Combine the two rows into a final layout
  combined_plot <- plot_grid(
    top_row,
    bottom_row,
    ncol = 1,
    rel_heights = c(2, 1)
  ) +
    theme(plot.margin = margin(20, 20, 20, 20))  # Adjust overall margins for the combined plot

  # Add a title using draw_plot_label
  final_plot <- ggdraw(combined_plot) +
    draw_label("Covariate Scrutinization", x = 0.5, y = 1, hjust = 0.5, vjust = 1, size = 16)

  # Print the final plot with title
  print(final_plot)
  
    if(expanded) {
      
  # Initialize P-value independent plot
  df <- object@Results
  
  # Group by model and summarize P.values into a list
  df_grouped <- df %>%
    group_by(Model) %>%
    summarize(P.values = list(P.value))
  
  # Count the number of occurrences of each model
  model_counts <- df %>%
    count(Model, sort = TRUE)
  
  # Select the top 12 models based on their counts
  top_12_models <- model_counts %>%
    top_n(12, n) %>%
    pull(Model)
  
  # Filter the dataset to include only the top 12 models
  df_top_12 <- df %>%
    filter(Model %in% top_12_models)
  
  # Add model count information to the dataset
  df_top_12 <- df_top_12 %>%
    left_join(model_counts, by = "Model") %>%
    mutate(Model_with_count = paste(Model, "(n =", n, ")"))
  
  # Process the data to create bins and calculate the percentage within each bin for each model
  df_top_12_processed <- df_top_12 %>%
    group_by(Model_with_count) %>%
    mutate(bin = cut(as.numeric(P.value), breaks = seq(0, 1, by = 0.1), include.lowest = TRUE)) %>%
    group_by(Model_with_count, bin) %>%
    summarise(count = n(), .groups = "drop") %>%
    group_by(Model_with_count) %>%
    mutate(percentage = (count / sum(count)) * 100)
  
  # Add a line break before the count in the model labels for better visualization
  df_top_12_processed <- df_top_12_processed %>%
    mutate(Model_with_count = gsub(" \\(", "\n(", Model_with_count))
  
  # Plot the histogram with percentages on the y-axis and a continuous x-axis from 0 to 1
  ggplot(df_top_12_processed, aes(x = bin, y = percentage)) +
    geom_bar(stat = "identity", fill = "skyblue", color = "black", alpha = 0.7) +
    facet_wrap(~Model_with_count, scales = "free_x") +
    labs(
      title = "P-Value Distribution for Top 12 Most Frequent Models",
      x = "P-Value Bin", y = "Percentage of P-values"
    ) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    scale_x_discrete(labels = seq(0, 1, by = 0.1)) +  # Keep bins as discrete values, but label them from 0 to 1
    geom_hline(yintercept = 10, linetype = "dashed", color = "red", size = 1) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold")
    )
  }
  
}

