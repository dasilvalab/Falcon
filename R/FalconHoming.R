#' FalconHoming: Validation and reproducibility of Falcon analysis.
#'
#' This function conducts all the basic steps of Falcon pipeline from data
#' normalization to differential expression analysis using FalconOutput log file.
#'
#' @param counts A numeric matrix of gene expression counts where rows are genes
#' and columns are samples.
#' @param colData A DataFrame containing sample metadata
#' @param FalconOutput A character string indicating a path to a
#' directory where the '.Falcon' is stored. If not provided,
#' the function will search for existing `.Falcon` files in the current directory.
#'
#' @importFrom stats aov step var deviance
#' @importFrom progress progress_bar
#' @importFrom car leveneTest Anova
#' @importFrom effects effect
#' @importFrom multcomp glht mcp
#' @importFrom edgeR cpm
#' @importFrom limma normalizeBetweenArrays
#' @importFrom stats sd lowess setNames rnorm model.matrix lm coef
#' @importFrom ggplot2 ggplot aes geom_point geom_boxplot labs theme_minimal theme ggtitle
#' @importFrom grid textGrob gpar
#' @importFrom reshape2 melt
#' @import RColorBrewer
#' @importFrom stats reorder setNames princomp loadings
#' @import methods
#'
#' @return A FalconHoming object containing results of the analysis.
#'
#' @examples
#' # Example usage of FalconHoming
#' result <- FalconHoming(counts = my_counts, colData = my_colData)
#'
#' @export
#'
FalconHoming <- function(counts, colData=NULL, FalconOutput=NULL) {
  
  utils::globalVariables(c("FalconInfo", "FalconDesign", "GroupOfInterest", "FalconCount", "FinalSampleList", 
                           "unique_level", "factor.level", "flooring.option", "random.covar", "mean.pick", 
                           "shuffle.group", "group"))

  # Check if FalconOutput is provided, if not, load the most recent .Falcon file
  if(is.null(FalconOutput)) {
    falcon_files <- list.files(pattern = "\\.Falcon$", full.names = TRUE)
    if(length(falcon_files) == 0) {
      stop("No '.Falcon' files were found. Please input the file using the FalconOutput option")
    }
    if(length(falcon_files) > 1) {
      message("Multiple '.Falcon' files detected. The most recent file will be used.")
      output <- dget(falcon_files[1])  # Load the first file if multiple
    } else {
      output <- dget(falcon_files)  # Load the single .Falcon file
    }
  } else {
    output <- FalconOutput  # Load user-provided FalconOutput
  }

  
  # Define a union class for Matrix or DataFrame
  setClassUnion("MatrixOrDf", c("matrix", "data.frame"))

  # Create a new class "FalconHoming" with specified slots
  setClass(
    "FalconHoming",
    slots = list(
      FalconCount = "MatrixOrDf",        # Counts data
      colData = "MatrixOrDf",         # Sample metadata
      NormalizedCounts = "MatrixOrDf",   # Normalized counts
      AdjustedCounts = "MatrixOrDf",     # Adjusted counts
      Results = "MatrixOrDf",            # Results data
      .InternalNormList = "list"                # Internal normalization list
    )
  )

  # Extract variable names from FalconOutput
  variable_names <- names(output)
  
  # Assign elements from FalconOutput to the FalconHoming object using @
  for (name in variable_names) {
    assign(name, output[[name]], envir = .GlobalEnv)  # This assigns the value of output[[name]] to a variable named `name`
  }
  
  if(is.null(colData)){
    colData<-FalconInfo
  } else{
    # Check for matching column names in colData
    colDataOriginal <- colnames(FalconInfo)
    col_names <- names(colData)
    if (!any(col_names %in% colDataOriginal)) {
      message("Warning: The provided colData does not match the colData used originally. Ignoring the inputted colData.")
      object@colData<-colData  # Fall back to original colData
    }
  }
  

  # Initialize the FalconHoming object
  object <- new("FalconHoming",
                FalconCount = counts,
                colData = colData,
                NormalizedCounts = matrix(),
                AdjustedCounts = matrix(),
                Results = matrix(),
                .InternalNormList=list())

  # Define a method to display the FalconHoming object
  setMethod("show", "FalconHoming", function(object) {
    cat("An object of class 'FalconHoming'\n")
    cat("Counts matrix dimensions: ", dim(object@FalconCount), "\n")
    cat("Number of samples in colData: ", nrow(object@colData), "\n")
  })

  # Initialize variables for analysis
  variables <- all.vars(FalconDesign)  # Extract variables from the design
  group_levels <- unique(colData[[GroupOfInterest]])  # Unique group levels
  reference.level <- levels(colData[[GroupOfInterest]])[1]  # Reference level
  group <- GroupOfInterest  # Group of interest
  partial_R2 <- list()  # List to store partial R2 values


  # Data validation checks
  if (any(FalconCount < 0)) {
    stop("Some values in the count data are negative. Provide raw count data.")
  }

  if (any(FalconCount %% 1 != 0)) {
    stop("Some values in the count data are non-integer. Provide raw count data.")
  }

  # Check if variables specified are in the data frame
  if (!all(variables %in% colnames(colData))) {
    stop("Some specified variables are not present in the data frame")
  }

  # Ensure the input is of the correct class
  if (!inherits(object, "FalconHoming")) {
    stop("Input must be of class 'FalconHoming'")
  }

  if (!(group %in% variables)) {
    stop("Variable ", group, " not found in design matrix")
  }

  # Normalize counts data
  temp1 <- apply(counts, 1, max) > 0  # Filter out rows with all zeros
  data.ancova <- counts[temp1, FinalSampleList]  # Subset counts data
  within.sample.normalized.data <- cpm(data.ancova)  # Calculate counts per million (CPM)
  within.sample.normalized.data <- log(within.sample.normalized.data + 2, 2)  # Pseudo log2 normalization
  running.normalized <- normalizeBetweenArrays(within.sample.normalized.data, "cyclicloess")  # Normalize between samples

  combined_means <- apply(running.normalized, 1, mean)  # Calculate mean across all data
  estimated_sd <- apply(running.normalized, 1, sd)  # Estimate SD across all data

  # Initialize lists to store results
  means <- list()
  sds <- list()
  cvs <- list()
  lowess_results <- list()
  percentage_data <- list()  # To store percentage of samples above mean.pick

  # Loop over each unique group to calculate mean, SD, CV, and perform LOWESS smoothing
  for (index in unique_level) {
    group_data <- running.normalized[, factor.level == index]  # Subset data for the current group
    group_mean <- apply(group_data, 1, mean)  # Calculate mean
    group_sd <- apply(group_data, 1, sd)  # Calculate standard deviation
    group_cv <- (group_sd / group_mean) * 100  # Calculate coefficient of variation
    lowess_result <- lowess(group_cv ~ group_mean, f = output$lowess.f)  # Perform LOWESS smoothing

    # Store the results
    means[[as.character(index)]] <- group_mean
    sds[[as.character(index)]] <- group_sd
    cvs[[as.character(index)]] <- group_cv
    lowess_results[[as.character(index)]] <- lowess_result

    # Calculate the percentage of samples above mean.pick
    percentage_above <- apply(group_data, 1, function(row) mean(row >= mean.pick))
    percentage_data[[index]] <- percentage_above
  }

  # Combine percentage data for all groups
  percentage_df <- do.call(cbind, percentage_data)
  colnames(percentage_df) <- unique_level

  # Filter genes based on the percentage cutoff
  keep_genes <- apply(percentage_df, 1, function(row) all(row >= flooring.ratio))

  # Apply the filtering to the working data
  Final.working.data <- running.normalized[keep_genes, ]
  # Floor values below mean.pick
  Final.working.data <- ifelse(Final.working.data < mean.pick, mean.pick, Final.working.data)

  # Filter CV values for each group and apply CV threshold
  cv_filtered <- lapply(unique_level, function(factor.level) {
    cv <- cvs[[as.character(factor.level)]]
    cv_f <- cv[rownames(Final.working.data)]  # Subset CV based on filtered data
    cv_f < cv.pick  # Apply CV threshold
  })

  cv_x <- do.call(cbind, cv_filtered)  # Combine filtered CVs
  cv_x <- apply(cv_x, 1, sum)  # Sum CVs across groups
  keepers <- names(cv_x[cv_x == length(unique_level)])  # Identify rows to keep

  # Filter Final.working.data based on CV
  Final.working.data <- Final.working.data[keepers, ]

  # Apply jitter to floored values based on flooring option
  means.post <- list()
  if(flooring.option == "soft") {
    Final.working.data <- add_jitter_to_mean_pick(Final.working.data, mean.pick, means, sds, verbose=FALSE)

    for (index in unique_level) {
      group_data.post <- Final.working.data[, factor.level == index]  # Subset data for the current group
      group_mean.post <- apply(group_data.post, 1, mean)  # Calculate mean
      means.post[[as.character(index)]] <- group_mean.post  # Store the results
    }

  } else if (flooring.option == "hard") {
    for (index in unique_level) {
      group_data.post <- Final.working.data[, factor.level == index]  # Subset data for the current group
      group_mean.post <- apply(group_data.post, 1, mean)  # Calculate mean
      means.post[[as.character(index)]] <- group_mean.post  # Store the results
    }
  }

  # Differential expression analysis
  # Apply binary transformation to factor covariates
  diff <- variables[variables != GroupOfInterest]
  colData.binary <- colData  # Initialize binary version of colData
  colData.binary <- apply_conversion(colData.binary, diff, verbose=FALSE)
  names(colData.binary)[names(colData.binary) == "binary"] <- GroupOfInterest

  names(colData.binary)[names(colData.binary) == GroupOfInterest] <- "GROUP"

  # Clean up column names
  colnames(colData.binary) <- gsub("[/ -]", "_", colnames(colData.binary))
  binary.names <- setdiff(as.vector(colnames(colData.binary)), as.vector(colnames(colData)))

  # Initialize progress bar
  pb <- progress_bar$new(
    format = "  Processing [:bar] :percent ETA: :eta",
    total = nrow(Final.working.data), clear = FALSE, width = 60
  )

  # Prepare ANCOVA formula
  components <- c(variables[variables %in% names(colData.binary)], binary.names)

  running.posthoc <- NULL
  running.adj.data <- NULL

  components[components == GroupOfInterest] <- "GROUP"
  fixed_data <- colData.binary[, components, drop = FALSE]

  if (!(is.null(random.covar))) {

    ANCOVA.formula <- paste("DATA ~", paste(c(components, "random.covar"), collapse = " + "))
    fixed_data <- fixed_data %>%
      mutate(random.covar = random.covar.list)

  } else {
    ANCOVA.formula <- paste("DATA ~", paste(components, collapse = " + "))
  }

  if (!(is.null(shuffle.group))) {

    ANCOVA.formula <- paste(c(ANCOVA.formula,"shuffled.group"), collapse = " + ")
    fixed_data <- fixed_data %>%
      mutate(shuffled.group = shuffle.group.list)

  }

  ANCOVA.formula.transformed <- as.formula(ANCOVA.formula)

  # Process each gene
  for (i in seq_len(nrow(Final.working.data))) {
    pb$tick()
    perm.temp <- fixed_data
    perm.temp$DATA <- as.numeric(t(Final.working.data[i, ]))

    if (var(perm.temp$DATA) > 0) {
      c1.init <- aov(ANCOVA.formula.transformed, data = perm.temp)
      c1.go <- try(stats::step(c1.init, trace = 0))

      if (inherits(c1.go, "try-error")) next

      c1.step <- try(stats::step(c1.init, trace = 0))
      if (inherits(c1.step, "try-error")) next

      ccc <- dimnames(summary(c1.step)[[1]])[[1]]
      ccc <- sub("\\s+", "", ccc)
      if (!"GROUP" %in% ccc) {
        ccc <- c("GROUP", ccc)
      }
      ccc.check <- ccc == "GROUP"
      if (sum(ccc.check) == 1) {
        ccc.report <- ccc

        if (deviance(c1.step) > sqrt(.Machine$double.eps)) {
          ccc.report <- ccc.report[-length(ccc.report)]
          ccc.report <- paste(ccc.report, collapse = "+")
          ccc.report.final <- gsub("\\+random.covar$", "", ccc.report)
          ccc <- paste("c1.final <- aov(DATA ~", ccc.report.final, ", data = perm.temp)")
          eval(parse(text = ccc))
          temp <- car::Anova(c1.final, type = "III", singular.ok = TRUE)
          l.test <- leveneTest(perm.temp$DATA, perm.temp$GROUP, center = median)
          selected.covariates<-unlist(str_split(ccc.report.final, "\\+"))

          # Extract sums of squares for covariates and residuals
          ss_covariates <- temp$`Sum Sq`[2:(nrow(temp))]  # Exclude the intercept and residuals

          # Calculate partial R² for each covariate
          gene_partial_r2 <- ss_covariates / sum(ss_covariates)

          # Create a data frame to display covariate names and their partial R²
          partial_r2_df <- data.frame(
            Covariate = rownames(temp)[-1],  # Exclude the intercept from the row names
            Partial_R2 = c(gene_partial_r2)
          )

          # Function to find matching root name
          partial_r2_df$Root <- sapply(partial_r2_df$Covariate, function(covariate) {
            # Find the original variable that matches the covariate's root
            match_index <- which(sapply(variables, function(root_var) {
              grepl(paste0("^", root_var), covariate)  # Check if the root_var is a prefix of covariate
            }))

            # Return the matched root variable if found
            if (length(match_index) > 0) {
              return(variables[match_index])
            } else {
              return(covariate)  # Return the covariate if no match found
            }
          })

          # Sum Partial_R2 values for covariates with the same root name
          summed_partial_r2_df <- aggregate(Partial_R2 ~ Root, data = partial_r2_df, sum)

          partial_R2[[i]]<-summed_partial_r2_df

          # VIF calculation
          if(length(selected.covariates) > 1) {

            ccc.VIF <- paste("VIF.final <- lm(DATA ~", ccc.report.final, ", data = perm.temp)")
            eval(parse(text = ccc.VIF)) # Evaluate the formula
            vif_results <- vif(VIF.final) # Compute VIF
            high_vif_terms <- names(vif_results)[vif_results > 8]  # Check for colinearity

            high_vif_terms <- sapply(high_vif_terms, function(suffix_var) {
              # Look for the original variable that is a substring of the suffixed variable
              match_index <- which(sapply(variables, function(orig_var) {
                grepl(paste0("^", orig_var), suffix_var)
              }))

              # If a match is found, return the original variable; otherwise return NA
              if (length(match_index) > 0) {
                return(variables[match_index])
              } else {
                return(suffix_var)  # Return the suffixed variable itself
              }
            })

            if(length(high_vif_terms) > 0) {
              colinear_covariate <- paste(high_vif_terms, collapse = "+")
            } else {
              colinear_covariate <- "FALSE"
            }

          } else {
            colinear_covariate <- NA
          }

          vif_test<-colinear_covariate

          # Collect results
          out2 <- c(ccc.report, l.test$'Pr(>F)'[1], vif_test, temp[[1]][2], temp[[2]][2], temp[[3]][2], temp[[4]][2])
          names(out2) <- c("Model", "Levene_P.value", "Colinearity", "SumSq", "DF", "Fstat", "P.value")

          # If the model has more coefficients than levels of GROUP, adjust the coefficients
          if (length(c1.final[[1]]) > length(levels(perm.temp$GROUP))) {
            # Get the names of the coefficients
            coef_names <- names(c1.final[[1]])  # Identify the coefficients related to the GROUP factor
            group_factors <- grep("^GROUP", coef_names, value = TRUE)
            adj.coeficients <- c1.final[[1]][!coef_names %in% c("(Intercept)", group_factors)] * -1 # Exclude the intercept (first coefficient) and the GROUP factor coefficients
            temp2 <- perm.temp[, names(adj.coeficients)]
            temp3 <- as.numeric(adj.coeficients)

            if (length(temp3) > 1) {
              for (j in seq_along(temp3)) {
                temp2[, j] <- (as.numeric(temp2[, j]) - 1) * temp3[j]
              }
              dimnames(temp2)[[2]] <- names(adj.coeficients)
            } else {
              temp2 <- as.numeric(temp2) * temp3
              temp2 <- matrix(temp2, , 1)
              dimnames(temp2)[[2]] <- names(adj.coeficients)
            }

            temp3 <- cbind(as.numeric(perm.temp$DATA), temp2)
            temp4 <- apply(temp3, 1, sum)
            temp5 <- temp4 - (mean(temp4) - mean(perm.temp$DATA))

            running.adj.data <- rbind(running.adj.data, c(dimnames(Final.working.data)[[1]][i], temp5))

          } else {
            temp5 <- Final.working.data[i, ]
            running.adj.data <- rbind(running.adj.data, c(dimnames(Final.working.data)[[1]][i], temp5))
          }

          # Compute means
          out1 <- sapply(group_levels, function(group.index) {
            mean(perm.temp$DATA[colData[[GroupOfInterest]] == group.index], na.rm = TRUE)
          })
          names(out1) <- paste("UnadjMean", group_levels, sep = "_")

          m.adj <- effect("GROUP", c1.final)
          out3 <- as.numeric(unlist(m.adj[[5]]))
          names(out3) <- paste("AdjMean", group_levels, sep = "_")

          # Post-hoc testing
          if (length(levels(colData[[GroupOfInterest]])) > 2) {
            if (identical(contrasts, "all")) {
              # Perform Tukey's HSD test for all pairwise comparisons
              postHocs <- glht(c1.final, linfct = mcp(GROUP = "Tukey"))
            } else if (identical(contrasts, "reference")) {
              # Perform Dunnett's test to compare Control vs. A and Control vs. B
              postHocs <- glht(c1.final, linfct = mcp(GROUP = "Dunnett"))
              message("Comparing all groups against", reference.level)
            } else if (any(grepl(":", contrasts))) {
              # Validate custom contrasts
              if (!all(sapply(contrasts, function(x) {
                parts <- strsplit(x, ":")[[1]]
                length(parts) == 2
              }))) {
                stop("Invalid custom contrast specification. Each custom contrast should be in the form 'level1:level2'.")
              }

              # Check if custom contrast levels are valid
              custom_levels <- unlist(sapply(contrasts, function(x) strsplit(x, ":")[[1]]))
              invalid_levels <- setdiff(custom_levels, group_levels)
              if (length(invalid_levels) > 0) {
                stop(paste("Invalid levels in custom contrasts:", paste(invalid_levels, collapse = ", ")))
              }

              # Create the contrast matrix
              contrast_matrix <- create_contrast_matrix(c1.final, contrasts, reference.level, n.levels)

              # Perform the post-hoc test with the custom contrast matrix
              postHocs <- glht(c1.final, linfct = mcp(GROUP = contrast_matrix))

            } else {
              stop("Invalid contrast specification. Please provide 'all', 'reference', or custom contrasts in the form 'level1:level2'.")
            }

            postHocs.s <- summary(postHocs)
            postHocs.c <- confint(postHocs)
            ph.p <- postHocs.s$test$pvalues
            ph.t <- postHocs.s$test$tstat
            ph.s <- postHocs.s$test$sigma
            ph.ci <- postHocs.c$confint
            compterms <- dimnames(ph.ci)[[1]]

            compterms <- sub(" - ", "_vs_", compterms)
            compterms <- paste("PostHoc_", compterms, sep = "")
            out4 <- sapply(seq_along(compterms), function(j) {
              c(ph.ci[j, ], ph.s[j], ph.t[j], ph.p[j])
            })
            out4 <- as.vector(unlist(out4))

            # Define the suffixes for each term
            suffixes <- c("_AdjMeanDiff", "_AdjMeanDiff_lwr", "_AdjMeanDiff_upr", "_stder", "_Fstat", "_Pvalue")
            # Generate the names directly
            names(out4) <- paste(rep(compterms, each = length(suffixes)), rep(suffixes, times = length(compterms)), sep = "")

            # Combine results
            ttt <- c(dimnames(Final.working.data)[[1]][i], c(out1, out2, out3, out4))
            names(ttt)[1] <- "Index"
            running.posthoc <- rbind(running.posthoc, ttt)
          } else {
            out4 <- out3[1] - out3[2]
            names(out4) <- "AdjMeanDiff"
            ttt <- c(dimnames(Final.working.data)[[1]][i], c(out1, out2, out3, out4))
            names(ttt)[1] <- "Index"
            running.posthoc <- rbind(running.posthoc, ttt)
          }
        }
      }
    }
  }

  # Adjust column names and perform multiple testing correction
  pvals <- as.numeric(running.posthoc[, "P.value"])
  fdr <- p.adjust(pvals, method = "BH")
  running.posthoc <- cbind(running.posthoc, FDR = fdr)
  rownames(running.posthoc) <- NULL

  # Rearrange columns
  posthoc_cols <- c("Index", "Fstat", "P.value", "FDR", setdiff(colnames(running.posthoc), c("Index", "P.value", "FDR", "Fstat")))
  running.posthoc <- running.posthoc[, posthoc_cols]

  InternalList <- list(
    means = means,
    sds = sds,
    cvs = cvs,
    lowess = lowess_results,
    keepers = keepers,
    cv.pick = cv.pick,
    mean.pick = mean.pick,
    PostNoiseSoftMeans = means.post
  )

  object@NormalizedCounts <- Final.working.data
  object@.InternalNormList <- InternalList
  object@AdjustedCounts <- as.data.frame(running.adj.data)
  object@Results <- as.data.frame(running.posthoc)
  names(partial_R2) <- rownames(Final.working.data)
  object@.InternalNormList$Partial_R2<-partial_R2

  if (random.covar) {
    flagged.genes <- object@Results[grepl("random", object@Results$Model, ignore.case = TRUE), "Index"]
    object@.InternalNormList[["CovarFlaggedGenes"]] <- flagged.genes
    object@Results$CovarFlaggedGenes <- ifelse(object@Results$Index %in% flagged.genes, TRUE, FALSE)
    if(shuffle.group==F){
      percentage_flagged <- (length(flagged.genes) / dim(running.posthoc)[1]) * 100

    }
  } else {
    object@.InternalNormList[["CovarFlaggedGenes"]] <- NULL
  }

  if (shuffle.group) {
    flagged.genes.2 <- object@Results[grepl("shuffled", object@Results$Model, ignore.case = TRUE), "Index"]
    object@.InternalNormList[["ShuffledFlaggedGenes"]] <- flagged.genes.2
    object@Results$ShuffledFlaggedGenes <- ifelse(object@Results$Index %in% flagged.genes.2, TRUE, FALSE)
    if(random.covar==F){
      percentage_flagged <- (length(flagged.genes) / dim(running.posthoc)[1]) * 100

    }
  } else {
    object@.InternalNormList[["ShuffledFlaggedGenes"]] <- NULL
  }
  if(random.covar && shuffle.group == T){
    unique.flagged.genes<-unique(c(flagged.genes.2, flagged.genes))
    percentage_flagged <- (length(unique.flagged.genes) / dim(running.posthoc)[1]) * 100

  }
  return(object)

}
