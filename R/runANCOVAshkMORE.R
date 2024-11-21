#' @name FalconGlide
#' @title Falcon Function for ANCOVA Analysis
#'
#' @description This function performs ANCOVA analysis on the normalized counts data from a
#' `FalconNest` object.
#'
#' @param object An object of class `FalconNest`, containing counts data, sample metadata, and normalized counts.
#' @param contrasts A character string specifying the type of contrasts for post-hoc
#' testing. Defaults to "all" for Tukey's HSD test. Other options include "reference" for Dunnett's test or custom contrasts.
#' @param random.covar Logical, indicating whether to apply random covariate control. Defaults to TRUE.
#' @param shuffle.group Logical, indicating whether to shuffle the group variable for control. Defaults to TRUE.
#' @param FalconOutput Logical, whether to capture the input parameters. Defaults to TRUE.
#'
#' @return An updated `FalconNest` object with two additional slots:
#' \item{AdjustedCounts}{A data frame of adjusted counts after ANCOVA.}
#' \item{Results}{A data frame of ANCOVA results, including post-hoc analysis, if applicable.}
#'
#' @details
#' This function runs ANCOVA on each gene in the `NormalizedCounts` slot of the `FalconNest` object.
#' It builds the ANCOVA model based on the provided design and group variables, with options for
#' post-hoc analysis using Tukey's HSD, Dunnett's test, or custom contrasts. The function also
#' adjusts for multiple comparisons using the Benjamini-Hochberg method and can perform controls
#' for random covariates or shuffled group variables.
#'
#' @importFrom stats aov step var deviance as.formula aggregate confint p.adjust
#' @importFrom progress progress_bar
#' @importFrom car leveneTest Anova vif
#' @importFrom effects effect
#' @importFrom multcomp glht mcp
#' @importFrom dplyr mutate
#' @import magrittr
#' @import stringr
#' @keywords internal

utils::globalVariables(c("progress_bar", "reference.level", "median", "n.levels", "c1.final", "VIF.final"))

runANCOVAshkMORE <- function(object, contrasts = "all", p.adjust="BH", random.covar = TRUE, shuffle.group = TRUE, FalconOutput=TRUE, verbose=TRUE) {

  set.seed(07142016)  # Set a seed for reproducibility

  # Ensure the input is of the correct class 'FalconNest'
  if (!inherits(object, "FalconNest")) {
    stop("Input must be of class 'FalconNest'")
  }

  # Check if contrasts are provided, stop if none specified
  if (missing(contrasts) || length(contrasts) == 0) {
    stop("No contrasts specified.")
  }

  # Extract relevant data from the input object
  FalconInfo <- object@colData
  FalconDesign <- object@.design
  Final.working.data <- object@NormalizedCounts
  GroupOfInterest <- object@.GroupVariable
  variables <- all.vars(FalconDesign)
  group_levels <- unique(FalconInfo[[GroupOfInterest]])
  reference.level <- levels(FalconInfo[[GroupOfInterest]])[1]
  n.levels <- levels(FalconInfo[[GroupOfInterest]])
  group<-object@.GroupVariable
  partial_R2<-list()
  smoothed_variance<-object@.InternalNormList$smoothed_variance
  
  if (contrasts == "all") {
    message("Comparing all groups against ", reference.level)
  }

  # Convert GroupOfInterest to factor if it's neither numeric nor factor
  if (!is.numeric(FalconInfo[[GroupOfInterest]]) && !is.factor(FalconInfo[[GroupOfInterest]])) {
    FalconInfo[[GroupOfInterest]] <- as.factor(FalconInfo[[GroupOfInterest]])
  }

  # Check if results already exist and handle flagged genes if present
  if (!is.null(dim(object@Results[1]))) {

    if (is.null(object@.InternalNormList$CovarFlaggedGenes) &&
        is.null(object@.InternalNormList$ShuffledFlaggedGenes)) {
      stop("Falcon analysis has already been completed. No flagged genes were identified for removal.")
    }

    # Remove CovarFlaggedGenes from the analysis
    if (!is.null(object@.InternalNormList$CovarFlaggedGenes)) {
      covar.flagged.genes <- object@.InternalNormList$CovarFlaggedGenes
      Final.working.data <- Final.working.data[!(row.names(Final.working.data) %in% covar.flagged.genes), ]
      if(verbose) {
      message(paste("Warning:", length(covar.flagged.genes), "'CovarFlaggedGenes' will be removed from the analysis."))
   }
      # Disable random.covar if applied already
      if(random.covar){
        random.covar=FALSE
        message(paste("Randomization control option cannot be applied more than once. 'random.covar' has been set to FALSE"))
      }
    }

    # Remove ShuffledFlaggedGenes from the analysis
    if (!is.null(object@.InternalNormList$ShuffledFlaggedGenes)) {
      shuffled.flagged.genes <- object@.InternalNormList$ShuffledFlaggedGenes
      Final.working.data <- Final.working.data[!(row.names(Final.working.data) %in% shuffled.flagged.genes), ]
      if(verbose) {
      message(paste("Warning:", length(shuffled.flagged.genes), " 'ShuffledFlaggedGenes' will be removed  from the analysis"))
      }
      # Disable shuffle.group if applied already
      if(shuffle.group){
        shuffle.group=FALSE
        message(paste("Warning: Randomization control option cannot be applied more than once. 'shuffle.group' has been set to FALSE."))
      }
    }
  }
  
  
  FalconInfo.binary <- FalconInfo
  FalconInfo.binary <- model.matrix(FalconDesign, data = FalconInfo.binary)
  names(FalconInfo.binary)[names(FalconInfo.binary) == GroupOfInterest] <- "GROUP"

  # Clean up column names
  colnames(FalconInfo.binary) <- gsub("[/ -]", "_", colnames(FalconInfo.binary))
  binary.names <- setdiff(as.vector(colnames(FalconInfo.binary)), as.vector(colnames(FalconInfo)))

  # Initialize progress bar
  pb <- progress_bar$new(
    format = "  Processing [:bar] :percent ETA: :eta",
    total = nrow(Final.working.data), clear = FALSE, width = 60
  )

  # Prepare ANCOVA formula
  components <- c(variables[variables %in% names(FalconInfo.binary)], binary.names)

  running.posthoc <- NULL
  running.adj.data <- NULL

  components[components == GroupOfInterest] <- "GROUP"
  fixed_data <- FalconInfo.binary[, components, drop = FALSE]

  # Add random covariate if specified
  if (random.covar) {
    ANCOVA.formula <- paste("DATA ~", paste(c(components, "random.covar"), collapse = " + "))
    fixed_data <- fixed_data %>%
      mutate(random.covar = rnorm(nrow(fixed_data), 0, 1))  # Add random covariate
    object@.InternalNormList[["random.covar"]] <- fixed_data$random.covar
  } else {
    ANCOVA.formula <- paste("DATA ~", paste(components, collapse = " + "))
  }

  # Shuffle group if specified
  if (shuffle.group) {
    shuffled.group<-data.frame(shuffled.group=as.factor(sample(FalconInfo[[group]])))
    if(length(levels(shuffled.group$shuffled.group)) == 2){
      ANCOVA.formula <- paste(c(ANCOVA.formula,"shuffled.group"), collapse = " + ")
    fixed_data <- fixed_data %>%
      mutate(shuffled.group = unlist(apply_conversion(shuffled.group, "shuffled.group", verbose=F)))
    fixed_data$shuffled.group <- as.numeric(fixed_data$shuffled.group)
    }else{
      conversion_result <- apply_conversion(shuffled.group, "shuffled.group", verbose=FALSE)
      fixed_data <- fixed_data %>%
        bind_cols(as.data.frame(conversion_result))
      ANCOVA.formula <- paste(c(ANCOVA.formula, names(conversion_result)), collapse = " + ")
    }
    object@.InternalNormList[["shuffled.group"]] <- fixed_data$shuffled.group
  }

  ANCOVA.formula.transformed <- as.formula(ANCOVA.formula)
  design.limma <- model.matrix(FalconDesign, data = FalconInfo)
  fit <- lmFit(Final.working.data, design.limma)
  eb <- eBayes(fit)
  mod.coefs <- eb$coefficients

  # Process each gene
  for (i in seq_len(nrow(Final.working.data))) {
    pb$tick()  # Update progress bar
    perm.temp <- fixed_data
    perm.temp <- model.matrix(~ ., data = fixed_data)
    perm.temp$DATA <- as.numeric(t(Final.working.data[i, ]))

    # Proceed if there is variation in the data
    if (var(perm.temp$DATA) > 0) {
      # Fit the initial ANCOVA model
      c1.init <- aov(ANCOVA.formula.transformed, data = perm.temp)

      # Stepwise model selection
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

        # Check for model significance
        if (deviance(c1.step) > sqrt(.Machine$double.eps)) {
          ccc.report <- ccc.report[-length(ccc.report)]
          ccc.report <- paste(ccc.report, collapse = "+")
          ccc.report.final <- gsub("\\+random.covar$", "", ccc.report)

          # Fit the final model
          weights <- rep(1 / smoothed_variance[i], dim(perm.temp)[1])
          
          ccc <- paste("c1.final <- lm(DATA ~", ccc.report.final, ", data = perm.temp", ", weights = weights)")
          eval(parse(text = ccc))
          
          mod.coefs.i <- mod.coefs[i,]
          c1.final$coefficients[grep(GroupOfInterest, names(c1.final$coefficients), value = F, ignore.case = TRUE)]<-mod.coefs.i[grep(GroupOfInterest, names(mod.coefs.i), value = F, ignore.case = TRUE)]
      
          
          # Perform ANOVA
          temp <- car::Anova(c1.final, type = "III", singular.ok = TRUE)

          # Conduct Levene's Test
          l.test <- leveneTest(perm.temp$DATA, perm.temp$GROUP, center = median)
          selected.covariates <- unlist(str_split(ccc.report.final, "\\+"))

          # Extract sums of squares
          ss_covariates <- temp$`Sum Sq`[2:(nrow(temp))]  # Exclude intercept and residuals

          # Calculate partial R²
          gene_partial_r2 <- ss_covariates / sum(ss_covariates)

          # Create a data frame for partial R²
          partial_r2_df <- data.frame(
            Covariate = rownames(temp)[-1],  # Exclude intercept
            Partial_R2 = c(gene_partial_r2)
          )

          # Find matching root names
          partial_r2_df$Root <- sapply(partial_r2_df$Covariate, function(covariate) {
            match_index <- which(sapply(variables, function(root_var) {
              grepl(paste0("^", root_var), covariate)
            }))
            if (length(match_index) > 0) {
              return(variables[match_index])
            } else {
              return(covariate)  # No match found
            }
          })

          # Sum Partial_R2 values for covariates with the same root name
          summed_partial_r2_df <- aggregate(Partial_R2 ~ Root, data = partial_r2_df, sum)

          # Store results
          partial_R2[[i]] <- summed_partial_r2_df

          # VIF Calculation
          
          if (length(selected.covariates) > 1) {
            ccc.VIF <- paste("VIF.final <- lm(DATA ~", ccc.report.final, ", data = perm.temp)")
            eval(parse(text = ccc.VIF))  # Evaluate the VIF formula
            vif_results <- vif(VIF.final)  # Compute VIF
            high_vif_terms <- names(vif_results)[vif_results > 8]  # Check for collinearity

            # Find original variables for high VIF terms
            high_vif_terms <- sapply(high_vif_terms, function(suffix_var) {
              match_index <- which(sapply(variables, function(orig_var) {
                grepl(paste0("^", orig_var), suffix_var)
              }))
              if (length(match_index) > 0) {
                return(variables[match_index])
              } else {
                return(suffix_var)  # Return the suffixed variable
              }
            })

            if (length(high_vif_terms) > 0) {
              colinear_covariate <- paste(high_vif_terms, collapse = "+")
            } else {
              colinear_covariate <- "FALSE"
            }
          } else {
            colinear_covariate <- "NA"
          }
          

          # Collect results
          group.position<-which(rownames(temp) == "GROUP")
          out2 <- c(ccc.report, l.test$'Pr(>F)'[1], colinear_covariate, temp[[1]][group.position], temp[[2]][group.position], temp[[3]][group.position], temp[[4]][group.position])
          names(out2) <- c("Model", "Levene_P.value", "Colinearity", "SumSq", "DF", "Fstat", "P.value")

          # Adjust coefficients if necessary
          if (length(c1.final[[1]]) > length(levels(perm.temp$GROUP))) {
            coef_names <- names(c1.final[[1]])  # Coefficients
            group_factors <- grep("^GROUP", coef_names, value = TRUE)
            adj.coefficients <- c1.final[[1]][!coef_names %in% c("(Intercept)", group_factors)] * -1  # Exclude intercept and group factor coefficients

            temp2 <- perm.temp[, names(adj.coefficients)]
            temp3 <- as.numeric(adj.coefficients)

            if (length(temp3) > 1) {
              for (j in seq_along(temp3)) {
                temp2[, j] <- (as.numeric(temp2[, j]) - 1) * temp3[j]
              }
              dimnames(temp2)[[2]] <- names(adj.coefficients)
            } else {
              temp2 <- as.numeric(temp2) * temp3
              temp2 <- matrix(temp2, , 1)
              dimnames(temp2)[[2]] <- names(adj.coefficients)
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
            mean(perm.temp$DATA[FalconInfo[[GroupOfInterest]] == group.index], na.rm = TRUE)
          })
          names(out1) <- paste("UnadjMean", group_levels, sep = "_")

          m.adj <- effect("GROUP", c1.final)
          out3 <- as.numeric(unlist(m.adj[[5]]))
          names(out3) <- paste("AdjMean", group_levels, sep = "_")

          # Post-hoc testing
          if (length(levels(FalconInfo[[GroupOfInterest]])) > 2) {
            if (identical(contrasts, "all")) {
              postHocs <- glht(c1.final, linfct = mcp(GROUP = "Tukey"))  # Tukey's HSD
            } else if (identical(contrasts, "reference")) {
              postHocs <- glht(c1.final, linfct = mcp(GROUP = "Dunnett"))
            } else if (any(grepl(":", contrasts))) {
              # Validate custom contrasts
              if (!all(sapply(contrasts, function(x) {
                length(strsplit(x, ":")[[1]]) == 2
              }))) {
                stop("Invalid custom contrast specification. Each custom contrast should be in the form 'level1:level2'.")
              }

              # Check for valid contrast levels
              custom_levels <- unlist(sapply(contrasts, function(x) strsplit(x, ":")[[1]]))
              invalid_levels <- setdiff(custom_levels, group_levels)
              if (length(invalid_levels) > 0) {
                stop(paste("Invalid levels in custom contrasts:", paste(invalid_levels, collapse = ", ")))
              }

              # Create contrast matrix
              contrast_matrix <- create_contrast_matrix(c1.final, contrasts, reference.level, n.levels)

              # Perform post-hoc test
              postHocs <- glht(c1.final, linfct = mcp(GROUP = contrast_matrix))
            } else {
              stop("Invalid contrast specification. Please provide 'all', 'reference', or custom contrasts in the form 'level1:level2'.")
            }

            # Extract and format post-hoc results
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

            # Define suffixes for output
            suffixes <- c("_AdjMeanDiff", "_AdjMeanDiff_lwr", "_AdjMeanDiff_upr", "_stder", "_Fstat", "_Pvalue")
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
  fdr <- p.adjust(pvals, method = p.adjust)
  running.posthoc <- cbind(running.posthoc, FDR = fdr)
  rownames(running.posthoc) <- NULL

  # Rearrange columns for better organization
  posthoc_cols <- c("Index", "Fstat", "P.value", "FDR", setdiff(colnames(running.posthoc), c("Index", "P.value", "FDR", "Fstat")))
  running.posthoc <- running.posthoc[, posthoc_cols]

  # Store results in the FalconNest object
  object@AdjustedCounts <- as.data.frame(running.adj.data)
  object@Results <- as.data.frame(running.posthoc)
  names(partial_R2) <- rownames(Final.working.data)
  object@.InternalNormList$Partial_R2 <- partial_R2
  object@.InternalNormList$fixed_data <- fixed_data

  # Handle flagged genes based on random covariates
  if (random.covar && (!shuffle.group) ) {
    flagged.genes <- object@Results[grepl("random", object@Results$Model, ignore.case = TRUE), "Index"]
    object@.InternalNormList[["CovarFlaggedGenes"]] <- flagged.genes
    object@Results$CovarFlaggedGenes <- object@Results$Index %in% flagged.genes
      percentage_flagged <- (length(flagged.genes) / nrow(running.posthoc)) * 100
      if(verbose) {
      message(sprintf("Warning: %d (%.2f%%) of values have been flagged using random.covar selection.",
                      length(flagged.genes), percentage_flagged))
      }
  # Handle flagged genes based on shuffled groups
  }else if (shuffle.group && (!random.covar)) {
    flagged.genes.2 <- object@Results[grepl("shuffled", object@Results$Model, ignore.case = TRUE), "Index"]
    object@.InternalNormList[["ShuffledFlaggedGenes"]] <- flagged.genes.2
    object@Results$ShuffledFlaggedGenes <- object@Results$Index %in% flagged.genes.2
      percentage_flagged <- (length(flagged.genes) / nrow(running.posthoc)) * 100
      if(verbose) {
      message(sprintf("Warning: %d (%.2f%%) of values have been flagged using shuffle.group selection.",
                      length(flagged.genes), percentage_flagged))
      }


  }else if (random.covar && shuffle.group) {
    flagged.genes.2 <- object@Results[grepl("shuffled", object@Results$Model, ignore.case = TRUE), "Index"]
    flagged.genes <- object@Results[grepl("random", object@Results$Model, ignore.case = TRUE), "Index"]
    unique.flagged.genes <- unique(c(flagged.genes.2, flagged.genes))
    percentage_flagged <- (length(unique.flagged.genes) / nrow(running.posthoc)) * 100
    object@.InternalNormList[["ShuffledFlaggedGenes"]] <- flagged.genes.2
    object@.InternalNormList[["CovarFlaggedGenes"]] <- flagged.genes
    object@Results$CovarFlaggedGenes <- object@Results$Index %in% flagged.genes
    object@Results$ShuffledFlaggedGenes <- object@Results$Index %in% flagged.genes.2

    if(verbose) {
    message(sprintf("Warning: %d (%.2f%%) unique genes have been flagged using shuffled group and random covariate selection",
                    length(unique.flagged.genes), percentage_flagged))
    }
  } else {
    object@.InternalNormList[["ShuffledFlaggedGenes"]] <- NULL
  }


  # Store metadata in the FalconHatch object
  object@.FalconHatch <- c(object@.FalconHatch ,list(
    contrasts = contrasts,
    random.covar = random.covar,
    shuffle.group = shuffle.group,
    FinalGeneList = rownames(Final.working.data),
    FalconDesign = FalconDesign,
    FalconInfo = FalconInfo,
    GroupOfInterest = object@.GroupVariable,
    random.covar.list = object@.InternalNormList[["random.covar"]],
    shuffle.group.list = object@.InternalNormList[["shuffled.group"]]
  ))

  # Output the FalconHatch object if required
  if (FalconOutput) {
    dput(object@.FalconHatch, file = paste0(Sys.time(),".Falcon"))
  }

  return(object)
}
