% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotNoiseCorrection.R
\name{plotNoiseCorrection}
\alias{plotNoiseCorrection}
\title{Plot Pre and Post Noise Correction Data
This function generates plots to visually inspect the quality of noise correction in the data.
It creates two plots: one showing the data post-noise filtering and one pre-noise filtering,
with LOESS fits for the different groups. The plots help in assessing the effectiveness of the noise
correction and identifying any remaining issues.}
\arguments{
\item{object}{An object of class \code{FalconNest} which contains slots with normalized counts and other relevant data.}
}
\value{
A grid of plots is displayed, including:
\itemize{
\item A plot showing post noise filtering LOESS with points colored by group.
\item A plot showing pre noise filtering LOESS with points and LOESS smooth lines.
}
The function does not return a value but displays the plots directly.
}
\description{
Plot Pre and Post Noise Correction Data
This function generates plots to visually inspect the quality of noise correction in the data.
It creates two plots: one showing the data post-noise filtering and one pre-noise filtering,
with LOESS fits for the different groups. The plots help in assessing the effectiveness of the noise
correction and identifying any remaining issues.
}
\details{
The function checks that the input object is of class \code{FalconNest} and contains the required slots.
It then extracts relevant data such as normalized counts, internal normalization lists, and group variables.
It creates two types of plots: one showing post-noise filtering data and another showing pre-noise filtering data.
Both plots include dashed lines for the mean and coefficient of variation thresholds. The plots are arranged
side-by-side with a title indicating that they are noise correction plots.
}
\examples{
colData <- data.frame(sample = factor(rep(c("A", "B"), each = 5)))
counts <- matrix(data = abs(rpois(50 * 10, lambda = 5)), nrow = 50, ncol = 10)
colnames(counts) <- rownames(colData)
rownames(counts) <- 1:50

# Initialize FalconNest object
falcon_obj <- FalconNest(counts = counts, colData = colData, design = ~sample)

# Normalize counts
normalizedObject <- normCounts(object = falcon_obj,
group = "sample", mean.pick = 1, cv.pick = 100)
plotNoiseCorrection(normalizedObject)

}
