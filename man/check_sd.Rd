% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/check_sd.R
\name{check_sd}
\alias{check_sd}
\title{Check Standard Deviation and Display Warnings}
\usage{
check_sd(data_frame, threshold = 5)
}
\arguments{
\item{data_frame}{A data frame containing the variables to be checked.}

\item{threshold}{A numeric value specifying the threshold for SD. Columns with SD greater than this value
will trigger a warning. Default is 5.}
}
\value{
The original data frame. The function returns the data frame invisibly and does not modify it.
}
\description{
This function checks the standard deviation (SD) of numeric columns in a data frame and displays warnings
if any column has an SD greater than a specified threshold. High SD values can affect the fit of statistical
models like ANCOVA, and the function suggests scaling these variables if needed.
}
\details{
The function first ensures that the input is a data frame. It then identifies numeric columns and
calculates their standard deviations. If any numeric column has an SD greater than the specified threshold,
a warning message is printed, listing these variables. If no columns exceed the threshold, the function does
nothing.
}
\keyword{internal}
