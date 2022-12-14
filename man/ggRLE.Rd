% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/PlotFunc_R.R
\name{ggRLE}
\alias{ggRLE}
\title{ggRLE}
\usage{
ggRLE(
  dat_x,
  annot,
  col_str,
  col_pal = NULL,
  isLog = TRUE,
  isLarge = FALSE,
  ylim = c(-2, 2),
  zero_line = TRUE,
  zero_col = "skyblue",
  medPoint = FALSE,
  whisk = 1.5
)
}
\arguments{
\item{dat_x}{DGEList object or subsettable data.frame or matrix object}

\item{annot}{DGEList$samples or data.frame object containing sample annotations, including col_str}

\item{col_str}{Character string indicating column of annot to colour boxplots by}

\item{col_pal}{Named character vector to use as colour palette, where values are hex colours and names match col_str}

\item{isLog}{Defaults to TRUE, set to FALSE if dat_x is not log-transformed}

\item{isLarge}{Defaults to TRUE, if FALSE removes boxplot whiskers}

\item{ylim}{Numeric vector indicating y-axis limits to plot, defaults to c(-2,2)}

\item{zero_line}{Defaults to TRUE, if FALSE does not plot line at y=0}

\item{zero_col}{Defaults to "skyblue", indicating colour for zero_line}

\item{medPoint}{Defaults to FALSE, if TRUE plots the median value as a point}

\item{whisk}{Defaults to 1.5, numeric indicating length of whiskers as multiple of IQR}
}
\value{
A ggplot object
}
\description{
Calculates RLE plot and returns ggplot object
}
\examples{
ggRLE(dge_x, dge_x$samples, "Treatment", isLog = FALSE)
ggRLE(dge_x$logCPM, dge_x$samples, "Patient", isLog = TRUE)
}
