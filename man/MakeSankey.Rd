% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Sankey.R
\name{MakeSankey}
\alias{MakeSankey}
\title{MakeSankey}
\usage{
MakeSankey(
  ListOfLabels_1,
  ListOfLabels_2,
  ListOfLabels_3 = NULL,
  loops = FALSE
)
}
\arguments{
\item{ListOfLabels_1}{Character vector of labels}

\item{ListOfLabels_2}{Character vector of labels}

\item{ListOfLabels_3}{Character vector of labels}

\item{loops}{TRUE or FALSE, whether or not to allow self-loops in a label common to both lists}
}
\value{
A sankeyNetwork() plot
}
\description{
Draws a Sankey diagram using networkD3::sankeyNetwork().
Input vectors must be same length.
}
