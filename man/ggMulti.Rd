% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/PlotFunc_R.R
\name{ggMulti}
\alias{ggMulti}
\title{ggMulti}
\usage{
ggMulti(
  plot_lst,
  leg_to_plot = NULL,
  leg_lst = NULL,
  n_col = 2,
  leg_title_pos = "top",
  leg_size = 4,
  leg_nrow = 2,
  plot_ratio = 6,
  ...
)
}
\arguments{
\item{plot_lst}{List of ggplots to plot}

\item{leg_to_plot}{Character vector or vector of indexes indicating the plots in `plot_lst` to include legends for. If leg_lst=NULL, leg_to_plot is not used}

\item{leg_lst}{List of legends (gtable objects) from cowplot::get_legend() to include}

\item{n_col}{Integer indicating the number of columns that `plot_lst` should be organised into}

\item{leg_title_pos}{Character string indicating position of legend titles within each legend gtable object, one of "top", "bottom" etc.}

\item{leg_size}{Integer size of legend keys}

\item{leg_nrow}{Integer, number of rows of keys in each legend, passed on to `ExtractLegends()`}

\item{plot_ratio}{Integer value that will be the height of `cowplot::plot_grid(plot_lst)` relative to the single row of legends at the bottom}

\item{...}{}
}
\value{
Plots with cowplot::plot_grid()
}
\description{
ggMulti
}
\examples{

## Only plot legend for first plot
ggMulti(list(A,B),leg_to_plot=c(1),leg_nrow=1)

## Include two custom legends
ggMulti(list(A, B), leg_lst = list(GetLegend(Colour_values,"CellType", n_row=2),
                                   GetShapeLegend(unique(dge_x$samples$Sample), "Sample")))

}
