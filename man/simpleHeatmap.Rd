% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/simpleHeatmap.R
\name{simpleHeatmap}
\alias{simpleHeatmap}
\title{Simple Heatmap}
\usage{
simpleHeatmap(
  heatmap_dat,
  row_annot_df,
  col_annot_df,
  value_str = "value",
  annot_cols = NULL,
  order_by = c(NA, NA),
  leg_text_size = 14,
  leg_side = "left",
  row_names = FALSE,
  col_names = FALSE,
  cmap = mako(100),
  ...
)
}
\arguments{
\item{heatmap_dat}{Dataframe containing the data values to plot. Should include rowname and columns matching row_annot_df and col_annot_df}

\item{row_annot_df}{Dataframe of row annotations according to `dict2Annotdf`, rownames are `heatmap_dat` row names}

\item{col_annot_df}{Dataframe of column annotations according to `dict2Annotdf`, rownames are `heatmap_dat` row names}

\item{value_str}{Character string, name of the value being plotted}

\item{annot_cols}{Named list of all the colour maps for annotations in `row_annot_df` and `col_annot_df`}

\item{order_by}{A vector of two elements, either character strings or NA, describing if rows and/or columns should be ordering according to the input `heatmap_dat` i.e. c(row,col). If NA, clustering is performed. Assumes annotation is a factor as per `GetOrdering`}

\item{leg_text_size}{Integer, size of legend text}

\item{leg_side}{Character string, side to put heatmap legend, either "right", "left", "bottom" or "top"}

\item{row_names}{Boolean, if TRUE displays row names}

\item{col_names}{Boolean, if TRUE displays column names}

\item{cmap}{Character vector of hex colours used to colour heatmap}
}
\value{
Function returns a HeatmapList-class object (ComplexHeatmap), a plot is created
}
\description{
A wrapper for ComplexHeatmap with the intention of
simplifying the coordination of row and column annotations
}
\examples{
## -- Run dict2Annotdf(dict_lst, data_names) to get row/col annot_dfs. Must include row/column annotations
col_annot <- dict2Annotdf(list(SigType=sig_group, Author=SigAuthor), colnames(heatmap_dat))
row_annot <- dict2Annotdf(list(CellType=setNames(tmp_dge$samples$CellType,rownames(tmp_dge$samples)),
                              Cluster=setNames(tmp_dge$samples$Cluster,rownames(tmp_dge$samples)),
                              Stage= setNames(tmp_dge$samples$Stage, rownames(tmp_dge$samples))),
                        rownames(heatmap_dat))
simpleHeatmap(heatmap_dat,row_annot,col_annot)

## -- Define colours
col_lst <- list(SigType=sig_cols,Author=Author_cols,
                   CellType=GA_cols,Cluster=NCC_clust_cols,Stage=TS_cols)
simpleHeatmap(heatmap_dat,row_annot,col_annot, annot_cols=col_lst)

## -- Different ordering
simpleHeatmap(heatmap_dat,row_annot,col_annot,order_by=c(NA, "SigType"))
simpleHeatmap(heatmap_dat,row_annot,col_annot,order_by=c("CellType", NA))
simpleHeatmap(heatmap_dat,row_annot,col_annot,order_by=c("CellType", "SigType"))

}
