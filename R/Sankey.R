


#' MakeSankey
#'
#' Draws a Sankey diagram using networkD3::sankeyNetwork().
#' Input vectors must be same length.
#'
#' @param ListOfLabels_1 Character vector of labels
#' @param ListOfLabels_2 Character vector of labels
#' @param ListOfLabels_3 Character vector of labels
#' @param loops TRUE or FALSE, whether or not to allow self-loops in a label common to both lists
#'
#' @return A sankeyNetwork() plot
#' @export
#'
#' @importFrom networkD3 sankeyNetwork
#' @import dplyr
#'
MakeSankey <- function(ListOfLabels_1, ListOfLabels_2,
                       ListOfLabels_3=NULL, loops=FALSE){

  if (!(is.null(ListOfLabels_3))){
    MakeSankey3d(ListOfLabels_1, ListOfLabels_2, ListOfLabels_3)
  } else {

    if (!(loops)){
      ## Get unique names
      int <- intersect(unique(ListOfLabels_1), unique(ListOfLabels_2))
      if (length(int) > 0){
        dict_1 <- setNames(paste0(int,"_1"),int)
        dict_2 <- setNames(paste0(int,"_2"),int)
        dict_1[setdiff(unique(ListOfLabels_1), int)] <- setdiff(unique(ListOfLabels_1), int)
        dict_2[setdiff(unique(ListOfLabels_2), int)] <- setdiff(unique(ListOfLabels_2), int)

        ListOfLabels_1 <- as.vector(dict_1[ListOfLabels_1])
        ListOfLabels_2 <- as.vector(dict_2[ListOfLabels_2])
      }
    }

    ## Build weighted table
    Sankey_df <- data.frame(Labels_1=ListOfLabels_1,
                            Labels_2=ListOfLabels_2)
    Weighted_table <- table(Sankey_df[, c("Labels_1", "Labels_2")])

    source_vector <- c()
    target_vector <- c()
    count_vector <- c()

    for (i in rownames(Weighted_table)){
      for (j in colnames(Weighted_table)){
        iWeight <- Weighted_table[i, j]
        if (iWeight != 0){
          source_vector <- append(source_vector, i)
          target_vector <- append(target_vector, j)
          count_vector <- append(count_vector, iWeight)
        }
      }
    }

    SankeyLinks_df <- data.frame(Source=source_vector,
                                 Target=target_vector,
                                 Count=count_vector)
    nodes <- data.frame(name=c(as.character(SankeyLinks_df$Source), as.character(SankeyLinks_df$Target)) %>% unique())

    SankeyLinks_df$IDsource <- match(SankeyLinks_df$Source, nodes$name)-1
    SankeyLinks_df$IDtarget <- match(SankeyLinks_df$Target, nodes$name)-1

    sankeyNetwork(Links = SankeyLinks_df, Nodes = nodes,
                  Source = "IDsource", Target = "IDtarget",
                  Value = "Count", NodeID = "name",
                  sinksRight=FALSE, fontSize=14, fontFamily="Helvetica")
  }

}

#' MakeSankey3d
#'
#' @param ListOfLabels_1 Character vector of labels
#' @param ListOfLabels_2 Character vector of labels
#' @param ListOfLabels_3 Character vector of labels
#'
#' @return A sankeyNetwork() plot
#' @export
#'

MakeSankey3d <- function(ListOfLabels_1, ListOfLabels_2, ListOfLabels_3){

  Sankey_df <- data.frame(Labels_1=as.character(ListOfLabels_1),
                          Labels_2=as.character(ListOfLabels_2),
                          Labels_3=as.character(ListOfLabels_3))
  Sankey_df <- Sankey_df[complete.cases(Sankey_df), ]

  Weighted_table_1to2 <- table(Sankey_df[, c("Labels_1", "Labels_2")])
  Weighted_table_2to3 <- table(Sankey_df[, c("Labels_2", "Labels_3")])

  source_vector <- c()
  target_vector <- c()
  count_vector <- c()

  for (i in rownames(Weighted_table_1to2)){
    for (j in colnames(Weighted_table_1to2)){
      iWeight <- Weighted_table_1to2[i, j]
      if (iWeight != 0){
        source_vector <- append(source_vector, i)
        target_vector <- append(target_vector, j)
        count_vector <- append(count_vector, iWeight)
      }
    }
  }

  for (i in rownames(Weighted_table_2to3)){
    for (j in colnames(Weighted_table_2to3)){
      iWeight <- Weighted_table_2to3[i, j]
      if (iWeight != 0){
        source_vector <- append(source_vector, i)
        target_vector <- append(target_vector, j)
        count_vector <- append(count_vector, iWeight)
      }
    }
  }

  SankeyLinks_df <- data.frame(Source=source_vector,
                               Target=target_vector,
                               Count=count_vector)
  nodes <- data.frame(name=c(as.character(SankeyLinks_df$Source), as.character(SankeyLinks_df$Target)) %>% unique())

  SankeyLinks_df$IDsource <- match(SankeyLinks_df$Source, nodes$name)-1
  SankeyLinks_df$IDtarget <- match(SankeyLinks_df$Target, nodes$name)-1

  networkD3::sankeyNetwork(Links = SankeyLinks_df, Nodes = nodes,
                     Source = "IDsource", Target = "IDtarget",
                     Value = "Count", NodeID = "name",
                     sinksRight=FALSE,fontSize =14, fontFamily="Helvetica")


}


