#' @title CircularPlot
#' @description a function to draw CircularPlot
#' @author Shine Shen
#' \email{qq951633542@@163.com}
#' @param label_size 5
#' @return  All the results can be got form other functions and instruction.
#' @export
#' @examples
#' \donttest{
#' ##---- Be sure the format of data and sample.info is correct!! ----
#' }
CircularPlot <- function(label_size=5){

  #library(tidyverse)

  # Create dataset
  data <- read.csv("cicular data.csv")
  data$value <- -log2(data$value)
  # Set a number of 'empty bar' to add at the end of each group
  empty_bar <- 3
  to_add <- data.frame( matrix(NA, empty_bar*nlevels(data$group), ncol(data)) )
  colnames(to_add) <- colnames(data)
  to_add$group <- rep(levels(data$group), each=empty_bar)
  data <- rbind(data, to_add)
  data <- data %>% arrange(group)
  data$id <- seq(1, nrow(data))

  # Get the name and the y position of each label
  label_data <- data
  number_of_bar <- nrow(label_data)
  angle <- 90 - 360 * (label_data$id-0.5) /number_of_bar     # I substract 0.5 because the letter must have the angle of the center of the bars. Not extreme right(1) or extreme left (0)
  label_data$hjust <- ifelse( angle < -90, 1, 0)
  label_data$angle <- ifelse(angle < -90, angle+180, angle)

  # prepare a data frame for base lines
  base_data <- data %>%
    group_by(group) %>%
    summarize(start=min(id), end=max(id) - empty_bar) %>%
    rowwise() %>%
    mutate(title=mean(c(start, end)))

  # prepare a data frame for grid (scales)
  grid_data <- base_data
  grid_data$end <- grid_data$end[ c( nrow(grid_data), 1:nrow(grid_data)-1)] + 1
  grid_data$start <- grid_data$start - 1
  grid_data <- grid_data[-1,]

  # Make the plot
  p <- ggplot(data, aes(x=as.factor(id), y=value, fill=individual))+     # Note that id is a factor. If x is numeric, there is some space between the first bar

    geom_bar(aes(x=as.factor(id), y=value, fill=individual), stat="identity")+

    # Add a val=100/75/50/25 lines. I do it at the beginning to make sur barplots are OVER it.
    #geom_segment(data=grid_data, aes(x = end, y = 20, xend = start, yend = 20), colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
    #geom_segment(data=grid_data, aes(x = end, y = 15, xend = start, yend = 15), colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
    geom_segment(data=grid_data, aes(x = end, y = 10, xend = start, yend = 10), colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
    geom_segment(data=grid_data, aes(x = end, y = 5, xend = start, yend = 5), colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +

    # Add text showing the value of each 100/75/50/25 lines
    #annotate("text", x = rep(max(data$id),4), y = c(20, 40, 60, 80), label = c("20", "40", "60", "80") ,
    #color="grey", size=3 , angle=0, fontface="bold", hjust=1) +
    annotate("text", x = rep(max(data$id),2), y = c(5, 10), label = c("5", "10") ,color="grey",
             size= 2.5 , angle=0, fontface="bold", hjust=1) +
    #geom_bar(aes(x=as.factor(id), y=value, fill=individual), stat="identity", alpha=0.5) +
    ylim(-100,120) +
    theme_minimal() +
    theme(legend.position = 'none',
     #legend.text = element_text(size = 14),
      axis.text = element_blank(),
      axis.title = element_blank(),
      panel.grid = element_blank(),
      plot.margin = unit(rep(-1,4), "cm")
    ) +
    coord_polar() +
     geom_text(data=label_data, aes(x=id, y=value, label=individual, hjust=hjust), color="black",
              fontface="bold",alpha=0.6, size=label_size, angle= label_data$angle, inherit.aes = FALSE ) +

    # Add base line information
    geom_segment(data=base_data, aes(x = start, y = -5, xend = end, yend = -5), colour = "black", alpha=0.8, size=0.6 , inherit.aes = FALSE )  +
    geom_text(data=base_data, aes(x = title, y = -18, label=group), hjust=c(1,1,0), colour = "black", alpha=0.8, size=4, fontface="bold", inherit.aes = FALSE)

  export::graph2ppt(x=p,file='circular.pptx',height=8,width=10)
}

