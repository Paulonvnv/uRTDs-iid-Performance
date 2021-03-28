# function to create a clustered heatmap plot ----

fx_HClust.Heatmap <-function(data = NULL,
                             names.variables = NULL,
                             names.rows = NULL,
                             dendrogram = NULL,
                             clusters = NULL
){# Step One: Data preparation
  data_df<-as.data.frame(data[names.variables])
  row.names(data_df)<-data[[names.rows]]
  
  #dendogram data
  dendro <- dendro_data(dendrogram)
  # add cluster info to dendro object
  dendro[["labels"]][["clusters"]]<-NA
  for(i in names(clusters)){
    dendro[["labels"]][dendro[["labels"]][["label"]]==i,"clusters"] <- clusters[i]
  }
  
  # define cluster dimensions 
  clusters_dim <- NULL
  
  for(i in levels(as.factor(dendro[["labels"]][["clusters"]]))){
    temp <- tibble(cluster = i, n = nrow(dendro[["labels"]][dendro[["labels"]][["clusters"]]==i,]))
    temp$x_min <- min(dendro[["labels"]][dendro[["labels"]][["clusters"]]==i,"x"])
    temp$x_max <- max(dendro[["labels"]][dendro[["labels"]][["clusters"]]==i,"x"])
    if(i == "1"){
      temp$y_max <- as.numeric(levels(as.factor(dendro$segments$y)))[nrow(dendro[["labels"]]) - as.integer(as.numeric(i)/2)-1]
      temp$yend_max <- as.numeric(levels(as.factor(dendro$segments$yend)))[nrow(dendro[["labels"]]) - as.integer(as.numeric(i)/2)-1]
    }else{
      temp$y_max <- as.numeric(levels(as.factor(dendro$segments$y)))[nrow(dendro[["labels"]]) - as.integer(as.numeric(i)/2)-1]
      temp$yend_max <- as.numeric(levels(as.factor(dendro$segments$yend)))[nrow(dendro[["labels"]]) - as.integer(as.numeric(i)/2)-1]
    }
    clusters_dim<-rbind(clusters_dim, temp)
  }
  
  # create a function to plot the dendrogram in ggplot2 format
  ggdend <- function(df) {
    ggplot() +
      geom_segment(data = df, aes(x=x, y=y, xend=xend, yend=yend, color = cluster)) +
      labs(x = "", y = "") + theme_minimal() +
      theme(axis.text = element_blank(), axis.ticks = element_blank(),
            panel.grid = element_blank(),
            legend.position = "none")
  }
  
  # add clusters to segments
  dendro$segments$cluster <- NA
  for(i in clusters_dim[["cluster"]]){
    dendro[["segments"]][dendro$segments$x >= clusters_dim[clusters_dim[["cluster"]]==i,][["x_min"]] &
                           dendro$segments$xend >= clusters_dim[clusters_dim[["cluster"]]==i,][["x_min"]] &
                           dendro$segments$x <= clusters_dim[clusters_dim[["cluster"]]==i,][["x_max"]] &
                           dendro$segments$xend <= clusters_dim[clusters_dim[["cluster"]]==i,][["x_max"]] &
                           (dendro$segments$y <= clusters_dim[clusters_dim[["cluster"]]==i,][["y_max"]] |
                              dendro$segments$yend <= clusters_dim[clusters_dim[["cluster"]]==i,][["yend_max"]]),][["cluster"]]<-i
    
  }
  
  # x/y dendograms
  dendro_plot <- ggdend(dendro$segments) + coord_flip()
  # heatmap
  col.ord <- order.dendrogram(dendrogram)
  data_matrix <- as.matrix(scale(data_df[col.ord,]))
  
  data_names <- attr(data_matrix, "dimnames")
  data_df2 <- as.data.frame(data_matrix)
  colnames(data_df2) <- data_names[[2]]
  data_df2$id <- data_names[[1]]
  data_df2$id <- with(data_df2, factor(id, levels=id, ordered=TRUE))
  melt_df <- reshape2::melt(data_df2, id.vars="id")
  Heatmap_plot <- ggplot(melt_df, aes(x = variable, y = id)) +
    geom_tile(aes(fill = value))+
    theme_bw()+
    scale_fill_gradient(low="white", high="red")+
    theme(axis.text.y= element_text(size=8))+
    theme(axis.text.x= element_text(angle = 0, 
                                    vjust = 1, 
                                    size = 10.5, 
                                    hjust = .5))+
    labs(y="eaid",
         x="Metrics",
         fill = "Scale Value")+
    theme(plot.title = element_text(hjust=0.5, size=12))+
    theme(axis.title.x = element_text(size=10))+
    theme(axis.title.y = element_text(size=10))+
    theme(legend.title = element_text(size=10))+
    theme(legend.position="bottom")
  
  Clust.Heatmap.Plots <- list(dendro = dendro,
                              clusters_dim = clusters_dim,
                              dendro_plot = dendro_plot,
                              Heatmap_plot = Heatmap_plot)
  
  return(Clust.Heatmap.Plots)
}