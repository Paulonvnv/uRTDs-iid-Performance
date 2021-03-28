# function to perform the heirarchical clustering method ----

fx_HClust.Analysis <- function(data = NULL,
                               names.variables = NULL,
                               names.rows = NULL,
                               distance = "auto",
                               method = "auto",
                               NbC.index ="alllong"){
  
  # Step One: Data preparation
  data_df<-as.data.frame(data[names.variables])
  row.names(data_df)<-data[[names.rows]]
  data_matrix <- as.matrix(scale(data_df))
  
  # step two: Measure pairwise distances
  # Define the method: auto, euclidean, manhattan,
  # canberra, binary, minkowski
  
  # Step three: define the clustering algorithm (two possibilities)
  ## Agglomrative clustering (hclust)
  # and the method (hc.method) to compute dissimilarities
  #between clusters:
  #hclust
  # "ward.D", "ward.D2", "single",
  #"complete", "average" (= UPGMA),
  #"mcquitty" (= WPGMA),
  #"median" (= WPGMC) or "centroid" (= UPGMC).
  
  # Step four: define the number of clusters
  
  if(distance == "auto" & method == "auto"){
    d <- c("euclidean", "manhattan",
           "canberra", "binary", "minkowski")
    m <- c("ward.D", "ward.D2", "single",
           "complete", "average",
           "mcquitty")
    ac <- NULL
    for (i in d){
      temp1 <- get_dist(data_matrix, method = i)
      for(j in m){
        temp2 <- hclust(temp1, method = j)
        temp3 <- tibble(distance = i, method = j, ac = coef.hclust(temp2))
        ac <- rbind(ac, temp3)
      }
    }
    dist.matrix <- get_dist(data_matrix, method = ac[which.max(ac$ac),][["distance"]])
    h.clust <- hclust(dist.matrix, method = ac[which.max(ac$ac),][["method"]])
    Nb.Clust<-NbClust(data_matrix,
                      distance = ac[which.max(ac$ac),][["distance"]],
                      method = ac[which.max(ac$ac),][["method"]],
                      index = NbC.index)
    rm(list = c("i", "j", "temp1", "temp2", "temp3"))
  }else if(distance == "auto" & method != "auto"){
    d <- c("euclidean", "manhattan",
           "canberra", "binary", "minkowski")
    ac <- NULL
    for (i in d){
      temp1<-get_dist(data_matrix, method = i)
      temp2 <- hclust(temp1, method = method)
      temp3 <- tibble(distance = i, method = method, ac = coef.hclust(temp2))
      ac <- rbind(ac, temp3)
    }
    dist.matrix <- get_dist(data_matrix, method = ac[which.max(ac$ac),][["distance"]])
    h.clust <- hclust(dist.matrix, method = method)
    Nb.Clust<-NbClust(data_matrix,
                      distance = ac[which.max(ac$ac),][["distance"]],
                      method = method,
                      index = NbC.index)
    rm(list = c("i", "temp1", "temp2", "temp3"))
  }else if(distance != "auto" & method == "auto"){
    m <- c("ward.D", "ward.D2", "single",
           "complete", "average",
           "mcquitty")
    ac <- NULL
    temp1<-get_dist(data_matrix, method = distance)
    for(j in m){
      temp2 <- hclust(temp1, method = j)
      temp3 <- tibble(distance = distance, method = j, ac = coef.hclust(temp2))
      ac <- rbind(ac, temp3)
    }
    dist.matrix <- get_dist(data_matrix, method = distance)
    h.clust <- hclust(dist.matrix, method = ac[which.max(ac$ac),][["method"]])
    Nb.Clust<-NbClust(data_matrix,
                      distance = distance,
                      method = ac[which.max(ac$ac),][["method"]],
                      index = NbC.index)
    rm(list = c("j", "temp1", "temp2", "temp3"))
  }else{
    dist.matrix <- get_dist(data_matrix, method = distance)
    h.clust <- hclust(dist.matrix, method = method)
    ac <- tibble(distance = distance, method = method, ac = coef.hclust(h.clust))
    Nb.Clust<-NbClust(data_matrix,
                      distance = distance,
                      method = method,
                      index = NbC.index)
  }
  
  HClust.Analysis <- list(dist.matrix = dist.matrix,
                          h.clust = h.clust,
                          ac = ac,
                          Nb.Clust = Nb.Clust)
  return(HClust.Analysis)
}
