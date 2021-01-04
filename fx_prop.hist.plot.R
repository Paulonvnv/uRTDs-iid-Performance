

fx_prop.hist.plot <- function(data = NULL,
                              bins = 60,
                              quant.var = NULL,
                              qual.var = NULL,
                              group.var = NULL,
                              log.scale = TRUE,
                              colors = c("dodgerblue4", "firebrick4"),
                              y.title = "Proportion of samples",
                              x.title = "Concentration in log10 scale",
                              legend.title = "Result",
                              legend.position = "right"){
  if (log.scale){
    data[data[[quant.var]]!=0,][[quant.var]] <- log10(data[data[[quant.var]]!=0,][[quant.var]])
    data1 <- tibble(lower = seq(0,
                                max(data[[quant.var]], na.rm = T) - max(data[[quant.var]], na.rm = T)/bins,
                                max(data[[quant.var]], na.rm = T)/bins),
                    upper = seq(max(data[[quant.var]], na.rm = T)/bins,
                                max(data[[quant.var]], na.rm = T),
                                max(data[[quant.var]], na.rm = T)/bins))
  }else{
    data1 <- tibble(lower = seq(0,
                                max(data[[quant.var]], na.rm = T) - max(data[[quant.var]], na.rm = T)/bins,
                                max(data[[quant.var]], na.rm = T)/bins),
                    upper = seq(max(data[[quant.var]], na.rm = T)/bins,
                                max(data[[quant.var]], na.rm = T),
                                max(data[[quant.var]], na.rm = T)/bins))    
  }
  
  
  
  
  data2 <- NULL
  for (i in levels(as.factor(as.character(data[[group.var]])))){
    for (j in levels(as.factor(as.character(data[[qual.var]])))){
      for (n in 1:bins){
        temp <- tibble(group = i,
                       result = j,
                       bin  = n,
                       lower = data1[n,][["lower"]],
                       upper = data1[n,][["upper"]])
        temp$prop <- nrow(data[data[[group.var]] == i&
                                 data[[qual.var]] == j&
                                 (data[[quant.var]]>=data1[n,][["lower"]]&data[[quant.var]]<=data1[n,][["upper"]]),])/
          nrow(data[data[[group.var]] == i&
                      (data[[quant.var]]>=data1[n,][["lower"]]&data[[quant.var]]<=data1[n,][["upper"]]),])
        
        data2 <- rbind(data2, temp)
      }
    }
  }
  
  
  prop.plot <- data2 %>% 
    ggplot(aes(x = upper, y = prop, fill = result, color = result)) +
    geom_bar(stat = "identity") + 
    facet_wrap(~group, nrow = 3) + 
    theme_classic()+
    scale_color_manual(values = colors)+
    scale_fill_manual(values = colors)+
    labs(y = y.title,
         x = x.title,
         color = legend.title,
         fill = legend.title)+
    theme(axis.title.x = element_text(size = 10)) +
    theme(axis.title.y = element_text(size = 10)) +
    theme(axis.text = element_text(size = 10), legend.position = legend.position)
  
  if (log.scale){
    prop.plot <- prop.plot +
      scale_x_continuous(breaks = 0:floor(max(data2$upper)), labels = 10^(0:floor(max(data2$upper))))
  }
  
  prop.hist <- list(prop.hist = data2, prop.plot = prop.plot)
  
  return(prop.hist)
   
}






