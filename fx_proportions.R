fx_proportions <- function(data = NULL,
                           strata = NULL,
                           var = NULL,
                           obs = NULL,
                           total = TRUE){
  data0 <- data[c(strata,var)]
  names(data0) <- c("strata", "var")
  
  
  proportions<-data0 %>% as_tibble() %>%
    group_by(strata) %>% summarise(count(var), n = length(strata))%>%
    filter(x == obs) %>% select(-x) %>% mutate(prop = binconf(freq,
                                                              n,
                                                              alpha = 0.05,
                                                              method = "wilson")[1],
                                               lower = binconf(freq,
                                                               n,
                                                               alpha = 0.05,
                                                               method = "wilson")[2],
                                               upper = binconf(freq,
                                                               n,
                                                               alpha = 0.05,
                                                               method = "wilson")[3])
  if (total == T){
    propTot <- tibble("strata" = "Total",
                      data0 %>% as_tibble() %>% summarise(count(var), n = length(strata))%>%
                        filter(x == obs) %>% select(-x) %>% mutate(prop = binconf(freq,
                                                                                  n,
                                                                                  alpha = 0.05,
                                                                                  method = "wilson")[1],
                                                                   lower = binconf(freq,
                                                                                   n,
                                                                                   alpha = 0.05,
                                                                                   method = "wilson")[2],
                                                                   upper = binconf(freq,
                                                                                   n,
                                                                                   alpha = 0.05,
                                                                                   method = "wilson")[3]))
    proportions <- rbind(proportions,propTot)
    return(proportions)
  }else{
    return(proportions)    
  }
}
