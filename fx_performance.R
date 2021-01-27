

fx_performance <- function (data = NULL,
                          reference = NULL,
                          test = NULL,
                          pos = NULL,
                          neg = NULL){
  
  data3 <- NULL
  data5 <- NULL
  
  for (i in test){
    data1 <- data[,c(reference,i)]
    
    names(data1)<-c("reference", "test")
    
    data1 %<>% mutate(kind_result = case_when(
      reference == pos & test == pos ~ "tp",
      reference == neg & test == pos ~ "fp",
      reference == pos & test == neg ~ "fn",
      reference == neg & test == neg ~ "tn",
    ))
    
    data2 <- data1 %>% group_by(kind_result) %>% summarise(count = n())
    
    data2[["test"]] <- i
    
    sens <- binconf(data2[data2[["kind_result"]]=="tp",][["count"]],
                     (data2[data2[["kind_result"]]=="tp",][["count"]]+data2[data2[["kind_result"]]=="fn",][["count"]]),
                     alpha = 0.05,
                     method = "wilson")
    spec <- binconf(data2[data2[["kind_result"]]=="tn",][["count"]],
                    (data2[data2[["kind_result"]]=="tn",][["count"]]+data2[data2[["kind_result"]]=="fp",][["count"]]),
                    alpha = 0.05,
                    method = "wilson")
    ppv <- binconf(data2[data2[["kind_result"]]=="tp",][["count"]],
                   (data2[data2[["kind_result"]]=="tp",][["count"]]+data2[data2[["kind_result"]]=="fp",][["count"]]),
                   alpha = 0.05,
                   method = "wilson")
    npv <- binconf(data2[data2[["kind_result"]]=="tn",][["count"]],
                   (data2[data2[["kind_result"]]=="tn",][["count"]]+data2[data2[["kind_result"]]=="fn",][["count"]]),
                   alpha = 0.05,
                   method = "wilson")
    accu <- binconf(data2[data2[["kind_result"]]=="tp",][["count"]] + data2[data2[["kind_result"]]=="tn",][["count"]],
                    sum(data2[["count"]]),
                    alpha = 0.05,
                    method = "wilson")
    
    data3 <- rbind(data3, data2)
    
    data4 <- rbind(tibble(test=i, metric = "sens",  value = sens[1]*100, lower = sens[2]*100, upper = sens[3]*100),
                   tibble(test=i, metric = "spec",  value = spec[1]*100, lower = spec[2]*100, upper = spec[3]*100),
                   tibble(test=i, metric = "ppv",  value = ppv[1]*100, lower = ppv[2]*100, upper = ppv[3]*100),
                   tibble(test=i, metric = "npv",  value = npv[1]*100, lower = npv[2]*100, upper = npv[3]*100),
                   tibble(test=i, metric = "accu",  value = accu[1]*100, lower = accu[2]*100, upper = accu[3]*100))
    
    data5 <- rbind(data5, data4)
    }
  
  data6 <- list(conf_table = data3, metrics = data5)
  return(data6)
  
}





