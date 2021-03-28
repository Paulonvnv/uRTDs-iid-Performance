fx_pROC_analysis <- function(response = NULL,
                             predictors = NULL,
                             data = NULL,
                             ci.method = "bootstrap",
                             boot.n = 10000,
                             boot.stratified = T,
                             best.method = "youden",
                             conf.level = .95,
                             x = "best",
                             input = NULL
){
  # Constructing ROC curve
  roc_analysis <- NULL
  temp2 <- NULL
  
  for(i in response){
    temp1 <- pROC::roc(as.formula(paste0(i, sep = "~", paste(predictors, collapse = "+"))),
                       data = data,
                       ci = T,
                       ci.method = ci.method,
                       boot.stratified = boot.stratified,
                       boot.n = boot.n)
    temp2[i] <- list(temp1)
    roc_analysis <- c(roc_analysis, temp2)
  }
  
  rm(list=c("temp1", "temp2"))
  # Calculating metrics pluss ci: Sens, Spec, ppv, npv, Auccuracy, threshold 
  
  fx_sensitivity <- function(tp, fn){
    sens <- as_tibble(binconf(tp,
                              (tp + fn),
                              alpha = 0.05,
                              method = "wilson"))
    colnames(sens) <- c("value", "lower", "upper")
    return(sens)
  }
  
  
  fx_specificity <- function(tn, fp){
    spec <- as_tibble(binconf(tn,
                              (tn + fp),
                              alpha = 0.05,
                              method = "wilson"))
    colnames(spec) <- c("value", "lower", "upper")
    return(spec)
  }
  
  fx_ppv <- function(tp, fp){
    ppv <- as_tibble(binconf(tp,
                             (tp + fp),
                             alpha = 0.05,
                             method = "wilson"))
    colnames(ppv) <- c("value", "lower", "upper")
    return(ppv)
  }
  
  fx_npv <- function(tn, fn){
    npv <- as_tibble(binconf(tn,
                             (tn + fn),
                             alpha = 0.05,
                             method = "wilson"))
    colnames(npv) <- c("value", "lower", "upper")
    return(npv)
  }
  
  fx_accuracy <- function(tp, fp, tn, fn){
    accuracy <- as_tibble(binconf((tp + tn),
                                  (tp + fp + tn + fn),
                                  alpha = 0.05,
                                  method = "wilson"))
    colnames(accuracy) <- c("value", "lower", "upper")
    return(accuracy)
  }
  
  
  fx_ci.auc <- function(boot.roc){
    auc <- tibble(value = boot.roc$AUC)
    auc$lower <- boot_ci(boot.roc, AUC, in_bag = TRUE, alpha = 0.05)[1,][["values"]]
    auc$upper <- boot_ci(boot.roc, AUC, in_bag = TRUE, alpha = 0.05)[2,][["values"]]
    return(auc)
  }
  
  
  metrics <-NULL
  for(j in levels(as.factor(names(roc_analysis)))){
    for(i in names(roc_analysis[[j]])){
      temp1 <- coords(roc_analysis[[j]][[i]],
                      x = x,
                      input = input,
                      best.method = best.method,
                      best.policy = "random",
                      ret = c("tp", "tn", "fp", "fn"),
                      transpose = F)
      
      tp<-temp1[["tp"]]
      fp<-temp1[["fp"]]
      tn<-temp1[["tn"]]
      fn<-temp1[["fn"]]
      
      temp2 <- ci.coords(roc_analysis[[j]][[i]],
                         x = "best",
                         ret = "threshold",
                         best.method = best.method,
                         best.weights = c(1, 0.5),
                         best.policy = "random",
                         boot.n = boot.n,
                         conf.level = conf.level,
                         transpose = T)
      
      temp3 <- rbind(tibble(outcome = j, predictor = i, metric = "sensitivity", fx_sensitivity(tp = tp, fn = fn)),
                     tibble(outcome = j, predictor = i, metric = "specificity", fx_specificity(tn = tn, fp = fp)),
                     tibble(outcome = j, predictor = i, metric = "ppv", fx_ppv(tp = tp, fp = fp)),
                     tibble(outcome = j, predictor = i, metric = "npv", fx_npv(tn = tn, fn = fn)),
                     tibble(outcome = j, predictor = i, metric = "accuracy", fx_accuracy(tp = tp, fp = fp, tn = tn, fn = fn)),
                     tibble(outcome = j, predictor = i, metric = "auc", value = roc_analysis[[j]][[i]][["ci"]][[2]],
                            lower = roc_analysis[[j]][[i]][["ci"]][[1]],
                            upper = roc_analysis[[j]][[i]][["ci"]][[3]]),
                     tibble(outcome = j, predictor = i, metric = "threshold", value = temp2[["threshold"]][[2]],
                            lower = temp2[["threshold"]][[1]],
                            upper = temp2[["threshold"]][[3]])
      )  
      
      metrics<-rbind(metrics, temp3)
      rm(list = c("temp1", "temp2", "temp3", "tp", "fp", "tn", "fn"))
    }}
  
  
  
  auc.test <- tibble(response = as.character(NULL),
                     comparison = as.character(NULL),
                     roc1 = as.character(NULL),
                     roc2 = as.character(NULL),
                     auc1 = as.double(NULL),
                     auc2 = as.double(NULL),
                     difference = as.double(NULL),
                     p.value = as.double(NULL),
                     power = as.double(NULL))
  
  for (h in levels(as.factor(names(roc_analysis)))){
    for(i in names(roc_analysis[[h]])){
      for(j in names(roc_analysis[[h]])){
        if (i != j & !(paste(j, i, sep = "-") %in% auc.test[auc.test[["response"]] == h,][["comparison"]])){
          temp <- tibble(response = h, comparison = paste(i, j, sep = "-"), roc1 = i, roc2 = j, auc1 = pROC::auc(roc_analysis[[h]][[i]]), auc2 = pROC::auc(roc_analysis[[h]][[j]]))
          temp$difference <- pROC::auc(roc_analysis[[h]][[i]]) - pROC::auc(roc_analysis[[h]][[j]])
          temp1 <- roc.test(roc_analysis[[h]][[i]],roc_analysis[[h]][[j]], method = ci.method)
          temp$p.value <- temp1$p.value
          temp2 <- power.roc.test(roc_analysis[[h]][[i]], roc_analysis[[h]][[j]], method = "delong")
          temp$power <- temp2$power
          auc.test <- rbind(auc.test, temp) 
        }
      }
    }  
  }
  
  
  rm(list = c("i", "j", "temp", "temp1", "temp2"))  
  
  
  roc_curves <- NULL
  
  for (h in levels(as.factor(names(roc_analysis)))){
    for (i in names(roc_analysis[[h]])){
      temp <- tibble(response = h,
                     predictor = i,
                     sensitivities = roc_analysis[[h]][[i]][["sensitivities"]],
                     specificities = roc_analysis[[h]][[i]][["specificities"]],
                     thresholds = roc_analysis[[h]][[i]][["thresholds"]])
      roc_curves <- rbind(roc_curves,temp)
      rm(temp)
    }  
  }
  
  roc_analysis <- list(roc_analysis = roc_analysis,
                       metrics = metrics,
                       roc_curves = roc_curves,
                       auc.test = auc.test)
  
  return(roc_analysis)
  
}