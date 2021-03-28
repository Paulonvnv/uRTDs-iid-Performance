# function to calculate the roc curve, the cut point and a summary table with cutpointer ----
fx_roc_analysis <- function(data = NULL,
                            predictor.vars = NULL,
                            outcome.vars = NULL,
                            method = maximize_metric, # method to determine the cutpoints
                            metric = youden, # function for computing a metric when maximize_metric or minimize_metric is used
                            pos_class = NULL, # value indicatin the positive class
                            neg_class = NULL,
                            direction = NULL,
                            boot_runs = 1000, # number of runs to assess the variability
                            boot_stratify = FALSE,
                            use_midpoints = FALSE,
                            break_ties = mean,
                            na.rm = FALSE,
                            allowParallel = TRUE,
                            silent = FALSE,
                            tol_metric = 1e-06,
                            ...
){
  
  # first run
  
  i = predictor.vars[1]
  j = outcome.vars[1]
  
  data1 <- data[,c(i,j)]
  
  names(data1) <- c("predictor.var", "outcome.var")
  
  roc_temp <- cutpointr(data = data1, # Data
                        x = predictor.var, # predictor
                        class = outcome.var, # class membership
                        method = method, # method to determine the cutpoints
                        metric = metric, # function for computing a metric when maximize_metric or minimize_metric is used
                        pos_class = pos_class, # value indicatin the positive class
                        neg_class = neg_class,
                        direction = direction,
                        boot_runs = boot_runs, # number of runs to assess the variability
                        boot_stratify = boot_stratify,
                        use_midpoints = use_midpoints,
                        break_ties = break_ties,
                        na.rm = na.rm,
                        allowParallel = allowParallel,
                        silent = silent,
                        tol_metric = tol_metric
  )
  rm(data1)
  
  sum_temp <- summary(roc_temp)
  
  tp <- sum_temp[["confusion_matrix"]][[1]][["tp"]]
  fn <- sum_temp[["confusion_matrix"]][[1]][["fn"]]
  fp <- sum_temp[["confusion_matrix"]][[1]][["fp"]]
  tn <- sum_temp[["confusion_matrix"]][[1]][["tn"]]
  
  temp <- rbind(tibble(outcome = j, predictor = i, metric = "sensitivity", fx_sensitivity(tp = tp, fn = fn)),
                tibble(outcome = j, predictor = i, metric = "specificity", fx_specificity(tn = tn, fp = fp)),
                tibble(outcome = j, predictor = i, metric = "ppv", fx_ppv(tp = tp, fp = fp)),
                tibble(outcome = j, predictor = i, metric = "npv", fx_npv(tn = tn, fn = fn)),
                tibble(outcome = j, predictor = i, metric = "accuracy", fx_accuracy(tp = tp, fp = fp, tn = tn, fn = fn)),
                tibble(outcome = j, predictor = i, metric = "auc", fx_ci.auc(roc_temp))
  )
  rm(list = c("tp","fp","tn","fn"))
  
  metrics <- NULL
  metrics <- rbind(metrics, temp)
  rm(temp)
  
  roc_temp$metrics<-list(metrics)
  rm(metrics)
  
  roc_temp[["outcome"]] <- j
  roc_temp[["predictor"]] <- i
  
  roc_temp[["data"]][[1]][["outcome"]] <- j
  roc_temp[["data"]][[1]][["predictor"]] <- i
  
  roc_temp[["roc_curve"]][[1]][["outcome"]] <- j 
  roc_temp[["roc_curve"]][[1]][["predictor"]] <- i 
  
  names(roc_temp[["data"]])<-paste(j,i,sep = "_")
  names(roc_temp[["roc_curve"]])<-paste(j,i,sep = "_")
  names(roc_temp[["boot"]])<-paste(j,i,sep = "_")
  names(roc_temp[["metrics"]])<-paste(j,i,sep = "_")
  
  roc_analysis<-roc_temp
  rm(list = c("roc_temp","i","j"))
  
  if (length(predictor.vars) > 1){
    for (j in outcome.vars){
      for (i in predictor.vars[-1]){
        
        data1 <- data[,c(i,j)]
        
        names(data1) <- c("predictor.var", "outcome.var")
        
        roc_temp <- cutpointr(data = data1, # Data
                              x = predictor.var, # predictor
                              class = outcome.var, # class membership
                              method = method, # method to determine the cutpoints
                              metric = metric, # function for computing a metric when maximize_metric or minimize_metric is used
                              pos_class = pos_class, # value indicatin the positive class
                              neg_class = neg_class,
                              direction = direction,
                              boot_runs = boot_runs, # number of runs to assess the variability
                              boot_stratify = boot_stratify,
                              use_midpoints = use_midpoints,
                              break_ties = break_ties,
                              na.rm = na.rm,
                              allowParallel = allowParallel,
                              silent = silent,
                              tol_metric = tol_metric)
        rm(data1)
        
        sum_temp <- summary(roc_temp)
        
        tp <- sum_temp[["confusion_matrix"]][[1]][["tp"]]
        fn <- sum_temp[["confusion_matrix"]][[1]][["fn"]]
        fp <- sum_temp[["confusion_matrix"]][[1]][["fp"]]
        tn <- sum_temp[["confusion_matrix"]][[1]][["tn"]]
        
        temp <- rbind(tibble(outcome = j, predictor = i, metric = "sensitivity", fx_sensitivity(tp = tp, fn = fn)),
                      tibble(outcome = j, predictor = i, metric = "specificity", fx_specificity(tn = tn, fp = fp)),
                      tibble(outcome = j, predictor = i, metric = "ppv", fx_ppv(tp = tp, fp = fp)),
                      tibble(outcome = j, predictor = i, metric = "npv", fx_npv(tn = tn, fn = fn)),
                      tibble(outcome = j, predictor = i, metric = "accuracy", fx_accuracy(tp = tp, fp = fp, tn = tn, fn = fn)),
                      tibble(outcome = j, predictor = i, metric = "auc", fx_ci.auc(roc_temp))
        )
        rm(list = c("tp","fp","tn","fn"))
        
        metrics <- NULL
        metrics <- rbind(metrics, temp)
        rm(temp)
        
        roc_temp$metrics<-list(metrics)
        rm(metrics)
        
        roc_temp[["outcome"]] <- j
        roc_temp[["predictor"]] <- i
        
        roc_temp[["data"]][[1]][["outcome"]] <- j
        roc_temp[["data"]][[1]][["predictor"]] <- i
        
        roc_temp[["roc_curve"]][[1]][["outcome"]] <- j 
        roc_temp[["roc_curve"]][[1]][["predictor"]] <- i 
        
        names(roc_temp[["data"]])<-paste(j,i,sep = "_")
        names(roc_temp[["roc_curve"]])<-paste(j,i,sep = "_")
        names(roc_temp[["boot"]])<-paste(j,i,sep = "_")
        names(roc_temp[["metrics"]])<-paste(j,i,sep = "_")
        
        roc_analysis <- rbind(roc_analysis,roc_temp)
        
      }
    }
    rm(list = c("roc_temp","i","j"))
  }
  
  metrics_comb <- NULL
  
  for (i in names(roc_analysis[["metrics"]])){
    temp<-roc_analysis[["metrics"]][[i]]
    metrics_comb<-rbind(metrics_comb,temp)
    rm(temp)
  }
  
  data_comb <- NULL
  
  for (i in names(roc_analysis[["data"]])){
    temp<-roc_analysis[["data"]][[i]]
    data_comb<-rbind(data_comb,temp)
    rm(temp)
  }
  
  
  rocc_comb <- NULL
  
  for (i in names(roc_analysis[["roc_curve"]])){
    temp<-roc_analysis[["roc_curve"]][[i]]
    rocc_comb<-rbind(rocc_comb,temp)
    rm(temp)
  }
  
  sens_points <- c(roc_analysis[["sensitivity"]], rep(NA, times = (nrow(rocc_comb) - nrow(roc_analysis))))
  spec_points <- c((1 - roc_analysis[["specificity"]]), rep(NA, times = (nrow(rocc_comb) - nrow(roc_analysis))))
  
  roc_plot <- ggplot(rocc_comb, aes(x=fpr, y=tpr, group=predictor)) +
    geom_line(aes(linetype = predictor, color=predictor), size = 1.25) +
    geom_point(aes(x = spec_points,
                   y = sens_points), size = 2) +
    labs(title = "ROC Curve",
         y = "Sensitivity",
         x = "1 - Specificity") +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5, size = 14)) +
    theme(axis.title.x = element_text(size = 12)) +
    theme(axis.title.y = element_text(size = 12)) +
    theme(axis.text = element_text(size = 11)) +
    theme(legend.position="bottom")
  
  metrics_plot <- ggplot(metrics_comb, aes(x = predictor, y = value, ymin = lower, ymax = upper)) +
    geom_pointrange() + 
    geom_hline(yintercept = 0.8, lty = 2) +  # add a dotted line at x=1 after flip
    coord_flip() +  # flip coordinates (puts labels on y axis)
    xlab("Predictor") + ylab("Mean (95% CI)") +
    theme_bw() +  # use a white background
    facet_wrap( ~ metric, ncol = 3)
  
  roc_analysis <- tibble(roc_analysis = list(roc_analysis),
                         data = list(data_comb),
                         roc_curve = list(rocc_comb),
                         metrics = list(metrics_comb),
                         roc_plot = list(roc_plot),
                         metrics_plot = list(metrics_plot)
  )
  
  return(roc_analysis)
  
}