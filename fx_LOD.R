

fx_lod <- function(data = NULL,
                  concentration = NULL,
                  result = NULL,
                  lower = TRUE,
                  upper = TRUE,
                  precision = 0.1,
                  from = 0,
                  to = 20000,
                  by = 1,
                  xlab = "Concentration",
                  ylab = "Proportion of positive results"){
  data0 <- data[c(concentration, result)]
  names(data0) <- c("concentration", "result")
  
  data0$result<-as.factor(data0$result)
  
  model = glm(result ~ concentration, family = binomial, data = data0)
  
  data1 <- tibble(concentration = seq(from = from, to = to, by = by))
  pred.data <- predict(model, newdata = data1, se.fit = TRUE)
  
  data1$fit <- pred.data$fit
  data1$se.fit <- pred.data$se.fit
  
  data1 %<>% mutate(prob = exp(fit)/(1 + exp(fit)))
  data1 %<>% mutate(upper = exp(fit + 1.96*se.fit)/(1 + exp(fit + 1.96*se.fit)))
  data1 %<>% mutate(lower = exp(fit - 1.96*se.fit)/(1 + exp(fit - 1.96*se.fit)))
  
  
  coef <- data.frame(t(coefficients(model)))
  names(coef)<-c("Intercept", "Concentration")
  
  # LOD determination
  
  lod <- (log(0.95/(1 - 0.95)) - coef$Intercept)/coef$Concentration
  pred <- predict(model, newdata = data.frame(concentration = lod), se.fit = TRUE)
  
  data2 <- tibble(lod = lod, lower = as.double(NA), upper = as.double(NA))
  
  # LOD lower limit determination
  if (lower){
    conc = lod
    pred.sup.n <- predict(model, newdata = data.frame(concentration=conc), se.fit = TRUE)
    pro.sup.n <- exp(pred.sup.n$fit + 1.96*pred.sup.n$se.fit)/(1 + exp(pred.sup.n$fit + 1.96*pred.sup.n$se.fit))
    while (pro.sup.n > 0.95) {
      pred.sup.n <- predict(model, newdata = data.frame(concentration = conc), se.fit = TRUE)
      pro.sup.n <- exp(pred.sup.n$fit + 1.96*pred.sup.n$se.fit)/(1 + exp(pred.sup.n$fit + 1.96*pred.sup.n$se.fit))
      conc = conc - precision
    }
    conc.inf.n <- (log(0.95/(1 - 0.95)) - 1.96*pred.sup.n$se.fit - coef$Intercept)/coef$Concentration  
    data2$lower <- conc.inf.n  
    }
  
  
  # LOD upper limit determination
  if (upper){
    conc = lod
    pred.inf.n <- predict(model, newdata = data.frame(concentration = conc), se.fit = TRUE)
    pro.inf.n <- exp(pred.inf.n$fit - 1.96*pred.inf.n$se.fit)/(1 + exp(pred.inf.n$fit - 1.96*pred.inf.n$se.fit))
    while (pro.inf.n < 0.95) {
      pred.inf.n <- predict(model, newdata = data.frame(concentration=conc), se.fit=TRUE)
      pro.inf.n <- exp(pred.inf.n$fit - 1.96*pred.inf.n$se.fit)/(1 + exp(pred.inf.n$fit - 1.96*pred.inf.n$se.fit))
      conc = conc + precision
    }
    conc.sup.n <- (log(0.95/(1 - 0.95)) + 1.96*pred.inf.n$se.fit - coef$Intercept)/coef$Concentration  
    data2$upper <- conc.sup.n
    }
  
  
  lod.plot <- data1 %>% ggplot(aes(x = concentration, y = prob)) +
    geom_line()+
    geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2) +
    geom_vline(xintercept = lod, color = "firebrick3", size = 1) +
    geom_point(x = lod, y = 0.95, color = "firebrick3") +
    geom_text(aes(x = lod*1.065, y = 0.95*1.05), label = round(lod, 2), color = "firebrick3", size = 3.2 )+
    theme_bw() +
    labs(y = ylab,
         x = xlab) +
    theme(plot.title = element_text(hjust = 0, size = 11),
          axis.title = element_text(size = 10),
          axis.text = element_text(size = 10)) +
    scale_x_continuous(expand = c(0, 0.02*max(data1$concentration)), breaks = seq(from, to, 2000)) +
    scale_y_continuous(expand = c(0, 0.02*max(data1$upper)), breaks = seq(0, 1, 0.1))
  
  if (lower){
    lod.plot <- lod.plot + 
      geom_vline(xintercept = conc.inf.n, color = "gray50", size = 1) +
      geom_point(x = conc.inf.n, y = 0.95, color = "gray40") +
      geom_text(aes(x = conc.inf.n*0.92, y = 0.95), label = round(conc.inf.n, 2), color = "gray40", size = 3.2)
  }
  
  if (upper){
    lod.plot <- lod.plot + 
      geom_vline(xintercept = conc.sup.n, color = "gray50", size = 1) + 
      geom_point(x = conc.sup.n, y = 0.95, color = "gray40") +
      geom_text(aes(x = conc.sup.n*1.06, y = 0.95), label = round(conc.sup.n, 2), color = "gray40", size = 3.2)
  }
    
  
  lod.analysis <- list(lod = data2, glm.lod = data1, lod.plot = lod.plot)
  
  return(lod.analysis)
}
