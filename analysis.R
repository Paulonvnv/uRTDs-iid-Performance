# 1. Data Managment----
# 1.1. Charge packages and functions ----

setwd("C:/Users/paulo/OneDrive/Github/uRDTs-iid-performance/")

source("Required_packages.R")
source("fx_LOD.R")
source("fx_prop.hist.plot.R")
source("fx_performance.R")
source("fx_hmean.R")
source("fx_hclust_analysis.R")
source("fx_hclustheatmap.R")
source("fx_pROC_analysis.R")

# 1.2. Upload iid and eaid data with read.dta13----

iid_wdf<-read.dta13("Data/uRDT_analysis_BW_4146.dta", nonint.factors= T) # import uRDT data base, nonint.factors = T allows to import variable labels (some variables were stored as floats in stata)

# 1.3. tibble data.frame by iid and eaid---- 

# 1.3.1. iid_df: ----
# long format tibble data.frame for iid`s`, just some variables will be used in order to facilitate its processing
iid_wdf <- iid_wdf[c("iid", "hhid", "eaid",
                            "age", "gender",
                            "body_temp", "fever48",
                            "hrp2_posneg", "hrp2_pg_ml",
                            "qPCRposneg", "pden_mean",
                            "pf_result", "pm_result",
                            "po_result", "pv_result",
                            "rdtresult", "rdtresult_posneg",
                            "hsrdtresult"
)]

# Rename variables
names(iid_wdf) <- c(names(iid_wdf[1:7]),"HRP2", "hrp2_pg_ml",
                    "qPCR", "pden_mean",
                    "pf_result", "pm_result",
                    "po_result", "pv_result",
                    "RDT_cat", "RDT",
                    "uRDT") # rename tests


iid_wdf <- as_tibble(iid_wdf) # convert to tibble

iid_wdf[iid_wdf[["hhid"]] == "XSHH0855" & iid_wdf[["eaid"]] == "10699004",][["hhid"]]<-"XSHH0856" # Change the name of a repeated hhid code

# removing NA's
for (i in c("pf_result", "pm_result", "po_result", "pv_result")){
  iid_wdf[[i]] <- as.character(iid_wdf[[i]])
  iid_wdf[is.na(iid_wdf[[i]]),][[i]] <- "No tested"
  iid_wdf[[i]] <- as.factor(iid_wdf[[i]])
  rm(i)
}

iid_wdf[is.na(iid_wdf$hrp2_pg_ml),][["hrp2_pg_ml"]]<-0
iid_wdf[is.na(iid_wdf$pden_mean),][["pden_mean"]]<-0

# 1.3.2. Collapse the qualitative PCR results ----

iid_wdf %<>% mutate(qual_diag = case_when(
  pf_result == "Positive" & pm_result == "Negative" ~ "P. falciparum",
  pf_result == "Positive" & pm_result == "Positive" ~ "Pf Mixed infection",
  pf_result == "Negative" & pm_result == "Positive" ~ "P. malariae",
  pf_result == "Negative" & pm_result == "Negative" ~ "Negative",
  pf_result == "No tested" & pm_result == "Negative" ~ "Negative",
  pf_result == "No tested" & pm_result == "No tested" ~ "No tested"
))


# 1.3.3. convert to long format----

test<-c("HRP2", "qPCR", "RDT", "uRDT", "qual_diag")

iid_ldf<- iid_wdf %>% pivot_longer(all_of(test),# define the variables that will be converted to long format
                                   names_to = "test",# assign a new variable named "test" where the categorical information is going to be stored
                                   values_to = "result")# assign a new variables named "result" where the result of each test is going to be stored

# 2. Descriptive Analysis----

# 2.1. Summary information of the sampling process----
stbl1_summ_sampling_iid <-tibble(
  # Number of eaids
  neaids = nlevels(as.factor(iid_wdf$eaid)),
  # Number of iids                  
  niids = nlevels(as.factor(iid_wdf$iid)),
  # Number of hhids                   
  nhhids = nlevels(as.factor(iid_wdf$hhid)),
  # Number of performed tests                  
  ntests = nlevels(as.factor(iid_ldf$test)),
  # List of al performed tests                   
  tests_list = paste0(paste(levels(as.factor(iid_ldf$test))[-length(levels(as.factor(iid_ldf$test)))], collapse= ", "), sep = ", and ", levels(as.factor(iid_ldf$test))[length(levels(as.factor(iid_ldf$test)))]))

# 2.2. Consistency among tests----  

stbl2_consistency<-tibble(  
  # Number of positive samples by any test
  npos_any4 = iid_wdf %>% filter(HRP2 == "Positive"|
                                   qPCR == "Positive"|
                                   RDT == "Positive"|
                                   uRDT == "Positive") %>% nrow(),
  
  # Number of all positive samples detected by HRP2                   
  npos_HRP2 = iid_ldf %>% filter(test == "HRP2", result == "Positive")%>%nrow(),
  
  # Number of all positive samples detected by qPCR
  npos_qPCR = iid_ldf %>% filter(test == "qPCR", result == "Positive")%>%nrow(),
  
  # Number of all positive samples detected by RDTs
  npos_RDT = iid_ldf %>% filter(test == "RDT", result == "Positive")%>%nrow(),
  
  # Number of all positive samples detected by uRDTs
  npos_uRDT = iid_ldf %>% filter(test == "uRDT", result == "Positive")%>%nrow(),
  
  
  # Number of only HRP2 positive samples
  only_HRP2 = iid_wdf %>% filter(HRP2 =="Positive",
                                 qPCR !="Positive",
                                 RDT !="Positive",
                                 uRDT !="Positive")%>%nrow(),
  
  # Number of only qPCR positive samples
  only_qPCR = iid_wdf %>% filter(HRP2 !="Positive",
                                 qPCR =="Positive",
                                 RDT !="Positive",
                                 uRDT !="Positive")%>%nrow(),
  
  # Number of only RDT positive samples
  only_RDT = iid_wdf %>% filter(HRP2 !="Positive",
                                qPCR !="Positive",
                                RDT =="Positive",
                                uRDT !="Positive")%>%nrow(),
  # Number of only uRDTs positive samples
  only_uRDT = iid_wdf %>% filter(HRP2 !="Positive",
                                 qPCR !="Positive",
                                 RDT !="Positive",
                                 uRDT =="Positive")%>%nrow(),
  
  # Qualitative qPCR of only qPCR and positive samples
  pos_qPCR_qualdiag = list(iid_wdf %>% filter(HRP2 !="Positive",
                                              qPCR =="Positive",
                                              RDT !="Positive",
                                              uRDT !="Positive") %>%
                             group_by(qual_diag) %>% summarise(nSamples= n())))

# 2.2.1.VennD representation----
positive_iids <- list(HRP2 = iid_ldf[iid_ldf[["test"]] == "HRP2" & iid_ldf[["result"]] == "Positive",][["iid"]],
                      qPCR = iid_ldf[iid_ldf[["test"]] == "qPCR" & iid_ldf[["result"]] == "Positive",][["iid"]],
                      RDT = iid_ldf[iid_ldf[["test"]] == "RDT" & iid_ldf[["result"]] == "Positive",][["iid"]],
                      uRDT = iid_ldf[iid_ldf[["test"]] == "uRDT" & iid_ldf[["result"]] == "Positive",][["iid"]])

iid_VennD_4tests <- ggVennDiagram(positive_iids, label = "count")


qPCR_vennD<- ggVennDiagram(list(qPCR = iid_ldf[iid_ldf[["test"]] == "qPCR" & iid_ldf[["result"]] == "Positive",][["iid"]],
                                only_qPCR = iid_wdf %>% filter(HRP2 !="Positive",
                                                               qPCR =="Positive",
                                                               RDT !="Positive",
                                                               uRDT !="Positive")%>%select(iid)%>%unlist,
                                
                                "P. falciparum" = iid_wdf %>% filter(pf_result == "Positive")%>%
                                  select(iid)%>%unlist,
                                "P malariae" = iid_wdf %>% filter(pm_result == "Positive")%>%
                                  select(iid)%>%unlist), label = "count")

plot1.consistency <-plot_grid(iid_VennD_4tests, qPCR_vennD, ncol = 2, labels = c("A", "B"), label_size = 11, label_fontface = "plain")

# 2.3. HRP2 concentration and parasite density by test ----  

summ_quant<-NULL
for (i in test[-5]){
  for (j in c("hrp2_pg_ml","pden_mean")){
    temp <- tibble(test = i,
                   measurement = j,
                   median = iid_wdf %>%
                     filter(.[i]=="Positive")%>%
                     select(as.name(j))%>%
                     unlist()%>%
                     median(),
                   q25 = iid_wdf %>%
                     filter(.[i]=="Positive")%>%
                     select(as.name(j))%>%
                     unlist()%>%
                     quantile(probs = .25),
                   q75 = iid_wdf %>%
                     filter(.[i]=="Positive")%>%
                     select(as.name(j))%>%
                     unlist()%>%
                     quantile(probs = .75))
    summ_quant<- rbind(summ_quant, temp)
  }
}

rm(list=c("i","j","temp"))

# 3. Limit of detection----

# 3.1. HRP2 as reference----

# 3.1.1 Distribution of pos and neg RDTs per range of hrp2 concentration----
plot.HRP2_hist<-iid_ldf %>% filter(test == "RDT" | test == "uRDT") %>%
  mutate(result2 = case_when(result=="Negative"~"g1",
                             test == "RDT"&result=="Positive"~"g2",
                             test == "uRDT"&result=="Positive"~"g3"))%>%
  ggplot(aes(x = log10(hrp2_pg_ml), color = result2, fill = result2)) +
  geom_histogram(position = "stack", bins = 60) +
  labs(y = "# of samples or iids",
       x = "HRP2 pg/mL in log10 scale") +
  facet_wrap(.~test, nrow = 2)+
  scale_color_manual(values = c("dodgerblue4", "coral", "firebrick3"))+
  scale_fill_manual(values = c("dodgerblue4", "coral", "firebrick3"))+
  theme_classic()+
  theme(plot.title = element_text(hjust = 0, size = 11)) +
  theme(axis.title.x = element_text(size = 10)) +
  theme(axis.title.y = element_text(size = 10)) +
  theme(axis.text = element_text(size = 10), legend.position = "none")+
  scale_x_continuous(breaks = 0:4, labels = c(1, 10, 100, 1000, 10000))

# 3.1.2. Proportion of pos and neg RDTs per range of hrp2 concentration----
plot.HRP2_prop<-iid_ldf %>% filter(test == "RDT" | test == "uRDT") %>%
  mutate(result2 = case_when(result=="Negative"~"g1",
                             test == "RDT"&result=="Positive"~"g2",
                             test == "uRDT"&result=="Positive"~"g3"))%>%
  fx_prop.hist.plot(bins = 60,
                    quant.var = "hrp2_pg_ml",
                    qual.var = "result2",
                    group.var = "test",
                    log.scale = TRUE,
                    colors = c("dodgerblue4","coral", "firebrick3"),
                    y.title = "Proportion of samples or iids",
                    x.title = "HRP2 pg/mL in log10 scale",
                    legend.title = "Result",
                    legend.position = "none")

# 3.1.3. LOD using HRP2 as reference----
lod.RDT.hrp2 <- iid_ldf %>% filter(test == "RDT") %>% select (hrp2_pg_ml, result)%>%
  fx_lod(concentration = "hrp2_pg_ml",
         result = "result",
         precision = 0.1,
         from = 0,
         to = 20000,
         by = 100,
         xlab = "HRP2 pg/mL",
         ylab = "Proportion of RDT positive results")


lod.uRDT.hrp2 <- iid_ldf %>% filter(test == "uRDT") %>% select (hrp2_pg_ml, result)%>%
  fx_lod(concentration = "hrp2_pg_ml",
         result = "result",
         precision = 0.1,
         from = 0,
         to = 20000,
         by = 100,
         xlab = "HRP2 pg/mL",
         ylab = "Proportion of uRDT positive results")


lod.RDT.hrp2$glm.lod$test <- "RDT"
lod.uRDT.hrp2$glm.lod$test <- "uRDT"

lod.hrp2 <- rbind(lod.RDT.hrp2$glm.lod,lod.uRDT.hrp2$glm.lod)

lod.hrp2.plot <- lod.hrp2 %>% ggplot(aes(x = concentration, y = prob, color = test)) +
  geom_line()+
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = test), alpha = 0.2, linetype = 3) +
  geom_vline(xintercept = c(lod.RDT.hrp2$lod$lod, lod.RDT.hrp2$lod$lower, lod.RDT.hrp2$lod$upper), color = c("coral", "coral", "coral"), size = 1, linetype = c(1,3,3)) +
  geom_vline(xintercept = c(lod.uRDT.hrp2$lod$lod, lod.uRDT.hrp2$lod$lower, lod.uRDT.hrp2$lod$upper), color = c("firebrick3", "firebrick3", "firebrick3"), size = 1, linetype = c(1,3,3)) +
  geom_point(x = lod.RDT.hrp2$lod$lod, y = 0.95, color = "black") +
  geom_point(x = lod.uRDT.hrp2$lod$lod, y = 0.95, color = "black") +
  geom_text(aes(x = lod.RDT.hrp2$lod$lod*1.1, y = 0.9), label = round(lod.RDT.hrp2$lod$lod, 2), color = "black", size = 3.2 )+
  geom_text(aes(x = lod.uRDT.hrp2$lod$lod*0.8, y = 0.95), label = round(lod.uRDT.hrp2$lod$lod, 2), color = "black", size = 3.2 )+
  theme_bw() +
  scale_color_manual(values = c("coral", "firebrick4"))+
  scale_fill_manual(values = c("coral", "firebrick4"))+
  labs(y = "Proportion of positive results",
       x = "HRP2 pg/mL",
       color = "Test",
       fill = "Test") +
  theme(plot.title = element_text(hjust = 0, size = 11),
        axis.title = element_text(size = 10),
        axis.text = element_text(size = 10),
        legend.position = "none") +
  scale_x_continuous(expand = c(0, 0.02*max(lod.hrp2$concentration)), breaks = seq(0, 20000, 2000)) +
  scale_y_continuous(expand = c(0, 0.02*max(lod.hrp2$upper)), breaks = seq(0, 1, 0.1))+
  theme(axis.text.x = element_text(angle = 45,vjust = 0.5))


# 3.2. qPCR as reference----

# 3.2.1. Distribution of pos and neg RDTs per range of parasite density----
plot.qPCR_hist<-iid_ldf %>% filter(test == "RDT" | test == "uRDT") %>%
  mutate(result2 = case_when(result=="Negative"~"g1",
                             test == "RDT"&result=="Positive"~"g2",
                             test == "uRDT"&result=="Positive"~"g3"))%>%
  ggplot(aes(x = log10(pden_mean), color = result2, fill = result2)) +
  geom_histogram(position = "stack", bins = 60) +
  labs(y = "# of samples or iids",
       x = "Parasites/uL in log10 scale") +
  facet_wrap(.~test, nrow = 2)+
  scale_color_manual(values = c("dodgerblue4", "coral", "firebrick3"))+
  scale_fill_manual(values = c("dodgerblue4", "coral", "firebrick3"))+
  theme_classic()+
  theme(plot.title = element_text(hjust = 0, size = 11)) +
  theme(axis.title.x = element_text(size = 10)) +
  theme(axis.title.y = element_text(size = 10)) +
  theme(axis.text = element_text(size = 10), legend.position = "none")+
  scale_x_continuous(breaks = 0:4, labels = c(1, 10, 100, 1000, 10000))

# 3.2.2. Proportion of pos and neg RDTs per range of parasite density----
plot.qPCR_prop<-iid_ldf %>% filter(test == "RDT" | test == "uRDT") %>%
  mutate(result2 = case_when(result=="Negative"~"g1",
                             test == "RDT"&result=="Positive"~"g2",
                             test == "uRDT"&result=="Positive"~"g3"))%>%
  fx_prop.hist.plot(bins = 20,
                    quant.var = "pden_mean",
                    qual.var = "result2",
                    group.var = "test",
                    log.scale = TRUE,
                    colors = c("dodgerblue4", "coral", "firebrick3"),
                    y.title = "Proportion of samples or iids",
                    x.title = "Parasites/uL in log10 scale",
                    legend.title = "Result",
                    legend.position = "none")

# 3.2.3. LOD using qPCR as reference----
lod.RDT.qpcr <- iid_ldf %>% filter(test == "RDT") %>% select (pden_mean, result)%>%
  fx_lod(concentration = "pden_mean",
         result = "result",
         upper = FALSE,
         precision = 10,
         from = 0,
         to = 20000,
         by = 100,
         xlab = "Parasites/mL",
         ylab = "Proportion of RDT positive results")

lod.uRDT.qpcr <- iid_ldf %>% filter(test == "uRDT") %>% select (pden_mean, result)%>%
  fx_lod(concentration = "pden_mean",
         result = "result",
         upper = FALSE,
         precision = 10,
         from = 0,
         to = 20000,
         by = 100,
         xlab = "Parasites/mL",
         ylab = "Proportion of uRDT positive results")

lod.RDT.qpcr$glm.lod$test <- "RDT"
lod.uRDT.qpcr$glm.lod$test <- "uRDT"

lod.qpcr <- rbind(lod.RDT.qpcr$glm.lod,lod.uRDT.qpcr$glm.lod)

lod.qPCR.plot <- lod.qpcr %>% ggplot(aes(x = concentration, y = prob, color = test)) +
  geom_line()+
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = test), alpha = 0.2, linetype = 3) +
  geom_vline(xintercept = c(lod.RDT.qpcr$lod$lod, lod.RDT.qpcr$lod$lower), color = c("coral", "coral"), size = 1, linetype = c(1,3)) +
  geom_vline(xintercept = c(lod.uRDT.qpcr$lod$lod, lod.uRDT.qpcr$lod$lower), color = c("firebrick3", "firebrick3"), size = 1, linetype = c(1,3)) +
  geom_point(x = lod.RDT.qpcr$lod$lod, y = 0.95, color = "black") +
  geom_point(x = lod.uRDT.qpcr$lod$lod, y = 0.95, color = "black") +
  geom_text(aes(x = lod.RDT.qpcr$lod$lod*1.1, y = 0.9), label = round(lod.RDT.qpcr$lod$lod, 2), color = "black", size = 3.2 )+
  geom_text(aes(x = lod.uRDT.qpcr$lod$lod*0.8, y = 0.95), label = round(lod.uRDT.qpcr$lod$lod, 2), color = "black", size = 3.2 )+
  theme_bw() +
  scale_color_manual(values = c("coral", "firebrick4"))+
  scale_fill_manual(values = c("coral", "firebrick4"))+
  labs(y = "Proportion of positive results",
       x = "Parasites/uL",
       color = "Test",
       fill = "Test") +
  theme(plot.title = element_text(hjust = 0, size = 11),
        axis.title = element_text(size = 10),
        axis.text = element_text(size = 10),
        legend.position = "none") +
  scale_x_continuous(expand = c(0, 0.02*max(lod.qpcr$concentration)), breaks = seq(0, 20000, 2000)) +
  scale_y_continuous(expand = c(0, 0.02*max(lod.qpcr$upper)), breaks = seq(0, 1, 0.1))+
  theme(axis.text.x = element_text(angle = 45,vjust = 0.5))

# 3.3. qPCR without HRP2 negative samples as reference----

# 3.3.1. Distribution of pos and neg RDTs per range of parasite density----
plot.qPCR_hist2<- iid_wdf %>%
  filter(!(qPCR == "Positive"&
             HRP2 !="Positive"))%>%
  pivot_longer(all_of(test),
               names_to = "test",
               values_to = "result") %>%
  filter(test == "RDT" | test == "uRDT") %>%
  mutate(result2 = case_when(result=="Negative"~"g1",
                             test == "RDT"&result=="Positive"~"g2",
                             test == "uRDT"&result=="Positive"~"g3"))%>%
  ggplot(aes(x = log10(pden_mean), color = result2, fill = result2)) +
  geom_histogram(position = "stack", bins = 60) +
  labs(y = "# of samples or iids",
       x = "Parasites/uL in log10 scale") +
  facet_wrap(.~test, nrow = 2)+
  scale_color_manual(values = c("dodgerblue4", "coral", "firebrick3"))+
  scale_fill_manual(values = c("dodgerblue4", "coral", "firebrick3"))+
  theme_classic()+
  theme(plot.title = element_text(hjust = 0, size = 11)) +
  theme(axis.title.x = element_text(size = 10)) +
  theme(axis.title.y = element_text(size = 10)) +
  theme(axis.text = element_text(size = 10), legend.position = "none")+
  scale_x_continuous(breaks = 0:4, labels = c(1, 10, 100, 1000, 10000))

# 3.3.2. Proportion of pos and neg RDTs per range of parasite density----
plot.qPCR_prop2 <- iid_wdf %>%
  filter(!(qPCR == "Positive"&
             HRP2 !="Positive"))%>%
  pivot_longer(all_of(test),
               names_to = "test",
               values_to = "result") %>%
  filter(test == "RDT" | test == "uRDT") %>%
  mutate(result2 = case_when(result=="Negative"~"g1",
                             test == "RDT"&result=="Positive"~"g2",
                             test == "uRDT"&result=="Positive"~"g3"))%>%
  fx_prop.hist.plot(bins = 20,
                    quant.var = "pden_mean",
                    qual.var = "result2",
                    group.var = "test",
                    log.scale = TRUE,
                    colors = c("dodgerblue4", "coral", "firebrick3"),
                    y.title = "Proportion of samples or iids",
                    x.title = "Parasites/uL in log10 scale",
                    legend.title = "Result",
                    legend.position = "none")

# 3.3.3. LOD using qPCR as reference----

lod.RDT.qpcr2 <- iid_wdf %>%
  filter(!(qPCR == "Positive"&
             HRP2 !="Positive") & pden_mean <= 1000)%>%
  pivot_longer(all_of(test),
               names_to = "test",
               values_to = "result") %>%
  filter(test == "RDT") %>% select (pden_mean, result)%>%
  fx_lod(concentration = "pden_mean",
         result = "result",
         upper = FALSE,
         precision = 10,
         from = 0,
         to = 500,
         by = .1,
         xlab = "Parasites/mL",
         ylab = "Proportion of RDT positive results")


lod.uRDT.qpcr2 <- iid_wdf %>%
  filter(!(qPCR == "Positive"&
             HRP2 !="Positive") & pden_mean <= 1000)%>%
  pivot_longer(all_of(test),
               names_to = "test",
               values_to = "result") %>%
  filter(test == "uRDT") %>% select (pden_mean, result)%>%
  fx_lod(concentration = "pden_mean",
         result = "result",
         upper = FALSE,
         precision = 10,
         from = 0,
         to = 500,
         by = .1,
         xlab = "Parasites/mL",
         ylab = "Proportion of uRDT positive results")


lod.RDT.qpcr2$glm.lod$test <- "RDT"
lod.uRDT.qpcr2$glm.lod$test <- "uRDT"

lod.qpcr2 <- rbind(lod.RDT.qpcr2$glm.lod,lod.uRDT.qpcr2$glm.lod)

lod.qPCR.plot2 <- lod.qpcr2 %>% ggplot(aes(x = concentration, y = prob, color = test)) +
  geom_line()+
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = test), alpha = 0.2, linetype = 3) +
  geom_vline(xintercept = c(lod.RDT.qpcr2$lod$lod, lod.RDT.qpcr2$lod$lower), color = c("coral", "coral"), size = 1, linetype = c(1,3)) +
  geom_vline(xintercept = c(lod.uRDT.qpcr2$lod$lod, lod.uRDT.qpcr2$lod$lower), color = c("firebrick3", "firebrick3"), size = 1, linetype = c(1,3)) +
  geom_point(x = lod.RDT.qpcr2$lod$lod, y = 0.95, color = "black") +
  geom_point(x = lod.uRDT.qpcr2$lod$lod, y = 0.95, color = "black") +
  geom_text(aes(x = lod.RDT.qpcr2$lod$lod*1.1, y = 0.9), label = round(lod.RDT.qpcr2$lod$lod, 2), color = "black", size = 3.2 )+
  geom_text(aes(x = lod.uRDT.qpcr2$lod$lod*0.8, y = 0.95), label = round(lod.uRDT.qpcr2$lod$lod, 2), color = "black", size = 3.2 )+
  theme_bw() +
  scale_color_manual(values = c("coral", "firebrick4"))+
  scale_fill_manual(values = c("coral", "firebrick4"))+
  labs(y = "Proportion of positive results",
       x = "Parasites/uL",
       color = "Test",
       fill = "Test") +
  theme(plot.title = element_text(hjust = 0, size = 11),
        axis.title = element_text(size = 10),
        axis.text = element_text(size = 10),
        legend.position = "none") +
  scale_x_continuous(expand = c(0, 0.02*max(lod.qpcr2$concentration)), breaks = seq(0, 500, 50), limits = c(0,500)) +
  scale_y_continuous(expand = c(0, 0.02*max(lod.qpcr2$upper)), breaks = seq(0, 1, 0.1))+
  theme(axis.text.x = element_text(angle = 45,vjust = 0.5))

# 3.4. Plot LOD----

plot2.lod <- plot_grid(plot_grid(plot.HRP2_hist, lod.hrp2.plot,
                             ncol = 2),
                   plot_grid(plot.qPCR_hist, lod.qPCR.plot,
                             ncol = 2),
                   labels = c("A", "B"), label_size = 11, label_fontface = "plain", nrow = 2)

supp.plot1.lod <- plot_grid(plot.HRP2_prop$prop.plot,
                                 plot.qPCR_prop$prop.plot,
                                 nrow = 2, labels = c("A", "B"), label_size = 11, label_fontface = "plain")

supp.plot2.lod <- plot_grid(plot.qPCR_hist2,
                                      plot.qPCR_prop2$prop.plot, lod.qPCR.plot2,
                                      ncol = 3, labels = c("A", "B", "C"), label_size = 11, label_fontface = "plain")


# 4. Performance of RDTs and uRDTs in detecting infected or recently exposed individuals----
# 4.1. HRP2 as standard----
performance_HRP2_all <- fx_performance(iid_wdf,
                                       reference = "HRP2",
                                       test = c("RDT", "uRDT"),
                                       pos = "Positive",
                                       neg = "Negative")

performance_HRP2.800high <- iid_wdf %>% filter(!(HRP2 == "Positive" & hrp2_pg_ml < 800))%>%
  fx_performance(reference = "HRP2",
                 test = c("RDT", "uRDT"),
                 pos = "Positive",
                 neg = "Negative")

perf_hrp2_metrics <- rbind(cbind(performance_HRP2_all$metrics,tibble(group = "All")),
                           cbind(performance_HRP2.800high$metrics,tibble(group = ">800"))
                           )

perf_hrp2_metrics$metric<-factor(perf_hrp2_metrics$metric,
                                 levels = c("sens", "spec", "ppv", "npv"))

plot.metrics.hrp2<-perf_hrp2_metrics%>%
  filter(metric != "accu")%>%
  mutate(group2=paste(test,group,sep="_"))%>%
  ggplot(aes(x = group, y = value, ymin = lower, ymax = upper, color = group2)) +
  geom_pointrange(size = 0.8) + 
  geom_hline(yintercept = 80, lty = 2) +  # add a dotted line at x=1 after flip
  coord_flip() +  # flip coordinates (puts labels on y axis)\
  scale_color_manual(values = c("dodgerblue1",
                                "dodgerblue3",
                                "firebrick1",
                                "firebrick3"))+
  scale_y_continuous(limits = c(0,100))+
  labs(title = "A) HRP2 as standard",
       x = "Test",
       y = "Mean (95% CI)")+
  theme_bw() +  # use a white background
  facet_grid(test ~ metric)+
  theme(plot.title = element_text(hjust = 0, size = 11)) +
  theme(axis.title.x = element_text(size = 10)) +
  theme(axis.title.y = element_blank()) +
  theme(axis.text = element_text(size = 10),
        axis.text.x = element_text(angle = 45, vjust = .5)) +
  theme(legend.position="none")

# 4.2. qPCR as standard----
performance_qPCR_all <- fx_performance(iid_wdf,
                                       reference = "qPCR",
                                       test = c("RDT", "uRDT"),
                                       pos = "Positive",
                                       neg = "Negative")

performance_qPCR.100high <- iid_wdf %>% filter(!(qPCR == "Positive" & pden_mean < 100))%>%
  fx_performance(reference = "qPCR",
                 test = c("RDT", "uRDT"),
                 pos = "Positive",
                 neg = "Negative")

performance_qPCR.2 <- iid_wdf %>%
  filter(!(qPCR == "Positive"&
             HRP2 !="Positive") & pden_mean <= 1000)%>%
  fx_performance(reference = "qPCR",
                 test = c("RDT", "uRDT"),
                 pos = "Positive",
                 neg = "Negative")

perf_qpcr_metrics <- rbind(cbind(performance_qPCR_all$metrics,tibble(group = "All")),
                           cbind(performance_qPCR.100high$metrics,tibble(group = ">100")),
                           cbind(performance_qPCR.2$metrics,tibble(group = "wo HRP2 neg"))
                           )

perf_qpcr_metrics$metric<-factor(perf_qpcr_metrics$metric,
                                 levels = c("sens", "spec", "ppv", "npv"))

perf_qpcr_metrics$group<-factor(perf_qpcr_metrics$group,
                                 levels = c("wo HRP2 neg", ">100", "All"))

plot.metrics.qpcr<-perf_qpcr_metrics%>%
  filter(metric != "accu")%>%
  mutate(group2=paste(test,group,sep="_"))%>%
  ggplot(aes(x = group, y = value, ymin = lower, ymax = upper, color = group2)) +
  geom_pointrange(size = 0.8) + 
  geom_hline(yintercept = 80, lty = 2) +  # add a dotted line at x=1 after flip
  coord_flip() +  # flip coordinates (puts labels on y axis)\
  scale_color_manual(values = c("dodgerblue2",
                                "dodgerblue3",
                                "dodgerblue1",
                                "firebrick2",
                                "firebrick3",
                                "firebrick1"))+
  scale_y_continuous(limits = c(0,100))+
  labs(title = "B) qPCR as standard",
       y = "Mean (95% CI)")+
  theme_bw() +  # use a white background
  facet_grid(test ~ metric)+
  theme(plot.title = element_text(hjust = 0, size = 11)) +
  theme(axis.title.x = element_text(size = 10)) +
  theme(axis.title.y = element_blank()) +
  theme(axis.text = element_text(size = 10),
        axis.text.x = element_text(angle = 45, vjust = .5)) +
  theme(legend.position="none")

# 4.3. Plot performance----
plot3.performance <- plot_grid(plot.metrics.hrp2,
          plot.metrics.qpcr, nrow = 2)

# 5. Performance of RDTs and uRDTs in detecting at EA hotspots----

# 5.1. Data managment----
# 5.1.1. hhid_df:----
# Number of iid screened and iid infected per house by test

iid_ldf %<>% # data take it from iid_ldf
  filter(test!="qual_diag")%>%
  filter(eaid != "10499014",
         eaid != "10699016",
         eaid != "10599010",
         eaid != "10599056")

hhid_ldf <- iid_ldf %>% # data take it from iid_ldf
  group_by(hhid, test, result) %>% 
  summarise(
    eaid = levels(as.factor(as.character(eaid))),
    n.iid = nlevels(as.factor(as.character(iid)))) %>% # new variable "n.iid": number of iids per each category
  pivot_wider(names_from = result, values_from = n.iid) %>% # convert to wide format to sum Negative and Positive results per hhid
  mutate(scr.iid = sum(c(Negative, Positive), na.rm = T), # count the number of screened iid
         inf.iid = sum(Positive, na.rm = T))%>% # count the number of
  select(hhid,eaid,test, scr.iid, inf.iid)

hhid_wdf <- hhid_ldf %>% pivot_wider(names_from = test, values_from = inf.iid)

# 5.1.2. Calculate the number of of hhid with infected indv per eaid detected by Test----
eaid_c.hh_df<-hhid_ldf%>%
  group_by(eaid, test) %>%
  summarise(n.cases.hh = levels(as.factor(inf.iid)),
            nhhid.eaid = sum(summary(as.factor(inf.iid))),
            n.hh = summary(as.factor(inf.iid)),
            freq = summary(as.factor(inf.iid))/length(inf.iid)
  )

eaid_c.hh_df%<>%pivot_wider(names_from = n.cases.hh, values_from = c(freq,n.hh))

for (i in names(eaid_c.hh_df[-1:-2])){
  eaid_c.hh_df[is.na(eaid_c.hh_df[[i]]),][[i]]<-0
  rm(i)
  }

eaid_c.hh_df2<-eaid_c.hh_df%>%pivot_longer(names(eaid_c.hh_df)[grepl("n.hh",names(eaid_c.hh_df))], names_to = "n.cases.hh", values_to = "n.hh")
eaid_c.hh_df2$n.cases.hh<-gsub("n.hh_","",eaid_c.hh_df2$n.cases.hh)
eaid_c.hh_df3<-eaid_c.hh_df%>%pivot_longer(names(eaid_c.hh_df)[grepl("freq",names(eaid_c.hh_df))], names_to = "n.cases.hh", values_to = "freq")

eaid_c.hh_df<-cbind(eaid_c.hh_df2[-4:-11],eaid_c.hh_df3[13])
eaid_c.hh_df$eaid<-as.factor(eaid_c.hh_df$eaid)

rm(list=c("eaid_c.hh_df2","eaid_c.hh_df3"))

# 5.1.3. Proportion of hhid with more tham 1 iid inf.----
eaid_pr.hh.iniid_df = NULL

for (x in levels(as.factor(hhid_ldf$test))){
  temp<-tibble("test"=x)
  for (n in levels(as.factor(hhid_ldf$eaid))){
    temp[["eaid"]]<-n
    for(m in 1:4){
      temp[["n.inf.hh"]]<-m
      temp[["Proportion"]]<-(nrow(hhid_ldf[hhid_ldf[["eaid"]]==n &
                                            hhid_ldf[["test"]]==x &
                                            hhid_ldf[["inf.iid"]]>=m,])/nrow(hhid_ldf[hhid_ldf[["eaid"]]==n&
                                                                                      hhid_ldf[["test"]]==x,]))
      eaid_pr.hh.iniid_df<-rbind(eaid_pr.hh.iniid_df,temp)
    }
  }
  rm(list= c("x","n","m"))
}

eaid_pr.hh.iniid_df$test<-as.factor(eaid_pr.hh.iniid_df$test)
eaid_pr.hh.iniid_df$n.inf.hh<-as.factor(eaid_pr.hh.iniid_df$n.inf.hh)

# 5.1.4. eaid_df: ----
#long format tbbile data.frame for eaid`s, just some variables will be used in order to facilitate its processing

eaid_wdf<-read.dta13("Data/EA_level_dataset_56.dta") # import ea data base

eaid_wdf <- eaid_wdf[c("ea_no","ea_pop_update","ea_local_cases_post_8wks","ea_incidence_rate")]
eaid_wdf <- as_tibble(eaid_wdf) # convert to a tibble
names(eaid_wdf) <- c("eaid", "popsize", "cases_8wks", "Incidence")# variables are renamed in order to facilitate its processing

eaid_wdf %<>% filter(eaid != "10499014",
                     eaid != "10699016",
                     eaid != "10599010",
                     eaid != "10599056")

# calculate the number of iid`s screened by eaid 

iid_wdf %<>% # Pick information from ii_df and store it into eaid_df
  filter(eaid != "10499014",
         eaid != "10699016",
         eaid != "10599010",
         eaid != "10599056")

eaid_wdf <- iid_wdf %>% # Pick information from ii_df and store it into eaid_df
  group_by(eaid) %>% # group information by eaid
  summarise(samsize = nlevels(as.factor(as.character(iid))),
            nhhid = nlevels(as.factor(as.character(hhid)))) %>% # summarise the number of iid`s per eaid and store that information in a varibale named samsize
  full_join(eaid_wdf, by = "eaid") # joint the new variable samsize with the information in eaid_df

eaid_wdf[is.na(eaid_wdf[["Incidence"]]),][["Incidence"]] <-
  hmean(eaid_wdf[!is.na(eaid_wdf[["Incidence"]]) & eaid_wdf[["Incidence"]] != 0,][["Incidence"]])

# 5.1.5. Calculate prevalence ---- 

prev_df <- NULL

for (i in levels(as.factor(iid_ldf$eaid))){
  temp <- tibble("eaid" = i)
  for (j in levels(as.factor(as.character(iid_ldf[iid_ldf[["eaid"]] == i,][["test"]])))){
    temp[["test"]] <- j
    temp[["ncases"]] <-nlevels(
      as.factor(
        as.character(
          iid_ldf[iid_ldf[["eaid"]] == i &
                    iid_ldf[["test"]] == j &
                    iid_ldf[["result"]] == "Positive",][["iid"]])))
    temp[["prev"]] <-100*(nlevels(
      as.factor(
        as.character(
          iid_ldf[iid_ldf[["eaid"]] == i &
                    iid_ldf[["test"]] == j &
                    iid_ldf[["result"]] == "Positive",][["iid"]]))))/eaid_wdf[eaid_wdf[["eaid"]] == i,][["samsize"]]
    prev_df <- rbind(prev_df,temp)
    
  }
  rm(list = c("temp", "i", "j"))
}

for (i in paste0("cases", sep = ".", levels(as.factor(prev_df$test)))){
  eaid_wdf[[i]]<-as.numeric(NA)  
  rm(i)
}

for (i in paste0("prev", sep = ".", levels(as.factor(prev_df$test)))){
  eaid_wdf[[i]]<-as.numeric(NA)
  rm(i)
}

for (i in levels(as.factor(eaid_wdf$eaid))){
  for (j in  levels(as.factor(prev_df$test))){
    eaid_wdf[eaid_wdf[["eaid"]] == i,][[paste0("cases", sep = ".", j)]] <- prev_df[prev_df[["eaid"]] == i & prev_df[["test"]] == j,][["ncases"]]
    eaid_wdf[eaid_wdf[["eaid"]] == i,][[paste0("prev", sep = ".", j)]] <- prev_df[prev_df[["eaid"]] == i & prev_df[["test"]] == j,][["prev"]]
  }
  rm(list = c("i", "j"))
}

# 5.1.6. add proportions of hhids with inf. iids ----

prop.names <- as.factor(mapply(function (x, y) paste("prop",x,y , sep = "."),
                               crossing(levels(eaid_pr.hh.iniid_df$test),
                                        levels(eaid_pr.hh.iniid_df$n.inf.hh))[1],
                               crossing(levels(eaid_pr.hh.iniid_df$test),
                                        levels(eaid_pr.hh.iniid_df$n.inf.hh))[2]))

for (i in levels(prop.names)){
  eaid_wdf[[i]]<-as.numeric(NA)
  rm(i)
}

for (i in levels(as.factor(eaid_wdf$eaid))){
  for (j in  levels(as.factor(eaid_pr.hh.iniid_df$test))){
    for (k in levels(as.factor(eaid_pr.hh.iniid_df$n.inf.hh))){
      eaid_wdf[eaid_wdf[["eaid"]] == i,][[paste("prop", j, k, sep = ".")]] <- 
        eaid_pr.hh.iniid_df[eaid_pr.hh.iniid_df[["eaid"]] == i &
                              eaid_pr.hh.iniid_df[["test"]] == j &
                              eaid_pr.hh.iniid_df[["n.inf.hh"]] == k,][["Proportion"]]
    }
  }
  rm(list = c("i", "j", "k"))
}


# 5.2. Defining hotspots at eaid level----
# 5.2.1 Criterion 1----
eaid_wdf <- mutate(eaid_wdf, Criterion1 = case_when(
  prev.HRP2 < quantile(eaid_wdf$prev.HRP2, probs = .75, na.rm = T)[1] &
    Incidence < 50 &
    prev.qPCR < 10 ~ "no",
  prev.HRP2 >= quantile(eaid_wdf$prev.HRP2, probs = .75, na.rm = T)[1] |
    Incidence >= 50 |
    prev.qPCR >= 10 ~ "yes"
))

# 5.2.2 Criterion 2: Clustered dendrogram----

HClust.Analysis <- fx_HClust.Analysis(data = eaid_wdf,
                                      names.variables = c("Incidence","prev.HRP2","prev.qPCR"),
                                      names.rows = "eaid",
                                      distance = "auto",
                                      method = "auto",
                                      NbC.index = "sdindex")

par(mfcol = c(1,1), mfrow = c(1,1))

clusters <- cutree(HClust.Analysis$h.clust, k=3)
dd.samp <- as.dendrogram(HClust.Analysis$h.clust)
set.seed(2)
dd.samp2 <- reorder(dd.samp, (sample(names(clusters[clusters==2]),1): sample(names(clusters[clusters==3]),1)))

clusters <- cutree(HClust.Analysis$h.clust, k=6)
Clust.Heatmap.Plots <-fx_HClust.Heatmap(data = eaid_wdf,
                                        names.variables = c("Incidence","prev.HRP2","prev.qPCR"),
                                        names.rows = "eaid",
                                        dendrogram = dd.samp2,
                                        clusters = clusters)

# define hotspots
eaid_wdf %<>% mutate(Criterion2 = case_when(eaid %in% unlist(Clust.Heatmap.Plots$dendro$labels %>% filter(x<=15) %>% select(label)) ~ "yes",
                                            !(eaid %in% unlist(Clust.Heatmap.Plots$dendro$labels %>% filter(x<=15) %>% select(label))) ~ "no"))

rm(list=c("dd.samp", "dd.samp2", "clusters"))




# 5.3. Performance at eaid level all criteria----
eaid_roc_analysis <- fx_pROC_analysis(response = c("Criterion1",
                                                   "Criterion2"),
                                      predictors = c("prev.uRDT",
                                                     "prev.RDT"),
                                      data = eaid_wdf,
                                      ci.method = "delong",
                                      boot.n = 10000,
                                      boot.stratified = T,
                                      best.method = "youden",
                                      conf.level = .95,
                                      x = .85,
                                      input = "spec")




# 6. hhid level analysis----

hhid_ldf2<-hhid_ldf %>% filter(scr.iid >= 2 & scr.iid < 20)

hhid_wdf2 <- hhid_ldf2 %>% pivot_wider(names_from = test, values_from = inf.iid)
hhid_wdf2 %<>% mutate(hh_hotspots = case_when(
  HRP2 < 2 & qPCR <1 ~ "no",
  HRP2 >= 2 | qPCR >= 1 ~ "yes")
)

hhid_ldf2<-full_join(hhid_ldf2, hhid_wdf2[c("hhid", "hh_hotspots1", "hh_hotspots2", "hh_hotspots3")], by = "hhid")

# 6.1. Number of infected invd per hhid`s by eaid----

# Calculate the number of of hhid with infected indv per eaid detected by Test

eaid_c.hh_df2<-hhid_ldf2%>%
  group_by(eaid, test) %>%
  summarise(n.cases.hh = levels(as.factor(inf.iid)),
            nhhid.eaid = sum(summary(as.factor(inf.iid))),
            n.hh = summary(as.factor(inf.iid)),
            freq = summary(as.factor(inf.iid))/length(inf.iid)
  )

eaid_c.hh_df2%<>%pivot_wider(names_from = n.cases.hh, values_from = c(freq,n.hh))

for (i in names(eaid_c.hh_df2[-1:-2])){
  eaid_c.hh_df2[is.na(eaid_c.hh_df2[[i]]),][[i]]<-0
}

eaid_c.hh_df3<-eaid_c.hh_df2%>%pivot_longer(names(eaid_c.hh_df2)[grepl("n.hh",names(eaid_c.hh_df2))], names_to = "n.cases.hh", values_to = "n.hh")
eaid_c.hh_df3$n.cases.hh<-gsub("n.hh_","",eaid_c.hh_df3$n.cases.hh)
eaid_c.hh_df4<-eaid_c.hh_df2%>%pivot_longer(names(eaid_c.hh_df2)[grepl("freq",names(eaid_c.hh_df2))], names_to = "n.cases.hh", values_to = "freq")

eaid_c.hh_df2<-cbind(eaid_c.hh_df3[-4:-11],eaid_c.hh_df4[13])
eaid_c.hh_df2$eaid<-as.factor(eaid_c.hh_df2$eaid)


rm(list = c("eaid_c.hh_df3","eaid_c.hh_df4"))

# 6.2. Performance at hhid level----
hhid_roc_analysis<- fx_pROC_analysis(response = "hh_hotspots",
                                     predictors = c("uRDT", "RDT"),
                                     data = hhid_wdf2,
                                     ci.method = "delong",
                                     boot.n = 10000,
                                     boot.stratified = T,
                                     best.method = "youden",
                                     conf.level = .95)
