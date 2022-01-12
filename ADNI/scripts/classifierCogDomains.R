rm(list=ls())
library(openxlsx)
library(purrr)
library(rlist)
library(plyr)
library(stringr)
library(Hmisc)
library(readr)
##library to obtain longest common substring from multiple strings
library(PTXQC)
library(naniar)
library(rle)
library(missForest)
library(SGL)
load("~/Documents/Documents_IT/paper/Altoida_VAMBN_paper/createbnAltoida.RData")
rm(list=setdiff(ls(),  c("targetNodes", "dataCP", "demog", "diagnostics", "cogDomains")))
cogDomains_vis6 = gsub("_VIS1", "_VIS6", cogDomains)
cogDomains_vis12 = gsub("_VIS1", "_VIS12", cogDomains)
cogDomains_vis24 = gsub("_VIS1", "_VIS24", cogDomains)
cogDomains_vis36 = gsub("_VIS1", "_VIS36", cogDomains)
name<-'main'
data_out<-paste0('data/data_out/',name)
data_all<-readRDS('data/data_out/data_all_imp_dig_cog.rds')
orig<-data_all
orig<-orig %>% reduce(merge, by = 'SUBJID')

origDigCog = orig[,c("SUBJID",cogDomains,cogDomains_vis6, cogDomains_vis12, cogDomains_vis24, cogDomains_vis36)]
data_all_other<-readRDS('data/data_condensed.rds')
orig_other<-data_all_other %>% reduce(merge, by = 'SUBJID')
orig_other = orig_other[,c("SUBJID", grep("MMSE", colnames(orig_other), value = TRUE), grep("FAQ", colnames(orig_other), value = TRUE),
                           grep("DX", colnames(orig_other), value = TRUE))]
orig_other["MMSE_VIS1"] = rowSums(orig_other[,grep("MMSE.*\\_VIS1$", colnames(orig_other), value = TRUE)])
orig_other["MMSE_VIS6"] = rowSums(orig_other[,grep("MMSE.*\\_VIS6$", colnames(orig_other), value = TRUE)])
orig_other["MMSE_VIS12"] = rowSums(orig_other[,grep("MMSE.*\\_VIS12$", colnames(orig_other), value = TRUE)])
orig_other["MMSE_VIS24"] = rowSums(orig_other[,grep("MMSE.*\\_VIS24$", colnames(orig_other), value = TRUE)])
orig_other["MMSE_VIS36"] = rowSums(orig_other[,grep("MMSE.*\\_VIS36$", colnames(orig_other), value = TRUE)])
orig_all = merge(origDigCog, orig_other, all = TRUE)
#colnames(orig_all) = gsub("SA_", "", colnames(orig_all))

# colnames(orig_cog_mmse) = gsub('SA_', '',colnames(orig_cog_mmse))
# colnames(orig_cog_mmse) = gsub('_VIS1', '',colnames(orig_cog_mmse))
#set.seed(1234)
# orig_cog_mmse_imp = missForest(orig_cog_mmse, maxiter = 10, ntree = 100, variablewise = FALSE,
#                                decreasing = FALSE, verbose = FALSE,
#                                mtry = floor(sqrt(ncol(orig_cog_mmse))), replace = TRUE,
#                                classwt = NULL, cutoff = NULL, strata = NULL,
#                                sampsize = NULL, nodesize = NULL, maxnodes = NULL,
#                                xtrue = NA, parallelize = c('no', 'variables', 'forests'))

#orig_cog_mmse_imp_classifier = orig_cog_mmse_imp$ximp


func.processing.for.classifier = function(data,first_stage){
  data[,grep("\\_VIS36", colnames(data), value = TRUE)] = NULL
  names(data)[names(data)=="SUBJID"] = "RID"
  data$RID = as.character(data$RID)
  sel <- grep(".*\\_VIS", colnames(data), value = TRUE)
  data <- data %>% rename_at(sel,
                             funs(str_replace(., "\\_VIS", "\\.")))
  data <- data %>% rename_at(colnames(data),
                             funs(str_replace(., "SA\\_", "")))
  data = as.data.frame(data)
  data_long = reshape(data,varying = grep(".*\\.[0-9]+$", colnames(data), value = TRUE),
                      direction = "long",
                      #v.names = "Value",
                      idvar = "SUBJID",
                      timevar = "VISCODE")
  names(data_long)[names(data_long)=="SA_DX"] = "DX"
  aggDiagnosis <- ddply(data_long, .(RID), summarise, DX = toString(DX))
  patientidConv = c()
  for(i in 1:nrow(aggDiagnosis)){
    if(grepl("^MCI.*Dementia.*MCI$", aggDiagnosis[i,2])){
      patientidConv <- c(patientidConv, aggDiagnosis[i,"RID"])
    }
  }
  
  patientidConvCNtoMCI = c()
  for(i in 1:nrow(aggDiagnosis)){
    if(grepl("^CN.*MCI.*CN$", aggDiagnosis[i,2])){
      patientidConvCNtoMCI <- c(patientidConvCNtoMCI, aggDiagnosis[i,1])
    }
  }
  
  patientidConvMCItoCN = c()
  for(i in 1:nrow(aggDiagnosis)){
    if((grepl("^MCI.*CN.*MCI$", aggDiagnosis[i,2]))|(grepl("^MCI.*CN$", aggDiagnosis[i,2]))){
      patientidConvMCItoCN <- c(patientidConvMCItoCN, aggDiagnosis[i,1])
    }
  }
  
  patientidConvDementiaToMCI = c()
  for(i in 1:nrow(aggDiagnosis)){
    if((grepl("^Dementia.*MCI.*Dementia$", aggDiagnosis[i,2]))|(grepl("^Dementia.*MCI$", aggDiagnosis[i,2]))){
      patientidConvDementiaToMCI <- c(patientidConvDementiaToMCI, aggDiagnosis[i,1])
    }
  }
  
  patientidConv = c(patientidConv, patientidConvCNtoMCI, patientidConvMCItoCN, patientidConvDementiaToMCI)
  aggDiagnosis <- aggDiagnosis[which(!aggDiagnosis$RID %in% patientidConv),]
  
  nlSubjIDS = c()
  for(i in 1:nrow(aggDiagnosis)){
    if(grepl("^CN.*CN$", aggDiagnosis[i,2])){
      nlSubjIDS <- c(nlSubjIDS, aggDiagnosis[i,1])
    }
  }
  mciSubjIDS = c()
  for(i in 1:nrow(aggDiagnosis)){
    if(grepl("^MCI.*MCI$", aggDiagnosis[i,2])){
      mciSubjIDS <- c(mciSubjIDS, aggDiagnosis[i,1])
    }
  }
  
  dementiaSubjIDS = c()
  for(i in 1:nrow(aggDiagnosis)){
    if(grepl("^Dementia.*Dementia$", aggDiagnosis[i,2])){
      dementiaSubjIDS <- c(dementiaSubjIDS, aggDiagnosis[i,1])
    }
  }
  ##subset for patients remaining healthy and patients remaining mci
  if(first_stage == "CN"){
    data_long = data_long[data_long$RID %in% c(nlSubjIDS,mciSubjIDS),]
  }
  if(first_stage == "MCI"){
    data_long = data_long[data_long$RID %in% c(mciSubjIDS, dementiaSubjIDS),]
  }
  data_long$DX = factor(data_long$DX)
  rownames(data_long) = 1:nrow(data_long)
  data_long$VISCODE[data_long$VISCODE==1] <- 0
  data_long$ridViscode = paste(data_long$RID, data_long$VISCODE, sep = ",")
  rownames(data_long) = 1:nrow(data_long)
 
  adniDataForClassifier = data_long
  ids = table(adniDataForClassifier$RID)
  idsNotRepeating = c()
  for(i in 1:length(ids)){
    if(ids[i][[1]] == 1){
      idsNotRepeating = c(idsNotRepeating, as.integer(names(ids[i])))
    }
  }
  adniDataForClassifier = adniDataForClassifier[!adniDataForClassifier$RID %in% idsNotRepeating,]

  if(first_stage == "CN"){
    adniDataForClassifier$DX <- factor(adniDataForClassifier$DX,
                                       levels = c("CN","MCI"),
                                       labels = c(0,1))}
  if(first_stage == "MCI"){
    adniDataForClassifier$DX <- factor(adniDataForClassifier$DX,
                                       levels = c("MCI", "Dementia"),
                                       labels = c(0,1))}
 
  adniDataForClassifier$ridViscode = NULL
  
  rownames(adniDataForClassifier) = 1:nrow(adniDataForClassifier)
  
  return(adniDataForClassifier)
}
orig_cn_mci = func.processing.for.classifier(orig_all, first_stage = "CN")
orig_mci_dementia = func.processing.for.classifier(orig_all, first_stage = "MCI")
orig_cn_mciBL = orig_cn_mci[orig_cn_mci$VISCODE == 0,]
orig_mci_dementiaBL = orig_mci_dementia[orig_mci_dementia$VISCODE == 0,]


write.csv(orig_cn_mciBL, "adni_data_classifier_cn_mci.csv")
write.csv(orig_mci_dementiaBL, "adni_data_classifier_mci_dementia.csv")




