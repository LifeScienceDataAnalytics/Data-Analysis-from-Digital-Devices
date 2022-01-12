rm(list=ls())
library(tidyverse)
library(bnlearn)
load("~/Documents/Documents_IT/VAMBN_Version5/randomForestPred.RData")
rm(list=setdiff(ls(),  c("targetNodes", "regResults", "discdata",
                         "mmseMeasures")))
targetNodes = unique(targetNodes)
targetNodes = unique(targetNodes)
data = inputDataVIS1
# ##set the working directory to the ALTOIDA model working directory
# setwd("~/Documents/Documents_IT/VAMBN_Version4/")
# #load the model trained on ALTOIDA data#
# finalBN<-readRDS('data/data_out/main_finalBN.rds')# load finalBN
# real<-readRDS('data/data_out/main_RealPPts.rds')# load full real dataset 
# real$SUBJID<-NULL
# fitted.bnlearn<-readRDS('data/data_out/main_finalBN_fitted.rds')
# tmp<-readRDS('data/data_out/data_all_imp.rds')# load demo dataset

finalModel = readRDS("~/Documents/Documents_IT/VAMBN_Version5/finalModel.rds")

#simulate the digital nodes for ADNI##
func.pred.digital.nodes = function(data, target){
  extractVariables = grep(".AGE.*|.*MMSE.*|.*DX.*|.*PTGENDER.*|SUBJID", colnames(data), value = TRUE)
  data = data[,extractVariables]
  names(data)[names(data)%in%c("SA_PTGENDER_VIS1","SA_PTGENDER_VIS6","SA_PTGENDER_VIS12","SA_PTGENDER_VIS18","SA_PTGENDER_VIS24","SA_PTGENDER_VIS36","SA_PTGENDER_VIS48")] = "SA_gender_VIS1"
  names(data)[names(data)%in%c("SA_AGE_VIS1","SA_AGE_VIS6","SA_AGE_VIS12","SA_AGE_VIS18","SA_AGE_VIS24","SA_AGE_VIS36","SA_AGE_VIS48")] = "SA_age_VIS1"
  names(data)[names(data)%in%c("SA_DX_VIS1","SA_DX_VIS6","SA_DX_VIS12","SA_DX_VIS18","SA_DX_VIS24","SA_DX_VIS36","SA_DX_VIS48")] = "SA_MCI_VIS1"
  #data$SA_gender_VIS1 = as.factor(data$SA_gender_VIS1)
  data$SA_gender_VIS1 = as.character(data$SA_gender_VIS1)
  data$SA_gender_VIS1[data$SA_gender_VIS1=="Female"] <- "2"
  data$SA_gender_VIS1[data$SA_gender_VIS1=="Male"] <- "1"
  data$SA_gender_VIS1 = factor(data$SA_gender_VIS1)
  data$SA_MCI_VIS1 = factor(data$SA_MCI_VIS1)
  if("CN" %in% unique(data$SA_MCI_VIS1)){
    levels(data$SA_MCI_VIS1)[levels(data$SA_MCI_VIS1)=="CN"] <- "0"
  }
  if("MCI" %in% unique(data$SA_MCI_VIS1)){
    levels(data$SA_MCI_VIS1)[levels(data$SA_MCI_VIS1)=="MCI"] <- "1"
  }
  if("Dementia" %in% unique(data$SA_MCI_VIS1)){
    levels(data$SA_MCI_VIS1)[levels(data$SA_MCI_VIS1)=="Dementia"] <- "0"
  }
  data$SA_MCI_VIS1 = factor(data$SA_MCI_VIS1, levels = c("0", "1"))
  #digData = data.frame()
  
  # nrmse_bn = nrmseDf[nrmseDf$target == targetNodes[i], "Mean"]
  # print(nrmse_bn)
  # nrmse_rf = regResultsFinal[regResultsFinal$target == targetNodes[i], "nrmseMean"]
  # print(nrmse_rf)
  pt = data["SUBJID"]
  data = data[,c("SA_gender_VIS1","SA_age_VIS1", "SA_MCI_VIS1", grep("^zcode.*MMSE", colnames(data), value = TRUE))]
  #if(nrmse_bn<nrmse_rf){
  #print("bayesian")
  # set.seed(123)
  # data[target] = predict(fitted.bnlearn, data = data, node = target, method = "bayes-lw")
  # 
  #}
  # else{
  #   print("rf")
  #   set.seed(123)
  #   #print(colnames(data))
  #   discdata <- discdata %>% 
  #     select(cogDomainsAltoida, mmseMeasures, "age", "gender", "AB+",
  #            "MCI", targetNodes[1])
  #   print(colnames(discdata))
  #   regr <- randomForest(discdata[,targetNodes[1]] ~ ., data = discdata, mtry = 3, 
  #                        importance = TRUE, na.action = na.omit)
  #   predict(regr, data)
  # }
  target = "zcode_BIT_AR_o2_ATT_VIS1"
  testPred <- predict(finalModel , data[,setdiff(colnames(data),target)])
  return(list(data,pt))
}

func.all.time.points = function(visData, visit){
  digDataVis = c()
  for(target in targetNodes){
    data = func.pred.digital.nodes(visData, target)
    digDataVis = c(digDataVis, data[[1]][target])
  }
  digDataVis = as.data.frame(digDataVis)
  digDataVis = cbind.data.frame(digDataVis, data[[2]])
  #colnames to be changed
  x = setdiff(colnames(digDataVis), "SUBJID")
  digDataVis <- digDataVis %>%
    rename_with(~ gsub('VIS[0-9]+', paste("VIS",visit,sep =""), .x))
  #digDataVis = digDataVis %>% rename_at(vars(-SUBJID), ~ gsub("VIS[0-9]+", paste("VIS",visit,sep ="")))
  return(digDataVis)
}
impADNIdata = readRDS("~/Documents/Documents_IT/ADNIVAMBN/VAMBNForADNI/data/data_out/data_all_imp.rds")
inputDataVIS1 = impADNIdata$stalone_VIS1
#inputDataVIS1[,grep("^zcode_MMSE", colnames(impADNIdata$stalone_VIS1))] <- lapply(inputDataVIS1[,grep("^zcode_MMSE", colnames(inputDataVIS1))], as.numeric)
digDataVIS1 = func.all.time.points(inputDataVIS1,discdata,1)

