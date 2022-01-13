setwd("~/Altoida_VAMBN_paper")

############################
############################ Dependencies and helper functions
############################
rm(list=ls())
library(bnlearn)
library(tidyverse)
library(plyr)
library(grid)
library(gridExtra)
library(psych)
library(FSA)
library(bestNormalize)


library(ggplot2)

load("~/Documents/Documents_IT/paper/ALTOIDA_VAMBN_paper/createbnAltoida.RData")
load("~/Documents/Documents_IT/paper/Altoida_VAMBN_paper/auc_upd.RData")
rm(list=setdiff(ls(),  c("coefDf","coefDfNoZero","coefDfNoZeroTop15","targetNodes", "dataCP", "demog", "diagnostics", "cogDomains", "dataDigAUC")))
altData = read.csv("altoidaDf.csv")
data = altData[altData$DX %in% c(0,1,2),]
data$DX[data$DX == 2] = 1
#names(data)[names(data) == "DX"] = "Class"
data$X = NULL

proDromalADMCIdata = altData[altData$DX %in% c(1,2,3),]
proDromalADMCIdata$DX[proDromalADMCIdata$DX == 2] = 1
proDromalADMCIdata$X = NULL



demMCIdata = altData[altData$DX %in% c(1,2,4),]
demMCIdata$DX[demMCIdata$DX == 2] = 1
demMCIdata$X = NULL

##only cognition
func.cog = function(cogDomains, data){
  cogDomains = gsub("_VIS1", "", cogDomains)
  cogDomains = gsub("SA_", "", cogDomains)
  dataCog = data[,c(cogDomains, "DX", "age", "gender")]
  resultsCog <- list()
  dataCogResult = data.frame()
  #plot(dataCog$age, dataCog$CognitiveProcessingSpeed)
  #cor(dataCog$age, dataCog$CognitiveProcessingSpeed)
  for(i in 1:length(dataCog[,1:9])){  
    print(i)
    set.seed(123)
    normalize_obj  = bestNormalize(dataCog[,i])
    dataCog[,i] = normalize_obj$x.t
    model <- lm(formula(paste(names(dataCog[i]), "~ DX + age + gender")), data = dataCog)
    summary = summary(model)
    confInt = confint(model)
    print(summary(model))
    print(confint(model))
    #print(resultsCog[[names(dataCog[i])]])
    dataCogResult[i,"Digital cognitive domains"] = names(dataCog[i])
    dataCogResult[i,"p.value"] = summary$coefficients[2,4]
    dataCogResult[i,"low.int"] = round(confInt[2,1],2)
    dataCogResult[i,"high.int"] = round(confInt[2,2],2)
  }
  dataCogResult["adj.p.value"] = round(p.adjust(dataCogResult$p.value, method = "holm"), 4)
  return(dataCogResult)
}

dataCogResultCNMCI = func.cog(cogDomains, data)
dataCogResultMCIproDromalAD = func.cog(cogDomains, proDromalADMCIdata)

dataCogResultMCIDem = func.cog(cogDomains, demMCIdata)


##only MMSE data with diagnosis, significant differences in digital data based on diagnosis
func.mmse = function(data){
  dataMMSE = data[,c(grep("MMSE", colnames(data), value = TRUE), "DX", "age", "gender")]
   dataMMSE$MMSE = NULL
  resultsMMSE <- list()
  dataMMSEResult = data.frame()
  for(i in 1:length(dataMMSE[,1:5])){  
    print(i)
    ##change for paper review Decemeber 
    # resultsMMSE[[names(dataMMSE[i])]] <- wilcox.test(formula(paste(names(dataMMSE[i]), "~ DX")), data = dataMMSE)
    # #resultsMMSE[[names(dataMMSE[i])]] <- dunnTest(formula(paste(names(dataMMSE[i]), "~ DX")), data = dataMMSE)
    # print(resultsMMSE[[names(dataMMSE[i])]])
    # dataMMSEResult[i,"MMSE sub-domain"] = names(dataMMSE[i])
    # dataMMSEResult[i,"p.value"] = resultsMMSE[[names(dataMMSE[i])]]["p.value"]
    set.seed(123)
    normalize_obj  = bestNormalize(dataMMSE[,i])
    dataMMSE[,i] = normalize_obj$x.t
    model <- lm(formula(paste(names(dataMMSE[i]), "~ DX + age + gender")), data = dataMMSE)
    summary = summary(model)
    confInt = confint(model)
    print(summary(model))
    print(confint(model))
    #print(resultsCog[[names(dataCog[i])]])
    dataMMSEResult[i,"Digital cognitive domains"] = names(dataMMSE[i])
    dataMMSEResult[i,"p.value"] = summary$coefficients[2,4]
    dataMMSEResult[i,"low.int"] = round(confInt[2,1],2)
    dataMMSEResult[i,"high.int"] = round(confInt[2,2],2)
    
  }
  dataMMSEResult["adj.p.value"] = round(p.adjust(dataMMSEResult$p.value, method = "holm"),4)
  return(dataMMSEResult)
}

dataMMSEResultCNMCI = func.mmse(data)
dataMMSEResultMCIproDromalAD = func.mmse(proDromalADMCIdata)##eception for i =5
dataMMSEResultMCIDDem = func.mmse(demMCIdata)

func.dig = function(dataDig){
  dataDigTaskResult = data.frame()
  resultsDig <- list()
  for(i in 1:length(dataDig[,1:15])){  
    print(i)
    ##change for paper review Decemeber 
    # resultsDig[[names(dataDig[i])]] <- wilcox.test(formula(paste(names(dataDig[i]), "~ DX")), data = dataDig)
    # #resultsMMSE[[names(dataMMSE[i])]] <- dunnTest(formula(paste(names(dataMMSE[i]), "~ DX")), data = dataMMSE)
    # print(resultsDig[[names(dataDig[i])]])
    # dataDigTaskResult[i,"Digital Task"] = names(dataDig[i])
    # dataDigTaskResult[i,"p.value"] = resultsDig[[names(dataDig[i])]]["p.value"]
    set.seed(123)
    normalize_obj  = bestNormalize(dataDig[,i])
    dataDig[,i] = normalize_obj$x.t
    model <- lm(formula(paste(names(dataDig[i]), "~ DX + age + gender")), data = dataDig)
    summary = summary(model)
    confInt = confint(model, "DX")
    print(summary(model))
    print(confint(model))
    #print(resultsCog[[names(dataCog[i])]])
    dataDigTaskResult[i,"Digital cognitive domains"] = names(dataDig[i])
    dataDigTaskResult[i,"p.value"] = summary$coefficients[2,4]
    dataDigTaskResult[i,"low.int"] = round(confInt[2,1],2)
    dataDigTaskResult[i,"high.int"] = round(confInt[2,2],2)
    
  }
  dataDigTaskResult["adj.p.value"] = round(p.adjust(dataDigTaskResult$p.value, method = "holm"),4)
  return(dataDigTaskResult)
}

##read the froup desription file 
technical_grouping_anonymous <- read_csv(".....csv")
names(technical_grouping_anonymous)[names(technical_grouping_anonymous) == "feature_name"] = "DigitalCognitiveDomains"
dataDigCNMCI = data[,c(coefDfNoZero$names, "DX", "age", "gender")]
dataDigTaskResultCNMCI = func.dig(dataDigCNMCI)
names(dataDigTaskResultCNMCI)[names(dataDigTaskResultCNMCI) == "Digital cognitive domains"] = "DigitalCognitiveDomains"
mergingDf = merge(dataDigTaskResultCNMCI,technical_grouping_anonymous, all.x = TRUE)
mergingDfCNMCI =mergingDf[ order(match(mergingDf$DigitalCognitiveDomains, dataDigTaskResultCNMCI$DigitalCognitiveDomains)), ]
mergingDfCNMCI$...3 = NULL

dataDigMCIprodromalAD = proDromalADMCIdata[,c(coefDfNoZero$names, "DX", "age", "gender")]
dataDigTaskResultMCIProdromalAD = func.dig(dataDigMCIprodromalAD)
names(dataDigTaskResultMCIProdromalAD)[names(dataDigTaskResultMCIProdromalAD) == "Digital cognitive domains"] = "DigitalCognitiveDomains"
mergingDfMCIProdromalAD = merge(dataDigTaskResultMCIProdromalAD,fraunhofer_v5_technical_grouping_anonymous, all.x = TRUE)
mergingDfMCIProdromalAD =mergingDfMCIProdromalAD[ order(match(mergingDfMCIProdromalAD$DigitalCognitiveDomains, dataDigTaskResultMCIProdromalAD$DigitalCognitiveDomains)), ]
mergingDfMCIProdromalAD$...3 = NULL

dataDigMCIDem = demMCIdata[,c(coefDfNoZero$names, "DX", "age", "gender")]
dataDigTaskResultMCIDem = func.dig(dataDigMCIDem)
names(dataDigTaskResultMCIDem)[names(dataDigTaskResultMCIDem) == "Digital cognitive domains"] = "DigitalCognitiveDomains"
mergingDfMCIDem = merge(dataDigTaskResultMCIDem,fraunhofer_v5_technical_grouping_anonymous, all.x = TRUE)
mergingDfMCIDem =mergingDfMCIDem[ order(match(mergingDfMCIDem$DigitalCognitiveDomains, dataDigTaskResultMCIDem$DigitalCognitiveDomains)), ]
mergingDfMCIDem$...3 = NULL


save.image("~/Altoida_VAMBN_paper/workspaces/linearModel.RData")
