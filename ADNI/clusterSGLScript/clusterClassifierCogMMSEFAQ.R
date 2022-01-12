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
library(PRROC)
library(naniar)
library(rle)
library(missForest)
library(SGL)
library(dplyr)
library(caret)
library(ggplot2)
library(tidyverse)
library(gridExtra)
require(gridExtra)
library(cowplot)

load("~/Documents/Documents_IT/paper/Altoida_VAMBN_paper/createbnAltoida.RData")
rm(list=setdiff(ls(),  c("diagnostics", "cogDomains", "targetNodes")))
cogDomains = gsub("_VIS1", "", cogDomains)
cogDomains = gsub("SA_", "", cogDomains)
adniClassifierCNtoMCI = read.csv("data_DigCNtoMCI.csv")
mmseFeatures = grep("MMSE", colnames(adniClassifierCNtoMCI), value = TRUE)
faqFeatures = grep("FAQ", colnames(adniClassifierCNtoMCI), value = TRUE)
digFeatures = grep("AR|^Motor|BIT|DOT", colnames(adniClassifierCNtoMCI), value = TRUE)
adniClassifierCNtoMCI$X = NULL

dataDf = adniClassifierCNtoMCIDigMMSEFAQ
func.classifier.sgl.best.param = function(dataDf,index, type, Class){
  xtrain = subset(dataDf, select=-c(Class))
  ytrain = dataDf$Class
  ntrain=length(ytrain)
  set.seed(123)
  train.ext=createFolds(ytrain,k=5,returnTrain=TRUE)
  test.ext=lapply(train.ext,function(x) (1:ntrain)[-x])
  seedVec = c()
  alphas = seq(0,1,0.1)
  alphaVec = c()
  meanRocVec = c()
  lambdasMin = c()
  #alpha = 0
  for(alpha in alphas){
    print(alpha)
    rocVec = c()
    lambdas = c()
    for(i in 1:5){
      set.seed(i+121)
      val = i+121
      xtrain_imp = missForest(xtrain[train.ext[[i]],653:662], maxiter = 2, ntree = 100, variablewise = FALSE,
                              decreasing = FALSE, verbose = FALSE,
                              mtry = floor(sqrt(ncol(xtrain))), replace = TRUE,
                              classwt = NULL, cutoff = NULL, strata = NULL,
                              sampsize = NULL, nodesize = NULL, maxnodes = NULL,
                              xtrue = NA)$ximp
      xtest_imp = missForest(xtrain[test.ext[[i]],653:662], maxiter = 2, ntree = 100, variablewise = FALSE,
                             decreasing = FALSE, verbose = FALSE,
                             mtry = floor(sqrt(ncol(xtrain))), replace = TRUE,
                             classwt = NULL, cutoff = NULL, strata = NULL,
                             sampsize = NULL, nodesize = NULL, maxnodes = NULL,
                             xtrue = NA)$ximp
      #xtrain_imp = xtrain_imp$ximp
      #xtest_imp = xtest_imp$ximp
      data = list(x=cbind.data.frame(xtrain[train.ext[[i]],1:652],xtrain_imp)
                  , y=ytrain[train.ext[[i]]])
      fit = cvSGL(data, index = index, type = "logit", maxit = 1000, thresh = 0.001,
                  min.frac = 0.05, nlam =20, gamma = 0.8, nfold = 10, standardize = TRUE,
                  verbose = FALSE, step = 1, reset = 10, alpha = alpha, lambdas = NULL,
                  foldid = NULL)
      #save(fit, file = paste("sgl_R/modelsNew/",type, "_model","_",alpha, "_", i, ".RData", sep = ""))
      minLambda = min(fit$fit$lambdas)
      trainFit = SGL(data, index = index, type = "logit", maxit = 1000, thresh = 0.001,
                     min.frac = 0.1,nlam=2, gamma = 0.8, standardize = TRUE,
                     verbose = FALSE, step = 1, reset = 10, alpha = alpha, lambdas=c(min(fit$fit$lambdas),1))
      x_test = cbind.data.frame(xtrain[test.ext[[i]],1:652],xtest_imp)
      yh = predictSGL(trainFit, as.matrix(x_test), 1)
      roc = (roc.curve(yh, weights.class0 = ytrain[test.ext[[i]]], curve = TRUE)$auc)
      print(roc)
      rocVec = c(rocVec, roc)
      lambdas = c(lambdas, minLambda)
    }
    meanRoc = mean(rocVec)
    meanRocVec = c(meanRocVec, meanRoc)
    lambdasMin = c(lambdasMin, min(lambdas))
    alphaVec = c(alphaVec, alpha)
  }
  maxRoc = match(max(meanRocVec),meanRocVec)
  print(maxRoc)
  bestAlpha = alphaVec[maxRoc]
  bestLambda = lambdasMin[maxRoc]
  print(bestAlpha)
  return(list(bestAlpha, bestLambda, meanRocVec))
}


##for digital MMSE and FAQ##
adniClassifierCNtoMCIDigMMSEFAQ = adniClassifierCNtoMCI[,c(digFeatures,mmseFeatures,faqFeatures,"Class")]
targetNodes = gsub("zcode_", "", targetNodes)
targetNodes = gsub("_VIS1", "", targetNodes)
groups = targetNodes
#data_all<-readRDS('data/data_out/data_all_imp.rds')
featureDesc =read_delim("~/Documents/Documents_IT/paper/Altoida_VAMBN_paper/ALTOIDA/fraunhofer_v5_technical_grouping.csv", ";", escape_double = FALSE, col_names = TRUE, trim_ws = TRUE)
targetNodesFeatures = c()
for(group in groups){
  features = as.vector(featureDesc[featureDesc$technical_group_name == group, "feature_name"])
  feat = features$feature_name 
  targetNodesFeatures = c(targetNodesFeatures,feat)
}
targetNodesFeatures = intersect(targetNodesFeatures, colnames(adniClassifierCNtoMCIDigMMSE))
adniClassifierCNtoMCIDigMMSEFAQ = adniClassifierCNtoMCIDigMMSEFAQ[,c(targetNodesFeatures, mmseFeatures, faqFeatures,"Class")]
adniClassifierCNtoMCIDigMMSEFAQ$MMSE = NULL
featureDesc = featureDesc[featureDesc$technical_group_name%in% groups,]
index = c()

group_name_vec = c()
for(names in colnames(adniClassifierCNtoMCIDigMMSEFAQ[,targetNodesFeatures])){
  group_name = featureDesc[featureDesc$feature_name == names, "technical_group_name"]
  group_name = as.character(group_name)
  print(group_name)
  index = c(index,which(groups==group_name))
  group_name_vec = c(group_name_vec, group_name)
}
indexDf = cbind.data.frame(index,group_name_vec)
indexDf = unique(indexDf)
#indexDf <- indexDf[order(indexDf$index),] 
##index for Dig and MMSE
indexDf = rbind.data.frame(indexDf, c((max(index)+1), "MMSE"))
indexDf = rbind.data.frame(indexDf, c((max(index)+2), "FAQ"))
index = c(index, rep((max(index)+1),5))
index = c(index, rep((max(index)+1),10))
indexDigMMSEFAQ = index
Class = "Class"
type = "DigMMSEFAQ"
paramDigMMSE = func.classifier.sgl.best.param(adniClassifierCNtoMCIDigMMSEFAQ, index = indexDigMMSEFAQ, type = "DigMMSEFAQ", "Class")

save.image(paramDigMMSE)