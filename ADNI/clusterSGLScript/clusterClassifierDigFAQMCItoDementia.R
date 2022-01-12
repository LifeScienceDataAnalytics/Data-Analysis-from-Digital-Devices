library(openxlsx)
library(purrr)
library(rlist)
library(plyr)
library(stringr)
library(Hmisc)
library(readr)
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
adniClassifierMCItoDementia = read.csv("data_DigMCItoDementia.csv")
adniClassifierMCItoDementia = na.omit(adniClassifierMCItoDementia)
#mmseFeatures = grep("MMSE", colnames(adniClassifierCNtoMCI), value = TRUE)
faqFeatures = grep("FAQ", colnames(adniClassifierMCItoDementia), value = TRUE)
digFeatures = grep("AR|^Motor|BIT|DOT", colnames(adniClassifierMCItoDementia), value = TRUE)
adniClassifierMCItoDementia$X = NULL
dataDf = adniClassifierMCItoDementiaDigFAQ
func.classifier.sgl.best.param = function(dataDf,index, type, Class){
  xtrain = subset(dataDf, select=-c(Class))
  ytrain = dataDf$Class
  ntrain=length(ytrain)
  set.seed(123)
  train.ext=createFolds(ytrain,k=5,returnTrain=TRUE)
  test.ext=lapply(train.ext,function(x) (1:ntrain)[-x])
  seedVec = c()
  alphas = seq(0.1,1,0.1)
  alphaVec = c()
  meanRocVec = c()
  lambdasMin = c()
  for(alpha in alphas){
    print(alpha)
    alpha = 0.1
    rocVec = c()
    lambdas = c()
    for(i in 1:5){
      i = 1
      data = list(x=xtrain[train.ext[[i]],],y=ytrain[train.ext[[i]]])
      fit = cvSGL(data, index = index, type = "logit", maxit = 1000, thresh = 0.001,
                  min.frac = 0.05, nlam =20, gamma = 0.8, nfold = 10, standardize = TRUE,
                  verbose = FALSE, step = 1, reset = 10, alpha = alpha, lambdas = NULL,
                  foldid = NULL)
      minLambda = min(fit$fit$lambdas)
      trainFit = SGL(data, index = index, type = "logit", maxit = 1000, thresh = 0.001,
                     min.frac = 0.1,nlam=2, gamma = 0.8, standardize = TRUE,
                     verbose = FALSE, step = 1, reset = 10, alpha = alpha, lambdas=c(min(fit$fit$lambdas),1))
      yh = predictSGL(trainFit, as.matrix(xtrain[test.ext[[i]],]), 1)
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

adniClassifierMCItoDementiaDigFAQ = adniClassifierMCItoDementia[,c(digFeatures,faqFeatures,"Class")]
targetNodes = gsub("zcode_", "", targetNodes)
targetNodes = gsub("_VIS1", "", targetNodes)
groups = targetNodes
featureDesc =read_delim("~/Documents/Documents_IT/paper/Altoida_VAMBN_paper/ALTOIDA/fraunhofer_v5_technical_grouping.csv", ";", escape_double = FALSE, col_names = TRUE, trim_ws = TRUE)
targetNodesFeatures = c()
for(group in groups){
  features = as.vector(featureDesc[featureDesc$technical_group_name == group, "feature_name"])
  feat = features$feature_name 
  targetNodesFeatures = c(targetNodesFeatures,feat)
}
targetNodesFeatures = intersect(targetNodesFeatures, colnames(adniClassifierMCItoDementiaDigFAQ))
adniClassifierMCItoDementiaDigFAQ = adniClassifierMCItoDementiaDigFAQ[,c(targetNodesFeatures, faqFeatures,"Class")]
adniClassifierMCItoDementiaDigFAQ$MMSE = NULL
featureDesc = featureDesc[featureDesc$technical_group_name%in% groups,]
index = c()

group_name_vec = c()
for(names in colnames(adniClassifierMCItoDementiaDigFAQ[,targetNodesFeatures])){
  group_name = featureDesc[featureDesc$feature_name == names, "technical_group_name"]
  group_name = as.character(group_name)
  print(group_name)
  index = c(index,which(groups==group_name))
  group_name_vec = c(group_name_vec, group_name)
}
indexDf = cbind.data.frame(index,group_name_vec)
indexDf = unique(indexDf)

##index for Dig and MMSE
indexDf = rbind.data.frame(indexDf, c((max(index)+1), "FAQ"))
index = c(index, rep((max(index)+1),10))
indexDigFAQ = index
Class = "Class"
type = "DigFAQ"
paramDigFAQ = func.classifier.sgl.best.param(adniClassifierMCItoDementiaDigFAQ, index = indexDigFAQ, type = "DigFAQ", "Class")
save.image("paramDigFAQMCIDem.RData")

