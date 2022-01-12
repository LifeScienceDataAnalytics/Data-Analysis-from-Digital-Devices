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

load("~/Altoida_VAMBN_paper/createbnAltoida.RData")
rm(list=setdiff(ls(),  c("diagnostics", "cogDomains", "targetNodes")))
cogDomains = gsub("_VIS1", "", cogDomains)
cogDomains = gsub("SA_", "", cogDomains)
adniClassifierCNtoMCI = read.csv("data_DigCNtoMCI.csv")
adniClassifierCNtoMCI = na.omit(adniClassifierCNtoMCI)
adniClassifierCNtoMCI$X.1 = NULL
mmseFeatures = grep("MMSE", colnames(adniClassifierCNtoMCI), value = TRUE)
faqFeatures = grep("FAQ", colnames(adniClassifierCNtoMCI), value = TRUE)
digFeatures = grep("AR|^Motor|BIT|DOT", colnames(adniClassifierCNtoMCI), value = TRUE)
adniClassifierCNtoMCI$X = NULL

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
  for(alpha in alphas){
    print(alpha)
    rocVec = c()
    lambdas = c()
    for(i in 1:5){
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
adniClassifierCNtoMCICogMMSEFAQ_only = adniClassifierCNtoMCI[,c(cogDomains, "MMSE",faqFeatures,"Class")]

##index for Cog Domains, MMSE and FAQ
index = c(rep(1,9), rep(2,1), rep(3,10))
indexCogMMSEFAQ_only = index
Class = "Class"
type = "CogMMSEFAQ_only"
paramCogMMSEFAQ_only = func.classifier.sgl.best.param(adniClassifierCNtoMCICogMMSEFAQ_only, index = indexCogMMSEFAQ_only, type = type, "Class")
save.image("~/ADNI_VAMBN_paper/workspaces/paramCogMMSEFAQ_only.RData")
