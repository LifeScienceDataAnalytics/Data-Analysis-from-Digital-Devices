rm(list=ls())
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
##index, for Cognitive measures and MMSE measures, FAQ measures
adniClassifierCNtoMCICogMMSEFAQ = adniClassifierCNtoMCI[,c(cogDomains, mmseFeatures,faqFeatures,"Class")]
adniClassifierCNtoMCICogMMSEFAQ$MMSE = NULL
##index for Cog Domains, MMSE and FAQ
index = c(rep(1,9), rep(2,5), rep(3,10))
indexCogMMSEFAQ = index
Class = "Class"
type = "CogMMSEFAQ"
paramCogMMSEFAQ = func.classifier.sgl.best.param(adniClassifierCNtoMCICogMMSEFAQ, index = indexCogMMSEFAQ, type = type, "Class")
alphaCogMMSEFAQ = paramCogMMSEFAQ[[1]]
lambdaCogMMSEFAQ = paramCogMMSEFAQ[[2]]
load("~/ADNI_VAMBN_paper/sgl_R/resultsFromCluster/adniClassifier/paramCogMMSEFAQ.RData")


x = subset(adniClassifierCNtoMCICogMMSEFAQ, select=-c(Class))
y= adniClassifierCNtoMCICogMMSEFAQ$Class
finalData = list(x = x, y=y)
finalFit = SGL(finalData, index = index, type = "logit", maxit = 1000, thresh = 0.001,
               min.frac = 0.1, nlam =2, gamma = 0.8,standardize = TRUE,
               verbose = FALSE, step = 1, reset = 10, alpha = alphaCogMMSEFAQ, lambdas = c(lambdaCogMMSEFAQ,1)
)

save(finalFit, file = paste("sgl_R/modelsNew/","final_", type, "_model", ".RData", sep = ""))

coef = finalFit$beta[,match(min(finalFit$lambdas),finalFit$lambdas)]
absCoef = abs(coef)
coefDf = cbind.data.frame(colnames(x), coef, absCoef, index)
names(coefDf)[names(coefDf)=="colnames(x)"] = "names"
names(coefDf)[names(coefDf)=="index"] = "group"
coefDf <- coefDf[with(coefDf, order(-absCoef)), ]

library(ggplot2) 
coefDfNoZero = subset(coefDf[coefDf$absCoef>0,])
coefDfNoZeroCp = coefDfNoZero
png("sgl_R/plots/classifier_ADNI_Cognition_MMSE_FAQ.png", width = 300, height = 225, units='mm',res = 300) 
ggplot(coefDfNoZero, aes(x=reorder(names, absCoef), weight=absCoef, fill=as.factor(group))) + 
  geom_bar() +
  scale_fill_manual(labels = c("Digital Cognitive Domains","FAQ"), values = c("dark orange","orangered"), name = "Variable group") +
  ylab("absolute coefficients") +
  xlab("Feature Names")+
  coord_flip()
dev.off()

load("~/ADNI_VAMBN_paper/sgl_R/resultsFromCluster/adniClassifier/paramCogMMSEFAQ_only.RData")
alphaCogMMSEFAQ_only = paramCogMMSEFAQ_only[[1]]
lambdaCogMMSEFAQ_only = paramCogMMSEFAQ_only[[2]]
rocCogMMSEFAQ_only = paramCogMMSEFAQ_only[[3]]
type = "CogMMSEFAQ_only"
alphaCogMMSEFAQ_only = 0.1
x = subset(adniClassifierCNtoMCICogMMSEFAQ_only, select=-c(Class))
y= adniClassifierCNtoMCICogMMSEFAQ_only$Class
finalData = list(x = x, y=y)
finalFit = SGL(finalData, index = index, type = "logit", maxit = 1000, thresh = 0.001,
               min.frac = 0.1, nlam =2, gamma = 0.8,standardize = TRUE,
               verbose = FALSE, step = 1, reset = 10, alpha = alphaCogMMSEFAQ_only, lambdas = c(lambdaCogMMSEFAQ_only,1)
)

save(finalFit, file = paste("sgl_R/modelsNew/","final_", type, "_model", ".RData", sep = ""))
finalFit_only = finalFit
coefMMSE_only = finalFit_only$beta[,match(min(finalFit_only$lambdas),finalFit_only$lambdas)]
absCoefMMSE_only = abs(coefMMSE_only)
coefDfMMSE_only = cbind.data.frame(colnames(x), coefMMSE_only, absCoefMMSE_only, index)
names(coefDfMMSE_only)[names(coefDfMMSE_only)=="colnames(x)"] = "names"
names(coefDfMMSE_only)[names(coefDfMMSE_only)=="index"] = "group"
coefDfMMSE_only <- coefDfMMSE_only[with(coefDfMMSE_only, order(-absCoefMMSE_only)), ]


library(ggplot2) 
coefDfNoZeroMMSE_only = subset(coefDfMMSE_only[coefDfMMSE_only$absCoefMMSE_only>0,])
png("~/ADNI_VAMBN_paper/sgl_R/plots/classifier_ADNI_Cognition_MMSE_FAQ_only.png", width = 300, height = 225, units='mm',res = 300) 
ggplot(coefDfNoZeroMMSE_only, aes(x=reorder(names, absCoefMMSE_only), weight=absCoefMMSE_only, fill=as.factor(group))) + 
  geom_bar() +
  scale_fill_manual(labels = c("Cognitive Domains from DMs", "MMSE","FAQ"), values = c("dark orange","steelblue","orangered"), name = "Variable group") +
  ylab("absolute coefficients") +
  xlab("Feature Names")+
  coord_flip()
dev.off()

###overall feature importance of FAQ, Cog and MMSE##
overall_imp_mmse_only2 =  coefDfNoZeroMMSE_only[coefDfNoZeroMMSE_only$names %in% c(grep("MMSE", coefDfNoZeroMMSE_only$names, value = TRUE)),]
overall_imp_mmse_sum_only2 = colSums(overall_imp_mmse_only2['absCoefMMSE_only'])/colSums(coefDfNoZeroMMSE_only['absCoefMMSE_only'])
overall_imp_faq_only2 =  coefDfNoZeroMMSE_only[coefDfNoZeroMMSE_only$names %in% c(grep("FAQ", coefDfNoZeroMMSE_only$names, value = TRUE)),]
overall_imp_faq_sum_only2 = colSums(overall_imp_faq_only2['absCoefMMSE_only'])/colSums(coefDfNoZeroMMSE_only['absCoefMMSE_only'])
overall_imp_cog_only2 =  coefDfNoZeroMMSE_only[coefDfNoZeroMMSE_only$names %in% gsub("SA_", "",cogDomains),]
overall_imp_cog_sum_only2 = colSums(overall_imp_cog_only2['absCoefMMSE_only'])/colSums(coefDfNoZeroMMSE_only['absCoefMMSE_only'])

overalImpDf_only2 = cbind.data.frame(c("MMSE", "FAQ","Cognitive readouts"), c(overall_imp_mmse_sum_only2,overall_imp_faq_sum_only2, overall_imp_cog_sum_only2), c("","",""))
names(overalImpDf_only2) = c("groups", "importance", "feature")

plot1 = ggplot(coefDfNoZeroCp, aes(x=reorder(names, absCoef), weight=absCoef, fill=as.factor(group))) + 
  geom_bar(show.legend = FALSE) +
  scale_fill_manual(labels = c("Digital Cognitive Domains","FAQ"), values = c("orange", "orangered"))+ 
  ylab("absolute coefficients") +
  xlab("Feature Names")+
  coord_flip()
plot2 = ggplot(overalImpDf_only2, aes(fill = groups, y= importance, x = feature)) + geom_bar(position = "stack", stat = "identity", width = 0.1)+scale_fill_manual(values = c("orange","orangered", "steelblue"))+ylab("overall feature importance")
png("~/ADNI_VAMBN_paper/sgl_R/plots/classifier_ADNI_CogMMSEFAQ_overallImp_only.png", width = 300, height = 100, units='mm',res = 300) 
plot_grid(plot1, plot2, align = "h",labels = "AUTO", label_size = 12, ncol = 2)
dev.off()


auc <- function(data, index,bestAlpha, bestLambda, seed){
  xtrain = subset(data, select=-c(Class))
  ytrain = data$Class
  ntrain=length(ytrain)
  set.seed(123)
  train.ext=createFolds(ytrain,k=10,returnTrain=TRUE)
  test.ext=lapply(train.ext,function(x) (1:ntrain)[-x])
  rocVec = c()
  for(i in 1:10){
    #set.seed(seed)
    #val = i+121
    #i = 1
    data = list(x=xtrain[train.ext[[i]],], y=ytrain[train.ext[[i]]])
    fit = cvSGL(data, index = index, type = "logit", maxit = 1000, thresh = 0.001,
                min.frac = 0.05, nlam =2, gamma = 0.8, nfold = 10, standardize = TRUE,
                verbose = FALSE, step = 1, reset = 10, alpha = bestAlpha, lambdas = NULL,
                foldid = NULL)
    #save(fit, file = paste("sgl_R/modelsNew/",type, "_model","_",alpha, "_", i, ".RData", sep = ""))
    #minLambda = min(fit$fit$lambdas)
    trainFit = SGL(data, index = index, type = "logit", maxit = 1000, thresh = 0.001,
                   min.frac = 0.1,nlam=2, gamma = 0.8, standardize = TRUE,
                   verbose = FALSE, step = 1, reset = 10, alpha = bestAlpha, lambdas=c(bestLambda,1))
    yh = predictSGL(trainFit, as.matrix(xtrain[test.ext[[i]],]), 1)
    roc = (roc.curve(yh, weights.class0 = ytrain[test.ext[[i]]], curve = TRUE)$auc)
    print(roc)
    rocVec = c(rocVec, roc)
    #lambdas = c(lambdas, minLambda)
  }
  #maxRoc = match(max(rocVec),rocVec)
  #print(maxRoc)
  #bestLambda = lambdas[maxRoc]
  return(rocVec)
}

##cogMMSEFAQ
indexCogMMSEFAQ 
rocCogMMSEFAQ= auc(adniClassifierCNtoMCICogMMSEFAQ,indexCogMMSEFAQ, alphaCogMMSEFAQ, lambdaCogMMSEFAQ, seed = 123)
mean(rocCogMMSEFAQ)

boxplot(rocCogMMSEFAQ)

##dig, mmse, faq
load("~/ADNI_VAMBN_paper/sgl_R/resultsFromCluster/adniClassifier/paramDigFAQMMSE.RData")
Class = "Class"
type = "DigMMSEFAQ"
x = subset(adniClassifierCNtoMCIDigMMSEFAQ, select=-c(Class))
y= adniClassifierCNtoMCIDigMMSEFAQ$Class
alphaDigMMSEFAQ = paramDigFAQMMSE[[1]]
lambdaDigMMSEFAQ = paramDigFAQMMSE[[2]]
finalData = list(x = x, y=y)
finalFitDigMMSEFAQ = SGL(finalData, index = indexDigMMSEFAQ, type = "logit", maxit = 1000, thresh = 0.001,
                      min.frac = 0.1, nlam =2, gamma = 0.8, standardize = TRUE,
                      verbose = FALSE, step = 1, reset = 10, alpha = alphaDigMMSEFAQ, lambdas = c(lambdaDigMMSEFAQ,1)
)

save(finalFitDigMMSEFAQ, file = paste("sgl_R/modelsNew/","final_", type, "_model", ".RData", sep = ""))

##auc##
rocDigMMSEFAQ= auc(adniClassifierCNtoMCIDigMMSEFAQ,indexDigMMSEFAQ, alphaDigMMSEFAQ, lambdaDigMMSEFAQ, seed = 123)
mean(rocDigMMSEFAQ)


