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

load("~/Documents/Documents_IT/paper/Altoida_VAMBN_paper/createbnAltoida.RData")
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

dataDf = adniClassifierCNtoMCIDigMMSEFAQ
func.classifier.sgl.best.param = function(dataDf,index, type, Class){
  xtrain = subset(dataDf, select=-c(Class))
  ytrain = dataDf$Class
  ntrain=length(ytrain)
  set.seed(123)
  train.ext=createFolds(ytrain,k=5,returnTrain=TRUE)
  test.ext=lapply(train.ext,function(x) (1:ntrain)[-x])
  seedVec = c()
  alphas = seq(0.1,1,0.1)
  #alphas = setdiff(alphas,"0.6")
  alphaVec = c()
  
  meanRocVec = c()
  lambdasMin = c()
  #alpha = 0
  for(alpha in alphas){
    print(alpha)
    #alpha = 0
    #alpha = 0.6
    rocVec = c()
    lambdas = c()
    for(i in 1:5){
      data = list(x=xtrain[train.ext[[i]],],y=ytrain[train.ext[[i]]])
      fit = cvSGL(data, index = index, type = "logit", maxit = 1000, thresh = 0.001,
                  min.frac = 0.05, nlam =20, gamma = 0.8, nfold = 10, standardize = TRUE,
                  verbose = FALSE, step = 1, reset = 10, alpha = alpha, lambdas = NULL,
                  foldid = NULL)
      #save(fit, file = paste("sgl_R/modelsNew/",type, "_model","_",alpha, "_", i, ".RData", sep = ""))
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
    #print(max(meanRocVec))
    alphaVec = c(alphaVec, alpha)
  }
  maxRoc = match(max(meanRocVec),meanRocVec)
  print(maxRoc)
  bestAlpha = alphaVec[maxRoc]
  bestLambda = lambdasMin[maxRoc]
  print(bestAlpha)
  return(list(bestAlpha, bestLambda, meanRocVec))
}



# ##index, for Cognitive measures and MMSE measures
# adniClassifierCNtoMCICogMMSE = adniClassifierCNtoMCI[,c(cogDomains, mmseFeatures, "Class")]
# adniClassifierCNtoMCICogMMSE$MMSE = NULL
# ##index for cog Domains and MMSE
# index = c(rep(1,9), rep(2,5))
# Class = "Class"
# type = "CogMMSE"
# paramCogMMSE = func.classifier.sgl.best.param(adniClassifierCNtoMCICogMMSE, index = index, type = type, "Class")
# #boxplot(rocCogMMSE, las = 2, xlab = "", ylab = "AUC", ylim = c(0.5,1),col = "skyblue")
# alphaCogMMSE = paramCogMMSE[[1]]
# lambdaCogMMSE = paramCogMMSE[[2]]
# #print(mean(sum(rocCogMMSE)/10))
# #maxRocCogMMSE = match(max(rocCogMMSE), rocCogMMSE)
# 
# ##load the model having minimum roc
# #load(paste("sgl_R/models/",type, "_model_" ,maxRocCogMMSE, ".RData", sep = ""))
# #x = read.csv(paste("sgl_R/data/",type,"_", maxRocCogMMSE, ".csv", sep = ""))
# #y= adniClassifierCNtoMCICogMMSE$Class
# #x$X = NULL
# x = subset(adniClassifierCNtoMCICogMMSE, select=-c(Class))
# y= adniClassifierCNtoMCICogMMSE$Class
# finalData = list(x = x, y=y)
# finalFit = SGL(finalData, index = index, type = "logit", maxit = 1000, thresh = 0.001,
#                  min.frac = 0.1, nlam =2, gamma = 0.8,standardize = TRUE,
#                  verbose = FALSE, step = 1, reset = 10, alpha = alphaCogMMSE, lambdas = c(lambdaCogMMSE,1)
#                  )
# 
# save(finalFit, file = paste("sgl_R/modelsNew/","final_", type, "_model", ".RData", sep = ""))
# 
# coef = finalFit$beta[,match(min(finalFit$lambdas),finalFit$lambdas)]
# absCoef = abs(coef)
# coefDf = cbind.data.frame(colnames(x), coef, absCoef, index)
# names(coefDf)[names(coefDf)=="colnames(x)"] = "names"
# names(coefDf)[names(coefDf)=="index"] = "group"
# coefDf <- coefDf[with(coefDf, order(-absCoef)), ]
# 
# 
# library(ggplot2) 
# coefDfNoZero = subset(coefDf[coefDf$absCoef>0,])
# png("sgl_R/plots/classifier_ADNI_Cognition_MMSE.png", width = 300, height = 225, units='mm',res = 300) 
# ggplot(coefDfNoZero, aes(x=reorder(names, absCoef), weight=absCoef, fill=as.factor(group))) + 
#   geom_bar() +
#   scale_fill_manual(labels = c("Digital Cognitive Domains","MMSE"), values = c("dark orange","steelblue"), name = "Variable group") +
#   ylab("absolute coefficients") +
#   xlab("Feature Names")+
#   coord_flip()
# dev.off()
# ###overall feature importance of Cog and MMSE##
# overall_imp_mmse =  coefDfNoZero[coefDfNoZero$names %in% c(grep("MMSE", coefDfNoZero$names, value = TRUE)),]
# overall_imp_mmse_sum = colSums(overall_imp_mmse['absCoef'])/colSums(coefDfNoZero['absCoef'])
# overall_imp_cog =  coefDfNoZero[coefDfNoZero$names %in% gsub("SA_", "",cogDomains),]
# overall_imp_cog_sum = colSums(overall_imp_cog['absCoef'])/colSums(coefDfNoZero['absCoef'])
# 
# overalImpDf = cbind.data.frame(c("MMSE", "Cog from DMs"), c(overall_imp_mmse_sum, overall_imp_cog_sum), c("",""))
# names(overalImpDf) = c("groups", "importance", "feature")
# 
# coefDfNoZeroCp = coefDfNoZero
# plot1 = ggplot(coefDfNoZeroCp, aes(x=reorder(names, absCoef), weight=absCoef, fill=as.factor(group))) + 
#   geom_bar(show.legend = FALSE) +
#   scale_fill_manual(labels = c("Cog from DMs", "MMSE"), values = c("orange", "steelblue"))+ 
#   ylab("absolute coefficients") +
#   xlab("Feature Names")+
#   coord_flip()
# plot2 = ggplot(overalImpDf, aes(fill = groups, y= importance, x = feature)) + geom_bar(position = "stack", stat = "identity", width = 0.1)+scale_fill_manual(values = c("orange", "steelblue"))+ylab("overall feature importance")
# png("sgl_R/plots/classifier_ADNI_CogMMSE_overallImp.png", width = 300, height = 100, units='mm',res = 300) 
# #plot_grid(plot1, plot2, align = "h",labels = "AUTO", label_size = 12, rel_widths = c(1.5,1))
# grid.arrange(plot1, plot2,ncol=2)
# dev.off()

# ##index, for Cognitive measures and MMSE measures only
# adniClassifierCNtoMCICogMMSE_only = adniClassifierCNtoMCI[,c(cogDomains, "MMSE", "Class")]
# #adniClassifierCNtoMCICogMMSE$MMSE = NULL
# ##index for cog Domains and MMSE
# index = c(rep(1,9), 2)
# indexCogMMSE_only = index 
# Class = "Class"
# type = "CogMMSE_only"
# paramCogMMSE_only = func.classifier.sgl.best.param(adniClassifierCNtoMCICogMMSE_only, index = indexCogMMSE_only, type = type, "Class")
# alphaCogMMSE_only = paramCogMMSE_only[[1]]
# lambdaCogMMSE_only = paramCogMMSE_only[[2]]
# #boxplot(rocCogMMSE_only, las = 2, xlab = "", ylab = "AUC", ylim = c(0.5,1),col = "skyblue")
# 
# #print(mean(sum(rocCogMMSE_only)/10))
# #maxRocCogMMSE_only = match(max(rocCogMMSE_only), rocCogMMSE_only)
# 
# ##load the model having minimum roc
# #load(paste("sgl_R/models/",type, "_model_" ,maxRocCogMMSE_only, ".RData", sep = ""))
# #x = read.csv(paste("sgl_R/data/",type,"_", maxRocCogMMSE_only, ".csv", sep = ""))
# #y= adniClassifierCNtoMCICogMMSE_only$Class
# #x$X = NULL
# x = subset(adniClassifierCNtoMCICogMMSE_only, select=-c(Class))
# y= adniClassifierCNtoMCICogMMSE_only$Class
# finalData_only = list(x = x, y=y)
# finalFit_only = SGL(finalData_only, index = index, type = "logit", maxit = 1000, thresh = 0.001,
#                min.frac = 0.1, nlam =2, gamma = 0.8,standardize = TRUE,
#                verbose = FALSE, step = 1, reset = 10, alpha = alphaCogMMSE_only, lambdas = c(lambdaCogMMSE_only,1)
# )
# 
# save(finalData_only, file = paste("sgl_R/modelsNew/","final_", type, "_model", ".RData", sep = ""))
# 
# coef_only = finalFit_only$beta[,match(min(finalFit_only$lambdas),finalFit_only$lambdas)]
# absCoef_only = abs(coef_only)
# coefDf_only = cbind.data.frame(colnames(x), coef_only, absCoef_only, index)
# names(coefDf_only)[names(coefDf_only)=="colnames(x)"] = "names"
# names(coefDf_only)[names(coefDf_only)=="index"] = "group"
# coefDf_only <- coefDf_only[with(coefDf_only, order(absCoef_only)), ]
# 
# 
# library(ggplot2) 
# coefDfNoZero_only = subset(coefDf_only[coefDf_only$absCoef_only>0,])
# png("sgl_R/plots/classifier_ADNI_Cognition_MMSE_only.png", width = 300, height = 225, units='mm',res = 300) 
# ggplot(coefDfNoZero_only, aes(x=reorder(names, absCoef_only), weight=absCoef_only, fill=as.factor(group))) + 
#   geom_bar() +
#   scale_fill_manual(labels = c("Cognitive Domains from DMs","MMSE"), values = c("dark orange","steelblue"), name = "Variable group") +
#   ylab("absolute coefficients") +
#   xlab("Feature Names")+
#   coord_flip()
# dev.off()
# ###overall feature importance of Cog and MMSE##
# overall_imp_mmse_only =  coefDfNoZero_only[coefDfNoZero_only$names %in% c(grep("MMSE", coefDfNoZero_only$names, value = TRUE)),]
# overall_imp_mmse_sum_only = colSums(overall_imp_mmse_only['absCoef_only'])/colSums(coefDfNoZero_only['absCoef_only'])
# overall_imp_cog_only =  coefDfNoZero_only[coefDfNoZero_only$names %in% gsub("SA_", "",cogDomains),]
# overall_imp_cog_sum_only = colSums(overall_imp_cog_only['absCoef_only'])/colSums(coefDfNoZero_only['absCoef_only'])
# 
# overalImpDf_only = cbind.data.frame(c("MMSE", "Cog from DMs"), c(overall_imp_mmse_sum_only, overall_imp_cog_sum_only), c("",""))
# names(overalImpDf_only) = c("groups", "importance", "feature")
# 
# coefDfNoZeroCp_only = coefDfNoZero_only
# plot1 = ggplot(coefDfNoZeroCp, aes(x=reorder(names, absCoef), weight=absCoef, fill=as.factor(group))) + 
#   geom_bar(show.legend = FALSE) +
#   scale_fill_manual(labels = c("Cog from DMs", "MMSE"), values = c("orange", "steelblue"))+ 
#   ylab("absolute coefficients") +
#   xlab("Feature Names")+
#   coord_flip()
# plot2 = ggplot(overalImpDf_only, aes(fill = groups, y= importance, x = feature)) + geom_bar(position = "stack", stat = "identity", width = 0.1)+scale_fill_manual(values = c("orange", "steelblue"))+ylab("overall feature importance")
# png("sgl_R/plots/classifier_ADNI_CogMMSE_overallImp_only.png", width = 300, height = 100, units='mm',res = 300) 
# #plot_grid(plot1, plot2, align = "h",labels = "AUTO", label_size = 12, rel_widths = c(1.5,1))
# grid.arrange(plot1, plot2,ncol=2)
# dev.off()

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
lambdaCogMMSEFAQ =  paramCogMMSEFAQ[[2]]
#boxplot(rocCogMMSEFAQ, las = 2, xlab = "", ylab = "AUC", ylim = c(0.5,1),col = "skyblue")

#print(mean(sum(rocCogMMSEFAQ)/10))
#maxRocCogMMSEFAQ = match(max(rocCogMMSEFAQ), rocCogMMSEFAQ)

##load the model having minimum roc
#load(paste("sgl_R/models/",type, "_model_" ,maxRocCogMMSEFAQ, ".RData", sep = ""))
#x = read.csv(paste("sgl_R/data/",type,"_", maxRocCogMMSEFAQ, ".csv", sep = ""))
x = subset(adniClassifierCNtoMCICogMMSEFAQ, select=-c(Class))
x = missForest(x, maxiter = 10, ntree = 100, variablewise = FALSE,
                        decreasing = FALSE, verbose = FALSE,
                        mtry = floor(sqrt(ncol(x))), replace = TRUE,
                        classwt = NULL, cutoff = NULL, strata = NULL,
                        sampsize = NULL, nodesize = NULL, maxnodes = NULL,
                        xtrue = NA)
#x = x$ximp
x$X = NULL
x = x$ximp
y= adniClassifierCNtoMCICogMMSEFAQ$Class
finalData = list(x = x, y=y)
finalFit = SGL(finalData, index = index, type = "logit", maxit = 1000, thresh = 0.001,
                 min.frac = 0.1, nlam =2, gamma = 0.8, standardize = TRUE,
                 verbose = FALSE, step = 1, reset = 10, alpha = alphaCogMMSEFAQ, lambdas = c(lambdaCogMMSEFAQ,1)
                 )

save(finalFit, file = paste("sgl_R/modelsNew/","final_", type, "_model",".RData", sep = ""))

coef = finalFit$beta[,match(min(finalFit$lambdas),finalFit$lambdas)]
absCoef = abs(coef)
coefDf = cbind.data.frame(colnames(x), coef, absCoef, index)
names(coefDf)[names(coefDf)=="colnames(x)"] = "names"
names(coefDf)[names(coefDf)=="index"] = "group"
coefDf <- coefDf[with(coefDf, order(-absCoef)), ]
#coefDf["label"] <- indexDf$group_name_vec[match(coefDf$group, indexDf$index)]


library(ggplot2) 
coefDfNoZero = subset(coefDf[coefDf$absCoef>0,])
png("sgl_R/plots/classifier_ADNI_Cognition_MMSE_FAQ.png", width = 300, height = 225, units='mm',res = 300) 
ggplot(coefDfNoZero, aes(x=reorder(names, absCoef), weight=absCoef, fill=as.factor(group))) + 
  geom_bar() +
  scale_fill_manual(labels = c("Cognitive Domains from DMs","FAQ"), values = c("dark orange","sky blue"), name = "Variable group") +
  ylab("absolute coefficients") +
  xlab("Feature Names")+
  coord_flip()
dev.off()


###overall feature importance of FAQ, Cog and MMSE##
overall_imp_mmse =  coefDfNoZero[coefDfNoZero$names %in% c(grep("MMSE", coefDfNoZero$names, value = TRUE)),]
overall_imp_mmse_sum = colSums(overall_imp_mmse['absCoef'])/colSums(coefDfNoZero['absCoef'])
overall_imp_faq =  coefDfNoZero[coefDfNoZero$names %in% c(grep("FAQ", coefDfNoZero$names, value = TRUE)),]
overall_imp_faq_sum = colSums(overall_imp_faq['absCoef'])/colSums(coefDfNoZero['absCoef'])
overall_imp_cog =  coefDfNoZero[coefDfNoZero$names %in% gsub("SA_", "",cogDomains),]
overall_imp_cog_sum = colSums(overall_imp_cog['absCoef'])/colSums(coefDfNoZero['absCoef'])

overalImpDf = cbind.data.frame(c("MMSE", "FAQ","Cog from DMs"), c(overall_imp_mmse_sum,overall_imp_faq_sum, overall_imp_cog_sum), c("","",""))
names(overalImpDf) = c("groups", "importance", "feature")


##index, for Cognitive measures and MMSE measures only, FAQ measures
adniClassifierCNtoMCICogMMSEFAQ_only = adniClassifierCNtoMCI[,c(cogDomains, "MMSE",faqFeatures,"Class")]

##index for Cog Domains, MMSE and FAQ
index = c(rep(1,9), rep(2,1), rep(3,10))
indexCogMMSEFAQ_only = index
Class = "Class"
type = "CogMMSEFAQ_only"
paramCogMMSEFAQ_only = func.classifier.sgl.best.param(adniClassifierCNtoMCICogMMSEFAQ_only, index = indexCogMMSEFAQ_only, type = type, "Class")
alphaCogMMSEFAQ_only = paramCogMMSEFAQ_only[[1]]
lambdaCogMMSEFAQ_only = paramCogMMSEFAQ_only[[2]]
#boxplot(rocCogMMSEFAQ_only, las = 2, xlab = "", ylab = "AUC", ylim = c(0.5,1),col = "skyblue")

#print(mean(sum(rocCogMMSEFAQ_only)/10))
#maxRocCogMMSEFAQ_only = match(max(rocCogMMSEFAQ_only), rocCogMMSEFAQ_only)

##load the model having minimum roc
#load(paste("sgl_R/models/",type, "_model_" ,maxRocCogMMSEFAQ_only, ".RData", sep = ""))
#x = read.csv(paste("sgl_R/data/",type,"_", maxRocCogMMSEFAQ_only, ".csv", sep = ""))
x = subset(adniClassifierCNtoMCICogMMSEFAQ_only, select=-c(Class))
x = missForest(x, maxiter = 10, ntree = 100, variablewise = FALSE,
               decreasing = FALSE, verbose = FALSE,
               mtry = floor(sqrt(ncol(x))), replace = TRUE,
               classwt = NULL, cutoff = NULL, strata = NULL,
               sampsize = NULL, nodesize = NULL, maxnodes = NULL,
               xtrue = NA)
x = x$ximp
y= adniClassifierCNtoMCICogMMSEFAQ_only$Class

finalDataMMSE_only = list(x = x, y=y)
finalFitMMSE_only = SGL(finalDataMMSE_only, index = index, type = "logit", maxit = 1000, thresh = 0.001,
               min.frac = 0.1, nlam =2, gamma = 0.8, standardize = TRUE,
               verbose = FALSE, step = 1, reset = 10, alpha = alphaCogMMSEFAQ_only, lambdas = c(lambdaCogMMSEFAQ_only,1)
)

save(finalFitMMSE_only, file = paste("sgl_R/modelsNew/","final_", type, "_model", ".RData", sep = ""))

coefMMSE_only = finalFitMMSE_only$beta[,match(min(finalFitMMSE_only$lambdas),finalFitMMSE_only$lambdas)]
absCoefMMSE_only = abs(coefMMSE_only)
coefDfMMSE_only = cbind.data.frame(colnames(x), coefMMSE_only, absCoefMMSE_only, index)
names(coefDfMMSE_only)[names(coefDfMMSE_only)=="colnames(x)"] = "names"
names(coefDfMMSE_only)[names(coefDfMMSE_only)=="index"] = "group"
coefDfMMSE_only <- coefDfMMSE_only[with(coefDfMMSE_only, order(-absCoefMMSE_only)), ]
#coefDf["label"] <- indexDf$group_name_vec[match(coefDf$group, indexDf$index)]


library(ggplot2) 
coefDfNoZeroMMSE_only = subset(coefDfMMSE_only[coefDfMMSE_only$absCoefMMSE_only>0,])
png("sgl_R/plots/classifier_ADNI_Cognition_MMSE_FAQ_only.png", width = 300, height = 225, units='mm',res = 300) 
ggplot(coefDfNoZeroMMSE_only, aes(x=reorder(names, absCoefMMSE_only), weight=absCoefMMSE_only, fill=as.factor(group))) + 
  geom_bar() +
  scale_fill_manual(labels = c("Cognitive Domains from DMs", "MMSE","FAQ"), values = c("dark orange","steelblue","sky blue"), name = "Variable group") +
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

overalImpDf_only2 = cbind.data.frame(c("MMSE", "FAQ","Digital cognitive domains"), c(overall_imp_mmse_sum_only2,overall_imp_faq_sum_only2, overall_imp_cog_sum_only2), c("","",""))
names(overalImpDf_only2) = c("groups", "importance", "feature")


plot1 = ggplot(coefDfNoZero, aes(x=reorder(names, absCoef), weight=absCoef, fill=as.factor(group))) + 
  geom_bar(show.legend = FALSE) +
  scale_fill_manual(labels = c("Digital Cognitive Domains","FAQ"), values = c("dark orange", "sky blue"))+ 
  ylab("absolute coefficients") +
  xlab("Feature Names")+
  coord_flip()
plot2 = ggplot(overalImpDf_only2, aes(fill = groups, y= importance, x = feature)) + geom_bar(position = "stack", stat = "identity", width = 0.1)+scale_fill_manual(values = c("dark orange","sky blue", "steelblue"))+ylab("overall feature importance")
png("sgl_R/plots/classifier_ADNI_CogMMSEFAQ_overallImp_only.png", width = 300, height = 100, units='mm',res = 300) 
plot_grid(plot1, plot2, align = "h", label_size = 12, ncol = 2)
dev.off()


# ##for Cognitive measures and FAQ measures
# adniClassifierCNtoMCICogFAQ = adniClassifierCNtoMCI[,c(cogDomains,faqFeatures,"Class")]
# ##index for Cog domains and FAQ
# index = c(rep(1,9), rep(2,10))
# Class = "Class"
# type = "CogFAQ"
# rocCogFAQ = func.classifier.sgl(adniClassifierCNtoMCICogFAQ, index = index, type = type, "Class")
# boxplot(rocCogFAQ, las = 2, xlab = "", ylab = "AUC", ylim = c(0.8,1),col = "skyblue")
# 
# print(mean(sum(rocCogFAQ)/10))
# maxRocCogFAQ = match(max(rocCogFAQ), rocCogFAQ)

# ##load the model having minimum roc
# load(paste("sgl_R/models/",type, "_model_" ,maxRocCogFAQ, ".RData", sep = ""))
# x = read.csv(paste("sgl_R/data/",type,"_", maxRocCogFAQ, ".csv", sep = ""))
# y= adniClassifierCNtoMCICogFAQ$Class
# 
# x$X = NULL
# finalData = list(x = x, y=y)
# finalFit = SGL(finalData, index = index, type = "logit", maxit = 1000, thresh = 0.001,
#                  min.frac = 0.05, nlam =20, gamma = 0.8, standardize = TRUE,
#                  verbose = FALSE, step = 1, reset = 10, alpha = 0.95, lambdas = fit$fit$lambdas)
# 
# save(finalFit, file = paste("sgl_R/models/","final_", type, "_model","_", maxRocCogFAQ, ".RData", sep = ""))
# 
# coef = finalFit$beta[,match(min(finalFit$lambdas),finalFit$lambdas)]
# absCoef = abs(coef)
# coefDf = cbind.data.frame(colnames(x), coef, absCoef, index)
# names(coefDf)[names(coefDf)=="colnames(x)"] = "names"
# names(coefDf)[names(coefDf)=="index"] = "group"
# coefDf <- coefDf[with(coefDf, order(-absCoef)), ]
# 
# 
# library(ggplot2) 
# coefDfNoZero = subset(coefDf[coefDf$absCoef>0,])
# png("output_figures/classifier_ADNI_Cognition_FAQ.png", width = 250, height = 150, units='mm',res = 300) 
# ggplot(coefDfNoZero, aes(x=reorder(names, absCoef), weight=absCoef, fill=as.factor(group))) + 
#   geom_bar() +
#   scale_fill_manual(labels = c("Cognitive Domains from DMs","FAQ"), values = c( "dark orange","orangered"), name = "Variable group") +
#   ylab("absolute coefficients") +
#   xlab("Feature Names")+
#   coord_flip()
# dev.off()
# 

# ###overall feature importance of FAQ and Cog Domains##
# overall_imp_faq =  coefDfNoZero[coefDfNoZero$names %in% c(grep("FAQ", coefDfNoZero$names, value = TRUE)),]
# overall_imp_faq_sum = colSums(overall_imp_faq['absCoef'])/colSums(coefDfNoZero['absCoef'])
# overall_imp_cog =  coefDfNoZero[coefDfNoZero$names %in% gsub("SA_", "",cogDomains),]
# overall_imp_cog_sum = colSums(overall_imp_cog['absCoef'])/colSums(coefDfNoZero['absCoef'])
# 
# overalImpDf = cbind.data.frame(c("FAQ","Cog from DMs"), c(overall_imp_faq_sum, overall_imp_cog_sum), c("",""))
# names(overalImpDf) = c("groups", "importance", "feature")
# 
# coefDfNoZeroCp = coefDfNoZero
# plot1 = ggplot(coefDfNoZeroCp, aes(x=reorder(names, absCoef), weight=absCoef, fill=as.factor(group))) + 
#   geom_bar(show.legend = FALSE) +
#   scale_fill_manual(labels = c("Cog from DMs",  "FAQ"), values = c("dark orange","orangered"))+ 
#   ylab("absolute coefficients") +
#   xlab("Feature Names")+
#   coord_flip()
# plot2 = ggplot(overalImpDf, aes(fill = groups, y= importance, x = feature)) + geom_bar(position = "stack", stat = "identity", width = 0.1)+scale_fill_manual(values = c("dark orange","orangered", "steelblue"))+ylab("overall feature importance")
# png("output_figures/classifier_ADNI_CogFAQ_overallImp.png", width = 300, height = 100, units='mm',res = 300) 
# plot_grid(plot1, plot2, align = "h",labels = "AUTO", label_size = 12, rel_heights = c(0.1, 1), rel_widths = c(0.5,0.5))
# dev.off()

# ##for digital MMSE and FAQ##
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
targetNodesFeatures = intersect(targetNodesFeatures, colnames(adniClassifierCNtoMCIDigMMSEFAQ))
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
boxplot(rocDigMMSE, las = 2, xlab = "", ylab = "AUC", ylim = c(0.5,1),col = "skyblue")

alphaDigMMSE = paramDigMMSE[[1]]
lambdaDigMMSE = paramDigMMSE[[2]]
#print(mean(sum(rocDigMMSE)/10))
#maxRocDigMMSE = match(max(rocDigMMSE), rocDigMMSE)

##final model
#rm(fit)
##load the model having minimum roc
#load(paste("sgl_R/models/",type, "_model_" ,maxRocDigMMSE, ".RData", sep = ""))
#x = subset(data, select=-c(Class))
# x = subset(adniClassifierCNtoMCIDigMMSE, select=-c(Class))
# #x = read.csv(paste("sgl_R/data/",type,"_", maxRocDigMMSE, ".csv", sep = ""))
# #x$X = NULL
# y= adniClassifierCNtoMCIDigMMSE$Class
# # x = missForest(x, maxiter = 10, ntree = 100, variablewise = FALSE,
# #                decreasing = FALSE, verbose = FALSE,
# #                mtry = floor(sqrt(ncol(x))), replace = TRUE,
# #                classwt = NULL, cutoff = NULL, strata = NULL,
# #                sampsize = NULL, nodesize = NULL, maxnodes = NULL,
# #                xtrue = NA, parallelize = c('no', 'variables', 'forests'))
# # x = x$ximp
# finalData = list(x = x, y=y)
# finalFitDigMMSE = SGL(finalData, index = index, type = "logit", maxit = 1000, thresh = 0.001,
#                         min.frac = 0.1, nlam =2, gamma = 0.8, standardize = TRUE,
#                         verbose = FALSE, step = 1, reset = 10, alpha = alphaDigMMSE, lambdas = c(lambdaDigMMSE,1)
#                         )
# 
# save(finalFitDigMMSE, file = paste("sgl_R/modelsNew/","final_", type, "_model", ".RData", sep = ""))
# 
# coef = finalFitDigMMSE$beta[,match(min(finalFit$lambdas),finalFit$lambdas)]
# absCoef = abs(coef)
# coefDf = cbind.data.frame(colnames(x), coef, absCoef, index)
# names(coefDf)[names(coefDf)=="colnames(x)"] = "names"
# names(coefDf)[names(coefDf)=="index"] = "group"
# coefDf <- coefDf[with(coefDf, order(-absCoef)), ]
# coefDf["label"] <- indexDf$group_name_vec[match(coefDf$group, indexDf$index)]
# 
# library(ggplot2) 
# coefDfNoZero = subset(coefDf[coefDf$absCoef>0,])
# coefDfNoZeroTop15 = coefDfNoZero[1:15,]
# png("sgl_R/plots/classifier_ADNI_DigMMSE.png", width = 300, height = 100, units='mm',res = 300) 
# require(RColorBrewer)
# myColors <- brewer.pal(length(unique(coefDfNoZeroTop15$label))-1, "Oranges")
# ggplot(coefDfNoZeroTop15, aes(x=reorder(names, absCoef), weight=absCoef, fill=as.factor(group))) + 
#   geom_bar() +
#   scale_fill_manual(labels = c("ARGlobalTelemetryVariance","ARObjectPlacementFFT" ), values =myColors, name = "Variable group")+ 
#   ylab("absolute coefficients") +
#   xlab("Feature Names")+
#   coord_flip()
# dev.off()
# 
# ###overall feature importance of Dig and MMSE##
# overall_imp_mmse =  coefDfNoZero[coefDfNoZero$names %in% c(grep("MMSE", coefDfNoZero$names, value = TRUE)),]
# overall_imp_mmse_sum = colSums(overall_imp_mmse['absCoef'])/colSums(coefDfNoZero['absCoef'])
# overall_imp_dig =  coefDfNoZero[coefDfNoZero$names %in% c(grep("AR|BIT|DOT|Motor", coefDfNoZero$names, value = TRUE)),]
# overall_imp_dig_sum = colSums(overall_imp_dig['absCoef'])/colSums(coefDfNoZero['absCoef'])
# 
# overalImpDf = cbind.data.frame(c("MMSE", "DMs"), c(overall_imp_mmse_sum, overall_imp_dig_sum), c("",""))
# names(overalImpDf) = c("groups", "importance", "feature")
# 
# coefDfNoZeroCp = coefDfNoZeroTop15
# coefDfNoZeroCp$group[coefDfNoZeroCp$group < 10] =1
# plot1 = ggplot(coefDfNoZeroCp, aes(x=reorder(names, absCoef), weight=absCoef, fill=as.factor(group))) + 
#   geom_bar(show.legend = FALSE) +
#   scale_fill_manual(labels = c("Cog from DMs",  "FAQ"), values = c("dark orange","orangered"))+ 
#   ylab("absolute coefficients") +
#   xlab("Feature Names")+
#   coord_flip()
# plot2 = ggplot(overalImpDf, aes(fill = groups, y= importance, x = feature)) + geom_bar(position = "stack", stat = "identity", width = 0.1)+scale_fill_manual(values = c("orange", "steelblue"))+ylab("overall feature importance")
# png("sgl_R/plots/classifier_Altoida_DigMMSE_overallImp.png", width = 200, height = 125, units='mm',res = 300) 
# plot_grid(plot1, plot2, align = "h",labels = "AUTO", label_size = 12, ncol=2,nrow=1)
# dev.off()

#for MMSE
dataMMSE = adniClassifierCNtoMCI[,c(mmseFeatures, "Class")]
dataMMSE$MMSE = NULL
index = c(rep(1,5))
indexMMSE = index
Class = "Class"
type = "MMSE"
paramMMSE = func.classifier.sgl.best.param(dataMMSE, index = index, type = type, "Class")
alphaMMSE = paramMMSE[[1]]
lambdaMMSE = paramMMSE[[2]]
x = subset(dataMMSE, select=-c(Class))
y= dataMMSE$Class

finalData = list(x =x, y=y)
finalFitMMSE = SGL(finalData, index = indexMMSE, type = "logit", maxit = 1000, thresh = 0.001,
                         min.frac = 0.1, nlam =2, gamma = 0.8, standardize = TRUE,
                         verbose = FALSE, step = 1, reset = 10, alpha = alphaMMSE, lambdas = c(lambdaMMSE,1)
)

save(finalFitMMSE, file = paste("sgl_R/modelsNew/","final_", type, "_model", ".RData", sep = ""))

coef = finalFitMMSE$beta[,match(min(finalFitMMSE$lambdas),finalFitDig$lambdas)]
absCoef = abs(coef)
coefDf = cbind.data.frame(colnames(x), coef, absCoef, indexMMSE)
names(coefDf)[names(coefDf)=="colnames(x)"] = "names"
names(coefDf)[names(coefDf)=="index"] = "group"
coefDf <- coefDf[with(coefDf, order(-absCoef)), ]


library(ggplot2) 
coefDfNoZero = subset(coefDf[coefDf$absCoef>0,])
png("sgl_R/plots/classifier_ADNI_MMSE.png", width = 200, height = 100, units='mm',res = 300) 
ggplot(coefDfNoZero, aes(x=reorder(names, absCoef), weight=absCoef, fill=as.factor(group))) + 
  geom_bar() +
  scale_fill_manual(labels = "MMSE", values = "steelblue", name = "Variable group") +
  ylab("absolute coefficients") +
  xlab("Feature Names")+theme(axis.text.x = element_text(size=14),
                               axis.text.y = element_text(size=14),
                               axis.title=element_text(size=16))+
  theme(legend.text=element_text(size=12))+
  coord_flip()
dev.off()


#for Cognition
##index, for Cognitive measures 
dataCognition =  adniClassifierCNtoMCI[,c(cogDomains, "Class")]
index = c(rep(1,9))
indexCog = index
Class = "Class"
type = "Cognition"
paramCog = func.classifier.sgl.best.param(dataCognition, index = index, type = type, "Class")
alphaCog = paramCog[[1]]
lambdaCog = paramCog[[2]]
x = subset(dataCognition, select=-c(Class))
y= dataCognition$Class

finalData = list(x =x, y=y)
finalFitCog = SGL(finalData, index = indexCog, type = "logit", maxit = 1000, thresh = 0.001,
                   min.frac = 0.1, nlam =2, gamma = 0.8, standardize = TRUE,
                   verbose = FALSE, step = 1, reset = 10, alpha = alphaCog, lambdas = c(lambdaCog,1)
)

save(finalFitCog, file = paste("sgl_R/modelsNew/","final_", type, "_model", ".RData", sep = ""))

coef = finalFitCog$beta[,match(min(finalFitCog$lambdas),finalFitCog$lambdas)]
absCoef = abs(coef)
coefDf = cbind.data.frame(colnames(x), coef, absCoef, index)
names(coefDf)[names(coefDf)=="colnames(x)"] = "names"
names(coefDf)[names(coefDf)=="index"] = "group"
coefDf <- coefDf[with(coefDf, order(-absCoef)), ]


library(ggplot2) 
coefDfNoZero = subset(coefDf[coefDf$absCoef>0,])
png("sgl_R/plots/classifier_ADNI_Cog.png", width = 200, height = 100, units='mm',res = 300) 
ggplot(coefDfNoZero, aes(x=reorder(names, absCoef), weight=absCoef, fill=as.factor(group))) + 
  geom_bar() +
  scale_fill_manual(labels = "Digital cognitive domains", values = "orange", name = "Variable group") +
  ylab("absolute coefficients") +
  xlab("Feature Names")+theme(axis.text.x = element_text(size=14),
                               axis.text.y = element_text(size=14),
                               axis.title=element_text(size=16))+
  theme(legend.text=element_text(size=12))+
  coord_flip()
dev.off()



#for Digital Measure
##index, for digital measures
dataDig = data[,c(targetNodesFeatures, "Class")]

featureDesc = featureDesc[featureDesc$technical_group_name%in% groups,]
index = c()
##index for groups of digital features

group_name_vec = c()
for(names in colnames(dataDig[,targetNodesFeatures])){
  group_name = featureDesc[featureDesc$feature_name == names, "technical_group_name"]
  group_name = as.character(group_name)
  print(group_name)
  index = c(index,which(groups==group_name))
  group_name_vec = c(group_name_vec, group_name)
}

indexDf = cbind.data.frame(index,group_name_vec)
indexDf = unique(indexDf)
##index for 

##index for 
Class = "Class"
type = "Dig"
indexDig = index
load("~/Documents/Documents_IT/paper/ADNI_VAMBN_paper/sgl_R/resultsFromCluster/adniClassifier/paramDig.RData")
alphaDig = paramDig[[1]]
lambdaDig = paramDig[[2]]

x = subset(dataDig, select=-c(Class))
y= dataDig$Class

finalData = list(x = x, y=y)

finalFitDig = SGL(finalData, index = indexDig, type = "logit", maxit = 1000, thresh = 0.001,
                  min.frac = 0.1, nlam =2, gamma = 0.8, standardize = TRUE,
                  verbose = FALSE, step = 1, reset = 10, alpha =alphaDig, lambdas = c(min(lambdaDig),1))

save(finalFitDig, file = paste("sgl_R/modelsNew/","final_", type, "_model", ".RData", sep = ""))

coef = finalFitDig$beta[,match(min(finalFit$lambdas),finalFit$lambdas)]
absCoef = abs(coef)
coefDf = cbind.data.frame(colnames(x), coef, absCoef, index)
names(coefDf)[names(coefDf)=="colnames(x)"] = "names"
names(coefDf)[names(coefDf)=="index"] = "group"
coefDf <- coefDf[with(coefDf, order(-absCoef)), ]
coefDf["label"] <- indexDf$group_name_vec[match(coefDf$group, indexDf$index)]


library(ggplot2) 
coefDfNoZero = subset(coefDf[coefDf$absCoef>0,])
coefDfNoZero = coefDfNoZero[1:15,]
png("sgl_R/plots/classifier_ADNI_Dig.png", width = 300, height = 200, units='mm',res = 300) 
require(RColorBrewer)
myColors <- brewer.pal(length(unique(coefDfNoZero$label)))
ggplot(coefDfNoZero, aes(x=reorder(names, absCoef), weight=absCoef, fill=as.factor(group))) + 
  geom_bar() +
  scale_fill_manual(labels = c("ARObjectFinding", "MotorTappingFeatures"), values =myColors, name = "Variable group")+ 
  ylab("absolute coefficients") +
  xlab("Feature Names")+theme(axis.text.x = element_text(size=14),
                              axis.text.y = element_text(size=14),
                              axis.title=element_text(size=16))+
  theme(legend.text=element_text(size=12))+
  coord_flip()
dev.off()

##multiple boxplots together
fill <- "#4271AE"
line <- "#1F3552"
boxplotDf = cbind.data.frame(rocFAQ, rocDig, rocCogMMSE, rocDigMMSE)
names(boxplotDf) = c("MMSE", "DMs", "MMSE and Cog", "MMSE and DMs")
boxplotDf = melt(boxplotDf)
plt <- ggplot(data = boxplotDf, aes(x = variable, y = value))
png("output_figures/AUC_boxplot.png", width = 200, height = 155, units='mm',res = 300) 
plt + geom_boxplot(fill = fill, colour = line) + labs(x = "Features", y = "AUC") + ylim(0.8,1)
dev.off()


mmseWithCog = wilcox.test(rocCogMMSE, rocCogMMSEFAQ, paired = TRUE, alternative = "two.sided")
mmseWithDig = wilcox.test(rocMMSE, rocDigMMSE, paired = TRUE, alternative = "two.sided")
