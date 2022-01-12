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
  alphas = seq(0.1,1,0.1)
  alphas = setdiff(alphas,"0.6")
  alphaVec = c()
  
  meanRocVec = c()
  lambdasMin = c()
  #alpha = 0
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

func.classifier.sgl.best.param.single = function(dataDf,index, type, Class){
  xtrain = subset(dataDf, select=-c(Class))
  ytrain = dataDf$Class
  ntrain=length(ytrain)
  set.seed(123)
  train.ext=createFolds(ytrain,k=5,returnTrain=TRUE)
  test.ext=lapply(train.ext,function(x) (1:ntrain)[-x])
  seedVec = c()
  meanRocVec = c()
  lambdasMin = c()
   rocVec = c()
    lambdas = c()
    for(i in 1:5){
      data = list(x=xtrain[train.ext[[i]],],y=ytrain[train.ext[[i]]])
      fit = cvSGL(data, index = index, type = "logit", maxit = 1000, thresh = 0.001,
                  min.frac = 0.05, nlam =20, gamma = 0.8, nfold = 10, standardize = TRUE,
                  verbose = FALSE, step = 1, reset = 10, alpha = 1, lambdas = NULL,
                  foldid = NULL)
      #save(fit, file = paste("sgl_R/modelsNew/",type, "_model","_",alpha, "_", i, ".RData", sep = ""))
      minLambda = min(fit$fit$lambdas)
      trainFit = SGL(data, index = index, type = "logit", maxit = 1000, thresh = 0.001,
                     min.frac = 0.1,nlam=2, gamma = 0.8, standardize = TRUE,
                     verbose = FALSE, step = 1, reset = 10, alpha = 1, lambdas=c(min(fit$fit$lambdas),1))
      yh = predictSGL(trainFit, as.matrix(xtrain[test.ext[[i]],]), 1)
      roc = (roc.curve(yh, weights.class0 = ytrain[test.ext[[i]]], curve = TRUE)$auc)
      print(roc)
      rocVec = c(rocVec, roc)
      lambdas = c(lambdas, minLambda)
    }
    meanRoc = mean(rocVec)
    meanRocVec = c(meanRocVec, meanRoc)
    lambdasMin = c(lambdasMin, min(lambdas))
  maxRoc = match(max(meanRocVec),meanRocVec)
  print(maxRoc)
  bestLambda = lambdasMin[maxRoc]
  return(list(bestLambda, meanRocVec))
}
##index, for Cognitive measures and MMSE measures, FAQ measures
adniClassifierCNtoMCICogMMSEFAQ = adniClassifierCNtoMCI[,c(cogDomains, mmseFeatures,faqFeatures,"Class")]
adniClassifierCNtoMCICogMMSEFAQ$MMSE = NULL
##index for Cog Domains, MMSE and FAQ
index = c(rep(1,9), rep(2,5), rep(3,10))
indexCogMMSEFAQ = index
Class = "Class"
type = "CogMMSEFAQ"
#paramCogMMSEFAQ = func.classifier.sgl.best.param(adniClassifierCNtoMCICogMMSEFAQ, index = indexCogMMSEFAQ, type = type, "Class")
load("~/Documents/Documents_IT/paper/ADNI_VAMBN_paper/sgl_R/resultsFromCluster/adniClassifier/paramCogMMSEFAQ.RData")

alphaCogMMSEFAQ = paramCogMMSEFAQ[[1]]
lambdaCogMMSEFAQ =  paramCogMMSEFAQ[[2]]
x = subset(adniClassifierCNtoMCICogMMSEFAQ, select=-c(Class))
y= adniClassifierCNtoMCICogMMSEFAQ$Class
finalDataCogMMSEFAQ = list(x = x, y=y)
finalFit_CogMMSEFAQ = SGL(finalDataCogMMSEFAQ, index = indexCogMMSEFAQ, type = "logit", maxit = 1000, thresh = 0.001,
               min.frac = 0.1, nlam =2, gamma = 0.8, standardize = TRUE,
               verbose = FALSE, step = 1, reset = 10, alpha = alphaCogMMSEFAQ, lambdas = c(lambdaCogMMSEFAQ,1)
)

save(finalFit_CogMMSEFAQ, file = paste("sgl_R/modelsNew/","final_", type, "_model",".RData", sep = ""))

coef = finalFit_CogMMSEFAQ$beta[,match(min(finalFit_CogMMSEFAQ$lambdas),finalFit_CogMMSEFAQ$lambdas)]
absCoef = abs(coef)
coefDf = cbind.data.frame(colnames(x), coef, absCoef, index)
names(coefDf)[names(coefDf)=="colnames(x)"] = "names"
names(coefDf)[names(coefDf)=="index"] = "group"
coefDf <- coefDf[with(coefDf, order(-absCoef)), ]
coefDf_CogMMSEFAQ = coefDf
coefDf_CogMMSEFAQ = subset(coefDf_CogMMSEFAQ[coefDf_CogMMSEFAQ$absCoef>0,])
##index, for Cognitive measures and MMSE measures only, FAQ measures
adniClassifierCNtoMCICogMMSEFAQ_only = adniClassifierCNtoMCI[,c(cogDomains, "MMSE",faqFeatures,"Class")]

##index for Cog Domains, MMSE and FAQ
index = c(rep(1,9), rep(2,1), rep(3,10))
indexCogMMSEFAQ_only = index
Class = "Class"
type = "CogMMSEFAQ_only"
#paramCogMMSEFAQ_only = func.classifier.sgl.best.param(adniClassifierCNtoMCICogMMSEFAQ_only, index = indexCogMMSEFAQ_only, type = type, "Class")
load("~/Documents/Documents_IT/paper/ADNI_VAMBN_paper/sgl_R/resultsFromCluster/adniClassifier/paramCogMMSEFAQ_only.RData")
alphaCogMMSEFAQ_only = paramCogMMSEFAQ_only[[1]]
lambdaCogMMSEFAQ_only = paramCogMMSEFAQ_only[[2]]
x = subset(adniClassifierCNtoMCICogMMSEFAQ_only, select=-c(Class))
y= adniClassifierCNtoMCICogMMSEFAQ_only$Class

finalDataCogMMSEFAQ_only = list(x = x, y=y)
finalFit_CogMMSEFAQ_only = SGL(finalDataCogMMSEFAQ_only, index = index, type = "logit", maxit = 1000, thresh = 0.001,
                        min.frac = 0.1, nlam =2, gamma = 0.8, standardize = TRUE,
                        verbose = FALSE, step = 1, reset = 10, alpha = alphaCogMMSEFAQ_only, lambdas = c(lambdaCogMMSEFAQ_only,1)
)

save(finalFit_CogMMSEFAQ_only, file = paste("sgl_R/modelsNew/","final_", type, "_model", ".RData", sep = ""))

coefMMSE_only = finalFit_CogMMSEFAQ_only$beta[,match(min(finalFit_CogMMSEFAQ_only$lambdas),finalFit_CogMMSEFAQ_only$lambdas)]
absCoefMMSE_only = abs(coefMMSE_only)
coefDfMMSE_only = cbind.data.frame(colnames(x), coefMMSE_only, absCoefMMSE_only, index)
names(coefDfMMSE_only)[names(coefDfMMSE_only)=="colnames(x)"] = "names"
names(coefDfMMSE_only)[names(coefDfMMSE_only)=="index"] = "group"
coefDfMMSE_only <- coefDfMMSE_only[with(coefDfMMSE_only, order(-absCoefMMSE_only)), ]


library(ggplot2) 
coefDfNoZeroMMSE_only = subset(coefDfMMSE_only[coefDfMMSE_only$absCoefMMSE_only>0,])

###overall feature importance of FAQ, Cog and MMSE##
overall_imp_mmse_only2 =  coefDfNoZeroMMSE_only[coefDfNoZeroMMSE_only$names %in% c(grep("MMSE", coefDfNoZeroMMSE_only$names, value = TRUE)),]
overall_imp_mmse_sum_only2 = colSums(overall_imp_mmse_only2['absCoefMMSE_only'])/colSums(coefDfNoZeroMMSE_only['absCoefMMSE_only'])
overall_imp_faq_only2 =  coefDfNoZeroMMSE_only[coefDfNoZeroMMSE_only$names %in% c(grep("FAQ", coefDfNoZeroMMSE_only$names, value = TRUE)),]
overall_imp_faq_sum_only2 = colSums(overall_imp_faq_only2['absCoefMMSE_only'])/colSums(coefDfNoZeroMMSE_only['absCoefMMSE_only'])
overall_imp_cog_only2 =  coefDfNoZeroMMSE_only[coefDfNoZeroMMSE_only$names %in% gsub("SA_", "",cogDomains),]
overall_imp_cog_sum_only2 = colSums(overall_imp_cog_only2['absCoefMMSE_only'])/colSums(coefDfNoZeroMMSE_only['absCoefMMSE_only'])

overalImpDf_only2 = cbind.data.frame(c("MMSE", "FAQ","Digital Cognitive Domains"), c(overall_imp_mmse_sum_only2,overall_imp_faq_sum_only2, overall_imp_cog_sum_only2), c("","",""))
names(overalImpDf_only2) = c("groups", "importance", "feature")


plot1 = ggplot(coefDf_CogMMSEFAQ, aes(x=reorder(names, absCoef), weight=absCoef, fill=as.factor(group))) + 
  geom_bar(show.legend = FALSE) +
  scale_fill_manual(labels = c("Digital Cognitive Domains","FAQ", "MMSE"), values = c("dark orange","steel blue", "sky blue"))+ 
  ylab("absolute coefficients") +
  xlab("Feature Names")+
  coord_flip()
plot2 = ggplot(overalImpDf_only2, aes(fill = groups, y= importance, x = feature)) + geom_bar(position = "stack", stat = "identity", width = 0.1)+scale_fill_manual(values = c("dark orange","sky blue", "steelblue"))+ylab("overall feature importance")
png("sgl_R/plots/classifier_ADNI_CogMMSEFAQ_overallImp.png", width = 300, height = 100, units='mm',res = 300) 
plot_grid(plot1, plot2, align = "h", label_size = 12, ncol = 2)
dev.off()

##index, for digital tasks and MMSE measures, FAQ measures
adniClassifierCNtoMCIDigMMSEFAQ = adniClassifierCNtoMCI[,c(digFeatures,mmseFeatures,faqFeatures,"Class")]
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

##index for Dig, MMSE and FAQ
indexDf = rbind.data.frame(indexDf, c((max(index)+1), "MMSE"))
indexDf = rbind.data.frame(indexDf, c((max(index)+2), "FAQ"))
index = c(index, rep((max(index)+1),5))
index = c(index, rep((max(index)+1),10))
indexDigMMSEFAQ = index
Class = "Class"
type = "DigMMSEFAQ"

#paramDigMMSE = func.classifier.sgl.best.param(adniClassifierCNtoMCIDigMMSEFAQ, index = indexDigMMSEFAQ, type = "DigMMSEFAQ", "Class")
load("~/Documents/Documents_IT/paper/ADNI_VAMBN_paper/sgl_R/resultsFromCluster/adniClassifier/paramDigFAQMMSE.RData")
alphaDigMMSEFAQ = paramDigFAQMMSE[[1]]
lambdaDigMMSEFAQ = paramDigFAQMMSE[[2]]
x = subset(adniClassifierCNtoMCIDigMMSEFAQ, select=-c(Class))
y= adniClassifierCNtoMCIDigMMSEFAQ$Class

finalData_DigMMSEFAQ = list(x = x, y=y)
finalFitDigMMSEFAQ = SGL(finalData_DigMMSEFAQ, index = indexDigMMSEFAQ, type = "logit", maxit = 1000, thresh = 0.001,
                        min.frac = 0.1, nlam =2, gamma = 0.8, standardize = TRUE,
                        verbose = FALSE, step = 1, reset = 10, alpha = alphaDigMMSEFAQ, lambdas = c(lambdaDigMMSEFAQ,1)
                        )

save(finalFitDigMMSEFAQ, file = paste("sgl_R/modelsNew/","final_", type, "_model", ".RData", sep = ""))

coef = finalFitDigMMSEFAQ$beta[,match(min(finalFitDigMMSEFAQ$lambdas),finalFitDigMMSEFAQ$lambdas)]
absCoef = abs(coef)
coefDf = cbind.data.frame(colnames(x), coef, absCoef, index)
names(coefDf)[names(coefDf)=="colnames(x)"] = "names"
names(coefDf)[names(coefDf)=="index"] = "group"
coefDf <- coefDf[with(coefDf, order(-absCoef)), ]
coefDf["label"] <- indexDf$group_name_vec[match(coefDf$group, indexDf$index)]

coefDfNoZero = subset(coefDf[coefDf$absCoef>0,])
coefDfNoZeroTop15 = coefDfNoZero[1:15,]
coefDfNoZeroTopDigMMSEFAQ = coefDfNoZeroTop15

##index for Dig, MMSE only and FAQ
adniClassifierCNtoMCIDigMMSEFAQ_only = adniClassifierCNtoMCI[,c(digFeatures,"MMSE",faqFeatures,"Class")]
targetNodes = gsub("zcode_", "", targetNodes)
targetNodes = gsub("_VIS1", "", targetNodes)
groups = targetNodes
adniClassifierCNtoMCIDigMMSEFAQ_only = adniClassifierCNtoMCIDigMMSEFAQ_only[,c(targetNodesFeatures, "MMSE", faqFeatures,"Class")]
index = c()

group_name_vec = c()
for(names in colnames(adniClassifierCNtoMCIDigMMSEFAQ_only[,targetNodesFeatures])){
  group_name = featureDesc[featureDesc$feature_name == names, "technical_group_name"]
  group_name = as.character(group_name)
  print(group_name)
  index = c(index,which(groups==group_name))
  group_name_vec = c(group_name_vec, group_name)
}
indexDf = cbind.data.frame(index,group_name_vec)
indexDf = unique(indexDf)

##index for Dig and MMSE only and FAQ
indexDf = rbind.data.frame(indexDf, c((max(index)+1), "MMSE"))
indexDf = rbind.data.frame(indexDf, c((max(index)+2), "FAQ"))
index = c(index, rep((max(index)+1),1))
index = c(index, rep((max(index)+1),10))
indexDigMMSEFAQ_only  = index
Class = "Class"
type = "DigMMSEFAQ_only"
load("~/Documents/Documents_IT/paper/ADNI_VAMBN_paper/sgl_R/resultsFromCluster/adniClassifier/paramDigMMSEFAQ_only.RData")
paramDigMMSEFAQ_only = indexDigMMSEFAQ_only
indexDigMMSEFAQ_only  = index
alphaDigMMSEFAQ_only =  paramDigMMSEFAQ_only[[1]]
lambdaDigMMSEFAQ_only =  paramDigMMSEFAQ_only[[2]]

x = subset(adniClassifierCNtoMCIDigMMSEFAQ_only, select=-c(Class))
y= adniClassifierCNtoMCIDigMMSEFAQ_only$Class

finalData_DigMMSEFAQ_only = list(x = x, y=y)
finalFitDigMMSEFAQ_only = SGL(finalData_DigMMSEFAQ_only, index = indexDigMMSEFAQ_only, type = "logit", maxit = 1000, thresh = 0.001,
                         min.frac = 0.1, nlam =2, gamma = 0.8, standardize = TRUE,
                         verbose = FALSE, step = 1, reset = 10, alpha = alphaDigMMSEFAQ_only, lambdas = c(lambdaDigMMSEFAQ_only,1)
)

save(finalFitDigMMSEFAQ_only, file = paste("sgl_R/modelsNew/","final_", type, "_model", ".RData", sep = ""))

coef = finalFitDigMMSEFAQ_only$beta[,match(min(finalFitDigMMSEFAQ$lambdas),finalFitDigMMSEFAQ$lambdas)]
absCoef = abs(coef)
coefDf = cbind.data.frame(colnames(x), coef, absCoef, index)
names(coefDf)[names(coefDf)=="colnames(x)"] = "names"
names(coefDf)[names(coefDf)=="index"] = "group"
coefDf <- coefDf[with(coefDf, order(-absCoef)), ]
coefDf["label"] <- indexDf$group_name_vec[match(coefDf$group, indexDf$index)]

coefDfNoZero_only = subset(coefDf[coefDf$absCoef>0,])
coefDfNoZeroTop15 = coefDfNoZero_only[1:15,]
coefDfNoZeroTopDigMMSEFAQ_only = coefDfNoZeroTop15


###overall feature importance of FAQ, Dig and MMSE##
overall_imp_mmse_dig_only =  coefDfNoZero_only[coefDfNoZero_only$names %in% c(grep("MMSE", coefDfNoZero_only$names, value = TRUE)),]
overall_imp_mmse_sum_only2 = colSums(overall_imp_mmse_dig_only['absCoef'])/colSums(coefDfNoZero_only['absCoef'])
overall_imp_faq_dig_only =  coefDfNoZero_only[coefDfNoZero_only$names %in% c(grep("FAQ", coefDfNoZero_only$names, value = TRUE)),]
overall_imp_faq_dig_sum_only = colSums(overall_imp_faq_dig_only['absCoef'])/colSums(coefDfNoZero_only['absCoef'])
overall_imp_dig_only =  coefDfNoZero_only[coefDfNoZero_only$names %in% c(grep("AR|BIT|DOT|Motor", coefDfNoZero_only$names, value = TRUE)),]
overall_imp_dig_sum_only2 = colSums(overall_imp_dig_only['absCoef'])/colSums(coefDfNoZero_only['absCoef'])

overalImpDf_only2 = cbind.data.frame(c("MMSE", "FAQ","Digital Tasks"), c(overall_imp_mmse_sum_only2,overall_imp_faq_dig_sum_only, overall_imp_dig_sum_only2), c("","",""))
names(overalImpDf_only2) = c("groups", "importance", "feature")
coefDfNoZeroCp_only = coefDfNoZeroTopDigMMSEFAQ_only
coefDfNoZeroCp_only$group[coefDfNoZeroCp_only$group < unique(coefDfNoZeroCp_only[coefDfNoZeroCp_only$label=="MMSE", "group"])] =1
require(RColorBrewer)
library(RColorBrewer)
myColors <- brewer.pal(length(unique(coefDfNoZeroTopDigMMSEFAQ$label)), "Oranges")
plot1 =   ggplot(coefDfNoZeroTopDigMMSEFAQ, aes(x=reorder(names, absCoef), weight=absCoef, fill = as.factor(group))) + 
  geom_bar() +scale_fill_manual(labels = c("ARObjectPlacementFFT", "MMSE", "FAQ"), values = c(myColors[1],"steelblue", "sky blue"), name = "Variable group")+ 
  ylab("absolute coefficients") +
  xlab("Feature Names")+
  coord_flip()
plot2 = ggplot(overalImpDf_only2, aes(fill = groups, y= importance, x = feature)) + geom_bar(position = "stack", stat = "identity", width = 0.1)+scale_fill_manual(values = c("#FB6A4A", "sky blue", "steelblue"))+ylab("overall feature importance")
png("sgl_R/plots/classifier_ADNI_DigMMSEFAQ_overallImp.png", width = 300, height = 100, units='mm',res = 300) 
plot_grid(plot1, plot2, align = "h", label_size = 12, rel_widths = c(2,1), ncol = 2)
#grid.arrange(plot1, plot2, ncol=2)
dev.off()


#for MMSE
dataMMSE = adniClassifierCNtoMCI[,c(mmseFeatures, "Class")]
dataMMSE$MMSE = NULL
index = c(rep(1,5))
indexMMSE = index
Class = "Class"
type = "MMSE"
paramMMSE = func.classifier.sgl.best.param.single(dataMMSE, index = index, type = type, "Class")
bestLambda = paramMMSE[[1]]
roc = paramMMSE[[2]]
x = subset(dataMMSE, select=-c(Class))
y= dataMMSE$Class

finalData = list(x =x, y=y)
finalFitMMSE = SGL(finalData, index = indexMMSE, type = "logit", maxit = 1000, thresh = 0.001,
                   min.frac = 0.1, nlam =2, gamma = 0.8, standardize = TRUE,
                   verbose = FALSE, step = 1, reset = 10, alpha = 1, lambdas = c(bestLambda,1)
)

save(finalFitMMSE, file = paste("sgl_R/modelsNew/","final_", type, "_model", ".RData", sep = ""))

coef = finalFitMMSE$beta[,match(min(finalFitMMSE$lambdas),finalFitMMSE$lambdas)]
absCoef = abs(coef)
coefDf = cbind.data.frame(colnames(x), coef, absCoef, indexMMSE)
names(coefDf)[names(coefDf)=="colnames(x)"] = "names"
names(coefDf)[names(coefDf)=="indexMMSE"] = "group"
coefDf <- coefDf[with(coefDf, order(-absCoef)), ]


library(ggplot2) 
coefDfNoZeroMMSE = subset(coefDf[coefDf$absCoef>0,])
png("sgl_R/plots/classifier_ADNI_MMSE_upd.png", width = 200, height = 100, units='mm',res = 300) 
ggplot(coefDfNoZeroMMSE, aes(x=reorder(names, absCoef), weight=absCoef, fill=as.factor(group))) + 
  geom_bar() +
  scale_fill_manual(labels = "MMSE", values = "steelblue", name = "Variable group") +
  ylab("absolute coefficients") +
  xlab("Feature Names")+theme(axis.text.x = element_text(size=14),
                              axis.text.y = element_text(size=14),
                              axis.title=element_text(size=16))+
  theme(legend.text=element_text(size=12))+
  coord_flip()
dev.off()

#for FAQ
dataFAQ = adniClassifierCNtoMCI[,c(faqFeatures, "Class")]
index = c(rep(1,10))
indexFAQ = index
Class = "Class"
type = "FAQ"
paramFAQ = func.classifier.sgl.best.param.single(dataFAQ, index = index, type = type, "Class")
lambdaFAQ = paramFAQ[[1]]
roc = paramFAQ[[2]]
x = subset(dataFAQ, select=-c(Class))
y= dataFAQ$Class

finalData = list(x =x, y=y)
finalFitFAQ = SGL(finalData, index = indexFAQ, type = "logit", maxit = 1000, thresh = 0.001,
                   min.frac = 0.1, nlam =2, gamma = 0.8, standardize = TRUE,
                   verbose = FALSE, step = 1, reset = 10, alpha = 1, lambdas = c(lambdaFAQ,1)
)

save(finalFitFAQ, file = paste("sgl_R/modelsNew/","final_", type, "_model", ".RData", sep = ""))

coef = finalFitFAQ$beta[,match(min(finalFitFAQ$lambdas),finalFitFAQ$lambdas)]
absCoef = abs(coef)
coefDf = cbind.data.frame(colnames(x), coef, absCoef, indexFAQ)
names(coefDf)[names(coefDf)=="colnames(x)"] = "names"
names(coefDf)[names(coefDf)=="indexFAQ"] = "group"
coefDf <- coefDf[with(coefDf, order(-absCoef)), ]


library(ggplot2) 
coefDfNoZeroFAQ = subset(coefDf[coefDf$absCoef>0,])
png("sgl_R/plots/classifier_ADNI_FAQ_upd.png", width = 200, height = 100, units='mm',res = 300) 
ggplot(coefDfNoZeroFAQ, aes(x=reorder(names, absCoef), weight=absCoef, fill=as.factor(group))) + 
  geom_bar() +
  scale_fill_manual(labels = "FAQ", values = "sky blue", name = "Variable group") +
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
paramCog = func.classifier.sgl.best.param.single(dataCognition, index = indexCog, type = type, "Class")
lambdaCog = paramCog[[1]]
rocCog = paramCog[[2]]
x = subset(dataCognition, select=-c(Class))
y= dataCognition$Class

finalData = list(x =x, y=y)
finalFitCog = SGL(finalData, index = indexCog, type = "logit", maxit = 1000, thresh = 0.001,
                  min.frac = 0.1, nlam =2, gamma = 0.8, standardize = TRUE,
                  verbose = FALSE, step = 1, reset = 10, alpha = 1, lambdas = c(lambdaCog,1)
)

save(finalFitCog, file = paste("sgl_R/modelsNew/","final_", type, "_model", ".RData", sep = ""))

coef = finalFitCog$beta[,match(min(finalFitCog$lambdas),finalFitCog$lambdas)]
absCoef = abs(coef)
coefDf = cbind.data.frame(colnames(x), coef, absCoef, indexCog)
names(coefDf)[names(coefDf)=="colnames(x)"] = "names"
names(coefDf)[names(coefDf)=="indexCog"] = "group"

coefDfCog <- coefDf[with(coefDf, order(-absCoef)), ]


library(ggplot2) 
coefDfCogNoZero = subset(coefDfCog[coefDfCog$absCoef>0,])
png("sgl_R/plots/classifier_ADNI_Cog_upd.png", width = 200, height = 100, units='mm',res = 300) 
ggplot(coefDfCogNoZero, aes(x=reorder(names, absCoef), weight=absCoef, fill=as.factor(group))) + 
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
dataDig = adniClassifierCNtoMCI[,c(targetNodesFeatures, "Class")]

featureDesc = featureDesc[featureDesc$technical_group_name%in% groups,]
index = c()
##index for groups of digital features
indexDig = index
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
paramDig = func.classifier.sgl.best.param(dataDig, index = indexDig, type = type, "Class")

#load("~/Documents/Documents_IT/paper/ADNI_VAMBN_paper/sgl_R/resultsFromCluster/adniClassifier/paramDig.RData")
alphaDig = paramDig[[1]]
lambdaDig = paramDig[[2]]

x = subset(dataDig, select=-c(Class))
y= dataDig$Class

finalData = list(x = x, y=y)

finalFitDig = SGL(finalData, index = indexDig, type = "logit", maxit = 1000, thresh = 0.001,
                  min.frac = 0.1, nlam =2, gamma = 0.8, standardize = TRUE,
                  verbose = FALSE, step = 1, reset = 10, alpha =alphaDig, lambdas = c(min(lambdaDig),1))

save(finalFitDig, file = paste("sgl_R/modelsNew/","final_", type, "_model", ".RData", sep = ""))

coef = finalFitDig$beta[,match(min(finalFitDig$lambdas),finalFitDig$lambdas)]
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
myColors <- brewer.pal(length(unique(coefDfNoZero$label)), "Oranges")
ggplot(coefDfNoZero, aes(x=reorder(names, absCoef), weight=absCoef, fill=as.factor(group))) + 
  geom_bar() +
  scale_fill_manual(labels = c("ARObjectFinding", "MotorTestDurations", "MotorTappingFeatures"), values =myColors, name = "Variable group")+ 
  ylab("absolute coefficients") +
  xlab("Feature Names")+theme(axis.text.x = element_text(size=14),
                              axis.text.y = element_text(size=14),
                              axis.title=element_text(size=16))+
  theme(legend.text=element_text(size=12))+
  coord_flip()
dev.off()





############################################################################################################
##index, for Cognitive measures and MMSE measures
adniClassifierCNtoMCICogMMSE = adniClassifierCNtoMCI[,c(cogDomains, mmseFeatures,"Class")]
adniClassifierCNtoMCICogMMSE$MMSE = NULL
##index for Cog Domains, MMSE and FAQ
index = c(rep(1,9), rep(2,5))
indexCogMMSE = index
Class = "Class"
type = "CogMMSE"
paramCogMMSE = func.classifier.sgl.best.param(adniClassifierCNtoMCICogMMSE, index = indexCogMMSE, type = type, "Class")
saveRDS(paramCogMMSE,"~/Documents/Documents_IT/paper/ADNI_VAMBN_paper/sgl_R/resultsFromCluster/adniClassifier/paramCogMMSE.rds" )

#load("~/Documents/Documents_IT/paper/ADNI_VAMBN_paper/sgl_R/resultsFromCluster/adniClassifier/paramCogMMSEFAQ.RData")

alphaCogMMSE = paramCogMMSE[[1]]
lambdaCogMMSE =  paramCogMMSE[[2]]
x = subset(adniClassifierCNtoMCICogMMSE, select=-c(Class))
y= adniClassifierCNtoMCICogMMSE$Class
finalDataCogMMSE = list(x = x, y=y)
finalFit_CogMMSE = SGL(finalDataCogMMSE, index = indexCogMMSE, type = "logit", maxit = 1000, thresh = 0.001,
                          min.frac = 0.1, nlam =2, gamma = 0.8, standardize = TRUE,
                          verbose = FALSE, step = 1, reset = 10, alpha = alphaCogMMSE, lambdas = c(lambdaCogMMSE,1)
)

save(finalFit_CogMMSE, file = paste("sgl_R/modelsNew/","final_", type, "_model",".RData", sep = ""))

coef = finalFit_CogMMSE$beta[,match(min(finalFit_CogMMSE$lambdas),finalFit_CogMMSE$lambdas)]
absCoef = abs(coef)
coefDf = cbind.data.frame(colnames(x), coef, absCoef, index)
names(coefDf)[names(coefDf)=="colnames(x)"] = "names"
names(coefDf)[names(coefDf)=="index"] = "group"
coefDf <- coefDf[with(coefDf, order(-absCoef)), ]
coefDf_CogMMSE = coefDf
coefDf_CogMMSE = subset(coefDf_CogMMSE[coefDf_CogMMSE$absCoef>0,])
##index, for Cognitive measures and MMSE measures only, FAQ measures
adniClassifierCNtoMCICogMMSE_only = adniClassifierCNtoMCI[,c(cogDomains, "MMSE","Class")]

##index for Cog Domains, MMSE only
index = c(rep(1,9), rep(2,1))
indexCogMMSE_only = index
Class = "Class"
type = "CogMMSE_only"
paramCogMMSE_only = func.classifier.sgl.best.param(adniClassifierCNtoMCICogMMSE_only, index = indexCogMMSE_only, type = type, "Class")
saveRDS(paramCogMMSE_only,"~/Documents/Documents_IT/paper/ADNI_VAMBN_paper/sgl_R/resultsFromCluster/adniClassifier/paramCogMMSE_only.rds" )
#load("~/Documents/Documents_IT/paper/ADNI_VAMBN_paper/sgl_R/resultsFromCluster/adniClassifier/paramCogMMSEFAQ_only.RData")
alphaCogMMSE_only = paramCogMMSE_only[[1]]
lambdaCogMMSE_only = paramCogMMSE_only[[2]]
x = subset(adniClassifierCNtoMCICogMMSE_only, select=-c(Class))
y= adniClassifierCNtoMCICogMMSE_only$Class

finalDataCogMMSE_only = list(x = x, y=y)
finalFit_CogMMSE_only = SGL(finalDataCogMMSE_only, index = index, type = "logit", maxit = 1000, thresh = 0.001,
                               min.frac = 0.1, nlam =2, gamma = 0.8, standardize = TRUE,
                               verbose = FALSE, step = 1, reset = 10, alpha = alphaCogMMSE_only, lambdas = c(lambdaCogMMSE_only,1)
)

save(finalFit_CogMMSE_only, file = paste("sgl_R/modelsNew/","final_", type, "_model", ".RData", sep = ""))

coefMMSE_only = finalFit_CogMMSE_only$beta[,match(min(finalFit_CogMMSE_only$lambdas),finalFit_CogMMSE_only$lambdas)]
absCoefMMSE_only = abs(coefMMSE_only)
coefDfMMSE_only = cbind.data.frame(colnames(x), coefMMSE_only, absCoefMMSE_only, index)
names(coefDfMMSE_only)[names(coefDfMMSE_only)=="colnames(x)"] = "names"
names(coefDfMMSE_only)[names(coefDfMMSE_only)=="index"] = "group"
coefDfMMSE_only <- coefDfMMSE_only[with(coefDfMMSE_only, order(-absCoefMMSE_only)), ]


library(ggplot2) 
coefDfNoZeroMMSE_only = subset(coefDfMMSE_only[coefDfMMSE_only$absCoefMMSE_only>0,])

###overall feature importance of FAQ, Cog and MMSE##
overall_imp_mmse_only2 =  coefDfNoZeroMMSE_only[coefDfNoZeroMMSE_only$names %in% c(grep("MMSE", coefDfNoZeroMMSE_only$names, value = TRUE)),]
overall_imp_mmse_sum_only2 = colSums(overall_imp_mmse_only2['absCoefMMSE_only'])/colSums(coefDfNoZeroMMSE_only['absCoefMMSE_only'])
overall_imp_cog_only2 =  coefDfNoZeroMMSE_only[coefDfNoZeroMMSE_only$names %in% gsub("SA_", "",cogDomains),]
overall_imp_cog_sum_only2 = colSums(overall_imp_cog_only2['absCoefMMSE_only'])/colSums(coefDfNoZeroMMSE_only['absCoefMMSE_only'])

overalImpDf_only2 = cbind.data.frame(c("MMSE","Digital Cognitive Domains"), c(overall_imp_mmse_sum_only2, overall_imp_cog_sum_only2), c("",""))
names(overalImpDf_only2) = c("groups", "importance", "feature")


plot1 = ggplot(coefDf_CogMMSE, aes(x=reorder(names, absCoef), weight=absCoef, fill=as.factor(group))) + 
  geom_bar(show.legend = FALSE) +
  scale_fill_manual(labels = c("Digital Cognitive Domains","MMSE"), values = c("dark orange", "steel blue"))+ 
  ylab("absolute coefficients") +
  xlab("Feature Names")+
  coord_flip()
plot2 = ggplot(overalImpDf_only2, aes(fill = groups, y= importance, x = feature)) + geom_bar(position = "stack", stat = "identity", width = 0.1)+scale_fill_manual(values = c("dark orange", "steelblue"))+ylab("overall feature importance")
png("sgl_R/plots/classifier_ADNI_CogMMSE_overallImp.png", width = 300, height = 100, units='mm',res = 300) 
plot_grid(plot1, plot2, align = "h", label_size = 12, ncol = 2)
dev.off()

##index, for FAQ and MMSE measures
adniClassifierCNtoMCIFAQMMSE = adniClassifierCNtoMCI[,c(faqFeatures, mmseFeatures,"Class")]
adniClassifierCNtoMCIFAQMMSE$MMSE = NULL
##index for Cog Domains, MMSE and FAQ
index = c(rep(1,10), rep(2,5))
indexFAQMMSE = index
Class = "Class"
type = "FAQMMSE"
paramFAQMMSE = func.classifier.sgl.best.param(adniClassifierCNtoMCIFAQMMSE, index = indexFAQMMSE, type = type, "Class")
saveRDS(paramFAQMMSE,"~/Documents/Documents_IT/paper/ADNI_VAMBN_paper/sgl_R/resultsFromCluster/adniClassifier/paramFAQMMSE.rds" )

#load("~/Documents/Documents_IT/paper/ADNI_VAMBN_paper/sgl_R/resultsFromCluster/adniClassifier/paramCogMMSEFAQ.RData")

alphaFAQMMSE = paramFAQMMSE[[1]]
lambdaFAQMMSE =  paramFAQMMSE[[2]]
x = subset(adniClassifierCNtoMCIFAQMMSE, select=-c(Class))
y= adniClassifierCNtoMCIFAQMMSE$Class
finalDataFAQMMSE = list(x = x, y=y)
finalFit_FAQMMSE = SGL(finalDataFAQMMSE, index = indexFAQMMSE, type = "logit", maxit = 1000, thresh = 0.001,
                       min.frac = 0.1, nlam =2, gamma = 0.8, standardize = TRUE,
                       verbose = FALSE, step = 1, reset = 10, alpha = alphaFAQMMSE, lambdas = c(lambdaFAQMMSE,1)
)

save(finalFit_FAQMMSE, file = paste("sgl_R/modelsNew/","final_", type, "_model",".RData", sep = ""))

coef = finalFit_FAQMMSE$beta[,match(min(finalFit_FAQMMSE$lambdas),finalFit_FAQMMSE$lambdas)]
absCoef = abs(coef)
coefDf = cbind.data.frame(colnames(x), coef, absCoef, index)
names(coefDf)[names(coefDf)=="colnames(x)"] = "names"
names(coefDf)[names(coefDf)=="index"] = "group"
coefDf <- coefDf[with(coefDf, order(-absCoef)), ]
coefDf_FAQMMSE = coefDf
coefDf_FAQMMSE = subset(coefDf_FAQMMSE[coefDf_FAQMMSE$absCoef>0,])
##index, for Cognitive measures and MMSE measures only, FAQ measures
adniClassifierCNtoMCIFAQMMSE_only = adniClassifierCNtoMCI[,c(faqFeatures, "MMSE","Class")]

##index for Cog Domains, MMSE only
index = c(rep(1,10), rep(2,1))
indexFAQMMSE_only = index
Class = "Class"
type = "FAQMMSE_only"
paramFAQMMSE_only = func.classifier.sgl.best.param(adniClassifierCNtoMCIFAQMMSE_only, index = indexFAQMMSE_only, type = type, "Class")
saveRDS(paramFAQMMSE_only,"~/Documents/Documents_IT/paper/ADNI_VAMBN_paper/sgl_R/resultsFromCluster/adniClassifier/paramFAQMMSE_only.rds")
#load("~/Documents/Documents_IT/paper/ADNI_VAMBN_paper/sgl_R/resultsFromCluster/adniClassifier/paramCogMMSEFAQ_only.RData")
alphaFAQMMSE_only = paramFAQMMSE_only[[1]]
lambdaFAQMMSE_only = paramFAQMMSE_only[[2]]
x = subset(adniClassifierCNtoMCIFAQMMSE_only, select=-c(Class))
y= adniClassifierCNtoMCIFAQMMSE_only$Class

finalDataFAQMMSE_only = list(x = x, y=y)
finalFit_FAQMMSE_only = SGL(finalDataFAQMMSE_only, index = index, type = "logit", maxit = 1000, thresh = 0.001,
                            min.frac = 0.1, nlam =2, gamma = 0.8, standardize = TRUE,
                            verbose = FALSE, step = 1, reset = 10, alpha = alphaFAQMMSE_only, lambdas = c(lambdaFAQMMSE_only,1)
)

save(finalFit_FAQMMSE_only, file = paste("sgl_R/modelsNew/","final_", type, "_model", ".RData", sep = ""))

coefMMSE_only = finalFit_FAQMMSE_only$beta[,match(min(finalFit_FAQMMSE_only$lambdas),finalFit_FAQMMSE_only$lambdas)]
absCoefMMSE_only = abs(coefMMSE_only)
coefDfMMSE_only = cbind.data.frame(colnames(x), coefMMSE_only, absCoefMMSE_only, index)
names(coefDfMMSE_only)[names(coefDfMMSE_only)=="colnames(x)"] = "names"
names(coefDfMMSE_only)[names(coefDfMMSE_only)=="index"] = "group"
coefDfMMSE_only <- coefDfMMSE_only[with(coefDfMMSE_only, order(-absCoefMMSE_only)), ]


library(ggplot2) 
coefDfNoZeroMMSE_only = subset(coefDfMMSE_only[coefDfMMSE_only$absCoefMMSE_only>0,])

###overall feature importance of FAQ, Cog and MMSE##
overall_imp_mmse_only2 =  coefDfNoZeroMMSE_only[coefDfNoZeroMMSE_only$names %in% c(grep("MMSE", coefDfNoZeroMMSE_only$names, value = TRUE)),]
overall_imp_mmse_sum_only2 = colSums(overall_imp_mmse_only2['absCoefMMSE_only'])/colSums(coefDfNoZeroMMSE_only['absCoefMMSE_only'])
overall_imp_FAQ_only2 =  coefDfNoZeroMMSE_only[coefDfNoZeroMMSE_only$names %in% faqFeatures,]
overall_imp_FAQ_sum_only2 = colSums(overall_imp_FAQ_only2['absCoefMMSE_only'])/colSums(coefDfNoZeroMMSE_only['absCoefMMSE_only'])

overalImpDf_only2 = cbind.data.frame(c("MMSE", "FAQ"), c(overall_imp_mmse_sum_only2, overall_imp_FAQ_sum_only2), c("",""))
names(overalImpDf_only2) = c("groups", "importance", "feature")


plot1 = ggplot(coefDf_FAQMMSE, aes(x=reorder(names, absCoef), weight=absCoef, fill=as.factor(group))) + 
  geom_bar(show.legend = FALSE) +
  scale_fill_manual(labels = c("FAQ","MMSE"), values = c("sky blue", "steel blue"))+ 
  ylab("absolute coefficients") +
  xlab("Feature Names")+
  coord_flip()
plot2 = ggplot(overalImpDf_only2, aes(fill = groups, y= importance, x = feature)) + geom_bar(position = "stack", stat = "identity", width = 0.1)+scale_fill_manual(values = c("sky blue", "steelblue"))+ylab("overall feature importance")
png("sgl_R/plots/classifier_ADNI_FAQMMSE_overallImp.png", width = 300, height = 100, units='mm',res = 300) 
plot_grid(plot1, plot2, align = "h", label_size = 12, ncol = 2)
dev.off()


save.image("~/Documents/Documents_IT/paper/ADNI_VAMBN_paper_final/resultsSGLCNMCI_upd.RData")
