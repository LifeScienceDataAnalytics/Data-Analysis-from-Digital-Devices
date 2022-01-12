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
adniClassifierMCItoDementia = read.csv("data_DigMCItoDementia.csv")
##remove missing values in FAQ columns##
adniClassifierMCItoDementia = na.omit(adniClassifierMCItoDementia)
#mmseFeatures = grep("MMSE", colnames(adniClassifierCNtoMCI), value = TRUE)
faqFeatures = grep("FAQ", colnames(adniClassifierMCItoDementia), value = TRUE)
digFeatures = grep("AR|^Motor|BIT|DOT", colnames(adniClassifierMCItoDementia), value = TRUE)
adniClassifierMCItoDementia$X = NULL
#dataDf = adniClassifierMCItoDementiaDigFAQ
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
    #alpha = 0.1
    rocVec = c()
    lambdas = c()
    for(i in 1:5){
      #i = 1
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

##index, for Cognitive measures FAQ measures
adniClassifierMCItoDemCogFAQ = adniClassifierMCItoDementia[,c(cogDomains,faqFeatures,"Class")]

##index for Cog Domains and FAQ
index = c(rep(1,9),  rep(2,10))
indexCogFAQ = index
Class = "Class"
type = "CogFAQ"
#paramCogMMSEFAQ = func.classifier.sgl.best.param(adniClassifierCNtoMCICogMMSEFAQ, index = indexCogMMSEFAQ, type = type, "Class")
load("~/Documents/Documents_IT/paper/ADNI_VAMBN_paper/sgl_R/resultsFromCluster/adniClassifier/paramCogFAQMCIDem.RData")

alphaCogFAQ = paramCogFAQ[[1]]
lambdaCogFAQ =  paramCogFAQ[[2]]
x = subset(adniClassifierMCItoDemCogFAQ, select=-c(Class))
y= adniClassifierMCItoDemCogFAQ$Class
finalDataCogFAQ = list(x = x, y=y)
finalFit_CogFAQ = SGL(finalDataCogFAQ, index = indexCogFAQ, type = "logit", maxit = 1000, thresh = 0.001,
                       min.frac = 0.1, nlam =2, gamma = 0.8, standardize = TRUE,
                       verbose = FALSE, step = 1, reset = 10, alpha = alphaCogFAQ, lambdas = c(lambdaCogFAQ,1)
)

save(finalFit_CogFAQ, file = paste("sgl_R/modelsNew/","final_", type, "_model",".RData", sep = ""))

coef = finalFit_CogFAQ$beta[,match(min(finalFit_CogFAQ$lambdas),finalFit_CogFAQ$lambdas)]
absCoef = abs(coef)
coefDf = cbind.data.frame(colnames(x), coef, absCoef, index)
names(coefDf)[names(coefDf)=="colnames(x)"] = "names"
names(coefDf)[names(coefDf)=="index"] = "group"
coefDf <- coefDf[with(coefDf, order(-absCoef)), ]
coefDf_CogFAQ = coefDf
coefDf_CogFAQ = subset(coefDf_CogFAQ[coefDf_CogFAQ$absCoef>0,])

###overall feature importance of FAQ, Cog ##
overall_imp_faq =  coefDf_CogFAQ[coefDf_CogFAQ$names %in% c(grep("FAQ", coefDf_CogFAQ$names, value = TRUE)),]
overall_imp_faq_sum = colSums(overall_imp_faq['absCoef'])/colSums(coefDf_CogFAQ['absCoef'])
overall_imp_cog =  coefDf_CogFAQ[coefDf_CogFAQ$names %in% gsub("SA_", "",cogDomains),]
overall_imp_cog_sum = colSums(overall_imp_cog['absCoef'])/colSums(coefDf_CogFAQ['absCoef'])

overalImpDf= cbind.data.frame(c("FAQ","Digital Cognitive Domains"), c(overall_imp_faq_sum,overall_imp_cog_sum), c("",""))
names(overalImpDf) = c("groups", "importance", "feature")


plot1 = ggplot(coefDf_CogFAQ, aes(x=reorder(names, absCoef), weight=absCoef, fill=as.factor(group))) + 
  geom_bar(show.legend = FALSE) +
  scale_fill_manual(labels = c("Digital Cognitive Domains","FAQ"), values = c("dark orange", "sky blue"))+ 
  ylab("absolute coefficients") +
  xlab("Feature Names")+
  coord_flip()
plot2 = ggplot(overalImpDf, aes(fill = groups, y= importance, x = feature)) + geom_bar(position = "stack", stat = "identity", width = 0.1)+scale_fill_manual(values = c("dark orange","sky blue"))+ylab("overall feature importance")
png("sgl_R/plots/classifier_ADNI_CogFAQMCIDem_overallImp.png", width = 300, height = 100, units='mm',res = 300) 
plot_grid(plot1, plot2, align = "h", label_size = 12, ncol = 2)
dev.off()


##index, for Digital measures FAQ measures
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
load("~/Documents/Documents_IT/paper/ADNI_VAMBN_paper/sgl_R/resultsFromCluster/adniClassifier/paramDigFAQMCIDem.RData")
paramDigFAQMCIDem = paramDigFAQ
alphaDigFAQ = paramDigFAQMCIDem[[1]]
lambdaDigFAQ= paramDigFAQMCIDem[[2]]
rocDigFAQ = paramDigFAQMCIDem[[3]]
type = "DigFAQ"

x = subset(adniClassifierMCItoDementiaDigFAQ, select=-c(Class))
y= adniClassifierMCItoDementiaDigFAQ$Class
finalData_digFAQ = list(x = x, y=y)
finalFit_digFAQ = SGL(finalData_digFAQ, index = indexDigFAQ, type = "logit", maxit = 1000, thresh = 0.001,
               min.frac = 0.1, nlam =2, gamma = 0.8,standardize = TRUE,
               verbose = FALSE, step = 1, reset = 10, alpha = alphaDigFAQ, lambdas = c(lambdaDigFAQ,1)
)

save(finalFit_digFAQ, file = paste("sgl_R/modelsNew/","final_", type, "_model", ".RData", sep = ""))
coefDigFAQ= finalFit_digFAQ$beta[,match(min(finalFit_digFAQ$lambdas),finalFit_digFAQ$lambdas)]
absCoefDigFAQ= abs(coefDigFAQ)
coefDfDigFAQ = cbind.data.frame(colnames(x), absCoefDigFAQ, indexDigFAQ)
names(coefDfDigFAQ)[names(coefDfDigFAQ)=="colnames(x)"] = "names"
names(coefDfDigFAQ)[names(coefDfDigFAQ)=="indexDigFAQ"] = "group"
coefDfDigFAQ <- coefDfDigFAQ[with(coefDfDigFAQ, order(-absCoefDigFAQ)), ]
coefDfDigFAQ["label"] <- indexDf$group_name_vec[match(coefDfDigFAQ$group, indexDf$index)]

coefDfDigFAQ <- coefDfDigFAQ[with(coefDfDigFAQ, order(-absCoefDigFAQ)), ]


library(ggplot2) 
coefDfDigFAQ = subset(coefDfDigFAQ[coefDfDigFAQ$absCoefDigFAQ>0,])
coefDfDigFAQTop15 = coefDfDigFAQ[1:15,]


###overall feature importance of FAQ and Cog Domains##
overall_imp_faq =  coefDfDigFAQ[coefDfDigFAQ$names %in% c(grep("FAQ", coefDfDigFAQ$names, value = TRUE)),]
overall_imp_faq_sum = colSums(overall_imp_faq['absCoefDigFAQ'])/colSums(coefDfDigFAQ['absCoefDigFAQ'])
round(overall_imp_faq_sum,2)
overall_imp_dig = coefDfDigFAQ[coefDfDigFAQ$names %in% c(grep("AR|BIT|DOT|Motor", coefDfDigFAQ$names, value = TRUE)),]
overall_imp_dig_sum = colSums(overall_imp_dig['absCoefDigFAQ'])/colSums(coefDfDigFAQ['absCoefDigFAQ'])
round(overall_imp_dig_sum,2)
overalImpDf = cbind.data.frame(c("FAQ", "Digital tasks"), c(overall_imp_faq_sum, overall_imp_dig_sum), c("",""))
names(overalImpDf) = c("groups", "importance", "feature")

##plot for main paper
coefDfDigFAQTop15Cp = coefDfDigFAQTop15
#coefDfDigFAQTop15Cp$group[coefDfDigFAQTop15Cp$group < unique(coefDfDigFAQ[coefDfDigFAQ$label=="MMSE", "group"])] =1
require(RColorBrewer)
library(RColorBrewer)
myColors <- brewer.pal(length(unique(coefDfDigFAQTop15$label)), "Oranges")
plot1 =   ggplot(coefDfDigFAQTop15, aes(x=reorder(names, absCoefDigFAQ), weight=absCoefDigFAQ, fill = as.factor(group))) + 
  geom_bar() +scale_fill_manual(labels = c("BITDOTMotorInstructionReadingTimeRatios", "ARGlobalTelemetryVariance","MotorTappingFeatures", "FAQ"), values = c(myColors[1:3], "sky blue"), name = "Variable group")+ 
  ylab("absolute coefficients") +
  xlab("Feature Names")+
  coord_flip()
plot2 = ggplot(overalImpDf, aes(fill = groups, y= importance, x = feature)) + geom_bar(position = "stack", stat = "identity", width = 0.1)+scale_fill_manual(values = c("#FB6A4A", "sky blue"))+ylab("overall feature importance")
png("sgl_R/plots/classifier_ADNI_DigFAQMCIDem.png", width = 300, height = 100, units='mm',res = 300) 
plot_grid(plot1, plot2, align = "h", label_size = 12, rel_widths = c(2,1), ncol = 2)
#grid.arrange(plot1, plot2, ncol=2)
dev.off()

#for FAQ
dataFAQ = adniClassifierMCItoDementia[,c(faqFeatures, "Class")]
index = c(rep(1,10))
indexFAQ = index
Class = "Class"
type = "FAQ"
paramFAQ = func.classifier.sgl.best.param.single(dataFAQ, index = indexFAQ, type = type, "Class")
lambdaFAQ = paramFAQ[[1]]
roc = paramFAQ[[2]]
x = subset(dataFAQ, select=-c(Class))
y= dataFAQ$Class

finalData = list(x =x, y=y)
finalFitFAQ = SGL(finalData, index = indexFAQ, type = "logit", maxit = 1000, thresh = 0.001,
                  min.frac = 0.1, nlam =2, gamma = 0.8, standardize = TRUE,
                  verbose = FALSE, step = 1, reset = 10, alpha = 1, lambdas = c(lambdaFAQ,1)
)

save(finalFitFAQ, file = paste("sgl_R/modelsNew/","final_mci_dem_", type, "_model", ".RData", sep = ""))

coef = finalFitFAQ$beta[,match(min(finalFitFAQ$lambdas),finalFitFAQ$lambdas)]
absCoef = abs(coef)
coefDf = cbind.data.frame(colnames(x), coef, absCoef, indexFAQ)
names(coefDf)[names(coefDf)=="colnames(x)"] = "names"
names(coefDf)[names(coefDf)=="indexFAQ"] = "group"
coefDf <- coefDf[with(coefDf, order(-absCoef)), ]


library(ggplot2) 
coefDfNoZeroFAQ = subset(coefDf[coefDf$absCoef>0,])
png("sgl_R/plots/classifier_ADNI_FAQ_MCIDEM_upd.png", width = 200, height = 100, units='mm',res = 300) 
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
dataCognition =  adniClassifierMCItoDementia[,c(cogDomains, "Class")]
index = c(rep(1,9))
indexCog = index
Class = "Class"
type = "Cognition"
paramCog = func.classifier.sgl.best.param.single(dataCognition, index = indexCog, type = type, "Class")
lambdaCog = paramCog[[1]]
roc = paramCog[[2]]
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
coefDf = cbind.data.frame(colnames(x), coef, absCoef, index)
names(coefDf)[names(coefDf)=="colnames(x)"] = "names"
names(coefDf)[names(coefDf)=="index"] = "group"
coefDfCog <- coefDf[with(coefDf, order(-absCoef)), ]


library(ggplot2) 
coefDfCogNoZero = subset(coefDfCog[coefDfCog$absCoef>0,])
png("sgl_R/plots/classifier_ADNI_Cog_mciDem_upd.png", width = 200, height = 100, units='mm',res = 300) 
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
dataDig = adniClassifierMCItoDementia[,c(targetNodesFeatures, "Class")]

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

save(finalFitDig, file = paste("sgl_R/modelsNew/","final_mciDem_", type, "_model", ".RData", sep = ""))

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
png("sgl_R/plots/classifier_ADNI_Dig_mciDem_upd.png", width = 300, height = 200, units='mm',res = 300) 
require(RColorBrewer)
myColors <- brewer.pal(length(unique(coefDfNoZero$label)), "Oranges")
ggplot(coefDfNoZero, aes(x=reorder(names, absCoef), weight=absCoef, fill=as.factor(group))) + 
  geom_bar() +
  scale_fill_manual(labels = c("ARObjectFinding", "MotorTappingFeatures"), values =myColors[1:2], name = "Variable group")+ 
  ylab("absolute coefficients") +
  xlab("Feature Names")+theme(axis.text.x = element_text(size=14),
                              axis.text.y = element_text(size=14),
                              axis.title=element_text(size=16))+
  theme(legend.text=element_text(size=12))+
  coord_flip()
dev.off()

#save.image("~/Documents/Documents_IT/paper/ADNI_VAMBN_paper/resultsSGLMCIDem_upd.RData")

save.image("~/Documents/Documents_IT/paper/ADNI_VAMBN_paper/resultsSGLMCIDem.RData")
