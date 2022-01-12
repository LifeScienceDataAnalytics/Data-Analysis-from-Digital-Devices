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
rm(list=setdiff(ls(),  c("diagnostics", "cogDomains")))
cogDomains = gsub("_VIS1", "", cogDomains)
cogDomains = gsub("SA_", "", cogDomains)
adniClassifierMCItoDementia = read.csv("data_DigMCItoDementia.csv")
#mmseFeatures = grep("MMSE", colnames(adniClassifierCNtoMCI), value = TRUE)
faqFeatures = grep("FAQ", colnames(adniClassifierMCItoDementia), value = TRUE)
digFeatures = grep("AR|^Motor|BIT|DOT", colnames(adniClassifierMCItoDementia), value = TRUE)
adniClassifierMCItoDementia$X = NULL

func.classifier.sgl = function(dataDf,index, type, Class){
  xtrain = subset(dataDf, select=-c(Class))
  ytrain = dataDf$Class
  ntrain=length(ytrain)
  #train.ext=createFolds(ytrain,k=10,returnTrain=TRUE)
  #test.ext=lapply(train.ext,function(x) (1:ntrain)[-x])
  
  rocVec = c()
  for(i in 1:10){
    if(nrow(xtrain[complete.cases(xtrain), ]) < nrow(xtrain)){
      set.seed(i+121)
      xtrain_imp = missForest(xtrain, maxiter = 10, ntree = 100, variablewise = FALSE,
                              decreasing = FALSE, verbose = FALSE,
                              mtry = floor(sqrt(ncol(xtrain))), replace = TRUE,
                              classwt = NULL, cutoff = NULL, strata = NULL,
                              sampsize = NULL, nodesize = NULL, maxnodes = NULL,
                              xtrue = NA, parallelize = c('no', 'variables', 'forests'))
      xtrain_imp = xtrain_imp$ximp
      data = list(x=xtrain_imp, y=ytrain)
    }
    else{
      data = list(x=xtrain, y=ytrain)
    }
    fit = cvSGL(data, index = index, type = "logit", maxit = 1000, thresh = 0.001,
                min.frac = 0.05, nlam =20, gamma = 0.8, nfold = 10, standardize = TRUE,
                verbose = FALSE, step = 1, reset = 10, alpha = 0.95, lambdas = NULL,
                foldid = NULL)
    
    save(fit, file = paste("sgl_R/models_mciDem/",type, "_model","_", i, ".RData", sep = ""))
    minLamdaInd = match(min(fit$fit$lambdas),fit$fit$lambdas)
    roc = (roc.curve(fit$prevals[,minLamdaInd], weights.class0 = ytrain)$auc)
    print(roc)
    rocVec = c(rocVec, roc)
    if(nrow(xtrain[complete.cases(xtrain), ]) < nrow(xtrain)){
      write.csv(xtrain_imp, file = paste("sgl_R/data_mciDem/",type,"_", i, ".csv", sep = ""))
    }
    else{
      write.csv(xtrain, file = paste("sgl_R/data_mciDem/",type,"_", i, ".csv", sep = ""))
    }
    #print(roc.curve(fit$prevals[,minLamdaInd], weights.class1 = ytrain[train.ext[[i]]])$auc)
    # print(minLamdaInd)
    # print(str(xtest_imp))
  }
  return(rocVec)
}


##for Cognitive measures and FAQ measures
adniClassifierMCItoDementiaCogFAQ = adniClassifierMCItoDementia[,c(cogDomains,faqFeatures,"Class")]
##index for Cog domains and FAQ
index = c(rep(1,9), rep(2,10))
Class = "Class"
type = "CogFAQ"
rocCogFAQ = func.classifier.sgl(adniClassifierMCItoDementiaCogFAQ, index = index, type = type, "Class")
boxplot(rocCogFAQ, las = 2, xlab = "", ylab = "AUC", ylim = c(0.5,1),col = "skyblue")

print(mean(sum(rocCogFAQ)/10))
maxRocCogFAQ = match(max(rocCogFAQ), rocCogFAQ)

##load the model having minimum roc
load(paste("sgl_R/models_mciDem/",type, "_model_" ,maxRocCogFAQ, ".RData", sep = ""))
x = read.csv(paste("sgl_R/data_mciDem/",type,"_", maxRocCogFAQ, ".csv", sep = ""))
y= adniClassifierMCItoDementiaCogFAQ$Class

x$X = NULL
finalData = list(x = x, y=y)
finalFit = SGL(finalData, index = index, type = "logit", maxit = 1000, thresh = 0.001,
               min.frac = 0.05, nlam =20, gamma = 0.8, standardize = TRUE,
               verbose = FALSE, step = 1, reset = 10, alpha = 0.95, lambdas = fit$fit$lambdas)

save(finalFit, file = paste("sgl_R/models_mciDem/","final_", type, "_model","_", maxRocCogFAQ, ".RData", sep = ""))

coef = finalFit$beta[,match(min(finalFit$lambdas),finalFit$lambdas)]
absCoef = abs(coef)
coefDf = cbind.data.frame(colnames(x), coef, absCoef, index)
names(coefDf)[names(coefDf)=="colnames(x)"] = "names"
names(coefDf)[names(coefDf)=="index"] = "group"
coefDf <- coefDf[with(coefDf, order(-absCoef)), ]


library(ggplot2) 
coefDfNoZero = subset(coefDf[coefDf$absCoef>0,])
png("output_figures/classifier_ADNI_Cognition_FAQ_MCI_Dem.png", width = 300, height = 225, units='mm',res = 300) 
ggplot(coefDfNoZero, aes(x=reorder(names, absCoef), weight=absCoef, fill=as.factor(group))) + 
  geom_bar() +
  scale_fill_manual(labels = c("Cognitive Domains from DMs","FAQ"), values = c( "dark orange","orangered"), name = "Variable group") +
  ylab("absolute coefficients") +
  xlab("Feature Names")+
  coord_flip()
dev.off()


###overall feature importance of FAQ and Cog Domains##
overall_imp_faq =  coefDfNoZero[coefDfNoZero$names %in% c(grep("FAQ", coefDfNoZero$names, value = TRUE)),]
overall_imp_faq_sum = colSums(overall_imp_faq['absCoef'])/colSums(coefDfNoZero['absCoef'])
round(overall_imp_faq_sum,4)
overall_imp_cog =  coefDfNoZero[coefDfNoZero$names %in% gsub("SA_", "",cogDomains),]
overall_imp_cog_sum = colSums(overall_imp_cog['absCoef'])/colSums(coefDfNoZero['absCoef'])
round(overall_imp_cog_sum,4)
overalImpDf = cbind.data.frame(c("FAQ","Cog from DMs"), c(overall_imp_faq_sum, overall_imp_cog_sum), c("",""))
names(overalImpDf) = c("groups", "importance", "feature")

coefDfNoZeroCp = coefDfNoZero
plot1 = ggplot(coefDfNoZeroCp, aes(x=reorder(names, absCoef), weight=absCoef, fill=as.factor(group))) + 
  geom_bar(show.legend = FALSE) +
  scale_fill_manual(labels = c("Cog from DMs",  "FAQ"), values = c("dark orange","orangered"))+ 
  ylab("absolute coefficients")+
  xlab("Feature Names")+
  coord_flip()
plot2 = ggplot(overalImpDf, aes(fill = groups, y= importance, x = feature)) + geom_bar(position = "stack", stat = "identity", width = 0.1)+scale_fill_manual(values = c("dark orange","orangered", "steelblue"))+ylab("overall feature importance")
png("output_figures/classifier_ADNI_CogFAQ_overallImp_mci_dem.png", width = 300, height = 100, units='mm',res = 300) 
plot_grid(plot1, plot2, align = "h",labels = "AUTO", label_size = 12, rel_heights = c(0.1, 1), rel_widths = c(0.5,0.5))
dev.off()

#for Cognition
##index, for Cognitive measures 
dataCognition =  adniClassifierMCItoDementia[,c(cogDomains, "Class")]
index = c(rep(1,9))
Class = "Class"
type = "Cognition"
rocCog = func.classifier.sgl(dataCognition, index = index, type = type, "Class")
boxplot(rocCog, las = 2, xlab = "", ylab = "AUC", ylim = c(0.8,1),col = "skyblue")

#print(mean(sum(rocCogMMSE)/10))
maxRoc = match(max(rocCog), rocCog)

##final model
rm(fit)
##load the model having minimum roc
load(paste("sgl_R/models/",type, "_model_" ,maxRoc, ".RData", sep = ""))
x = read.csv(paste("sgl_R/data/",type,"_", maxRoc, ".csv", sep = ""))
y= dataCognition$Class
# x = missForest(x, maxiter = 10, ntree = 100, variablewise = FALSE,
#                         decreasing = FALSE, verbose = FALSE,
#                         mtry = floor(sqrt(ncol(x))), replace = TRUE,
#                         classwt = NULL, cutoff = NULL, strata = NULL,
#                         sampsize = NULL, nodesize = NULL, maxnodes = NULL,
#                         xtrue = NA, parallelize = c('no', 'variables', 'forests'))
#x = x$ximp
x$X = NULL
finalData = list(x = x, y=y)
finalFit = SGL(finalData, index = index, type = "logit", maxit = 1000, thresh = 0.001,
               min.frac = 0.05, nlam =20, gamma = 0.8, standardize = TRUE,
               verbose = FALSE, step = 1, reset = 10, alpha = 0.95, lambdas = fit$fit$lambdas
)

save(finalFit, file = paste("sgl_R/models/","final_", type, "_model","_", maxRoc, ".RData", sep = ""))

coef = finalFit$beta[,match(min(finalFit$lambdas),finalFit$lambdas)]
absCoef = abs(coef)
coefDf = cbind.data.frame(colnames(x), coef, absCoef, index)
names(coefDf)[names(coefDf)=="colnames(x)"] = "names"
names(coefDf)[names(coefDf)=="index"] = "group"
coefDf <- coefDf[with(coefDf, order(-absCoef)), ]


library(ggplot2) 
coefDfNoZero = subset(coefDf[coefDf$absCoef>0,])
png("output_figures/classifier_Altoida_Cog.png", width = 300, height = 225, units='mm',res = 300) 
ggplot(coefDfNoZero, aes(x=reorder(names, absCoef), weight=absCoef, fill=as.factor(group))) + 
  geom_bar() +
  scale_fill_manual(labels = "Cogntive domains from DMs", values = "orange", name = "Variable group") +
  ylab("absolute coefficients") +
  xlab("Feature Names")+
  coord_flip()
dev.off()


#for Digital Measure
##index, for digital measures
data = adniClassifierCNtoMCI[,c(digFeatures, "Class")]

featureDesc = featureDesc[featureDesc$technical_group_name%in% groups,]
index = c()
##index for groups of digital features

group_name_vec = c()
for(names in colnames(adniClassifierCNtoMCI[,digFeatures])){
  group_name = featureDesc[featureDesc$feature_name == names, "technical_group_name"]
  group_name = as.character(group_name)
  print(group_name)
  index = c(index,which(groups==group_name))
  group_name_vec = c(group_name_vec, group_name)
}
indexDf = cbind.data.frame(index,group_name_vec)
indexDf = unique(indexDf)

##index for 
Class = "Class"
type = "Dig"
rocDig = func.classifier.sgl(data, index = index, type = type, "Class")
boxplot(rocDig, las = 2, xlab = "", ylab = "AUC", ylim = c(0.8,1),col = "skyblue")

print(mean(sum(rocDig)/10))
maxRoc = match(max(rocDig), rocDig)

##final model
rm(fit)
##load the model having minimum roc
load(paste("sgl_R/models/",type, "_model_" ,maxRoc, ".RData", sep = ""))
#x = subset(data, select=-c(Class))
x = read.csv(paste("sgl_R/data/",type,"_", maxRoc, ".csv", sep = ""))
x$X = NULL
y= data$Class
# x = missForest(x, maxiter = 10, ntree = 100, variablewise = FALSE,
#                decreasing = FALSE, verbose = FALSE,
#                mtry = floor(sqrt(ncol(x))), replace = TRUE,
#                classwt = NULL, cutoff = NULL, strata = NULL,
#                sampsize = NULL, nodesize = NULL, maxnodes = NULL,
#                xtrue = NA, parallelize = c('no', 'variables', 'forests'))
# x = x$ximp
finalData = list(x = x, y=y)
finalFitDig = cvSGL(finalData, index = index, type = "logit", maxit = 1000, thresh = 0.001,
                    min.frac = 0.05, nlam =20, gamma = 0.8, nfold = 10, standardize = TRUE,
                    verbose = FALSE, step = 1, reset = 10, alpha = 0.95, lambdas = fit$fit$lambdas,
                    foldid = NULL)

save(finalFitDig, file = paste("sgl_R/models/","final_", type, "_model","_", i, ".RData", sep = ""))

coef = finalFitDig$fit$beta[,match(min(finalFit$fit$lambdas),finalFit$fit$lambdas)]
absCoef = abs(coef)
coefDf = cbind.data.frame(colnames(x), coef, absCoef, index)
names(coefDf)[names(coefDf)=="colnames(x)"] = "names"
names(coefDf)[names(coefDf)=="index"] = "group"
coefDf <- coefDf[with(coefDf, order(-absCoef)), ]


library(ggplot2) 
coefDfNoZero = subset(coefDf[coefDf$absCoef>0,])
coefDfNoZero = coefDfNoZero[1:10,]
png("output_figures/classifier_Altoida_Dig.png", width = 300, height = 225, units='mm',res = 300) 
require(RColorBrewer)
myColors <- brewer.pal(3, "Oranges")
ggplot(coefDfNoZero, aes(x=reorder(names, absCoef), weight=absCoef, fill=as.factor(group))) + 
  geom_bar() +
  scale_fill_manual(labels = c("ARObjectFinding", "ARObjectPlacement", "MotorDrawingFeatures"), values =myColors, name = "Variable group")+ 
  ylab("absolute coefficients") +
  xlab("Feature Names")+
  coord_flip()
dev.off()

##multiple boxplots together
fill <- "#4271AE"
line <- "#1F3552"
boxplotDf = cbind.data.frame(rocMMSE, rocDig, rocCogMMSE, rocDigMMSE)
names(boxplotDf) = c("MMSE", "DMs", "MMSE and Cog", "MMSE and DMs")
boxplotDf = melt(boxplotDf)
plt <- ggplot(data = boxplotDf, aes(x = variable, y = value))
png("output_figures/AUC_boxplot.png", width = 200, height = 155, units='mm',res = 300) 
plt + geom_boxplot(fill = fill, colour = line) + labs(x = "Features", y = "AUC") + ylim(0.8,1)
dev.off()
