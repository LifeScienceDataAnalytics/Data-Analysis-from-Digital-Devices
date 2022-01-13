setwd("~/Altoida_VAMBN_paper")
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
library(caret)
library(pdp)
library(groupdata2)


load("~/createbnAltoida.RData")
rm(list=setdiff(ls(),  c("targetNodes", "dataCP", "demog", "diagnostics", "cogDomains")))
cogDomains = gsub("_VIS1", "", cogDomains)
cogDomains = gsub("SA_", "", cogDomains)
data = read.csv("altoidaDf.csv")
data = data[data$DX %in% c(0,1,2),]
data$DX[data$DX == 2] = 1
names(data)[names(data) == "DX"] = "Class"
data$X = NULL
mmseFeatures = grep("MMSE", colnames(data), value = TRUE)

data$SUBJID = gsub("_[0-9]$", "", data$SUBJID)

func.classifier.sgl.best.param = function(dataDf,index, type, Class){
  xtrain = subset(dataDf, select=-c(Class))
  ytrain = dataDf$Class
  ntrain=length(ytrain)
  set.seed(123)
  ##changed for revision Dec 21###
  train.ext = groupKFold(xtrain$SUBJID, k = 5) 
  #train.ext=createFolds(ytrain,k=5,returnTrain=TRUE)
  ##changed for revision Dec 21###
  test.ext=lapply(train.ext,function(x) (1:ntrain)[-x])
  seedVec = c()
  alphas = seq(0.1,1,0.1)
  ##changed for revision Dec 21
  alphas = setdiff(alphas,"0.6")
  alphaVec = c()
  meanRocVec = c()
  lambdasMin = c()
  ##changed for revision Dec 21###
  rownames(xtrain) = paste(xtrain$SUBJID,"_", row.names(xtrain),sep = "")
  xtrain$SUBJID = NULL
  ##changed for revision Dec 21###
  for(alpha in alphas){
    print(alpha)
    rocVec = c()
    lambdas = c()
    for(i in 1:5){
      print(i)
      set.seed(i+121)
      val = i+121
      xtrain_imp = missForest(xtrain[train.ext[[i]],], maxiter = 10, ntree = 100, variablewise = FALSE,
                                                 decreasing = FALSE, verbose = FALSE,
                                                 mtry = floor(sqrt(ncol(xtrain))), replace = TRUE,
                                                 classwt = NULL, cutoff = NULL, strata = NULL,
                                                 sampsize = NULL, nodesize = NULL, maxnodes = NULL,
                                                 xtrue = NA)

      set.seed(i+121)
      xtest_imp = missForest(xtrain[test.ext[[i]],], maxiter = 10, ntree = 100, variablewise = FALSE,
                            decreasing = FALSE, verbose = FALSE,
                            mtry = floor(sqrt(ncol(xtrain))), replace = TRUE,
                            classwt = NULL, cutoff = NULL, strata = NULL,
                            sampsize = NULL, nodesize = NULL, maxnodes = NULL,
                            xtrue = NA)
      xtrain_imp = xtrain_imp$ximp
      ##changed for revision Dec 21###
      xtrain_imp["SUBJID"] = rownames(xtrain_imp)
      vec = rownames(xtrain_imp)
      xtrain_imp["ID"] = gsub("_[0-9]*$","",rownames(xtrain_imp))
      xtrain_imp$SUBJID = as.factor(xtrain_imp$SUBJID)
      xtrain_imp$ID = as.factor(xtrain_imp$ID)
      # ##changed for revision Dec 21###
      xtrain_imp <- fold(
        data = xtrain_imp,
        k = 10,
        id_col = "ID",
        method = "n_dist"
      )
      xtrain_imp = as.data.frame(xtrain_imp)
      ##to order with the same vector order##
      xtrain_imp = xtrain_imp[match(vec, xtrain_imp$SUBJID),]
      rownames(xtrain_imp) = xtrain_imp$SUBJID
      xtest_imp = xtest_imp$ximp
      if(type == "CogMMSE"|type == "DigMMSE"){
        xtrain_imp[,grep("MMSE", colnames(xtrain_imp), value = TRUE)]=lapply(xtrain_imp[,grep("MMSE", colnames(xtrain_imp), value = TRUE)],
                                                                            as.numeric)
        xtest_imp[,grep("MMSE", colnames(xtest_imp), value = TRUE)]=lapply(xtest_imp[,grep("MMSE", colnames(xtest_imp), value = TRUE)],
                                                                        as.numeric)
      }
      if((type == "CogMMSE_only")| (type == "DigMMSE_only")){
        xtrain_imp$MMSE = as.numeric(xtrain_imp$MMSE)
        xtest_imp$MMSE = as.numeric(xtest_imp$MMSE)
      }
      ##changed for revision Dec 21###
      foldid = as.integer(xtrain_imp$.folds)
      xtrain_imp$.folds = NULL
      data = list(x=xtrain_imp, y=ytrain[train.ext[[i]]])
      ##changed for revision Dec 21###
      data$x$SUBJID = NULL
      data$x$ID = NULL
      data$x$.folds = NULL
      fit = cvSGL(data, index = index, type = "logit", maxit = 1000, thresh = 0.001,
              min.frac = 0.05, nlam =10, gamma = 0.8, nfold = 10, standardize = TRUE,
              verbose = FALSE, step = 1, reset = 10, alpha = alpha, lambdas = NULL,
              foldid = foldid)
      minLambda = min(fit$fit$lambdas)
      trainFit = SGL(data, index = index, type = "logit", maxit = 1000, thresh = 0.001,
                     min.frac = 0.1,nlam=2, gamma = 0.8, standardize = TRUE,
                     verbose = FALSE, step = 1, reset = 10, alpha = alpha, lambdas=c(min(fit$fit$lambdas),1))
       yh = predictSGL(trainFit, as.matrix(xtest_imp), 1)
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
  bestAlpha = alphaVec[maxRoc]
  bestLambda = lambdasMin[maxRoc]
  print(bestAlpha)
  return(list(bestAlpha, bestLambda, meanRocVec))
}


##index, for Cognitive measures and MMSE measures
altoidaClassifierCogMMSE = data[,c(cogDomains, mmseFeatures, "Class", "SUBJID")]
altoidaClassifierCogMMSE$MMSE = NULL
index = c(rep(1,9), rep(2,5))
Class = "Class"
type = "CogMMSE"

altoidaClassifierCogMMSE[,grep("MMSE", colnames(altoidaClassifierCogMMSE), value = TRUE)]=lapply(altoidaClassifierCogMMSE[,grep("MMSE", colnames(altoidaClassifierCogMMSE), value = TRUE)],
                                                                 as.factor)

paramCogMMSE = func.classifier.sgl.best.param(altoidaClassifierCogMMSE, index = index, type = type, "Class")
alpha = paramCogMMSE[[1]]
lambda  = min(paramCogMMSE[[2]])
altoidaClassifierCogMMSE$SUBJID = NULL
x = subset(altoidaClassifierCogMMSE, select=-c(Class))
y= altoidaClassifierCogMMSE$Class

set.seed(123)
x = missForest(x, maxiter = 10, ntree = 100, variablewise = FALSE,
                        decreasing = FALSE, verbose = FALSE,
                        mtry = floor(sqrt(ncol(x))), replace = TRUE,
                        classwt = NULL, cutoff = NULL, strata = NULL,
                        sampsize = NULL, nodesize = NULL, maxnodes = NULL,
                        xtrue = NA)
x = x$ximp
x[,grep("MMSE", colnames(x), value = TRUE)]=lapply(x[,grep("MMSE", colnames(x), value = TRUE)],
                                                                   as.numeric)

finalData = list(x = x, y=y)
finalFit = SGL(finalData, index = index, type = "logit", maxit = 1000, thresh = 0.001,
            min.frac = 0.1, nlam =2, gamma = 0.8, standardize = TRUE,
            verbose = FALSE, step = 1, reset = 10, alpha = alpha, lambdas = c(lambda,1))

save(finalFit, file = paste("sgl_R/modelsNew/","final_", type, "_model", ".RData", sep = ""))


coef = finalFit$beta[,match(min(finalFit$lambdas),finalFit$lambdas)]
absCoef = abs(coef)
coefDf = cbind.data.frame(colnames(x), coef, absCoef, index)
names(coefDf)[names(coefDf)=="colnames(x)"] = "names"
names(coefDf)[names(coefDf)=="index"] = "group"
coefDf <- coefDf[with(coefDf, order(-absCoef)), ]


library(ggplot2) 
coefDfNoZero = subset(coefDf[coefDf$absCoef>0,])
png("sgl_R/plots_revised/classifier_Altoida_Cognition.png", width = 300, height = 225, units='mm',res = 300) 
ggplot(coefDfNoZero, aes(x=reorder(names, absCoef), weight=absCoef, fill=as.factor(group))) + 
  geom_bar() +
  scale_fill_manual(labels = c("MMSE", "Cognitive Domains from DMs"), values = c("dark orange", "steelblue"), name = "Variable group") +
  ylab("absolute coefficients") +
  xlab("Feature Names")+
  coord_flip()
dev.off()

##for total MMSE and Cognition
##index, for Cognitive measures and MMSE measures
altoidaClassifierCogMMSE_only = data[,c(cogDomains, "MMSE", "Class", "SUBJID")]
index = c(rep(1,9), rep(2,1))
Class = "Class"
type = "CogMMSE_only"
altoidaClassifierCogMMSE_only$MMSE = as.factor(altoidaClassifierCogMMSE_only$MMSE)

paramCogMMSE_only = func.classifier.sgl.best.param(altoidaClassifierCogMMSE_only, index = index, type = type, "Class")
alphaCogMMSE_only = paramCogMMSE_only[[1]]
lambdaCogMMSE_only = paramCogMMSE_only[[2]]

altoidaClassifierCogMMSE_only$SUBJID = NULL
x = subset(altoidaClassifierCogMMSE_only, select=-c(Class))
y= altoidaClassifierCogMMSE_only$Class
set.seed(123)
x = missForest(x, maxiter = 10, ntree = 500, variablewise = FALSE,
               decreasing = FALSE, verbose = FALSE,
               mtry = floor(sqrt(ncol(x))), replace = TRUE,
               classwt = NULL, cutoff = NULL, strata = NULL,
               sampsize = NULL, nodesize = NULL, maxnodes = NULL,
               xtrue = NA)
x = x$ximp
x$MMSE = as.numeric(x$MMSE)
finalData = list(x = x, y=y)
finalFit_only = SGL(finalData, index = index, type = "logit", maxit = 1000, thresh = 0.001,
               min.frac = 0.1, nlam =2, gamma = 0.8, standardize = TRUE,
               verbose = FALSE, step = 1, reset = 10, alpha = alphaCogMMSE_only, lambdas = c(lambdaCogMMSE_only,1))

save(finalFit_only, file = paste("sgl_R/models_revised/","final_", type, "_model", ".RData", sep = ""))

coef_only = finalFit_only$beta[,match(min(finalFit_only$lambdas),finalFit_only$lambdas)]
absCoef_only = abs(coef_only)
coefDf_only = cbind.data.frame(colnames(x), coef_only, absCoef_only, index)
names(coefDf_only)[names(coefDf_only)=="colnames(x)"] = "names"
names(coefDf_only)[names(coefDf_only)=="index"] = "group"
coefDf_only <- coefDf_only[with(coefDf_only, order(-absCoef_only)), ]


library(ggplot2) 
coefDfNoZero_only = subset(coefDf_only[coefDf_only$absCoef_only>0,])
png("sgl_R/plots_revised/classifier_Altoida_Cognition_only.png", width = 300, height = 225, units='mm',res = 300) 
ggplot(coefDfNoZero_only, aes(x=reorder(names, absCoef_only), weight=absCoef_only, fill=as.factor(group))) + 
  geom_bar() +
  scale_fill_manual(labels = c("MMSE", "Cognitive Domains from DMs"), values = c("dark orange", "steelblue"), name = "Variable group") +
  ylab("absolute coefficients") +
  xlab("Feature Names")+
  coord_flip()
dev.off()


###overall feature importance of Cog and MMSE##
overall_imp_mmse_only =  coefDfNoZero_only[coefDfNoZero_only$names %in% c(grep("MMSE", coefDfNoZero_only$names, value = TRUE)),]
overall_imp_mmse_only_sum = colSums(overall_imp_mmse_only['absCoef_only'])/colSums(coefDfNoZero_only['absCoef_only'])
round(overall_imp_mmse_only_sum,4)
overall_imp_cog_only =  coefDfNoZero_only[coefDfNoZero_only$names %in% gsub("SA_", "",cogDomains),]
overall_imp_cog_only_sum = colSums(overall_imp_cog_only['absCoef_only'])/colSums(coefDfNoZero_only['absCoef_only'])
round(overall_imp_cog_only_sum,4)
overalImpDf_only = cbind.data.frame(c("MMSE", "Digital Cognitive Domains"), c(overall_imp_mmse_only_sum, overall_imp_cog_only_sum), c("",""))
names(overalImpDf_only) = c("groups", "importance", "feature")


coefDfNoZeroCp_only = coefDfNoZero_only
plot1 = ggplot(coefDfNoZero, aes(x=reorder(names, absCoef), weight=absCoef, fill=as.factor(group))) + 
  geom_bar(show.legend = FALSE) +
  scale_fill_manual(labels = c("Digital Cognitive Domains", "MMSE"), values = c("orange", "steelblue"))+ 
  ylab("absolute coefficients") +
  xlab("feature Names")+
  coord_flip()
plot2 = ggplot(overalImpDf_only, aes(fill = groups, y= importance, x = feature)) + geom_bar(position = "stack", stat = "identity", width = 0.1)+scale_fill_manual(values = c("orange", "steelblue"))+ylab("overall feature importance")
png("sgl_R/plots_revised/classifier_Altoida_CogMMSE_overallImp_only.png", width = 300, height = 100, units='mm',res = 300) 
plot_grid(plot1, plot2, align = "h",labels = "AUTO", label_size = 12, rel_widths = c(1.5,1), ncol = 2)
grid.arrange(plot1, plot2,ncol=2)
dev.off()
save.image("~/Altoida_VAMBN_paper/sgl_R/cogMMSE.RData")


##for digital and MMSE##
##subset the digital features for the target nodes
targetNodes = gsub("zcode_", "", targetNodes)
targetNodes = gsub("_VIS1", "", targetNodes)
groups = targetNodes

featureDesc =read_delim("~/.....csv", ";", escape_double = FALSE, col_names = TRUE, trim_ws = TRUE)
targetNodesFeatures = c()
for(group in groups){
  features = as.vector(featureDesc[featureDesc$technical_group_name == group, "feature_name"])
  feat = features$feature_name 
  targetNodesFeatures = c(targetNodesFeatures,feat)
}
targetNodesFeatures = intersect(targetNodesFeatures, colnames(data))
altoidaClassifierDigMMSE = data[,c(targetNodesFeatures, mmseFeatures, "Class", "SUBJID")]
altoidaClassifierDigMMSE$MMSE = NULL
featureDesc = featureDesc[featureDesc$technical_group_name%in% groups,]
index = c()
##index for groups of digital features

group_name_vec = c()
for(names in colnames(altoidaClassifierDigMMSE[,targetNodesFeatures])){
  group_name = featureDesc[featureDesc$feature_name == names, "technical_group_name"]
  group_name = as.character(group_name)
  print(group_name)
  index = c(index,which(groups==group_name))
  group_name_vec = c(group_name_vec, group_name)
}

indexDf = cbind.data.frame(index,group_name_vec)
indexDf = unique(indexDf)
##index for 
indexDf = rbind.data.frame(indexDf, c((max(index)+1), "MMSE"))
index = c(index, rep((max(index)+1),5))
Class = "Class"
type = "DigMMSE"

altoidaClassifierDigMMSE[,grep("MMSE", colnames(altoidaClassifierDigMMSE), value = TRUE)]=lapply(altoidaClassifierDigMMSE[,grep("MMSE", colnames(altoidaClassifierDigMMSE), value = TRUE)],
                                                                                                 as.factor)

paramDigMMSE = func.classifier.sgl.best.param(altoidaClassifierDigMMSE, index = index, type = type, "Class")
alphaDigMMSE = paramDigMMSE[[1]]
lambaDigMMSE = paramDigMMSE[[2]]

##train the final Model##
altoidaClassifierDigMMSE$SUBJID = NULL
x = subset(altoidaClassifierDigMMSE, select=-c(Class))
y= altoidaClassifierDigMMSE$Class
set.seed(123)
x = missForest(x, maxiter = 10, ntree = 100, variablewise = FALSE,
               decreasing = FALSE, verbose = FALSE,
               mtry = floor(sqrt(ncol(x))), replace = TRUE,
               classwt = NULL, cutoff = NULL, strata = NULL,
               sampsize = NULL, nodesize = NULL, maxnodes = NULL,
               xtrue = NA)
x = x$ximp
x[,grep("MMSE", colnames(x), value = TRUE)]=lapply(x[,grep("MMSE", colnames(x), value = TRUE)],
                                                   as.numeric)
finalDataDigMMSE = list(x = x, y=y)
finalFitDigMMSE = SGL(finalDataDigMMSE, index = index, type = "logit", maxit = 1000, thresh = 0.001,
                    min.frac = 0.1, nlam =2, gamma = 0.8, standardize = TRUE,
                    verbose = FALSE, step = 1, reset = 10, alpha = alphaDigMMSE, lambdas = c(lambaDigMMSE,1))

save(finalFitDigMMSE, file = paste("sgl_R/models_revised/","final_", type, "_model", ".RData", sep = ""))


coef = finalFitDigMMSE$beta[,match(min(finalFitDigMMSE$lambdas),finalFitDigMMSE$lambdas)]
absCoef = abs(coef)
coefDf = cbind.data.frame(colnames(x), coef, absCoef, index)
names(coefDf)[names(coefDf)=="colnames(x)"] = "names"
names(coefDf)[names(coefDf)=="index"] = "group"
coefDf <- coefDf[with(coefDf, order(-absCoef)), ]
coefDf["label"] <- indexDf$group_name_vec[match(coefDf$group, indexDf$index)]
write.csv(coefDf, "featureImpCoefDfMMSEDigAlt.csv")
library(ggplot2) 
coefDfNoZero = subset(coefDf[coefDf$absCoef>0,])
write.csv(coefDfNoZero, "featureImpCoefDfMMSEDigAlt.csv")
coefDfNoZeroTop15 = coefDfNoZero[1:15,]

png("sgl_R/plots_revised/classifier_Altoida_DigMMSE.png", width = 300, height = 100, units='mm',res = 300) 
require(RColorBrewer)
if("MMSE" %in% coefDfNoZeroTop15$label){
  myColors <- brewer.pal(length(unique(coefDfNoZeroTop15$label))-1, "Oranges")
  ggplot(coefDfNoZeroTop15, aes(x=reorder(names, absCoef), weight=absCoef, fill = as.factor(group))) + 
    geom_bar() +
    scale_fill_manual(labels = c("ARObjectFinding", "MotorTestDurations", "MotorTappingFeatures",
                                 "ARObjectPlacement", "MMSE"), values = c(myColors, "steelblue"), name = "Variable group")+ 
    ylab("absolute coefficients") +
    xlab("Feature Names")+
    coord_flip()
}else{
  myColors <- brewer.pal(length(unique(coefDfNoZeroTop15$label)), "Oranges")
  ggplot(coefDfNoZeroTop15, aes(x=reorder(names, absCoef), weight=absCoef, fill=as.factor(group))) + 
    geom_bar() +
    scale_fill_manual(labels = unique(coefDfNoZeroTop15$label), values = c(myColors), name = "Variable group")+ 
    ylab("absolute coefficients") +
    xlab("Feature Names")+
    coord_flip()
}
dev.off()

##############for digital and only MMSE
altoidaClassifierDigMMSE_only = data[,c(targetNodesFeatures, "MMSE", "Class", "SUBJID")]
featureDesc = featureDesc[featureDesc$technical_group_name%in% groups,]
index = c()
##index for groups of digital features

group_name_vec = c()
for(names in colnames(altoidaClassifierDigMMSE[,targetNodesFeatures])){
  group_name = featureDesc[featureDesc$feature_name == names, "technical_group_name"]
  group_name = as.character(group_name)
  print(group_name)
  index = c(index,which(groups==group_name))
  group_name_vec = c(group_name_vec, group_name)
}

indexDf = cbind.data.frame(index,group_name_vec)
indexDf = unique(indexDf)
##index for 
indexDf = rbind.data.frame(indexDf, c((max(index)+1), "MMSE"))
index = c(index, rep((max(index)+1),1))
Class = "Class"
type = "DigMMSE_only"
paramDigMMSE_only = func.classifier.sgl.best.param(altoidaClassifierDigMMSE_only, index = index, type = type, "Class")
alphaDigMMSE_only = paramDigMMSE_only[[1]]
lambdaDigMMSE_only = paramDigMMSE_only[[2]]
altoidaClassifierDigMMSE_only$SUBJID = NULL
x = subset(altoidaClassifierDigMMSE_only, select=-c(Class))
y= altoidaClassifierDigMMSE_only$Class
set.seed(123)
x = missForest(x, maxiter = 10, ntree = 100, variablewise = FALSE,
               decreasing = FALSE, verbose = FALSE,
               mtry = floor(sqrt(ncol(x))), replace = TRUE,
               classwt = NULL, cutoff = NULL, strata = NULL,
               sampsize = NULL, nodesize = NULL, maxnodes = NULL,
               xtrue = NA)
x = x$ximp

finalData_only = list(x = x, y=y)
finalFitDigMMSE_only = SGL(finalData_only, index = index, type = "logit", maxit = 1000, thresh = 0.001,
                      min.frac = 0.1, nlam =2, gamma = 0.8, standardize = TRUE,
                      verbose = FALSE, step = 1, reset = 10, alpha = alphaDigMMSE_only, lambdas = c(lambdaDigMMSE_only,1))


save(finalFitDigMMSE_only, file = paste("sgl_R/models_revised/","final_", type, "_model", ".RData", sep = ""))

coef_only  = finalFitDigMMSE_only$beta[,match(min(finalFitDigMMSE_only$lambdas),finalFitDigMMSE_only$lambdas)]
absCoef_only = abs(coef_only)
coefDf_only = cbind.data.frame(colnames(x), coef_only, absCoef_only, index)
names(coefDf_only)[names(coefDf_only)=="colnames(x)"] = "names"
names(coefDf_only)[names(coefDf_only)=="index"] = "group"
coefDf_only <- coefDf_only[with(coefDf_only, order(-absCoef_only)), ]
coefDf_only["label"] <- indexDf$group_name_vec[match(coefDf_only$group, indexDf$index)]

library(ggplot2) 
coefDfNoZero_only = subset(coefDf_only[coefDf_only$absCoef_only>0,])
coefDfNoZeroTop15_only = coefDfNoZero_only[1:15,]

###overall feature importance of Dig and MMSE##
overall_imp_mmse_only =  coefDfNoZero_only[coefDfNoZero_only$names %in% c(grep("MMSE", coefDfNoZero_only$names, value = TRUE)),]
overall_imp_mmse_only_sum = colSums(overall_imp_mmse_only['absCoef_only'])/colSums(coefDfNoZero_only['absCoef_only'])
round(overall_imp_mmse_only_sum,4)
overall_imp_dig_only =  coefDfNoZero_only[coefDfNoZero_only$names %in% c(grep("AR|BIT|DOT|Motor", coefDfNoZero_only$names, value = TRUE)),]
overall_imp_dig_only_sum = colSums(overall_imp_dig_only['absCoef_only'])/colSums(coefDfNoZero_only['absCoef_only'])
round(overall_imp_dig_only_sum,4)
overalImpDf_only = cbind.data.frame(c("MMSE", "Digital tasks"), c(overall_imp_mmse_only_sum, overall_imp_dig_only_sum), c("",""))
names(overalImpDf_only) = c("groups", "importance", "feature")

coefDfNoZeroCp_only = coefDfNoZeroTop15_only
coefDfNoZeroCp_only$group[coefDfNoZeroCp_only$group < unique(coefDf_only[coefDf_only$label=="MMSE", "group"])] =1
plot1 =   ggplot(coefDfNoZeroTop15, aes(x=reorder(names, absCoef), weight=absCoef, fill = as.factor(group))) + 
  geom_bar() +scale_fill_manual(labels = c("ARObjectFinding", "MotorTestDurations", "MotorTappingFeatures",
                               "ARObjectPlacement", "MMSE"), values = c(myColors[1:4], "steelblue"), name = "Variable group")+ 
  ylab("absolute coefficients") +
  xlab("Feature Names")+
  coord_flip()
plot2 = ggplot(overalImpDf_only, aes(fill = groups, y= importance, x = feature)) + geom_bar(position = "stack", stat = "identity", width = 0.1)+scale_fill_manual(values = c("#FB6A4A", "steelblue"))+ylab("overall feature importance")
png("sgl_R/plots_revised/classifier_Altoida_DigMMSE_overallImp_only.png", width = 300, height = 100, units='mm',res = 300) 
plot_grid(plot1, plot2, align = "h", label_size = 12, rel_widths = c(2,1), ncol = 2)
dev.off()

save.image("~/Altoida_VAMBN_paper/classifiersgl_upd.RData")

#for Digital Measure
##index, for digital measures
dataDig = data[,c(targetNodesFeatures, "Class", "SUBJID")]

featureDesc = featureDesc[featureDesc$technical_group_name%in% groups,]
index = c()
##index for groups of digital features

group_name_vec = c()
for(names in colnames(altoidaClassifierDigMMSE[,targetNodesFeatures])){
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
alphaDig = paramDig[[1]]
lambdaDig = paramDig[[2]]

dataDig$SUBJID = NULL
x = subset(dataDig, select=-c(Class))
y= dataDig$Class
set.seed(123)
x = missForest(x, maxiter = 10, ntree = 100, variablewise = FALSE,
               decreasing = FALSE, verbose = FALSE,
               mtry = floor(sqrt(ncol(x))), replace = TRUE,
               classwt = NULL, cutoff = NULL, strata = NULL,
               sampsize = NULL, nodesize = NULL, maxnodes = NULL,
               xtrue = NA)

x = x$ximp
finalData = list(x = x, y=y)

finalFitDig = SGL(finalData, index = indexDig, type = "logit", maxit = 1000, thresh = 0.001,
                  min.frac = 0.1, nlam =2, gamma = 0.8, standardize = TRUE,
                  verbose = FALSE, step = 1, reset = 10, alpha =alphaDig, lambdas = c(min(lambdaDig),1))

save(finalFitDig, file = paste("sgl_R/models_revised/","final_", type, "_model", ".RData", sep = ""))

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
png("sgl_R/plots_revised/classifier_Altoida_Dig.png", width = 300, height = 200, units='mm',res = 300) 
require(RColorBrewer)
myColors <- brewer.pal(length(unique(coefDfNoZero$label)), "Oranges")
ggplot(coefDfNoZero, aes(x=reorder(names, absCoef), weight=absCoef, fill=as.factor(group))) + 
  geom_bar() +
  scale_fill_manual(labels = c("ARObjectFinding", "MotorTestDurations", "ARObjectPlacement"), values =myColors, name = "Variable group")+ 
  ylab("absolute coefficients") +
  xlab("Feature Names")+theme(axis.text.x = element_text(size=14),
                                 axis.text.y = element_text(size=14),
                              axis.title=element_text(size=16))+
  theme(legend.text=element_text(size=12))+
  coord_flip()
dev.off()
save.image("~/Altoida_VAMBN_paper/classifiersgl_upd.RData")


auc <- function(data, index,bestAlpha, bestLambda, seed){
  xtrain = subset(data, select=-c(Class))
  xtrainCp = xtrain
  ytrain = data$Class
  ntrain=length(ytrain)
  set.seed(123)
  ## change December 2021, not to keep same id in training and testing
  #selected_index = createMultiFolds(ytrain, k=5, times = 10) 
  rocVec = list()
  f1Vec = list()
  pltVec = c()
  predictionsDfVec = data.frame()
  for(rep in 1:10){
    #rep = 1
    print(paste("rep",rep))
    set.seed(rep + 123)
    # Data split 
    ## change December 2021, not to keep same id in training and testing
    selected_index = groupKFold(xtrainCp$SUBJID, k = 5) 
    print(paste("selected_index",selected_index))
    test.ext=lapply(selected_index,function(x) (1:ntrain)[-x])
    for(i in 1:length(selected_index)){
    #i = 1
      #if CogMMSE##
    if(rep==1&&i==2|rep==5&&i==2){
      next
    }
    print(paste("i",i))
    f_index = unlist(selected_index[[i]])
    xtrain$SUBJID = NULL
    rf_train = xtrain[f_index,]
    rf_test = xtrain[-f_index,]
    y_train=ytrain[f_index]
    if(sum(is.na(rf_train))>0){
    set.seed(seed+i)
    xtrain_imp = missForest(rf_train, maxiter = 2, ntree = 100, variablewise = FALSE,
                            decreasing = FALSE, verbose = FALSE,
                            mtry = floor(sqrt(ncol(xtrain))), replace = TRUE,
                            classwt = NULL, cutoff = NULL, strata = NULL,
                            sampsize = NULL, nodesize = NULL, maxnodes = NULL,
                            xtrue = NA)
    xtrain_imp = xtrain_imp$ximp
    }
    else{
      xtrain_imp = rf_train
    }
    if(sum(is.na(rf_test))>0){
    set.seed(seed+i)
    xtest_imp = missForest(rf_test, maxiter = 2, ntree = 100, variablewise = FALSE,
                           decreasing = FALSE, verbose = FALSE,
                           mtry = floor(sqrt(ncol(xtrain))), replace = TRUE,
                           classwt = NULL, cutoff = NULL, strata = NULL,
                           sampsize = NULL, nodesize = NULL, maxnodes = NULL,
                           xtrue = NA)

    xtest_imp = xtest_imp$ximp
    }
    else{
      xtest_imp = rf_test
    }
    y_train=ytrain[f_index]
    xtrain_imp[,grep("MMSE", colnames(xtrain_imp), value = TRUE)]=lapply(xtrain_imp[,grep("MMSE", colnames(xtrain_imp), value = TRUE)],
                                                                         as.numeric)
    xtest_imp[,grep("MMSE", colnames(xtest_imp), value = TRUE)]=lapply(xtest_imp[,grep("MMSE", colnames(xtest_imp), value = TRUE)],
                                                                       as.numeric)
    
    data = list(x=xtrain_imp, y = y_train)
    trainFit = SGL(data, index = index, type = "logit", maxit =1000, thresh = 0.001,
                   min.frac = 0.1,nlam=2, gamma = 0.8, standardize = TRUE,
                   verbose = FALSE, step = 1, reset = 10, alpha = bestAlpha, lambdas=c(bestLambda,1))
    yh = predictSGL(trainFit, as.matrix(xtest_imp), 1)
    actualLabel =  ytrain[-f_index]
    predictedLabel =  yh
    folds = rep(i,length(actualLabel) )
    predictionsDf = cbind.data.frame(folds, predictedLabel, actualLabel)
    predictionsDfVec = rbind.data.frame(predictionsDfVec, predictionsDf)
    pltRoc = roc.curve(yh, weights.class0 = ytrain[-f_index], curve = TRUE)
    roc = (roc.curve(yh, weights.class0 = ytrain[-f_index], curve = TRUE)$auc)
    rocVec = cbind(rocVec, roc)
    pltVec[[length(pltVec)+1]] = pltRoc
    #Rename columns
    last_index = ncol(rocVec)
    fold_name = paste("Fold" , as.character(i))
    colnames(rocVec)[last_index] <- fold_name
    ##calculating f1 score, change Jan 20222 for the review
    predictedLabelUpd = ifelse(predictedLabel <0.5,0,1)
    predictedLabelUpd = as.factor(as.vector(predictedLabelUpd))
    actualLabel = as.factor(actualLabel)
    #confusionMatrix(actualLabel,predictedLabelUpd)
    precision <- posPredValue(predictedLabelUpd, actualLabel, positive="0")
    recall <- sensitivity(predictedLabelUpd, actualLabel, positive="0")
    f1 = 2*((precision*recall)/(precision+recall))
    f1Vec = cbind(f1Vec, f1)
    }
  }
  return(list(rocVec,pltVec, predictionsDfVec, f1Vec))
}

                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                ##cogMMSE
##index, for Cognitive measures and MMSE measures
altoidaClassifierCogMMSE = data[,c(cogDomains, mmseFeatures, "Class", "SUBJID")]
altoidaClassifierCogMMSE$MMSE = NULL


altoidaClassifierCogMMSE[,grep("MMSE", colnames(altoidaClassifierCogMMSE), value = TRUE)]=lapply(altoidaClassifierCogMMSE[,grep("MMSE", colnames(altoidaClassifierCogMMSE), value = TRUE)],
                                                                                                 as.factor)


indexCogMMSE = c(rep(1,9), rep(2,5))
Class = "Class"
type = "CogMMSE"
alphaCogMMSE = paramCogMMSE[[1]]
lambdaCogMMSE  = min(paramCogMMSE[[2]])

rocCogMMSE= auc(altoidaClassifierCogMMSE,indexCogMMSE, alphaCogMMSE, lambdaCogMMSE, seed = 123)
rocCogMMSEAUC = rocCogMMSE[[1]]
rocCogMMSEdf_mean =  mean(unlist(rocCogMMSEAUC)) 
boxplot(unlist(rocCogMMSEAUC), col = "#00BFC4")

f1CogMMSE = rocCogMMSE[[4]]
f1CogMMSEdf_mean =  mean(unlist(f1CogMMSE)) 
boxplot(unlist(f1CogMMSE), col = "#00BFC4")


##digMMSE
targetNodes = gsub("zcode_", "", targetNodes)
targetNodes = gsub("_VIS1", "", targetNodes)
groups = targetNodes

featureDesc =read_delim("~/Altoida_VAMBN_paper/ALTOIDA/fraunhofer_v5_technical_grouping.csv", ";", escape_double = FALSE, col_names = TRUE, trim_ws = TRUE)
targetNodesFeatures = c()
for(group in groups){
  features = as.vector(featureDesc[featureDesc$technical_group_name == group, "feature_name"])
  feat = features$feature_name 
  targetNodesFeatures = c(targetNodesFeatures,feat)
}
targetNodesFeatures = intersect(targetNodesFeatures, colnames(data))
altoidaClassifierDigMMSE = data[,c(targetNodesFeatures, mmseFeatures, "Class","SUBJID")]
altoidaClassifierDigMMSE$MMSE = NULL
featureDesc = featureDesc[featureDesc$technical_group_name%in% groups,]
index = c()
##index for groups of digital features

group_name_vec = c()
for(names in colnames(altoidaClassifierDigMMSE[,targetNodesFeatures])){
  group_name = featureDesc[featureDesc$feature_name == names, "technical_group_name"]
  group_name = as.character(group_name)
  print(group_name)
  index = c(index,which(groups==group_name))
  group_name_vec = c(group_name_vec, group_name)
}

indexDf = cbind.data.frame(index,group_name_vec)
indexDf = unique(indexDf)
##index for 
indexDf = rbind.data.frame(indexDf, c((max(index)+1), "MMSE"))
indexDigMMSE = c(index, rep((max(index)+1),5))
altoidaClassifierDigMMSE[,grep("MMSE", colnames(altoidaClassifierDigMMSE), value = TRUE)]=lapply(altoidaClassifierDigMMSE[,grep("MMSE", colnames(altoidaClassifierDigMMSE), value = TRUE)],
                                                                                                 as.factor)

rocDigMMSE= auc(altoidaClassifierDigMMSE,indexDigMMSE, alphaDigMMSE, lambaDigMMSE, seed = 123)
rocDigMMSEAUC = rocDigMMSE[[1]]
rocDigMMSEdf_mean =  mean(unlist(rocDigMMSEAUC) ) 
boxplot(unlist(rocDigMMSEdf_mean), col = "#00BFC4")

f1DigMMSE = rocDigMMSE[[4]]
f1DigMMSEdf_mean =  mean(unlist(f1DigMMSE)) 
boxplot(unlist(f1DigMMSE), col = "#00BFC4")



fill <- "#4271AE"
line <- "#1F3552"

##Dig

rocDig = auc(dataDig, indexDig, bestAlpha = alphaDig, bestLambda = lambdaDig, seed = 123)
rocDigAUC = rocDig[[1]]
rocDigdf_mean =  mean(unlist(rocDigAUC) ) 
boxplot(unlist(rocDigdf_mean), col = "#00BFC4")

f1Dig = rocDig[[4]]
f1Digdf_mean =  mean(unlist(f1Dig)) 
boxplot(unlist(f1Dig), col = "#00BFC4")


predictionsDfVec  = rocDig[[3]]
write.csv(predictionsDfVec, "predictionsDfVec.csv")
prediction = list()
labels = list()

for(i in unique(predictionsDfVec$folds)){
  prediction[[length(prediction)+1]] = predictionsDfVec[predictionsDfVec$folds == i, "predictedLabel"]
}

for(i in unique(predictionsDfVec$folds)){
  labels[[length(labels)+1]] = predictionsDfVec[predictionsDfVec$folds == i, "actualLabel"]
}

combineData = list("prediction"= prediction, 
                   "labels" = labels)

library(ROCR)
pred <- prediction(combineData$prediction , combineData$labels)
perf <- performance(pred,"tpr","fpr")
plot(perf,colorize=TRUE,lty=3)
plot(perf,lwd=3,avg="vertical",
     spread.estimate="stderror",spread.scale=2,
     show.spread.at = c(),
     colorize.palette = rev(rainbow(256, start = 0, end = 4/6)),
     colorkey.relwidth = 0.25, colorkey.pos = "right", print.cutoffs.at = c(), cutoff.label.function = function(x) { round(x, 2) },
     downsampling = 0,
     add = FALSE )



save.image("~/Altoida_VAMBN_paper/workspaces/auc_upd.RData")

