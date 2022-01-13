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
library(pdp)


load("~/Documents/Documents_IT/paper/ALTOIDA_VAMBN_paper/createbnAltoida.RData")
rm(list=setdiff(ls(),  c("targetNodes", "dataCP", "demog", "diagnostics", "cogDomains")))
cogDomains = gsub("_VIS1", "", cogDomains)
cogDomains = gsub("SA_", "", cogDomains)
data = read.csv("altoidaDf.csv")
data = data[data$DX %in% c(0,1,2),]
data$DX[data$DX == 2] = 1
names(data)[names(data) == "DX"] = "Class"
data$X = NULL
mmseFeatures = grep("MMSE", colnames(data), value = TRUE)

func.classifier.sgl.best.param.single = function(dataDf,index, type, Class){
  xtrain = subset(dataDf, select=-c(Class))
  ytrain = dataDf$Class
  ntrain=length(ytrain)
  set.seed(123)
  ##changed for revision Dec 21###
  train.ext = groupKFold(xtrain$SUBJID, k = 5) 
  #train.ext=createFolds(ytrain,k=5,returnTrain=TRUE)
  test.ext=lapply(train.ext,function(x) (1:ntrain)[-x])
  seedVec = c()
  meanRocVec = c()
  lambdasMin = c()
  rownames(xtrain) = paste(xtrain$SUBJID,"_", row.names(xtrain),sep = "")
  xtrain$SUBJID = NULL
  #alpha = 0
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
      #print(xtrain_imp)
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
      if(type == "MMSE"){
        xtrain_imp[,grep("MMSE", colnames(xtrain_imp), value = TRUE)]=lapply(xtrain_imp[,grep("MMSE", colnames(xtrain_imp), value = TRUE)],
                                                                             as.numeric)
        xtest_imp[,grep("MMSE", colnames(xtest_imp), value = TRUE)]=lapply(xtest_imp[,grep("MMSE", colnames(xtest_imp), value = TRUE)],
                                                                           as.numeric)
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
                  min.frac = 0.05, nlam =20, gamma = 0.8, nfold = 10, standardize = TRUE,
                  verbose = FALSE, step = 1, reset = 10, alpha = 1, lambdas = NULL,
                  foldid = foldid)
      #save(fit, file = paste("sgl_R/modelsNew/",type, "_model","_",alpha, "_", i, ".RData", sep = ""))
      minLambda = min(fit$fit$lambdas)
      trainFit = SGL(data, index = index, type = "logit", maxit = 1000, thresh = 0.001,
                     min.frac = 0.1,nlam=2, gamma = 0.8, standardize = TRUE,
                     verbose = FALSE, step = 1, reset = 10, alpha = 1, lambdas=c(min(fit$fit$lambdas),1))
      yh = predictSGL(trainFit, as.matrix(xtest_imp), 1)
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

#for MMSE
##index, for MMSE measures
dataMMSE = data[,c(mmseFeatures, "Class", "SUBJID")]
dataMMSE$MMSE = NULL
index =rep(1,5)
indexMMSE = index
type = "MMSE"
dataMMSE[,grep("MMSE", colnames(dataMMSE), value = TRUE)]=lapply(dataMMSE[,grep("MMSE", colnames(dataMMSE), value = TRUE)],
                                                                                                 as.factor)

paramMMSE = func.classifier.sgl.best.param.single(dataMMSE, index = index, type = type, "Class")
lambdaMMSE = paramMMSE[[1]]
rocMMSE = paramMMSE[[2]]


dataMMSE$SUBJID = NULL
x = subset(dataMMSE, select=-c(Class))
y= dataMMSE$Class
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
               verbose = FALSE, step = 1, reset = 10, alpha = 1, lambdas = c(lambdaMMSE, 1)
)

save(finalFit, file = paste("sgl_R/models_revised/","final_", type, "_model", ".RData", sep = ""))

coef = finalFit$beta[,match(min(finalFit$lambdas),finalFit$lambdas)]
absCoef = abs(coef)
coefDf = cbind.data.frame(colnames(x), coef, absCoef, index)
names(coefDf)[names(coefDf)=="colnames(x)"] = "names"
names(coefDf)[names(coefDf)=="index"] = "group"
coefDf <- coefDf[with(coefDf, order(-absCoef)), ]


library(ggplot2) 
coefDfNoZero = subset(coefDf[coefDf$absCoef>0,])
png("sgl_R/plots_revised/classifier_Altoida_MMSE.png", width = 200, height = 100, units='mm',res = 300) 
ggplot(coefDfNoZero, aes(x=reorder(names, absCoef), weight=absCoef, fill=as.factor(group))) + 
  geom_bar() +
  scale_fill_manual(labels = "MMSE", values = "steelblue", name = "Variable group") +
  ylab("absolute coefficients") +
  xlab("Feature Names")+  theme(axis.text.x = element_text(size=14),
                                axis.text.y = element_text(size=14),
                                axis.title=element_text(size=16))+
  theme(legend.text=element_text(size=12))+
  coord_flip()
dev.off()

#for Cognition
##index, for Cognitive measures 
cogDomains = gsub("SA_", "", cogDomains)
cogDomains = gsub("_VIS1", "", cogDomains)
dataCognition = data[,c(cogDomains, "Class", "SUBJID")]
index = c(rep(1,9))
indexCog = index
Class = "Class"
type = "Cognition"
paramCog = func.classifier.sgl.best.param.single(dataCognition, index = index, type = type, "Class")
lambdaCog = paramCog[[1]]
rocCog = paramCog[[2]]

dataCognition$SUBJID = NULL
x = subset(dataCognition, select=-c(Class))
y= dataCognition$Class
set.seed(123)

x = missForest(x, maxiter = 10, ntree = 100, variablewise = FALSE,
                        decreasing = FALSE, verbose = FALSE,
                        mtry = floor(sqrt(ncol(x))), replace = TRUE,
                        classwt = NULL, cutoff = NULL, strata = NULL,
                        sampsize = NULL, nodesize = NULL, maxnodes = NULL,
                        xtrue = NA)
x = x$ximp
x$X = NULL
finalData = list(x = x, y=y)
finalFit = SGL(finalData, index = index, type = "logit", maxit = 1000, thresh = 0.001,
               min.frac = 0.1, nlam =2, gamma = 0.8,standardize = TRUE,
               verbose = FALSE, step = 1, reset = 10, alpha = 0.95, lambdas = c(lambdaCog,1)
)

save(finalFit, file = paste("sgl_R/models_revised/","final_", type, "_model", ".RData", sep = ""))

coef = finalFit$beta[,match(min(finalFit$lambdas),finalFit$lambdas)]
absCoef = abs(coef)
coefDf = cbind.data.frame(colnames(x), coef, absCoef, index)
names(coefDf)[names(coefDf)=="colnames(x)"] = "names"
names(coefDf)[names(coefDf)=="index"] = "group"
coefDf <- coefDf[with(coefDf, order(-absCoef)), ]


library(ggplot2) 
coefDfNoZero = subset(coefDf[coefDf$absCoef>0,])
png("sgl_R/plots_revised/classifier_Altoida_Cog.png", width = 200, height = 100, units='mm',res = 300) 
ggplot(coefDfNoZero, aes(x=reorder(names, absCoef), weight=absCoef, fill=as.factor(group))) + 
  geom_bar() +
  scale_fill_manual(labels = "Digital cognitive domains", values = "orange", name = "Variable group") +
  ylab("absolute coefficients") +
  xlab("Feature Names")+  theme(axis.text.x = element_text(size=14),
                                 axis.text.y = element_text(size=14))+
  theme(legend.text=element_text(size=12),axis.title=element_text(size=16))+
  coord_flip()
dev.off()

##boxplots auc for final model
auc <- function(data, index,bestAlpha,  bestLambda, seed){
  xtrain = subset(data, select=-c(Class))
  xtrainCp = xtrain
  ytrain = data$Class
  ntrain=length(ytrain)
  set.seed(123)
  # Data split 
  ## change December 2021, not to keep same id in training and testing
  #selected_index = createMultiFolds(ytrain, k=5, times = 10) 
  #train.ext=createFolds(ytrain,k=10, returnTrain=TRUE)
  #test.ext=lapply(train.ext,function(x) (1:ntrain)[-x])
  f1Vec = c()
  rocVec = c()
  for(rep in 1:10){
    print(paste("rep",rep))
    ##MMSE
    # if (rep == 7|rep == 9){
    #   next
    # }
    # Data split 
    ## change December 2021, not to keep same id in training and testing
    set.seed(rep + 123)
    selected_index = groupKFold(xtrainCp$SUBJID, k = 5) 
    print(paste("selected_index",selected_index))
    test.ext=lapply(selected_index,function(x) (1:ntrain)[-x])
  for(i in 1:length(selected_index)){
    print(i)
    #if MMSE##
    if(rep==9 && i ==3){
    next
    }
    ###
    f_index = unlist(selected_index[[i]])
    xtrain$SUBJID = NULL
    rf_train = xtrain[f_index,]
    rf_test = xtrain[-f_index,]
    y_train=ytrain[f_index]
    #print(sum(is.na(rf_train)))
    #print(sum(is.na(rf_test)))
    if(sum(is.na(rf_train))>0){
      set.seed(seed+i)
    xtrain_imp = missForest(rf_train, maxiter = 10, ntree = 100, variablewise = FALSE,
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
    xtest_imp = missForest(rf_test, maxiter = 10, ntree = 100, variablewise = FALSE,
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
    # fit = cvSGL(data, index = index, type = "logit", maxit = 1000, thresh = 0.001,
    #             min.frac = 0.05, nlam =2, gamma = 0.8, nfold = 10, standardize = TRUE,
    #             verbose = FALSE, step = 1, reset = 10, alpha = alphaDigMMSE, lambdas = NULL,
    set.seed(seed+i)
    trainFit = SGL(data, index = index, type = "logit", maxit =1000, thresh = 0.001,
                   min.frac = 0.1,nlam=2, gamma = 0.8, standardize = TRUE,
                   verbose = FALSE, step = 1, reset = 10, alpha = 1, lambdas=c(bestLambda,1))
    
    
    # if(is.na(trainFit) == FALSE){
    
    #   colnames(rocVec)[last_index] <- fold_name
    # }
    yh = predictSGL(trainFit, as.matrix(xtest_imp), 1)
    #print(yh)
    roc = (roc.curve(yh, weights.class0 = ytrain[-f_index], curve = TRUE)$auc)
    #print(roc)
    #rocVec = c(rocVec, roc)
    rocVec = cbind(rocVec, roc)
    #Rename columns
    last_index = ncol(rocVec)
    fold_name = paste("Fold" , as.character(i))
    colnames(rocVec)[last_index] <- fold_name
    #lambdas = c(lambdas, minLambda)
    ##calculating f1 score, change Jan 20222 for the review
    actualLabel =  ytrain[-f_index]
    predictedLabel =  yh
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
  #maxRoc = match(max(rocVec),rocVec)
  #print(maxRoc)
  #bestLambda = lambdas[maxRoc]
  return(list(rocVec, f1Vec))
}

#boxplots##
##MMSE
dataMMSE = data[,c(mmseFeatures, "Class", "SUBJID")]
dataMMSE$MMSE = NULL
index =rep(1,5)
indexMMSE = index
type = "MMSE"
dataMMSE[,grep("MMSE", colnames(dataMMSE), value = TRUE)]=lapply(dataMMSE[,grep("MMSE", colnames(dataMMSE), value = TRUE)],
                                                                 as.factor)


rocRESMMSE = auc(dataMMSE, indexMMSE, bestAlpha = 1, bestLambda = lambdaMMSE, seed = 123)
#rocRESMMSE$auc = NULL
rocRESMMSEdf_mean =  mean(unlist(rocRESMMSE[1]))
boxplot(rocRESMMSE[1], col = "#00BFC4")

f1MMSE = rocRESMMSE[[2]]
f1MMSEdf_mean =  mean(unlist(f1MMSE)) 
boxplot(unlist(rocRESMMSE[2]), col = "#00BFC4")



##Cog
cogDomains = gsub("SA_", "", cogDomains)
cogDomains = gsub("_VIS1", "", cogDomains)
dataCognition = data[,c(cogDomains, "Class", "SUBJID")]

rocRESCog = auc(dataCognition, indexCog, bestAlpha = 1, bestLambda = lambdaCog, seed = 123)
#rocRESCog$auc = NULL
rocCogdf_mean =  mean(unlist(rocRESCog[1]) )
boxplot(unlist(rocRESCog[1]), col = "#00BFC4")

f1Cog = rocRESCog[[2]]
f1Cogdf_mean =  mean(unlist(f1Cog)) 
boxplot(unlist(rocRESCog[2]), col = "#00BFC4")



save.image("~/Altoida_VAMBN_paper/workspaces/classifiersglInd_upd.RData")
