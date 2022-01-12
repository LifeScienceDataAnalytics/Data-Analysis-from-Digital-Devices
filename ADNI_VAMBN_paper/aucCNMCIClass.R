rm(list=ls())
library(caret)
library(SGL)
library(PRROC)
library(stringr)
##load the workspace saved from the scripts that are run on cluster 
load("~/ADNI_VAMBN_paper/workspaces/resultsSGLCNMCI_upd.RData")
auc <- function(data, index,bestAlpha, bestLambda, seed){
  xtrain = subset(data, select=-c(Class))
  ytrain = data$Class
  ntrain=length(ytrain)
  set.seed(123)
  selected_index = createMultiFolds(ytrain, k=5, times = 10) 
  rocVec = c()
  f1Vec = c()
  for(i in 1:length(selected_index)){
    print(i)
    f_index = unlist(selected_index[[i]])
    rf_train = xtrain[f_index,]
    rf_test = xtrain[-f_index,]
    y_train=ytrain[f_index]
    data = list(x=rf_train, y=y_train)
    trainFit = SGL(data, index = index, type = "logit", maxit = 1000, thresh = 0.001,
                   min.frac = 0.1,nlam=2, gamma = 0.8, standardize = TRUE,
                   verbose = FALSE, step = 1, reset = 10, alpha = bestAlpha, lambdas=c(bestLambda,1))
    yh = predictSGL(trainFit,as.matrix(rf_test), 1)
    print(yh)
    actualLabel =  ytrain[-f_index]
    predictedLabel =  yh
    roc = (roc.curve(yh, weights.class0 = ytrain[-f_index], curve = TRUE)$auc)
    print(roc)
    rocVec = cbind(rocVec, roc)
    ##calculating f1 score, change Jan 2022 for the review
    predictedLabelUpd = ifelse(predictedLabel <0.5,0,1)
    predictedLabelUpd = as.factor(as.vector(predictedLabelUpd))
    actualLabel = as.factor(actualLabel)
    precision <- posPredValue(predictedLabelUpd, actualLabel, positive="0")
    recall <- sensitivity(predictedLabelUpd, actualLabel, positive="0")
    f1 = 2*((precision*recall)/(precision+recall))
    f1Vec = cbind(f1Vec, f1)
    #Rename columns
    last_index = ncol(rocVec)
    fold_name = paste("Fold" , as.character(i))
    colnames(rocVec)[last_index] <- fold_name
  }
  return(list(rocVec, f1Vec))
}

##Calculating AUC and F1

##MMSE
rocMMSE = auc(dataMMSE, indexMMSE, bestAlpha = 1, bestLambda = bestLambda, seed = 123)
rocMMSE$auc = NULL
rocMMSEdf_mean =  mean(unlist(rocMMSE[[1]]) ) 
boxplot(unlist(rocMMSE[1]), col = "#00BFC4")
f1MMSEdf_mean =  mean(unlist(rocMMSE[[2]])) 
boxplot(unlist(rocMMSE[2]), col = "#00BFC4")


##FAQ
rocFAQ = auc(dataFAQ, indexFAQ, bestAlpha = 1, bestLambda = lambdaFAQ, seed = 123)
rocFAQ$auc = NULL
rocFAQdf_mean =  mean(unlist(rocFAQ[[1]]) ) 
boxplot(unlist(rocFAQ[1]), col = "#00BFC4")
f1FAQdf_mean =  mean(unlist(rocFAQ[[2]])) 
boxplot(unlist(rocFAQ[2]), col = "#00BFC4")

##Dig
rocDig = auc(dataDig, indexDig, bestAlpha = alphaDig, bestLambda = lambdaDig, seed = 123)
rocDig$auc = NULL
rocDigdf_mean =  mean(unlist(rocDig) ) 
boxplot(unlist(rocDig), col = "#00BFC4")

##Cog
rocCog = auc(dataCognition, indexCog, bestAlpha = 1, bestLambda = lambdaCog, seed = 123)
rocCog$auc = NULL
rocCogdf_mean =  mean(unlist(rocCog) ) 
boxplot(unlist(rocCog), col = "#00BFC4")

##CogMMSEFAQ
rocCogMMSEFAQ = auc(adniClassifierCNtoMCICogMMSEFAQ, indexCogMMSEFAQ, bestAlpha = alphaCogMMSEFAQ, 
                    bestLambda = lambdaCogMMSEFAQ, seed = 123)
rocCogMMSEFAQ$auc = NULL
rocCogMMSEFAQdf_mean =  mean(unlist(rocCogMMSEFAQ[[1]]) ) 
boxplot(unlist(rocCogMMSEFAQ[1]), col = "#00BFC4")

##CogMMSE
rocCogMMSE = auc(adniClassifierCNtoMCICogMMSE, indexCogMMSE, bestAlpha = alphaCogMMSE, 
                    bestLambda = lambdaCogMMSE, seed = 123)
rocCogMMSE$auc = NULL
rocCogMMSEdf_mean =  mean(unlist(rocCogMMSE) ) 
boxplot(unlist(rocCogMMSE), col = "#00BFC4")

##MMSEFAQ
rocFAQMMSE = auc(adniClassifierCNtoMCIFAQMMSE, indexFAQMMSE, bestAlpha = alphaFAQMMSE, 
                 bestLambda = lambdaFAQMMSE, seed = 123)
rocFAQMMSE$auc = NULL
rocFAQMMSEdf_mean =  mean(unlist(rocFAQMMSE[[1]]) ) 
boxplot(unlist(rocFAQMMSE[1]), col = "#00BFC4")
f1FAQMMSEdf_mean =  mean(unlist(rocFAQMMSE[[2]])) 
boxplot(unlist(rocFAQMMSE[2]), col = "#00BFC4")

##DigMMSEFAQ
rocDigMMSEFAQ = auc(adniClassifierCNtoMCIDigMMSEFAQ, indexDigMMSEFAQ, bestAlpha = alphaDigMMSEFAQ, 
                    bestLambda = lambdaDigMMSEFAQ, seed = 123)
rocDigMMSEFAQ$auc = NULL
rocDigMMSEFAQdf_mean =  mean(unlist(rocDigMMSEFAQ) ) 
boxplot(unlist(rocDigMMSEFAQ), col = "#00BFC4")

boxplotDf = cbind.data.frame(unlist(rocMMSE), unlist(rocFAQ),unlist(rocFAQMMSE), unlist(rocCog), unlist(rocCogMMSEFAQ), unlist(rocDigMMSEFAQ))

boxplotDfUpdAUC = cbind.data.frame(unlist(rocMMSE[1]), unlist(rocFAQ[1]),unlist(rocFAQMMSE[1]))

boxplotDfUpdF1 = cbind.data.frame(unlist(rocMMSE[2]), unlist(rocFAQ[2]),unlist(rocFAQMMSE[2]))


names(boxplotDfUpdAUC) = c( "MMSE", "FAQ",  "MMSE, FAQ")
names(boxplotDfUpdF1) = c( "MMSE", "FAQ",  "MMSE, FAQ")
fill <- "#00BFC4"
line <- "#1F3552"

##auc
boxplotDfUpdAUC = reshape2::melt(boxplotDfUpdAUC)
png("~/ADNI_VAMBN_paper/sgl_R/plots_revised/AUC_boxplot_adni_cnmci_upd.png", width = 200, height = 150, units='mm',res = 300) 
plt <- ggplot(data = boxplotDfUpdAUC, aes(x = variable, y = value))
plt + geom_boxplot(fill = fill, colour = line) + labs(x = "Features", y = "AUC") + ylim(0.5,1) + theme(axis.text = element_text(size = 14)) +theme(text = element_text(size = 18))+
scale_x_discrete(labels = function(x) str_wrap(x, width = 10)) + theme_classic(base_size = 15)
dev.off()

##f1
boxplotDfUpdF1 = reshape2::melt(boxplotDfUpdF1)
png("~/ADNI_VAMBN_paper/sgl_R/plots_revised/F1_boxplot_adni_cnmci_upd.png", width = 200, height = 150, units='mm',res = 300) 
plt <- ggplot(data = boxplotDfUpdF1, aes(x = variable, y = value))
plt + geom_boxplot(fill = fill, colour = line) + labs(x = "Features", y = "F1") + ylim(0.5,1) + theme(axis.text = element_text(size = 14)) +theme(text = element_text(size = 18))+
  scale_x_discrete(labels = function(x) str_wrap(x, width = 10)) + theme_classic(base_size = 15)
dev.off()

##Wilcoxon tests to test the significance ##two sided 
mmseWithCog1 = wilcox.test(unlist(rocMMSE), unlist(rocCogMMSEFAQ), paired = TRUE, alternative = "two.sided")
mmseWithCog2 = wilcox.test(unlist(rocFAQ), unlist(rocCogMMSEFAQ), paired = TRUE, alternative = "two.sided")
mmseWithDig1 = wilcox.test(unlist(rocMMSE), unlist(rocDigMMSEFAQ), paired = TRUE, alternative = "two.sided")
mmseWithDig2 = wilcox.test(unlist(rocFAQ), unlist(rocDigMMSEFAQ), paired = TRUE, alternative = "two.sided")
mmseAndCog = wilcox.test(unlist(rocMMSE), unlist(rocCog), paired = TRUE, alternative = "two.sided")
FAQAndCog = wilcox.test(unlist(rocFAQ), unlist(rocCog), paired = TRUE, alternative = "two.sided")
mmseAndDig = wilcox.test(unlist(rocMMSE), unlist(rocDig), paired = TRUE, alternative = "two.sided")
FAQAndDig = wilcox.test(unlist(rocFAQ), unlist(rocDig), paired = TRUE, alternative = "two.sided")
mmseAndFAQ = wilcox.test(unlist(rocFAQ), unlist(rocMMSE), paired = TRUE, alternative = "two.sided")
mmseAndmmseAndFAQ = wilcox.test(unlist(rocFAQ), unlist(rocFAQMMSE), paired = TRUE, alternative = "two.sided")

save.image("~/ADNI_VAMBN_paper_final/workspaces/aucSGLCNMCI.RData")

