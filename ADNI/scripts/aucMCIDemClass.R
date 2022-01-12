load("~/Documents/Documents_IT/paper/ADNI_VAMBN_paper/resultsSGLMCIDem.RData")
#load("~/Documents/Documents_IT/paper/ADNI_VAMBN_paper/resultsSGLMCIDem_upd.RData")
auc <- function(data, index,bestAlpha, bestLambda, seed){
  xtrain = subset(data, select=-c(Class))
  ytrain = data$Class
  ntrain=length(ytrain)
  set.seed(123)
  selected_index = createMultiFolds(ytrain, k=5, times = 10) 
  #train.ext=createFolds(ytrain,k=10,returnTrain=TRUE)
  #test.ext=lapply(train.ext,function(x) (1:ntrain)[-x])
  rocVec = c()
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
    roc = (roc.curve(yh, weights.class0 = ytrain[-f_index], curve = TRUE)$auc)
    print(roc)
    #rocVec = c(rocVec, roc)
    rocVec = cbind(rocVec, roc)
    #Rename columns
    last_index = ncol(rocVec)
    fold_name = paste("Fold" , as.character(i))
    colnames(rocVec)[last_index] <- fold_name
  }
  return(rocVec)
}



##FAQ
rocFAQ = auc(dataFAQ, indexFAQ, bestAlpha = 1, bestLambda = lambdaFAQ, seed = 123)
rocFAQ$auc = NULL
rocFAQdf_mean =  mean(unlist(rocFAQ) ) 
boxplot(unlist(rocFAQ), col = "#00BFC4")

##Cog
rocCog = auc(dataCognition, indexCog, bestAlpha = 1, bestLambda = lambdaCog, seed = 123)
rocCog$auc = NULL
rocCogdf_mean =  mean(unlist(rocCog) ) 
boxplot(unlist(rocCog), col = "#00BFC4")

##Dig
rocDig = auc(dataDig, indexDig, bestAlpha = alphaDig, bestLambda = lambdaDig, seed = 123)
rocDig_upd = rocDig
rocDig$auc = NULL
rocDigdf_mean =  mean(unlist(rocDig) ) 
boxplot(unlist(rocDig), col = "#00BFC4")

##CogFAQ
adniClassifierMCItoDemCogFAQ[,grep("FAQ", colnames(adniClassifierMCItoDemCogFAQ), value = TRUE)]=lapply(adniClassifierMCItoDemCogFAQ[,grep("FAQ", colnames(adniClassifierMCItoDemCogFAQ), value = TRUE)],
        as.numeric)
rocCogFAQ = auc(adniClassifierMCItoDemCogFAQ, indexCogFAQ, bestAlpha = alphaCogFAQ, 
                    bestLambda = lambdaCogFAQ, seed = 123)
rocCogFAQ $auc = NULL
rocCogFAQdf_mean =  mean(unlist(rocCogFAQ) ) 
boxplot(unlist(rocCogFAQ), col = "#00BFC4")

rocDigFAQ = auc(adniClassifierMCItoDementiaDigFAQ, indexDigFAQ, bestAlpha = alphaDigFAQ, 
                bestLambda = lambdaDigFAQ, seed = 123)
rocDigFAQ$auc = NULL
rocDigFAQdf_mean =  mean(unlist(rocCogFAQ) ) 
boxplot(unlist(rocDigFAQ), col = "#00BFC4")

boxplotDf = cbind.data.frame(unlist(rocFAQ), unlist(rocCog), unlist(rocCogFAQ), unlist(rocDigFAQ))


names(boxplotDf) = c("FAQ", "Digital cognitive domains", "FAQ and Digital cognitive domains",
                     "FAQ and Digital tasks")
fill <- "#00BFC4"
line <- "#1F3552"

boxplotDf = reshape2::melt(boxplotDf)
png("sgl_R/plots/AUC_boxplot_adni_mciDem.png", width = 300, height = 200, units='mm',res = 300) 
plt <- ggplot(data = boxplotDf, aes(x = variable, y = value))
plt + geom_boxplot(fill = fill, colour = line) + labs(x = "Features", y = "AUC") + ylim(0.5,1) + theme(axis.text = element_text(size = 14)) +theme(text = element_text(size = 18))+
scale_x_discrete(labels = function(x) str_wrap(x, width = 10))+ theme_classic(base_size = 15)
dev.off()
mmseWithCog1 = wilcox.test(unlist(rocFAQ), unlist(rocCogFAQ), paired = TRUE, alternative = "two.sided")
mmseWithDig1 = wilcox.test(unlist(rocFAQ), unlist(rocDigFAQ), paired = TRUE, alternative = "two.sided")
save.image("~/Documents/Documents_IT/paper/ADNI_VAMBN_paper/aucSGLMCIDem.RData")  

