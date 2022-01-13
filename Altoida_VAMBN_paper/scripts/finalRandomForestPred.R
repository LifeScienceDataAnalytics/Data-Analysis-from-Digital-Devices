rm(list=ls())
library(randomForest)
library(ranger)
library(caret)
library(e1071)
library(Metrics)
library(spm)

load("~/Altoida_VAMBN_paper/createbnAltoida.RData")
rm(list=setdiff(ls(),  c("foldListRun","cogDomains","targetNodes", "discdataCp", "discdata","mmseMeasures", "demog", "diagnostics", "fromNodes", "nrmseDf")))
targetNodes = unique(targetNodes)

func.prediction.performance=function(data, target){
  xtrain = data[,setdiff(colnames(data), target)]
  ytrain = data[,target]
  ntrain=length(ytrain)
  #train.ext=createFolds(ytrain,k=10,returnTrain=TRUE)
  #train.ext=groupKFold(xtrain$ID,k=10)
  #test.ext=lapply(train.ext,function(x) (1:ntrain)[-x])
  seed <-7
  #metric<-'RMSE'
  ##changed for keeping the same training and testing folds as in bnlearn, Jan 2021, for paper review
  for(i in 1:10){
    train.ext = unlist(foldListRun[[i]][-i])
    test.ext = unlist(foldListRun[[i]][i])
  }
# Set grid search parameters
  control <- trainControl(method="repeatedcv", number=10, repeats=3, search='grid')

# Outline the grid of parameters
  #tunegrid <- expand.grid(.maxnodes=c(70,80,90,100), .ntree=c(900, 1000, 1100))
  splitrule = c("variance", "extratrees","maxstat")
  tunegrid <- expand.grid(
    .mtry = 2:4,
    .splitrule = splitrule,
    .min.node.size = c(10, 20)
  )
  set.seed(seed)
  xtrain$ID = NULL
  rf_gridsearch <- train(x=xtrain, y=ytrain, method='ranger', tuneGrid=tunegrid, num.trees = 100,  trControl=control)

  finalModel = rf_gridsearch$finalModel
  data<-data[complete.cases(data),]
  rmseAll = c()
  nrmseAll = c()
  ##best hyperparameters, number of trees = 900, No. of variables tried at each split: 4
  for (i in 1:10){
    set.seed(121+i)
    model<-train(#x = xtrain[train.ext[[i]],],y = ytrain[train.ext[[i]]],
      ##changed for keeping the same training and testing folds as in bnlearn, Jan 2021, for paper review
                 x = xtrain[train.ext,],y = ytrain[train.ext],
               method = 'ranger',tuneGrid = data.frame(mtry=rf_gridsearch$bestTune$mtry, min.node.size = rf_gridsearch$bestTune$min.node.size,
                                                       splitrule = rf_gridsearch$bestTune$splitrule),
               trControl=trainControl(method="repeatedcv", number=10, repeats=3))
    #testPred <- predict(model , xtrain[test.ext[[i]],])
    testPred <- predict(model , xtrain[test.ext,])
    print(testPred)
    #rmse = caret::postResample(testPred ,  ytrain[test.ext[[i]]])['RMSE']^2
    #nrmse = rmse/(max(ytrain[test.ext[[i]]])-min(ytrain[test.ext[[i]]]))
    rmse = caret::postResample(testPred ,  ytrain[test.ext])['RMSE']^2
    nrmse = rmse/(max(ytrain[test.ext])-min(ytrain[test.ext]))
    print(nrmse)
    rmseAll = c(rmseAll, rmse)
    nrmseAll = c(nrmseAll, nrmse)
  }
  print(nrmseAll)
  rmseMean = mean(rmseAll)
  nrmseMean = mean(nrmseAll)
  nrmseMean = round(nrmseMean,4)
  sdNRMSE = round(sd(nrmseAll,na.rm = TRUE),4)
  seNRMSE = round(sdNRMSE/sqrt(length(nrmseAll)),4)
  saveRDS(finalModel, paste("~/Documents/Documents_IT/paper/Altoida_VAMBN_paper/finalModels/",target, ".rds", sep = ""))
  return(list(rmseMean,nrmseMean,sdNRMSE,seNRMSE))
}

##prediction of digital nodes
regResults = data.frame()
for(i in 1:length(targetNodes)){
  print(targetNodes[i])
  data = discdataCp[,c(mmseMeasures, "SA_age_VIS1", "SA_gender_VIS1", "SA_yearsOfEducation_VIS1", "SA_DX_VIS1", targetNodes[i], "ID")]
  results = func.prediction.performance(data,targetNodes[i])
  regResults[i,"target"] = targetNodes[i]
  #regResults[i,"maeMean"] = results[1]
  regResults[i,"rmseMean"] = results[1]
  regResults[i,"nrmseMean"] = results[2]
  regResults[i,"sdNRMSE"] = results[3]
  regResults[i,"seNRMSE"] = results[4]
}

write.csv(regResults, "regResults.csv")

##prediction of cognitive nodes
regResultsCog = data.frame()
for(i in 1:length(cogDomains)){
  print(targetNodes[i])
  data = discdataCp[,c(mmseMeasures, "SA_age_VIS1", "SA_gender_VIS1", "SA_yearsOfEducation_VIS1", "SA_DX_VIS1", cogDomains[i], "ID")]
  results = func.prediction.performance(data,cogDomains[i])
  regResultsCog[i,"target"] = cogDomains[i]
  #regResults[i,"maeMean"] = results[1]
  regResultsCog[i,"rmseMean"] = results[1]
  regResultsCog[i,"nrmseMean"] = results[2]
  regResultsCog[i,"sdNRMSE"] = results[3]
  regResultsCog[i,"seNRMSE"] = results[4]
}

write.csv(regResultsCog, "regResultsCog.csv")


##predict MMSE nodes from digital nodes
digitalNodes = targetNodes 
regMMSE = data.frame()
for(i in 1:length(mmseMeasures)){
  print(mmseMeasures[i])
  data = discdataCp[,c(digitalNodes,mmseMeasures[i], "ID")]
  results = func.prediction.performance(data,mmseMeasures[i])
  regMMSE[i,"from"] = paste(toString(targetNodes))
  regMMSE[i,"target"] = mmseMeasures[i]
  regMMSE[i,"rmseMean"] =results[1]
  regMMSE[i,"nrmseMean"] =results[2]
  regMMSE[i,"sdNRMSE"] =results[3]
  regMMSE[i,"seNRMSE"] =results[4]
}


write.csv(regMMSE, "regMMSE.csv")

