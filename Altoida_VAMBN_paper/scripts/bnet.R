############# README
# This is the main analysis file. 
# 1. run the R files clean_data->format_data->impute_aux (scripts with fixed settings)
# 2. run HI-VAE jupyter notebook (up until VP decoding)
# 3. run full script below for bayesian network and VirtualPatient validation
# 4. decode generated VPs in the HI-VAE notebook
# 5. Additional analyses:
#     - hi-vae_decoded.R (get all the comparison plots with confidence intervals)
#     - bnet_likelihoods.R (get likelihoods of VP/RP under fitted model)
#     - counterfactuals_bnlearn (interventions - decoded plots in HI-VAE notebook)
############# 

############################## Dependencies and helper functions (dependencies for helper functions are loaded here)
rm(list=ls())
library(tidyverse)
library(beepr)
library(arules)
library(mclust)
library(rpart)
library(bnlearn) # hc might be overwritten by arules or some such package "bnlearn::hc" if so; not currently used though
library(parallel)
library(randomForest)
library(randomForestSRC)
library(doMC)

library(binst)

# general helpers
source('helper/plot_bn.R')
source('helper/clean_help.R')
source('helper/simulate_VP.R')
source('helper/VP_missslist.R')
source('helper/save_VPmisslist.R')
# study specific helpers
source('helper/merge_data.R')
source('helper/addnoise.R')
source('helper/make_bl_wl_altoida.R')

############################## Settings and preprocessing

########## Name output files
name<-'main'
data_out<-paste0('data/data_out/',name)

########## Load data & remaining formatting of standalone
data<-merge_data() # could load just basic autoencoded with T (if so need to fully discretize)
baselineDX = data[grep("\\_1", data$SUBJID),]
#save data for Classifier
dataCP = data
##extract the data frame for the testing the logistic regression classifier on MMSE and Digital data 
##for Normal and MCI at Risk for AD


# remove subject variable
pt<-data$SUBJID
data$SUBJID<-NULL

colnames(data) = paste(colnames(data),"_VIS1", sep="")

######### Discretize & set score
discdata<-data
discdata<-addnoise(discdata,0.01)

##remove not required aux columns##
removeCols = c()
for(cols in grep("AUX", colnames(discdata), value= TRUE)){
  if(sum(discdata[,cols]==0)==nrow(discdata)){
    removeCols = c(removeCols, cols)
  }
}
#discdata = discdata[,-which(colnames(discdata) %in% removeCols)]
#discdata<-addnoise(discdata,0.01)## for aux
scr<-"bic-cg" # 'bic' for basic autoencoder and fully discretized
mth<-"mle" # 'bayes' for basic autoencoder and fully discretized

# ######### Add AUX superior
discdata$visitmiss_VIS1<-factor(ifelse(apply(discdata[,grepl('AUX_',colnames(discdata))&grepl('_VIS1',colnames(discdata))],1,function(x) (all(x==1))),1,0))
d1<-as.data.frame(t(discdata[,(grepl('AUX_',colnames(discdata))&grepl('_VIS1',colnames(discdata)))|grepl('visitmiss_VIS1',colnames(discdata))]))
rm1 =c()
for(i in 1:nrow(d1)){
  x = duplicated(d1[i,nrow(d1),])
  print(x)
  if(x == TRUE){
    print(rownames(d1[i,]))
    rm1 = c(rm1,rownames(d1[i,]))
  }
}
 
# # orphaned nodes are those whose immediate parent AUX is identical to the visitmiss at that visit
lowaux<-discdata[,grepl('AUX_',colnames(discdata))&!(colnames(discdata) %in% rm1)]
lowaux<-colnames(lowaux)[sapply(colnames(lowaux),function(x) sum(as.numeric(as.character(lowaux[,x])))<=5)]
discdata<-discdata[ , !(names(discdata) %in% rm1)]
discdata<-discdata[ , !(names(discdata) %in% lowaux)]
orphans<-gsub('AUX_','',rm1)
orphans<-unname(sapply(orphans,function(x) ifelse(!grepl('SA_',x),paste0('zcode_',x),x)))

########## Make bl/wl
blname<-paste0(data_out,'_bl.csv')
wlname<-paste0(data_out,'_wl.csv')
make_bl_wl_altoida(discdata,blname,wlname,F,orphans) # rm has info about orphaned nodes
bl<-read.csv(blname)
wl<-read.csv(wlname)
wlCogDomain = read.csv("wlCogDomain.csv")
wlCogDomain$X = NULL

wlCogDomain = wlCogDomain[wlCogDomain$from %in% colnames(discdata),]

cogDomains = grep("PerceptualMotorCoordination|Planning|ProspectiveMemory|SpatialMemory|CognitiveProcessingSpeed|ComplexAttention|EyeMovement|Flexibility|Inhibition|VisualPerception", colnames(discdata), value = TRUE)
##remove Sppech
wlCogDomain = wlCogDomain[wlCogDomain$to %in% cogDomains,]
wl = rbind.data.frame(wl, wlCogDomain)
############################## Bnet

######### Final bayesian network
#discdata[,grep("^SA_MMSE", colnames(discdata), value = TRUE)] = lapply(discdata[,grep("^SA_MMSE", colnames(discdata), value = TRUE)], as.numeric)
discdata[] <- lapply(discdata, function(x) if(is.integer(x)) as.numeric(x) else x)
saveRDS(discdata,paste0(data_out,'_discdata.rds'))

##test the best  algorithm
set.seed(1234)
registerDoMC(cores=10)
clGn =makeCluster(10)


##depending on the lowest average loss over runs
set.seed(123)
finalBN = hc(discdata,maxp = 5, blacklist = bl,whitelist = wl, score="bic-cg")

saveRDS(finalBN,paste0(data_out,'_finalBN.rds'))

######### Bootstrapped network
cores = detectCores()
cl =  makeCluster(cores)
set.seed(122)
boot.stren = boot.strength(discdata, algorithm="hc", R=1000, algorithm.args = list(maxp=5,blacklist=bl, whitelist=wl, score=scr), cluster=cl)

library(plyr)
stopCluster(cl)
boot.strenwithThreshold = boot.stren[boot.stren$strength >= 0.5 & boot.stren$direction >= 0.5,]
boot.strenwithThreshold["typeEdges"] = rep("solid", nrow(boot.strenwithThreshold))
boot.strenwithThreshold$strength = round(boot.strenwithThreshold$strength, 2)
matchDf1 = as.data.frame(match_df(boot.strenwithThreshold, wlCogDomain, on = c("from", "to")))
matchDf2 = as.data.frame(match_df(boot.strenwithThreshold, wl, on = c("from", "to")))
matchingIDs = rownames(matchDf1,matchDf2)
boot.strenwithThreshold[matchingIDs, "strength"] = NA
boot.strenwithThreshold[matchingIDs, "typeEdges"] = "dashed"
boot.strenwithThreshold$from = gsub('zcode_', '', boot.strenwithThreshold$from)
boot.strenwithThreshold$to = gsub('zcode_', '', boot.strenwithThreshold$to)
boot.strenwithThreshold$from = gsub('SA_', '', boot.strenwithThreshold$from)
boot.strenwithThreshold$to = gsub('SA_', '', boot.strenwithThreshold$to)
boot.strenwithThreshold$from = gsub('_VIS1', '', boot.strenwithThreshold$from)
boot.strenwithThreshold$to = gsub('_VIS1', '', boot.strenwithThreshold$to)
boot.strenwithThreshold = boot.strenwithThreshold[- grep("AUX", boot.strenwithThreshold$from),]


write.csv(boot.strenwithThreshold, "boot_stren_50_percent.csv")
stopCluster(cl)
saveRDS(boot.stren,paste0(data_out,'_bootBN.rds'))
beep()

boot.strenwithThresholdCP = boot.strenwithThreshold
boot.strenwithThresholdCP = boot.stren[boot.stren$strength >= 0.1 & boot.stren$direction >= 0.1,]

boot.strenwithThresholdMMSEtoDM = boot.strenwithThresholdCP[grepl("MMSE", boot.strenwithThresholdCP$from) & grepl("zcode\\_AR|Motor|BIT", boot.strenwithThresholdCP$to),]
boot.strenwithThresholdMMSEtoCogDomain = boot.strenwithThresholdCP[grepl("MMSE", boot.strenwithThresholdCP$from) & boot.strenwithThresholdCP$to %in% cogDomains,]


boot.stren.threshold.visual = rbind.data.frame(boot.strenwithThreshold,boot.strenwithThresholdMMSEtoDM,boot.strenwithThresholdMMSEtoCogDomain)
boot.stren.threshold.visual$from = gsub('zcode_', '', boot.stren.threshold.visual$from)
boot.stren.threshold.visual$to = gsub('zcode_', '', boot.stren.threshold.visual$to)
boot.stren.threshold.visual$from = gsub('SA_', '', boot.stren.threshold.visual$from)
boot.stren.threshold.visual$to = gsub('SA_', '', boot.stren.threshold.visual$to)
boot.stren.threshold.visual$from = gsub('_VIS1', '', boot.stren.threshold.visual$from)
boot.stren.threshold.visual$to = gsub('_VIS1', '', boot.stren.threshold.visual$to)
boot.stren.threshold.visual = boot.stren.threshold.visual[- grep("visitmiss", boot.stren.threshold.visual$from),]
boot.stren.threshold.visual$EdgeType =  ifelse(boot.stren.threshold.visual$strength >= 0.50,"highStrength", "LowStrength")
boot.stren.threshold.visual$strength = round(boot.stren.threshold.visual$strength, 2)
boot.stren.threshold.visual = unique(boot.stren.threshold.visual)
write.csv(boot.stren.threshold.visual, "boot_stren_visual.csv")



# save fitted network
real = discdata

#real$SUBJID<-NULL
finalBN<-readRDS(paste0(data_out,'_finalBN.rds'))
finalBNCp = finalBN
finalBNarcs = as.data.frame(finalBN$arcs)
finalBNarcs$from = gsub('zcode_', '', finalBNarcs$from)
finalBNarcs$to = gsub('zcode_', '', finalBNarcs$to)
finalBNarcs$from = gsub('SA_', '', finalBNarcs$from)
finalBNarcs$to = gsub('SA_', '', finalBNarcs$to)
finalBNarcs$from = gsub('_VIS1', '', finalBNarcs$from)
finalBNarcs$to = gsub('_VIS1', '', finalBNarcs$to)
fitted = bn.fit(finalBN, real, method=mth)
saveRDS(fitted,paste0(data_out,'_finalBN_fitted.rds'))

##assesing the stabiliity of the dynamic bayesina network
png("output_figures/dyn.str.png") 
plot(boot.stren)
dev.off()


#replace.unidentifiable: a boolean value. 
##If TRUE and method is mle, unidentifiable parameters 
##are replaced by zeroes (in the case of regression coefficients 
##and standard errors in Gaussian and conditional Gaussian nodes) 
##or by uniform conditional probabilities (in discrete nodes).
set.seed(234)
fitted = bn.fit(finalBN, discdata,replace.unidentifiable= TRUE)
saveRDS(fitted,paste0(data_out,'_finalBN_fitted.rds'))
set.seed(234)
kfold = bn.cv(discdata, finalBN, "hc", k = 10,loss = "logl-cg")
ensemble = lapply(kfold, `[[`, "fitted")
saveRDS(ensemble,paste0(data_out,'_ensemble.rds'))
finalBNarcs = as.data.frame(finalBN$arcs)
mmseMeasures = grep("^zcode_MMSE", colnames(discdata), value = TRUE)
demog = c("SA_age_VIS1", "SA_gender_VIS1", "SA_dominantHand_VIS1", "SA_yearsOfEducation_VIS1", "SA_Amyloid_VIS1")
diagnostics = "SA_DX_VIS1"

digNodes = setdiff(colnames(discdata), grep("SA|AUX|MMSE|miss|scode", colnames(discdata), value = TRUE))
staloneFromFinalBN = as.data.frame(finalBNarcs[finalBNarcs$from %in% c(diagnostics, mmseMeasures, demog) & 
                                   finalBNarcs$to %in% digNodes,])


staloneFromFinalBN$from = as.character(staloneFromFinalBN$from)
fromNodes = c(unique(unique(staloneFromFinalBN$from), setdiff(unique(staloneFromFinalBN$from), mmseMeasures)))
targetNodes = unique(c(grep("zcode", staloneFromFinalBN$to, value = TRUE)))
targetNodes = unique(targetNodes)
func.nrmse = function(bnCv,feat){
  mse.loss = loss(bnCv)
  nrmse.all = c()
  for(i in mse.loss){
    rmse = sqrt(i)
    print(rmse)
    nrmse = rmse/(max(feat)-min(feat))
    nrmse.all = c(nrmse.all, nrmse)
  }
  #print(nrmse.all)
  print("Mean NRMSE")
  meanNRMSE = round(mean(nrmse.all,na.rm = TRUE),4)
  print("Mean S.D")
  sdNRMSE = round(sd(nrmse.all,na.rm = TRUE),4)
  seNRMSE = round((round(sd(nrmse.all,na.rm = TRUE),4))/sqrt(length(nrmse.all)),4)
  #print(sd(mse.loss))
  return(list(meanNRMSE,sdNRMSE, seNRMSE))
}

library(caret)
bnCvVec = c()
nrmseDf = data.frame()
##not taking Amyloid, as DX and Amyloid are highly correlated
discdataCp = cbind.data.frame(dataCP$SUBJID, discdata)
names(discdataCp)[names(discdataCp) == "dataCP$SUBJID"] = "SUBJID"
vec = discdataCp$SUBJID
discdataCp["ID"] = gsub("_[0-9]*$","",discdataCp$SUBJID)
discdataCp$SUBJID = as.factor(discdataCp$SUBJID)
discdataCp$ID = as.factor(discdataCp$ID)
discdataCp <- fold(discdataCp, k = 10, id_col = 'ID', method = "n_dist")      
discdataCp = as.data.frame(discdataCp)
##to order with the same vector order##
discdataCp = discdataCp[match(vec, discdataCp$SUBJID),]
#rownames(discdataCp) = discdataCp$SUBJID
foldid = as.integer(discdataCp$.folds)
##list of indices for folds
foldListRun = list()
for(runs in 1:10){
  foldList = list()
  for(i in sort(unique(foldid))){
    print(i)
    set.seed(runs+121)
    folds = as.integer(rownames(discdataCp[discdataCp$.folds == i,]))
    foldList[[length(foldList)+1]] = folds
  }
  foldListRun[[length(foldListRun)+1]] = foldList
}

for(i in 1:length(targetNodes)){
  print(targetNodes[i])
  set.seed(123)
  bnCvList = list()
  ##change for paper review
  bnCv = bn.cv(discdata, finalBN,folds= foldListRun,method="custom-folds", loss = "mse-lw-cg",loss.args = list(from = c(mmseMeasures, "SA_age_VIS1", "SA_gender_VIS1", "SA_yearsOfEducation_VIS1", "SA_DX_VIS1"),target= targetNodes[[i]]))
  #bnCv = bn.cv(discdata, finalBN, runs=10,  k = 10, loss = "mse-lw-cg",loss.args = list(from = c(mmseMeasures, "SA_age_VIS1", "SA_gender_VIS1", "SA_yearsOfEducation_VIS1", "SA_DX_VIS1"),target= targetNodes[[i]]))
  nrmse = func.nrmse(bnCv, discdata[,targetNodes[i]])
  print(nrmse)
  #nrmseDf[i,"from"] = toString(fromNodes)
  nrmseDf[i,"target"] = targetNodes[i]
  nrmseDf[i,"Mean"] = nrmse[[1]]
  nrmseDf[i,"SD"] = nrmse[[2]]
  nrmseDf[i,"SE"] = nrmse[[3]]
  bnCvVec = c(bnCvVec, bnCv)
}
nrmseDf = unique(nrmseDf)

##predict cognitive nodes#

bnCvVecCog = c()
nrmseDfCog = data.frame()
for(i in 1:length(cogDomains)){
  print(cogDomains[i])
  set.seed(123)
  ##change for paper review
  bnCv = bn.cv(discdata, finalBN,folds=foldListRun,method="custom-folds", loss = "mse-lw-cg",loss.args = list(from = c(mmseMeasures, "SA_age_VIS1", "SA_gender_VIS1", "SA_yearsOfEducation_VIS1", "SA_DX_VIS1"),target= cogDomains[[i]]))
  #bnCv = bn.cv(discdata, fitted, runs=10,  loss = "mse-lw-cg",loss.args = list(target= targetNodes[i]))
  nrmse = func.nrmse(bnCv, discdata[,cogDomains[i]])
  print(nrmse)
  #nrmseDf[i,"from"] = toString(fromNodes)
  nrmseDfCog[i,"target"] = cogDomains[i]
  nrmseDfCog[i,"Mean"] = nrmse[[1]]
  nrmseDfCog[i,"SD"] = nrmse[[2]]
  nrmseDfCog[i,"SE"] = nrmse[[3]]
  bnCvVecCog = c(bnCvVecCog, bnCv)
}
nrmseDfCog = unique(nrmseDfCog)


#PREDICT MMSE from digital nodes
#target node much be a continous variable
bnCvVecMMSE = c()
nrmseDfMMSE = data.frame()
for(i in 1:length(mmseMeasures)){
  print(mmseMeasures[i])
  set.seed(123)
  ##change for paper review
  bnCvMMSE = bn.cv(discdata, finalBN,folds=foldListRun,method = "custom-folds", loss = "mse-lw-cg",loss.args = list(from = targetNodes,target= mmseMeasures[i]))
  #bnCv = bn.cv(discdata, fitted, runs=10,  loss = "mse-lw-cg",loss.args = list(target= targetNodes[i]))
  nrmseMMSE = func.nrmse(bnCvMMSE, discdata[,mmseMeasures[i]])
  print(nrmseMMSE)
  nrmseDfMMSE[i,"from"] = paste(toString(targetNodes))
  nrmseDfMMSE[i,"target"] = mmseMeasures[i]
  nrmseDfMMSE[i,"Mean"] = nrmseMMSE[[1]]
  nrmseDfMMSE[i,"SD"] = nrmseMMSE[[2]]
  nrmseDfMMSE[i,"SE"] = nrmseMMSE[[3]]
  bnCvVecMMSE = c(bnCvVecMMSE, bnCvMMSE)
}
nrmseDfMMSE = unique(nrmseDfMMSE)
  
save.image("createbnAltoida.RData")


############################## VP vs RP

############################
############################ VP vs RP
############################

# Virtual Patient Generation
#virtual<-simulate_VPs(real,finalBN,fitted,iterative=F,scr,mth,wl,bl)
set.seed(456)
virtual<-simulate_VPs(real,finalBN,iterative=FALSE,scr,mth,wl,bl, n=nrow(real))
############################
############################ save out all data
############################

# save out real and virtual patients
real$SUBJID<-pt
write.csv(real, "real.csv")
real$SUBJID<-NULL
real <- real %>% rename_at(vars(starts_with("zcode")), 
                                 funs(str_replace(., "_VIS1", "")))
real <- real %>% rename_at(vars(starts_with("scode")), 
                                 funs(str_replace(., "_VIS1", "")))
write.csv(real,paste0(data_out,'_RealPPts.csv'),row.names=FALSE)
saveRDS(real,paste0(data_out,'_RealPPts.rds'))

virtual <- virtual %>% rename_at(vars(starts_with("zcode")), 
                             funs(str_replace(., "_VIS1", "")))
virtual <- virtual %>% rename_at(vars(starts_with("scode")), 
                                 funs(str_replace(., "_VIS1", "")))
# save out VP misslist (for HIVAE decoding, tells HIVAE which zcodes the BN considers missing)
data_meta1<-read.csv('data/HI-VAE/metaenc.csv')
data_meta2<-read.csv('data/HI-VAE/metaenc_nm.csv')
#data_meta3<- read.csv('data/HI-VAE/metaenc_mmse.csv')
##missing columns
misCols = setdiff(grep("^zcode|scode", colnames(data_meta1), value = TRUE), "SUBJID")
#misCols3 = setdiff(grep("^zcode|scode", colnames(data_meta3), value = TRUE), "SUBJID")
#misCols = c(misCols1, misCols3)
misCols = intersect(misCols, colnames(virtual))
# for(i in grep("zcode",colnames(virtual),value = TRUE)){
#   names(virtual)[names(virtual)==i] = substr(i,1,nchar(i)-5)
# }
write.csv(virtual,paste0(data_out,'_VirtualPPts.csv'),row.names=FALSE)
saveRDS(virtual,paste0(data_out,'_VirtualPPts.rds'))
save_VPmisslist(virtual[,misCols],'data/HI-VAE/')

cols_nm = intersect(colnames(data_meta2), colnames(virtual))
write.csv(virtual[,cols_nm],paste0(data_out,'_VirtualPPts_nm.csv'),row.names=FALSE)
saveRDS(virtual[,cols_nm],paste0(data_out,'_VirtualPPts_nm.rds'))

