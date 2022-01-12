rm(list=ls())
library(tidyverse)
library(bnlearn)
##load the workspace from  the model trained on ALTOIDA data##
load("~/Documents/Documents_IT/paper/Altoida_VAMBN_paper/createbnAltoida.RData")
#load("~/Documents/Documents_IT/VAMBN/randomForestPred_Nov26.RData")
##keep the required target nodes that needs to be simulated in ADNI data#
rm(list=setdiff(ls(),  c("targetNodes", "nrmseDf", "regResultsFinal", "discdata", "cogDomains",
                         "mmseMeasures", "fromNodes", "nrmseDfMMSE", "finalBNarcs")))
#keep(targetNodes, sure=TRUE)
cogDomains = unique(cogDomains)

##set the working directory to the ALTOIDA model working directory
setwd("~/Documents/Documents_IT/paper/Altoida_VAMBN_paper/")
#load the model trained on ALTOIDA data#
finalBN<-readRDS('data/data_out/main_finalBN.rds')# load finalBN
#ensembleBN<-readRDS('data/data_out/main_ensemble.rds')# load finalBN

# decRP_m<-read.csv('data/HI-VAE/reconRP.csv')
# decRP_nm<-read.csv('data/HI-VAE/reconRP_nm.csv')
# decRP = merge(decRP_m, decRP_nm)
# real<-read.csv("real.csv")# load full real dataset 
# real <- merge(real[,grep("SA", colnames(real), value= TRUE)], decRP[,c("")])
# real$SUBJID<-NULL
fitted.bnlearn<-readRDS('data/data_out/main_finalBN_fitted.rds')
tmp<-readRDS('data/data_out/data_all_imp.rds')# load demo dataset

parentsNodestobeSubstituted = as.data.frame(finalBNarcs[finalBNarcs$from %in% c(targetNodes,mmseMeasures,"SA_Amyloid_VIS1", "SA_age_VIS1","SA_gender_VIS1", "SA_yearsOfEducation_VIS1", "SA_DX_VIS1") & 
                                                          finalBNarcs$to %in% cogDomains,])

parentsNodestobeSubstituted = unique(parentsNodestobeSubstituted)
rownames(parentsNodestobeSubstituted) = 1:nrow(parentsNodestobeSubstituted)
cogDomainsUpd = as.character(unique(parentsNodestobeSubstituted$to))

#digNodes = targetNodes
#simulate the digital nodes for ADNI##
#data = predData1
#target = cogDomains[1]
#visit = 1
func.pred.cog.nodes = function(data, target, visit){
  names(data)[names(data)%in%c("SA_PTGENDER_VIS1","SA_PTGENDER_VIS2","SA_PTGENDER_VIS6","SA_PTGENDER_VIS12","SA_PTGENDER_VIS18","SA_PTGENDER_VIS24","SA_PTGENDER_VIS36")] = "SA_gender_VIS1"
  names(data)[names(data)%in%c("SA_AGE_VIS1","SA_AGE_VIS2","SA_AGE_VIS6","SA_AGE_VIS12","SA_AGE_VIS18","SA_AGE_VIS24","SA_AGE_VIS36")] = "SA_age_VIS1"
  names(data)[names(data)%in%c("SA_PTEDUCAT_VIS1","SA_PTEDUCAT_VIS6","SA_PTEDUCAT_VIS12","SA_PTEDUCAT_VIS18","SA_PTEDUCAT_VIS24","SA_PTEDUCAT_VIS36")] = "SA_yearsOfEducation_VIS1"
  names(data)[names(data)%in%c("SA_Amyloid_VIS1","SA_Amyloid_VIS6","SA_Amyloid_VIS12","SA_Amyloid_VIS18","SA_Amyloid_VIS24","SA_Amyloid_VIS36")] = "SA_Amyloid_VIS1"
  
  #data$SA_gender_VIS1 = as.factor(data$SA_gender_VIS1)
  data$SA_gender_VIS1 = as.character(data$SA_gender_VIS1)
  data$SA_gender_VIS1[data$SA_gender_VIS1=="Female"] <- "2"
  data$SA_gender_VIS1[data$SA_gender_VIS1=="Male"] <- "1"
  data$SA_gender_VIS1 = factor(data$SA_gender_VIS1)
  #data$SA_Amyloid_VIS1 = factor(data$SA_Amyloid_VIS1)
  dxVisit = paste("SA_DX_VIS", visit, sep = "")
  print(dxVisit)
  if("CN" %in% unique(data[,dxVisit])){
    levels(data[,dxVisit])[levels(data[,dxVisit])=="CN"] <- "0"
  }
  if("MCI" %in% unique(data[,dxVisit])){
    levels(data[,dxVisit])[levels(data[,dxVisit])=="MCI"] <- "1"
  }
  if("Dementia" %in% unique(data[,dxVisit])){
    levels(data[,dxVisit])[levels(data[,dxVisit])=="Dementia"] <- "4"
  }
  data[,dxVisit] = factor(data[,dxVisit])
  if(visit >1){
    print(visit)
    data = data %>%
      setNames(gsub("[0-9]+", 1, names(.)))
    dxVisit = grep("DX", colnames(data), value = TRUE)
  }
  data$SA_gender_VIS1 = factor(data$SA_gender_VIS1,levels = levels(discdata$SA_gender_VIS1))
  data$SA_Amyloid_VIS1 = factor(data$SA_Amyloid_VIS1,levels = levels(discdata$SA_Amyloid_VIS1))
  data[,dxVisit] = factor(data[,dxVisit],levels = levels(discdata$SA_DX_VIS1))
  #mmseCommon = intersect(grep("MMSE",colnames(data), value=TRUE), grep("MMSE",fromNodes, value = TRUE))
  #common = intersect(colnames(data), c(mmseMeasures, "SA_age_VIS1","SA_gender_VIS1", "SA_yearsOfEducation_VIS1", "SA_DX_VIS1"))
  common = c(mmseMeasures,  "SA_age_VIS1", "SA_gender_VIS1", "SA_yearsOfEducation_VIS1", "SA_DX_VIS1")
  #print("mmse common")
  #print(mmseCommon)
  #parentNodes = parentsNodestobeSubstituted[parentsNodestobeSubstituted$to == target, "from"]
  #parentNodes = intersect(c("SA_gender_VIS1","SA_age_VIS1", dxVisit,mmseCommon), fromNodes)
  #print(parentNodes)
  pt = data["SUBJID"]
  data = data[,common]
  data = as.data.frame(data)
  set.seed(121)
  print("colnames data")
  print(colnames(data))
  if(ncol(data) ==1){
    names(data)[names(data) == "data"] = parentNodes
  }
  #pred.ensemble = sapply(ensembleBN, predict, node = target, data = data, method = "bayes-lw")
  #data[target]=rowMeans(pred.ensemble)
  data[target] = predict(fitted.bnlearn, data = data, node = target,method = "bayes-lw", debug = FALSE)
  return(list(data,pt))
}

func.all.time.points = function(visData, visit){
  digDataVis = c()
  for(target in cogDomainsUpd){
    print(target)
    data = func.pred.cog.nodes(visData, target, visit)
    #print(data)
    digDataVis = c(digDataVis, data[[1]][target])
  }
  digDataVis = as.data.frame(digDataVis)
  #print(digDataVis)
  digDataVis = cbind.data.frame(digDataVis, data[[2]])
  #colnames to be changed
  x = setdiff(colnames(digDataVis), "SUBJID")
  digDataVis <- digDataVis %>%
    rename_with(~ gsub('VIS[0-9]+', paste("VIS",visit,sep =""), .x))
  #digDataVis = digDataVis %>% rename_at(vars(-SUBJID), ~ gsub("VIS[0-9]+", paste("VIS",visit,sep ="")))
  return(digDataVis)
}



##simulate the digital nodes at each time point##
data_meta<-read.csv('~/Documents/Documents_IT/paper/ADNI_VAMBN_paper_final/data/HI-VAE/metaenc.csv')
data_meta_nm<-read.csv('~/Documents/Documents_IT/paper/ADNI_VAMBN_paper_final/data/HI-VAE/metaenc_nm.csv')
data_meta = merge(data_meta, data_meta_nm)
impADNIdata = readRDS("~/Documents/Documents_IT/paper/ADNI_VAMBN_paper_final/data/data_out/data_all_imp.rds")
inputDataVIS1_dx = impADNIdata$stalone_VIS1_dx
inputDataVIS1_demog = impADNIdata$stalone_VIS1_demog
inputDataVIS1 = merge(inputDataVIS1_dx, inputDataVIS1_demog)
inputDataVIS1 = merge(data_meta[,c(grep("zcode\\_MMSE.*VIS1$", colnames(data_meta), value=TRUE), "SUBJID")], inputDataVIS1)
#data_dig_meta_vis1 = read.csv("~/Documents/Documents_IT/paper/ADNI_VAMBN_paper_final/data/HI-VAE/digitalDataPredicted/data_dig_meta_vis1.csv")
#data_dig_meta_vis1$X = NULL
#predData1 = merge(inputDataVIS1, data_dig_meta_vis1)
cogVIS1 = func.all.time.points(inputDataVIS1,1)


inputDataVIS6_dx = impADNIdata$stalone_VIS6_dx
inputDataVIS6_demog = impADNIdata$stalone_VIS1_demog
inputDataVIS6 = merge(inputDataVIS6_dx, inputDataVIS1_demog)
inputDataVIS6 = merge(data_meta[,c(grep("zcode\\_MMSE.*VIS6", colnames(data_meta), value=TRUE), "SUBJID")], inputDataVIS6)
inputDataVIS6["SA_PTGENDER_VIS6"] = inputDataVIS1$SA_PTGENDER_VIS1 
inputDataVIS6["SA_AGE_VIS6"] = inputDataVIS1$SA_AGE_VIS1 
inputDataVIS6["SA_PTEDUCAT_VIS6"] = inputDataVIS1$SA_PTEDUCAT_VIS1
#data_dig_meta_vis6 = read.csv("~/Documents/Documents_IT/paper/ADNI_VAMBN_paper_final/data/HI-VAE/digitalDataPredicted/data_dig_meta_vis6.csv")
#data_dig_meta_vis6$X = NULL
#predData6 = merge(inputDataVIS6, data_dig_meta_vis6)
cogVIS6 = func.all.time.points(inputDataVIS6,6)

inputDataVIS12_dx = impADNIdata$stalone_VIS12_dx
inputDataVIS12_demog = impADNIdata$stalone_VIS1_demog
inputDataVIS12 = merge(inputDataVIS12_dx, inputDataVIS1_demog)
inputDataVIS12 = merge(data_meta, inputDataVIS12)
inputDataVIS12["SA_PTGENDER_VIS12"] = inputDataVIS1$SA_PTGENDER_VIS1 
inputDataVIS12["SA_AGE_VIS12"] = inputDataVIS1$SA_AGE_VIS1 
inputDataVIS12["SA_PTEDUCAT_VIS12"] = inputDataVIS1$SA_PTEDUCAT_VIS1
#data_dig_meta_vis12 = read.csv("~/Documents/Documents_IT/paper/ADNI_VAMBN_paper_final/data/HI-VAE/digitalDataPredicted/data_dig_meta_vis12.csv")
#data_dig_meta_vis12$X = NULL
#predData12 = merge(inputDataVIS12, data_dig_meta_vis12)
cogVIS12 = func.all.time.points(inputDataVIS12,12)


inputDataVIS24_dx = impADNIdata$stalone_VIS24_dx
inputDataVIS24_demog = impADNIdata$stalone_VIS1_demog
inputDataVIS24 = merge(inputDataVIS24_dx, inputDataVIS1_demog)
inputDataVIS24 = merge(data_meta, inputDataVIS24)
inputDataVIS24["SA_PTGENDER_VIS24"] = inputDataVIS1$SA_PTGENDER_VIS1 
inputDataVIS24["SA_AGE_VIS24"] = inputDataVIS1$SA_AGE_VIS1 
inputDataVIS24["SA_PTEDUCAT_VIS24"] = inputDataVIS1$SA_PTEDUCAT_VIS1
#data_dig_meta_vis24 = read.csv("~/Documents/Documents_IT/paper/ADNI_VAMBN_paper_final/data/HI-VAE/digitalDataPredicted/data_dig_meta_vis24.csv")
#data_dig_meta_vis24$X = NULL
#predData24 = merge(inputDataVIS24, data_dig_meta_vis24)
cogVIS24 = func.all.time.points(inputDataVIS24,24)

inputDataVIS36_dx = impADNIdata$stalone_VIS36_dx
inputDataVIS36_demog = impADNIdata$stalone_VIS1_demog
inputDataVIS36 = merge(inputDataVIS36_dx, inputDataVIS36_demog)
inputDataVIS36 = merge(data_meta, inputDataVIS36)
inputDataVIS36["SA_PTGENDER_VIS36"] = inputDataVIS1$SA_PTGENDER_VIS1
inputDataVIS36["SA_AGE_VIS36"] = inputDataVIS1$SA_AGE_VIS1 
inputDataVIS36["SA_PTEDUCAT_VIS36"] = inputDataVIS1$SA_PTEDUCAT_VIS1
#data_dig_meta_vis36 = read.csv("~/Documents/Documents_IT/paper/ADNI_VAMBN_paper_final/data/HI-VAE/digitalDataPredicted/data_dig_meta_vis36.csv")
#data_dig_meta_vis36$X = NULL
#predData36 = merge(inputDataVIS36, data_dig_meta_vis36)
cogVIS36 = func.all.time.points(inputDataVIS36,36)

cogDataAll = Reduce(function(x, y) merge(x, y, all=TRUE), list(cogVIS1, cogVIS6, cogVIS12, cogVIS24, 
                                                               cogVIS36))

cogDataAll= cogDataAll[!duplicated(as.list(cogDataAll))]
cogDataAll$SUBJID

##add the digital data to data_all
data_all = readRDS("~/Documents/Documents_IT/paper/ADNI_VAMBN_paper_final/data/data_out/data_all_imp.rds")
#cogVIS1 = cogVIS1 %>% rename_at(vars(-SUBJID), ~ paste0("SA_",.))
data_all[["stalone_cog_VIS1"]] = cogVIS1

#cogVIS6 = cogVIS6 %>% rename_at(vars(-SUBJID), ~ paste0("SA_",.))
data_all[["stalone_cog_VIS6"]] = cogVIS6

#cogVIS12 = cogVIS12 %>% rename_at(vars(-SUBJID), ~ paste0("SA_",.))
data_all[["stalone_cog_VIS12"]] = cogVIS12

#cogVIS24 = cogVIS24 %>% rename_at(vars(-SUBJID), ~ paste0("SA_",.))
data_all[["stalone_cog_VIS24"]] = cogVIS24

#cogVIS36 = cogVIS36 %>% rename_at(vars(-SUBJID), ~ paste0("SA_",.))
data_all[["stalone_cog_VIS36"]] = cogVIS36


save.image("~/Documents/Documents_IT/paper/ADNI_VAMBN_paper_final/cogDataPred.RData")
setwd("~/Documents/Documents_IT/paper/ADNI_VAMBN_paper_final")
data_out<-'data/data_out/'
saveRDS(data_all, file = paste0(data_out,'data_all_imp_cog.rds'))

