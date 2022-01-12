############################## Dependencies and helper functions (dependencies for helper functions are loaded here)
setwd("~/ADNI_VAMBN_paper")
rm(list=ls())
load("~/Altoida_VAMBN_paper/createbnAltoida.RData")

##keep the required target nodes that needs to be simulated in ADNI data#
rm(list=setdiff(ls(),  "discdata"))
discdataAltoida = discdata
library(tidyverse)
library(beepr)
library(arules)
library(mclust)
library(rpart)
library(bnlearn) # hc might be overwritten by arules or some such package "bnlearn::hc" if so; not currently used though
library(parallel)
# general helpers
source('helper/plot_bn.R')
source('helper/clean_help.R')
source('helper/simulate_VP.R')
source('helper/VP_missslist.R')
# study specific helpers
source('helper/merge_data.R')
source('helper/addnoise.R')
source('helper/make_bl_wl_adni_altoida_common.R')
source('helper/save_VPmisslist.R')
############################## Settings and preprocessing

########## Name output files
name<-'main'
data_out<-paste0('data/data_out/',name)

########## Load data & remaining formatting of standalone
data<-merge_data() 
colnames(data) = gsub("SA_zcode_", "zcode_", colnames(data))
##import the features from altoida data##
names(data)[names(data) == "SA_AGE_VIS1"] = "SA_age_VIS1"
names(data)[names(data) == "SA_PTGENDER_VIS1"] = "SA_gender_VIS1"
names(data)[names(data) == "SA_PTEDUCAT_VIS1"] = "SA_yearsOfEducation_VIS1"
commonFeatures = c("SA_age_VIS1", "SA_gender_VIS1", "SA_yearsOfEducation_VIS1", "SA_Amyloid_VIS1")
commonFeatures = c(commonFeatures, grep("MMSE|zcode\\_AR|BIT|DOT|Motor|DX|Com|Per|Inh|Fle|Vis|Pla|Pro|Spa", colnames(data), value = TRUE))
data = data[,commonFeatures]

# refactor all factor columns (so there are no empty levels)
data$SA_gender_VIS1 = as.factor(data$SA_gender_VIS1)
for(col in colnames(data)){
  if (is.factor(data[,col])|grepl('scode',col))
    data[,col]<-factor(data[,col])
}

# remove subject variable
pt<-data$SUBJID
data$SUBJID<-NULL

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
scr<-"bic-cg" # 'bic' for basic autoencoder and fully discretized
mth<-"mle" # 'bayes' for basic autoencoder and fully discretized

######### Add AUX superior
discdata$visitmiss_VIS6<-factor(ifelse(apply(discdata[,grepl('AUX_',colnames(discdata))&grepl('_VIS6',colnames(discdata))],1,function(x) (all(x==1))),1,0))
discdata$visitmiss_VIS12<-factor(ifelse(apply(discdata[,grepl('AUX_',colnames(discdata))&grepl('_VIS12',colnames(discdata))],1,function(x) (all(x==1))),1,0))
discdata$visitmiss_VIS18<-factor(ifelse(apply(discdata[,grepl('AUX_',colnames(discdata))&grepl('_VIS18',colnames(discdata))],1,function(x) (all(x==1))),1,0))
discdata$visitmiss_VIS24<-factor(ifelse(apply(discdata[,grepl('AUX_',colnames(discdata))&grepl('_VIS24',colnames(discdata))],1,function(x) (all(x==1))),1,0))
discdata$visitmiss_VIS36<-factor(ifelse(apply(discdata[,grepl('AUX_',colnames(discdata))&grepl('_VIS36',colnames(discdata))],1,function(x) (all(x==1))),1,0))

d6<-as.data.frame(t(discdata[,(grepl('AUX_',colnames(discdata))&grepl('_VIS6',colnames(discdata)))|grepl('visitmiss_VIS6',colnames(discdata))]))
d12<-as.data.frame(t(discdata[,(grepl('AUX_',colnames(discdata))&grepl('_VIS12',colnames(discdata)))|grepl('visitmiss_VIS12',colnames(discdata))]))
d18<-as.data.frame(t(discdata[,(grepl('AUX_',colnames(discdata))&grepl('_VIS18',colnames(discdata)))|grepl('visitmiss_VIS18',colnames(discdata))]))
d24<-as.data.frame(t(discdata[,(grepl('AUX_',colnames(discdata))&grepl('_VIS24',colnames(discdata)))|grepl('visitmiss_VIS24',colnames(discdata))]))
d36<-as.data.frame(t(discdata[,(grepl('AUX_',colnames(discdata))&grepl('_VIS36',colnames(discdata)))|grepl('visitmiss_VIS36',colnames(discdata))]))

func.rm = function(d){
  rm =c()
  for(i in 1:nrow(d)){
    x = duplicated(d[i,nrow(d),])
    print(x)
    if(x == TRUE){
      print(rownames(d[i,]))
      rm = c(rm,rownames(d[i,]))
    }
    #print(rownames(d1[i,]))
  }
}

rm6 = func.rm(d6)
rm12 = func.rm(d12)
rm18 = func.rm(d18)
rm24 = func.rm(d24)
rm36 = func.rm(d36)


rm<-c(rm6, rm12,rm24,rm36) # orphaned nodes are those whose immediate parent AUX is identical to the visitmiss at that visit
lowaux<-discdata[,grepl('AUX_',colnames(discdata))&!(colnames(discdata) %in% rm)]
lowaux<-colnames(lowaux)[sapply(colnames(lowaux),function(x) sum(as.numeric(as.character(lowaux[,x])))<=5)]
discdata<-discdata[ , !(names(discdata) %in% rm)]
discdata<-discdata[ , !(names(discdata) %in% lowaux)]
orphans<-gsub('AUX_','',rm)
orphans<-unname(sapply(orphans,function(x) ifelse(!grepl('SA_',x),paste0('zcode_',x),x)))

########## Make bl/wl
blname<-paste0(data_out,'_common_bl.csv')
wlname<-paste0(data_out,'_common_wl.csv')


make_bl_wl_adni_altoida_common(discdata,blname,wlname,F,orphans) # rm has info about orphaned nodes
bl<-read.csv(blname)
bl$from = as.character(bl$from)
bl$to = as.character(bl$to)
wl<-read.csv(wlname)
wlCogDomain = read.csv("wlCogDomain.csv")
wlCogDomain$X = NULL

wlCogDomain = wlCogDomain[wlCogDomain$from %in% colnames(discdata),]

wlCogDomain_vis6 = data.frame(lapply(wlCogDomain, function(x) {
  gsub("_VIS1", "_VIS6", x)
}))
wlCogDomain_vis12 = data.frame(lapply(wlCogDomain, function(x) {
  gsub("_VIS1", "_VIS12", x)
}))
wlCogDomain_vis24 = data.frame(lapply(wlCogDomain, function(x) {
  gsub("_VIS1", "_VIS24", x)
}))
wlCogDomain_vis36 = data.frame(lapply(wlCogDomain, function(x) {
  gsub("_VIS1", "_VIS36", x)
}))
cogDomains = grep("PerceptualMotorCoordination|Planning|ProspectiveMemory|SpatialMemory|CognitiveProcessingSpeed|ComplexAttention|EyeMovement|Flexibility|Inhibition|VisualPerception", colnames(discdata), value = TRUE)
##remove Speech
wlCogDomainAll = rbind.data.frame(wlCogDomain,wlCogDomain_vis6, wlCogDomain_vis12, wlCogDomain_vis24,
                                  wlCogDomain_vis36)
wlCogDomainAll = wlCogDomainAll[wlCogDomainAll$to %in% cogDomains,]
wl = rbind.data.frame(wl, wlCogDomainAll)  
wl$from = as.character(wl$from)
wl$to = as.character(wl$to)
wl$from = as.character(wl$from)
wl$to = as.character(wl$to)
bl = unique(bl)
wl = unique(wl)
bl = anti_join(bl,wl)
############################## Bnet

######### Final bayesian network
set.seed(123)
finalBNCommon =hc(discdata,maxp = 5,blacklist = bl, whitelist = wl, score="bic-cg")
save.image("finalBN_common.RData")
##switch to hc if tabu doesn't work,
##try anther scoring algorithms

##extract the common connections between the graph##

######### Bootstrapped network
cores = detectCores()
set.seed(234)
boot.stren = boot.strength(discdata, algorithm="hc", R=1000, algorithm.args = list(maxp=5, blacklist=bl, whitelist = wl, score=scr))
#boot.stren = read.csv("boot_stren_maxp5_common.csv")
#stopCluster(cl)
##results from cluster
boot.strenwithThreshold = boot.stren[boot.stren$strength >=0.50&boot.stren$direction>=0.5, ]
boot.strenwithThreshold$from = gsub('zcode_', '', boot.strenwithThreshold$from)
boot.strenwithThreshold$to = gsub('zcode_', '', boot.strenwithThreshold$to)
boot.strenwithThreshold$from = gsub('SA_', '', boot.strenwithThreshold$from)
boot.strenwithThreshold$to = gsub('SA_', '', boot.strenwithThreshold$to)

boot.strenwithThreshold = boot.strenwithThreshold[- grep("scode", boot.strenwithThreshold$from),]
boot.strenwithThreshold = boot.strenwithThreshold[- grep("AUX", boot.strenwithThreshold$from),]
boot.strenwithThreshold = boot.strenwithThreshold[- grep("AUX", boot.strenwithThreshold$to),]
boot.strenwithThreshold = boot.strenwithThreshold[- grep("visitmiss", boot.strenwithThreshold$from),]
boot.strenwithThreshold$strength = round(boot.strenwithThreshold$strength,2)
#boot.strenwithThreshold = boot.strenwithThreshold[- grep("visitmiss", boot.strenwithThreshold$to),]
write.csv(boot.strenwithThreshold, "boot.strenwithThresholdCommon_paper.csv")
#bootstrenPaper = read.csv("boot.strenwithThreshold_paper.csv")
boot.strenwithThresholdCP = boot.strenwithThreshold
boot.strenwithThresholdDirectToTasks = boot.strenwithThresholdCP[grepl("MMSE|PT|age|DX|csf|APOE|ABETA|gender|Amyloid", boot.strenwithThresholdCP$from) & grepl("AR|Motor|BIT", boot.strenwithThresholdCP$to),]
cogDomains = gsub("SA_", "", cogDomains)
boot.strenwithThresholdDirectToDigCogDomains = boot.strenwithThresholdCP[grepl("MMSE|PT|age|DX|csf|APOE|ABETA|gender|Amyloid", boot.strenwithThresholdCP$from) & (boot.strenwithThresholdCP$to%in% cogDomains),]
boot.strenwithThresholdDirectFromDigtoCog = boot.strenwithThresholdCP[grepl("AR|Motor|BIT", boot.strenwithThresholdCP$from) & (boot.strenwithThresholdCP$to%in% cogDomains),]

boot.strenwithThresholdDirect = rbind.data.frame(boot.strenwithThresholdDirectToTasks, boot.strenwithThresholdDirectToDigCogDomains)
##visual bootstrapped network

write.csv(boot.strenwithThresholdDirect, "boot.strenwithThresholdDirectCommon_visual.csv")


saveRDS(boot.stren,paste0(data_out,'_bootBN_common.rds'))
