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
setwd("~/Documents/Documents_IT/paper/ADNI_VAMBN_paper_final")
rm(list=ls())
load("~/Documents/Documents_IT/paper/Altoida_VAMBN_paper/createbnAltoida.RData")
#load("~/Documents/Documents_IT/VAMBN/randomForestPred_Nov26.RData")
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

#data <- data %>% select(-contains("VIS48"))
#colnames(data) = paste(colnames(data),"_VIS1", sep="")
# could load just basic autoencoded with T (if so need to fully discretize)
#remCorticalBl = grep("SA_cortical", colnames(data), value = TRUE)
#data<-data[,!grepl('SA_cortical',colnames(data))]
#data<-data[,!grepl('DX|APOE',colnames(data))]

# refactor all factor columns (so there are no empty levels)
data$SA_gender_VIS1 = as.factor(data$SA_gender_VIS1)
#data$SA_PTMARRY_VIS1 = as.factor(data$SA_PTMARRY_VIS1)
#data[,grep("^SA.*MMSE", colnames(data))] <- lapply(data[,grep("^SA.*MMSE", colnames(data))], as.numeric)
for(col in colnames(data)){
  if (is.factor(data[,col])|grepl('scode',col))
    data[,col]<-factor(data[,col])
}

# remove subject variable
pt<-data$SUBJID
data$SUBJID<-NULL

######### Discretize & set score
discdata<-data
#dat = discdata

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
#discdata$visitmiss_VIS48<-factor(ifelse(apply(discdata[,grepl('AUX_',colnames(discdata))&grepl('_VIS48',colnames(discdata))],1,function(x) (all(x==1))),1,0))
#discdata$visitmiss_VIS60<-factor(ifelse(apply(discdata[,grepl('AUX_',colnames(discdata))&grepl('_VIS60',colnames(discdata))],1,function(x) (all(x==1))),1,0))
#discdata$visitmiss_VIS72<-factor(ifelse(apply(discdata[,grepl('AUX_',colnames(discdata))&grepl('_VIS72',colnames(discdata))],1,function(x) (all(x==1))),1,0))


d6<-as.data.frame(t(discdata[,(grepl('AUX_',colnames(discdata))&grepl('_VIS6',colnames(discdata)))|grepl('visitmiss_VIS6',colnames(discdata))]))
d12<-as.data.frame(t(discdata[,(grepl('AUX_',colnames(discdata))&grepl('_VIS12',colnames(discdata)))|grepl('visitmiss_VIS12',colnames(discdata))]))
d18<-as.data.frame(t(discdata[,(grepl('AUX_',colnames(discdata))&grepl('_VIS18',colnames(discdata)))|grepl('visitmiss_VIS18',colnames(discdata))]))
d24<-as.data.frame(t(discdata[,(grepl('AUX_',colnames(discdata))&grepl('_VIS24',colnames(discdata)))|grepl('visitmiss_VIS24',colnames(discdata))]))
d36<-as.data.frame(t(discdata[,(grepl('AUX_',colnames(discdata))&grepl('_VIS36',colnames(discdata)))|grepl('visitmiss_VIS36',colnames(discdata))]))
#d48<-as.data.frame(t(discdata[,(grepl('AUX_',colnames(discdata))&grepl('_VIS48',colnames(discdata)))|grepl('visitmiss_VIS48',colnames(discdata))]))
#d60<-as.data.frame(t(discdata[,(grepl('AUX_',colnames(discdata))&grepl('_VIS60',colnames(discdata)))|grepl('visitmiss_VIS60',colnames(discdata))]))
#d72<-as.data.frame(t(discdata[,(grepl('AUX_',colnames(discdata))&grepl('_VIS72',colnames(discdata)))|grepl('visitmiss_VIS72',colnames(discdata))]))

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
# rm6<-rownames(d6)[duplicated(d6,fromLast=TRUE)]
# rm12<-rownames(d12)[duplicated(d12,fromLast=TRUE)]
# rm18<-rownames(d18)[duplicated(d18,fromLast=TRUE)]
# rm24<-rownames(d24)[duplicated(d24,fromLast=TRUE)]
# rm36<-rownames(d36)[duplicated(d36,fromLast=TRUE)]
#rm48<-rownames(d36)[duplicated(d48,fromLast=TRUE)]
#rm60<-rownames(d36)[duplicated(d60,fromLast=TRUE)]
#rm72<-rownames(d36)[duplicated(d72,fromLast=TRUE)]

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
#wlCogDomain$from = paste("SA_", wlCogDomain$from, sep = "")
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
##remove Sppech
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


#ds = discdata[!grepl("AUX|visitmiss", colnames(discdata))]
#blUpd = bl[!grepl("AUX|visitmiss", bl$from),]
#blUpd = blUpd[!grepl("AUX|visitmiss", blUpd$to),]
#cvres1 = bn.cv(discdata, "hc", runs=10, algorithm.args = list(maxp=5, blacklist=bl,  whitelist=wl, score="bic-cg")) 
#cvres2 = bn.cv(discdata, "tabu", runs=10, algorithm.args = list(maxp=5, blacklist=bl,  whitelist=wl, score="bic-cg")) 
set.seed(123)
finalBNCommon =hc(discdata,maxp = 5,blacklist = bl, whitelist = wl, score="bic-cg")
save.image("finalBN_common.RData")
##switch to hc if tabu doesn't work,
##try anther scoring algorithms

##extract the common connections between the graph##

######### Bootstrapped network
cores = detectCores()
#cl =  makeCluster(cores)
set.seed(234)
boot.stren = boot.strength(discdata, algorithm="hc", R=1000, algorithm.args = list(maxp=5, blacklist=bl, whitelist = wl, score=scr))
boot.stren = read.csv("boot_stren_maxp5_common.csv")
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


#bootStrenVisualZoom = boot.strenwithThreshold[grep("FAQ|MMSE|subgraph|rs", boot.strenwithThreshold$from),]
#write.csv(bootStrenVisualZoom, "bootStrenVisualZoom.csv")
saveRDS(boot.stren,paste0(data_out,'_bootBN_common.rds'))
#beep()

# save fitted network
# #real = discdata
# real = discdata
# finalBNarcs = as.data.frame(finalBN$arcs)
# finalBNarcsCp = finalBNarcs
# finalBNarcsCp$from = gsub('zcode_', '', finalBNarcsCp$from)
# finalBNarcsCp$to = gsub('zcode_', '', finalBNarcsCp$to)
# finalBNarcsCp$from = gsub('SA_', '', finalBNarcsCp$from)
# finalBNarcsCp$to = gsub('SA_', '', finalBNarcsCp$to)
# finalBNarcsCp = finalBNarcsCp[- grep("AUX", finalBNarcsCp$from),]
# finalBNarcsCp = finalBNarcsCp[- grep("AUX", finalBNarcsCp$to),]
# finalBNarcsCp = finalBNarcsCp[- grep("visitmiss", finalBNarcsCp$from),]
# finalBNarcsCp = finalBNarcsCp[- grep("visitmiss", finalBNarcsCp$to),]
# #finalBNarcsCp = finalBNarcsCp[- grep("scode", finalBNarcsCp$to),]
# #finalBNarcsCpMMSEFAQ = with(finalBNarcsCp, finalBNarcsCp[ grepl( 'FAQ|PT|cog|volume|brain68|cortical|snp', finalBNarcsCp$from) & grepl( 'MMSE|PT|cog|volume|brain68|cortical|snp', finalBNarcsCp$to),])
# #finalBNarcsCpMMSEFAQUpd = with(finalBNarcsCp, finalBNarcsCp[ grepl( 'MMSE|PT|cog|volume|brain68|cortical|snp', finalBNarcsCp$from) & grepl( 'FAQ|PT|cog|volume|brain68|cortical|snp', finalBNarcsCp$to),])
# #finalBNarcsCpMMSEFAQ = rbind.data.frame(finalBNarcsCpMMSEFAQ, finalBNarcsCpMMSEFAQUpd)
# #finalBNarcsCpMMSEFAQ = unique(finalBNarcsCpMMSEFAQ)
# 
# write.csv(as.data.frame(finalBNarcsCp), "finalBN.csv")
# 
# #boot_stren_adni_100 <- read_csv("boot.stren.adni.100.csv")
# #final_BayesianNetwork = write.csv(as.data.frame(finalBNarcsCp), "/Users/msood/Documents/Documents – IT-Admin’s MacBook Pro/ADNIVAMBN/VAMBNForADNI/finalBN.csv")
# #finalBN = readRDS(paste0(data_out,'_finalBN.rds'))
# fitted = bn.fit(finalBN, real, method=mth)
# saveRDS(fitted,paste0(data_out,'_finalBN_fitted.rds'))
# 
# 
# 
# #real$SUBJID<-NULL
# #finalBN<-readRDS(paste0(data_out,'_finalBN.rds'))
# 
# ##visualization overlap between boot strength and final bayesian network
# 
# finalBNarcsVisual <- finalBNarcs
# finalBNarcsVisual$strength <- ""
# finalBNarcsVisual$direction <- ""
# 
# ##finding the overlap between two dataframes##
# boot.strenwithThresholdVisual = boot.strenwithThreshold
# boot.strenwithThresholdVisual$Overlap <- ifelse(is.na(match(paste0(boot.strenwithThresholdVisual$from, boot.strenwithThresholdVisual$to), paste0(finalBNarcsVisual$from, finalBNarcsVisual$to))),"OnlyInBootStrapped", "Overlap")
# #boot.strenwithThresholdVisual = boot.strenwithThresholdVisual[boot.strenwithThresholdVisual$Overlap=="Overlap",]
# boot.strenwithThresholdVisual$from = gsub('zcode_', '', boot.strenwithThresholdVisual$from)
# boot.strenwithThresholdVisual$to = gsub('zcode_', '', boot.strenwithThresholdVisual$to)
# boot.strenwithThresholdVisual$from = gsub('SA_', '', boot.strenwithThresholdVisual$from)
# boot.strenwithThresholdVisual$to = gsub('SA_', '', boot.strenwithThresholdVisual$to)
# #boot.strenwithThresholdVisual = boot.strenwithThresholdVisual[- grep("AUX|visitmiss", boot.strenwithThresholdVisual$from),]
# #boot.strenwithThresholdVisual = boot.strenwithThresholdVisual[- grep("AUX|visitmiss", boot.strenwithThresholdVisual$to),]
# boot.strenwithThreshold = boot.strenwithThreshold[grepl("MMSE|FAQ|subgraph|csf|volume",boot.strenwithThreshold$from)&
#                                                     grepl("BIT|DOT", boot.strenwithThreshold$to),]
# 
# 
# write.csv(boot.strenwithThresholdVisual,"boot.strenwithThresholdVisual1000.csv")
# 
# finalBNarcsVisual$Overlap <- ifelse(is.na(match(paste0(finalBNarcsVisual$from, finalBNarcsVisual$to), paste0(boot.strenwithThresholdVisual$from, boot.strenwithThresholdVisual$to))),"OnlyInfinalBN", "Overlap")
# finalBNarcsVisual = finalBNarcsVisual[- grep("AUX", finalBNarcsVisual$from),]
# finalBNarcsVisual = finalBNarcsVisual[- grep("AUX", finalBNarcsVisual$to),]
# finalBNarcsVisual = finalBNarcsVisual[- grep("visitmiss", finalBNarcsVisual$from),]
# finalBNarcsVisual = finalBNarcsVisual[- grep("visitmiss", finalBNarcsVisual$to),]
# finalBNarcsVisual = finalBNarcsVisual[finalBNarcsVisual$Overlap=="Overlap",]
# combiOverlap <- rbind.data.frame(boot.strenwithThresholdVisual, finalBNarcsVisual[which(finalBNarcsVisual$Overlap=="OnlyInfinalBN"),])
# combiOverlap$strength <- as.numeric(combiOverlap$strength )
# combiOverlapNew <- combiOverlap
# combiOverlapWithOverLabel <- combiOverlap[combiOverlap$Overlap == "Overlap",]
# combiOverlapWithOverLabel$from = gsub('zcode_', '', combiOverlapWithOverLabel$from)
# combiOverlapWithOverLabel$to = gsub('zcode_', '', combiOverlapWithOverLabel$to)
# write.csv(combiOverlapWithOverLabel,"combiOverlapAltoidaFinalBNOOverlap.csv")
# 
# 
# ############################## VP vs RP
# 
# ############################
# ############################ VP vs RP
# ############################
# 
# # Virtual Patient Generation
# set.seed(123)
# virtual<-simulate_VPs(real,finalBN,iterative=FALSE,scr,mth,wl,bl, n= nrow(real))
# 
# ############################
# ############################ save out all data
# ############################
# 
# # save out real and virtual patients
# #virtual = virtual[1:1000,]
# real$SUBJID<-pt
# #real$SUBJID<-NULL
# # real <- real %>% rename_at(vars(starts_with("zcode")), 
# #                            funs(str_replace(., "_VIS1", "")))
# # real <- real %>% rename_at(vars(starts_with("scode")), 
# #                            funs(str_replace(., "_VIS1", "")))
# write.csv(real,paste0(data_out,'_RealPPts.csv'),row.names=FALSE)
# saveRDS(real,paste0(data_out,'_RealPPts.rds'))
# # virtual <- virtual %>% rename_at(vars(starts_with("zcode")), 
# #                                  funs(str_replace(., "_VIS1", "")))
# # virtual <- virtual %>% rename_at(vars(starts_with("scode")), 
# #                                  funs(str_replace(., "_VIS1", "")))
# 
# # save out VP misslist (for HIVAE decoding, tells HIVAE which zcodes the BN considers missing)
# data_meta1<-read.csv('data/HI-VAE/metaenc.csv')
# data_meta2<-read.csv('data/HI-VAE/metaenc_nm.csv')
# names(data_meta1)[names(data_meta1) == "zcode_Bcl.2_subgraph" ] = "zcode_Bcl-2_subgraph"
# names(data_meta1)[names(data_meta1) == "zcode_Calcium.dependent_signal_transduction" ] = "zcode_Calcium-dependent_signal_transduction"
# names(data_meta1)[names(data_meta1) == "scode_Bcl.2_subgraph" ] = "scode_Bcl_2_subgraph"
# names(data_meta1)[names(data_meta1) == "scode_Calcium.dependent_signal_transduction" ] = "scode_Calcium-dependent_signal_transduction"
# names(data_meta1)[names(data_meta1) == "AUX_Bcl.2_subgraph" ] = "AUX_Bcl-2_subgraph"
# names(data_meta1)[names(data_meta1) == "AUX_Calcium.dependent_signal_transduction" ] = "AUX_Calcium-dependent_signal_transduction"
# 
# names(virtual)[names(virtual) == "scode_Bcl.2_subgraph" ] = "scode_Bcl-2_subgraph"
# names(virtual)[names(virtual) == "zcode_Bcl.2_subgraph" ] = "zcode_Bcl-2_subgraph"
# names(virtual)[names(virtual) == "scode_Calcium.dependent_signal_transduction" ] = "scode_Calcium-dependent_signal_transduction"
# names(virtual)[names(virtual) == "zcode_Calcium.dependent_signal_transduction" ] = "zcode_Calcium-dependent_signal_transduction"
# names(virtual)[names(virtual) == "AUX_Bcl.2_subgraph" ] = "AUX_Bcl-2_subgraph"
# names(virtual)[names(virtual) == "AUX_Calcium.dependent_signal_transduction" ] = "AUX_Calcium-dependent_signal_transduction"
# misCols = setdiff(grep("^zcode|scode", colnames(data_meta1), value = TRUE), "SUBJID")
# 
# misCols = intersect(misCols, colnames(virtual))
# #digitalCols = grep("BIT|DOT", colnames(virtual), value = TRUE)
# #virtualWithoutDigital = virtual[,setdiff(colnames(virtual), digitalCols)]
# write.csv(virtual,paste0(data_out,'_VirtualPPts.csv'),row.names=FALSE)
# saveRDS(virtual,paste0(data_out,'_VirtualPPts.rds'))
# save_VPmisslist(virtual[,misCols],'data/HI-VAE/')
# 
# 
# cols_nm = intersect(colnames(data_meta2), colnames(virtual))
# write.csv(virtual[,cols_nm],paste0(data_out,'_VirtualPPts_nm.csv'),row.names=FALSE)
# saveRDS(virtual[,cols_nm],paste0(data_out,'_VirtualPPts_nm.rds'))
# ####### Virtual Patient Validation
# #roc<-validate_VP(real=real,virtual=virtual,proc=F) # full AUC rather than partial AUC
# #proc<-validate_VP(real=real,virtual=virtual,proc=T) # partial AUC with focus on sensitivity
# #beep()
# # if validation of ARAE/VAE, saveout is in jupyter notebooks (VAE_decoded and mclust_decode(R notebook! could make it a script if easier)/ARAE_decoded)
# 
# 
# 
# 
# 
# 
# #save.image("~/Documents/Documents_IT/ADNIVAMBN/VAMBNForADNI/ADNIVAMBN.RData")
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
