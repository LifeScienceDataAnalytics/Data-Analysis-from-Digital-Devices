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
library(tidyverse)
library(beepr)
library(arules)
library(mclust)
library(rpart)
library(bnlearn) # hc might be overwritten by arules or some such package "bnlearn::hc" if so; not currently used though
library(parallel)
library(igraph)
# general helpers
source('helper/plot_bn.R')
source('helper/clean_help.R')
source('helper/simulate_VP.R')
source('helper/VP_missslist.R')
# study specific helpers
source('helper/merge_data.R')
source('helper/addnoise.R')
source('helper/make_bl_wl_adni.R')
source('helper/save_VPmisslist.R')
############################## Settings and preprocessing

########## Name output files
name<-'main'
data_out<-paste0('data/data_out/',name)

########## Load data & remaining formatting of standalone
data<-merge_data() 

##add VIS1 as suffix to the mechanisms
pathways = c(grep("[subgraph|Receptors|transduction]$", colnames(data), value = TRUE))
data = data %>% rename_with(~paste0(., "_VIS1"), pathways)

# refactor all factor columns (so there are no empty levels)
data$SA_PTGENDER_VIS1 = as.factor(data$SA_PTGENDER_VIS1)
data$SA_PTMARRY_VIS1 = as.factor(data$SA_PTMARRY_VIS1)
data$SA_PTEDUCAT_VIS1 = as.factor(data$SA_PTEDUCAT_VIS1)
#data$SA_PTEDUCAT_VIS1 = as.factor(data$SA_PTEDUCAT_VIS1)
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

# remove subject variable
#pt<-data$SUBJID
#data$SUBJID<-NULL
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

rm<-c(rm6,rm12,rm24,rm36) # orphaned nodes are those whose immediate parent AUX is identical to the visitmiss at that visit
lowaux<-discdata[,grepl('AUX_',colnames(discdata))&!(colnames(discdata) %in% rm)]
lowaux<-colnames(lowaux)[sapply(colnames(lowaux),function(x) sum(as.numeric(as.character(lowaux[,x])))<=5)]
discdata<-discdata[ , !(names(discdata) %in% rm)]
discdata<-discdata[ , !(names(discdata) %in% lowaux)]
orphans<-gsub('AUX_','',rm)
orphans<-unname(sapply(orphans,function(x) ifelse(!grepl('SA_',x),paste0('zcode_',x),x)))

########## Make bl/wl
blname<-paste0(data_out,'_bl.csv')
wlname<-paste0(data_out,'_wl.csv')
make_bl_wl_adni(discdata,blname,wlname,F,orphans) # rm has info about orphaned nodes
bl<-read.csv(blname)
bl$from = as.character(bl$from)
bl$to = as.character(bl$to)
wl<-read.csv(wlname)
wlCogDomain = read.csv("wlCogDomain.csv")
wlCogDomain$X = NULL
wlCogDomain$from = paste("SA_", wlCogDomain$from, sep = "")
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
bl = unique(bl)
wl = unique(wl)
bl = anti_join(bl,wl)
bl = unique(bl)
wl = unique(wl)
############################## Bnet

######### Final bayesian network


#ds = discdata[!grepl("AUX|visitmiss", colnames(discdata))]
#blUpd = bl[!grepl("AUX|visitmiss", bl$from),]
#blUpd = blUpd[!grepl("AUX|visitmiss", blUpd$to),]
#cvres1 = bn.cv(discdata, "hc", runs=10, algorithm.args = list(maxp=5, blacklist=bl,  whitelist=wl, score="bic-cg")) 
#cvres2 = bn.cv(discdata, "tabu", runs=10, algorithm.args = list(maxp=5, blacklist=bl,  whitelist=wl, score="bic-cg")) 
discdata[] <- lapply(discdata, function(x) if(is.integer(x)) as.numeric(x) else x)
saveRDS(discdata,paste0(data_out,'_discdata.rds'))

set.seed(123)
#cvres3Gn = bn.cv(discdata, "hc", runs=10, loss="logl-cg", algorithm.args = list(maxp=5,blacklist=bl,  whitelist=wl, score="bic-cg"), debug = TRUE) 
finalBNFive =hc(discdata,maxp = 5, blacklist = bl, whitelist = wl, score="bic-cg")
set.seed(121)
#cvres3Gn = bn.cv(discdata, "hc", runs=10, loss="logl-cg", algorithm.args = list(maxp=5,blacklist=bl,  whitelist=wl, score="bic-cg"), debug = TRUE) 
finalBNFifteen =hc(discdata,maxp = 15, blacklist = bl, whitelist = wl, score="bic-cg")
set.seed(123)
finalBNAll =hc(discdata, blacklist = bl, whitelist = wl, score="bic-cg")
save.image("finalBN.RData")



##switch to hc if tabu doesn't work,
##try anther scoring algorithms
saveRDS(finalBNAll,paste0(data_out,'_finalBN.rds'))

######### Bootstrapped network
cores = detectCores()
#cl =  makeCluster(cores)
set.seed(234)
boot.stren = boot.strength(discdata, algorithm="hc", R=100, algorithm.args = list(blacklist=bl, whitelist = wl, score=scr))
#stopCluster(cl)

##results from cluster IN ADNIVAMBN_BootStrenpaper
booot.stren.1000 = read.csv("boot_stren_all.csv")
boot.strenwithThreshold = booot.stren.1000[booot.stren.1000$strength >=0.1&booot.stren.1000$direction>=0.5, ]

##extract nodes connected from mechanism to digital nodes
boot.stren.subgraphToDigTasks = booot.stren.1000[grepl("zcode\\_.*subgraph", booot.stren.1000$from) & grepl("AR|Motor|BIT", booot.stren.1000$to)&booot.stren.1000$strength >=0.3&booot.stren.1000$direction>=0.5,]
boot.stren.subgraphToCog = booot.stren.1000[grepl("zcode\\_.*subgraph", booot.stren.1000$from) & booot.stren.1000$to%in%cogDomains & booot.stren.1000$strength>0.3&booot.stren.1000$direction>=0.5,]
boot.stren.faqToDigTasks = booot.stren.1000[grepl("FAQ", booot.stren.1000$from) & grepl("AR|Motor|BIT", booot.stren.1000$to) & booot.stren.1000$strength>0.3&booot.stren.1000$direction>=0.5,]
boot.stren.faqToCog = booot.stren.1000[grepl("FAQ", booot.stren.1000$from) & booot.stren.1000$to%in%cogDomains & booot.stren.1000$strength>0.3&booot.stren.1000$direction>=0.5,]
boot.stren.volumecsfimagingPHStoDigTasks = booot.stren.1000[grepl("volume|csf|imaging|PHS", booot.stren.1000$from) & grepl("AR|Motor|BIT", booot.stren.1000$to) & booot.stren.1000$strength>0.3&booot.stren.1000$direction>=0.5,]
boot.stren.volumecsfimagingPHSToCog = booot.stren.1000[grepl("volume|csf|imaging|PHS", booot.stren.1000$from) & booot.stren.1000$to%in%cogDomains & booot.stren.1000$strength>0.3&booot.stren.1000$direction>=0.5,]


#boot.strenwithThresholdCP = boot.strenwithThreshold

boot.strenwithThreshold$from = gsub('zcode_', '', boot.strenwithThreshold$from)
boot.strenwithThreshold$to = gsub('zcode_', '', boot.strenwithThreshold$to)
boot.strenwithThreshold$from = gsub('SA_', '', boot.strenwithThreshold$from)
boot.strenwithThreshold$to = gsub('SA_', '', boot.strenwithThreshold$to)

boot.strenwithThreshold = boot.strenwithThreshold[- grep("scode", boot.strenwithThreshold$from),]
boot.strenwithThreshold = boot.strenwithThreshold[- grep("AUX", boot.strenwithThreshold$from),]
boot.strenwithThreshold = boot.strenwithThreshold[- grep("AUX", boot.strenwithThreshold$to),]
boot.strenwithThreshold = boot.strenwithThreshold[- grep("visitmiss", boot.strenwithThreshold$from),]


write.csv(boot.strenwithThreshold, "boot.strenwithThreshold_all.csv")

##visual bootstrapped network
bootStrenVisual = rbind.data.frame(boot.stren.subgraphToDigTasks,boot.stren.subgraphToCog, boot.stren.faqToDigTasks, 
                                   boot.stren.faqToCog, boot.stren.volumecsfimagingPHStoDigTasks, boot.stren.volumecsfimagingPHSToCog,
                                   boot.strenwithThreshold)
cogDomains = gsub("SA_", "", cogDomains)
bootStrenVisual$from = gsub('zcode_', '', bootStrenVisual$from)
bootStrenVisual$to = gsub('zcode_', '', bootStrenVisual$to)
bootStrenVisual$from = gsub('SA_', '', bootStrenVisual$from)
bootStrenVisual$to = gsub('SA_', '', bootStrenVisual$to)

bootStrenVisualDig = bootStrenVisual[grepl("AR|BIT|Motor", bootStrenVisual$from)&grepl("AR|BIT|Motor", bootStrenVisual$to),]
bootStrenVisualDigtoCog = bootStrenVisual[grepl("AR|BIT|Motor", bootStrenVisual$from)&bootStrenVisual$to%in%cogDomains,]
bootStrenVisualCog = bootStrenVisual[bootStrenVisual$from%in%cogDomains&bootStrenVisual$to%in%cogDomains,]
bootStrenVisualFAQ = bootStrenVisual[grepl("FAQ", bootStrenVisual$from)&grepl("FAQ", bootStrenVisual$to),]
bootStrenVisualsubgraph = bootStrenVisual[grepl("subgraph", bootStrenVisual$from)&grepl("subgraph", bootStrenVisual$to),]
bootStrenVisualMMSE = bootStrenVisual[grepl("MMSE", bootStrenVisual$from)&grepl("MMSE", bootStrenVisual$to),]
bootStrenVisualDX = bootStrenVisual[grepl("DX", bootStrenVisual$from),]
bootStrenVisualage = bootStrenVisual[grepl("AGE", bootStrenVisual$from),]

bootStrenVisual = bootStrenVisual[!bootStrenVisual$X %in% c(bootStrenVisualDig$X,bootStrenVisualFAQ$X,bootStrenVisualsubgraph$X,
                                                            bootStrenVisualMMSE$X, bootStrenVisualDigtoCog$X, bootStrenVisualCog$X,
                                                            bootStrenVisualDX$X, bootStrenVisualage$X),]
bootStrenVisual$strength = round(bootStrenVisual$strength, 2)
bootStrenVisual$EdgeType =  ifelse(bootStrenVisual$strength >= 0.50,"highStrength", "LowStrength")
bootStrenVisual = unique(bootStrenVisual)
write.csv(bootStrenVisual, "boot.strenwithThreshold_visual.csv")
saveRDS(booot.stren.1000,paste0(data_out,'_bootBN.rds'))
#beep()

# save fitted network
#real = discdata
real = discdata
finalBNarcs = as.data.frame(finalBNAll$arcs)
finalBNarcsCp = finalBNarcs
finalBNarcsCp$from = gsub('zcode_', '', finalBNarcsCp$from)
finalBNarcsCp$to = gsub('zcode_', '', finalBNarcsCp$to)
finalBNarcsCp$from = gsub('SA_', '', finalBNarcsCp$from)
finalBNarcsCp$to = gsub('SA_', '', finalBNarcsCp$to)
finalBNarcsCp = finalBNarcsCp[- grep("AUX", finalBNarcsCp$from),]
finalBNarcsCp = finalBNarcsCp[- grep("AUX", finalBNarcsCp$to),]
finalBNarcsCp = finalBNarcsCp[- grep("visitmiss", finalBNarcsCp$from),]

write.csv(as.data.frame(finalBNarcsCp), "finalBN.csv")

fitted = bn.fit(finalBNAll, real, method=mth)
saveRDS(fitted,paste0(data_out,'_finalBN_fitted.rds'))




############################## VP vs RP
# Virtual Patient Generation
set.seed(123)
virtual<-simulate_VPs(real,finalBNAll,iterative=FALSE,scr,mth,wl,bl, n= nrow(real))

############################
############################ save out all data
############################

# save out real and virtual patients
#virtual = virtual[1:1000,]
real$SUBJID<-pt
write.csv(real,paste0(data_out,'_RealPPts.csv'),row.names=FALSE)
saveRDS(real,paste0(data_out,'_RealPPts.rds'))

# save out VP misslist (for HIVAE decoding, tells HIVAE which zcodes the BN considers missing)
data_meta1<-read.csv('data/HI-VAE/metaenc.csv')
pathways = c(grep("[subgraph|Receptors|transduction]$", colnames(data_meta1), value = TRUE))
data_meta1 = data_meta1 %>% rename_with(~paste0(., "_VIS1"), pathways)

data_meta2<-read.csv('data/HI-VAE/metaenc_nm.csv')

misCols = setdiff(grep("^zcode|scode", colnames(data_meta1), value = TRUE), "SUBJID")

misCols = intersect(misCols, colnames(virtual))
#digitalCols = grep("BIT|DOT", colnames(virtual), value = TRUE)
#virtualWithoutDigital = virtual[,setdiff(colnames(virtual), digitalCols)]
write.csv(virtual,paste0(data_out,'_VirtualPPts.csv'),row.names=FALSE)
saveRDS(virtual,paste0(data_out,'_VirtualPPts.rds'))
save_VPmisslist(virtual[,misCols],'data/HI-VAE/')


cols_nm = intersect(colnames(data_meta2), colnames(virtual))
write.csv(virtual[,cols_nm],paste0(data_out,'_VirtualPPts_nm.csv'),row.names=FALSE)
saveRDS(virtual[,cols_nm],paste0(data_out,'_VirtualPPts_nm.rds'))


##add suffix to subgraph columns VIS1
metaenc <- read_csv("data/HI-VAE/metaenc.csv")
pathways = c(grep("[subgraph|Receptors|transduction]$", colnames(metaenc), value = TRUE))
metaencUpd = metaenc %>% rename_with(~paste0(., "_VIS1"), pathways)
write.csv(metaenc,"data/HI-VAE/metaencUpd.csv")

##change the subgraph names in hyperparameters file
best_hyper_ADNI_processed_upd <- read_csv("data/HI-VAE/best_hyper_ADNI_processed_upd.csv")
best_hyper_ADNI_processed_upd_rename <- 
  best_hyper_ADNI_processed_upd %>% 
  mutate_at("files", str_replace, "subgraph", "subgraph_VIS1")
best_hyper_ADNI_processed_upd_rename <- 
  best_hyper_ADNI_processed_upd_rename %>% 
  mutate_at("files", str_replace, "transduction", "transduction_VIS1")

write.csv(best_hyper_ADNI_processed_upd_rename, "data/HI-VAE/best_hyper_ADNI_processed_upd_rename.csv")

###### Virtual Patient Validation
#roc<-validate_VP(real=real,virtual=virtual,proc=F) # full AUC rather than partial AUC
#proc<-validate_VP(real=real,virtual=virtual,proc=T) # partial AUC with focus on sensitivity
#beep()
# if validation of ARAE/VAE, saveout is in jupyter notebooks (VAE_decoded and mclust_decode(R notebook! could make it a script if easier)/ARAE_decoded)













