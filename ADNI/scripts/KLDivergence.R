############# README
# This is the analysis file for computing KL divergence
# Run after bnet.R
############# 

rm(list=ls())
library(tidyverse)
library(plyr)
library(grid)
library(gridExtra)
library(Matching)
library(FNN)
source('helper/plot_bn.R')
source('helper/clean_help.R')
source('helper/compare_dists.R')

############ data prep
name<-'main'
data_out<-paste0('data/data_out/',name)
#data_all<-readRDS('data/data_condensed.rds')
data_all<-readRDS('data/data_out/data_all_imp_dig.rds')
virtual<-readRDS(paste0(data_out,'_VirtualPPts.rds'))
virtual_nm <-  readRDS(paste0(data_out,'_VirtualPPts_nm.rds'))
virtual = merge(virtual, virtual_nm)
#colnames(virtual)<-gsub('_VIS1_','_',colnames(virtual))
real<-readRDS(paste0(data_out,'_RealPPts.rds'))
#colnames(real)<-gsub('_VIS1_','_',colnames(real))
orig<-data_all[!grepl('stalone',names(data_all))]
orig<-orig %>% reduce(merge, by = 'SUBJID')
#orig <- reduce(orig, left_join)
colnames(orig)<-gsub(' ','\\.',colnames(orig))
colnames(orig)<-gsub('[()]','\\.',colnames(orig))
decRP<-read_csv('data/HI-VAE/reconRP.csv')
decRP_nm<-read_csv('data/HI-VAE/reconRP_nm.csv')
decRP = merge(decRP, decRP_nm)
decVP<-read.csv('data/HI-VAE/decodedVP.csv')
decVP_nm<-read.csv('data/HI-VAE/decodedVP_nm.csv')
decVP = cbind.data.frame(decVP, decVP_nm)
saveVirtual = cbind(virtual[,grep("SA", colnames(virtual), value= TRUE)], decVP)
saveVirtual = saveVirtual[,!grepl("AUX",colnames(saveVirtual)),]
saveVirtual$SUBJID = NULL
saveVirtual$ID <- seq.int(nrow(saveVirtual))
#write.csv(saveVirtual, "virtualADNI.csv")
orig<-orig[,(colnames(orig) %in% colnames(decRP))]
for (col in colnames(orig)){
  if (is.factor(orig[,col])){#|grepl('COGT_RAVLT.learning|COGT_RAVLT.forgetting|COGT_FAQ|COGT_MMSE|COGT_ADAS11|COGT_ADAS13',col)
    orig[,col]<-factor(orig[,col])
    lvs<-levels(orig[,col])
    names(lvs)<-as.character(0:(length(lvs)-1))
    decRP[,col]<-factor(decRP[,col],labels=unlist(lvs[levels(factor(decRP[,col]))]))
    decVP[,col]<-factor(decVP[,col],labels=unlist(lvs[levels(factor(decVP[,col]))]))
  }
}


virtual<-virtual[!grepl("AUX|visitmiss|scode", colnames(virtual))]
real<-real[!grepl("AUX|visitmiss|scode", colnames(real))]
names(real)[names(real) == "zcode_Bcl.2_subgraph" ] = "zcode_Bcl-2_subgraph"
names(real)[names(real) == "zcode_Calcium.dependent_signal_transduction" ] = "zcode_Calcium-dependent_signal_transduction"
real <- real %>% rename_at(vars(starts_with("zcode_Bcl")), 
                           funs(str_replace(., "Bcl-2_", "Bcl_")))
real <- real %>% rename_at(vars(starts_with("zcode_Calcium")), 
                           funs(str_replace(., "Calcium-", "Calcium_")))
virtual <- virtual %>% rename_at(vars(starts_with("zcode_Bcl")), 
                                 funs(str_replace(., "Bcl-2_", "Bcl_")))
virtual <- virtual %>% rename_at(vars(starts_with("zcode_Calcium")), 
                                 funs(str_replace(., "Calcium-", "Calcium_")))
real <- real %>% rename_at(vars(starts_with("SA_Calpastatin")), 
                           funs(str_replace(., "SA_Calpastatin-", "SA_Calpastatin_")))
real <- real %>% rename_at(vars(starts_with("SA_Beta")), 
                           funs(str_replace(., "SA_Beta-", "SA_Beta_")))
virtual <- virtual %>% rename_at(vars(starts_with("SA_Calpastatin")), 
                                 funs(str_replace(., "SA_Calpastatin-", "SA_Calpastatin_")))
virtual <- virtual %>% rename_at(vars(starts_with("SA_Beta")), 
                                 funs(str_replace(., "SA_Beta-", "SA_Beta_")))


## make violin plots of all other comparisons
orig<-orig[,sort(colnames(orig))]
decRP<-decRP[,sort(colnames(decRP))]
#decRP[is.na(orig)]<-NA
decVP$SUBJID<-NULL
decVP<-decVP[,sort(colnames(decVP))]
decVP <- decVP %>% rename_at(vars(starts_with("Calcium")), 
                             funs(str_replace(., "Calcium.", "Calcium-")))
decVP <- decVP %>% rename_at(vars(starts_with("Bcl")), 
                             funs(str_replace(., "Bcl.", "Bcl-")))
KL.divergence(real, virtual)

write.csv(orig, "data/HI-VAE/origForKD.csv")
write.csv(decRP, "data/HI-VAE/reconRpForKD.csv")
write.csv(decVP, "data/HI-VAE/decodedVPForKD.csv")



decRP$SUBJID<-NULL
decVP$SUBJID<-NULL
orig<-orig[,sort(colnames(orig))]
decRP<-decRP[,sort(colnames(decRP))]



