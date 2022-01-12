############# README
# Correlation structure
############# 

rm(list=ls())
library(tidyverse)
library(beepr)
library(arules)
library(mclust)
library(rpart)
library(bnlearn)
library(parallel)
library(randomForestSRC)
library(coin)
library(Matching)
library(plyr)
library(corrplot)
library(Hmisc)
source('helper/plot_bn.R')
source('helper/clean_help.R')
source('helper/merge_data.R')
source('helper/discr.R')
source('helper/addnoise.R')
source('helper/make_bl_wl_adni.R')

source('helper/simulate_VP.R')
source('helper/validate_VP.R')
source('helper/make_dummy.R')
source('helper/compare_dists.R')

############ data prep
name<-'main'
data_out<-paste0('data/data_out/',name)
#data_all<-readRDS('data/data_condensed.rds')
data_all<-readRDS('data/data_out/data_all_imp_dig_cog.rds')
virtual_m<-readRDS(paste0(data_out,'_VirtualPPts.rds'))
virtual_nm <-  readRDS(paste0(data_out,'_VirtualPPts_nm.rds'))
virtual = merge(virtual_m, virtual_nm)
#colnames(virtual)<-gsub('_VIS1_','_',colnames(virtual))
real<-readRDS(paste0(data_out,'_RealPPts.rds'))
#colnames(real)<-gsub('_VIS1_','_',colnames(real))
orig<-data_all[!grepl('SA',names(data_all))]
orig<-orig %>% reduce(merge, by = 'SUBJID')
#orig <- reduce(orig, left_join)
colnames(orig)<-gsub(' ','\\.',colnames(orig))
colnames(orig)<-gsub('[()]','\\.',colnames(orig))
names(orig)[names(orig)=="Insulin_signal_transduction_VIS1"] = "Insulin_signal_transduction_rs1999763_VIS1"
names(orig)[names(orig)=="Nerve_growth_factor_subgraph_VIS1"] = "Nerve_growth_factor_subgraph_rs3775256_VIS1"

decRP<-read_csv('data/HI-VAE/reconRP.csv')
decRP_nm<-read_csv('data/HI-VAE/reconRP_nm.csv')
decRP = merge(decRP, decRP_nm)
names(decRP)[names(decRP)=="Insulin_signal_transduction"] = "Insulin_signal_transduction_rs1999763_VIS1"
names(decRP)[names(decRP)=="Nerve_growth_factor_subgraph"] = "Nerve_growth_factor_subgraph_rs3775256_VIS1"

decVP<-read.csv('data/HI-VAE/decodedVP.csv')
decVP_nm<-read.csv('data/HI-VAE/decodedVP_nm.csv')
decVP = cbind.data.frame(decVP, decVP_nm)
names(decVP)[names(decVP)=="Insulin_signal_transduction_VIS1"] = "Insulin_signal_transduction_rs1999763_VIS1"
names(decVP)[names(decVP)=="Nerve_growth_factor_subgraph_VIS1"] = "Nerve_growth_factor_subgraph_rs3775256_VIS1"



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
orig<-orig[,colnames(decRP)]
orig<-orig[,sort(colnames(orig))]
decRP<-decRP[,sort(colnames(decRP))]
decVP<-decVP[,sort(colnames(decVP))]
orig<-orig[,sapply(colnames(orig),function(x) !is.factor(orig[,x]))]
decRP<-decRP[,sapply(colnames(decRP),function(x) !is.factor(decRP[,x]))]
decVP<-decVP[,sapply(colnames(decVP),function(x) !is.factor(decVP[,x]))]

if (!((all((colnames(decRP))==(colnames(orig))))&(all((colnames(decRP))==(colnames(decVP))))&(all((colnames(orig))==(colnames(decVP))))))
  warning('columns not matching!')

oc<-rcorr(as.matrix(orig))
rc<-rcorr(as.matrix(decRP))
vc<-rcorr(as.matrix(decVP))

png("output_figures/Adni_corr_orig.png") 
corrplot(oc$r, type = "upper",method='color',main=paste0('Original'),tl.cex=0.25, tl.pos='n',mar=c(0,0,1.5,0))
dev.off() 

png("output_figures/Adni_corr_decRP.png") 
corrplot(rc$r, type = "upper",method='color',main=paste0('Decoded real'),tl.cex=0.25, tl.pos='n',mar=c(0,0,1.5,0))
dev.off() 

png("output_figures/Adni_corr_decVP.png") 
corrplot(vc$r, type = "upper",method='color',main=paste0('Decoded virtual'),tl.cex=0.25, tl.pos='n',mar=c(0,0,1.5,0))
dev.off() 

pdf("output_figures/Adni_corr_orig.pdf") 
corrplot(oc$r, type = "upper",method='color',main=paste0('Original (Norm: ',round(norm(oc$r,type='F'),2),')'),tl.cex=0.25, tl.pos='n',mar=c(0,0,1.5,0))
dev.off() 

pdf("output_figures/Adni_corr_decRP.pdf") 
corrplot(rc$r, type = "upper",method='color',main=paste0('Decoded real (Norm: ',round(norm(rc$r,type='F'),2),')'),tl.cex=0.25, tl.pos='n',mar=c(0,0,1.5,0))
dev.off() 

pdf("output_figures/Adni_corr_decVP.pdf") 
corrplot(vc$r, type = "upper",method='color',main=paste0('Decoded virtual (Norm: ',round(norm(vc$r,type='F'),2),')'),tl.cex=0.25, tl.pos='n',mar=c(0,0,1.5,0))
dev.off() 


