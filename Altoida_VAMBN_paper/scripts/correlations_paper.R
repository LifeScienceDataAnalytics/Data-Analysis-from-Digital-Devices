############# README
# Correlation structure
############# 

############# 
setwd("~/Altoida_VAMBN_paper")
rm(list=ls())
library(tidyverse)
library(plyr)
library(grid)
library(gridExtra)
library(Matching)
library(Hmisc)
library(corrplot)
library(ggpubr)
source('helper/plot_bn.R')
source('helper/clean_help.R')
source('helper/compare_dists.R')

############ data prep
name<-'main'
data_out<-paste0('data/data_out/',name)
data_all<-readRDS('data/data_condensed.rds')
virtual<-readRDS(paste0(data_out,'_VirtualPPts.rds'))
real<-readRDS(paste0(data_out,'_RealPPts.rds'))
data_all<-readRDS('data/data_out/data_all_imp.rds')
orig<-data_all[!grepl('stalone',names(data_all))]
orig<-orig %>% reduce(merge, by = 'SUBJID')
decRP_m<-read.csv('data/HI-VAE/reconRP.csv')
decRP_nm<-read.csv('data/HI-VAE/reconRP_nm.csv')
decRP = merge(decRP_m, decRP_nm)
decVP_m<-read.csv('data/HI-VAE/decodedVP.csv')
decVP_nm<-read.csv('data/HI-VAE/decodedVP_nm.csv')
decVP = merge(decVP_m, decVP_nm)

orig<-orig[,(colnames(orig) %in% colnames(decRP))]
for (col in colnames(orig)){
  if (is.factor(orig[,col])){#|grepl('COGT_RAVLT.learning|COGT_RAVLT.forgetting|COGT_FAQ|COGT_MMSE',col)
    orig[,col]<-factor(orig[,col])
    lvs<-levels(orig[,col])
    names(lvs)<-as.character(0:(length(lvs)-1))
    decRP[,col]<-factor(decRP[,col],labels=unlist(lvs[levels(factor(decRP[,col]))]))
    decVP[,col]<-factor(decVP[,col],labels=unlist(lvs[levels(factor(decVP[,col]))]))
  }
}

orig$SUBJID<-NULL
real$SUBJID<-NULL
decRP$SUBJID<-NULL
decVP$SUBJID<-NULL

orig<-orig[,colnames(decRP)]
orig<-orig[,sapply(colnames(orig),function(x) !is.factor(orig[,x]))]
decRP<-decRP[,sapply(colnames(decRP),function(x) !is.factor(decRP[,x]))]
decVP<-decVP[,sapply(colnames(decVP),function(x) !is.factor(decVP[,x]))]

if (!((all((colnames(decRP))==(colnames(orig))))&(all((colnames(decRP))==(colnames(decVP))))&(all((colnames(orig))==(colnames(decVP))))))
  warning('columns not matching!')


oc <- rcorr(as.matrix(orig))

rc <- rcorr(as.matrix(decRP))

vc <- rcorr(as.matrix(decVP))
library(RColorBrewer)

png("output_figures/Altoida_corr_orig.png") 
#tiff("output_figures/Altoida_corr_orig.tiff", width = 4, height = 4, units = 'in', res = 300)
corrplot(oc$r, type = "upper",method='color',main=paste0('Original'),tl.cex=0.25, tl.pos='n',mar=c(0,0,1.5,0))
dev.off() 

png("output_figures/Altoida_corr_decRP.png") 
corrplot(rc$r, type = "upper",method='color',main=paste0('Decoded real'),tl.cex=0.25, tl.pos='n',mar=c(0,0,1.5,0))
dev.off() 

png("output_figures/Altoida_corr_decVP.png") 
corrplot(vc$r, type = "upper",method='color',main=paste0('Decoded virtual'),tl.cex=0.25, tl.pos='n',mar=c(0,0,1.5,0))
dev.off() 


png("output_figures/Altoida_corr_orig.pdf") 
corrplot(oc$r, type = "upper",method='color',main=paste0('Original (Norm: ',round(norm(oc$r,type='F'),2),')'),tl.cex=0.25, tl.pos='n',mar=c(0,0,1.5,0))
dev.off() 

pdf("output_figures/Altoida_corr_decRP.pdf") 
corrplot(rc$r, type = "upper",method='color',main=paste0('Decoded real (Norm: ',round(norm(rc$r,type='F'),2),')'),tl.cex=0.25, tl.pos='n',mar=c(0,0,1.5,0))
dev.off() 

pdf("output_figures/Altoida_corr_decVP.pdf") 
corrplot(vc$r, type = "upper",method='color',main=paste0('Decoded virtual (Norm: ',round(norm(vc$r,type='F'),2),')'),tl.cex=0.25, tl.pos='n',mar=c(0,0,1.5,0))
dev.off() 




