############# README
# This is the analysis file for the correlations in the real, decoded VP & RP data.
setwd("~/Altoida_VAMBN_paper")
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
library(grid)
library(gridExtra)
source('helper/plot_bn.R')
source('helper/clean_help.R')
source('helper/merge_data.R')
source('helper/discr.R')
source('helper/addnoise.R')
source('helper/make_bl_wl_altoida.R')
source('helper/simulate_VP.R')
source('helper/validate_VP.R')
source('helper/make_dummy.R')
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
ss = colnames(oc$r)[colSums(is.na(oc$r)) > 0]
oc<-rcorr(as.matrix(orig[,!names(orig)%in% ss]))
rc<-rcorr(as.matrix(decRP[,!names(orig)%in% ss]))
vc<-rcorr(as.matrix(decVP[,!names(orig)%in% ss]))
n_OP<-norm(oc$r,type='F')
n_RP<-norm(rc$r,type='F')
n_VP<-norm(vc$r,type='F')
relER_VP<-norm(oc$r-vc$r,type='F')/norm(oc$r,type='F')
relER_RP<-norm(oc$r-rc$r,type='F')/norm(oc$r,type='F')


# testing
data_corrs=data.frame(
  corr=c(as.vector(oc$r[upper.tri(oc$r)]),as.vector(rc$r[upper.tri(rc$r)]),as.vector(vc$r[upper.tri(vc$r)])),
  type=c(rep('real',length(as.vector(oc$r[upper.tri(oc$r)]))),
         rep('decoded real',length(as.vector(rc$r[upper.tri(rc$r)]))),
         rep('decoded virtual',length(as.vector(vc$r[upper.tri(vc$r)]))))
)
stats<-matrix(c('real',round(n_OP,2),'-',
                'decoded real',round(n_RP,2),round(relER_RP,3),
                'decoded virtual',round(n_VP,2),round(relER_VP,3)
)
,ncol=3,byrow=TRUE)
colnames(stats) <- c("Type","Norm","rel. error")
data_corrs$type<-factor(data_corrs$type,levels=c('real','decoded real','decoded virtual'))
pl<-ggplot(data_corrs, aes(x=corr,fill=type))+geom_density(alpha=.2)+
  xlab('Correlation')+ylab('density')+ggtitle('ALTOIDA')+
  annotation_custom(tableGrob(stats), xmin=0.15, ymax=10)
pl
ggsave('output_figures/Altoida_normcorrs.png', pl, device = "png", width=7.26, height = 4.35)
ggsave('output_figures/Altoida_normcorrs.eps', pl, device = "eps", width=7.26, height = 4.35)
pl<-ggplot(data_corrs, aes(x=type,y=corr,fill=type))+
  geom_violin(position=position_dodge(1)) + 
  geom_boxplot(position=position_dodge(1),width=0.05) + 
  xlab('type')+ylab('density')+ggtitle('ALTOIDA')+coord_cartesian(ylim=c(-1,1))+
  annotation_custom(tableGrob(stats), xmin=1.5, ymax=-0.35)
pl
ggsave(paste0('output_figures/ALTOIDA_normcorrs_violin.png'), pl, device = "png", width=7.26, height = 4.35)
ggsave(paste0('output_figures/ALTOIDA_normcorrs_violin.eps'), pl, device = "eps", width=7.26, height = 4.35)



pdf("data/data_out/corr_plots/orig.pdf") 
corrplot(oc$r, type = "upper",main=paste0('Original (Norm: ',round(norm(oc$r,type='F'),2),')'),tl.cex=0.25)
dev.off() 

pdf("data/data_out/corr_plots/decRP.pdf") 
corrplot(rc$r, type = "upper",main=paste0('Decoded real (Norm: ',round(norm(rc$r,type='F'),2),')'),tl.cex=0.25)
dev.off() 

pdf("data/data_out/corr_plots/decVP.pdf") 
corrplot(vc$r, type = "upper",main=paste0('Decoded virtual (Norm: ',round(norm(vc$r,type='F'),2),')'),tl.cex=0.25)
dev.off() 
