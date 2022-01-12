############# README
# This is to analyse the likelihood under the model in terms of VP/RP
############# 
rm(list=ls())
library(tidyverse)
library(beepr)
library(arules)
library(mclust)
library(igraph)
library(rpart)
library(bnlearn) # hc might be overwritten by arules or some such package "bnlearn::hc" if so; not currently used though
library(parallel)
source('helper/plot_bn.R')
source('helper/simulate_VP.R')

############################## Settings and preprocessing

########## Name output files
name<-'main'
data_out<-paste0('data/data_out/',name)
scr<-"bic-cg" # 'bic' for basic autoencoder and fully discretized
mth<-"mle" # 'bayes' for basic autoencoder and fully discretized

# load model and data
bl<-read.csv(paste0(data_out,'_bl.csv'))
wl<-read.csv(paste0(data_out,'_wl.csv'))
finalBN<-readRDS(paste0(data_out,'_finalBN.rds'))
#bootBN<-readRDS(paste0(data_out,'_bootBN.rds'))
virtual<-readRDS(paste0(data_out,'_VirtualPPts.rds'))
real<-readRDS(paste0(data_out,'_RealPPts.rds'))
real$SUBJID<-NULL
fitted<-readRDS(paste0(data_out,'_finalBN_fitted.rds'))


# ######### plot graphs in cytoscape (careful about dashed line maps!)
# cyt_graph(finalBN,'final', F)
# cyt_graph(bootBN,'boot', T)
# source('network_plot_paper.R')

########### Loglikelihood of RP and VP data
# set.seed(123)
# virtual2 = simulate_VPs(real,finalBN,iterative=FALSE,scr,mth,wl,bl,n=250)
# virtual2 = virtual2[1:250,]
#fittedRP<-bn.fit(finalBN, real, method = mth)
#virtual<-virtual[!grepl("AUX|visitmiss|scode", colnames(virtual))]
#real<-real[!grepl("AUX|visitmiss|scode", colnames(real))]

virtual <- virtual %>% rename_at(vars(starts_with("zcode_Bcl")), 
                           funs(str_replace(., "Bcl-2_", "Bcl.2_")))
virtual <- virtual %>% rename_at(vars(starts_with("scode_Bcl")), 
                                 funs(str_replace(., "Bcl-2_", "Bcl.2_")))
virtual <- virtual %>% rename_at(vars(starts_with("AUX_Bcl")), 
                                 funs(str_replace(., "Bcl-2_", "Bcl.2_")))
virtual <- virtual %>% rename_at(vars(starts_with("zcode_Calcium")), 
                                 funs(str_replace(., "Calcium-", "Calcium.")))
virtual <- virtual %>% rename_at(vars(starts_with("AUX_Calcium")), 
                                 funs(str_replace(., "Calcium-", "Calcium.")))
virtual <- virtual %>% rename_at(vars(starts_with("scode_Calcium")), 
                                 funs(str_replace(., "Calcium-", "Calcium.")))


RP=logLik(fitted, real, by.sample=TRUE)
VP=logLik(fitted, virtual, by.sample=TRUE)

data_lik<-data.frame(
  likelihood=c(RP,VP),
  type=c(rep('real',length(RP)),rep('virtual',length(VP)))
)
print(paste('Mean RP:',mean(RP,rm.na=T),'Mean VP:',mean(VP[is.finite(VP)],rm.na=T)))
pl<-ggplot(data_lik, aes(x=likelihood,fill=type))+geom_density(alpha=.2)+xlab('Likelihood')+ylab('density')+ggtitle('ADNI')
pl
ggsave('output_figures/ADNI_likelihoods.png', pl, device = "png")
 