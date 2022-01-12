############# README
# This is the analysis file for the decoded VP data.
# Run after bnet.R
############# 
setwd("~/ADNI_VAMBN_paper")
rm(list=ls())
library(tidyverse)
library(plyr)
library(grid)
library(gridExtra)
library(Matching)
library(LaplacesDemon)
source('helper/plot_bn.R')
source('helper/clean_help.R')
source('helper/compare_dists.R')

############ data prep
name<-'main'
data_out<-paste0('data/data_out/',name)
data_all<-readRDS('data/data_out/data_all_imp_dig_cog.rds')
virtual<-readRDS(paste0(data_out,'_VirtualPPts.rds'))
real<-readRDS(paste0(data_out,'_RealPPts.rds'))
orig<-data_all[!grepl('stalone',names(data_all))]
orig<-orig %>% reduce(merge, by = 'SUBJID')
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

saveVirtual = cbind(virtual[,grep("SA", colnames(virtual), value= TRUE)], decVP)
saveVirtual = saveVirtual[,!grepl("AUX",colnames(saveVirtual)),]
saveVirtual$SUBJID = NULL
saveVirtual$ID <- seq.int(nrow(saveVirtual))
write.csv(saveVirtual, "virtualADNI.csv")
orig<-orig[,(colnames(orig) %in% colnames(decRP))]
for (col in colnames(orig)){
  if (is.factor(orig[,col])){#|grepl('COGT_RAVLT.learning|COGT_RAVLT.forgetting|COGT_FAQ|COGT_MMSE|COGT_ADAS11|COGT_ADAS13',col)
    orig[,col]<-factor(orig[,col])
    lvs<-levels(orig[,col])
    names(lvs)<-as.character(0:(length(lvs)-1))
    decRP[,col]<-factor(decRP[,col],labels=unlist(lvs[levels(factor(decRP[,col]))]))
  }
}

orig$SUBJID<-NULL
real$SUBJID<-NULL
decRP$SUBJID<-NULL
decVP$SUBJID<-NULL

orig$type<-'real'
virtual$type<-'virtual'
real$type<-'real'
decVP$type<-'decoded virtual'
decRP$type<-'decoded real'



########## Bnet generated codes
virtual<-virtual[!grepl("AUX|visitmiss|scode", colnames(virtual))]
real<-real[!grepl("AUX|visitmiss|scode", colnames(real))]

codes<-rbind(virtual,real)
if (!all(sort(colnames(virtual))==sort(colnames(real))))
  error('column names!')
codes$type<-factor(codes$type,levels=c('real','virtual'))

compare_dists(codes,'type','data/data_out/dist_plots_codes/',0.95)

## make violin plots of all other comparisons
orig<-orig[,sort(colnames(orig))]
decRP<-decRP[,sort(colnames(decRP))]
decVP$SUBJID<-NULL
decVP<-decVP[,sort(colnames(decVP))]
all<-rbind(orig,decRP,decVP)
  
compare_dists_paper(all,'type','data/data_out/dist_violin_joint/')

all = all[,grepl("MMSE|type", colnames(all))]
data = all
typecol = 'type'
folder = 'data/data_out/dist_violin_joint/'
gnames<-unique(data[,typecol])
vnames<-colnames(data[,!grepl('type',colnames(data))])
data[,typecol] = factor(data[,typecol], levels = c("real", "decoded real", "decoded virtual"))
lv<-levels(factor(data[,typecol]))
group.1<-subset(data,eval(parse(text = typecol))==lv[1])
group.2<-subset(data,eval(parse(text = typecol))==lv[2])
group.3<-subset(data,eval(parse(text = typecol))==lv[3])

alpha<-0.05/length(vnames)


for (col in vnames){
  set.seed(123)
  px =  dnorm(group.2[,col])
  py =  dnorm(group.3[,col])
  kld = KLD(px, py)
  klValue = round(kld$sum.KLD.py.px,4)
  dat<-data[,grepl(paste0(col,'|type'),colnames(data))]
  gd <- dat %>% group_by(type) %>% count
  plot<-ggplot(gd, aes(x=eval(parse(text = col)),y=freq,fill=type))+geom_bar(position="dodge", stat="identity")+scale_fill_brewer(palette="Dark2")+xlab(col)+ylab('Count')#+ggtitle(paste('Permutation Pval:',pval,'Alpha:',signif(alpha,digits = 3),'Sig:',sig))
  ggsave(paste0(folder,col,'.png'), plot, device = "png", width=7.26, height = 4.35)
  ggsave(paste0(folder,col,'.eps'), plot, device = "eps", width=7.26, height = 4.35)
}

