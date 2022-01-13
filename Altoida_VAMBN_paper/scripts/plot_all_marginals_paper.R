############# README
# This is the analysis file for the decoded VP data.
# Run after bnet.R
############# 
setwd("~/Altoida_VAMBN_paper")
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
data_all<-readRDS('data/data_out/data_all_imp.rds')

#data_stalone = data_all[[2]]
virtual<-readRDS(paste0(data_out,'_VirtualPPts.rds'))
virtual_nm<-readRDS(paste0(data_out,'_VirtualPPts_nm.rds'))
real<-readRDS(paste0(data_out,'_RealPPts.rds'))
orig<-data_all[!grepl('stalone',names(data_all))]
orig<-orig %>% reduce(merge, by = 'SUBJID')
decRP_m<-read.csv('data/HI-VAE/reconRP.csv')
decRP_nm<-read.csv('data/HI-VAE/reconRP_nm.csv')
decRP = merge(decRP_m, decRP_nm)
##mmse for classifier##
decRPMMSE = decRP[,grep("MMSE\\_|SUBJID", colnames(decRP), value = TRUE)]
write.csv(decRPMMSE, "decRPMMSE.csv")
decVP_m<-read.csv('data/HI-VAE/decodedVP.csv')
decVP_nm<-read.csv('data/HI-VAE/decodedVP_nm.csv')
decVP = merge(decVP_m, decVP_nm)
saveVirtual = cbind(virtual[,grep("SA", colnames(virtual), value= TRUE)], decVP)
saveVirtual = saveVirtual[,!grepl("AUX",colnames(saveVirtual)),]
saveVirtual$SUBJID = NULL
saveVirtual$ID <- seq.int(nrow(saveVirtual))
saveVirtual = saveVirtual %>% rename_at(vars(starts_with("SA")), 
                                        funs(str_replace(., "_VIS1", "")))
saveVirtual = saveVirtual %>% rename_at(vars(starts_with("SA")), 
                                        funs(str_replace(., "SA_", "")))

write.csv(saveVirtual, "virtualALTOIDA.csv")
orig<-orig[,(colnames(orig) %in% colnames(decRP))]
for (col in colnames(orig)){
  if (is.factor(orig[,col])){
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

orig$type<-'real'
virtual$type<-'virtual'
real$type<-'real'
decVP$type<-'decoded virtual'
decRP$type<-'decoded real'

real <- real %>% rename_at(vars(starts_with("zcode")), 
                           funs(str_replace(., "_VIS1", "")))
real <- real %>% rename_at(vars(starts_with("scode")), 
                                 funs(str_replace(., "_VIS1", "")))

########## Bnet generated codes
codes<-rbind.data.frame(virtual,real)
if (!all(sort(colnames(virtual))==sort(colnames(real)))){
  error('column names!')
}
codes$type<-factor(codes$type,levels=c('real','virtual'))

compare_dists(codes,'type','data/data_out/dist_plots_codes/',0.95)


## make violin plots of all other comparisons
#orig<-orig[,sort(colnames(orig))]
decRP<-decRP[,sort(colnames(decRP))]
orig<-orig[,sort(colnames(decVP))]
decVP<-decVP[,sort(colnames(decVP))]
all<-rbind(orig,decRP,decVP)
all$type<-factor(all$type,levels=c('real','decoded real','decoded virtual'))



compare_dists_paper(all,'type','data/data_out/dist_violin_joint/')
all = all[,grepl("MMSE|type", colnames(all))]
data = all
typecol = 'type'
folder = 'data/data_out/dist_violin_joint/'
gnames<-unique(data[,typecol])
vnames<-colnames(data[,!grepl('type',colnames(data))])
lv<-levels(factor(data[,typecol]))
group.1<-subset(data,eval(parse(text = typecol))==lv[1])
group.2<-subset(data,eval(parse(text = typecol))==lv[2])
group.3<-subset(data,eval(parse(text = typecol))==lv[3])

gnames<-unique(data[,typecol])
vnames<-colnames(data[,!grepl('type',colnames(data))])
alpha<-0.05/length(vnames)


for (col in vnames){
  set.seed(123)
  px =  dnorm(group.2[,col])
  py =  dnorm(group.3[,col])
  kld = KLD(px, py)
  klValue = round(kld$sum.KLD.py.px,4)
  dat<-data[,grepl(paste0(col,'|type'),colnames(data))]
  gd <- dat %>% group_by(type) %>% count
  plot<-ggplot(gd, aes(x=eval(parse(text = col)),y=freq,fill=type))+geom_bar(position="dodge", stat="identity")+scale_fill_brewer(palette="Dark2")+xlab(col)+ylab('Count')+ggtitle(paste('kldivergence:', klValue))
  ggsave(paste0(folder,col,'.png'), plot, device = "png", width=7.26, height = 4.35)
  ggsave(paste0(folder,col,'.eps'), plot, device = "eps", width=7.26, height = 4.35)
}

