############# README
# This is an additional analysis file for after bnet run. 
############# 

rm(list=ls())
library(bnlearn)
library(Rgraphviz)
library(ggplot2)
library(tidyverse)

################################# load discdata & bnet, fit bnet
name<-'main'
data_out<-paste0('data/data_out/',name)
finalBN<-readRDS(paste0(data_out,'_finalBN.rds'))
discdata<-readRDS(paste0(data_out,'_RealPPts.rds'))
discdata$SUBJID<-NULL
any(is.na(discdata))

fitted.bnlearn = bn.fit(finalBN,discdata)

######################### Counterfactuals
#https://stackoverflow.com/questions/41441812/prediction-with-cpdist-using-probabilities-as-evidence

n=round(10000/326) #(n*326 samples)

low<-fitted.bnlearn
low$SA_AGE_VIS1 = list(coef = fitted.bnlearn$SA_AGE_VIS1$coefficients-20, sd = fitted.bnlearn$SA_AGE_VIS1$sd)
high<-fitted.bnlearn
high$SA_AGE_VIS1 = list(coef = fitted.bnlearn$SA_AGE_VIS1$coefficients+20, sd = fitted.bnlearn$SA_AGE_VIS1$sd)

# age & UPDRS
BL<-cpdist(fitted.bnlearn, nodes=c('zcode_cogtest_VIS1'),evidence=T, method = "lw",n=length(discdata[,1])*n)
gend1<-cpdist(low, nodes=c('zcode_cogtest_VIS1'),evidence=T, method = "lw",n=length(discdata[,1])*n)
gend2<-cpdist(high, nodes=c('zcode_cogtest_VIS1'),evidence=T, method = "lw",n=length(discdata[,1])*n)
dat<-data.frame(dv=c(BL$zcode_cogtest_VIS1,gend1$zcode_cogtest_VIS1,gend2$zcode_cogtest_VIS1),level=factor(c(rep('No intervention',dim(gend1)[1]),rep('Age -20yrs',dim(gend1)[1]),rep('Age +20yrs',dim(gend2)[1]))))
pl<-ggplot(dat, aes(x=dv,fill=level))+geom_density(alpha=.2)+ggtitle('Age on CogTest (BL Visit)')
pl
write.csv(dat,paste0('data/data_out/counter_cog_age.csv'),row.names=FALSE)
ggsave('/Users/msood/Documents/ADNIVAMBN/VAMBNForADNI/output_figures/ADNI_CF_COG_code.png', pl, device = "png")
ggsave('/Users/msood/Documents/ADNIVAMBN/VAMBNForADNI/output_figures/ADNI_CF_COG_code.eps', pl, device = "eps")

# decoded (after running jupyter notebook)
dat<-read.csv('data/HI-VAE/CF_output.csv')
data_all<-readRDS('data/data_condensed.rds')
orig<-data_all[!grepl('stalone',names(data_all))]
orig<-orig %>% reduce(merge, by = 'SUBJID')
colnames(orig)<-gsub(' ','\\.',colnames(orig))
colnames(orig)<-gsub('[()]','\\.',colnames(orig))
orig<-orig[,grepl('UPDRS',colnames(orig))&grepl('_VIS00',colnames(orig))]
orig$Intervention<-'No intervention'
dat<-rbind(subset(dat,Intervention!='No intervention'),orig)
pl<-ggplot(dat, aes(x=UPDRS_UPDRS_VIS00,fill=Intervention))+geom_density(alpha=.2)+ggtitle('PPMI: Age on UPDRS (BL Visit)')+xlab('UPDRS total score')
#pl<-pl+theme(
#  plot.title = element_text(size=18),
#  axis.text = element_text(size=12),
#  axis.title=element_text(size=14)
#)
pl
ggsave('../output_figures/PPMI_CF_UPDRS.png', pl, device = "png")
ggsave('../output_figures/PPMI_CF_UPDRS.eps', pl, device = "eps")

