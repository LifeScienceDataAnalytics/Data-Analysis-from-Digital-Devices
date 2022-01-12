rm(list=ls())
library(tidyverse)
library(bnlearn)
library(PPMI)
#data(ppmimerge) #PD_MED_USE.post

finalBN<-readRDS('data/data_out/main_finalBN.rds')# load finalBN
real<-readRDS('data/data_out/main_RealPPts.rds')# load full real dataset 
real$SUBJID<-NULL
fitted.bnlearn<-readRDS('data/data_out/main_finalBN_fitted.rds')
tmp<-readRDS('data/data_out/data_all_imp.rds')# load demo dataset
demo_adni<-tmp[['demo_VIS1']]
demo_adni$SUBJID<-NULL

################# Extract PPMI data and save out jupyter data
# input (BL visit)
# - age,gender,UPDRS2,UPDRS3,ESS,[height,weight]* | *separate group
data_VIS1<- extractVariables(patients=getCohort('PD'),
  variables=c("SimpleGender", "ENROLL_AGE","UPDRS2","UPDRS3"), events='SC')
colnames(data_VIS1)<-c('SA_gender_VIS1','SA_age_VIS1','SA_UPDRS2_VIS1','SA_UPDRS3_VIS1')
data_VIS2<- extractVariables(patients=getCohort('PD'),
  variables=c("UPDRS2","UPDRS3","ESS"), events='BL')
colnames(data_VIS2)<-c('SA_UPDRS2_VIS2','SA_UPDRS3_VIS2','SA_EPWORTH_VIS2')
data<-cbind(data_VIS1,data_VIS2)

# demo data
# - note! BL visit here but SC in SP513, I think
demo<- extractVariables(patients=getCohort('PD'),variables=c("WGTKG", "HTCM"),events='BL')
colnames(demo)<-c('demo_WEIGHKG_VIS1','demo_X_HEIGHCM_VIS1')
demo[,colnames(demo_SP513)[!(colnames(demo_SP513) %in% colnames(demo))]] <- NA
demo<-demo[,colnames(demo_SP513)]
demo<-rbind(demo,matrix(NA,(NROW(demo_SP513)-NROW(demo)),length(colnames(demo)),dimnames=list(rep(NA,(NROW(demo_SP513)-NROW(demo))),colnames(demo))))
# save out demo data
write.table(which(is.na(demo), arr.ind=TRUE),'data/HI-VAE/PPMI/PPMI_demo_VIS1_missing.csv',sep=',',row.names = F,col.names = F,quote=F)
write.table(demo,'data/HI-VAE/PPMI/PPMI_demo_VIS1.csv',sep=',',row.names = F,col.names = F,quote=F, na = "NaN")
write.table(as.character(rownames(demo)),'data/HI-VAE/PPMI/PPMI_demo_VIS1_subj.csv',sep=',',row.names = F,col.names = T,quote=T, na = "NaN")

################# Run Jupyter Notebook for Demographic data (PPMI)
codes<-read.csv('data/HI-VAE/PPMI/PPMI_metaenc.csv')
codes<-codes[!grepl('NA',codes$SUBJID),]
codes$SUBJID<-NULL
if (!(all(codes$SUBJID==rownames(demo))))
  warning('Something went wrong!')
colnames(codes)<-gsub('PPMI_','',colnames(codes))
data<-cbind(data,codes)

################# Formatting
data$SA_gender_VIS1<-factor(ifelse(data$SA_gender_VIS1=='Female',2,1))
data$SA_age_VIS1<-factor(cut(floor(data$SA_age_VIS1), 
    breaks = c(-Inf, 53, 62, 72, Inf), 
    labels = c("<=52", "53-61", "62-70", ">=71"), 
    right = F))
#data<-data[complete.cases(data),]

################# sample conditional on evidence

# for every observation (evidence), get a single sample conditioned on it
pred<-rep(NA,NROW(data))
pred_trt<-rep(NA,NROW(data))
for (obs in 1:NROW(data)){
  # NO TREATMENT
  liste<-list(
    SA_gender_VIS1 = as.character(data[obs,'SA_gender_VIS1']),
    SA_age_VIS1 = as.character(data[obs,'SA_age_VIS1']),
    SA_UPDRS2_VIS1 = as.numeric(data[obs,'SA_UPDRS2_VIS1']),
    SA_UPDRS3_VIS1 = as.numeric(data[obs,'SA_UPDRS3_VIS1']),
    SA_UPDRS2_VIS2 = as.numeric(data[obs,'SA_UPDRS2_VIS2']),
    SA_UPDRS3_VIS2 = as.numeric(data[obs,'SA_UPDRS3_VIS2']),
    SA_EPWORTH_VIS2 = as.numeric(data[obs,'SA_EPWORTH_VIS2']),
    #scode_demo_VIS1 = as.character(data[obs,'scode_demo_VIS1']),
    #zcode_demo_VIS1 = as.numeric(data[obs,'zcode_demo_VIS1']),
    visitmiss_VIS15 = '0',
    SA_TREATMNT_VIS1 = 'Placebo'
  )
  if (!any(is.na(data[obs,])))
    pred[obs]<-mean(unlist(cpdist(fitted.bnlearn, 
         nodes=c('SA_UPDRS3_VIS15'),
         evidence=liste,method = "lw")),na.rm=T)
  
  # TREATMENT
  liste['SA_TREATMNT_VIS1']<-'Ropinirole'
  if (!any(is.na(data[obs,])))
    pred_trt[obs]<-mean(unlist(cpdist(fitted.bnlearn, nodes=c('SA_UPDRS3_VIS15'),
         evidence=liste,method = "lw")),na.rm=T)
}


################# Plots
# outcome - UPDRS 3 total post treatment
outcome<-extractVariables(patients=getCohort('PD'),
                          variables=c("UPDRS3.post","UPDRS3"),
                          events='V04') # 04/06/08/10
outcome1<-outcome$UPDRS3
outcome2<-outcome$UPDRS3.post
outcome<-outcome$UPDRS3.post
length(outcome[!is.na(outcome)])

dat<-data.frame(UPDRS3=c(pred,pred_trt),level=factor(c(rep('No Treatment',length(pred)),rep('Treatment',length(pred_trt)))))
pl<-ggplot(dat, aes(x=UPDRS3,fill=level))+geom_density(alpha=.2)+ggtitle('Simulated effect of Ropinirole in PPMI')
pl
#ggsave('../output_figures/CF_PPMI_pred.png', pl, device = "png")

dat<-data.frame(UPDRS3=c(outcome1,outcome2),level=factor(c(rep('No Treatment',length(outcome1)),rep('Treatment',length(outcome2)))))
pl<-ggplot(dat, aes(x=UPDRS3,fill=level))+geom_density(alpha=.2)+ggtitle('Observed L-DOPA effect in PPMI')
pl
#ggsave('../output_figures/CF_PPMI_real.png', pl, device = "png")

################# Joint figure
# outcome - UPDRS 3 total post treatment
df_PPMI<-extractVariables(patients=getCohort('PD'),
                          variables=c("UPDRS3.post","UPDRS3"),
                          events='V04') # 04/06/08/10
df_SP513<-data.frame(UPDRS3=pred,UPDRS3.post=pred_trt)
df_SP513$delt1<-df_SP513$UPDRS3.post-df_SP513$UPDRS3
df_PPMI$pred1<-df_PPMI$UPDRS3+df_SP513$delt1
dat<-data.frame(UPDRS3=c(df_PPMI$UPDRS3,df_PPMI$UPDRS3.post,df_PPMI$pred1),
                level=factor(c(rep('Baseline',length(df_PPMI$UPDRS3)),rep('Treatment',length(df_PPMI$UPDRS3.post)),rep('BL + Ropinirole effect',length(df_PPMI$pred1)))))
pl<-ggplot(dat, aes(x=UPDRS3,fill=level))+geom_density(alpha=.2)+ggtitle('Simulated effect of Ropinirole in PPMI')
pl
ggsave('../output_figures/CF_PPMI_SP513_cross.png', pl, device = "png")
ggsave('../output_figures/CF_PPMI_SP513_cross.eps', pl, device = "eps")

####################
####################
vis<-'V06'
# outcome - UPDRS 3 total post treatment
df_PPMI<-extractVariables(patients=getCohort('PD'),variables=c("UPDRS3.post","UPDRS3"),
                          events=vis) # 04/06/08/10
df_SP513<-data.frame(UPDRS3=pred,UPDRS3.post=pred_trt)
df_SP513$delt1<-df_SP513$UPDRS3.post-df_SP513$UPDRS3
df_PPMI$pred1<-df_PPMI$UPDRS3+df_SP513$delt1
dat<-data.frame(UPDRS3=c(df_PPMI$UPDRS3,df_PPMI$UPDRS3.post,df_PPMI$pred1),
                level=factor(c(rep('Baseline',length(df_PPMI$UPDRS3)),rep('Treatment',length(df_PPMI$UPDRS3.post)),rep('BL + Ropinirole effect',length(df_PPMI$pred1)))))
pl<-ggplot(dat, aes(x=UPDRS3,fill=level))+geom_density(alpha=.2)+ggtitle(paste0('Simulated effect of Ropinirole in PPMI (',vis,')'))
pl
ggsave(paste0('data/data_out/CF_PPMI/CF_PPMI_SP513_cross_',vis,'.png'), pl, device = "png")
ggsave(paste0('data/data_out/CF_PPMI/CF_PPMI_SP513_cross_',vis,'.eps'), pl, device = "eps")

