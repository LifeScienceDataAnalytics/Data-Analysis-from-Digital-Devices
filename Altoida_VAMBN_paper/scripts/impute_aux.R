############# README
# This is the data imputation file. For HI-VAE there is no actual imputation, just saveout.
# Before this, run the R files clean_data->format_data->impute_aux (scripts with fixed settings).
# After this, run the autoencoder jupyter notebook (HI-VAE).
############## 
setwd("~/Altoida_VAMBN_paper")
rm(list=ls())
library(missForest)
source('helper/make_dummy.R') # create dummies for categorical variables
source('helper/clean_help.R') # check for constant variables
source('helper/fill_na.R') # fill na with mean or most frequent cat
source('helper/save_py.R') # fill na with mean or most frequent cat
data_out<-'data/data_out/'
data_out_py<-'data/HI-VAE/data_python/'
###################### Imputation & AUX
data_all<-readRDS(file = paste0("data/data_condensed.rds"))
dxDf = data_all$stalone_demog_dx

colnames(data_all$stalone_demog_dx) = paste("SA_", colnames(data_all$stalone_demog_dx), sep ="")
colnames(data_all$stalone_cog) = paste("SA_", colnames(data_all$stalone_cog), sep ="")
names(data_all$stalone_demog_dx)[names(data_all$stalone_demog_dx)=="SA_SUBJID"] = "SUBJID"
names(data_all$stalone_cog)[names(data_all$stalone_cog)=="SA_SUBJID"] = "SUBJID"
library(gtsummary)



data_aux=list()
for (datan in names(data_all)){ # for every variable group
  # load data & remove SUBJID
   data<-data_all[[datan]]
   pt <- data$SUBJID
   data$SUBJID<-NULL

  if (!grepl('stalone', datan)){ 
    print(datan)
    # remove bad data
    data=data[,includeVar(data)]
    data=data[,rmMiss(data)]
 
  ###################### AUX variables

  # make AUX columns and save in separate list (with SUBJID)
  nms<-colnames(data)
  if (grepl('stalone', datan)){
    dataux<-as.data.frame(sapply(as.data.frame(is.na(data)), as.numeric))
    dataux<-as.data.frame(sapply(dataux,factor))
    colnames(dataux)<-paste('AUX',nms,sep='_')
  }else{
    dataux<-data.frame(factor(apply(data,1,function(x) as.numeric(all(is.na(x))))))
    colnames(dataux)<-paste('AUX',datan,sep='_')
  }
   # update AUX list
   dataux$SUBJID<-pt
   data_aux[[datan]]<-dataux

   ###################### Imputation
  if (grepl('stalone', datan))
    data<-fillna(data) # if standalone data, mean and most frequent class imputation

  if (!grepl('stalone', datan)){
    # remove bad data
    data=data[,includeVar(data)]
    data=data[,rmMiss(data)]
  }

   # add ppt variable and update data list
   data$SUBJID <- pt
   data_all[[datan]]<-data

#   # save out csv's of scaled continous and dummy coded categorical data for autoencoders
   pt<-data$SUBJID
   ##changed Dec 2021##
   data$SUBJID<-NULL

#   #missing write
  if (!grepl('stalone', datan) & length(which(is.na(data), arr.ind = TRUE)) >0){
    print("memememe")
    print(datan)
    write.table(which(is.na(data), arr.ind=TRUE),paste0(data_out_py,datan,'_missing.csv'),sep=',',row.names = F,col.names = F,quote=F)
    write.table(data,paste0(data_out_py,datan,'.csv'),sep=',',row.names = F,col.names = F,quote=F, na = "NaN")
    write.table(which(is.na(data), arr.ind=TRUE),paste0("GridSearch/data_python/",datan,'_missing.csv'),sep=',',row.names = F,col.names = F,quote=F)
    write.table(data,paste0("GridSearch/data_python/",datan,'.csv'),sep=',',row.names = F,col.names = F,quote=F, na = "NaN")
    write.table(pt,paste0('data/HI-VAE/python_names/',datan,'_subj.csv'),sep=',',row.names = F,col.names = T,quote=T, na = "NaN")
    write.table(colnames(data),paste0('data/HI-VAE/python_names/',datan,'_cols.csv'),sep=',',row.names = F,col.names = T,quote=T, na = "NaN")
  }
   else{
     data_out_py_nm = 'data/HI-VAE/data_python_Notmissing/'
     write.table(data,paste0("GridSearch/data_python_Notmissing/",datan, '.csv'),sep=',',row.names = F,col.names = F,quote=F, na = "NaN")
     write.table(data,paste0(data_out_py_nm,datan, '.csv'),sep=',',row.names = F,col.names = F,quote=F, na = "NaN")
     write.table(pt,paste0('data/HI-VAE/python_names_Notmissing/',datan,'_subj.csv'),sep=',',row.names = F,col.names = T,quote=T, na = "NaN")
     write.table(colnames(data),paste0('data/HI-VAE/python_names_Notmissing/',datan,'_cols.csv'),sep=',',row.names = F,col.names = T,quote=T, na = "NaN")
   }
}

# save all
saveRDS(data_all, file = paste0(data_out,'data_all_imp.rds'))
saveRDS(data_aux, file = paste0(data_out,'data_aux.rds'))



