############# README
# This is the data imputation file. For HI-VAE there is no actual imputation, just saveout.
# Before this, run the R files clean_data->format_data->impute_aux (scripts with fixed settings).
# After this, run the autoencoder jupyter notebook (HI-VAE).
############## 

rm(list=ls())
library(missForest)
#setwd("~/Documents/PhDWork/BayesianNetworkAD/luise work/VirtualPatients-master/paper/")
source('helper/make_dummy.R') # create dummies for categorical variables
source('helper/clean_help.R') # check for constant variables
source('helper/fill_na.R') # fill na with mean or most frequent cat
source('helper/save_py.R') # fill na with mean or most frequent cat
data_out<-'data/data_out/'
data_out_py<-'data/HI-VAE/data_python_Park/'
load("data_types_Park.RData")
###################### Imputation & AUX

data_all<-readRDS(file = paste0("data/data_condensed_ParkSNP.rds"))
names(data_all) = dataFrameNamesCombine
data_aux=list()
for (datan in names(data_all)){ # for every variable group
  
  # load data & remove SUBJID
  data<-data_all[[datan]]
  #data1<-data_all[["stalone"]]
  #pt <- data1$SUBJID
  pt<-data$SUBJID
  data$SUBJID<-NULL
  # 
  #   if (!grepl('stalone_VIS6|stalone_VIS12|stalone_VIS24|snp_VIS1', datan)){
  #     # remove bad data
  #     data=data[,includeVar(data)]
  #     data=data[,rmMiss(data)]
  #   }
  
  ###################### AUX variables
  
  # make AUX columns and save in separate list (with SUBJID)
  # nms<-colnames(data)
  # if (grepl('stalone', datan)){
  #   dataux<-as.data.frame(sapply(as.data.frame(is.na(data)), as.numeric))
  #   dataux<-as.data.frame(sapply(dataux,factor))
  #   colnames(dataux)<-paste('AUX',nms,sep='_')
  # }else{
  #   dataux<-data.frame(factor(apply(data,1,function(x) as.numeric(all(is.na(x))))))
  #   colnames(dataux)<-paste('AUX',datan,sep='_')
  # }
  # 
  # # update AUX list
  # dataux$SUBJID<-pt
  # data_aux[[datan]]<-dataux
  
  ###################### Imputation
  #   print(datan)
  #   if (grepl('stalone', datan))
  #     data<-fillna(data) # if standalone data, mean and most frequent class imputation
  # # 
  #   if (!grepl('stalone_VIS6|stalone_VIS12|stalone_VIS24|snp_VIS1', datan)){
  #     # remove bad data
  #     data=data[,includeVar(data)]
  #     data=data[,rmMiss(data)]
  #   }
  
  # add ppt variable and update data list
  data$SUBJID <- pt
  data_all[[datan]]<-data
  
  # save out csv's of scaled continous and dummy coded categorical data for autoencoders
  pt<-data$SUBJID
  data$SUBJID<-NULL
  
  #missing write
  # if (!grepl('stalone', datan))
  #   write.table(which(is.na(data), arr.ind=TRUE),paste0(data_out_py,datan,'_missing.csv'),sep=',',row.names = F,col.names = F,quote=F)
  
  #data write
 
  #write.table(data,paste0(data_out_py,datan),sep=',',row.names = F,col.names = F,quote=F, na = "NaN")
  datan = unlist(strsplit(datan, split='.', fixed=TRUE))[1]
  write.table(as.character(pt),paste0('data/HI-VAE/python_names_park_snp/',datan,'_subj.csv'),sep=',',row.names = F,col.names = T,quote=T, na = "NaN")
  write.table(colnames(data),paste0('data/HI-VAE/python_names_park_snp/',datan,'_cols.csv'),sep=',',row.names = F,col.names = T,quote=T, na = "NaN")
}

# save all
saveRDS(data_all, file = paste0(data_out,'data_all_imp_park_snp.rds'))
#saveRDS(data_aux, file = paste0(data_out,'data_aux.rds'))

library(beepr)
beep()











###################### Imputation & AUX
#setwd("~/Documents/PhDWork/BayesianNetworkAD/luise work/VirtualPatients-master/paper/ADNI")
data_all<-readRDS(file = paste0("data/data_condensed_ParkSNP.rds"))

#for (datan in names(data_all)){ # for every variable group

  # load data & remove SUBJID
data<-data_all
pt <- data$SUBJID
data$SUBJID<-NULL
data$SUBJID <- pt
data_all<-data

write.table(as.character(pt),paste0('data/HI-VAE/python_names_park_snp/','ParkSnp_subj.csv'),sep=',',row.names = F,col.names = T,quote=T, na = "NaN")
write.table(colnames(data),paste0('data/HI-VAE/python_names_park_snp/','ParkSnp_cols.csv'),sep=',',row.names = F,col.names = T,quote=T, na = "NaN")


# save all
saveRDS(data_all, file = paste0(data_out,'data_all_imp_park_snp.rds'))
#saveRDS(data_aux, file = paste0(data_out,'data_aux.rds'))

library(beepr)
beep()


