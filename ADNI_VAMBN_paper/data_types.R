rm(list=ls())
data_all = readRDS("~/ADNI_VAMBN_paper/data/data_out/data_all_imp.rds")

library(ggdendro)
library(ade4)
library(PTXQC)
library(gdata)
for(i in 1:length(data_all)){
  data_all[[i]]$SUBJID = NULL
  require(sqldf)
  dataFrameNames = names(data_all[i])
  sv = data.frame()
  for(j in 1:ncol(data_all[[i]])){
      if(length(class(data_all[[i]][[j]]))>1){
        class(data_all[[i]][[j]]) = class(data_all[[i]][[j]])[2]
      }
    if(class(data_all[[i]][[j]]) == "factor"){
        sv[j,"type"] = "cat"
        sv[j,"dim"] = nlevels(data_all[[i]][[j]])
        sv[j,"nclass"] = nlevels(data_all[[i]][[j]])
    }
    if(class(data_all[[i]][[j]]) == "integer"){
        print(names(data_all[[i]][j]))
        sv[j,"type"]= "ordinal"
        sv[j,"dim"] = length(unique(na.omit(data_all[[i]][[j]])))
        print(sv[j,"dim"])
        sv[j,"nclass"] = length(unique(na.omit(data_all[[i]][[j]])))
        print(sv[j,"nclass"])
    }
   if(class(data_all[[i]][[j]]) == "numeric"){
        #print("bvbvbv")
        sv[j,"type"] = "real"
        sv[j,"dim"] = 1
        sv[j,"nclass"] = ""
        if(all(data_all[[i]][[j]] >= 0, na.rm = TRUE) == TRUE){
          #print("bvbvbv")
          sv[j,"type"] = "pos"
          sv[j,"dim"] = 1
          sv[j,"nclass"] = ""
        }
   }
  }
  sv = sv[complete.cases(sv), ]
  if (length(which(is.na(data_all[[i]]), arr.ind = TRUE)) >0){
    write.table(sv,paste0('data/HI-VAE/data_python/',names(data_all[i]),'_types.csv'),sep=',',row.names = F,col.names = T,quote=T, na = "NaN")
    write.table(sv, paste0("GridSearch/data_python/",names(data_all[i]),'_types.csv'),  sep=',',row.names = F,col.names = T,quote=T, na = "NaN")
  }
  else{
    write.table(sv,paste0('data/HI-VAE/data_python_Notmissing/',names(data_all[i]),'_types.csv'),sep=',',row.names = F,col.names = T,quote=T, na = "NaN")
    write.table(sv, paste0("GridSearch/data_python_Notmissing/",names(data_all[i]),'_types.csv'),  sep=',',row.names = F,col.names = T,quote=T, na = "NaN")
  }
}

save.image("~/ADNI_VAMBN_paper/data_types.RData")

