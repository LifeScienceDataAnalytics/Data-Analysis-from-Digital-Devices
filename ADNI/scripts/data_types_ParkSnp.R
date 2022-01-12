# data_all = readRDS("data/data_condensed_ParkSNP.rds")
# data_all[] <- lapply(data_all[], factor)
# library(ggdendro)
# library(ade4)
# library(PTXQC)
# dataFrameNamesCombine = c()
# data_all_cp = data_all
# data_all_cp$SUBJID = NULL
# sv = data.frame()
# for(j in 1:ncol(data_all)){
#   sv[j,"type"] = "ordinal"
#   sv[j,"dim"] =
#   sv[j,"nclass"] = 
# }
# sv = sv[complete.cases(sv), ]

rm(list=ls())
data_all = readRDS("data/data_condensed_ParkSNP.rds")
library(ggdendro)
library(ade4)
library(PTXQC)
dataFrameNamesCombine = c()
for(i in 1:length(data_all)){
  data_all[[i]]$SUBJID = NULL
  #print(i)
  #dataFrameNames = LCSn(colnames(data_all[[1]]), min_LCS_length = 0)
  dataFrameNames = names(data_all[i])
  dataFrameNames = unlist(strsplit(dataFrameNames, split='.', fixed=TRUE))[1]
  print(dataFrameNames)
  sv = data.frame()
  #print(sv)
  #print(str(data_all[[i]]))
  for(j in 1:ncol(data_all[[i]])){
    #   print(j)
    if(length(class(data_all[[i]][[j]]))>1){
      #print(class(data_all[[i]][[j]]))
      class(data_all[[i]][[j]]) = class(data_all[[i]][[j]])[2]
      #print(class(data_all[[i]][[j]]))
    }
    #print(class(data_all[[i]][[j]])[2])
    #print(length(class(data_all[[i]][[j]])))
    if(class(data_all[[i]][[j]]) == "factor"){
      sv[j,"type"] = "cat"
      sv[j,"dim"] = nlevels(data_all[[i]][[j]])
      sv[j,"nclass"] = nlevels(data_all[[i]][[j]])
    }
    if(class(data_all[[i]][[j]]) == "integer"){
      sv[j,"type"]= "ordinal"
      sv[j,"dim"] = length(unique(data_all[[i]][[j]]))
      sv[j,"nclass"] = length(unique(data_all[[i]][[j]]))
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
  #print(sv)
  write.table(sv,paste0('data/HI-VAE/data_python_Park/',dataFrameNames,'_types.csv'),sep=',',row.names = F,col.names = T,quote=T, na = "NaN")
  write.table(sv, paste0("GridSearch/data_python_Park/",dataFrameNames,'_types.csv'),  sep=',',row.names = F,col.names = T,quote=T, na = "NaN")
  dataFrameNamesCombine = c(dataFrameNamesCombine, dataFrameNames)
}


#dataFrameNamesCombine = replace(dataFrameNamesCombine, dataFrameNamesCombine="", "stalone")
names(data_all) = dataFrameNamesCombine
save.image("~/Documents/Documents – IT-Admin’s MacBook Pro/ADNIVAMBN/VAMBNForADNI/data_types_Park.RData")
  

