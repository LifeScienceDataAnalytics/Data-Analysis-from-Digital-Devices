setwd("~/Altoida_VAMBN_paper")
rm(list=ls())
library(ggdendro)
library(ade4)
library(PTXQC)
data_all = readRDS("data/data_out/data_all_imp.rds")
##calculate the baseline DX
dxDf = data_all$stalone_demog_dx
baselineDX = dxDf[grep("\\_1", dxDf$SUBJID),]
altoidaDfSummary_vs2 = merge(baselineDX[grep("\\_1", baselineDX$SUBJID),], data_all$MMSE_Attention_Concentration,
                             data_all$MMSE_Language, data_all$MMSE_Memory_Recall, data_all$MMSE_Orientation,
                             data_all$MMSE_Working_Memory_Registration)
altoidaDfSummary_vs2 = list(baselineDX[grep("\\_1", baselineDX$SUBJID),], data_all$MMSE_Attention_Concentration,
                            data_all$MMSE_Language, data_all$MMSE_Memory_Recall, data_all$MMSE_Orientation,
                            data_all$MMSE_Working_Memory_Registration) %>% reduce(left_join, by = "SUBJID")
write.csv(altoidaDfSummary_vs2, "altoidaDfSummary_vs2.csv")
allDX = table(dxDf$SA_DX)
data_all$ARObjectFinding$BIT_AR_countdownFail = as.factor(data_all$ARObjectFinding$BIT_AR_countdownFail)
data_all$ARObjectFinding$DOT_AR_countdownFail = as.factor(data_all$ARObjectFinding$DOT_AR_countdownFail)
for(i in 1:length(data_all)){
  print(i)
  data_all[[i]]['SUBJID'] = NULL
  sv = data.frame()
  for(j in 1:ncol(data_all[[i]])){
   if(length(class(data_all[[i]][[j]]))>1){
      class(data_all[[i]][[j]]) = class(data_all[[i]][[j]])[2]
    }
    if(class(data_all[[i]][[j]]) == "factor"){
        sv[j,"type"] = "cat"
        sv[j,"dim"] = nlevels(data_all[[i]][[j]])
        sv[j,"nclass"] = nlevels(data_all[[i]][[j]])
        if(grepl("count", names(data_all[[i]][j]))==TRUE){
          print(names(data_all[[i]][j]))
          sv[j,"type"] = "cat"
          sv[j,"dim"] = nlevels(data_all[[i]][[j]])
          sv[j,"nclass"] = nlevels(data_all[[i]][[j]])
        }
    }
    if(class(data_all[[i]][[j]]) == "integer"){
        sv[j,"type"]= "ordinal"
        sv[j,"dim"] = length(unique(na.omit(data_all[[i]][[j]])))
        sv[j,"nclass"] = length(unique(na.omit(data_all[[i]][[j]])))
        
    }
   if(class(data_all[[i]][[j]]) == "numeric"){
        #print("bvbvbv")
        sv[j,"type"] = "real"
        sv[j,"dim"] = 1
        sv[j,"nclass"] = ""
        if(all(data_all[[i]][[j]] >= 0, na.rm = TRUE) == TRUE){
          #print("bvbkjljvbv")
          sv[j,"type"] = "pos"
          sv[j,"dim"] = 1
          sv[j,"nclass"] = ""
        }
      
    }
  }
   sv = sv[complete.cases(sv), ]
   #print(sv)
   if (length(which(is.na(data_all[[i]]), arr.ind = TRUE)) >0){
     write.table(sv,paste0('data/HI-VAE/data_python/',names(data_all[i]),'_types.csv'),sep=',',row.names = F,col.names = T,quote=T, na = "NaN")
     write.table(sv, paste0("GridSearch/data_python/",names(data_all[i]),'_types.csv'),  sep=',',row.names = F,col.names = T,quote=T, na = "NaN")
   }
  else{
    write.table(sv,paste0('data/HI-VAE/data_python_Notmissing/',names(data_all[i]),'_types.csv'),sep=',',row.names = F,col.names = T,quote=T, na = "NaN")
    write.table(sv, paste0("GridSearch/data_python_Notmissing/",names(data_all[i]),'_types.csv'),  sep=',',row.names = F,col.names = T,quote=T, na = "NaN")
  }
}

save.image("data_types.RData")

