setwd("~/Documents/Documents_IT/ADNIVAMBN/VAMBNForADNI")
rm(list=ls())
load("test.RData")
test = readRDS(file = paste0("~/Documents/Documents_IT/ADNIVAMBN/VAMBNForADNI/data/testData.rds"))
visitData_complete=test
pt<-factor(visitData_complete$SUBJID)

##IDS taked for digital data##

colnames(visitData_complete)<-gsub('.bl','_VIS1',colnames(visitData_complete))
colnames(visitData_complete)<-gsub('.m06','_VIS6',colnames(visitData_complete))
colnames(visitData_complete)<-gsub('.m12','_VIS12',colnames(visitData_complete))
colnames(visitData_complete)<-gsub('.m18','_VIS18',colnames(visitData_complete))
colnames(visitData_complete)<-gsub('.m24','_VIS24',colnames(visitData_complete))
colnames(visitData_complete)<-gsub('.m36','_VIS36',colnames(visitData_complete))
colnames(visitData_complete)<-gsub('.m48','_VIS48',colnames(visitData_complete))
#colnames(visitData_complete)<-gsub('.m60','_VIS60',colnames(visitData_complete))
#colnames(visitData_complete)<-gsub('.m72','_VIS72',colnames(visitData_complete))
######## Variable groups
csf <- visitData_complete[,grep("ABETA|TAU|PTAU", colnames(visitData_complete), value = TRUE)]
#colnames(csf_VIS1)<-paste0('CSF_',colnames(csf_VIS1))
colnames(csf) = gsub("csf", "CSF", colnames(csf))
csf$SUBJID<-pt
csf_VIS1<-csf[,grepl('SUBJID|_VIS1$',colnames(csf))]
#csf_VIS12<-csf[,grepl('SUBJID|_VIS12',colnames(csf))]
#csf_VIS24<-csf[,grepl('SUBJID|_VIS24',colnames(csf))]


volume <- visitData_complete[,grep("Entorhinal|Fusiform|MidTemp|Hippocampus|Ventricles|WholeBrain", colnames(visitData_complete), value = TRUE)]
#colnames(volume)<-paste0('VOL_',colnames(volume))
colnames(volume) = gsub("brain", "VOL", colnames(volume))
volume$SUBJID<-pt
volume_VIS1<-volume[,grepl('SUBJID|_VIS1$',colnames(volume))]
volume_VIS6<-volume[,grepl('SUBJID|_VIS6',colnames(volume))]
volume_VIS12<-volume[,grepl('SUBJID|_VIS12',colnames(volume))]
#volume_VIS18<-volume[,grepl('SUBJID|_VIS18',colnames(volume))]
volume_VIS24<-volume[,grepl('SUBJID|_VIS24',colnames(volume))]
#volume_VIS36<-volume[,grepl('SUBJID|_VIS36',colnames(volume))]
#volume_VIS48<-volume[,grepl('SUBJID|_VIS48',colnames(volume))]



#brain68_VIS1<-visitData_complete[,c("SUBJID",grep('Right|Left',colnames(visitData_complete), value = TRUE))]

#cogtest <- visitData_complete[,grep("^CDRSB|^ADASQ4|^MMSE\\.|^RAVLT|^MOCA|^LDELTOTAL|^DIGITSCOR|^TRABSCOR|
#                                   ^mPACCdigit|^mPACCtrailsB|^FAQ\\.|ADAS[^0-9|^Q4]|EcogPt|FAQTOTAL", 
#                                  colnames(visitData_complete), value = TRUE)]


# colnames(cogtest)<-gsub("Cog", "COGT", colnames(cogtest))
# cogtest$SUBJID<-pt
# cogtest_VIS1<-cogtest[,grepl('SUBJID|_VIS1$',colnames(cogtest))]
# cogtest_VIS6<-cogtest[,grepl('SUBJID|_VIS6$',colnames(cogtest))]
# cogtest_VIS12<-cogtest[,grepl('SUBJID|_VIS12',colnames(cogtest))]
# cogtest_VIS18<-cogtest[,grepl('SUBJID|_VIS18',colnames(cogtest))]
# cogtest_VIS24<-cogtest[,grepl('SUBJID|_VIS24',colnames(cogtest))]
# cogtest_VIS36<-cogtest[,grepl('SUBJID|_VIS36',colnames(cogtest))]
# cogtest_VIS48<-cogtest[,grepl('SUBJID|_VIS48',colnames(cogtest))]
# cogtest_VIS60<-cogtest[,grepl('SUBJID|_VIS60',colnames(cogtest))]
# cogtest_VIS72<-cogtest[,grepl('SUBJID|_VIS72',colnames(cogtest))]

# imagingPET <- visitData_complete[,grep("FDG|PET|AV45", colnames(visitData_complete), value = TRUE)]
# imagingPET$SUBJID<-pt
# imagingPET_VIS1<-imagingPET[,grepl('SUBJID|_VIS1$',colnames(imagingPET))]
# imagingPET_VIS6<-imagingPET[,grepl('SUBJID|_VIS6$',colnames(imagingPET))]
# imagingPET_VIS12<-imagingPET[,grepl('SUBJID|_VIS12',colnames(imagingPET))]
# imagingPET_VIS24<-imagingPET[,grepl('SUBJID|_VIS24',colnames(imagingPET))]
# 


# SNPS
# snp_VIS1 <- visitData_complete[,grep("rs\\d+", colnames(visitData_complete), value = TRUE)]
# #colnames(snp_VIS1)<-paste0('SNP_',colnames(snp_VIS1))
# colnames(snp_VIS1)<-gsub("snp", "SNP", colnames(snp_VIS1))
# snp_VIS1$SUBJID<-pt
# snp_VIS1$PHS<-NULL

# pathways
# path_VIS1 <- visitData_complete[,grep("path", colnames(visitData_complete), value = TRUE)]
# #colnames(snp_VIS1)<-paste0('SNP_',colnames(snp_VIS1))
# colnames(path_VIS1)<-gsub("path", "PATH", colnames(path_VIS1))
# path_VIS1$SUBJID<-pt

## Standalone
snpsToBeIncludedAsStalone = reduce(snpsToBeIncludedAsStalone, full_join, by = "SUBJID")
snpsToBeIncludedAsStalone = snpsToBeIncludedAsStalone %>% rename_at(vars(-SUBJID), ~ paste0(.,"_VIS1"))
demogs <- visitData_complete[,grep("SUBJID|PT[^AU]|APOE4|AGE|PHS", colnames(visitData_complete), value = TRUE)]
#fdg <- visitData_complete[,grep("SUBJID|FDG", colnames(visitData_complete), value = TRUE)]
dx <- visitData_complete[,grep("SUBJID|DX", colnames(visitData_complete), value = TRUE)]
mmse <- visitData_complete[,grep("SUBJID|^MMSE\\_", colnames(visitData_complete), value = TRUE)]
faq <- visitData_complete[,grep("SUBJID|^FAQ", colnames(visitData_complete), value = TRUE)]
stalone<-visitData_complete[,c(colnames(demogs), colnames(dx), colnames(mmse), colnames(faq))]
stalone<-stalone[,grepl('_VIS|SUBJID$',colnames(stalone))]
#stalone = stalone[!grepl("FAQTOTAL.*", colnames(stalone))]
stalone$SUBJID<-pt
stalone = merge(stalone,snpsToBeIncludedAsStalone)
stalone$SUBJID = NULL
colnames(stalone)<-paste0('SA_',colnames(stalone))
stalone["SUBJID"] = pt
stalone$SA_APOE4_VIS1<-factor(stalone$SA_APOE4_VIS1)
stalone$SA_PTMARRY_VIS1  <-factor(stalone$SA_PTMARRY_VIS1)
stalone$SA_PTGENDER_VIS1  <-factor(stalone$SA_PTGENDER_VIS1)
#stalone$SA_PTEDUCAT_VIS1<-factor(ifelse(stalone$SA_PTEDUCAT_VIS1<16,'low','high'))

stalone_VIS1<-stalone[,grepl('SUBJID|_VIS1$',colnames(stalone))]
stalone_VIS6<-stalone[,grepl('SUBJID|_VIS6$',colnames(stalone))]
stalone_VIS12<-stalone[,grepl('SUBJID|_VIS12',colnames(stalone))]
stalone_VIS18<-stalone[,grepl('SUBJID|_VIS18',colnames(stalone))]
stalone_VIS24<-stalone[,grepl('SUBJID|_VIS24',colnames(stalone))]
stalone_VIS36<-stalone[,grepl('SUBJID|_VIS36',colnames(stalone))]
stalone_VIS48<-stalone[,grepl('SUBJID|_VIS48',colnames(stalone))]
#stalone_VIS60<-stalone[,grepl('SUBJID|_VIS60',colnames(stalone))]
#stalone_VIS72<-stalone[,grepl('SUBJID|_VIS72',colnames(stalone))]

snps_list = visitData_complete[,grep("SUBJID|\\_rs[0-9]+", colnames(visitData_complete), value=TRUE)]
library(ggdendro)
library(ade4)
library(PTXQC)
#snpFrameNames = LCSn(colnames(snps_list), min_LCS_length = 0)
snps_list$SUBJID = NULL
snps_list = split.default(snps_list, sub("\\_rs.*", "", names(snps_list)))
for(i in 1:length(snps_list)){
  snps_list[[i]]["SUBJID"] = pt
}

data_all_test<-list('csf_VIS1'=csf_VIS1,
               'volume_VIS1'=volume_VIS1,
               'volume_VIS6'=volume_VIS6,
               'volume_VIS12'=volume_VIS12,
               'volume_VIS24'=volume_VIS24,
               'stalone_VIS1'=stalone_VIS1,
               'stalone_VIS6'=stalone_VIS6,
               'stalone_VIS12'=stalone_VIS12,
               'stalone_VIS18'=stalone_VIS18,
               'stalone_VIS24'=stalone_VIS24,
               'stalone_VIS36'=stalone_VIS36,
               'stalone_VIS48'=stalone_VIS48)

data_all_with_snps_test = c(data_all_test, snps_list)

##add the snp mechanism data##


#merge the snps with the main file

saveRDS(data_all_with_snps_test, file = paste0("~/Documents/Documents_IT/ADNIVAMBN/VAMBNForADNI/data/data_condensed_testData.rds"))

############# README
# This is the data imputation file. For HI-VAE there is no actual imputation, just saveout.
# Before this, run the R files clean_data->format_data->impute_aux (scripts with fixed settings).
# After this, run the autoencoder jupyter notebook (HI-VAE).
############## 

rm(list=ls())
setwd("~/Documents/Documents_IT/ADNIVAMBN/VAMBNForADNI")
library(missForest)
source('helper/make_dummy.R') # create dummies for categorical variables
source('helper/clean_help.R') # check for constant variables
source('helper/fill_na.R') # fill na with mean or most frequent cat
source('helper/save_py.R') # fill na with mean or most frequent cat
data_out<-'data/data_out/'
data_out_py<-'data/HI-VAE/data_python/'

###################### Imputation & AUX
data_all<-readRDS(file = paste0("data/data_condensed_testData.rds"))
adniDf = data_all$stalone_VIS1

library(gtsummary)
# make dataset with a few variables to summarize
summary_adni <- adniDf %>% select("SA_PTGENDER_VIS1","SA_PTEDUCAT_VIS1", "SA_PTETHCAT_VIS1", 
                                  "SA_PTRACCAT_VIS1", "SA_PTMARRY_VIS1", "SA_APOE4_VIS1", "SA_PHS_VIS1",
                                  "SA_AGE_VIS1","SA_DX_VIS1")

# summarize the data with our package
summary_adni <- tbl_summary(summary_adni)

data_all$stalone_VIS1[,grep("MMSE|FAQ", colnames(data_all$stalone_VIS1))] <- lapply(data_all$stalone_VIS1[,grep("MMSE|FAQ", colnames(data_all$stalone_VIS1))], as.factor)
data_all$stalone_VIS6[,grep("MMSE|FAQ", colnames(data_all$stalone_VIS6))] <- lapply(data_all$stalone_VIS6[,grep("MMSE|FAQ", colnames(data_all$stalone_VIS6))], as.factor)
data_all$stalone_VIS12[,grep("MMSE|FAQ", colnames(data_all$stalone_VIS12))] <- lapply(data_all$stalone_VIS12[,grep("MMSE|FAQ", colnames(data_all$stalone_VIS12))], as.factor)
data_all$stalone_VIS18[,grep("MMSE|FAQ", colnames(data_all$stalone_VIS18))] <- lapply(data_all$stalone_VIS18[,grep("MMSE|FAQ", colnames(data_all$stalone_VIS18))], as.factor)
data_all$stalone_VIS24[,grep("MMSE|FAQ", colnames(data_all$stalone_VIS24))] <- lapply(data_all$stalone_VIS24[,grep("MMSE|FAQ", colnames(data_all$stalone_VIS24))], as.factor)
data_all$stalone_VIS36[,grep("MMSE|FAQ", colnames(data_all$stalone_VIS36))] <- lapply(data_all$stalone_VIS36[,grep("MMSE|FAQ", colnames(data_all$stalone_VIS36))], as.factor)
data_all$stalone_VIS48[,grep("MMSE|FAQ", colnames(data_all$stalone_VIS48))] <- lapply(data_all$stalone_VIS48[,grep("MMSE|FAQ", colnames(data_all$stalone_VIS48))], as.factor)
data_all$stalone_VIS1$`SA_Calpastatin-calpain_subgraph_rs10053056_VIS1` = as.factor(data_all$stalone_VIS1$`SA_Calpastatin-calpain_subgraph_rs10053056_VIS1`)
data_all$stalone_VIS1$`SA_Beta-Catenin_subgraph_rs1543332_VIS1` = as.factor(data_all$stalone_VIS1$`SA_Beta-Catenin_subgraph_rs1543332_VIS1`)
data_all$stalone_VIS1$`SA_Beta-Oxidation_of_Fatty_Acids_pathway_rs6444175_VIS1` = as.factor(data_all$stalone_VIS1$`SA_Beta-Oxidation_of_Fatty_Acids_pathway_rs6444175_VIS1`)
data_aux=list()
for (datan in names(data_all)){ # for every variable group
  # load data & remove SUBJID
  print(datan)
  data<-data_all[[datan]]
  pt<-data$SUBJID
  data$SUBJID<-NULL
  # 
  if (!grepl('stalone_VIS1|stalone_VIS6|stalone_VIS12|stalone_VIS18|stalone_VIS24|stalone_VIS36|
             stalone_VIS48', datan)){
    #print("meemansa1")
    # remove bad data
    data=data[,includeVar(data)]
    data=data[,rmMiss(data)]
  }
  if(ncol(data) == 0){
    data = NULL
    data_all[[datan]] = NULL
  }
  # 
  # ###################### AUX variables
  # 
  # # make AUX columns and save in separate list (with SUBJID)
  nms<-colnames(data)
  if(is.null(data) == FALSE){
    if (grepl('stalone', datan)){
      dataux<-as.data.frame(sapply(as.data.frame(is.na(data)), as.numeric))
      dataux<-as.data.frame(sapply(dataux,factor))
      colnames(dataux)<-paste('AUX',nms,sep='_')
    }else{
      dataux<-data.frame(factor(apply(data,1,function(x) as.numeric(all(is.na(x))))))
      colnames(dataux)<-paste('AUX',datan,sep='_')
    }
    
    # 
    # # update AUX list
    dataux$SUBJID<-pt
    data_aux[[datan]]<-dataux
    
    # ###################### Imputation
    if (grepl('stalone', datan))
      data<-fillna(data) # if standalone data, mean and most frequent class imputation
    
    if (!grepl('stalone_VIS1|stalone_VIS6|stalone_VIS12|stalone_VIS18|stalone_VIS24|stalone_VIS36|
             stalone_VIS48', datan)){
      # remove bad data
      #print("meemansa5")
      data=data[,includeVar(data)]
      data=data[,rmMiss(data)]
    }
    # 
    # # add ppt variable and update data list
    data$SUBJID <- pt
    data_all[[datan]]<-data
    # 
    # # save out csv's of scaled continuous and dummy coded categorical data for autoencoders
    pt<-data[,"SUBJID"]
    data$SUBJID<-NULL
    # 
    # #missing write
    # 
    # if (!grepl('stalone', datan))
    #   #print(which(is.na(data), arr.ind=TRUE))
    #   #print("meemansa6")
    #   write.table(which(is.na(data), arr.ind=TRUE),paste0(data_out_py,datan,'_missing.csv'),sep=',',row.names = F,col.names = F,quote=F)
    # 
    # # #data write
    # if (!grepl('stalone', datan))
    #   #print("meemansa7")
    #   write.table(data,paste0(data_out_py,datan,'.csv'),sep=',',row.names = F,col.names = F,quote=F, na = "NaN")
    # 
    # write.table(as.character(pt),paste0('data/HI-VAE/python_names/',datan,'_subj.csv'),sep=',',row.names = F,col.names = T,quote=T, na = "NaN")
    # write.table(colnames(data),paste0('data/HI-VAE/python_names/',datan,'_cols.csv'),sep=',',row.names = F,col.names = T,quote=T, na = "NaN")
  }
}

# save all
saveRDS(data_all, file = paste0(data_out,'data_all_imp_testData.rds'))
saveRDS(data_aux, file = paste0(data_out,'data_aux_testData.rds'))

library(beepr) 
beep()




rm(list=ls())
library(tidyverse)
library(bnlearn)
##load the workspace from  the model trained on ALTOIDA data##
load("~/Documents/Documents_IT/VAMBN/createbnAltoida_Nov26.RData")
load("~/Documents/Documents_IT/VAMBN/randomForestPred_Nov26.RData")
##keep the required target nodes that needs to be simulated in ADNI data#
rm(list=setdiff(ls(),  c("targetNodes", "nrmseDf", "regResultsFinal", "discdata", "cogDomainsAltoida",
                         "mmseMeasures")))
#keep(targetNodes, sure=TRUE)
targetNodes = unique(targetNodes)

##set the working directory to the ALTOIDA model working directory
setwd("~/Documents/Documents_IT/VAMBN/")
#load the model trained on ALTOIDA data#
finalBN<-readRDS('data/data_out/main_finalBN.rds')# load finalBN
real<-readRDS('data/data_out/main_RealPPts.rds')# load full real dataset 
real$SUBJID<-NULL
fitted.bnlearn<-readRDS('data/data_out/main_finalBN_fitted.rds')
tmp<-readRDS('data/data_out/data_all_imp.rds')# load demo dataset



#simulate the digital nodes for ADNI##
func.pred.digital.nodes = function(data, target){
  extractVariables = grep(".AGE.*|.*MMSE.*|.*DX.*|.*PTGENDER.*|SUBJID", colnames(data), value = TRUE)
  data = data[,extractVariables]
  names(data)[names(data)%in%c("SA_PTGENDER_VIS1","SA_PTGENDER_VIS6","SA_PTGENDER_VIS12","SA_PTGENDER_VIS18","SA_PTGENDER_VIS24","SA_PTGENDER_VIS36","SA_PTGENDER_VIS48")] = "SA_gender_VIS1"
  names(data)[names(data)%in%c("SA_AGE_VIS1","SA_AGE_VIS6","SA_AGE_VIS12","SA_AGE_VIS18","SA_AGE_VIS24","SA_AGE_VIS36","SA_AGE_VIS48")] = "SA_age_VIS1"
  names(data)[names(data)%in%c("SA_DX_VIS1","SA_DX_VIS6","SA_DX_VIS12","SA_DX_VIS18","SA_DX_VIS24","SA_DX_VIS36","SA_DX_VIS48")] = "SA_MCI_VIS1"
  #data$SA_gender_VIS1 = as.factor(data$SA_gender_VIS1)
  data$SA_gender_VIS1 = as.character(data$SA_gender_VIS1)
  data$SA_gender_VIS1[data$SA_gender_VIS1=="Female"] <- "2"
  data$SA_gender_VIS1[data$SA_gender_VIS1=="Male"] <- "1"
  data$SA_gender_VIS1 = factor(data$SA_gender_VIS1)
  data$SA_MCI_VIS1 = factor(data$SA_MCI_VIS1)
  if("CN" %in% unique(data$SA_MCI_VIS1)){
    levels(data$SA_MCI_VIS1)[levels(data$SA_MCI_VIS1)=="CN"] <- "0"
  }
  if("MCI" %in% unique(data$SA_MCI_VIS1)){
    levels(data$SA_MCI_VIS1)[levels(data$SA_MCI_VIS1)=="MCI"] <- "1"
  }
  if("Dementia" %in% unique(data$SA_MCI_VIS1)){
    levels(data$SA_MCI_VIS1)[levels(data$SA_MCI_VIS1)=="Dementia"] <- "0"
  }
  data$SA_MCI_VIS1 = factor(data$SA_MCI_VIS1, levels = c("0", "1"))
  #digData = data.frame()
  
  # nrmse_bn = nrmseDf[nrmseDf$target == targetNodes[i], "Mean"]
  # print(nrmse_bn)
  # nrmse_rf = regResultsFinal[regResultsFinal$target == targetNodes[i], "nrmseMean"]
  # print(nrmse_rf)
  pt = data["SUBJID"]
  data = data[,c("SA_gender_VIS1","SA_age_VIS1", "SA_MCI_VIS1", grep("^MMSE", colnames(data), value = TRUE))]
  #if(nrmse_bn<nrmse_rf){
  #print("bayesian")
  set.seed(123)
  data[target] = predict(fitted.bnlearn, data = data, node = target, method = "bayes-lw")
  
  #}
  # else{
  #   print("rf")
  #   set.seed(123)
  #   #print(colnames(data))
  #   discdata <- discdata %>% 
  #     select(cogDomainsAltoida, mmseMeasures, "age", "gender", "AB+",
  #            "MCI", targetNodes[1])
  #   print(colnames(discdata))
  #   regr <- randomForest(discdata[,targetNodes[1]] ~ ., data = discdata, mtry = 3, 
  #                        importance = TRUE, na.action = na.omit)
  #   predict(regr, data)
  # }
  return(list(data,pt))
}

func.all.time.points = function(visData, visit){
  digDataVis = c()
  for(target in targetNodes){
    data = func.pred.digital.nodes(visData, target)
    digDataVis = c(digDataVis, data[[1]][target])
  }
  digDataVis = as.data.frame(digDataVis)
  digDataVis = cbind.data.frame(digDataVis, data[[2]])
  digDataVis = digDataVis %>% rename_at(vars(-SUBJID), ~ paste0(.,"_VIS", visit))
  return(digDataVis)
}

##simulate the digital nodes at each time point##
impADNIdata = readRDS("~/Documents/Documents_IT/ADNIVAMBN/VAMBNForADNI/data/data_out/data_all_imp_testData.rds")
inputDataVIS1 = impADNIdata$stalone_VIS1
inputDataVIS1[,grep("^SA_MMSE", colnames(impADNIdata$stalone_VIS1))] <- lapply(inputDataVIS1[,grep("^SA_MMSE", colnames(inputDataVIS1))], as.numeric)
digDataVIS1 = func.all.time.points(inputDataVIS1,1)


inputDataVIS6 = impADNIdata$stalone_VIS6
inputDataVIS6["SA_PTGENDER_VIS6"] = inputDataVIS1$SA_PTGENDER_VIS1 
inputDataVIS6["SA_AGE_VIS6"] = inputDataVIS1$SA_AGE_VIS1 
inputDataVIS6[,grep("^SA_MMSE", colnames(impADNIdata$stalone_VIS6))] <- lapply(inputDataVIS6[,grep("^SA_MMSE", colnames(inputDataVIS6))], as.numeric)
digDataVIS6 = func.all.time.points(inputDataVIS6,6)

inputDataVIS12 = impADNIdata$stalone_VIS12
inputDataVIS12["SA_PTGENDER_VIS12"] = inputDataVIS1$SA_PTGENDER_VIS1 
inputDataVIS12["SA_AGE_VIS12"] = inputDataVIS1$SA_AGE_VIS1 
inputDataVIS12[,grep("^SA_MMSE", colnames(impADNIdata$stalone_VIS12))] <- lapply(inputDataVIS12[,grep("^SA_MMSE", colnames(inputDataVIS12))], as.numeric)
digDataVIS12 = func.all.time.points(inputDataVIS12,12)


inputDataVIS18 = impADNIdata$stalone_VIS18
inputDataVIS18["SA_PTGENDER_VIS18"] = inputDataVIS1$SA_PTGENDER_VIS1 
inputDataVIS18["SA_AGE_VIS18"] = inputDataVIS1$SA_AGE_VIS1 
inputDataVIS18[,grep("^SA_MMSE", colnames(impADNIdata$stalone_VIS18))] <- lapply(inputDataVIS18[,grep("^SA_MMSE", colnames(inputDataVIS18))], as.numeric)
digDataVIS18 = func.all.time.points(inputDataVIS18,18)


inputDataVIS24 = impADNIdata$stalone_VIS24
inputDataVIS24["SA_PTGENDER_VIS24"] = inputDataVIS1$SA_PTGENDER_VIS1 
inputDataVIS24["SA_AGE_VIS24"] = inputDataVIS1$SA_AGE_VIS1 
inputDataVIS24[,grep("^SA_MMSE", colnames(impADNIdata$stalone_VIS24))] <- lapply(inputDataVIS24[,grep("^SA_MMSE", colnames(inputDataVIS24))], as.numeric)
digDataVIS24 = func.all.time.points(inputDataVIS24,24)

inputDataVIS36 = impADNIdata$stalone_VIS36
inputDataVIS36["SA_PTGENDER_VIS36"] = inputDataVIS1$SA_PTGENDER_VIS1
inputDataVIS36["SA_AGE_VIS36"] = inputDataVIS1$SA_AGE_VIS1 
inputDataVIS36[,grep("^SA_MMSE", colnames(impADNIdata$stalone_VIS36))] <- lapply(inputDataVIS36[,grep("^SA_MMSE", colnames(inputDataVIS36))], as.numeric)
digDataVIS36 = func.all.time.points(inputDataVIS36,36)

# inputDataVIS48 = impADNIdata$stalone_VIS48
# inputDataVIS48["SA_PTGENDER_VIS48"] = inputDataVIS1$SA_PTGENDER_VIS1
# inputDataVIS48["SA_AGE_VIS48"] = inputDataVIS1$SA_AGE_VIS1
# inputDataVIS48[,grep("^SA_MMSE", colnames(impADNIdata$stalone_VIS48))] <- lapply(inputDataVIS48[,grep("^SA_MMSE", colnames(inputDataVIS48))], as.numeric)
# digDataVIS48 = func.all.time.points(inputDataVIS48,48)

# inputDataVIS60 = impADNIdata$stalone_VIS60
# inputDataVIS60["SA_PTGENDER_VIS60"] = inputDataVIS1$SA_PTGENDER_VIS1
# inputDataVIS60["SA_AGE_VIS6"] = inputDataVIS1$SA_AGE_VIS1 
# digDataVIS60 = func.all.time.points(inputDataVIS60,60)
# 
# inputDataVIS72 = impADNIdata$stalone_VIS72
# inputDataVIS72["SA_PTGENDER_VIS72"] = inputDataVIS1$SA_PTGENDER_VIS1
# inputDataVIS72["SA_AGE_VIS72"] = inputDataVIS1$SA_AGE_VIS1 
# digDataVIS72 = func.all.time.points(inputDataVIS72,72)

digDataAll = cbind.data.frame(digDataVIS1, digDataVIS6, digDataVIS12, digDataVIS18, digDataVIS24, 
                              digDataVIS36)

digDataAll= digDataAll[!duplicated(as.list(digDataAll))]
digDataAll$SUBJID

##add the digital data to data_all
data_all = readRDS("~/Documents/Documents_IT/ADNIVAMBN/VAMBNForADNI/data/data_out/data_all_imp_testData.rds")
digDataVIS1 = digDataVIS1 %>% rename_at(vars(-SUBJID), ~ paste0("SA_",.))
data_all$stalone_VIS1 = merge(data_all$stalone_VIS1, digDataVIS1)

digDataVIS6 = digDataVIS6 %>% rename_at(vars(-SUBJID), ~ paste0("SA_",.))
data_all$stalone_VIS6 = merge(data_all$stalone_VIS6, digDataVIS6)

digDataVIS12 = digDataVIS12 %>% rename_at(vars(-SUBJID), ~ paste0("SA_",.))
data_all$stalone_VIS12 = merge(data_all$stalone_VIS12, digDataVIS12)

digDataVIS18 = digDataVIS18 %>% rename_at(vars(-SUBJID), ~ paste0("SA_",.))
data_all$stalone_VIS18 = merge(data_all$stalone_VIS18, digDataVIS18)

digDataVIS24 = digDataVIS24 %>% rename_at(vars(-SUBJID), ~ paste0("SA_",.))
data_all$stalone_VIS24 = merge(data_all$stalone_VIS24, digDataVIS24)

digDataVIS36 = digDataVIS36 %>% rename_at(vars(-SUBJID), ~ paste0("SA_",.))
data_all$stalone_VIS36 = merge(data_all$stalone_VIS36, digDataVIS36)

#digDataVIS48 = digDataVIS48 %>% rename_at(vars(-SUBJID), ~ paste0("SA_",.))
#data_all$stalone_VIS48 = merge(data_all$stalone_VIS48, digDataVIS48)

# digDataVIS60 = digDataVIS60 %>% rename_at(vars(-SUBJID), ~ paste0("SA_",.))
# data_all$stalone_VIS60 = merge(data_all$stalone_VIS60, digDataVIS60)
# 
# digDataVIS72 = digDataVIS72 %>% rename_at(vars(-SUBJID), ~ paste0("SA_",.))
# data_all$stalone_VIS72 = merge(data_all$stalone_VIS72, digDataVIS72)
save.image("~/Documents/Documents_IT/ADNIVAMBN/VAMBNForADNI/digitalDataPred_testData.RData")
setwd("~/Documents/Documents_IT/ADNIVAMBN/VAMBNForADNI")
data_out<-'data/data_out/'
saveRDS(data_all, file = paste0(data_out,'data_all_imp_dig_testData.rds'))
######### return column i's where variance is not 0
includeVar_testData <- function(dat) {
  out <- lapply(dat, function(x) length(unique(x[!is.na(x)])))
  want <- which(out > 1)
  unlist(want)
}

merge_data_test<-function(){
  # Dir
  data_out<-'data/data_out/'
  #(standalone)
  data_all<-readRDS(file = paste0(data_out,'data_all_imp_dig_testData.rds'))
  data_stalone<-list(data_all[['stalone_VIS1']],
                     data_all[['stalone_VIS6']],
                     data_all[['stalone_VIS12']],
                     data_all[['stalone_VIS18']],
                     data_all[['stalone_VIS24']],
                     data_all[['stalone_VIS36']])
                     #data_all[['stalone_VIS48']])
  #data_all[['stalone_VIS60']],
  #data_all[['stalone_VIS72']])
  
  data_stalone<-data_stalone %>% reduce(merge, by = 'SUBJID')
  #(aux)
  data_aux<-readRDS('data/data_out/data_aux_testData.rds')
  data_aux<-data_aux %>% reduce(merge, by = 'SUBJID')
  data_aux<-as.data.frame(lapply(data_aux,factor))
  
  #(meta)
  data_meta<-read.csv('data/HI-VAE/metaenc.csv')
  
  # merge all
  
  data<-list(data_meta,data_aux,data_stalone) %>% reduce(merge, by = 'SUBJID')
  #data_aux = as.data.frame(data_aux)
  #data_meta = as.data.frame(data_meta)
  #data_stalone = as.data.frame(data_stalone)
  #data<-cbind.data.frame(data_meta,data_aux,data_stalone)
  #flag 0 var cols
  print(colnames(data)[-includeVar_testData(data)])
  data<-data[includeVar_testData(data)]
  
  name<-'data_final'
  #saveRDS(data,paste0(data_out,name,'.rds'))
  return(data)
}

data<-merge_data_test() 
data = write.csv(data, "testData.csv")  
  
  
rm(list=ls())
testDataADNI <- read_csv("testData.csv")
testDataADNI$X1 = NULL
dig = grep("BIT|DOT", colnames(testDataADNI), value = TRUE)
testDataADNI <- testDataADNI %>% rename_at(dig, 
                                           funs(str_replace(., "\\_VIS1", "")))
testDataADNI[,grep("VOL", colnames(testDataADNI), value = TRUE)] = NULL
testDataADNI[,grep("CSF", colnames(testDataADNI), value = TRUE)] = NULL
demogs = grep("SA_PT|SA_PH|SA_AGE|SA_APOE4|rs.*VIS", colnames(testDataADNI), value = TRUE)
testDataADNI <- testDataADNI %>% rename_at(demogs, 
                                         funs(str_replace(., "\\_VIS1", "")))
demogs = grep("SA_PT|SA_PH|SA_AGE|SA_APOE4|rs", colnames(testDataADNI), value = TRUE)
testDataADNI <- testDataADNI %>% rename_at(demogs, 
                                         funs(str_replace(., "SA_", "")))
names(testDataADNI)[names(testDataADNI)=="ID"] = "RID"
testDataADNI = testDataADNI[,!grepl("AUX|scode",colnames(testDataADNI)),]
sel <- grep(".*\\_VIS", colnames(testDataADNI), value = TRUE)
testDataADNI <- testDataADNI %>% rename_at(sel, 
                                         funs(str_replace(., "\\_VIS", "\\.")))
sel <- grep("SA", colnames(testDataADNI), value = TRUE)
testDataADNI <- testDataADNI %>% rename_at(sel, 
                                         funs(str_replace(., "SA_", "")))
sel <- grep("zcode", colnames(testDataADNI), value = TRUE)
testDataADNI <- testDataADNI %>% rename_at(sel, 
                                         funs(str_replace(., "zcode_", "")))
testDataADNI = testDataADNI[,!grepl("vol|csf",colnames(testDataADNI)),]
testDataADNI = as.data.frame(testDataADNI)
testDataADNI$MMSE_Working_Memory_Registration.36 = NULL
testDataADNI$MMSE_Working_Memory_Registration.6 = NULL
testDataADNI_long = reshape(testDataADNI,
                           varying = grep(".*\\.[0-9]+$", colnames(testDataADNI)),
                           direction = "long",
                           #v.names = "Value",
                           idvar = "RID",
                           timevar = "VISCODE")
rownames(testDataADNI_long) = 1:nrow(testDataADNI_long)
modify.educat <- function(value){
  ifelse(value <= 9, value<-"Low",value<-"High")
}
baseline.adni <- testDataADNI_long[,grep("PT|PH|AGE|APOE4|AGE|rs.*", colnames(testDataADNI_long), value = TRUE)]

### Grouping of patients by progression patterns ###
aggDiagnosis <- ddply(testDataADNI_long, .(RID), summarize, DX = toString(DX))

### Filtering last MCI/first Dementia points ###
mci <- testDataADNI_long[testDataADNI_long$DX == "MCI",]
dementia <- testDataADNI_long[testDataADNI_long$DX == "Dementia",]

last_mci <- mci %>% group_by(RID) %>% slice(which.max(VISCODE))
first_dementia <- dementia %>% group_by(RID) %>% slice(which.min(VISCODE))

threshold_df <- rbind(last_mci,first_dementia)
table(threshold_df$DX)
adni.split.count <- ddply(threshold_df, .(RID), summarize, DX = toString(DX))
converters <- adni.split.count[adni.split.count$DX == "MCI, Dementia",]
mciToADConverters <- threshold_df[which(threshold_df$RID %in% converters$RID),]

### Label Encoding ###
mciToADConverters <- threshold_df[which(threshold_df$RID %in% converters$RID),]

mciToADConverters$DX <- factor(mciToADConverters$DX,
                               levels = c("Dementia","MCI"),
                               labels = c(0,1))

mciToADConverters$PTGENDER <- factor(mciToADConverters$PTGENDER,
                                     levels = c("Female","Male"),
                                     labels = c(0,1))
mciToADConverters$PTGENDER <- as.factor(mciToADConverters$PTGENDER)
mciToADConverters$PTETHCAT <- factor(mciToADConverters$PTETHCAT,
                                     levels = c("Hisp/Latino","Not Hisp/Latino","Unknown"),
                                     labels = c(0,1,2))
mciToADConverters$PTETHCAT <- as.factor(mciToADConverters$PTETHCAT)
mciToADConverters$PTMARRY <- factor(mciToADConverters$PTMARRY,
                                    levels = c("Divorced","Married","Never married","Widowed"),
                                    labels = c(0,1,2,3))
mciToADConverters$PTMARRY <- as.factor(mciToADConverters$PTMARRY)

set.seed(1)
rows <- sample(nrow(mciToADConverters))
shuffleTest <- mciToADConverters[rows,]
shuffleTest = shuffleTest[,!grepl("subgraph|Receptors|transduction|Beta-Oxidation|AD_T2DM_SNPs",colnames(shuffleTest)),]
names(shuffleTest)[names(shuffleTest)=="SUBJID"] = "RID"
write_csv(shuffleTest,"testData_processed.csv")
sadniVirtual = read.csv("ADNI_virtual.csv")
adniVirtual= adniVirtual[,colnames(shuffleTest)]
write.csv(adniVirtual, "ADNI_virtual_training.csv")
  