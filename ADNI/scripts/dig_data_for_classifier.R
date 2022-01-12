###
# Run nested cross-validation for sparse group lasso
####
rm(list=ls())
setwd("~/Documents/Documents_IT/paper/ADNI_VAMBN_paper_final")
library(caret)
library(PRROC)
library(SGL)
library(missForest)
library(parallel)
library(readr)
library(plyr)


data1 = read.csv("data/HI-VAE/reconRP_all_Dig.csv")
data2 = read.csv("data/HI-VAE/reconRP_all_Dig_last178.csv")
data2 = data2[!data2$SUBJID%in%data1$SUBJID,]
data = rbind.data.frame(data1, data2[!data2$SUBJID%in%data1$SUBJID,])
data_nm1 = read.csv("data/HI-VAE/reconRP_all_Dig_nm.csv")
data_nm2 = read.csv("data/HI-VAE/reconRP_all_Dig_las_nm_178.csv")
data_nm2 = data_nm2[!data_nm2$SUBJID%in%data_nm1$SUBJID,]
data_nm = rbind.data.frame(data_nm1, data_nm2)

#data_nm = read.csv("data/HI-VAE/reconRP_all_Dig_nm.csv")
data = as.data.frame(data)
data_nm = as.data.frame(data_nm)
data_merge = join(data_nm,data)
data_classifier_cn_mci = read.csv("adni_data_classifier_cn_mci.csv")
data_classifier_cn_mci$SUBJID = NULL
names(data_classifier_cn_mci)[names(data_classifier_cn_mci)=="RID"] ="SUBJID"
data_Dig_cn_mci = merge(data_merge, data_classifier_cn_mci)
data_Dig_cn_mci$SUBJID = NULL
data_Dig_cn_mci$ID <- seq.int(nrow(data_Dig_cn_mci))
names(data_Dig_cn_mci)[names(data_Dig_cn_mci) == "DX"] = "Class"
write.csv(data_Dig_cn_mci, "data_DigCNtoMCI.csv")

data_classifier_mci_dementia = read.csv("adni_data_classifier_mci_dementia.csv")
data_classifier_mci_dementia$SUBJID = NULL
names(data_classifier_mci_dementia)[names(data_classifier_mci_dementia)=="RID"] ="SUBJID"
data_Dig_mci_dementia = merge(data_merge, data_classifier_mci_dementia)
data_Dig_mci_dementia$SUBJID = NULL
data_Dig_mci_dementia$ID <- seq.int(nrow(data_Dig_mci_dementia))
names(data_Dig_mci_dementia)[names(data_Dig_mci_dementia) == "DX"] = "Class"
write.csv(data_Dig_mci_dementia, "data_DigMCItoDementia.csv")

