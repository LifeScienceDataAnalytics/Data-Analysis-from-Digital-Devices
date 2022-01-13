setwd("~/Altoida_VAMBN_paper")
rm(list=ls())
library(openxlsx)
library(purrr)
library(rlist)
library(plyr)
library(stringr)
library(Hmisc)
library(readr)
##library to obtain longest common substring from multiple strings
library(PTXQC)
library(naniar)
library(rle)
library(missForest)
library(SGL)
load("~/ALTOIDA_VAMBN_paper/createbnAltoida.RData")
rm(list=setdiff(ls(),  c("targetNodes", "dataCP", "demog", "diagnostics", "cogDomains")))
cogDomains = gsub("_VIS1", "", cogDomains)
name<-'main'
data_out<-paste0('data/data_out/',name)
data_all<-readRDS('data/data_condensed.rds')
data_stalone<-list(data_all$stalone_cog,
                   data_all$stalone_demog_dx)
data_stalone<-data_stalone %>% reduce(merge, by = 'SUBJID')
orig<-data_all[!grepl('stalone_cog',names(data_all))]
orig<-orig %>% reduce(merge, by = 'SUBJID')
orig = merge(orig, data_stalone[,c("SUBJID",cogDomains)])
orig$SUBJID = NULL
orig_cog_mmse = orig[,c(cogDomains, grep("MMSE", colnames(orig), value = TRUE), "SA_DX")]
colnames(orig_cog_mmse) = gsub('SA_', '',colnames(orig_cog_mmse))
orig_cog_mmse = orig_cog_mmse[orig_cog_mmse$DX %in% c(0,1,2),]
orig_cog_mmse$DX[orig_cog_mmse$DX==2] <- 1
orig_cog_mmse$DX = factor(orig_cog_mmse$DX)
orig_cog_mmse[,grep("MMSE", colnames(orig_cog_mmse), value = TRUE)] = lapply(orig_cog_mmse[,grep("MMSE", colnames(orig_cog_mmse), value = TRUE)], as.factor)
#set.seed(1234)
write.csv(orig_cog_mmse, "altoidaDataClassifierCogDomains_LReg.csv")
