setwd("~/Altoida_VAMBN_paper")


rm(list=ls())
library(bnlearn)
library(tidyverse)
library(plyr)
library(grid)
library(gridExtra)
library(psych)



library(ggplot2)
name<-'main'
data_out<-paste0('data/data_out/',name)
finalBN<-readRDS(paste0(data_out,'_finalBN.rds'))
fitted = readRDS(paste0(data_out,'_finalBN_fitted.rds'))
data = readRDS(paste0(data_out, '_discdata.rds'))
dataCor = data
dataCor[,grep("SUBJID|visitmiss|AUX|scode|DX|GENDER|ETHCAT|RACCAT|MARRY|APOE4|Amyloid", colnames(dataCor), value = TRUE)] = NULL
dataCor = lapply(dataCor, as.numeric)
dataCor = as.data.frame(dataCor)
corData =corr.test(dataCor, dataCor, use = "pairwise",method="spearman",adjust="holm", 
                   alpha=.05,ci=TRUE,minlength=5)

correlations = corData$r
corDataSub = as.data.frame(correlations)
corDataSub$from = rownames(corDataSub)

spearmanPadj = pivot_longer(corDataSub, cols=1:(ncol(corDataSub)-1), names_to = "to", values_to = "r")

pAdj = corData$p.adj
pAdjSub = as.data.frame(pAdj)
pAdjSub$from = rownames(pAdjSub)
padj = pivot_longer(pAdjSub, cols=1:(ncol(pAdjSub)-1), names_to = "to", values_to = "p.adj")


ciAdj = corData$ci.adj
ciAdjSub = as.data.frame(ciAdj)


spearmanPadj$p.Adj = padj$p.adj
spearmanPadj =cbind.data.frame(spearmanPadj, ciAdjSub)
spearmanPadj[,6] = NULL
spearmanPadj$r = round(spearmanPadj$r, 2)
spearmanPadj$lower.adj = round(spearmanPadj$lower.adj, 2)
spearmanPadj$upper.adj = round(spearmanPadj$upper.adj, 2)

spearmanPadj$from = gsub("zcode_", "", spearmanPadj$from)
spearmanPadj$from = gsub("SA_", "", spearmanPadj$from)
spearmanPadj$from = gsub("_VIS1", "", spearmanPadj$from)


spearmanPadj$to = gsub("zcode_", "", spearmanPadj$to)
spearmanPadj$to = gsub("SA_", "", spearmanPadj$to)
spearmanPadj$to = gsub("_VIS1", "", spearmanPadj$to)
boot_stren_50_percent <- read_delim("boot_stren_50_percent.csv",  ",", escape_double = FALSE, trim_ws = TRUE)

boot_stren_50_percent = merge(boot_stren_50_percent, spearmanPadj, all.x = TRUE)
boot_stren_50_percent$X1 = NULL

write.csv(boot_stren_50_percent, "network_Altoida_r_CI.csv")
network_Altoida_r_CI_upd = network_Altoida_r_CI[!network_Altoida_r_CI$from == "visitmiss",]
