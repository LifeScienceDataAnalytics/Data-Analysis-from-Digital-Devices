
setwd("~/Documents/Documents_IT/paper/ADNI_VAMBN_paper_final")

############################
############################ Dependencies and helper functions
############################
rm(list=ls())
library(bnlearn)
library(tidyverse)
library(plyr)
library(grid)
library(gridExtra)
library(Matching)
library(ggpubr)
library("pcalg")
library(psych)
library(ggplot2)
name<-'main'
data_out<-paste0('data/data_out/',name)
finalBN<-readRDS(paste0(data_out,'_finalBN.rds'))
fitted = readRDS(paste0(data_out,'_finalBN_fitted.rds'))
data = readRDS(paste0(data_out, '_discdata.rds'))
data[,grep("FAQ", colnames(data), value = TRUE)] = lapply(data[,grep("FAQ", colnames(data), value = TRUE)] , as.integer)
dataCor = data
dataCor[,grep("SUBJID|visitmiss|AUX|scode|DX|GENDER|ETHCAT|RACCAT|MARRY|APOE4|Amyloid", colnames(dataCor), value = TRUE)] = NULL
dataCor = lapply(dataCor, as.numeric)
dataCor = as.data.frame(dataCor)
corData =corr.test(dataCor, dataCor, use = "pairwise",method="spearman",adjust="holm", 
          alpha=.05,ci=TRUE,minlength=5)
cogDomains= grep("PerceptualMotorCoordination|Planning|ProspectiveMemory|SpatialMemory|CognitiveProcessingSpeed|ComplexAttention|EyeMovement|Flexibility|Inhibition|VisualPerception", colnames(data), value = TRUE)


correlations = corData$r
corDataSub = as.data.frame(correlations)
corDataSub = as.data.frame(corDataSub)
corDataSub$from = rownames(corDataSub)
library(tidyr)
spearmanRankPadj = pivot_longer(corDataSub, cols=1:(ncol(corDataSub)-1), names_to = "to", values_to = "r")


pAdj = corData$p.adj
pAdjSub = as.data.frame(pAdj)
pAdjSub$from = rownames(pAdjSub)
padj = pivot_longer(pAdjSub, cols=1:(ncol(pAdjSub)-1), names_to = "to", values_to = "p.adj")

ciAdj = corData$ci.adj
ciAdjSub = as.data.frame(ciAdj)

spearmanRankPadj$p.Adj = padj$p.adj
spearmanRankPadj =cbind.data.frame(spearmanRankPadj, ciAdjSub)
spearmanRankPadj$lower.adj = round(spearmanRankPadj$lower.adj,2)
spearmanRankPadj$upper.adj = round(spearmanRankPadj$upper.adj,2)
spearmanRankPadj$r = round(spearmanRankPadj$r,2)

spearmanRankPadjGoodCor = spearmanRankPadj[(spearmanRankPadj$r>0.20)|(spearmanRankPadj$r< -0.20),]
spearmanRankPadjGoodCor$lower.adj = round(spearmanRankPadjGoodCor$lower.adj,2)
spearmanRankPadjGoodCor$upper.adj = round(spearmanRankPadjGoodCor$upper.adj,2)
spearmanRankPadjGoodCor$r = round(spearmanRankPadjGoodCor$r,2)


######pearson correlations###
dataCorPs = dataCor[,!grepl("scode|DX|FAQ|MMSE|GENDER|ETHCAT|RACCAT|MARRY|APOE4|Amyloid", colnames(dataCor))]
corDataPearson=corr.test(dataCorPs, dataCorPs, use = "pairwise",method="pearson",adjust="holm", 
                         alpha=.05,ci=TRUE,minlength=5)

correlationsPs = corDataPearson$r
corDataSubPs = as.data.frame(correlationsPs)
corDataSubPs = as.data.frame(corDataSubPs)
corDataSubPs$from = rownames(corDataSubPs)

library(tidyr)
pearsonPadj = pivot_longer(corDataSubPs, cols=1:(ncol(corDataSubPs)-1), names_to = "to", values_to = "r")


pAdjPs = corDataPearson$p.adj
pAdjSubPs = as.data.frame(pAdjPs)
pAdjSubPs$from = rownames(pAdjSubPs)
padjPs = pivot_longer(pAdjSubPs, cols=1:(ncol(pAdjSubPs)-1), names_to = "to", values_to = "p.adj")

ciAdjPs = corDataPearson$ci.adj
ciAdjSubPs = as.data.frame(ciAdjPs)

pearsonPadj$p.Adj =  padjPs$p.adj
pearsonPadj =cbind.data.frame(pearsonPadj, ciAdjSubPs)
pearsonPadj$lower.adj = round(pearsonPadj$lower.adj,2)
pearsonPadj$upper.adj = round(pearsonPadj$upper.adj,2)
pearsonPadj$r = round(pearsonPadj$r,2)

pearsonPadjGoodCor = pearsonPadj[(pearsonPadj$r>0.20)|(pearsonPadj$r< -0.20),]
pearsonPadjGoodCor$lower.adj = round(pearsonPadjGoodCor$lower.adj,2)
pearsonPadjGoodCor$upper.adj = round(pearsonPadjGoodCor$upper.adj,2)
pearsonPadjGoodCor$r = round(pearsonPadjGoodCor$r,2)
########################################################################################################

res0 <-cor.test(data$zcode_imagingPET_VIS1, data$zcode_csf_VIS1,  method = "spearman")
res0
r0 = paste(spearmanRankPadjGoodCor[spearmanRankPadjGoodCor$from=="zcode_imagingPET_VIS1"&spearmanRankPadjGoodCor$to=="zcode_csf_VIS1","r"])
r0_tx = paste("R =", r0, sep = " ")
low_int0= spearmanRankPadjGoodCor[spearmanRankPadjGoodCor$from=="zcode_imagingPET_VIS1"&spearmanRankPadjGoodCor$to=="zcode_csf_VIS1","lower.adj"]
upper_int0= spearmanRankPadjGoodCor[spearmanRankPadjGoodCor$from=="zcode_imagingPET_VIS1"&spearmanRankPadjGoodCor$to=="zcode_csf_VIS1","upper.adj"]
conf_int_tx0 = paste("95% CI [", low_int0, ",", upper_int0, "]", sep = "")
p_adj0 = spearmanRankPadjGoodCor[spearmanRankPadjGoodCor$from=="zcode_imagingPET_VIS1"&spearmanRankPadjGoodCor$to=="zcode_csf_VIS1","p.Adj"]
p_adj_tx0 = paste("p.Adj =", "<0.0001", sep = " ")
label  = paste(r0_tx, conf_int_tx0, p_adj_tx0, sep = ", ")
png(filename="scatterPlots/scatterimagingcsf.png", width = 7, height = 7, units = 'in',res=300)
p0 = ggscatter(data, x = "zcode_imagingPET_VIS1", y = "zcode_csf_VIS1", 
               add = "reg.line", conf.int = TRUE, conf.int.level = 0.95, fullrange = TRUE,
               cor.coef = FALSE,  cor.method = "spearman",
               cor.coeff.args = list(method = "spearman", label.x.npc = "center", label.y.npc = "top"),
               xlab = "imagingPET_VIS1", ylab = "csf_VIS1", font.x = 16,
               font.y = 16) + 
  annotate("text", x = 1, y = 3, label = label,  size=5) 
p0
dev.off()

res0 <-cor.test(data$SA_PHS_VIS1, data$zcode_imagingPET_VIS1,  method = "spearman")
res0
r0 = paste(spearmanRankPadjGoodCor[spearmanRankPadjGoodCor$from=="SA_PHS_VIS1"&spearmanRankPadjGoodCor$to=="zcode_imagingPET_VIS1","r"])
r0_tx = paste("R =", r0, sep = " ")
low_int0= spearmanRankPadjGoodCor[spearmanRankPadjGoodCor$from=="SA_PHS_VIS1"&spearmanRankPadjGoodCor$to=="zcode_imagingPET_VIS1","lower.adj"]
upper_int0= spearmanRankPadjGoodCor[spearmanRankPadjGoodCor$from=="SA_PHS_VIS1"&spearmanRankPadjGoodCor$to=="zcode_imagingPET_VIS1","upper.adj"]
conf_int_tx0 = paste("95% CI [", low_int0, ",", upper_int0, "]", sep = "")
p_adj0 = spearmanRankPadjGoodCor[spearmanRankPadjGoodCor$from=="SA_PHS_VIS1"&spearmanRankPadjGoodCor$to=="zcode_imagingPET_VIS1","p.Adj"]
p_adj_tx0 = paste("p.Adj =", "<0.0001", sep = " ")
label  = paste(r0_tx, conf_int_tx0, p_adj_tx0, sep = ", ")
png(filename="scatterPlots/scatterPHSimaging.png", width = 7, height = 7, units = 'in',res=300)
p0 = ggscatter(data, x = "SA_PHS_VIS1", y = "zcode_imagingPET_VIS1", 
               add = "reg.line", conf.int = TRUE, conf.int.level = 0.95, fullrange = TRUE,
               cor.coef = FALSE,  cor.method = "spearman",
               cor.coeff.args = list(method = "spearman", label.x.npc = "center", label.y.npc = "top"),
               xlab = "PHS_VIS1", ylab = "imagingPET_VIS1", font.x = 16,
               font.y = 16) + 
  annotate("text", x = 1, y =3, label = label,  size=5) 
p0
dev.off()

res0 <-cor.test(data$zcode_csf_VIS1, data$SA_Planning_VIS1,  method = "spearman")
res0
r0 = paste(spearmanRankPadj[spearmanRankPadj$from=="zcode_csf_VIS1"&spearmanRankPadj$to=="SA_Planning_VIS1","r"])
r0_tx = paste("R =", r0, sep = " ")
low_int0= spearmanRankPadj[spearmanRankPadj$from=="zcode_csf_VIS1"&spearmanRankPadj$to=="SA_Planning_VIS1","lower.adj"]
upper_int0= spearmanRankPadj[spearmanRankPadj$from=="zcode_csf_VIS1"&spearmanRankPadj$to=="SA_Planning_VIS1","upper.adj"]
conf_int_tx0 = paste("95% CI [", low_int0, ",", upper_int0, "]", sep = "")
p_adj0 = spearmanRankPadj[spearmanRankPadj$from=="zcode_csf_VIS1"&spearmanRankPadj$to=="SA_Planning_VIS1","p.Adj"]
p_adj_tx0 = paste("p.Adj =", "<0.0001", sep = " ")
label  = paste(r0_tx, conf_int_tx0, p_adj_tx0, sep = ", ")
png(filename="scatterPlots/scattercsfPlanning.png", width = 7, height = 7, units = 'in',res=300)
p01 = ggscatter(data, x = "zcode_csf_VIS1", y = "SA_Planning_VIS1", 
               add = "reg.line", conf.int = TRUE, conf.int.level = 0.95, fullrange = TRUE,
               cor.coef = FALSE,  cor.method = "spearman",
               cor.coeff.args = list(method = "spearman", label.x.npc = "center", label.y.npc = "top"),
               xlab = "csf_VIS1", ylab = "Planning_VIS1", font.x = 16,
               font.y = 16) + 
  annotate("text", x = -0.5, y = 67, label = label,  size=5) 
p01
dev.off()


##not significant
res0 <-cor.test(data$zcode_csf_VIS1, data$SA_ComplexAttention_VIS1,  method = "spearman")
res0
r0 = paste(spearmanRankPadj[spearmanRankPadj$from=="zcode_csf_VIS1"&spearmanRankPadj$to=="SA_ComplexAttention_VIS1","r"])
r0_tx = paste("R =", r0, sep = " ")
low_int0= spearmanRankPadj[spearmanRankPadj$from=="zcode_csf_VIS1"&spearmanRankPadj$to=="SA_ComplexAttention_VIS1","lower.adj"]
upper_int0= spearmanRankPadj[spearmanRankPadj$from=="zcode_csf_VIS1"&spearmanRankPadj$to=="SA_ComplexAttention_VIS1","upper.adj"]
conf_int_tx0 = paste("95% CI [", low_int0, ",", upper_int0, "]", sep = "")
p_adj0 = spearmanRankPadj[spearmanRankPadj$from=="zcode_csf_VIS1"&spearmanRankPadj$to=="SA_ComplexAttention_VIS1","p.Adj"]
p_adj_tx0 = paste("p.Adj =", "<2.2e-16", sep = " ")
label  = paste(r0_tx, conf_int_tx0, p_adj_tx0, sep = ", ")
png(filename="scatterPlots/scattercsfComplexAttention.png", width = 7, height = 7, units = 'in',res=300)
p02 = ggscatter(data, x = "zcode_csf_VIS1", y = "SA_ComplexAttention_VIS1", 
               add = "reg.line", conf.int = TRUE, conf.int.level = 0.95, fullrange = TRUE,
               cor.coef = FALSE,  cor.method = "spearman",
               cor.coeff.args = list(method = "spearman", label.x.npc = "center", label.y.npc = "top"),
               xlab = "csf_VIS1", ylab = "ComplexAttention_VIS1", font.x = 16,
               font.y = 16) + 
  annotate("text", x = -0.5, y = 53, label = label,  size=5) 
p02
dev.off()

res0 <-cor.test(data$SA_PHS_VIS1, data$SA_zcode_ARObjectPlacement_VIS1,  method = "spearman")
res0
r0 = paste(spearmanRankPadjGoodCor[spearmanRankPadjGoodCor$from=="SA_PHS_VIS1"&spearmanRankPadjGoodCor$to=="SA_zcode_ARObjectPlacement_VIS1","r"])
r0_tx = paste("R =", r0, sep = " ")
low_int0= spearmanRankPadjGoodCor[spearmanRankPadjGoodCor$from=="SA_PHS_VIS1"&spearmanRankPadjGoodCor$to=="SA_zcode_ARObjectPlacement_VIS1","lower.adj"]
upper_int0= spearmanRankPadjGoodCor[spearmanRankPadjGoodCor$from=="SA_PHS_VIS1"&spearmanRankPadjGoodCor$to=="SA_zcode_ARObjectPlacement_VIS1","upper.adj"]
conf_int_tx0 = paste("95% CI [", low_int0, ",", upper_int0, "]", sep = "")
p_adj0 = spearmanRankPadjGoodCor[spearmanRankPadjGoodCor$from=="SA_PHS_VIS1"&spearmanRankPadjGoodCor$to=="SA_zcode_ARObjectPlacement_VIS1","p.Adj"]
p_adj_tx0 = paste("p.Adj =", "<2.2e-16", sep = " ")
label  = paste(r0_tx, conf_int_tx0, p_adj_tx0, sep = ", ")
png(filename="scatterPlots/scatterPHSObjectPlacement.png", width = 7, height = 7, units = 'in',res=300)
p03 = ggscatter(data, x = "SA_PHS_VIS1", y = "SA_zcode_ARObjectPlacement_VIS1", 
               add = "reg.line", conf.int = TRUE, conf.int.level = 0.95, fullrange = TRUE,
               cor.coef = FALSE,  cor.method = "spearman",
               cor.coeff.args = list(method = "spearman", label.x.npc = "center", label.y.npc = "top"),
               xlab = "SA_PHS_VIS1", ylab = "ARObjectPlacement_VIS1", font.x = 16,
               font.y = 16) + 
  annotate("text", x = -0.5, y = 53, label = label,  size=5) 
p03
dev.off()


res0 <-cor.test(data$SA_PHS_VIS1, data$SA_zcode_ARObjectFinding_VIS1,  method = "spearman")
res0
r0 = paste(spearmanRankPadj[spearmanRankPadj$from=="SA_PHS_VIS1"&spearmanRankPadj$to=="SA_zcode_ARObjectFinding_VIS1","r"])
r0_tx = paste("R =", r0, sep = " ")
low_int0= spearmanRankPadj[spearmanRankPadj$from=="SA_PHS_VIS1"&spearmanRankPadj$to=="SA_zcode_ARObjectFinding_VIS1","lower.adj"]
upper_int0= spearmanRankPadj[spearmanRankPadj$from=="SA_PHS_VIS1"&spearmanRankPadj$to=="SA_zcode_ARObjectFinding_VIS1","upper.adj"]
conf_int_tx0 = paste("95% CI [", low_int0, ",", upper_int0, "]", sep = "")
p_adj0 = spearmanRankPadj[spearmanRankPadj$from=="SA_PHS_VIS1"&spearmanRankPadj$to=="SA_zcode_ARObjectFinding_VIS1","p.Adj"]
p_adj_tx0 = paste("p.Adj =", "<0.0001", sep = " ")
label  = paste(r0_tx, conf_int_tx0, p_adj_tx0, sep = ", ")
png(filename="scatterPlots/scatterPHSObjectFinding.png", width = 7, height = 7, units = 'in',res=300)
p04 = ggscatter(data, x = "SA_PHS_VIS1", y = "SA_zcode_ARObjectFinding_VIS1", 
               add = "reg.line", conf.int = TRUE, conf.int.level = 0.95, fullrange = TRUE,
               cor.coef = FALSE,  cor.method = "spearman",
               cor.coeff.args = list(method = "spearman", label.x.npc = "center", label.y.npc = "top"),
               xlab = "SA_PHS_VIS1", ylab = "ARObjectFinding_VIS1", font.x = 16,
               font.y = 16) + 
  annotate("text", x = 1, y = 5, label = label,  size=5) 
p04
dev.off()

##MMSE Orientation and AR Screen Button Presses
library("ggpubr")
res1 <-cor.test(data$zcode_MMSE_Orientation_VIS1, data$SA_zcode_ARScreenButtonPresses_VIS1,  method = "spearman")
res1
r = paste(spearmanRankPadjGoodCor[spearmanRankPadjGoodCor$from=="zcode_MMSE_Orientation_VIS1"&spearmanRankPadjGoodCor$to=="SA_zcode_ARScreenButtonPresses_VIS1","r"])
r_tx = paste("R =", r, sep = " ")
low_int= spearmanRankPadjGoodCor[spearmanRankPadjGoodCor$from=="zcode_MMSE_Orientation_VIS1"&spearmanRankPadjGoodCor$to=="SA_zcode_ARScreenButtonPresses_VIS1","lower.adj"]
upper_int= spearmanRankPadjGoodCor[spearmanRankPadjGoodCor$from=="zcode_MMSE_Orientation_VIS1"&spearmanRankPadjGoodCor$to=="SA_zcode_ARScreenButtonPresses_VIS1","upper.adj"]
conf_int_tx = paste("95% CI [", low_int, ",", upper_int, "]", sep = "")
p_adj = spearmanRankPadjGoodCor[spearmanRankPadjGoodCor$from=="zcode_MMSE_Orientation_VIS1"&spearmanRankPadjGoodCor$to=="SA_zcode_ARScreenButtonPresses_VIS1","p.Adj"]
p_adj_tx = paste("p.Adj =", "<0.0001", sep = " ")
label  = paste(r_tx, conf_int_tx, p_adj_tx, sep = ", ")
png(filename="scatterPlots/scatterOriScreenButton.png", width = 7, height = 7, units = 'in',res=300)
p1 = ggscatter(data, x = "zcode_MMSE_Orientation_VIS1", y = "SA_zcode_ARScreenButtonPresses_VIS1", 
               add = "reg.line", conf.int = TRUE, conf.int.level = 0.95, fullrange = TRUE,
               cor.coef = FALSE,  cor.method = "spearman",
               cor.coeff.args = list(method = "spearman", label.x.npc = "center", label.y.npc = "top"),
               xlab = "MMSE_Orientation_VIS1", ylab = "ARScreenButtonPresses_VIS1", font.x = 16,
               font.y = 16) + 
  annotate("text", x = 4, y = 3.5, label = label,  size=5) 
p1
dev.off()

res2 <-cor.test(data$zcode_MMSE_Orientation_VIS6, data$SA_zcode_ARScreenButtonPresses_VIS6,  method = "spearman")
res2
r2 = paste(spearmanRankPadjGoodCor[spearmanRankPadjGoodCor$from=="zcode_MMSE_Orientation_VIS6"&spearmanRankPadjGoodCor$to=="SA_zcode_ARScreenButtonPresses_VIS6","r"])
r2_tx = paste("R =", r2, sep = " ")
low_int2= spearmanRankPadjGoodCor[spearmanRankPadjGoodCor$from=="zcode_MMSE_Orientation_VIS6"&spearmanRankPadjGoodCor$to=="SA_zcode_ARScreenButtonPresses_VIS6","lower.adj"]
upper_int2= spearmanRankPadjGoodCor[spearmanRankPadjGoodCor$from=="zcode_MMSE_Orientation_VIS6"&spearmanRankPadjGoodCor$to=="SA_zcode_ARScreenButtonPresses_VIS6","upper.adj"]
conf_int_tx2 = paste("95% CI [", low_int2,",", upper_int2, "]", sep = "")
p_adj2 = spearmanRankPadjGoodCor[spearmanRankPadjGoodCor$from=="zcode_MMSE_Orientation_VIS6"&spearmanRankPadjGoodCor$to=="SA_zcode_ARScreenButtonPresses_VIS6","p.Adj"]
p_adj_tx2 = paste("p.Adj =", "<0.0001", sep = " ")
label2  = paste(r2_tx, conf_int_tx2, p_adj_tx2, sep = ", ")
png(filename="scatterPlots/scatterOriScreenButton_vis6.png", width = 7, height = 7, units = 'in',res=300)
p2 = ggscatter(data, x = "zcode_MMSE_Orientation_VIS6", y = "SA_zcode_ARScreenButtonPresses_VIS6", 
               add = "reg.line", conf.int = TRUE, conf.int.level = 0.95, fullrange = TRUE,
               cor.coef = FALSE,  cor.method = "spearman",
               cor.coeff.args = list(method = "spearman", label.x.npc = "center", label.y.npc = "top"),
               xlab = "MMSE_Orientation_VIS6", ylab = "ARScreenButtonPresses_VIS6", font.x = 16,
               font.y = 16) + 
  annotate("text", x = 5, y = 3.5, label = label2,  size=5) 
p2
dev.off()


res3 <-cor.test(data$zcode_MMSE_Orientation_VIS12, data$SA_zcode_ARScreenButtonPresses_VIS12,  method = "spearman")
res3
r3 = paste(spearmanRankPadjGoodCor[spearmanRankPadjGoodCor$from=="zcode_MMSE_Orientation_VIS12"&spearmanRankPadjGoodCor$to=="SA_zcode_ARScreenButtonPresses_VIS12","r"])
r3_tx = paste("R =", r3, sep = " ")
low_int3= spearmanRankPadjGoodCor[spearmanRankPadjGoodCor$from=="zcode_MMSE_Orientation_VIS12"&spearmanRankPadjGoodCor$to=="SA_zcode_ARScreenButtonPresses_VIS12","lower.adj"]
upper_int3= spearmanRankPadjGoodCor[spearmanRankPadjGoodCor$from=="zcode_MMSE_Orientation_VIS12"&spearmanRankPadjGoodCor$to=="SA_zcode_ARScreenButtonPresses_VIS12","upper.adj"]
conf_int_tx3 = paste("95% CI [", low_int3,",", upper_int3, "]", sep = "")
p_adj3 = spearmanRankPadjGoodCor[spearmanRankPadjGoodCor$from=="zcode_MMSE_Orientation_VIS6"&spearmanRankPadjGoodCor$to=="SA_zcode_ARScreenButtonPresses_VIS6","p.Adj"]
p_adj_tx3 = paste("p.Adj =", "<0.0001", sep = " ")
label3  = paste(r3_tx, conf_int_tx3, p_adj_tx3, sep = ", ")
png(filename="scatterPlots/scatterOriScreenButton_vis12.png", width = 7, height = 7, units = 'in',res=300)
p3 = ggscatter(data, x = "zcode_MMSE_Orientation_VIS12", y = "SA_zcode_ARScreenButtonPresses_VIS12", 
               add = "reg.line", conf.int = TRUE, conf.int.level = 0.95, fullrange = TRUE,
               cor.coef = FALSE,  cor.method = "spearman",
               cor.coeff.args = list(method = "spearman", label.x.npc = "center", label.y.npc = "top"),
               xlab = "MMSE_Orientation_VIS12", ylab = "ARScreenButtonPresses_VIS12", font.x = 16,
               font.y = 16) + 
  annotate("text", x = 3, y = 3.5, label = label3,  size=5) 
p3
dev.off()


##MMSE Language and zcode_ARPlaceAndFindTelemetryVariance_VIS1##
res4 <-cor.test(data$zcode_MMSE_Language_VIS1, data$SA_zcode_ARPlaceAndFindTelemetryVariance_VIS1,  method = "spearman")
res4
r4 = paste(spearmanRankPadjGoodCor[spearmanRankPadjGoodCor$from=="zcode_MMSE_Language_VIS1"&spearmanRankPadjGoodCor$to=="SA_zcode_ARPlaceAndFindTelemetryVariance_VIS1","r"])
r4_tx = paste("R =", r4, sep = " ")
low_int4= spearmanRankPadjGoodCor[spearmanRankPadjGoodCor$from=="zcode_MMSE_Language_VIS1"&spearmanRankPadjGoodCor$to=="SA_zcode_ARPlaceAndFindTelemetryVariance_VIS1","lower.adj"]
upper_int4= spearmanRankPadjGoodCor[spearmanRankPadjGoodCor$from=="zcode_MMSE_Language_VIS1"&spearmanRankPadjGoodCor$to=="SA_zcode_ARPlaceAndFindTelemetryVariance_VIS1","upper.adj"]
conf_int_tx4 = paste("95% CI [", low_int4,",", upper_int4, "]", sep = "")
p_adj4 = spearmanRankPadjGoodCor[spearmanRankPadjGoodCor$from=="zcode_MMSE_Language_VIS1"&spearmanRankPadjGoodCor$to=="SA_zcode_ARPlaceAndFindTelemetryVariance_VIS1","p.Adj"]
p_adj_tx4 = paste("p.Adj =", "<0.0001", sep = " ")
label4  = paste(r4_tx, conf_int_tx4, p_adj_tx4, sep = ", ")
png(filename="scatterPlots/scatterLanTelVar_vis1.png", width = 7, height = 7, units = 'in',res=300)
p4 = ggscatter(data, x = "zcode_MMSE_Language_VIS1", y = "SA_zcode_ARPlaceAndFindTelemetryVariance_VIS1", 
               add = "reg.line", conf.int = TRUE, conf.int.level = 0.95, fullrange = TRUE,
               cor.coef = FALSE,  cor.method = "spearman",
               cor.coeff.args = list(method = "spearman", label.x.npc = "center", label.y.npc = "top"),
               xlab = "MMSE_Language_VIS1", ylab = "ARPlaceAndFindTelemetryVariance_VIS1", font.x = 16,
               font.y = 16) + annotate("text", x = -5, y = 7, label = label4,  size=5) 
p4
dev.off()


##month 6
res5 <-cor.test(data$zcode_MMSE_Language_VIS6, data$SA_zcode_ARPlaceAndFindTelemetryVariance_VIS6,  method = "spearman")
res5
r5 = paste(spearmanRankPadjGoodCor[spearmanRankPadjGoodCor$from=="zcode_MMSE_Language_VIS6"&spearmanRankPadjGoodCor$to=="SA_zcode_ARPlaceAndFindTelemetryVariance_VIS6","r"])
r5_tx = paste("R =", r5, sep = " ")
low_int5= spearmanRankPadjGoodCor[spearmanRankPadjGoodCor$from=="zcode_MMSE_Language_VIS6"&spearmanRankPadjGoodCor$to=="SA_zcode_ARPlaceAndFindTelemetryVariance_VIS6","lower.adj"]
upper_int5= spearmanRankPadjGoodCor[spearmanRankPadjGoodCor$from=="zcode_MMSE_Language_VIS6"&spearmanRankPadjGoodCor$to=="SA_zcode_ARPlaceAndFindTelemetryVariance_VIS6","upper.adj"]
conf_int_tx5 = paste("95% CI [", low_int5, ",",upper_int5, "]", sep = "")
p_adj5 = spearmanRankPadjGoodCor[spearmanRankPadjGoodCor$from=="zcode_MMSE_Language_VIS6"&spearmanRankPadjGoodCor$to=="SA_zcode_ARPlaceAndFindTelemetryVariance_VIS6","p.Adj"]
p_adj_tx5 = paste("p.Adj =", "<0.0001", sep = " ")
label5  = paste(r5_tx, conf_int_tx5, p_adj_tx5, sep = ", ")
png(filename="scatterPlots/scatterLanTelVar_vis6.png", width = 7, height = 7, units = 'in',res=300)
p5 = ggscatter(data, x = "zcode_MMSE_Language_VIS6", y = "SA_zcode_ARPlaceAndFindTelemetryVariance_VIS6", 
               add = "reg.line", conf.int = TRUE, conf.int.level = 0.95, fullrange = TRUE,
               cor.coef = FALSE,  cor.method = "spearman",
               cor.coeff.args = list(method = "spearman", label.x.npc = "center", label.y.npc = "top"),
               xlab = "MMSE_Language_VIS6", ylab = "ARPlaceAndFindTelemetryVariance_VIS6", font.x = 16,
               font.y = 16) + annotate("text", x = 4, y = 5, label = label5,  size=5) 
p5
dev.off()

##month 12
res6 <-cor.test(data$zcode_MMSE_Language_VIS12, data$SA_zcode_ARPlaceAndFindTelemetryVariance_VIS12,  method = "spearman")
res6
r6 = paste(spearmanRankPadjGoodCor[spearmanRankPadjGoodCor$from=="zcode_MMSE_Language_VIS12"&spearmanRankPadjGoodCor$to=="SA_zcode_ARPlaceAndFindTelemetryVariance_VIS12","r"])
r6_tx = paste("R =", r6, sep = " ")
low_int6= spearmanRankPadjGoodCor[spearmanRankPadjGoodCor$from=="zcode_MMSE_Language_VIS12"&spearmanRankPadjGoodCor$to=="SA_zcode_ARPlaceAndFindTelemetryVariance_VIS12","lower.adj"]
upper_int6= spearmanRankPadjGoodCor[spearmanRankPadjGoodCor$from=="zcode_MMSE_Language_VIS12"&spearmanRankPadjGoodCor$to=="SA_zcode_ARPlaceAndFindTelemetryVariance_VIS12","upper.adj"]
conf_int_tx6 = paste("95% CI [", low_int6, ",",upper_int6, "]", sep = "")
p_adj6 = spearmanRankPadjGoodCor[spearmanRankPadjGoodCor$from=="zcode_MMSE_Language_VIS12"&spearmanRankPadjGoodCor$to=="SA_zcode_ARPlaceAndFindTelemetryVariance_VIS12","p.Adj"]
p_adj_tx6 = paste("p.Adj =", "<0.0001", sep = " ")
label6  = paste(r6_tx, conf_int_tx6, p_adj_tx6, sep = ", ")
png(filename="scatterPlots/scatterLanTelVar_vis12.png", width = 7, height = 7, units = 'in',res=300)
p6 = ggscatter(data, x = "zcode_MMSE_Language_VIS12", y = "SA_zcode_ARPlaceAndFindTelemetryVariance_VIS12", 
               add = "reg.line", conf.int = TRUE, conf.int.level = 0.95, fullrange = TRUE,
               cor.coef = FALSE,  cor.method = "spearman",
               cor.coeff.args = list(method = "spearman", label.x.npc = "center", label.y.npc = "top"),
               xlab = "MMSE_Language_VIS12", ylab = "ARPlaceAndFindTelemetryVariance_VIS12", font.x = 16,
               font.y = 16) + annotate("text", x = 3, y = 6, label = label6,  size=5) 
p6
dev.off()

##MMSE Language and zcode_ARObjectPlacement##
res7 <-cor.test(data$zcode_MMSE_Language_VIS1, data$SA_zcode_ARObjectPlacement_VIS1,  method = "spearman")
res7
r7 = paste(spearmanRankPadjGoodCor[spearmanRankPadjGoodCor$from=="zcode_MMSE_Language_VIS1"&spearmanRankPadjGoodCor$to=="SA_zcode_ARObjectPlacement_VIS1","r"])
r7_tx = paste("R =", r7, sep = " ")
low_int7= spearmanRankPadjGoodCor[spearmanRankPadjGoodCor$from=="zcode_MMSE_Language_VIS1"&spearmanRankPadjGoodCor$to=="SA_zcode_ARObjectPlacement_VIS1","lower.adj"]
upper_int7= spearmanRankPadjGoodCor[spearmanRankPadjGoodCor$from=="zcode_MMSE_Language_VIS1"&spearmanRankPadjGoodCor$to=="SA_zcode_ARObjectPlacement_VIS1","upper.adj"]
conf_int_tx7 = paste("95% CI [", low_int7,",", upper_int7, "]", sep = "")
p_adj7 = spearmanRankPadjGoodCor[spearmanRankPadjGoodCor$from=="zcode_MMSE_Language_VIS1"&spearmanRankPadjGoodCor$to=="SA_zcode_ARObjectPlacement_VIS1","p.Adj"]
p_adj_tx7 = paste("p.Adj =", "<0.0001", sep = " ")
label7  = paste(r7_tx, conf_int_tx7, p_adj_tx7, sep = ", ")
png(filename="scatterPlots/scatterLanObjPlac_vis1.png", width = 7, height = 7, units = 'in',res=300)
p7 = ggscatter(data, x = "zcode_MMSE_Language_VIS1", y = "SA_zcode_ARObjectPlacement_VIS1", 
               add = "reg.line", conf.int = TRUE, conf.int.level = 0.95, fullrange = TRUE,
               cor.coef = FALSE,  cor.method = "spearman",
               cor.coeff.args = list(method = "spearman", label.x.npc = "center", label.y.npc = "top"),
               xlab = "MMSE_Language_VIS1", ylab = "ARObjectPlacement_VIS1", font.x = 16,
               font.y = 16) + annotate("text", x = -5, y = 4, label = label7,  size=5) 
p7
dev.off()

##month 6
res8 <-cor.test(data$zcode_MMSE_Language_VIS6, data$SA_zcode_ARObjectPlacement_VIS6,  method = "spearman")
res8
r8 = paste(spearmanRankPadj[spearmanRankPadj$from=="zcode_MMSE_Language_VIS6"&spearmanRankPadj$to=="SA_zcode_ARObjectPlacement_VIS6","r"])
r8_tx = paste("R =", r8, sep = " ")
low_int8= spearmanRankPadj[spearmanRankPadj$from=="zcode_MMSE_Language_VIS6"&spearmanRankPadj$to=="SA_zcode_ARObjectPlacement_VIS6","lower.adj"]
upper_int8= spearmanRankPadj[spearmanRankPadj$from=="zcode_MMSE_Language_VIS6"&spearmanRankPadj$to=="SA_zcode_ARObjectPlacement_VIS6","upper.adj"]
conf_int_tx8 = paste("95% CI [", low_int8,",", upper_int8, "]", sep = "")
p_adj8 = spearmanRankPadj[spearmanRankPadj$from=="zcode_MMSE_Language_VIS6"&spearmanRankPadj$to=="SA_zcode_ARObjectPlacement_VIS6","p.Adj"]
p_adj_tx8 = paste("p.Adj =", "<0.0001", sep = " ")
label8  = paste(r8_tx, conf_int_tx8, p_adj_tx8, sep = ", ")
png(filename="scatterPlots/scatterLanObjPlac_vis6.png", width = 7, height = 7, units = 'in',res=300)
p8 = ggscatter(data, x = "zcode_MMSE_Language_VIS6", y = "SA_zcode_ARObjectPlacement_VIS6", 
               add = "reg.line", conf.int = TRUE, conf.int.level = 0.95, fullrange = TRUE,
               cor.coef = FALSE,  cor.method = "spearman",
               cor.coeff.args = list(method = "spearman", label.x.npc = "center", label.y.npc = "top"),
               xlab = "MMSE_Language_VIS6", ylab = "ARObjectPlacement_VIS6", font.x = 16,
               font.y = 16) + annotate("text", x = 4, y = 4, label = label8,  size=5) 
p8
dev.off()

##month 12
res9<-cor.test(data$zcode_MMSE_Language_VIS12, data$SA_zcode_ARObjectPlacement_VIS12,  method = "spearman")
res9
r9 = paste(spearmanRankPadj[spearmanRankPadj$from=="zcode_MMSE_Language_VIS12"&spearmanRankPadj$to=="SA_zcode_ARObjectPlacement_VIS12","r"])
r9_tx = paste("R =", r9, sep = " ")
low_int9= spearmanRankPadj[spearmanRankPadj$from=="zcode_MMSE_Language_VIS12"&spearmanRankPadj$to=="SA_zcode_ARObjectPlacement_VIS12","lower.adj"]
upper_int9= spearmanRankPadj[spearmanRankPadj$from=="zcode_MMSE_Language_VIS12"&spearmanRankPadj$to=="SA_zcode_ARObjectPlacement_VIS12","upper.adj"]
conf_int_tx9 = paste("95% CI [", low_int9,",", upper_int9, "]", sep = "")
p_adj9 = spearmanRankPadj[spearmanRankPadj$from=="zcode_MMSE_Language_VIS12"&spearmanRankPadj$to=="SA_zcode_ARObjectPlacement_VIS12","p.Adj"]
p_adj_tx9 = paste("p.Adj =", "<0.0001", sep = " ")
label9  = paste(r9_tx, conf_int_tx9, p_adj_tx9, sep = ", ")
png(filename="scatterPlots/scatterLanObjPlac_vis12.png", width = 7, height = 7, units = 'in',res=300)
p9 = ggscatter(data, x = "zcode_MMSE_Language_VIS12", y = "SA_zcode_ARObjectPlacement_VIS12", 
               add = "reg.line", conf.int = TRUE, conf.int.level = 0.95, fullrange = TRUE,
               cor.coef = FALSE,  cor.method = "spearman",
               cor.coeff.args = list(method = "spearman", label.x.npc = "center", label.y.npc = "top"),
               xlab = "MMSE_Language_VIS12", ylab = "ARObjectPlacement_VIS12", font.x = 16,
               font.y = 16) + annotate("text", x = 3.5, y = 4, label = label9,  size=5) 
p9
dev.off()

##MMSE MEMORY RECALL and Flexibility at vis1##
res10 <-cor.test(data$zcode_MMSE_Memory_Recall_VIS1, data$SA_Flexibility_VIS1,  method = "spearman")
res10
r10 = paste(spearmanRankPadjGoodCor[spearmanRankPadjGoodCor$from=="zcode_MMSE_Memory_Recall_VIS1"&spearmanRankPadjGoodCor$to=="SA_Flexibility_VIS1","r"])
r10_tx = paste("R =", r10, sep = " ")
low_int10= spearmanRankPadjGoodCor[spearmanRankPadjGoodCor$from=="zcode_MMSE_Memory_Recall_VIS1"&spearmanRankPadjGoodCor$to=="SA_Flexibility_VIS1","lower.adj"]
upper_int10= spearmanRankPadjGoodCor[spearmanRankPadjGoodCor$from=="zcode_MMSE_Memory_Recall_VIS1"&spearmanRankPadjGoodCor$to=="SA_Flexibility_VIS1","upper.adj"]
conf_int_tx10 = paste("95% CI [", low_int10,",", upper_int10, "]", sep = "")
p_adj10 = spearmanRankPadjGoodCor[spearmanRankPadjGoodCor$from=="zcode_MMSE_Memory_Recall_VIS1"&spearmanRankPadjGoodCor$to=="SA_Flexibility_VIS1","p.Adj"]
p_adj_tx10 = paste("p.Adj =", "<0.0001", sep = " ")
label10  = paste(r10_tx, conf_int_tx10, p_adj_tx10, sep = ", ")
png(filename="scatterPlots/scatterMemRecFlex_vis1.png", width = 7, height = 7, units = 'in',res=300)
p10 = ggscatter(data, x = "zcode_MMSE_Memory_Recall_VIS1", y = "SA_Flexibility_VIS1", 
               add = "reg.line", conf.int = TRUE, conf.int.level = 0.95, fullrange = TRUE,
               cor.coef = FALSE,  cor.method = "spearman",
               cor.coeff.args = list(method = "spearman", label.x.npc = "center", label.y.npc = "top"),
               xlab = "MMSE_Recall_VIS1", ylab = "Flexibility_VIS1", font.x = 16,
               font.y = 16) + annotate("text", x = -7.5, y = 84, label = label10,  size=5) 
p10
dev.off()

##MMSE MEMORY RECALL and Flexibility at vis6##
res11 <-cor.test(data$zcode_MMSE_Memory_Recall_VIS6, data$SA_Flexibility_VIS6,  method = "spearman")
res11
r11 = paste(spearmanRankPadjGoodCor[spearmanRankPadjGoodCor$from=="zcode_MMSE_Memory_Recall_VIS6"&spearmanRankPadjGoodCor$to=="SA_Flexibility_VIS6","r"])
r11_tx = paste("R =", r11, sep = " ")
low_int11= spearmanRankPadjGoodCor[spearmanRankPadjGoodCor$from=="zcode_MMSE_Memory_Recall_VIS6"&spearmanRankPadjGoodCor$to=="SA_Flexibility_VIS6","lower.adj"]
upper_int11= spearmanRankPadjGoodCor[spearmanRankPadjGoodCor$from=="zcode_MMSE_Memory_Recall_VIS6"&spearmanRankPadjGoodCor$to=="SA_Flexibility_VIS6","upper.adj"]
conf_int_tx11 = paste("95% CI [", low_int11,",", upper_int11, "]", sep = "")
p_adj11 = spearmanRankPadjGoodCor[spearmanRankPadjGoodCor$from=="zcode_MMSE_Memory_Recall_VIS6"&spearmanRankPadjGoodCor$to=="SA_Flexibility_VIS6","p.Adj"]
p_adj_tx11 = paste("p.Adj =", "<0.0001", sep = " ")
label11  = paste(r11_tx, conf_int_tx11, p_adj_tx11, sep = ", ")
png(filename="scatterPlots/scatterMemRecFlex_vis6.png", width = 7, height = 7, units = 'in',res=300)
p11 = ggscatter(data, x = "zcode_MMSE_Memory_Recall_VIS6", y = "SA_Flexibility_VIS6", 
                add = "reg.line", conf.int = TRUE, conf.int.level = 0.95, fullrange = TRUE,
                cor.coef = FALSE,  cor.method = "spearman",
                cor.coeff.args = list(method = "spearman", label.x.npc = "center", label.y.npc = "top"),
                xlab = "MMSE_Recall_VIS6", ylab = "Flexibility_VIS6", font.x = 16,
                font.y = 16) + annotate("text", x = -5, y = 80, label = label11,  size=5) 
p11
dev.off()

##MMSE MEMORY RECALL and Flexibility at vis12##
res12 <-cor.test(data$zcode_MMSE_Memory_Recall_VIS12, data$SA_Flexibility_VIS12,  method = "spearman")
res12
r12 = paste(spearmanRankPadjGoodCor[spearmanRankPadjGoodCor$from=="zcode_MMSE_Memory_Recall_VIS12"&spearmanRankPadjGoodCor$to=="SA_Flexibility_VIS12","r"])
r12_tx = paste("R =", r12, sep = " ")
low_int12= spearmanRankPadjGoodCor[spearmanRankPadjGoodCor$from=="zcode_MMSE_Memory_Recall_VIS12"&spearmanRankPadjGoodCor$to=="SA_Flexibility_VIS12","lower.adj"]
upper_int12= spearmanRankPadjGoodCor[spearmanRankPadjGoodCor$from=="zcode_MMSE_Memory_Recall_VIS12"&spearmanRankPadjGoodCor$to=="SA_Flexibility_VIS12","upper.adj"]
conf_int_tx12 = paste("95% CI [", low_int12,",", upper_int12, "]", sep = "")
p_adj12 = spearmanRankPadjGoodCor[spearmanRankPadjGoodCor$from=="zcode_MMSE_Memory_Recall_VIS12"&spearmanRankPadjGoodCor$to=="SA_Flexibility_VIS12","p.Adj"]
p_adj_tx12 = paste("p.Adj =", "<0.0001", sep = " ")
label12  = paste(r12_tx, conf_int_tx12, p_adj_tx12, sep = ", ")
png(filename="scatterPlots/scatterMemRecFlex_vis12.png", width = 7, height = 7, units = 'in',res=300)
p12 = ggscatter(data, x = "zcode_MMSE_Memory_Recall_VIS12", y = "SA_Flexibility_VIS12", 
                add = "reg.line", conf.int = TRUE, conf.int.level = 0.95, fullrange = TRUE,
                cor.coef = FALSE,  cor.method = "spearman",
                cor.coeff.args = list(method = "spearman", label.x.npc = "center", label.y.npc = "top"),
                xlab = "MMSE_Recall_VIS12", ylab = "Flexibility_VIS12", font.x = 16,
                font.y = 16) + annotate("text", x = 2, y = 88, label = label12,  size=5) 
p12
dev.off()

##confirmend scatter plots with Altoida
library(cowplot)
##nested plot grids
p1_plt = plot_grid(p1, p2, p3, align = "hv",rel_widths = c(1, 1, 1), ncol = 3, nrow = 1)+theme(plot.margin = unit(c(1, 1, 1, 1), "cm"))
p2_plt = plot_grid(p4, p5, p6,align = "hv",rel_widths = c(1, 1, 1), ncol = 3, nrow = 1)+theme(plot.margin = unit(c(1, 1, 1, 1), "cm"))
p3_plt = plot_grid(p7, p8, p9, align = "hv",rel_widths = c(1, 1, 1), ncol = 3, nrow = 1)+theme(plot.margin = unit(c(1, 1, 1, 1), "cm"))
p4_plt = plot_grid(p10, p11,p12,align = "hv",rel_widths = c(1, 1,1), ncol = 3, nrow = 1)+theme(plot.margin = unit(c(1, 1, 1, 1), "cm"))
png("scatterPlots/scatPlotAllconfirmed.png", width = 500, height = 420, units='mm',res = 300) 
plot_grid(p1_plt, p2_plt, p3_plt, p4_plt, labels = c('A', 'B', 'C', 'D'), label_size = 20, ncol = 1)
dev.off()

##other scatter plots not common to Altoida
res12 <-cor.test(data$zcode_MMSE_Orientation_VIS1, data$SA_zcode_BITDOTMotorInstructionReadingTimeRatios_VIS1,  method = "spearman")
res12
r12 = paste(spearmanRankPadj[spearmanRankPadj$from=="zcode_MMSE_Orientation_VIS1"&spearmanRankPadj$to=="SA_zcode_BITDOTMotorInstructionReadingTimeRatios_VIS1","r"])
r12_tx = paste("R =", r12, sep = " ")
low_int12= spearmanRankPadj[spearmanRankPadj$from=="zcode_MMSE_Orientation_VIS1"&spearmanRankPadj$to=="SA_zcode_BITDOTMotorInstructionReadingTimeRatios_VIS1","lower.adj"]
upper_int12= spearmanRankPadj[spearmanRankPadj$from=="zcode_MMSE_Orientation_VIS1"&spearmanRankPadj$to=="SA_zcode_BITDOTMotorInstructionReadingTimeRatios_VIS1","upper.adj"]
conf_int_tx12 = paste("95% CI [", low_int12, ",",upper_int12, "]", sep = "")
p_adj12 = spearmanRankPadj[spearmanRankPadj$from=="zcode_MMSE_Orientation_VIS1"&spearmanRankPadj$to=="SA_zcode_BITDOTMotorInstructionReadingTimeRatios_VIS1","p.Adj"]
p_adj_tx12 = paste("p.Adj =", "<0.0001", sep = " ")
label12  = paste(r12_tx, conf_int_tx12, p_adj_tx12, sep = ", ")
png(filename="scatterPlots/scatterOrienBITDOT_vis1.png", width = 7, height = 7, units = 'in',res=300)
p12 = ggscatter(data, x = "zcode_MMSE_Orientation_VIS1", y = "SA_zcode_BITDOTMotorInstructionReadingTimeRatios_VIS1", 
                add = "reg.line", conf.int = TRUE, conf.int.level = 0.95, fullrange = TRUE,
                cor.coef = FALSE,  cor.method = "spearman",
                cor.coeff.args = list(method = "spearman", label.x.npc = "center", label.y.npc = "top"),
                xlab = "MMSE_Orientation_VIS1", ylab = "BITDOTMotorInstructionReadingTimeRatios_VIS1", font.x = 16,
                font.y = 16) + annotate("text", x = 4, y = 3, label = label12,  size=5) 
p12
dev.off()

res13 <-cor.test(data$zcode_MMSE_Orientation_VIS6, data$SA_zcode_BITDOTMotorInstructionReadingTimeRatios_VIS6,  method = "spearman")
res13
r13 = paste(spearmanRankPadjGoodCor[spearmanRankPadjGoodCor$from=="zcode_MMSE_Orientation_VIS6"&spearmanRankPadjGoodCor$to=="SA_zcode_BITDOTMotorInstructionReadingTimeRatios_VIS6","r"])
r13_tx = paste("R =", r13, sep = " ")
low_int13= spearmanRankPadjGoodCor[spearmanRankPadjGoodCor$from=="zcode_MMSE_Orientation_VIS6"&spearmanRankPadjGoodCor$to=="SA_zcode_BITDOTMotorInstructionReadingTimeRatios_VIS6","lower.adj"]
upper_int13= spearmanRankPadjGoodCor[spearmanRankPadjGoodCor$from=="zcode_MMSE_Orientation_VIS6"&spearmanRankPadjGoodCor$to=="SA_zcode_BITDOTMotorInstructionReadingTimeRatios_VIS6","upper.adj"]
conf_int_tx13 = paste("95% CI [", low_int13, ",",upper_int13, "]", sep = "")
p_adj13 = spearmanRankPadjGoodCor[spearmanRankPadjGoodCor$from=="zcode_MMSE_Orientation_VIS6"&spearmanRankPadjGoodCor$to=="SA_zcode_BITDOTMotorInstructionReadingTimeRatios_VIS6","p.Adj"]
p_adj_tx13 = paste("p.Adj =", "<0.0001", sep = " ")
label13  = paste(r13_tx, conf_int_tx13, p_adj_tx13, sep = ", ")
png(filename="scatterPlots/scatterOrienBITDOT_vis6.png", width = 7, height = 7, units = 'in',res=300)
p13 = ggscatter(data, x = "zcode_MMSE_Orientation_VIS6", y = "SA_zcode_BITDOTMotorInstructionReadingTimeRatios_VIS6", 
                add = "reg.line", conf.int = TRUE, conf.int.level = 0.95, fullrange = TRUE,
                cor.coef = FALSE,  cor.method = "spearman",
                cor.coeff.args = list(method = "spearman", label.x.npc = "center", label.y.npc = "top"),
                xlab = "MMSE_Orientation_VIS6", ylab = "BITDOTMotorInstructionReadingTimeRatios_VIS6", font.x = 16,
                font.y = 16) + annotate("text", x = 5, y = 4, label = label13,  size=5) 
p13
dev.off()

res14 <-cor.test(data$zcode_MMSE_Orientation_VIS12, data$SA_zcode_BITDOTMotorInstructionReadingTimeRatios_VIS12,  method = "spearman")
res14
r14 = paste(spearmanRankPadj[spearmanRankPadj$from=="zcode_MMSE_Orientation_VIS12"&spearmanRankPadj$to=="SA_zcode_BITDOTMotorInstructionReadingTimeRatios_VIS12","r"])
r14_tx = paste("R =", r14, sep = " ")
low_int14= spearmanRankPadj[spearmanRankPadj$from=="zcode_MMSE_Orientation_VIS12"&spearmanRankPadj$to=="SA_zcode_BITDOTMotorInstructionReadingTimeRatios_VIS12","lower.adj"]
upper_int14= spearmanRankPadj[spearmanRankPadj$from=="zcode_MMSE_Orientation_VIS12"&spearmanRankPadj$to=="SA_zcode_BITDOTMotorInstructionReadingTimeRatios_VIS12","upper.adj"]
conf_int_tx14 = paste("95% CI [", low_int14, ",",upper_int14, "]", sep = "")
p_adj14 = spearmanRankPadj[spearmanRankPadj$from=="zcode_MMSE_Orientation_VIS12"&spearmanRankPadj$to=="SA_zcode_BITDOTMotorInstructionReadingTimeRatios_VIS12","p.Adj"]
p_adj_tx14 = paste("p.Adj =", "<0.0001", sep = " ")
label14  = paste(r14_tx, conf_int_tx14, p_adj_tx14, sep = ", ")
png(filename="scatterPlots/scatterOrienBITDOT_vis12.png", width = 7, height = 7, units = 'in',res=300)
p14 = ggscatter(data, x = "zcode_MMSE_Orientation_VIS12", y = "SA_zcode_BITDOTMotorInstructionReadingTimeRatios_VIS12", 
                add = "reg.line", conf.int = TRUE, conf.int.level = 0.95, fullrange = TRUE,
                cor.coef = FALSE,  cor.method = "spearman",
                cor.coeff.args = list(method = "spearman", label.x.npc = "center", label.y.npc = "top"),
                xlab = "MMSE_Orientation_VIS12", ylab = "BITDOTMotorInstructionReadingTimeRatios_VIS12", font.x = 16,
                font.y = 16) + annotate("text", x = 3.5, y = 3, label = label14,  size=5) 

p14
dev.off()

res15 <-cor.test(data$zcode_MMSE_Language_VIS6, data$SA_SpatialMemory_VIS6,  method = "spearman")
res15
r15 = paste(spearmanRankPadj[spearmanRankPadj$from=="zcode_MMSE_Language_VIS6"&spearmanRankPadj$to=="SA_SpatialMemory_VIS6","r"])
r15_tx = paste("R =", r15, sep = " ")
low_int15= spearmanRankPadj[spearmanRankPadj$from=="zcode_MMSE_Language_VIS6"&spearmanRankPadj$to=="SA_SpatialMemory_VIS6","lower.adj"]
upper_int15= spearmanRankPadj[spearmanRankPadj$from=="zcode_MMSE_Language_VIS6"&spearmanRankPadj$to=="SA_SpatialMemory_VIS6","upper.adj"]
conf_int_tx15 = paste("95% CI [", low_int15, ",",upper_int15, "]", sep = "")
p_adj15 = spearmanRankPadj[spearmanRankPadj$from=="zcode_MMSE_Language_VIS6"&spearmanRankPadj$to=="SA_SpatialMemory_VIS6","p.Adj"]
p_adj_tx15 = paste("p.Adj =", "<0.0001", sep = " ")
label15  = paste(r15_tx, conf_int_tx15, p_adj_tx15, sep = ", ")
png(filename="scatterPlots/scatterLangSptMem_vis6.png", width = 7, height = 7, units = 'in',res=300)
p15 = ggscatter(data, x = "zcode_MMSE_Language_VIS6", y = "SA_SpatialMemory_VIS6", 
                add = "reg.line", conf.int = TRUE, conf.int.level = 0.95, fullrange = TRUE,
                cor.coef = FALSE,  cor.method = "spearman",
                cor.coeff.args = list(method = "spearman", label.x.npc = "center", label.y.npc = "top"),
                xlab = "MMSE_Language_VIS6", ylab = "SpatialMemory_VIS6", font.x = 16,
                font.y = 16) + annotate("text", x = 4, y = 70, label = label15,  size=5) 

p15
dev.off()

res16 <-cor.test(data$zcode_MMSE_Language_VIS1, data$SA_ProspectiveMemory_VIS1,  method = "spearman")
res16
r16 = paste(spearmanRankPadj[spearmanRankPadj$from=="zcode_MMSE_Language_VIS1"&spearmanRankPadj$to=="SA_ProspectiveMemory_VIS1","r"])
r16_tx = paste("R =", r16, sep = " ")
low_int16= spearmanRankPadj[spearmanRankPadj$from=="zcode_MMSE_Language_VIS1"&spearmanRankPadj$to=="SA_ProspectiveMemory_VIS1","lower.adj"]
upper_int16= spearmanRankPadj[spearmanRankPadj$from=="zcode_MMSE_Language_VIS1"&spearmanRankPadj$to=="SA_ProspectiveMemory_VIS1","upper.adj"]
conf_int_tx16 = paste("95% CI [", low_int16, ",",upper_int16, "]", sep = "")
p_adj16 = spearmanRankPadj[spearmanRankPadj$from=="zcode_MMSE_Language_VIS1"&spearmanRankPadj$to=="SA_ProspectiveMemory_VIS1","p.Adj"]
p_adj_tx16 = paste("p.Adj =", "<0.005", sep = " ")
label16  = paste(r16_tx, conf_int_tx16, p_adj_tx16, sep = ", ")
png(filename="scatterPlots/scatterLangPerspMem_vis1.png", width = 7, height = 7, units = 'in',res=300)
p16 = ggscatter(data, x = "zcode_MMSE_Language_VIS1", y = "SA_ProspectiveMemory_VIS1", 
                add = "reg.line", conf.int = TRUE, conf.int.level = 0.95, fullrange = TRUE,
                cor.coef = FALSE,  cor.method = "spearman",
                cor.coeff.args = list(method = "spearman", label.x.npc = "center", label.y.npc = "top"),
                xlab = "MMSE_Language_VIS1", ylab = "ProspectiveMemory_VIS1", font.x = 16,
                font.y = 16) + annotate("text", x = -5, y =70, label = label16,  size=5) 

p16
dev.off()

res16_2 <-cor.test(data$zcode_MMSE_Language_VIS1, data$SA_PerceptualMotorCoordination_VIS1,  method = "spearman")
res16_2
r16_2 = paste(spearmanRankPadj[spearmanRankPadj$from=="zcode_MMSE_Language_VIS1"&spearmanRankPadj$to=="SA_PerceptualMotorCoordination_VIS1","r"])
r16_2_tx = paste("R =", r16_2, sep = " ")
low_int16_2= spearmanRankPadj[spearmanRankPadj$from=="zcode_MMSE_Language_VIS1"&spearmanRankPadj$to=="SA_PerceptualMotorCoordination_VIS1","lower.adj"]
upper_int16_2= spearmanRankPadj[spearmanRankPadj$from=="zcode_MMSE_Language_VIS1"&spearmanRankPadj$to=="SA_PerceptualMotorCoordination_VIS1","upper.adj"]
conf_int_tx16_2 = paste("95% CI [", low_int16_2, ",",upper_int16_2, "]", sep = "")
p_adj16_2 = spearmanRankPadj[spearmanRankPadj$from=="zcode_MMSE_Language_VIS1"&spearmanRankPadj$to=="SA_PerceptualMotorCoordination_VIS1","p.Adj"]
p_adj_tx16_2 = paste("p.Adj =", "<0.0001", sep = " ")
label16_2  = paste(r16_2_tx, conf_int_tx16_2, p_adj_tx16_2, sep = ", ")
png(filename="scatterPlots/scatterLangPercMotorCoord_vis1.png", width = 7, height = 7, units = 'in',res=300)
p16_2 = ggscatter(data, x = "zcode_MMSE_Language_VIS1", y = "SA_PerceptualMotorCoordination_VIS1", 
                  add = "reg.line", conf.int = TRUE, conf.int.level = 0.95, fullrange = TRUE,
                  cor.coef = FALSE,  cor.method = "spearman",
                  cor.coeff.args = list(method = "spearman", label.x.npc = "center", label.y.npc = "top"),
                  xlab = "MMSE_Language_VIS1", ylab = "PerceptualMotorCoordination_VIS1", font.x = 16,
                  font.y = 16) + annotate("text", x = -5, y =70, label = label16_2,  size=5) 

p16_2
dev.off()

res16_3 <-cor.test(data$zcode_MMSE_Language_VIS12, data$SA_SpatialMemory_VIS12,  method = "spearman")
res16_3
r16_3 = paste(spearmanRankPadj[spearmanRankPadj$from=="zcode_MMSE_Language_VIS12"&spearmanRankPadj$to=="SA_SpatialMemory_VIS12","r"])
r16_3_tx = paste("R =", r16_3, sep = " ")
low_int16_3= spearmanRankPadj[spearmanRankPadj$from=="zcode_MMSE_Language_VIS12"&spearmanRankPadj$to=="SA_SpatialMemory_VIS12","lower.adj"]
upper_int16_3= spearmanRankPadj[spearmanRankPadj$from=="zcode_MMSE_Language_VIS12"&spearmanRankPadj$to=="SA_SpatialMemory_VIS12","upper.adj"]
conf_int_tx16_3 = paste("95% CI [", low_int16_3, ",",upper_int16_3, "]", sep = "")
p_adj16_3 = spearmanRankPadj[spearmanRankPadj$from=="zcode_MMSE_Language_VIS12"&spearmanRankPadj$to=="SA_SpatialMemory_VIS12","p.Adj"]
p_adj_tx16_3 = paste("p.Adj =", "<0.0001", sep = " ")
label16_3  = paste(r16_3_tx, conf_int_tx16_3, p_adj_tx16_3, sep = ", ")
png(filename="scatterPlots/scatterLangSpatMem_vis12.png", width = 7, height = 7, units = 'in',res=300)
p16_3 = ggscatter(data, x = "zcode_MMSE_Language_VIS12", y = "SA_SpatialMemory_VIS12", 
                  add = "reg.line", conf.int = TRUE, conf.int.level = 0.95, fullrange = TRUE,
                  cor.coef = FALSE,  cor.method = "spearman",
                  cor.coeff.args = list(method = "spearman", label.x.npc = "center", label.y.npc = "top"),
                  xlab = "MMSE_Language_VIS12", ylab = "SpatialMemory_VIS12", font.x = 16,
                  font.y = 16) + annotate("text", x = 3.5, y =70, label = label16_3,  size=5) 

p16_3
dev.off()


###faq and dig and cog##
##to include in the paper 
#1.yes
res17 <-cor.test(data$SA_FAQFORM_VIS1, data$SA_Flexibility_VIS1,  method = "spearman")
res17
r17 = paste(spearmanRankPadj[spearmanRankPadj$from=="SA_FAQFORM_VIS1"&spearmanRankPadj$to=="SA_Flexibility_VIS1","r"])
r17_tx = paste("R =", r17, sep = " ")
low_int17= spearmanRankPadj[spearmanRankPadj$from=="SA_FAQFORM_VIS1"&spearmanRankPadj$to=="SA_Flexibility_VIS1","lower.adj"]
upper_int17= spearmanRankPadj[spearmanRankPadj$from=="SA_FAQFORM_VIS1"&spearmanRankPadj$to=="SA_Flexibility_VIS1","upper.adj"]
conf_int_tx17 = paste("95% CI [", low_int17, ",",upper_int17, "]", sep = "")
p_adj17 = spearmanRankPadj[spearmanRankPadj$from=="SA_FAQFORM_VIS1"&spearmanRankPadj$to=="SA_Flexibility_VIS1","p.Adj"]
p_adj_tx17 = paste("p.Adj =",  "<0.0001", sep = " ")
label17  = paste(r17_tx, conf_int_tx17, p_adj_tx17, sep = ", ")
png(filename="scatterPlots/scatterFAQFORMFlex_vis1.png", width = 7, height = 7, units = 'in',res=300)
p17 = ggscatter(data, x = "SA_FAQFORM_VIS1", y = "SA_Flexibility_VIS1", 
                add = "reg.line", conf.int = TRUE, conf.int.level = 0.95, fullrange = TRUE,
                cor.coef = FALSE,  cor.method = "spearman",
                cor.coeff.args = list(method = "spearman", label.x.npc = "center", label.y.npc = "top"),
                xlab = "FAQFORM_VIS1", ylab = "Flexibility_VIS1", font.x = 16,
                font.y = 16) + annotate("text", x = 2.5, y = 90, label = label17,  size=5) 

p17
dev.off()

##yes
res17_1 <-cor.test(data$SA_FAQFORM_VIS1, data$SA_Flexibility_VIS6,  method = "spearman")
res17_1
r17_1 = paste(spearmanRankPadj[spearmanRankPadj$from=="SA_FAQFORM_VIS1"&spearmanRankPadj$to=="SA_Flexibility_VIS6","r"])
r17_1_tx = paste("R =", r17_1, sep = " ")
low_int17_1= spearmanRankPadj[spearmanRankPadj$from=="SA_FAQFORM_VIS1"&spearmanRankPadj$to=="SA_Flexibility_VIS6","lower.adj"]
upper_int17_1= spearmanRankPadj[spearmanRankPadj$from=="SA_FAQFORM_VIS1"&spearmanRankPadj$to=="SA_Flexibility_VIS6","upper.adj"]
conf_int_tx17_1 = paste("95% CI [", low_int17_1, ",",upper_int17_1, "]", sep = "")
p_adj17_1 = spearmanRankPadj[spearmanRankPadj$from=="SA_FAQFORM_VIS1"&spearmanRankPadj$to=="SA_Flexibility_VIS6","p.Adj"]
p_adj_tx17_1 = paste("p.Adj =",  "<0.0001", sep = " ")
label17_1  = paste(r17_1_tx, conf_int_tx17_1, p_adj_tx17_1, sep = ", ")
png(filename="scatterPlots/scatterFAQFORMFlex_vis6.png", width = 7, height = 7, units = 'in',res=300)
p17_1 = ggscatter(data, x = "SA_FAQFORM_VIS1", y = "SA_Flexibility_VIS6", 
                  add = "reg.line", conf.int = TRUE, conf.int.level = 0.95, fullrange = TRUE,
                  cor.coef = FALSE,  cor.method = "spearman",
                  cor.coeff.args = list(method = "spearman", label.x.npc = "center", label.y.npc = "top"),
                  xlab = "FAQFORM_VIS1", ylab = "Flexibility_VIS6", font.x = 16,
                  font.y = 16) + annotate("text", x = 2.5, y = 90, label = label17_1,  size=5) 


dev.off()

##yes
res17_2 <-cor.test(data$SA_FAQFORM_VIS1, data$SA_zcode_ARScreenButtonPresses_VIS1,  method = "spearman")
res17_2
r17_2 = paste(spearmanRankPadj[spearmanRankPadj$from=="SA_FAQFORM_VIS1"&spearmanRankPadj$to=="SA_zcode_ARScreenButtonPresses_VIS1","r"])
r17_2_tx = paste("R =", r17_2, sep = " ")
low_int17_2= spearmanRankPadj[spearmanRankPadj$from=="SA_FAQFORM_VIS1"&spearmanRankPadj$to=="SA_zcode_ARScreenButtonPresses_VIS1","lower.adj"]
upper_int17_2= spearmanRankPadj[spearmanRankPadj$from=="SA_FAQFORM_VIS1"&spearmanRankPadj$to=="SA_zcode_ARScreenButtonPresses_VIS1","upper.adj"]
conf_int_tx17_2 = paste("95% CI [", low_int17_2, ",",upper_int17_2, "]", sep = "")
p_adj17_2 = spearmanRankPadj[spearmanRankPadj$from=="SA_FAQFORM_VIS1"&spearmanRankPadj$to=="SA_zcode_ARScreenButtonPresses_VIS1","p.Adj"]
p_adj_tx17_2 = paste("p.Adj =",  "<0.0001", sep = " ")
label17_2  = paste(r17_2_tx, conf_int_tx17_2, p_adj_tx17_2, sep = ", ")
png(filename="scatterPlots/scatterFAQFORMARScreenButton_vis1.png", width = 7, height = 7, units = 'in',res=300)
p17_2 = ggscatter(data, x = "SA_FAQFORM_VIS1", y = "SA_zcode_ARScreenButtonPresses_VIS1", 
                  add = "reg.line", conf.int = TRUE, conf.int.level = 0.95, fullrange = TRUE,
                  cor.coef = FALSE,  cor.method = "spearman",
                  cor.coeff.args = list(method = "spearman", label.x.npc = "center", label.y.npc = "top"),
                  xlab = "FAQFORM_VIS1", ylab = "ARScreenButtonPresses_VIS1", font.x = 16,
                  font.y = 16) + annotate("text", x = 2.5, y = 5, label = label17_2,  size=5) 

p17_2
dev.off()



##5. Yes
res23 <-cor.test(data$SA_FAQFORM_VIS24, data$SA_PerceptualMotorCoordination_VIS36,  method = "spearman")
res23
r23 = paste(spearmanRankPadj[spearmanRankPadj$from=="SA_FAQFORM_VIS24"&spearmanRankPadj$to=="SA_PerceptualMotorCoordination_VIS36","r"])
r23_tx = paste("R =", r23, sep = " ")
low_int23= spearmanRankPadj[spearmanRankPadj$from=="SA_FAQFORM_VIS24"&spearmanRankPadj$to=="SA_PerceptualMotorCoordination_VIS36","lower.adj"]
upper_int23= spearmanRankPadj[spearmanRankPadj$from=="SA_FAQFORM_VIS24"&spearmanRankPadj$to=="SA_PerceptualMotorCoordination_VIS36","upper.adj"]
conf_int_tx23 = paste("95% CI [", low_int23, ",",upper_int23, "]", sep = "")
p_adj23 = spearmanRankPadj[spearmanRankPadj$from=="SA_FAQFORM_VIS24"&spearmanRankPadj$to=="SA_PerceptualMotorCoordination_VIS36","p.Adj"]
p_adj_tx23 = paste("p.Adj =",  "<0.005", sep = " ")
label23  = paste(r23_tx, conf_int_tx23, p_adj_tx23, sep = ", ")
png(filename="scatterPlots/scatterFAQFORMPerceptualMotorCoord_VIS36.png", width = 7, height = 7, units = 'in',res=300)
p23 = ggscatter(data, x = "SA_FAQFORM_VIS24", y = "SA_PerceptualMotorCoordination_VIS36", 
                add = "reg.line", conf.int = TRUE, conf.int.level = 0.95, fullrange = TRUE,
                cor.coef = FALSE,  cor.method = "spearman",
                cor.coeff.args = list(method = "spearman", label.x.npc = "center", label.y.npc = "top"),
                xlab = "FAQFORM_VIS24", ylab = "PerceptualMotorCoordination_VIS36", font.x = 16,
                font.y = 16) + annotate("text", x = 2.5, y = 63, label = label23,  size=5) 

p23
dev.off()




##6. yes
res24 <-cor.test(data$SA_FAQFORM_VIS24, data$SA_Planning_VIS36,  method = "spearman")
res24
r24 = paste(spearmanRankPadj[spearmanRankPadj$from=="SA_FAQFORM_VIS24"&spearmanRankPadj$to=="SA_Planning_VIS36","r"])
r24_tx = paste("R =", r24, sep = " ")
low_int24= spearmanRankPadj[spearmanRankPadj$from=="SA_FAQFORM_VIS24"&spearmanRankPadj$to=="SA_Planning_VIS36","lower.adj"]
upper_int24= spearmanRankPadj[spearmanRankPadj$from=="SA_FAQFORM_VIS24"&spearmanRankPadj$to=="SA_Planning_VIS36","upper.adj"]
conf_int_tx24 = paste("95% CI [", low_int24, ",",upper_int24, "]", sep = "")
p_adj24 = spearmanRankPadj[spearmanRankPadj$from=="SA_FAQFORM_VIS24"&spearmanRankPadj$to=="SA_Planning_VIS36","p.Adj"]
p_adj_tx24 = paste("p.Adj =",  "<0.0001", sep = " ")
label24  = paste(r24_tx, conf_int_tx24, p_adj_tx24, sep = ", ")
png(filename="scatterPlots/scatterFAQFORMPlan_VIS36.png", width = 7, height = 7, units = 'in',res=300)
p24 = ggscatter(data, x = "SA_FAQFORM_VIS24", y = "SA_Planning_VIS36", 
                add = "reg.line", conf.int = TRUE, conf.int.level = 0.95, fullrange = TRUE,
                cor.coef = FALSE,  cor.method = "spearman",
                cor.coeff.args = list(method = "spearman", label.x.npc = "center", label.y.npc = "top"),
                xlab = "FAQFORM_VIS24", ylab = "Planning_VIS36", font.x = 16,
                font.y = 16) + annotate("text", x = 2.5, y = 65, label = label24,  size=5) 

p24
dev.off()

#7. yes
res25 <-cor.test(data$SA_FAQREM_VIS24, data$SA_zcode_MotorTestDurations_VIS36,  method = "spearman")
res25
r25 = paste(spearmanRankPadj[spearmanRankPadj$from=="SA_FAQREM_VIS24"&spearmanRankPadj$to=="SA_zcode_MotorTestDurations_VIS36","r"])
r25_tx = paste("R =", r25, sep = " ")
low_int25= spearmanRankPadj[spearmanRankPadj$from=="SA_FAQREM_VIS24"&spearmanRankPadj$to=="SA_zcode_MotorTestDurations_VIS36","lower.adj"]
upper_int25= spearmanRankPadj[spearmanRankPadj$from=="SA_FAQREM_VIS24"&spearmanRankPadj$to=="SA_zcode_MotorTestDurations_VIS36","upper.adj"]
conf_int_tx25 = paste("95% CI [", low_int25, ",",upper_int25, "]", sep = "")
p_adj25 = spearmanRankPadj[spearmanRankPadj$from=="SA_FAQREM_VIS24"&spearmanRankPadj$to=="SA_zcode_MotorTestDurations_VIS36","p.Adj"]
p_adj_tx25 = paste("p.Adj =",  "<0.005", sep = " ")
label25  = paste(r25_tx, conf_int_tx25, p_adj_tx25, sep = ", ")
png(filename="scatterPlots/scatterFAQREM24MotorTest_VIS36.png", width = 7, height = 7, units = 'in',res=300)
p25 = ggscatter(data, x = "SA_FAQREM_VIS24", y = "SA_zcode_MotorTestDurations_VIS36", 
                add = "reg.line", conf.int = TRUE, conf.int.level = 0.95, fullrange = TRUE,
                cor.coef = FALSE,  cor.method = "spearman",
                cor.coeff.args = list(method = "spearman", label.x.npc = "center", label.y.npc = "top"),
                xlab = "FAQREM_VIS24", ylab = "MotorTestDurations_VIS36", font.x = 16,
                font.y = 16) + annotate("text", x = 2.5, y = 5, label = label25,  size=5) 

p25
dev.off()

##7_1 yes
res25_1 <-cor.test(data$SA_FAQREM_VIS24, data$SA_zcode_ARObjectFinding_VIS24,  method = "spearman")
res25_1
r25_1 = paste(spearmanRankPadj[spearmanRankPadj$from=="SA_FAQREM_VIS24"&spearmanRankPadj$to=="SA_zcode_ARObjectFinding_VIS24","r"])
r25_1_tx = paste("R =", r25_1, sep = " ")
low_int25_1= spearmanRankPadj[spearmanRankPadj$from=="SA_FAQREM_VIS24"&spearmanRankPadj$to=="SA_zcode_ARObjectFinding_VIS24","lower.adj"]
upper_int25_1= spearmanRankPadj[spearmanRankPadj$from=="SA_FAQREM_VIS24"&spearmanRankPadj$to=="SA_zcode_ARObjectFinding_VIS24","upper.adj"]
conf_int_tx25_1 = paste("95% CI [", low_int25_1, ",",upper_int25_1, "]", sep = "")
p_adj25_1 = spearmanRankPadj[spearmanRankPadj$from=="SA_FAQREM_VIS24"&spearmanRankPadj$to=="SA_zcode_ARObjectFinding_VIS24","p.Adj"]
p_adj_tx25_1 = paste("p.Adj =",  "<0.0001", sep = " ")
label25_1  = paste(r25_1_tx, conf_int_tx25_1, p_adj_tx25_1, sep = ", ")
png(filename="scatterPlots/scatterFAQREM24ARObjectFinding_VIS24.png", width = 7, height = 7, units = 'in',res=300)
p25_1 = ggscatter(data, x = "SA_FAQREM_VIS24", y = "SA_zcode_ARObjectFinding_VIS24", 
                  add = "reg.line", conf.int = TRUE, conf.int.level = 0.95, fullrange = TRUE,
                  cor.coef = FALSE,  cor.method = "spearman",
                  cor.coeff.args = list(method = "spearman", label.x.npc = "center", label.y.npc = "top"),
                  xlab = "FAQREM_VIS24", ylab = "ARObjectFinding_VIS24", font.x = 16,
                  font.y = 16) + annotate("text", x = 2.5, y = 5, label = label25_1,  size=5) 

p25_1
dev.off()


#8. yes
res26 <-cor.test(data$SA_FAQFINAN_VIS12, data$SA_zcode_ARObjectFinding_VIS24,  method = "spearman")
res26
r26 = paste(spearmanRankPadj[spearmanRankPadj$from=="SA_FAQFINAN_VIS12"&spearmanRankPadj$to=="SA_zcode_ARObjectFinding_VIS24","r"])
r26_tx = paste("R =", r26, sep = " ")
low_int26= spearmanRankPadj[spearmanRankPadj$from=="SA_FAQFINAN_VIS12"&spearmanRankPadj$to=="SA_zcode_ARObjectFinding_VIS24","lower.adj"]
upper_int26= spearmanRankPadj[spearmanRankPadj$from=="SA_FAQFINAN_VIS12"&spearmanRankPadj$to=="SA_zcode_ARObjectFinding_VIS24","upper.adj"]
conf_int_tx26 = paste("95% CI [", low_int26, ",",upper_int26, "]", sep = "")
p_adj26 = spearmanRankPadj[spearmanRankPadj$from=="SA_FAQFINAN_VIS12"&spearmanRankPadj$to=="SA_zcode_ARObjectFinding_VIS24","p.Adj"]
p_adj_tx26 = paste("p.Adj =",  "<0.0001", sep = " ")
label26  = paste(r26_tx, conf_int_tx26, p_adj_tx26, sep = ", ")
png(filename="scatterPlots/scatterFAQFINAN12ARObjFind24.png", width = 7, height = 7, units = 'in',res=300)
p26 = ggscatter(data, x = "SA_FAQFINAN_VIS12", y = "SA_zcode_ARObjectFinding_VIS24", 
                add = "reg.line", conf.int = TRUE, conf.int.level = 0.95, fullrange = TRUE,
                cor.coef = FALSE,  cor.method = "spearman",
                cor.coeff.args = list(method = "spearman", label.x.npc = "center", label.y.npc = "top"),
                xlab = "FAQFINAN_VIS12", ylab = "ARObjectFinding_VIS24", font.x = 16,
                font.y = 16) + annotate("text", x = 2.5, y = 4.5, label = label26,  size=5) 

p26
dev.off()

#9. yes
res26_1 <-cor.test(data$SA_FAQREM_VIS12, data$SA_Flexibility_VIS12,  method = "spearman")
res26_1
r26_1 = paste(spearmanRankPadj[spearmanRankPadj$from=="SA_FAQREM_VIS12"&spearmanRankPadj$to=="SA_Flexibility_VIS12","r"])
r26_1_tx = paste("R =", r26_1, sep = " ")
low_int26_1= spearmanRankPadj[spearmanRankPadj$from=="SA_FAQREM_VIS12"&spearmanRankPadj$to=="SA_Flexibility_VIS12","lower.adj"]
upper_int26_1= spearmanRankPadj[spearmanRankPadj$from=="SA_FAQREM_VIS12"&spearmanRankPadj$to=="SA_Flexibility_VIS12","upper.adj"]
conf_int_tx26_1 = paste("95% CI [", low_int26_1, ",",upper_int26_1, "]", sep = "")
p_adj26_1 = spearmanRankPadj[spearmanRankPadj$from=="SA_FAQREM_VIS12"&spearmanRankPadj$to=="SA_Flexibility_VIS12","p.Adj"]
p_adj_tx26_1 = paste("p.Adj =",  "<0.0001", sep = " ")
label26_1  = paste(r26_1_tx, conf_int_tx26_1, p_adj_tx26_1, sep = ", ")
png(filename="scatterPlots/scatterFAQREMFlexiVIS12.png", width = 7, height = 7, units = 'in',res=300)
p26_1 = ggscatter(data, x = "SA_FAQREM_VIS12", y = "SA_Flexibility_VIS12", 
                  add = "reg.line", conf.int = TRUE, conf.int.level = 0.95, fullrange = TRUE,
                  cor.coef = FALSE,  cor.method = "spearman",
                  cor.coeff.args = list(method = "spearman", label.x.npc = "center", label.y.npc = "top"),
                  xlab = "FAQREM_VIS12", ylab = "Flexibility_VIS12", font.x = 16,
                  font.y = 16) + annotate("text", x = 2.5, y = 90, label = label26_1,  size=5) 

p26_1
dev.off()


#10. yes
res26_2 <-cor.test(data$SA_FAQREM_VIS12, data$SA_Flexibility_VIS36,  method = "spearman")
res26_2
r26_2 = paste(spearmanRankPadj[spearmanRankPadj$from=="SA_FAQREM_VIS12"&spearmanRankPadj$to=="SA_Flexibility_VIS36","r"])
r26_2_tx = paste("R =", r26_2, sep = " ")
low_int26_2= spearmanRankPadj[spearmanRankPadj$from=="SA_FAQREM_VIS12"&spearmanRankPadj$to=="SA_Flexibility_VIS36","lower.adj"]
upper_int26_2= spearmanRankPadj[spearmanRankPadj$from=="SA_FAQREM_VIS12"&spearmanRankPadj$to=="SA_Flexibility_VIS36","upper.adj"]
conf_int_tx26_2 = paste("95% CI [", low_int26_2, ",",upper_int26_2, "]", sep = "")
p_adj26_2 = spearmanRankPadj[spearmanRankPadj$from=="SA_FAQREM_VIS12"&spearmanRankPadj$to=="SA_Flexibility_VIS36","p.Adj"]
p_adj_tx26_2 = paste("p.Adj =",  "<0.0001", sep = " ")
label26_2  = paste(r26_2_tx, conf_int_tx26_2, p_adj_tx26_2, sep = ", ")
png(filename="scatterPlots/scatterFAQREM12FlexiVIS36.png", width = 7, height = 7, units = 'in',res=300)
p26_2 = ggscatter(data, x = "SA_FAQREM_VIS12", y = "SA_Flexibility_VIS36", 
                  add = "reg.line", conf.int = TRUE, conf.int.level = 0.95, fullrange = TRUE,
                  cor.coef = FALSE,  cor.method = "spearman",
                  cor.coeff.args = list(method = "spearman", label.x.npc = "center", label.y.npc = "top"),
                  xlab = "FAQREM_VIS12", ylab = "Flexibility_VIS36", font.x = 16,
                  font.y = 16) + annotate("text", x = 2.5, y = 100, label = label26_2,  size=5) 

p26_2
dev.off()

#10. yes
res26_3 <-cor.test(data$SA_FAQFORM_VIS24, data$SA_zcode_ARObjectFinding_VIS36,  method = "spearman")
res26_3
r26_3 = paste(spearmanRankPadj[spearmanRankPadj$from=="SA_FAQFORM_VIS24"&spearmanRankPadj$to=="SA_zcode_ARObjectFinding_VIS36","r"])
r26_3_tx = paste("R =", r26_3, sep = " ")
low_int26_3= spearmanRankPadj[spearmanRankPadj$from=="SA_FAQFORM_VIS24"&spearmanRankPadj$to=="SA_zcode_ARObjectFinding_VIS36","lower.adj"]
upper_int26_3= spearmanRankPadj[spearmanRankPadj$from=="SA_FAQFORM_VIS24"&spearmanRankPadj$to=="SA_zcode_ARObjectFinding_VIS36","upper.adj"]
conf_int_tx26_3 = paste("95% CI [", low_int26_3, ",",upper_int26_3, "]", sep = "")
p_adj26_3 = spearmanRankPadj[spearmanRankPadj$from=="SA_FAQFORM_VIS24"&spearmanRankPadj$to=="SA_zcode_ARObjectFinding_VIS36","p.Adj"]
p_adj_tx26_3 = paste("p.Adj =",  "<0.0001", sep = " ")
label26_3  = paste(r26_3_tx, conf_int_tx26_3, p_adj_tx26_3, sep = ", ")
png(filename="scatterPlots/scatterFAQFORM2ARObjectFinding_VIS3636.png", width = 7, height = 7, units = 'in',res=300)
p26_3 = ggscatter(data, x = "SA_FAQREM_VIS12", y = "SA_zcode_ARObjectFinding_VIS36", 
                  add = "reg.line", conf.int = TRUE, conf.int.level = 0.95, fullrange = TRUE,
                  cor.coef = FALSE,  cor.method = "spearman",
                  cor.coeff.args = list(method = "spearman", label.x.npc = "center", label.y.npc = "top"),
                  xlab = "FAQFORM_VIS24", ylab = "ARObjectFinding_VIS36", font.x = 16,
                  font.y = 16) + annotate("text", x = 2.5, y = 5, label = label26_3,  size=5) 

p26_3
dev.off()




##yes
res26_5 <-cor.test(data$SA_FAQREM_VIS12, data$SA_zcode_ARObjectFinding_VIS24,method = "spearman")
res26_5
r26_5 = paste(spearmanRankPadj[spearmanRankPadj$from=="SA_FAQREM_VIS12"&spearmanRankPadj$to=="SA_zcode_ARObjectFinding_VIS24","r"])
r26_5_tx = paste("R =", r26_5, sep = " ")
low_int26_5= spearmanRankPadj[spearmanRankPadj$from=="SA_FAQREM_VIS12"&spearmanRankPadj$to=="SA_zcode_ARObjectFinding_VIS24","lower.adj"]
upper_int26_5= spearmanRankPadj[spearmanRankPadj$from=="SA_FAQREM_VIS12"&spearmanRankPadj$to=="SA_zcode_ARObjectFinding_VIS24","upper.adj"]
conf_int_tx26_5 = paste("95% CI [", low_int26_5, ",",upper_int26_5, "]", sep = "")
p_adj26_5 = spearmanRankPadj[spearmanRankPadj$from=="SA_FAQREM_VIS12"&spearmanRankPadj$to=="SA_zcode_ARObjectFinding_VIS24","p.Adj"]
p_adj_tx26_5 = paste("p.Adj =",  "<0.0001", sep = " ")
label26_5  = paste(r26_5_tx, conf_int_tx26_5, p_adj_tx26_5, sep = ", ")
png(filename="scatterPlots/scatterFAQREM12ARObjectFinding_VIS24.png", width = 7, height = 7, units = 'in',res=300)
p26_5 = ggscatter(data, x = "SA_FAQREM_VIS24", y = "SA_zcode_ARObjectFinding_VIS24", 
                  add = "reg.line", conf.int = TRUE, conf.int.level = 0.95, fullrange = TRUE,
                  cor.coef = FALSE,  cor.method = "spearman",
                  cor.coeff.args = list(method = "spearman", label.x.npc = "center", label.y.npc = "top"),
                  xlab = "SA_FAQREM_VIS12", ylab = "ARObjectFinding_VIS24", font.x = 16,
                  font.y = 16) + annotate("text", x = 2.5, y = 5, label = label26_5,  size=5) 

p26_5
dev.off()

##nested plot grids, FAQ
p7_plt = plot_grid(p17,p17_1,p17_2,align = "hv",rel_widths = c(1, 1, 1), ncol = 3, nrow = 1)+theme(plot.margin = unit(c(2, 1, 1, 2), "cm"))
p8_plt = plot_grid(p23,p23_1,p26_3, align = "hv",rel_widths = c(1, 1, 1), ncol = 3, nrow = 1)+theme(plot.margin = unit(c(2, 1, 1, 2), "cm"))
p9_plt = plot_grid(p26_1,p26_2, align = "hv",rel_widths = c(1, 1, 1), ncol = 3, nrow = 1)+theme(plot.margin = unit(c(2, 1, 1, 2), "cm"))
p10_plt = plot_grid(p25, p25_1,align = "hv",rel_widths = c(1, 1), ncol = 3, nrow = 1)+theme(plot.margin = unit(c(2, 1, 1, 2), "cm"))
p11_plt = plot_grid(p26,align = "hv", ncol = 3, nrow = 1)+theme(plot.margin = unit(c(2, 1, 1, 2), "cm"))

png("scatterPlots/scatPlotFAQ.png", width =500, height = 700, units='mm',res = 300) 
plot_grid(p7_plt,p8_plt,p9_plt,p10_plt,p11_plt, rel_heights = c(3, 3, 3,3,3),labels = c('A', 'B', 'C','D','E'), label_size = 30, ncol = 1)
dev.off()


##digital with imaging##
res27 <-cor.test(data$zcode_volume_VIS1, data$SA_zcode_ARScreenButtonPresses_VIS1,  method = "spearman")
res27
r27 = paste(spearmanRankPadjGoodCor[spearmanRankPadjGoodCor$from=="zcode_volume_VIS1"&spearmanRankPadjGoodCor$to=="SA_zcode_ARScreenButtonPresses_VIS1","r"])
r27_tx = paste("R =", r27, sep = " ")
low_int27= spearmanRankPadjGoodCor[spearmanRankPadjGoodCor$from=="zcode_volume_VIS1"&spearmanRankPadjGoodCor$to=="SA_zcode_ARScreenButtonPresses_VIS1","lower.adj"]
upper_int27= spearmanRankPadjGoodCor[spearmanRankPadjGoodCor$from=="zcode_volume_VIS1"&spearmanRankPadjGoodCor$to=="SA_zcode_ARScreenButtonPresses_VIS1","upper.adj"]
conf_int_tx27 = paste("95% CI [", low_int27, ",",upper_int27, "]", sep = "")
p_adj27 = spearmanRankPadjGoodCor[spearmanRankPadjGoodCor$from=="zcode_volume_VIS1"&spearmanRankPadjGoodCor$to=="SA_zcode_ARScreenButtonPresses_VIS1","p.Adj"]
p_adj_tx27 = paste("p.Adj =", "<0.0001", sep = " ")
label27  = paste(r27_tx, conf_int_tx27, p_adj_tx27, sep = ", ")
png(filename="scatterPlots/volumeARScreenButtonPresses.png", width = 7, height = 7, units = 'in',res=300)
p27 = ggscatter(data, x = "zcode_volume_VIS1", y = "SA_zcode_ARScreenButtonPresses_VIS1", 
                add = "reg.line", conf.int = TRUE, conf.int.level = 0.95, fullrange = TRUE,
                cor.coef = FALSE,  cor.method = "spearman",
                cor.coeff.args = list(method = "spearman", label.x.npc = "center", label.y.npc = "top"),
                xlab = "volume_VIS1", ylab = "ARScreenButtonPresses_VIS1", font.x = 16,
                font.y = 16) + annotate("text", x = 0, y =3.5, label = label27,  size=5) 

p27
dev.off()



res28_1 <-cor.test(data$zcode_volume_VIS1, data$SA_ComplexAttention_VIS1,  method = "spearman")
res28_1
r28_1 = paste(spearmanRankPadj[spearmanRankPadj$from=="zcode_volume_VIS1"&spearmanRankPadj$to=="SA_ComplexAttention_VIS1","r"])
r28_1_tx = paste("R =", r28_1, sep = " ")
low_int28_1= spearmanRankPadj[spearmanRankPadj$from=="zcode_volume_VIS1"&spearmanRankPadj$to=="SA_ComplexAttention_VIS1","lower.adj"]
upper_int28_1= spearmanRankPadj[spearmanRankPadj$from=="zcode_volume_VIS1"&spearmanRankPadj$to=="SA_ComplexAttention_VIS1","upper.adj"]
conf_int_tx28_1 = paste("95% CI [", low_int28_1, ",",upper_int28_1, "]", sep = "")
p_adj28_1 = spearmanRankPadj[spearmanRankPadj$from=="zcode_volume_VIS1"&spearmanRankPadj$to=="SA_ComplexAttention_VIS1","p.Adj"]
p_adj_tx28_1 = paste("p.Adj =", "<0.0001", sep = " ")
label28_1  = paste(r28_1_tx, conf_int_tx28_1, p_adj_tx28_1, sep = ", ")
png(filename="scatterPlots/volumePercComplexAtt.png", width = 7, height = 7, units = 'in',res=300)
p28_1 = ggscatter(data, x = "zcode_volume_VIS1", y = "SA_ComplexAttention_VIS1", 
                  add = "reg.line", conf.int = TRUE, conf.int.level = 0.95, fullrange = TRUE,
                  cor.coef = FALSE,  cor.method = "spearman",
                  cor.coeff.args = list(method = "spearman", label.x.npc = "center", label.y.npc = "top"),
                  xlab = "volume_VIS1", ylab = "ComplexAttention_VIS1", font.x = 16,
                  font.y = 16) + annotate("text", x = 0, y =66, label = label28_1,  size=5) 

p28_1
dev.off()

##no
##PHS connected to object finding visit via Object Placement visit 1 (to Object finding 1-6-..)
res29 <-cor.test(data$SA_PHS_VIS1, data$SA_zcode_ARObjectFinding_VIS6,  method = "pearson")
res29
r29 = paste(spearmanRankPadj[spearmanRankPadj$from=="SA_PHS_VIS1"&spearmanRankPadj$to=="SA_zcode_ARObjectFinding_VIS6","r"])
r29_tx = paste("R =", r29, sep = " ")
low_int29= spearmanRankPadj[spearmanRankPadj$from=="SA_PHS_VIS1"&spearmanRankPadj$to=="SA_zcode_ARObjectFinding_VIS6","lower.adj"]
upper_int29= spearmanRankPadj[spearmanRankPadj$from=="SA_PHS_VIS1"&spearmanRankPadj$to=="SA_zcode_ARObjectFinding_VIS6","upper.adj"]
conf_int_tx29 = paste("95% CI [", low_int29, ",",upper_int29, "]", sep = "")
p_adj29 = spearmanRankPadj[spearmanRankPadj$from=="SA_PHS_VIS1"&spearmanRankPadj$to=="SA_zcode_ARObjectFinding_VIS6","p.Adj"]
p_adj_tx29 = paste("p.Adj =", "<0.0001", sep = " ")
label29  = paste(r29_tx, conf_int_tx29, p_adj_tx29, sep = ", ")
png(filename="scatterPlots/scatterPHSArObFin.png", width = 7, height = 7, units = 'in',res=300)
p29 = ggscatter(data, x = "SA_PHS_VIS1", y = "SA_zcode_ARObjectFinding_VIS6", 
                add = "reg.line", conf.int = TRUE, conf.int.level = 0.95, fullrange = TRUE,
                cor.coef = FALSE,  cor.method = "spearman",
                cor.coeff.args = list(method = "spearman", label.x.npc = "center", label.y.npc = "top"),
                xlab = "PHS_VIS1", ylab = "ARObjectFinding_VIS6", font.x = 16,
                font.y = 16) + annotate("text", x = 1, y =5, label = label29,  size=5) 


p29
dev.off()


##yes
res30 <-cor.test(as.numeric(data$zcode_csf_VIS1), data$SA_zcode_ARObjectPlacement_VIS1,  method = "spearman")
res30
r30 = paste(spearmanRankPadj[spearmanRankPadj$from=="zcode_csf_VIS1"&spearmanRankPadj$to=="SA_zcode_ARObjectPlacement_VIS1","r"])
r30_tx = paste("R =", r30, sep = " ")
low_int30= spearmanRankPadj[spearmanRankPadj$from=="zcode_csf_VIS1"&spearmanRankPadj$to=="SA_zcode_ARObjectPlacement_VIS1","lower.adj"]
upper_int30= spearmanRankPadj[spearmanRankPadj$from=="zcode_csf_VIS1"&spearmanRankPadj$to=="SA_zcode_ARObjectPlacement_VIS1","upper.adj"]
conf_int_tx30 = paste("95% CI [", low_int30, ",",upper_int30, "]", sep = "")
p_adj30 = spearmanRankPadj[spearmanRankPadj$from=="zcode_csf_VIS1"&spearmanRankPadj$to=="SA_zcode_ARObjectPlacement_VIS1","p.Adj"]
p_adj_tx30 = paste("p.Adj =", "<0.005", sep = " ")
label30  = paste(r30_tx, conf_int_tx30, p_adj_tx30, sep = ", ")
png(filename="scatterPlots/scattercsfARObjPlacement.png", width = 7, height = 7, units = 'in',res=300)
p30 = ggscatter(data, x = "zcode_csf_VIS1", y = "SA_zcode_ARObjectPlacement_VIS1", 
                add = "reg.line", conf.int = TRUE, conf.int.level = 0.95, fullrange = TRUE,
                cor.coef = FALSE,  cor.method = "spearman",
                cor.coeff.args = list(method = "spearman", label.x.npc = "center", label.y.npc = "top"),
                xlab = "csf_VIS1", ylab = "ARObjectPlacement_VIS1", font.x = 16,
                font.y = 16) + annotate("text", x = -0.5, y =4, label = label30,  size=5) 

p30
dev.off()

##yes
res31 <-cor.test(as.numeric(data$zcode_csf_VIS1), data$zcode_imagingPET_VIS1,  method = "spearman")
res31
r31 = paste(spearmanRankPadjGoodCor[spearmanRankPadjGoodCor$from=="zcode_csf_VIS1"&spearmanRankPadjGoodCor$to=="zcode_imagingPET_VIS1","r"])
r31_tx = paste("R =", r31, sep = " ")
low_int31= spearmanRankPadjGoodCor[spearmanRankPadjGoodCor$from=="zcode_csf_VIS1"&spearmanRankPadjGoodCor$to=="zcode_imagingPET_VIS1","lower.adj"]
upper_int31= spearmanRankPadjGoodCor[spearmanRankPadjGoodCor$from=="zcode_csf_VIS1"&spearmanRankPadjGoodCor$to=="zcode_imagingPET_VIS1","upper.adj"]
conf_int_tx31 = paste("95% CI [", low_int31, ",",upper_int31, "]", sep = "")
p_adj31 = spearmanRankPadjGoodCor[spearmanRankPadjGoodCor$from=="zcode_csf_VIS1"&spearmanRankPadjGoodCor$to=="zcode_imagingPET_VIS1","p.Adj"]
p_adj_tx31 = paste("p.Adj =", "<0.0001", sep = " ")
label31  = paste(r31_tx, conf_int_tx31, p_adj_tx31, sep = ", ")
png(filename="scatterPlots/scattercsfimaging.png", width = 7, height = 7, units = 'in',res=310)
p31 = ggscatter(data, x = "zcode_csf_VIS1", y = "zcode_imagingPET_VIS1", 
                add = "reg.line", conf.int = TRUE, conf.int.level = 0.95, fullrange = TRUE,
                cor.coef = FALSE,  cor.method = "spearman",
                cor.coeff.args = list(method = "spearman", label.x.npc = "center", label.y.npc = "top"),
                xlab = "csf_VIS1", ylab = "imagingPET_VIS1", font.x = 16,
                font.y = 16) + annotate("text", x = -0.5, y =4, label = label31,  size=5) 

p31
dev.off()


res32 <-cor.test(as.numeric(data$zcode_imagingPET_VIS1), data$zcode_volume_VIS1,  method = "spearman")
res32
r32 = paste(spearmanRankPadj[spearmanRankPadj$from=="zcode_imagingPET_VIS1"&spearmanRankPadj$to=="zcode_volume_VIS1","r"])
r32_tx = paste("R =", r32, sep = " ")
low_int32= spearmanRankPadj[spearmanRankPadj$from=="zcode_imagingPET_VIS1"&spearmanRankPadj$to=="zcode_volume_VIS1","lower.adj"]
upper_int32= spearmanRankPadj[spearmanRankPadj$from=="zcode_imagingPET_VIS1"&spearmanRankPadj$to=="zcode_volume_VIS1","upper.adj"]
conf_int_tx32 = paste("95% CI [", low_int32, ",",upper_int32, "]", sep = "")
p_adj32 = spearmanRankPadj[spearmanRankPadj$from=="zcode_imagingPET_VIS1"&spearmanRankPadj$to=="zcode_volume_VIS1","p.Adj"]
p_adj_tx32 = paste("p.Adj =", "<0.0001", sep = " ")
label32  = paste(r32_tx, conf_int_tx32, p_adj_tx32, sep = ", ")
png(filename="scatterPlots/scattervolimaging.png", width = 7, height = 7, units = 'in',res=320)
p32 = ggscatter(data, x = "SA_PHS_VIS1", y = "zcode_volume_VIS1", 
                add = "reg.line", conf.int = TRUE, conf.int.level = 0.95, fullrange = TRUE,
                cor.coef = FALSE,  cor.method = "spearman",
                cor.coeff.args = list(method = "spearman", label.x.npc = "center", label.y.npc = "top"),
                xlab = "PHS_VIS1", ylab = "zcode_volume_VIS1", font.x = 16,
                font.y = 16) + annotate("text", x = -0.5, y =4, label = label32,  size=5) 

p32
dev.off()

##not included as lower ci is lower than 0
# res31 <-cor.test(as.numeric(data$zcode_csf_VIS1), data$SA_CognitiveProcessingSpeed_VIS6,  method = "spearman")
# res31
# r31 = paste(spearmanRankPadj[spearmanRankPadj$from=="zcode_csf_VIS1"&spearmanRankPadj$to=="SA_CognitiveProcessingSpeed_VIS6","r"])
# r31_tx = paste("R =", r31, sep = " ")
# low_int31= spearmanRankPadj[spearmanRankPadj$from=="zcode_csf_VIS1"&spearmanRankPadj$to=="SA_CognitiveProcessingSpeed_VIS6","lower.adj"]
# upper_int31= spearmanRankPadj[spearmanRankPadj$from=="zcode_csf_VIS1"&spearmanRankPadj$to=="SA_CognitiveProcessingSpeed_VIS6","upper.adj"]
# conf_int_tx31 = paste("95% CI [", low_int31, ",",upper_int31, "]", sep = "")
# p_adj31 = spearmanRankPadj[spearmanRankPadj$from=="zcode_csf_VIS1"&spearmanRankPadj$to=="SA_CognitiveProcessingSpeed_VIS6","p.Adj"]
# p_adj_tx31 = paste("p.Adj =", "<0.005", sep = " ")
# label31  = paste(r31_tx, conf_int_tx31, p_adj_tx31, sep = ", ")
# png(filename="scatterPlots/scatterPHSArObFin.png", width = 7, height = 7, units = 'in',res=310)
# p31 = ggscatter(data, x = "zcode_csf_VIS1", y = "SA_CognitiveProcessingSpeed_VIS6", 
#                 add = "reg.line", conf.int = TRUE, conf.int.level = 0.95, fullrange = TRUE,
#                 cor.coef = FALSE,  cor.method = "spearman",
#                 cor.coeff.args = list(method = "spearman", label.x.npc = "center", label.y.npc = "top"),
#                 xlab = "csf_VIS1", ylab = "CognitiveProcessingSpeed_VIS6", font.x = 16,
#                 font.y = 16) + annotate("text", x = -5, y =70, label = label31,  size=5) 
# 
# p31
# dev.off()
##apoe4 to SA_zcode_ARObjectPlacementFFT_VIS6##(strength 0.37)

##new scatter plots with volume, csf and PHS
library(cowplot)
pear1_plt = plot_grid(p27, p28_1,p30, align = "hv",rel_widths = c(1, 1), ncol = 2, nrow = 1)+theme(plot.margin = unit(c(1,1,1,1), "cm"))
#pear2_plt = plot_grid(p29, p30, align = "hv",rel_widths = c(1, 1), ncol = 2, nrow = 1)+theme(plot.margin = unit(c(1,1,1,1), "cm"))
png("scatterPlots/scatPlotAllcsfVolumePHS.png", width =300, height = 300, units='mm',res = 300) 
plot_grid(pear1_plt,pear2_plt, rel_heights = c(1,1),labels = c('A', 'B'), label_size = 20, ncol = 1)
dev.off()


##SA_DX_VIS1 to Dig measures
booot.stren.1000 = read.csv("boot_stren_all.csv")
boot.strenwithThreshold = booot.stren.1000[booot.stren.1000$strength >=0.5&booot.stren.1000$direction>=0.5, ]


func.dx.to.DM = function(dat, visit){
  dx = paste("SA_DX_VIS", visit, sep = "")
  print(dx)
  toNodesForDXVis =  dat[boot.strenwithThreshold$from ==dx & dat$to %in% c(grep("AR|BIT|DOT|Motor",dat$to, value = TRUE), cogDomains),]
  toNodesForDXVis = toNodesForDXVis[,c("from", "to")]
  rownames(toNodesForDXVis) = 1:nrow(toNodesForDXVis)
  dataVis = spearmanRankPadjGoodCor[spearmanRankPadjGoodCor$from%in% toNodesForDXVis$from & spearmanRankPadjGoodCor$to%in% toNodesForDXVis$to,]
  dataVis$from = gsub("SA_", "", dataVis$from)
  #dataVis$from = gsub("_VIS1", "", dataVis$from)
  dataVis$to = gsub("SA_zcode_", "", dataVis$to)
  dataVis[,6] = NULL
  # for(i in 1:nrow(dataVis)){
  #   if(dataVis[i,"p.Adj"] < 2.2e-16){
  #     dataVis[i,"p.Adj"] = "<2.2e-16"
  #   }
  #   else if(dataVis[i,"p.Adj"] <0.0001 & dataVis[i,"p.Adj"] >2.2e-16){
  #     dataVis[i,"p.Adj"] = 0.0001
  #   }
  #   else if(dataVis[i,"p.Adj"] <0.05 & dataVis[i,"p.Adj"] >0.0001){
  #     dataVis[i,"p.Adj"] = 0.05
  #   }
  # }
  return(dataVis)
}

dxtoDM_vis1 = func.dx.to.DM(boot.strenwithThreshold,1)
dxtoDM_vis6 = func.dx.to.DM(boot.strenwithThreshold,6)
dxtoDM_vis12 = func.dx.to.DM(boot.strenwithThreshold,12)
dxtoDM_vis24 = func.dx.to.DM(boot.strenwithThreshold,24)
dxtoDM_vis36 = func.dx.to.DM(boot.strenwithThreshold,36)
dxtoDM_vis36_all = rbind.data.frame(dxtoDM_vis1,dxtoDM_vis6,dxtoDM_vis12,dxtoDM_vis24,dxtoDM_vis36)



func.age.to.DM = function(dat){
  toNodesForAge =  dat[boot.strenwithThreshold$from =="SA_AGE_VIS1" & dat$to %in% c(grep("AR|BIT|DOT|Motor",dat$to, value = TRUE), cogDomains),]
  toNodesForAge = toNodesForAge[,c("from", "to")]
  dataVis = spearmanRankPadj[spearmanRankPadj$from%in% toNodesForAge$from & spearmanRankPadj$to%in% toNodesForAge$to,]
  dataVis$from = gsub("SA_", "", dataVis$from)
  dataVis$from = gsub("_VIS1", "", dataVis$from)
  dataVis$to = gsub("SA_zcode_", "", dataVis$to)
  dataVis[,6] = NULL
  # for(i in 1:nrow(dataVis)){
  # if(dataVis[i,"p.Adj"] < 2.2e-16){
  #   dataVis[i,"p.Adj"] = "<2.2e-16"
  # }
  # else if(dataVis[i,"p.Adj"] <0.0001 & dataVis[i,"p.Adj"] >2.2e-16){
  #   dataVis[i,"p.Adj"] = 0.0001
  # }
  # else if(dataVis[i,"p.Adj"] <0.05 & dataVis[i,"p.Adj"] >0.0001){
  #   dataVis[i,"p.Adj"] = 0.05
  # }
  # else if(dataVis[i,"p.Adj"] == 0){
  #     dataVis[i,"p.Adj"] = 0
  # }
  # }
  return(dataVis)
}
ageToDM = func.age.to.DM(boot.strenwithThreshold)

flattenCorrMatrix <- function(cormat, pmat) {
  ut <- upper.tri(cormat)
  data.frame(
    row = rownames(cormat)[row(cormat)[ut]],
    column = rownames(cormat)[col(cormat)[ut]],
    cor  =(cormat)[ut],
    p = pmat[ut]
  )
}
library(Hmisc)
res<-rcorr(as.matrix(data[,2:7]))
flattenCorrMatrix(res$r, res$P)
model.lm <- lm(SA_zcode_MotorTappingFeatures_VIS1  ~ SA_DX_VIS1, data = data)
save.image("~/Documents/Documents_IT/paper/ADNI_VAMBN_paper_final/scatterPlots.RData")
