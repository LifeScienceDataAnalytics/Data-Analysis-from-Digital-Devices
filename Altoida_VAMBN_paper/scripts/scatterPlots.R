
setwd("~/Altoida_VAMBN_paper")


rm(list=ls())
library(bnlearn)
library(tidyverse)
library(plyr)
library(grid)
library(gridExtra)
library(psych)
library(ggplot2)
library(corrplot)
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

##dx to digital measures##
boot_stren_50_percent <- read_csv("boot_stren_50_percent.csv")
dxCor = (boot_stren_50_percent[boot_stren_50_percent$from == "DX", "to"])
dxCor = as_tibble(dxCor)%>% 
  unlist(., use.names=FALSE)

dxCorDf = spearmanRankPadj[spearmanRankPadj$from == "SA_DX_VIS1"&   spearmanRankPadj$to%in% grep("ARObjectFinding|MotorTappingFeatures|scode_ARIntroReadTimes|Flexibility|Planning|ProspectiveMemory", spearmanRankPadj$to,value=TRUE),]
##age to digital measures
ageCor = (boot_stren_50_percent[boot_stren_50_percent$from == "age", "to"])
ageCorDf = spearmanRankPadj[spearmanRankPadj$from == "SA_age_VIS1"&   spearmanRankPadj$to%in% grep("ARScreenButtonPresses|ARGlobalTelemetryVariance|ARObjectPlacementFFT", spearmanRankPadj$to,value=TRUE),]
ageCorDf[,6] = NULL

spearmanRankPadjGoodCor = spearmanRankPadj[(spearmanRankPadj$r>0.20)|(spearmanRankPadj$r< -0.20),]
spearmanRankPadjGoodCor$lower.adj = round(spearmanRankPadjGoodCor$lower.adj,2)
spearmanRankPadjGoodCor$upper.adj = round(spearmanRankPadjGoodCor$upper.adj,2)
spearmanRankPadjGoodCor$r = round(spearmanRankPadjGoodCor$r,2)

##MMSE Orientation and AR Screen Button Presses
library("ggpubr")
r = paste(spearmanRankPadjGoodCor[spearmanRankPadjGoodCor$from=="zcode_MMSE_Orientation_VIS1"&spearmanRankPadjGoodCor$to=="zcode_ARScreenButtonPresses_VIS1","r"])
r_tx = paste("R =", r, sep = " ")
low_int= spearmanRankPadjGoodCor[spearmanRankPadjGoodCor$from=="zcode_MMSE_Orientation_VIS1"&spearmanRankPadjGoodCor$to=="zcode_ARScreenButtonPresses_VIS1","lower.adj"]
upper_int= spearmanRankPadjGoodCor[spearmanRankPadjGoodCor$from=="zcode_MMSE_Orientation_VIS1"&spearmanRankPadjGoodCor$to=="zcode_ARScreenButtonPresses_VIS1","upper.adj"]
conf_int_tx = paste("95% CI [", low_int, upper_int, "]", sep = "")
p_adj = spearmanRankPadjGoodCor[spearmanRankPadjGoodCor$from=="zcode_MMSE_Orientation_VIS1"&spearmanRankPadjGoodCor$to=="zcode_ARScreenButtonPresses_VIS1","p.Adj"]
p_adj_tx = paste("p.Adj =", "<0.0001", sep = " ")
label  = paste(r_tx, conf_int_tx, p_adj_tx, sep = ", ")

png(filename="scatterPlots/scatterOriScreenButton.png", width = 5, height = 5, units = 'in',res=300)
p1 = ggscatter(data, x = "zcode_MMSE_Orientation_VIS1", y = "zcode_ARScreenButtonPresses_VIS1", 
          add = "reg.line", conf.int = TRUE, conf.int.level = 0.95, fullrange = TRUE,
          cor.coef = FALSE,  cor.method = "spearman",
          cor.coeff.args = list(method = "spearman", label.x.npc = "center", label.y.npc = "top"),
          xlab = "MMSE_Orientation", ylab = "ARScreenButtonPresses", font.x = 14,
          font.y = 14) + 
          annotate("text", x = 4, y = 2.5, label = label,  size=5) 

dev.off()

##MMSE Language and zcode_ARPlaceAndFindTelemetryVariance_VIS1##
r1 = paste(spearmanRankPadjGoodCor[spearmanRankPadjGoodCor$from=="zcode_MMSE_Language_VIS1"&spearmanRankPadjGoodCor$to=="zcode_ARPlaceAndFindTelemetryVariance_VIS1","r"])
r1_tx = paste("R =", r1, sep = " ")
low_int1= spearmanRankPadjGoodCor[spearmanRankPadjGoodCor$from=="zcode_MMSE_Language_VIS1"&spearmanRankPadjGoodCor$to=="zcode_ARPlaceAndFindTelemetryVariance_VIS1","lower.adj"]
upper_int1= spearmanRankPadjGoodCor[spearmanRankPadjGoodCor$from=="zcode_MMSE_Language_VIS1"&spearmanRankPadjGoodCor$to=="zcode_ARPlaceAndFindTelemetryVariance_VIS1","upper.adj"]
conf_int_tx1 = paste("95% CI [", low_int1, upper_int1, "]", sep = "")
p_adj1 = spearmanRankPadjGoodCor[spearmanRankPadjGoodCor$from=="zcode_MMSE_Language_VIS1"&spearmanRankPadjGoodCor$to=="zcode_ARPlaceAndFindTelemetryVariance_VIS1","p.Adj"]
p_adj_tx1 = paste("p.Adj =", "<0.0001", sep = " ")
label1  = paste(r1_tx, conf_int_tx1, p_adj_tx1, sep = ", ")
png(filename="scatterPlots/scatterLanARPlaceandFind.png", width = 5, height = 5, units = 'in',res=300)
#tiff("scatterOriScreenButton.tif",width = 4, height = 4, units = 'in', res = 300)
p2=ggscatter(data, x = "zcode_MMSE_Language_VIS1", y = "zcode_ARPlaceAndFindTelemetryVariance_VIS1", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef.coord= c(-4,5),cor.method = "spearman",
          xlab = "MMSE_Language", ylab = "ARPlaceAndFindTelemetryVariance",font.x = 14,
          font.y = 14) + annotate("text", x = -3, y = 5, label = label1, size=5) 
dev.off()


#MMSE Language and zcode_ARObjectPlacement_VIS1##
r2 = paste(spearmanRankPadjGoodCor[spearmanRankPadjGoodCor$from=="zcode_MMSE_Language_VIS1"&spearmanRankPadjGoodCor$to=="zcode_ARObjectPlacement_VIS1","r"])
r2_tx = paste("R =", r2, sep = " ")
low_int2= spearmanRankPadjGoodCor[spearmanRankPadjGoodCor$from=="zcode_MMSE_Language_VIS1"&spearmanRankPadjGoodCor$to=="zcode_ARObjectPlacement_VIS1","lower.adj"]
upper_int2= spearmanRankPadjGoodCor[spearmanRankPadjGoodCor$from=="zcode_MMSE_Language_VIS1"&spearmanRankPadjGoodCor$to=="zcode_ARObjectPlacement_VIS1","upper.adj"]
conf_int_tx2 = paste("95% CI [", low_int2, "-",upper_int2, "]", sep = "")
p_adj2 = spearmanRankPadjGoodCor[spearmanRankPadjGoodCor$from=="zcode_MMSE_Language_VIS1"&spearmanRankPadjGoodCor$to=="zcode_ARObjectPlacement_VIS1","p.Adj"]
p_adj_tx2 = paste("p.Adj =", "<0.05", sep = " ")
label2  = paste(r2_tx, conf_int_tx2, p_adj_tx2, sep = ", ")
png(filename="scatterPlots/scatterLanARPlace.png", width = 5, height = 5, units = 'in',res=300)
#tiff("scatterOriScreenButton.tif",width = 4, height = 4, units = 'in', res = 300)
p3=ggscatter(data, x = "zcode_MMSE_Language_VIS1", y = "zcode_ARObjectPlacement_VIS1", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef.coord= c(-4,5),cor.method = "spearman",
          xlab = "MMSE_Language", ylab = "ARObjectPlacement", font.x = 14,
          font.y = 14) + annotate("text", x = -3, y = 5, label = label2, size=5) 
dev.off()


#MMSE Recall and Flexibility##
r3 = paste(spearmanRankPadjGoodCor[spearmanRankPadjGoodCor$from=="zcode_MMSE_Memory_Recall_VIS1"&spearmanRankPadjGoodCor$to=="SA_Flexibility_VIS1","r"])
r3_tx = paste("R =", r3, sep = " ")
low_int3= spearmanRankPadjGoodCor[spearmanRankPadjGoodCor$from=="zcode_MMSE_Memory_Recall_VIS1"&spearmanRankPadjGoodCor$to=="SA_Flexibility_VIS1","lower.adj"]
upper_int3= spearmanRankPadjGoodCor[spearmanRankPadjGoodCor$from=="zcode_MMSE_Memory_Recall_VIS1"&spearmanRankPadjGoodCor$to=="SA_Flexibility_VIS1","upper.adj"]
conf_int_tx3 = paste("95% CI [", low_int3, "-",upper_int3, "]", sep = "")
p_adj3 = spearmanRankPadjGoodCor[spearmanRankPadjGoodCor$from=="zcode_MMSE_Memory_Recall_VIS1"&spearmanRankPadjGoodCor$to=="SA_Flexibility_VIS1","p.Adj"]
p_adj_tx3 = paste("p.Adj =", "<0.05", sep = " ")
label3 = paste(r3_tx, conf_int_tx3, p_adj_tx3, sep = ", ")
png(filename="scatterPlots/scatterRecallFlex.png", width = 5, height = 5, units = 'in',res=300)
#tiff("scatterOriScreenButton.tif",width = 4, height = 4, units = 'in', res = 300)
p4=ggscatter(data, x = "zcode_MMSE_Memory_Recall_VIS1", y = "SA_Flexibility_VIS1", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef.coord= c(-4,5),cor.method = "spearman",
          xlab = "MMSE_Memory_Recall_VIS1", ylab = "Flexibility", font.x = 14,
          font.y = 14) + annotate("text", x = 1.5, y = 85, label = label3,  size=5) 
dev.off()
library(cowplot)
png("scatterPlots/scatPlotAll.png", width = 300, height = 200, units='mm',res = 300) 
plot_grid(p1, p2, p3, p4, align = "h", label_size = 12, rel_widths = c(1,1), ncol = 2, nrow = 2)
dev.off()


#MMSE attention concentraation and AR place and find telemetry variance##
r4 = paste(spearmanRankPadj[spearmanRankPadj$from=="zcode_MMSE_Attention_Concentration_VIS1"&spearmanRankPadj$to=="zcode_ARPlaceAndFindTelemetryVariance_VIS1","r"])
r4_tx = paste("R =", r4, sep = " ")
low_int4= spearmanRankPadjGoodCor[spearmanRankPadjGoodCor$from=="zcode_MMSE_Attention_Concentration_VIS1"&spearmanRankPadjGoodCor$to=="zcode_ARPlaceAndFindTelemetryVariance_VIS1","lower.adj"]
upper_int4= spearmanRankPadjGoodCor[spearmanRankPadjGoodCor$from=="zcode_MMSE_Attention_Concentration_VIS1"&spearmanRankPadjGoodCor$to=="zcode_ARPlaceAndFindTelemetryVariance_VIS1","upper.adj"]
conf_int_tx4 = paste("95% CI [", low_int4, "-",upper_int4, "]", sep = "")
p_adj4 = spearmanRankPadjGoodCor[spearmanRankPadjGoodCor$from=="zcode_MMSE_Attention_Concentration_VIS1"&spearmanRankPadjGoodCor$to=="zcode_ARPlaceAndFindTelemetryVariance_VIS1","p.Adj"]
p_adj_tx4 = paste("p.Adj =", "<0.05", sep = " ")
label4 = paste(r4_tx, conf_int_tx4, p_adj_tx4, sep = ", ")
png(filename="scatterPlots/scatterMMSEAttConcARTelVar.png", width = 5, height = 5, units = 'in',res=300)
#tiff("scatterOriScreenButton.tif",width = 4, height = 4, units = 'in', res = 300)
ggscatter(data, x = "zcode_MMSE_Attention_Concentration_VIS1", y = "zcode_ARPlaceAndFindTelemetryVariance_VIS1", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef.coord= c(-4,5),cor.method = "spearman",
          xlab = "MMSE_Attention_Concentration_VIS1", ylab = "ARPlaceAndFindTelemetryVariance_VIS1", font.x = 14,
          font.y = 14) + annotate("text", x = 1.5, y = 85, label = label3,  size=5) 
dev.off()
library(cowplot)
png("scatterPlots/scatPlotAll.png", width = 300, height = 200, units='mm',res = 300) 
plot_grid(p1, p2, p3, p4, align = "h", label_size = 12, rel_widths = c(1,1), ncol = 2, nrow = 2)
dev.off()

##MMSE Language and zcode_ARObjectPlacement##
res5 <-cor.test(data$zcode_MMSE_Working_Memory_Registration_VIS1, data$zcode_ARObjectFinding_VIS1,  method = "spearman")
res5
png(filename="scatterWorkMemARObjFing.png", width = 4, height = 4, units = 'in',res=300)
#tiff("scatterOriScreenButton.tif",width = 4, height = 4, units = 'in', res = 300)
ggscatter(data, x = "zcode_MMSE_Working_Memory_Registration_VIS1", y = "zcode_ARObjectFinding_VIS1", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE,cor.coef.coord= c(1.2,4),cor.method = "spearman",
          xlab = "MMSE_Working_Memory_Registration", ylab = "ARObjectFinding")
dev.off()

