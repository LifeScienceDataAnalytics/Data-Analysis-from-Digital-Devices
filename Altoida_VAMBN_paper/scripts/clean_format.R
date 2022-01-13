setwd("~/Altoida_VAMBN_paper")
rm(list=ls())
library(openxlsx)
library(rlist)
library(plyr)
library(stringr)
library(Hmisc)
library(readr)
##library to obtain longest common substring from multiple strings
library(PTXQC)
library(naniar)
library(rle)
print(getwd())


##load the altoida data 
##load the file for altoida data
altoidaDf = read.csv("ALTOIDA/.....csv")

#remove the not required features, Visuospatial just had 1 value and increase the total score of MMSE more than 30
altoidaDf$MMSE_Visuospatial = NULL
##MMSE_Working_Memory_Registration can have max value as 3, remove the patient having value 9
altoidaDf = altoidaDf[!altoidaDf$MMSE_Working_Memory_Registration == 9,] ##Patient ID 009-20_047

##change the data types of required variables
altoidaDf$age = as.numeric(altoidaDf$age)
altoidaDf$gender = as.factor(altoidaDf$gender)
altoidaDf$unique_subject_id = as.character(as.factor(altoidaDf$unique_subject_id))

##only take preferred session 1 for the data as mentioned by Max
altoidaDf = altoidaDf[altoidaDf$PreferredSession == 1,]
id.rle = rle(altoidaDf$unique_subject_id)
###Change December 2021 to generate data fro Grid search (comment) and change it back later###
altoidaDf$unique_subject_id <- paste0(rep(id.rle$values, times = id.rle$lengths), "_",
                                      unlist(lapply(id.rle$lengths, seq_len)))
##Change December 2021###
# write out all the offending strings, replacing all -1 with NA's 
na_strings <- c("-1", "-1 ", "-1.0000", "-1.0000 ")
altoidaDf = altoidaDf %>% replace_with_na_all(condition = ~.x  %in% na_strings)

##removing columns that only contain NA values
altoidaDf = altoidaDf[, colSums(is.na(altoidaDf)) != nrow(altoidaDf)]

altoidaDf = altoidaDf[which(rowMeans(!is.na(altoidaDf)) > 0.5), which(colMeans(!is.na(altoidaDf)) > 0.5)]

##testing the missing data for all these features
speechPressureEye = altoidaDf[,grep("Speech|Eye|pressure", value = TRUE, colnames(altoidaDf))]
##testing if the data is more than 50% missing for all these columns##
speechPressureEye = speechPressureEye[which(rowMeans(!is.na(speechPressureEye)) > 0.5), which(colMeans(!is.na(speechPressureEye)) > 0.5)]

##also delete BITDOT features as they are not of much significance, and other irrelevant features,
##remove speech and eye as it has just few values
excludeColNames = c(grep("Speech|Eye", colnames(altoidaDf), value = TRUE))
excludeColNames = c(excludeColNames, "PreferredSession")

##change the name of the ID column to fit to the VAMBN model 
names(altoidaDf)[which(names(altoidaDf) == "unique_subject_id")] = "SUBJID"

##remove the features that you are not sure about##
altoidaDf = altoidaDf[,which(!colnames(altoidaDf)%in% excludeColNames)]

featureDesc =read_delim("~/Documents/Documents_IT/VAMBN_Version5_mmseEncoded_final/ALTOIDA/fraunhofer_v5_technical_grouping.csv", ";", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)
featuresRemovedFromAltoida = setdiff(featureDesc$X1,colnames(altoidaDf))
featureDesc = featureDesc[which(!featureDesc$X1 %in% featuresRemovedFromAltoida),]
names(featureDesc)[names(featureDesc) == "X1"] = "Features"
names(featureDesc)[names(featureDesc) == "X2"] = "GroupOrNot"

##summary statistics table##
#install.packages("gtsummary")
library(gtsummary)

##change the data type for ones that is required
##add a new diagnostics column and remove MCI,Healthy and AD column , better for the bayesian network 
altoidaDf["DX"] = ifelse(altoidaDf$MCI==1&altoidaDf$Amyloid == 0, "MCI",
                   ifelse(altoidaDf$MCI==1&altoidaDf$Amyloid == 1, "MCIatRiskforAD",
                        ifelse(altoidaDf$Healthy==1, "CN",
                               ifelse(altoidaDf$AD==1, "Dementia",
                                 ifelse(altoidaDf$ProAD==1, "ProAD",
                                      NA)))))

altoidaDfSummary = altoidaDf[grep("\\_1", altoidaDf$SUBJID),]
write.csv(altoidaDfSummary, "altoidaDfSummary.csv")

altoidaDf$DX[altoidaDf$DX == "CN"] = 0
altoidaDf$DX[altoidaDf$DX == "MCI"] = 1
altoidaDf$DX[altoidaDf$DX == "MCIatRiskforAD"] = 2
altoidaDf$DX[altoidaDf$DX == "ProAD"] = 3
altoidaDf$DX[altoidaDf$DX == "Dementia"] = 4

altoidaDf$MCI = NULL
altoidaDf$ProAD = NULL
altoidaDf$AD = NULL
altoidaDf$Healthy = NULL
altoidaDf$ADRisk = NULL

altoidaDf$DX = as.factor(altoidaDf$DX)
altoidaDf$Amyloid = as.factor(altoidaDf$Amyloid)
altoidaDf$gender = as.factor(altoidaDf$gender)
altoidaDf$dominantHand = as.factor(altoidaDf$dominantHand)

altoidaDf$SUBJID = as.factor(altoidaDf$SUBJID)
write.csv(altoidaDf, "altoidaDf.csv")
pt<-factor(altoidaDf$SUBJID)
######## Variable groups, load the variable groups to be autoencoded in each module 
groups <- ddply(featureDesc, .(GroupOrNot), summarise, Features = toString(Features))
##don't need MMSE total score,  delete it 
write.csv(altoidaDf, "altoidaDf.csv")
altoidaDf$MMSE = NULL 

mmseCols = as.vector(grep("MMSE", colnames(altoidaDf), value = TRUE))
mmseGroupDf = data.frame()
mmseGroupDf[1:5,"GroupOrNot"] = mmseCols
mmseGroupDf[1:5,"Features"] = mmseCols
groups = rbind.data.frame(groups, mmseGroupDf)

groupsList = c()
for(i in 1:nrow(groups)){
  #print(i)
  features =  str_trim(as.vector(unlist(strsplit(groups[i, "Features"], ","))))
  #print(features)
  groupsList[[i]] = altoidaDf[,c(features)]
  groupsList[[i]]['SUBJID']  = pt
}
names(groupsList) = groups$GroupOrNot
##featuresNotAdded yet
additionalFeatures = setdiff(colnames(altoidaDf), featureDesc$Features)
additionalDXFeatures = c("DX", "Amyloid" )
groupsList[["stalone_demog_dx"]] = cbind.data.frame(groupsList[["stalone_demog_dx"]], altoidaDf[,c(additionalDXFeatures)])
mmseStalone = grep("MMSE", additionalFeatures, value = TRUE)
cogDomainStalone = setdiff(additionalFeatures, c(additionalDXFeatures, mmseStalone))
groupsList[["stalone_cog"]] = altoidaDf[,c(cogDomainStalone)] 

saveRDS(groupsList,"data/data_condensed.rds")


