rm(list=ls())
library(openxlsx)
library(caTools)
library(dplyr)
setwd("~/snpADNIVAMBN/data_paths")
relevantPathwaysADTop = read.csv("/snpWork/relevantPathwaysADTop.csv")

list.files(pattern=".csv$") # use the pattern argument to define a common pattern  for import files with regex. Here: .csv

# create a list from these files
list.filenames<-list.files(pattern=".csv$")
#list.filenames = gsub(" ", "_", list.filenames)

# create an empty list that will serve as a container to receive the incoming files
list.data<-list()

# create a loop to read in your data
for (i in 1:length(list.filenames))
{
  list.data[[i]]<-read.csv(list.filenames[i])
}

list.filenames<- gsub(".csv", "", list.filenames)
# add the names of your data to the list
names(list.data)<- gsub(" ", "_", list.filenames)
relevantPathwaysADTop$Pathways = gsub(" ", "_",relevantPathwaysADTop$Pathways)
data_snp = list.data[relevantPathwaysADTop$Pathways]
toBeIncludedAsStalone= c()
for(i in 1:length(data_snp)){
  data_snp[[i]]["X"] = NULL
  if(ncol(data_snp[[i]]) == 2){
    print(ncol(data_snp[[i]]))
    toBeIncludedAsStalone = c(toBeIncludedAsStalone,names(data_snp[i]))
  }
}
print(toBeIncludedAsStalone)




for(i in 1:length(data_snp)){
  ##adding specific pathway to each snp file##
  data_snp[[i]] = data_snp[[i]] %>% rename_at(vars(-SUBJID), ~ paste0(names(data_snp[i]),"_",.))
}

snpsToBeIncludedAsStalone = data_snp[toBeIncludedAsStalone]

for(sugrph in toBeIncludedAsStalone){
  print(sugrph)
  data_snp[[sugrph]] = NULL
}
load("~/ADNI_VAMBN_paper/ADNIVAMBN.RData")

rm(list=setdiff(ls(), c("data_snp", "visitData", "snpsToBeIncludedAsStalone")))
names(data_snp) <- paste(names(data_snp), "VIS1", sep = "_")

visitData_complete = visitData
#visitData_complete = merge(visitData_complete, digDataAll)
require(tidyverse)
data_snp =data_snp  %>% reduce(merge, by = 'SUBJID')
data_snp = data_snp[data_snp$SUBJID%in%visitData_complete$PTID,]
names(visitData_complete)[names(visitData_complete) == "PTID"] = "SUBJID"
visitData_complete = merge(visitData_complete, data_snp)

pt<-factor(visitData_complete$SUBJID)

##IDS taked for digital data##

colnames(visitData_complete)<-gsub('.bl','_VIS1',colnames(visitData_complete))
colnames(visitData_complete)<-gsub('.m06','_VIS6',colnames(visitData_complete))
colnames(visitData_complete)<-gsub('.m12','_VIS12',colnames(visitData_complete))
colnames(visitData_complete)<-gsub('.m18','_VIS18',colnames(visitData_complete))
colnames(visitData_complete)<-gsub('.m24','_VIS24',colnames(visitData_complete))
colnames(visitData_complete)<-gsub('.m36','_VIS36',colnames(visitData_complete))
colnames(visitData_complete)<-gsub('.m48','_VIS48',colnames(visitData_complete))

######## Variable groups
csf <- visitData_complete[,grep("ABETA|TAU|PTAU", colnames(visitData_complete), value = TRUE)]
colnames(csf) = gsub("csf", "CSF", colnames(csf))
csf$SUBJID<-pt
csf_VIS1<-csf[,grepl('SUBJID|_VIS1$',colnames(csf))]


volume <- visitData_complete[,grep("Entorhinal|Fusiform|MidTemp|Hippocampus|Ventricles|WholeBrain", colnames(visitData_complete), value = TRUE)]
colnames(volume) = gsub("brain", "VOL", colnames(volume))
volume$SUBJID<-pt
volume_VIS1<-volume[,grepl('SUBJID|_VIS1$',colnames(volume))]
volume_VIS6<-volume[,grepl('SUBJID|_VIS6',colnames(volume))]
volume_VIS12<-volume[,grepl('SUBJID|_VIS12',colnames(volume))]
volume_VIS24<-volume[,grepl('SUBJID|_VIS24',colnames(volume))]



imagingPET <- visitData_complete[,grep("FDG|AV45", colnames(visitData_complete), value = TRUE)]
imagingPET$SUBJID<-pt
imagingPET_VIS1<-imagingPET[,grepl('SUBJID|_VIS1$',colnames(imagingPET))]



## Standalone
snpsToBeIncludedAsStalone = reduce(snpsToBeIncludedAsStalone, full_join, by = "SUBJID")
snpsToBeIncludedAsStalone = snpsToBeIncludedAsStalone %>% rename_at(vars(-SUBJID), ~ paste0(.,"_VIS1"))
##add abeta column in the data frame based on csf abeta present or not
#Amyloid positivity was defined as Florbetapir PET (AV-45) SUVRs > 1.11 and amyloid negativity as SUVRs ≤1.11 
visitData_complete$Amyloid_VIS1 <- ifelse(visitData_complete$AV45_VIS1 > 1.11, 1,0)
visitData_complete$Amyloid_VIS1[is.na(visitData_complete$Amyloid_VIS1)] <- 0
visitData_complete$APOE4_VIS1[is.na(visitData_complete$APOE4_VIS1)] <- 0
demogs <- visitData_complete[,grep("SUBJID|PT[^AU]|APOE4|AGE|PHS|Amyloid_VIS1", colnames(visitData_complete), value = TRUE)]
dx <- visitData_complete[,grep("SUBJID|DX", colnames(visitData_complete), value = TRUE)]
mmseOrig <- visitData_complete[,grep("SUBJID|^MMSE\\_", colnames(visitData_complete), value = TRUE)]
mmseOrig[,grep("MMSE", colnames(mmseOrig), value = TRUE)] = lapply(mmseOrig[,grep("MMSE", colnames(mmseOrig), value = TRUE)], as.integer)

########################################################################################
faq <- visitData_complete[,grep("SUBJID|^FAQ", colnames(visitData_complete), value = TRUE)]
faq[,grep("FAQ", colnames(faq), value = TRUE)] = lapply(faq[,grep("FAQ", colnames(faq), value = TRUE)], as.integer)


stalone<-visitData_complete[,c(colnames(demogs), colnames(dx), colnames(faq))]
stalone<-stalone[,grepl('_VIS|SUBJID$',colnames(stalone))]
#stalone = stalone[!grepl("FAQTOTAL.*", colnames(stalone))]
stalone$SUBJID<-pt
stalone = merge(stalone,snpsToBeIncludedAsStalone)
stalone$SUBJID = NULL
colnames(stalone)<-paste0('SA_',colnames(stalone))
stalone["SUBJID"] = pt
stalone$SA_APOE4_VIS1<-factor(stalone$SA_APOE4_VIS1)
stalone$SA_Amyloid_VIS1<-factor(stalone$SA_Amyloid_VIS1)
stalone$SA_PTMARRY_VIS1  <-factor(stalone$SA_PTMARRY_VIS1)
stalone$SA_PTGENDER_VIS1  <-factor(stalone$SA_PTGENDER_VIS1)
#stalone$SA_PTEDUCAT_VIS1<-factor(ifelse(stalone$SA_PTEDUCAT_VIS1<16,'low','high'))


stalone_VIS1_demog<-stalone[,grepl('SUBJID|PT|APOE|PH|AGE|Amyloid\\_VIS1$',colnames(stalone))]
stalone_VIS1_dx<-stalone[,grepl('SUBJID|DX\\_VIS1$',colnames(stalone))]
stalone_baseline = merge(stalone_VIS1_demog, stalone_VIS1_dx)
saveRDS(stalone_baseline, file = paste0("~/Documents/Documents_IT/paper/ADNI_VAMBN_paper/data/stalone_VIS1.rds"))

stalone_VIS1_faq<-stalone[,grepl('SUBJID|FAQ.*\\_VIS1$',colnames(stalone))]
stalone_VIS1_rs<-stalone[,grepl('SUBJID|.*rs.*\\_VIS1$',colnames(stalone))]
stalone_VIS6_dx<-stalone[,grepl('SUBJID|DX\\_VIS6$',colnames(stalone))]
stalone_VIS6_faq<-stalone[,grepl('SUBJID|FAQ.*\\_VIS6$',colnames(stalone))]
stalone_VIS12_dx<-stalone[,grepl('SUBJID|DX\\_VIS12$',colnames(stalone))]
stalone_VIS12_faq<-stalone[,grepl('SUBJID|FAQ.*\\_VIS12$',colnames(stalone))]
stalone_VIS24_dx<-stalone[,grepl('SUBJID|DX\\_VIS24$',colnames(stalone))]
stalone_VIS24_faq<-stalone[,grepl('SUBJID|FAQ.*\\_VIS24$',colnames(stalone))]
stalone_VIS36_dx<-stalone[,grepl('SUBJID|DX\\_VIS36$',colnames(stalone))]

##for mmse individual encoding##
mmse_VIS1 = mmseOrig[,grepl('SUBJID|_VIS1$',colnames(mmseOrig))]
mmse_VIS6 = mmseOrig[,grepl('SUBJID|_VIS6$',colnames(mmseOrig))]
mmse_VIS12 = mmseOrig[,grepl('SUBJID|_VIS12$',colnames(mmseOrig))]
mmse_VIS24 = mmseOrig[,grepl('SUBJID|_VIS24$',colnames(mmseOrig))]
mmse_VIS36 = mmseOrig[,grepl('SUBJID|_VIS36$',colnames(mmseOrig))]

func.mmse.encoding<-function(data){
  mmseList = c()
  mmseList$SUBJID = NULL
  listNames = c()
  for(i in 2:ncol(data)){
    print(i)
    mmseList[[i-1]] = data[i]
    data$SUBJID = factor(data$SUBJID)
    mmseList[[i-1]]["SUBJID"] = data$SUBJID
    listNames = c(listNames, names(data[i]))
  }
  names(mmseList) = listNames
  return (mmseList)
}
mmseList1 = func.mmse.encoding(mmse_VIS1)
mmseList2 = func.mmse.encoding(mmse_VIS6)
mmseList3 = func.mmse.encoding(mmse_VIS12)
mmseList4 = func.mmse.encoding(mmse_VIS24)
mmseList5 = func.mmse.encoding(mmse_VIS36)

snps_list = visitData_complete[,grep("SUBJID|\\_rs[0-9]+", colnames(visitData_complete), value=TRUE)]
library(ggdendro)
library(ade4)
library(PTXQC)

snps_list$SUBJID = NULL
snps_list = split.default(snps_list, sub("\\_rs.*", "", names(snps_list)))
for(i in 1:length(snps_list)){
  snps_list[[i]]["SUBJID"] = pt
}
names(snps_list) <- paste(names(snps_list), "VIS1", sep = "_")


data_all<-list('csf_VIS1'=csf_VIS1,
               'imagingPET_VIS1' = imagingPET_VIS1,
               'volume_VIS1'=volume_VIS1,
               'volume_VIS6'=volume_VIS6,
               'volume_VIS12'=volume_VIS12,
               'volume_VIS24'=volume_VIS24,
               'stalone_VIS1_demog'=stalone_VIS1_demog,
               'stalone_VIS1_dx'=stalone_VIS1_dx,
               'stalone_VIS1_faq'=stalone_VIS1_faq,
               'stalone_VIS1_rs'=stalone_VIS1_rs,
               'stalone_VIS6_dx'=stalone_VIS6_dx,
               'stalone_VIS6_faq'=stalone_VIS6_faq,
               'stalone_VIS12_dx'=stalone_VIS12_dx,
               'stalone_VIS12_faq'=stalone_VIS12_faq,
               'stalone_VIS24_dx'=stalone_VIS24_dx,
               'stalone_VIS24_faq'=stalone_VIS24_faq,
               'stalone_VIS36_dx'=stalone_VIS36_dx
)

data_all_with_snps_mmse = c(data_all, snps_list,mmseList1,mmseList2,mmseList3,
                            mmseList4, mmseList5)


saveRDS(data_all_with_snps_mmse, file = paste0("~/ADNI_VAMBN_paper/data/data_condensed.rds"))

