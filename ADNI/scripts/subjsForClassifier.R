rm(list=ls())
library(dplyr)
library(dominanceanalysis)
library(VGAM)
library(caret)
library(readr)
library(plyr)
name<-'main'
data_out<-paste0('data/data_out/',name)
#virtual_m<-readRDS(paste0(data_out,'_VirtualPPts.rds'))
#virtual_nm <-  readRDS(paste0(data_out,'_VirtualPPts_nm.rds'))
#virtual = merge(virtual_m, virtual_nm)
data <- readRDS(paste0(data_out,'_RealPPts.rds'))
#data <- data %>% select(-contains("VIS36"))
#colnames(data)<-gsub('_VIS1_','_',colnames(data))
decRP<-read_csv('data/HI-VAE/reconRP.csv')
decRP_nm<-read_csv('data/HI-VAE/reconRP_nm.csv')
decRP = merge(decRP, decRP_nm)
#decRP <- decRP %>% select(-contains("VIS36"))
decVP<-read.csv('data/HI-VAE/decodedVP.csv')
decVP_nm<-read.csv('data/HI-VAE/decodedVP_nm.csv')
decVP = merge(decVP, decVP_nm)
#decVP <- decVP %>% select(-contains("VIS36"))
#saveVirtual = cbind(virtual[,grep("SA", colnames(virtual), value= TRUE)], decVP)
#saveVirtual = saveVirtual[,!grepl("AUX",colnames(saveVirtual)),]
#saveVirtual$SUBJID = NULL
#saveVirtual$ID <- seq.int(nrow(saveVirtual))


#data = saveReal
data = data_Dig
func.processing.for.classifier = function(data,first_stage){
        data$X1 = NULL
        data[,grep("VOL|ICV|FDG|AV45|\\_VIS36", colnames(data), value = TRUE)] = NULL
        data[,grep("CSF", colnames(data), value = TRUE)] = NULL
        data[,grep("rs", colnames(data), value = TRUE)] = NULL
        demogs = grep("SA_PT|SA_PH|SA_AGE|SA_APOE4|Amyloid", colnames(data), value = TRUE)
  
        data <- data %>% rename_at(demogs,
                                         funs(str_replace(., "\\_VIS1", "")))
        demogs = grep("SA_PT|SA_PH|SA_AGE|SA_APOE4|SA_Amyloid", colnames(data), value = TRUE)
        data <- data %>% rename_at(demogs, 
                                                 funs(str_replace(., "SA_", "")))
        names(data)[names(data)=="SUBJID"] = "RID"
        sel <- grep(".*\\_VIS", colnames(data), value = TRUE)
        data <- data %>% rename_at(sel,
                                                 funs(str_replace(., "\\_VIS", "\\.")))
        sel <- grep("SA", colnames(data), value = TRUE)
        data <- data %>% rename_at(sel, 
                                                 funs(str_replace(., "SA_", "")))
        sel <- grep("zcode", colnames(data), value = TRUE)
        data <- data %>% rename_at(sel, 
                                                 funs(str_replace(., "zcode_", "")))
        data = as.data.frame(data)
        data_long = reshape(data,varying = grep(".*\\.[0-9]+$", colnames(data), value = TRUE),
                                   direction = "long",
                                   #v.names = "Value",
                                   idvar = "RID",
                                  timevar = "VISCODE")
        aggDiagnosis <- ddply(data_long, .(RID), summarise, DX = toString(DX))
        patientidConv = c()
        for(i in 1:nrow(aggDiagnosis)){
          if(grepl("^MCI.*Dementia.*MCI$", aggDiagnosis[i,2])){
            patientidConv <- c(patientidConv, aggDiagnosis[i,1])
          }
        }
        
        patientidConvCNtoMCI = c()
        for(i in 1:nrow(aggDiagnosis)){
          if(grepl("^CN.*MCI.*CN$", aggDiagnosis[i,2])){
            patientidConvCNtoMCI <- c(patientidConvCNtoMCI, aggDiagnosis[i,1])
          }
        }
        
        patientidConvMCItoCN = c()
        for(i in 1:nrow(aggDiagnosis)){
          if(grepl("^MCI.*CN.*MCI$", aggDiagnosis[i,2])){
            patientidConvMCItoCN <- c(patientidConvMCItoCN, aggDiagnosis[i,1])
          }
        }
        
        patientidConvDementiaToMCI = c()
        for(i in 1:nrow(aggDiagnosis)){
          if(grepl("^Dementia.*MCI.*Dementia$", aggDiagnosis[i,2])){
            patientidConvDementiaToMCI <- c(patientidConvDementiaToMCI, aggDiagnosis[i,1])
          }
        }
        
        patientidConv = c(patientidConv, patientidConvCNtoMCI, patientidConvMCItoCN, patientidConvDementiaToMCI)
        aggDiagnosis <- aggDiagnosis[which(!aggDiagnosis$RID %in% patientidConv),]
        
        nlSubjIDS = c()
        for(i in 1:nrow(aggDiagnosis)){
          if(grepl("^CN.*CN$", aggDiagnosis[i,2])){
            nlSubjIDS <- c(nlSubjIDS, aggDiagnosis[i,1])
          }
        }
        mciSubjIDS = c()
        for(i in 1:nrow(aggDiagnosis)){
          if(grepl("^MCI.*MCI$", aggDiagnosis[i,2])){
            mciSubjIDS <- c(mciSubjIDS, aggDiagnosis[i,1])
          }
        }
        
        dementiaSubjIDS = c()
        for(i in 1:nrow(aggDiagnosis)){
          if(grepl("^Dementia.*Dementia$", aggDiagnosis[i,2])){
            dementiaSubjIDS <- c(dementiaSubjIDS, aggDiagnosis[i,1])
          }
        }
        ##subset for patients remaining healthy and patients remaining mci
        if(first_stage == "CN"){
          data_long = data_long[data_long$RID %in% c(nlSubjIDS,mciSubjIDS),]
        }
        if(first_stage == "MCI"){
          data_long = data_long[data_long$RID %in% c(mciSubjIDS, dementiaSubjIDS),]
        }
        data_long$DX = factor(data_long$DX)
        rownames(data_long) = 1:nrow(data_long)
        data_long$VISCODE[data_long$VISCODE==1] <- 0
        data_long$ridViscode = paste(data_long$RID, data_long$VISCODE, sep = ",")
        rownames(data_long) = 1:nrow(data_long)
        # # ##mohamed's data#
        # if(originalClassifier == TRUE){
        #   final_data_classifier = read_csv("ADNI_shuffle_08_10_2020.csv")
        #   final_data_classifier$ridViscode = paste(final_data_classifier$RID, final_data_classifier$VISCODE, sep = ",")
        #   adniDataForClassifier = data_long[data_long$ridViscode %in% final_data_classifier$ridViscode,]
        # }
        #else{
          adniDataForClassifier = data_long
        #}
        ids = table(adniDataForClassifier$RID)
        idsNotRepeating = c()
        for(i in 1:length(ids)){
        if(ids[i][[1]] == 1){
            idsNotRepeating = c(idsNotRepeating, as.integer(names(ids[i])))
          }
        }
      adniDataForClassifier = adniDataForClassifier[!adniDataForClassifier$RID %in% idsNotRepeating,]
      # modify.educat <- function(value){
      # ifelse(value <= 9, value<-"Low",value<-"High")
      #   }
      # adniDataForClassifier$DX <- factor(adniDataForClassifier$DX,
      #                                levels = c("Dementia","MCI"),
      #                                labels = c(0,1))
      if(first_stage == "CN"){
      adniDataForClassifier$DX <- factor(adniDataForClassifier$DX,
                                     levels = c("CN","MCI"),
                                     labels = c(0,1))}
      if(first_stage == "MCI"){
        adniDataForClassifier$DX <- factor(adniDataForClassifier$DX,
                                           levels = c("MCI", "Dementia"),
                                           labels = c(0,1))}
      adniDataForClassifier$PTGENDER <- factor(adniDataForClassifier$PTGENDER,
                                           levels = c("Female","Male"),
                                           labels = c(0,1))
      adniDataForClassifier$PTGENDER <- as.factor(adniDataForClassifier$PTGENDER)
      adniDataForClassifier$PTETHCAT <- factor(adniDataForClassifier$PTETHCAT,
                                           levels = c("Hisp/Latino","Not Hisp/Latino","Unknown"),
                                           labels = c(0,1,2))
      adniDataForClassifier$PTETHCAT <- as.factor(adniDataForClassifier$PTETHCAT)
      adniDataForClassifier$PTMARRY <- factor(adniDataForClassifier$PTMARRY,
                                          levels = c("Divorced","Married","Never married","Widowed", "Unknown"),
                                          labels = c(0,1,2,3,4))
      adniDataForClassifier$PTEDUCAT = as.factor(adniDataForClassifier$PTEDUCAT)
      adniDataForClassifier$AGE = round(adniDataForClassifier$AGE, 1)
      adniDataForClassifier$PTMARRY <- as.factor(adniDataForClassifier$PTMARRY)
      adniDataForClassifier$ridViscode = NULL
      
      rownames(adniDataForClassifier) = 1:nrow(adniDataForClassifier)

      return(adniDataForClassifier)
}

saveReal = merge(data[,c(grep("SA", colnames(data), value= TRUE), "SUBJID")], decRP)
saveReal = saveReal[,!grepl("AUX",colnames(saveReal)),]
saveReal$RID = as.integer(sub('.*\\_', '', saveReal$SUBJID))
saveReal$SUBJID = NULL
saveReal$ID <- seq.int(nrow(saveReal))

realDataClassifierCNMCI = func.processing.for.classifier(saveReal, first_stage = "CN")
realDataClassifierMCIDementia = func.processing.for.classifier(saveReal, first_stage = "MCI")
realDataClassifierCNMCIBL = realDataClassifierCNMCI[realDataClassifierCNMCI$VISCODE == 0,]
realDataClassifierMCIDementiaBL = realDataClassifierMCIDementia[realDataClassifierMCIDementia$VISCODE == 0,]

write.csv(realDataClassifierCNMCI, "real_decoded_adni_data_classifier_all.csv")
write.csv(realDataClassifierCNMCIBL, "real_decoded_adni_data_classifier_cn_mci_bl.csv")


  
write.csv(realDataClassifierMCIDementia, "real_decoded_adni_data_classifier_mci_dementia.csv")
write.csv(realDataClassifierMCIDementiaBL, "real_decoded_adni_data_classifier_mci_dementia_bl.csv")
#virtualDataClassifier = func.processing.for.clssifier(saveVirtual, originalClassifier=FALSE)
#write.csv(virtualDataClassifier, "virtual_decoded_adni_data_classifier.csv")



####IDSN virtual data##
# idsnData = saveVirtual[,grep("PT|FAQ|MMSE|VOL|AGE|APOE4", colnames(saveVirtual), value = TRUE)]
# sel <- grep("SA", colnames(idsnData), value = TRUE)
# idsnData <- idsnData %>% rename_at(sel, 
#                            funs(str_replace(., "SA_", "")))
# demogs = grep("PT|PH|AGE|APOE4|rs.*VIS", colnames(idsnData), value = TRUE)
# idsnData <- idsnData %>% rename_at(demogs,
#                            funs(str_replace(., "\\_VIS1", "")))
# idsnData$PTEDUCAT = round(idsnData$PTEDUCAT, 1)
# idsnData$AGE = round(idsnData$AGE, 1)
# idsnData = idsnData[1:1000,]
# write.csv(idsnData, "idsnData_1000.csv")
