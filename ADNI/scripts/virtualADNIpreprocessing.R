rm(list=ls())
virtualADNI <- read_csv("virtualADNI.csv")
virtualADNI$X1 = NULL
virtualADNI[,grep("VOL", colnames(virtualADNI), value = TRUE)] = NULL
virtualADNI[,grep("CSF", colnames(virtualADNI), value = TRUE)] = NULL
demogs = grep("SA_PT|SA_PH|SA_AGE|SA_APOE4|rs.*VIS", colnames(virtualADNI), value = TRUE)
virtualADNI <- virtualADNI %>% rename_at(demogs, 
                                         funs(str_replace(., "\\_VIS1", "")))
demogs = grep("SA_PT|SA_PH|SA_AGE|SA_APOE4|rs", colnames(virtualADNI), value = TRUE)
virtualADNI <- virtualADNI %>% rename_at(demogs, 
                                         funs(str_replace(., "SA_", "")))
names(virtualADNI)[names(virtualADNI)=="ID"] = "RID"
sel <- grep(".*\\_VIS", colnames(virtualADNI), value = TRUE)
virtualADNI <- virtualADNI %>% rename_at(sel, 
                            funs(str_replace(., "\\_VIS", "\\.")))
sel <- grep("SA", colnames(virtualADNI), value = TRUE)
virtualADNI <- virtualADNI %>% rename_at(sel, 
                                         funs(str_replace(., "SA_", "")))
sel <- grep("zcode", colnames(virtualADNI), value = TRUE)
virtualADNI <- virtualADNI %>% rename_at(sel, 
                                         funs(str_replace(., "zcode_", "")))
virtualADNI = as.data.frame(virtualADNI)
virtualADNI_long = reshape(virtualADNI,
        varying = grep(".*\\.[0-9]+$", colnames(virtualADNI)),
        direction = "long",
        #v.names = "Value",
        idvar = "RID",
        timevar = "VISCODE")
rownames(virtualADNI_long) = 1:nrow(virtualADNI_long)
modify.educat <- function(value){
  ifelse(value <= 9, value<-"Low",value<-"High")
}
baseline.adni <- virtualADNI_long[,grep("PT|PH|AGE|APOE4|AGE|rs.*", colnames(virtualADNI_long), value = TRUE)]

### Grouping of patients by progression patterns ###
aggDiagnosis <- ddply(virtualADNI_long, .(RID), summarize, DX = toString(DX))

### Filtering last MCI/first Dementia points ###
mci <- virtualADNI_long[virtualADNI_long$DX == "MCI",]
dementia <- virtualADNI_long[virtualADNI_long$DX == "Dementia",]

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
write_csv(shuffleTest,"ADNI_virtual.csv")


