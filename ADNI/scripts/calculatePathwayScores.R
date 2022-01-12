combined_ADNI_1_csv <- read_csv("~/Documents/Documents – IT-Admin’s MacBook Pro/snpWork/processed_ADNI_1_GWAS/combined_ADNI_1_csv.csv")
colnames(combined_ADNI_1_csv )[colnames(combined_ADNI_1_csv )=='Unnamed: 0']<-'SUBJID'
combined_ADNI_1_csv = combined_ADNI_1_csv[!grepl("Unnamed", colnames(combined_ADNI_1_csv))]

combined_ADNIGO_2_csv <- read_csv("~/Documents/Documents – IT-Admin’s MacBook Pro/snpWork/processed_ADNI_GO_2_GWAS/combined_ADNIGO_2_csv.csv")
colnames(combined_ADNIGO_2_csv )[colnames(combined_ADNIGO_2_csv )=='Unnamed: 0']<-'SUBJID'
combined_ADNIGO_2_csv = combined_ADNIGO_2_csv[!grepl("Unnamed", colnames(combined_ADNIGO_2_csv))]

setdiff(colnames(combined_ADNI_1_csv), colnames(combined_ADNIGO_2_csv))
common <- intersect(names(combined_ADNI_1_csv), names(combined_ADNIGO_2_csv))
allSNPs <- rbind(combined_ADNI_1_csv[,common], combined_ADNIGO_2_csv[,common])
SNPsToGenesToMechansims_mapping_final <- read_csv("~/Documents/Documents – IT-Admin’s MacBook Pro/snpWork/SNPsToGenesToMechansims_mapping_final.csv")
relevantPathwaysADTop = read.csv("/Users/msood/Documents/Documents – IT-Admin’s MacBook Pro/snpWork/relevantPathwaysADTop.csv")
relevantPathwaysADTop = as.data.frame(relevantPathwaysADTop[1:10,])
names(relevantPathwaysADTop)[1] = "Pathways"
SNPsToGenesToMechansims_mapping_final_top10 <- SNPsToGenesToMechansims_mapping_final[SNPsToGenesToMechansims_mapping_final$`Subgraph Name` %in% relevantPathwaysADTop$Pathways,]
pathways = unique(SNPsToGenesToMechansims_mapping_final_top10$`Subgraph Name`)

dataSnpAll = list()
for(subGraph in pathways){
  #print(subGraph)
  snps = SNPsToGenesToMechansims_mapping_final_top10[SNPsToGenesToMechansims_mapping_final_top10$`Subgraph Name`
                                                     == subGraph,]
  snps = snps$rsID
  snpForPath = intersect(colnames(allSNPs), snps)
  snpDf = allSNPs[,snpForPath]
  dfName = paste(subGraph, "_VIS1", sep = "")
  print(dfName)
  print(snpDf)
  dataSnpAll[[dfName]] = snpDf
}




write.csv(allSNPs, "/Users/msood/Documents/Documents – IT-Admin’s MacBook Pro/snpWork/allSNPsAndBurdenScores.csv")









