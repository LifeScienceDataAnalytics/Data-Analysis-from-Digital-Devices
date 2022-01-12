data_all = readRDS("data/data_out/data_all_imp_dig.rds")

featureNames = unlist(lapply(data_all, colnames))

featureNames = setdiff(featureNames, "SUBJID")
featureNames = as.data.frame(featureNames)
write.csv(featureNames, "featureNamesADNIforClassifier.csv")
