data_dig_meta_vis1 <- read_csv("data/HI-VAE/digitalDataPredicted/data_dig_meta_vis1.csv")

for(names in colnames(data_dig_meta_vis1[,3:ncol(data_dig_meta_vis1)])){
  print(names)
  newColname = gsub("zcode", "scode", names)
  data_dig_meta_vis1[newColname] = rep(0, nrow(data_dig_meta_vis1))
  data = data_dig_meta_vis1[,c(newColname,"SUBJID",names)]
  csvName = paste(gsub("zcode_", "",names), "_meta.csv", sep = "")
  write.csv(data, paste("data/HI-VAE/Saved_Networks/", csvName, sep = ""))
}

write.csv(data_dig_meta_vis1, "data/HI-VAE/digitalDataPredicted/data_dig_meta_vis1_scode.csv")
