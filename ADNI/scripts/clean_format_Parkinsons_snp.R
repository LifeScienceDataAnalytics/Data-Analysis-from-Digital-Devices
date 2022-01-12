rm(list=ls())
library(openxlsx)
library(rlist)
library(plyr)
library(stringr)
library(Hmisc)
##library to obtain longest common substring from multiple strings
library(PTXQC)

print(getwd())

# fileList = list.files(path = paste(getwd(), "/Mechanism_Genotypes", sep = ""), pattern = ".*csv")
# 
# dataFiles <- lapply(Sys.glob("Mechanism_Genotypes/*.csv"), read.csv)
# ## import_multiple_csv_files_to_R
# Purpose: Import multiple csv files to the Global Environment in R
setwd("~/Documents/Documents – IT-Admin’s MacBook Pro/ADNIVAMBN/VAMBNForADNI")
# set working directory
setwd(paste(getwd(), "/Mechanism_Genotypes", sep = ""))

# list all csv files from the current directory
list.files(pattern=".csv$") # use the pattern argument to define a common pattern  for import files with regex. Here: .csv

# create a list from these files
list.filenames<-list.files(pattern=".csv$")
list.filenames

# create an empty list that will serve as a container to receive the incoming files
list.data<-list()

# create a loop to read in your data
for (i in 1:length(list.filenames))
{
  list.data[[i]]<-read.csv(list.filenames[i])
}

# add the names of your data to the list
names(list.data)<-list.filenames

toBeRemoved = c()
dataList = list()
for(i in 1:length(list.data)){
  fileDf = as.data.frame(list.data[[i]])
  print(names(fileDf)[1])
  names(fileDf)[1] = "SUBJID"
  if(ncol(fileDf) == 2){
    print(ncol(fileDf))
    toBeRemoved = c(toBeRemoved,names(list.data[i]))
  }
  dataList[[length(dataList)+1]] = fileDf
}
print(toBeRemoved)

names(dataList)<-list.filenames
for(sugrph in toBeRemoved){
  dataList[sugrph] = NULL
}
data_all = dataList
setwd("~/Documents/Documents – IT-Admin’s MacBook Pro/ADNIVAMBN/VAMBNForADNI")
saveRDS(data_all, file = paste0("data/data_condensed_ParkSNP.rds"))


