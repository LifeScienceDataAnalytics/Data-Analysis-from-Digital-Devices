rm(list=ls()) 
library(igraph)
library(bnlearn)
library(stringr)

### load final BN from ADNI
finalBNADNI = readRDS("/Users/meemansasood/Documents/Documents_IT/ADNIVAMBN_paper_all/VAMBNForADNI/data/data_out/main_finalBN.rds")

## load final BN from Altoida
finalBNAlt = readRDS("/Users/meemansasood/Documents/Documents_IT/VAMBN_Version5_mmseEncoded_final/data/data_out/main_finalBN.rds")

#converted the bn's to edgelists
#edgesBNadni = finalBNADNI$arcs



edgesBNadni = finalBNADNI$arcs
edgesBNadni = edgesBNadni[- grep("scode", edgesBNadni[,1]),]
edgesBNadni = edgesBNadni[- grep("scode", edgesBNadni[,2]),]
edgesBNadni = edgesBNadni[- grep("AUX", edgesBNadni[,1]),]
edgesBNadni = edgesBNadni[- grep("visitmiss", edgesBNadni[,1]),]
edgesBNadni[,1] = gsub("^SA_zcode", "zcode", edgesBNadni[,1])
edgesBNadni[,2] = gsub("^SA_zcode", "zcode", edgesBNadni[,2])

edgesBNadniCp = edgesBNadni

edgesBNalt = finalBNAlt$arcs
edgesBNalt = edgesBNalt[- grep("scode",  edgesBNalt[,1]),]
edgesBNalt = edgesBNalt[- grep("scode",  edgesBNalt[,2]),]
edgesBNalt = edgesBNalt[- grep("AUX", edgesBNalt[,1]),]
edgesBNalt = edgesBNalt[- grep("visitmiss", edgesBNalt[,1]),]
#then to igraph objects
adni_igraph = graph_from_edgelist(edgesBNadni)
alt_igraph = graph_from_edgelist(edgesBNalt)

#checked the overlap and unique edges of each net
##change the adni_graph names
edgesBNadniCp[,1] = gsub("[[:digit:]]+", "1", edgesBNadniCp[,1])
edgesBNadniCp[,2] = gsub("[[:digit:]]+", "1", edgesBNadniCp[,2])
adni_igraphCp = graph_from_edgelist(edgesBNadniCp)
int_edges = as_edgelist(intersection(alt_igraph, adni_igraphCp)) # 29 nodes, 32 edges
write.csv(int_edges, "int_edges.csv")
onlyAdniBN = as_edgelist(difference(adni_igraph, alt_igraph)) # in bn but not in kg -> 124 edges
onlyAltBN = as_edgelist(difference(alt_igraph, adni_igraph)) # in kg but not in bn -> 0 edges



distMatrixAlt = distances(alt_igraph, v = V(alt_igraph), to = V(alt_igraph), mode = c("all", "out",
                                                                                 "in"), weights = NULL, algorithm = c("automatic", "unweighted",
                                                                                                                           "dijkstra", "bellman-ford", "johnson"))

##subset the distance matrix for Altoida
digFeatures = grep("AR|BIT|DOT|zcode\\_Motor", colnames(distMatrixAlt), value = TRUE)
mmseAndDemogs = grep("MMSE|DX|Education|age|Amyloid", colnames(distMatrixAlt), value = TRUE)
distMatrixAlt = distMatrixAlt[rownames(distMatrixAlt)%in%digFeatures,colnames(distMatrixAlt)%in%mmseAndDemogs]




write.csv(distMatrixAlt, "distMatrixAlt.csv")

##shortes path in adni network
distance_table(adni_igraph, directed = TRUE)


distMatrixADNI = distances(adni_igraph, v = V(adni_igraph), 
                           to = V(adni_igraph), mode = c("all", "out",
                                                       "in"), weights = NULL, algorithm = c("automatic", "unweighted","dijkstra", "bellman-ford", "johnson"))

digFeaturesADNI = grep("zcode\\_AR|BIT|DOT|zcode\\_Motor", colnames(distMatrixADNI), value = TRUE)
mmseFAQPathsVolDemogs = grep("MMSE|DX|PTEDUCAT|AGE|Amyloid", colnames(distMatrixADNI), value = TRUE)

distMatrixADNI = distMatrixADNI[rownames(distMatrixADNI)%in%digFeaturesADNI,colnames(distMatrixADNI)%in%mmseFAQPathsVolDemogs]
#rownames(distMatrixADNI) <- gsub("[[:digit:]]+", "1", rownames(distMatrixADNI))
names(distMatrixADNI)[names(distMatrixADNI)=="SA_AGE_VIS1"] = "SA_age_VIS1"

#threshold = 2
#data = distMatrixADNI
func.threshold = function(data, threshold){
  thresholdVec = c()
  for(i in 1:nrow(data)){
    for(j in 1:ncol(data)){
      if(data[i,j] <= threshold){
        from = (rownames(data)[i])
        to = (colnames(data)[j])
        fromUpd = gsub("[[:digit:]]+", "1", from)
        #print(from)
        toUpd = gsub("[[:digit:]]+", "1", to)
        #print(to)
        thresholdVec = c(thresholdVec, paste(fromUpd, fromUpd, data[from, to], sep = ","))
      }
    }
  }
  return(thresholdVec)
}

uniqueADNI = sort(unique(as.numeric(as.character(unique(unlist(distMatrixADNI))))))
uniqueAlt = sort(unique(as.numeric(as.character(unique(unlist(distMatrixAlt))))))
uniqueThreshold = intersect(uniqueADNI, uniqueAlt)
uniqueThreshold = c(1,2,3,4)
thresholdAlt1 = func.threshold(distMatrixAlt, uniqueADNI[1])
thresholdAlt2 = func.threshold(distMatrixAlt, uniqueADNI[2])
thresholdAlt3 = func.threshold(distMatrixAlt, uniqueADNI[3])
thresholdAlt4 = func.threshold(distMatrixAlt, uniqueADNI[4])
thresholdAlt5 = func.threshold(distMatrixAlt, uniqueADNI[5])
thresholdAlt6 = func.threshold(distMatrixAlt, uniqueADNI[6])
thresholdAlt7 = func.threshold(distMatrixAlt, uniqueADNI[7])
allAltCombs = c(thresholdAlt1,thresholdAlt2,thresholdAlt3,thresholdAlt4,thresholdAlt5,thresholdAlt6,thresholdAlt7)
allAltCombs = unique(allAltCombs)
fractionDf = data.frame()
for(i in 1:length(uniqueThreshold)){
  #print(i)
  thresholdADNI = func.threshold(distMatrixADNI, uniqueThreshold[i])
  thresholdADNI = unique(thresholdADNI)
  #print(thresholdADNI)
  thresholdAlt = func.threshold(distMatrixAlt, uniqueThreshold[i])
  thresholdAlt = unique(thresholdAlt)
  #print(thresholdAlt)
 #print(length(intersect(thresholdAlt, thresholdADNI)))
  #print(length(thresholdAlt))
  fraction = (length(intersect(thresholdAlt, thresholdADNI)))/length(allAltCombs)
  print(fraction)
  fractionDf[i, "threshold"] = uniqueADNI[i]
  fractionDf[i, "fraction"] = fraction
}

fractionDf = na.omit(fractionDf)
plot(fractionDf$threshold, fractionDf$fraction, type = "b", pch = 19, 
     col = "black", xlab = "threshold", ylab = "fraction of validated edges", ylim = c(0,1))

DistGraph(adni_igraph, dist.method="shortest.paths")
write.csv(distMatrixADNI, "distanceMatrixADNI.csv")
# Classical MDS
# N rows (objects) x p columns (variables)
# each row identified by a unique row name
path.length.hist (adni_igraph, directed = TRUE) 
heatmap(distMatrixADNI, cexRow=0.3, cexCol = 0.22)
heatmap(distMatrixAlt, cexRow=0.3, cexCol = 0.22)
#distance_table(adni_igraph)
s = all_shortest_paths(adni_igraph, from="zcode_ARPlaceAndFindTelemetryVariance_VIS12", to = V(adni_igraph), mode = c("out", "all",
                                                        "in"), weights = NULL)


dist.from.NYT <- distances(adni_igraph, v=V(adni_igraph)[from=="zcode_ARPlaceAndFindTelemetryVariance_VIS12"], to=V(adni_igraph), weights=NA)





