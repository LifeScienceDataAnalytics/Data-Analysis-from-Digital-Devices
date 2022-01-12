  ############# README
# This is to plot the network for the paper
############# 
rm(list=ls())
library(tidyverse)
library(igraph)
library(bnlearn) # hc might be overwritten by arules or some such package "bnlearn::hc" if so; not currently used though
library(parallel)
library(stringr)
library(igraph)
library(RCy3)

########## Get models
name<-'Final Bootstrap'
data_out<-'data/data_out/main'
finalBN<-readRDS(paste0(data_out,'_finalBN.rds'))
final_graph<-igraph.from.graphNEL(as.graphNEL(finalBN))
edg.final<-as.data.frame(as_edgelist(final_graph, names = TRUE),stringsAsFactors=F)
colnames(edg.final)<-c('from','to')
final_graph<-set_vertex_attr(final_graph,'visit',value = as.numeric(gsub('[a-zA-Z0-9_]{1,}_VIS','',V(final_graph)$name)))
V(final_graph)$name<-gsub('SA_|VIS','',V(final_graph)$name)

bootBN<-readRDS(paste0(data_out,'_bootBN.rds'))
#boot_subgraph = bootBN[bootBN$strength >0.4 & bootBN$direction >0.4 , ]
boot_subgraph =merge(x=bootBN,y=edg.final,by=c('from','to'))
boot_subgraph$strength = as.character(round(boot_subgraph$strength, 2 ))
boot_graph<-graph_from_data_frame(boot_subgraph)
boot_graph<-set_vertex_attr(boot_graph,'visit',value = as.numeric(gsub('[a-zA-Z0-9_]{1,}_VIS','',V(boot_graph)$name)))
V(boot_graph)$name<-gsub('SA_|VIS','',V(boot_graph)$name)

#overlap<-intersection(boot_graph,final_graph)
#E(boot_graph)$in_final<-ifelse(E(boot_graph) %in% E(final_graph),1,0)

# set attributes
E(boot_graph)$strength<-as.character(E(boot_graph)$strength)
E(boot_graph)$direction<-as.character(E(boot_graph)$direction)

boot_graph<-delete_vertices(boot_graph, V(boot_graph)$name[grepl('AUX_|visitmiss_',V(boot_graph)$name)])

netcont<-createNetworkFromIgraph(boot_graph,name,collection="Paper Plot")
layoutNetwork('hierarchical')
setVisualStyle('default')