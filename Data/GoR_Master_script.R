####Gulf of Riga Robustness master script####
#load packages
library(igraph)
library(MASS)
library(NetIndices)
library(fluxweb)
library(NetworkExtinction)
library(network)
library(dplyr)
library(patchwork)
library(tidyverse)

#load data
suffix<-"R1000"
insuffix<-"R1000"
load("GoRoff_v2.Rdata")
load(paste0("ConstSW5_GoRoff_", suffix, ".Rdata"))
load(paste0("GoRoff_FWTSConstS_",insuffix, ".Rdata")) # food web metrics for the 1000 food webs


##load Food web meta web data
load("metawebs_GoR.RData")
load("uw_fw_GoR.RData")
load("w_fw_GoR.RData")


####Networks####

for (i in 1:1000){
  
  #i<-1
  sp.pa.year<-1*((t(biomC[,,i])[, match(row.names(un_mw), row.names(biomC[,,1]))])>0) #reorder to use for making networks (1000 per year)
  
  webs<-list()
  for (j in 1:34){
    fw.sel<-which(sp.pa.year[j,]==1)  
    fw<-un_mw[fw.sel, fw.sel]
    webs[[j]]<-fw
  }
  index.arr[[i]]<-webs
}
index.arr[[10]][[34]]
index.arr[[1000]][[34]]# unweighted networks in matrix
saveRDS(index.arr, file = "index.arr.RDS")
readRDS("index.arr.RDS")



####fluxweb####

w.webs.1000<-list() #list with 1000 foodwebs of which each contains 34, one per year

for (i in 1:1000){
  
  #i<-1
  sp.pa.year<-1*((t(biomC[,,i])[, match(row.names(un_mw), row.names(biomC[,,i]))])>0) #reorder to use for making networks (1000 per year)
  
  w.webs<-list()  
  for (j in 1:34){
    #j<-1
    fw.sel<-which(sp.pa.year[j,]==1)  
    fw<-un_mw[fw.sel, fw.sel]
    
    biom.id<-which(1*((biomC[,,i][,j])>0)==1)
    #length(biom.id)  
    id.order<-match(row.names(fw), names(biom.id))
    ord.mat<-fw[order(id.order),order(id.order)]# web now in same order as biomC
    
    # print(match(colnames(ord.mat), names(biom.id))) #check that they match!
    
    #
    info.id<-which(rownames(info)%in%names(biom.id))
    
    #  print(match(rownames(info[info.id,]), names(biom.id)))#check that info.id and biom.id match
    
    fluxes <- fluxing(ord.mat, biomasses = biomC[,,i][biom.id,j], losses = info[info.id,]$losses,
                      efficiencies = info[info.id,]$efficiencies, ef.level="prey")
    #conversion from J/m2/sec to kJ/m2/day 
    #1 J/m2/sec = 86.4 kJ/m2/day (t here are 86400 sec/day)
    fluxes <- fluxes*86.4
    w.webs[[j]]<-fluxes # store weighted webs when fluxing function runs?
  }
  w.webs.1000[[i]]<-w.webs
}
saveRDS(w.webs.1000, file="w.webs.1000.RDS")
w.webs.1000<-readRDS("w.webs.1000.RDS")