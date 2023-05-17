
setwd("C:/LocalData/susakort/GitHub/GoR-Food-webs/Data")


#load packages
library(igraph)
library(MASS)
library(NetIndices)
library(fluxweb)
library(NetworkExtinction)
library(network)

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


#reduced iterations to 10 for testing purposes
index<-matrix(NA, 34, 10)
deg<-matrix(NA, 34, 10)
R50<-matrix(NA, 34, 10)
n.sp.extR50<-matrix(NA, 34, 10)
extinctionoutput<-list()
w.webs<-list()
#index.arr <- list() #could be list of output for the 1000 food webs per year

for (i in 1:10){

  #i<-1
  sp.pa.year<-1*((t(biomC[,,i])[, match(row.names(un_mw), row.names(biomC[,,1]))])>0) #reorder to use for making networks (1000 per year)
  

  for (j in 1:34){
  fw.sel<-which(sp.pa.year[j,]==1)  
  fw<-un_mw[fw.sel, fw.sel]
  
  #WHAT METRIC DO YOU WANT TO CALCULATE ADD HERE
  #Example with number of species per network (unweighted)
  
  index[j,i]<-dim(fw)[1] #number of species
  deg[j,i]<-c(mean(degree(graph.adjacency(fw)))) # mean degree test, will add more later when weights are up and running
  #links
  #generality
  #vulnrability
  #connectance
  #Weighted metrics
  
  #R50 test run with un-weighted extinction by degree, when fluxes are integrated can change it to/ add extinction by sum of link weights
  fw_i<-fw
  colnames(fw_i)<-c(1:dim(fw)[1])# remove node names so extinction function can run
  rownames(fw_i)<-c(1:dim(fw)[1])# remove node names so extinction function can run
  
  #degree of nodes
  fw1<-graph_from_adjacency_matrix(
    fw_i,
    mode = c("directed"),
    weighted = T
  )
  deg.sp<-degree(fw1, mode = c("all"),loops = TRUE,normalized = FALSE)
  order_degree<-as.numeric(names(sort(deg.sp, decreasing = T))) #order of nodes by degree
  
  #remove basal nodes from extinction order: autotrophs, mixotrophs, detritus
  order_degree<-setdiff(order_degree, c(which(colnames(fw)=="Autotroph"),
                                        which(colnames(fw)=="Mixotroph"),
                                        which(colnames(fw)=="Detritus")))
  
  #ext. by highest to lowest degree
  degree.extinc.mw<-ExtinctionOrder(fw_i, Order = order_degree, NetworkType = "Trophic",verbose = F)
  extinctionoutput[[i]]<-degree.extinc.mw
  #R50
  R50[j,i]<-degree.extinc.mw$R50 # the R50 index, between 1 and 0
  n.sp.extR50[j,i]<-length(which(degree.extinc.mw$sims$TotalExt<= dim(fw_i)[1]/2)) # number of species extinct at R50
  

  #YOU NEED TO DO THE SAME FOR WEIGHTED NETWORKS
  #YOU CAN USE THE BIOMASS INFO FOR THIS AND THE UNEWEIGHTED fw
  
 #se work in progress below
  
  #INSTEAD OF MATRICES YOU CAN ASLO STORE YOUR OUTPUT IN LISTS
  #THESE LISTS WE CAN ADD TO ARRAYS AND THEN CALCULATE THE MEAN VALUES
  
  }
}

index.median<- apply(index, 1, median)
round(apply(deg, 1, median), 2)
apply(R50, 1, median)
apply(n.sp.extR50,1,median)

extinctionoutput

#### Work in progress:

#fluxweb

w.webs.1000<-list() #list with 1000 foodwebs of which each contains 34, one per year

for (i in 1:10){
 
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

    print(match(colnames(ord.mat), names(biom.id))) #check that they match!
    
    #
    info.id<-which(rownames(info)%in%names(biom.id))
    
    print(match(rownames(info[info.id,]), names(biom.id)))#check that info.id and biom.id match
    
    fluxes <- fluxing(ord.mat, biomasses = biomC[,,i][biom.id,j], losses = info[info.id,]$losses,
                  efficiencies = info[info.id,]$efficiencies, ef.level="prey")
    #conversion from J/m2/sec to kJ/m2/day 
    #1 J/m2/sec = 86.4 kJ/m2/day (t here are 86400 sec/day)
    fluxes <- fluxes*86.4
    w.webs[[j]]<-fluxes # store weighted webs when fluxing function runs?
  }
  w.webs.1000[[i]]<-w.webs
  }



#examples from Susanne

#Example ordering & matching 
id.ord<-match(colnames(bm),rownames(traits_all))
bm.mat<-bm[,order(id.ord)]
      
#To make the array example
years<-34
dim.index<-5
no.fw<-1000

arr <- array(unlist(index), dim = c(years, dim.index, no.fw))
dimnames(arr) <- list(dimnames(biomC)[[2]],
                      colnames(dim.index), 
                      1:dim(biomC)[3])

#calc median value of array
apply(arr, 1:2, median, na.rm = TRUE)




