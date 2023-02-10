setwd("D:/R_Git/GoR-Food-webs/Rscripts")
w_fw<-readRDS("../Data/weighted_foodwebs.rds")
uw_fw<-readRDS("../Data/unweighted_foodwebs.rds")
uw.top.roles<-readRDS("../Data/uw.top.roles.RDS")
w.top.roles<-readRDS("../Data/w.top.roles.RDS")

library(NetworkExtinction)
library(network)
library(igraph)

#test/ proof of concept code
fw.1981<-uw_fw[[1]]

fw.1981
colnames(fw.1981)<-c(1:25)
rownames(fw.1981)<-c(1:25)
fw.1981
fw<-graph_from_adjacency_matrix(
  fw.1981,
  mode = c("directed"),
  weighted = TRUE,
)
deg.sp<-degree(fw, mode = c("all"),loops = TRUE,normalized = FALSE)
order_degree<-as.numeric(names(sort(deg.sp, decreasing = T)))

#remove basal nodes from extinction order: autotrophs=2, mixotrophs=3, detritus=15
remove<-c(2,3,15)
order_degree<-setdiff(order_degree, remove)

#degree högsta till lägsta 
degree.extinc<-ExtinctionOrder(fw.1981, Order = order_degree, NetworkType = "Trophic")
ExtinctionPlot(History = degree.extinc$sims)

#########################################
#Loop code for extinction by degree
extinc.plots<-list()
par(mfrow=c(4,2))
for(i in 1:length(uw_fw)){
  fw_i<-uw_fw[[i]]
  a<-dim(fw_i)
  colnames(fw_i)<-c(1:a[1])# remove node names so extinction function can run
  rownames(fw_i)<-c(1:a[1])
   
  #degree of nodes
    fw<-graph_from_adjacency_matrix(
      fw_i,
      mode = c("directed"),
      weighted = T
    )
    deg.sp<-degree(fw,mode = c("all"),loops = TRUE,normalized = FALSE)
    
  order_degree<-as.numeric(names(sort(deg.sp, decreasing = T)))

  #remove basal nodes from extinction order: autotrophs=2, mixotrophs=3, detritus=15
  remove<-c(2,3,15)
  order_degree<-setdiff(order_degree, remove)

  #degree högsta till lägsta 
  degree.extinc<-ExtinctionOrder(fw_i, Order = order_degree, NetworkType = "Trophic")
  extinc.plots[i]<-degree.extinc
  #ExtinctionPlot(History = extinc.plots) # can work in loop with list 
  year<-c("1981 degree", "1986 degree", "1991 degree", "1996 degree", "2001 degree", "2006 degree", "2011 degree")
  plot(extinc.plots[[i]]$NumExt, extinc.plots[[i]]$AccSecExt, type = "l",main = year[i], xlab="Number of extinctions", ylab="Accumulated Secondary Extinctions")
}


##### loop for extinctions with among module degree#
extinc.plots<-list()
par(mfrow=c(4,2))
for(i in 1:length(uw_fw)){
  fw_i<-uw_fw[[i]]
  a<-dim(fw_i)
  colnames(fw_i)<-c(1:a[1]) # remove node names so extinction function can run
  rownames(fw_i)<-c(1:a[1])
  
  #among module degree in order
  roles<-uw.top.roles[[i]]
  roles[,3]<-cbind(1:a[1])
  sort<-sort(roles[2:3], decreasing = T)
  order_among.mod<-as.numeric(unlist(sort[,2]))
  #remove basal nodes from extinction order: autotrophs=2, mixotrophs=3, detritus=15
  remove<-c(2,3,15)
  order_among.mod<-setdiff(order_among.mod, remove)
  
  #extinction högsta till lägsta among module degree
  among.mod.extinc<-ExtinctionOrder(fw_i, Order = order_among.mod, NetworkType = "Trophic")
  extinc.plots[i]<-among.mod.extinc
  #ExtinctionPlot(History = extinc.plots) # did not work in loop with list, something about the format
  
  #workaround ploting
  year<-c("1981 among module deg", "1986 among module deg", "1991 among module deg", "1996 among module deg",
          "2001 among module deg", "2006 among module deg", "2011 among module deg")
  plot(extinc.plots[[i]]$NumExt, extinc.plots[[i]]$AccSecExt, type = "l",main = year[i], xlab="Number of extinctions", ylab="Accumulated Secondary Extinctions")
}

##### loop for extinctions with within module degree#
extinc.plots<-list()
par(mfrow=c(4,2))
for(i in 1:length(uw_fw)){
  fw_i<-uw_fw[[i]]
  a<-dim(fw_i)
  colnames(fw_i)<-c(1:a[1]) # remove node names so extinction function can run
  rownames(fw_i)<-c(1:a[1])
  
  #within module degree in order
  roles<-uw.top.roles[[i]]
  roles[,3]<-cbind(1:a[1])
  sort<-sort(roles[1:3], decreasing = T)
  order_within.mod<-as.numeric(unlist(sort[,3]))
  #remove basal nodes from extinction order: autotrophs=2, mixotrophs=3, detritus=15
  remove<-c(2,3,15)
  order_within.mod<-setdiff(order_within.mod, remove)
  
  #extinction högsta till lägsta within module degree
  within.mod.extinc<-ExtinctionOrder(fw_i, Order = order_within.mod, NetworkType = "Trophic")
  extinc.plots[i]<-within.mod.extinc
  #ExtinctionPlot(History = extinc.plots) # did not work in loop with list, something about the format
  
  #workaround ploting
  year<-c("1981 within module deg", "1986 within module deg", "1991 within module deg", "1996 within module deg",
          "2001 within module deg", "2006 within module deg", "2011 within module deg")
  plot(extinc.plots[[i]]$NumExt, extinc.plots[[i]]$AccSecExt, type = "l",main = year[i], xlab="Number of extinctions", ylab="Accumulated Secondary Extinctions")
}



