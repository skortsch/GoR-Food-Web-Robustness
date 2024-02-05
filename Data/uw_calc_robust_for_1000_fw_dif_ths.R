webs<-readRDS("index.arr.RDS")

uw.R50th10<-matrix(NA, 34, 1000)
#n.sp.extR50<-matrix(NA, 34, 10)
#uw.ext.ind.th10<-list()

for (i in 1:1000) {
  #i<-1
  a<-webs[[i]]
  
  #
  for (j in 1:34) {
    #j<-1
    fw<-a[[j]]
    
    fw_i<-fw
    colnames(fw_i)<-c(1:dim(fw)[1])# remove node names so extinction function can run
    rownames(fw_i)<-c(1:dim(fw)[1])# remove node names so extinction function can run
    
    #degree of nodes
    fw1<-graph_from_adjacency_matrix(
      fw_i,
      mode = c("directed"),
      weighted = T
    )
    sum.link.weights.sp<-apply(fw_i, 1, sum)+apply(fw_i, 2, sum)
    order_sum_link<-as.numeric(names(sort(sum.link.weights.sp, decreasing = T)))
    
    
    #remove basal nodes from extinction order: autotrophs, mixotrophs, detritus
    order_sum_link<-setdiff(order_sum_link, c(which(colnames(fw)=="Autotroph"),
                                              which(colnames(fw)=="Mixotroph"),
                                              which(colnames(fw)=="Detritus")))
    
    #ext. by highest to lowest sum of link weights
    sumlink.extinc.mw<-ExtinctionOrder(fw_i, Order = order_sum_link, NetworkType = "Trophic",verbose = F, IS=0.1)
    #ext.output[[j]]<-sumlink.extinc.mw
    #R50
   uw.R50th10[j,i]<-sumlink.extinc.mw$R50 ## R50 value extracted from the ext.output, R50 between 1 and 0
    #n.sp.extR50[i,j]<-length(which(degree.extinc.mw$sims$TotalExt<= dim(fw_i)[1]/2)) # number of species extinct at R50
    
  }
 #uw.ext.ind.th10[[i]]<-ext.output
}
#uw.ext.ind.th10[[1000]][[34]]
uw.R50th10
#saveRDS(uw.ext.ind.th10,file="uw.ext.ind.th10.RDS")
saveRDS(uw.R50th10,file="uw.R50th10.RDS")
apply(uw.R50th10,1,median)
readRDS("uw.R50th10.RDS")



uw.R50th20<-matrix(NA, 34, 1000)
#n.sp.extR50<-matrix(NA, 34, 10)
#uw.ext.ind.th20<-list()

for (i in 1:1000) {
  #i<-1
  a<-webs[[i]]
  
  ##
  for (j in 1:34) {
    #j<-1
    fw<-a[[j]]
    
    fw_i<-fw
    colnames(fw_i)<-c(1:dim(fw)[1])# remove node names so extinction function can run
    rownames(fw_i)<-c(1:dim(fw)[1])# remove node names so extinction function can run
    
    #degree of nodes
    fw1<-graph_from_adjacency_matrix(
      fw_i,
      mode = c("directed"),
      weighted = T
    )
    sum.link.weights.sp<-apply(fw_i, 1, sum)+apply(fw_i, 2, sum)
    order_sum_link<-as.numeric(names(sort(sum.link.weights.sp, decreasing = T)))
    
    
    #remove basal nodes from extinction order: autotrophs, mixotrophs, detritus
    order_sum_link<-setdiff(order_sum_link, c(which(colnames(fw)=="Autotroph"),
                                              which(colnames(fw)=="Mixotroph"),
                                              which(colnames(fw)=="Detritus")))
    
    #ext. by highest to lowest sum of link weights
    sumlink.extinc.mw<-ExtinctionOrder(fw_i, Order = order_sum_link, NetworkType = "Trophic",verbose = F, IS=0.2)
    #ext.output[[j]]<-sumlink.extinc.mw
    #R50
   uw.R50th20[j,i]<-sumlink.extinc.mw$R50 ## R50 value extracted from the ext.output, R50 between 1 and 0
    #n.sp.extR50[i,j]<-length(which(degree.extinc.mw$sims$TotalExt<= dim(fw_i)[1]/2)) # number of species extinct at R50
    
  }
 #uw.ext.ind.th20[[i]]<-ext.output
}
#saveRDS(uw.ext.ind.th20, file="uw.ext.ind.th20.RDS")
saveRDS(uw.R50th20, file="uw.R50th20.RDS")
apply(uw.R50th20, 1, median)

readRDS("uw.R50th20.RDS")


uw.R50th30<-matrix(NA, 34, 1000)
#n.sp.extR50<-matrix(NA, 34, 10)
#uw.ext.ind.th30<-list()

for (i in 1:1000) {
  #i<-1
  a<-webs[[i]]
  
  ##
  for (j in 1:34) {
    #j<-1
    fw<-a[[j]]
    
    fw_i<-fw
    colnames(fw_i)<-c(1:dim(fw)[1])# remove node names so extinction function can run
    rownames(fw_i)<-c(1:dim(fw)[1])# remove node names so extinction function can run
    
    #degree of nodes
    fw1<-graph_from_adjacency_matrix(
      fw_i,
      mode = c("directed"),
      weighted = T
    )
    sum.link.weights.sp<-apply(fw_i, 1, sum)+apply(fw_i, 2, sum)
    order_sum_link<-as.numeric(names(sort(sum.link.weights.sp, decreasing = T)))
    
    
    #remove basal nodes from extinction order: autotrophs, mixotrophs, detritus
    order_sum_link<-setdiff(order_sum_link, c(which(colnames(fw)=="Autotroph"),
                                              which(colnames(fw)=="Mixotroph"),
                                              which(colnames(fw)=="Detritus")))
    
    #ext. by highest to lowest sum of link weights
    sumlink.extinc.mw<-ExtinctionOrder(fw_i, Order = order_sum_link, NetworkType = "Trophic",verbose = F, IS=0.3)
    #ext.output[[j]]<-sumlink.extinc.mw
    #R50
    uw.R50th30[j,i]<-sumlink.extinc.mw$R50 ## R50 value extracted from the ext.output, R50 between 1 and 0
    #n.sp.extR50[i,j]<-length(which(degree.extinc.mw$sims$TotalExt<= dim(fw_i)[1]/2)) # number of species extinct at R50
    
  }
  #uw.ext.ind.th30[[i]]<-ext.output
}
#saveRDS(uw.ext.ind.th30, file="uw.ext.ind.th30.RDS")
saveRDS(uw.R50th30, file="uw.R50th30.RDS")
apply(uw.R50th30, 1, median)

readRDS("uw.R50th30.RDS")


uw.R50th40<-matrix(NA, 34, 1000)
#n.sp.extR50<-matrix(NA, 34, 10)
#uw.ext.ind.th40<-list()

for (i in 1:1000) {
  #i<-1
  a<-webs[[i]]
  
  #
  for (j in 1:34) {
    #j<-1
    fw<-a[[j]]
    
    fw_i<-fw
    colnames(fw_i)<-c(1:dim(fw)[1])# remove node names so extinction function can run
    rownames(fw_i)<-c(1:dim(fw)[1])# remove node names so extinction function can run
    
    #degree of nodes
    fw1<-graph_from_adjacency_matrix(
      fw_i,
      mode = c("directed"),
      weighted = T
    )
    sum.link.weights.sp<-apply(fw_i, 1, sum)+apply(fw_i, 2, sum)
    order_sum_link<-as.numeric(names(sort(sum.link.weights.sp, decreasing = T)))
    
    
    #remove basal nodes from extinction order: autotrophs, mixotrophs, detritus
    order_sum_link<-setdiff(order_sum_link, c(which(colnames(fw)=="Autotroph"),
                                              which(colnames(fw)=="Mixotroph"),
                                              which(colnames(fw)=="Detritus")))
    
    #ext. by highest to lowest sum of link weights
    sumlink.extinc.mw<-ExtinctionOrder(fw_i, Order = order_sum_link, NetworkType = "Trophic",verbose = F, IS=0.4)
    #ext.output[[j]]<-sumlink.extinc.mw
    #R50
   uw.R50th40[j,i]<-sumlink.extinc.mw$R50 ## R50 value extracted from the ext.output, R50 between 1 and 0
    #n.sp.extR50[i,j]<-length(which(degree.extinc.mw$sims$TotalExt<= dim(fw_i)[1]/2)) # number of species extinct at R50
    
  }
 #uw.ext.ind.th40[[i]]<-ext.output
}
#saveRDS(uw.ext.ind.th40, file="uw.ext.ind.th40.RDS")
saveRDS(uw.R50th40, file="uw.R50th40.RDS")
apply(uw.R50th40, 1, median)

readRDS("uw.R50th40.RDS")





uw.R50th50<-matrix(NA, 34, 1000)
#n.sp.extR50<-matrix(NA, 34, 10)
#uw.ext.ind.th50<-list()

for (i in 1:1000) {
  #i<-1
  a<-webs[[i]]
  
  #
  for (j in 1:34) {
    #j<-1
    fw<-a[[j]]
    
    fw_i<-fw
    colnames(fw_i)<-c(1:dim(fw)[1])# remove node names so extinction function can run
    rownames(fw_i)<-c(1:dim(fw)[1])# remove node names so extinction function can run
    
    #degree of nodes
    fw1<-graph_from_adjacency_matrix(
      fw_i,
      mode = c("directed"),
      weighted = T
    )
    sum.link.weights.sp<-apply(fw_i, 1, sum)+apply(fw_i, 2, sum)
    order_sum_link<-as.numeric(names(sort(sum.link.weights.sp, decreasing = T)))
    
    
    #remove basal nodes from extinction order: autotrophs, mixotrophs, detritus
    order_sum_link<-setdiff(order_sum_link, c(which(colnames(fw)=="Autotroph"),
                                              which(colnames(fw)=="Mixotroph"),
                                              which(colnames(fw)=="Detritus")))
    
    #ext. by highest to lowest sum of link weights
    sumlink.extinc.mw<-ExtinctionOrder(fw_i, Order = order_sum_link, NetworkType = "Trophic",verbose = F, IS=0.5)
   # ext.output[[j]]<-sumlink.extinc.mw
    #R50
   uw.R50th50[j,i]<-sumlink.extinc.mw$R50 ## R50 value extracted from the ext.output, R50 between 1 and 0
    #n.sp.extR50[i,j]<-length(which(degree.extinc.mw$sims$TotalExt<= dim(fw_i)[1]/2)) # number of species extinct at R50
    
  }
 #uw.ext.ind.th50[[i]]<-ext.output
}
#saveRDS(uw.ext.ind.th50, file="uw.ext.ind.th50.RDS")
saveRDS(uw.R50th50, file="uw.R50th50.RDS")
apply(uw.R50th50, 1, median)

readRDS("uw.R50th50.RDS")






uw.R50th60<-matrix(NA, 34, 1000)
#n.sp.extR50<-matrix(NA, 34, 10)
#uw.ext.ind.th60<-list()

for (i in 1:1000) {
  #i<-1
  a<-webs[[i]]
  
  #
  for (j in 1:34) {
    #j<-1
    fw<-a[[j]]
    
    fw_i<-fw
    colnames(fw_i)<-c(1:dim(fw)[1])# remove node names so extinction function can run
    rownames(fw_i)<-c(1:dim(fw)[1])# remove node names so extinction function can run
    
    #degree of nodes
    fw1<-graph_from_adjacency_matrix(
      fw_i,
      mode = c("directed"),
      weighted = T
    )
    sum.link.weights.sp<-apply(fw_i, 1, sum)+apply(fw_i, 2, sum)
    order_sum_link<-as.numeric(names(sort(sum.link.weights.sp, decreasing = T)))
    
    
    #remove basal nodes from extinction order: autotrophs, mixotrophs, detritus
    order_sum_link<-setdiff(order_sum_link, c(which(colnames(fw)=="Autotroph"),
                                              which(colnames(fw)=="Mixotroph"),
                                              which(colnames(fw)=="Detritus")))
    
    #ext. by highest to lowest sum of link weights
    sumlink.extinc.mw<-ExtinctionOrder(fw_i, Order = order_sum_link, NetworkType = "Trophic",verbose = F, IS=0.6)
  #  ext.output[[j]]<-sumlink.extinc.mw
    #R50
   uw.R50th60[j,i]<-sumlink.extinc.mw$R50 ## R50 value extracted from the ext.output, R50 between 1 and 0
    #n.sp.extR50[i,j]<-length(which(degree.extinc.mw$sims$TotalExt<= dim(fw_i)[1]/2)) # number of species extinct at R50
    
  }
 #uw.ext.ind.th60[[i]]<-ext.output
}
#saveRDS(uw.ext.ind.th60, file="uw.ext.ind.th60.RDS")
saveRDS(uw.R50th60, file="uw.R50th60.RDS")
apply(uw.R50th60, 1, median)

readRDS("uw.R50th60.RDS")





uw.R50th70<-matrix(NA, 34, 1000)
#n.sp.extR50<-matrix(NA, 34, 10)
uw.ext.ind.th70<-list()

for (i in 1:1000) {
  #i<-1
  a<-webs[[i]]
  
  #
  for (j in 1:34) {
    #j<-1
    fw<-a[[j]]
    
    fw_i<-fw
    colnames(fw_i)<-c(1:dim(fw)[1])# remove node names so extinction function can run
    rownames(fw_i)<-c(1:dim(fw)[1])# remove node names so extinction function can run
    
    #degree of nodes
    fw1<-graph_from_adjacency_matrix(
      fw_i,
      mode = c("directed"),
      weighted = T
    )
    sum.link.weights.sp<-apply(fw_i, 1, sum)+apply(fw_i, 2, sum)
    order_sum_link<-as.numeric(names(sort(sum.link.weights.sp, decreasing = T)))
    
    
    #remove basal nodes from extinction order: autotrophs, mixotrophs, detritus
    order_sum_link<-setdiff(order_sum_link, c(which(colnames(fw)=="Autotroph"),
                                              which(colnames(fw)=="Mixotroph"),
                                              which(colnames(fw)=="Detritus")))
    
    #ext. by highest to lowest sum of link weights
    sumlink.extinc.mw<-ExtinctionOrder(fw_i, Order = order_sum_link, NetworkType = "Trophic",verbose = F, IS=0.7)
    ext.output[[j]]<-sumlink.extinc.mw
    #R50
   uw.R50th70[j,i]<-sumlink.extinc.mw$R50 ## R50 value extracted from the ext.output, R50 between 1 and 0
    #n.sp.extR50[i,j]<-length(which(degree.extinc.mw$sims$TotalExt<= dim(fw_i)[1]/2)) # number of species extinct at R50
    
  }
 uw.ext.ind.th70[[i]]<-ext.output
}
#saveRDS(uw.ext.ind.th70, file="uw.ext.ind.th70.RDS")
saveRDS(uw.R50th70, file="uw.R50th70.RDS")
apply(uw.R50th70, 1, median)

readRDS("uw.R50th70.RDS")




uw.R50th80<-matrix(NA, 34, 1000)
#n.sp.extR50<-matrix(NA, 34, 10)
uw.ext.ind.th80<-list()

for (i in 1:1000) {
  #i<-1
  a<-webs[[i]]
  
  #
  for (j in 1:34) {
    #j<-1
    fw<-a[[j]]
    
    fw_i<-fw
    colnames(fw_i)<-c(1:dim(fw)[1])# remove node names so extinction function can run
    rownames(fw_i)<-c(1:dim(fw)[1])# remove node names so extinction function can run
    
    #degree of nodes
    fw1<-graph_from_adjacency_matrix(
      fw_i,
      mode = c("directed"),
      weighted = T
    )
    sum.link.weights.sp<-apply(fw_i, 1, sum)+apply(fw_i, 2, sum)
    order_sum_link<-as.numeric(names(sort(sum.link.weights.sp, decreasing = T)))
    
    
    #remove basal nodes from extinction order: autotrophs, mixotrophs, detritus
    order_sum_link<-setdiff(order_sum_link, c(which(colnames(fw)=="Autotroph"),
                                              which(colnames(fw)=="Mixotroph"),
                                              which(colnames(fw)=="Detritus")))
    
    #ext. by highest to lowest sum of link weights
    sumlink.extinc.mw<-ExtinctionOrder(fw_i, Order = order_sum_link, NetworkType = "Trophic",verbose = F, IS=0.8)
    ext.output[[j]]<-sumlink.extinc.mw
    #R50
   uw.R50th80[j,i]<-sumlink.extinc.mw$R50 ## R50 value extracted from the ext.output, R50 between 1 and 0
    #n.sp.extR50[i,j]<-length(which(degree.extinc.mw$sims$TotalExt<= dim(fw_i)[1]/2)) # number of species extinct at R50
    
  }
 uw.ext.ind.th80[[i]]<-ext.output
}
#saveRDS(uw.ext.ind.th80, file="uw.ext.ind.th80.RDS")
saveRDS(uw.R50th80, file="uw.R50th80.RDS")
apply(uw.R50th80, 1, median)

readRDS("uw.R50th80.RDS")



uw.R50th90<-matrix(NA, 34, 1000)
#n.sp.extR50<-matrix(NA, 34, 10)
uw.ext.ind.th90<-list()

for (i in 1:1000) {
  #i<-1
  a<-webs[[i]]
  
  #
  for (j in 1:34) {
    #j<-1
    fw<-a[[j]]
    
    fw_i<-fw
    colnames(fw_i)<-c(1:dim(fw)[1])# remove node names so extinction function can run
    rownames(fw_i)<-c(1:dim(fw)[1])# remove node names so extinction function can run
    
    #degree of nodes
    fw1<-graph_from_adjacency_matrix(
      fw_i,
      mode = c("directed"),
      weighted = T
    )
    sum.link.weights.sp<-apply(fw_i, 1, sum)+apply(fw_i, 2, sum)
    order_sum_link<-as.numeric(names(sort(sum.link.weights.sp, decreasing = T)))
    
    
    #remove basal nodes from extinction order: autotrophs, mixotrophs, detritus
    order_sum_link<-setdiff(order_sum_link, c(which(colnames(fw)=="Autotroph"),
                                              which(colnames(fw)=="Mixotroph"),
                                              which(colnames(fw)=="Detritus")))
    
    #ext. by highest to lowest sum of link weights
    sumlink.extinc.mw<-ExtinctionOrder(fw_i, Order = order_sum_link, NetworkType = "Trophic",verbose = F, IS=0.9)
    ext.output[[j]]<-sumlink.extinc.mw
    #R50
   uw.R50th90[j,i]<-sumlink.extinc.mw$R50 ## R50 value extracted from the ext.output, R50 between 1 and 0
    #n.sp.extR50[i,j]<-length(which(degree.extinc.mw$sims$TotalExt<= dim(fw_i)[1]/2)) # number of species extinct at R50
    
  }
 uw.ext.ind.th90[[i]]<-ext.output
}
#saveRDS(uw.ext.ind.th90, file="uw.ext.ind.th90.RDS")
saveRDS(uw.R50th90, file="uw.R50th90.RDS")
apply(uw.R50th90, 1, median)

readRDS("uw.R50th90.RDS")




uw.R50Binary<-matrix(NA, 34, 1000)
#n.sp.extR50<-matrix(NA, 34, 10)
uw.ext.ind.Binary<-list()

for (i in 1:1000) {
  #i<-1
  a<-webs[[i]]
  
  #
  for (j in 1:34) {
    #j<-1
    fw<-a[[j]]
    
    fw_i<-fw
    colnames(fw_i)<-c(1:dim(fw)[1])# remove node names so extinction function can run
    rownames(fw_i)<-c(1:dim(fw)[1])# remove node names so extinction function can run
    
    #order of extinction, sum of linkweights
    fw1<-graph_from_adjacency_matrix(
      fw_i,
      mode = c("directed"),
      weighted = T
    )
    sum.link.weights.sp<-apply(fw_i, 1, sum)+apply(fw_i, 2, sum)
    order_sum_link<-as.numeric(names(sort(sum.link.weights.sp, decreasing = T)))
    
    
    #remove basal nodes from extinction order: autotrophs, mixotrophs, detritus
    order_sum_link<-setdiff(order_sum_link, c(which(colnames(fw)=="Autotroph"),
                                              which(colnames(fw)=="Mixotroph"),
                                              which(colnames(fw)=="Detritus")))
    
    #ext. by highest to lowest sum of link weights
    sumlink.extinc.mw<-ExtinctionOrder(fw_i, Order = order_sum_link, NetworkType = "Trophic",verbose = F)
    ext.output[[j]]<-sumlink.extinc.mw
    #R50
   uw.R50Binary[j,i]<-sumlink.extinc.mw$R50 ## R50 value extracted from the ext.output, R50 between 1 and 0
    #n.sp.extR50[i,j]<-length(which(degree.extinc.mw$sims$TotalExt<= dim(fw_i)[1]/2)) # number of species extinct at R50
    
  }
 uw.ext.ind.Binary[[i]]<-ext.output
}
saveRDS(uw.ext.ind.Binary, file="uw.ext.ind.Binary.RDS")
saveRDS(uw.R50Binary, file="uw.R50Binary.RDS")
apply(uw.R50Binary, 1, median)

readRDS("uw.R50Binary.RDS")
