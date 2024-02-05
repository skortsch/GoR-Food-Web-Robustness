### metrics loop 1000 webs per year###

#
deg<-matrix(NA, 34, 1000)
no.spec<-matrix(NA, 34, 1000)
links<-matrix(NA, 34, 1000)
link.dens<-matrix(NA, 34, 1000)
deg<-matrix(NA, 34, 1000)
Con<-matrix(NA, 34, 1000)
V<-matrix(NA,34,1000)
G<-matrix(NA,34,1000)
#index.arr <- list() #could be list of output for the 1000 food webs per year

for (i in 1:1000){
  
  #i<-1
  sp.pa.year<-1*((t(biomC[,,i])[, match(row.names(un_mw), row.names(biomC[,,1]))])>0) #reorder to use for making networks (1000 per year)
  
  
  for (j in 1:34){
    fw.sel<-which(sp.pa.year[j,]==1)  
    fw<-un_mw[fw.sel, fw.sel]
    fw_g<-graph.adjacency(fw, mode = "directed")
    
    no.spec[j,i]<-dim(fw)[1] #number of species
    links[j,i]<-ecount(fw_g)#number of links
    link.dens[j,i]<-ecount(fw_g)/vcount(fw_g)#link density
    deg[j,i]<-c(mean(degree(fw_g))) # mean degree
    #generality:
    pred <- degree(fw_g, mode="in")>0
    G[j,i]<- sum(degree(fw_g, mode="in")[pred])/sum(pred)
    
    #vulnerability:
    prey <- degree(fw_g, mode="out")>0
    V[j,i] <- sum(degree(fw_g, mode="out")[prey])/sum(prey)
    
    #connectance 
    Con[j,i]<-edge_density(fw_g)

  }
}
no.spec
link.dens
links
deg
V
G
Con


w.webs.1000<-readRDS("w.webs.1000.RDS")
w.metrics<-list()

for (i in 1:1000){
  #i<-1
  a<-w.webs.1000[[i]]
  met.list<-list()
  
  for (j in 1:34){
    fw<-a[[j]]
    
    met.list[[j]]<-fluxind(fw)
    
  }
  w.metrics[[i]]<-met.list
}

w.metrics[[1000]][[34]]

w.con<-matrix(NA, 34, 1000)
qcon<-matrix(NA,34,1000)
qL<-matrix(NA,34,1000)
qLD.w<-matrix(NA,34,1000)
qLD.uw<-matrix(NA,34,1000)
qV.w<-matrix(NA,34,1000)
qG.w<-matrix(NA,34,1000)


for(i in 1:1000){
  #i<-1
  a<-w.metrics[[i]]
  
  for( j in 1:34){
    #j<-1
    b<-a[[j]]$qC.w
    w.con[j,i]<-b
    c<-a[[j]]$qC.uw
    qcon[j,i]<-c
    d<-a[[j]]$qL
    qL[j,i]<-d
    e<-a[[j]]$qLD.w
    qLD.w[j,i]<-e
    f<-a[[j]]$qLD.uw
    qLD.uw[j,i]<-f
    
    qG.w[j,i]<-a[[j]]$qG.w
    qV.w[j,i]<-a[[j]]$qV.w
  }
}

w.con
apply(w.con,1,median)

saveRDS(w.con, file="W.Con_34000.RDS")

