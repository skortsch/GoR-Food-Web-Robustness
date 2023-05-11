setwd("D:/R_Git/GoR-Food-webs/Rscripts")
w_fw<-readRDS("../Data/weighted_foodwebs.rds")
uw_fw<-readRDS("../Data/unweighted_foodwebs.rds")
uw.top.roles<-readRDS("../Data/uw.top.roles.RDS")
w.top.roles<-readRDS("../Data/w.top.roles.RDS")

library(NetworkExtinction)
library(network)
library(igraph)
library(RColorBrewer)

#test/ proof of concept code
#fw.1981<-uw_fw[[1]]

#fw.1981
#colnames(fw.1981)<-c(1:25)
#rownames(fw.1981)<-c(1:25)
#fw.1981
#fw<-graph_from_adjacency_matrix(
#  fw.1981,
#  mode = c("directed"),
#  weighted = TRUE,
#)
#deg.sp<-degree(fw, mode = c("all"),loops = TRUE,normalized = FALSE)
#order_degree<-as.numeric(names(sort(deg.sp, decreasing = T)))

#remove basal nodes from extinction order: autotrophs=2, mixotrophs=3, detritus=15
#remove<-c(2,3,which(colnames(w_fw[[i]])=="Detritus"))
#order_degree<-setdiff(order_degree, remove)

#degree högsta till lägsta 
#degree.extinc<-ExtinctionOrder(fw.1981, Order = order_degree, NetworkType = "Trophic")
#ExtinctionPlot(History = degree.extinc$sims)



###Loop code for extinction by degree###
{R50degr<-c()
extinc.plots<-list()
par(mfrow=c(4,2))
for(i in 1:length(uw_fw)){
  fw_i<-uw_fw[[i]]
  sp<-colnames(uw_fw[[i]])
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
  
  #what species go extinct? 
  b<-sort(`names<-`(deg.sp, sp), decreasing = T)
  b<-names(b)
  remove2<-c("Autotroph", "Mixotroph", "Detritus")
  b<-setdiff(b, remove2)

  #remove basal nodes from extinction order: autotrophs=2, mixotrophs=3, detritus=15
  remove<-c(2,3,which(colnames(w_fw[[i]])=="Detritus"))
  order_degree<-setdiff(order_degree, remove)

  #degree högsta till lägsta 
  degree.extinc<-ExtinctionOrder(fw_i, order_degree, NetworkType = "Trophic")
  extinc.plots[i]<-degree.extinc
  #R50
  dim(fw_i)[1]/2
  R50degr[i]<-length(which(degree.extinc$sims$TotalExt<= dim(fw_i)[1]/2))
  #R50degr[i]<-R50degr[i]/dim(fw_i)[1]
  #ExtinctionPlot(History = extinc.plots) # can work in loop with list 
  year<-c("1981 degree", "1986 degree", "1991 degree", "1996 degree", "2001 degree", "2006 degree", "2011 degree")
  plot(extinc.plots[[i]]$NumExt, extinc.plots[[i]]$AccSecExt, type = "l",main = year[i],
       xlab="Number of primary extinctions", ylab="Accumulated Secondary Extinctions", ylim = c(0,34), xlim = c(0,34))
   abline(v=R50degr[[i]])
   legend("topleft", inset = 0.1, legend= R50degr[[i]], title = "R50=", box.lty = 0, cex = 0.8)
   #legend("right", legend = b[1:R50degr[[i]]], title = "Extinction order",box.lty = 0, cex = 0.8)
   
}
plot(R50degr, type = "l", ylab= "R50", xlab = "year",xaxt = "n", main = "R50")
year2<-c("1981", "1986", "1991", "1996", "2001", "2006", "2011")
axis(1, at=1:7, labels= year2[1:7])
}

##### loop for extinctions with among module degree#
{
R50among.mod<-c()
extinc.plots<-list()
par(mfrow=c(4,2))
for(i in 1:length(uw_fw)){
  fw_i<-uw_fw[[i]]
  sp<-colnames(uw_fw[[i]])
  a<-dim(fw_i)
  colnames(fw_i)<-c(1:a[1]) # remove node names so extinction function can run
  rownames(fw_i)<-c(1:a[1])
  
  #among module degree in order
  roles<-uw.top.roles[[i]]
  roles2<-roles
  roles[,3]<-cbind(1:a[1])
  sort<-sort(roles[2:3], decreasing = T)
  order_among.mod<-as.numeric(unlist(sort[,2]))
  #remove basal nodes from extinction order: autotrophs=2, mixotrophs=3, detritus=15
  remove<-c(2,3,which(colnames(w_fw[[i]])=="Detritus"))
  order_among.mod<-setdiff(order_among.mod, remove)
  
  #what species go extinct?
  b<-sort(roles2[2:3], decreasing = T)
  remove2<-c("Autotroph", "Mixotroph", "Detritus")
  b<-as.character(unlist(b[,2]))
  b<-setdiff(b, remove2)
  
  #extinction högsta till lägsta among module degree
  among.mod.extinc<-ExtinctionOrder(fw_i, order_among.mod, NetworkType = "Trophic")
  extinc.plots[i]<-among.mod.extinc
  #ExtinctionPlot(History = extinc.plots) # did not work in loop with list, something about the format
  
  #R50
  dim(fw_i)[1]/2
  R50among.mod[i]<-length(which(among.mod.extinc$sims$TotalExt<= dim(fw_i)[1]/2))
  #R50among.mod[i]<-R50among.mod[i]/dim(fw_i)[1]
  #workaround ploting
  year<-c("1981 among module deg", "1986 among module deg", "1991 among module deg", "1996 among module deg",
          "2001 among module deg", "2006 among module deg", "2011 among module deg")
  plot(extinc.plots[[i]]$NumExt, extinc.plots[[i]]$AccSecExt, type = "l",main = year[i],
       xlab="Number of primary extinctions", ylab="Accumulated Secondary Extinctions", ylim = c(0,34), xlim = c(0,34))
  abline(v=R50among.mod[[i]])
  legend("left", legend= R50among.mod[[i]], title = "R50=", box.lty = 0, cex = 0.8)
  #legend("right", legend = b[1:R50among.mod[[i]]], box.lty = 0, cex = 0.8)
 }
plot(R50among.mod, type = "l", ylab= "R50", xlab = "year",xaxt = "n", main = "R50")
year2<-c("1981", "1986", "1991", "1996", "2001", "2006", "2011")
axis(1, at=1:7, labels= year2[1:7])
}

##### loop for extinctions with within module degree#
{R50within.mod<-c()
extinc.plots<-list()
par(mfrow=c(4,2))
for(i in 1:length(uw_fw)){
  fw_i<-uw_fw[[i]]
  sp<-colnames(uw_fw[[i]])
  a<-dim(fw_i)
  colnames(fw_i)<-c(1:a[1]) # remove node names so extinction function can run
  rownames(fw_i)<-c(1:a[1])
  
  #within module degree in order
  roles<-uw.top.roles[[i]]
  roles2<-roles
  roles[,3]<-cbind(1:a[1])
  sort<-sort(roles[1:3], decreasing = T)
  order_within.mod<-as.numeric(unlist(sort[,3]))
  #remove basal nodes from extinction order: autotrophs=2, mixotrophs=3, detritus=15
  remove<-c(2,3,which(colnames(w_fw[[i]])=="Detritus"))
  order_within.mod<-setdiff(order_within.mod, remove)
  
  #extinction högsta till lägsta within module degree
  within.mod.extinc<-ExtinctionOrder(fw_i, order_within.mod, NetworkType = "Trophic")
  extinc.plots[i]<-within.mod.extinc
  
  #what species go extinct?
  b<-sort(roles2[1:3], decreasing = T)
  remove2<-c("Autotroph", "Mixotroph", "Detritus")
  b<-as.character(unlist(b[,3]))
  b<-setdiff(b, remove2)
  
  #R50
  dim(fw_i)[1]/2
  R50within.mod[i]<-length(which(within.mod.extinc$sims$TotalExt<= dim(fw_i)[1]/2))
  #R50within.mod[i]<-R50within.mod[i]/dim(fw_i)[1]
  
  #ExtinctionPlot(History = extinc.plots) # did not work in loop with list, something about the format
  
  #workaround ploting
  year<-c("1981 within module deg", "1986 within module deg", "1991 within module deg", "1996 within module deg",
          "2001 within module deg", "2006 within module deg", "2011 within module deg")
  plot(extinc.plots[[i]]$NumExt, extinc.plots[[i]]$AccSecExt, type = "l",main = year[i],
       xlab="Number of primary extinctions", ylab="Accumulated Secondary Extinctions", ylim = c(0,34), xlim = c(0,34))
  abline(v=R50within.mod[[i]])
  legend("left", legend= R50within.mod[[i]], title = "R50=", box.lty = 0, cex = 0.8)
  #legend("right", legend = sp[order_within.mod[1:R50within.mod[[i]]]], box.lty = 0, cex = 0.8)
  }
plot(R50within.mod, type = "l", ylab= "R50", xlab = "year",xaxt = "n", main = "R50")
year2<-c("1981", "1986", "1991", "1996", "2001", "2006", "2011")
axis(1, at=1:7, labels= year2[1:7])
}



############## weighted versions ################
#fw_81<-w_fw[[1]]
#fw_81
#colnames(fw_81)<-c(1:25)
#rownames(fw_81)<-c(1:25)
#fw_81

#fw81<-graph_from_adjacency_matrix(
#  fw_81,
#  mode = c("directed"),
#  weighted = TRUE,
#)
#sum.link.weights.sp<-apply(fw_81, 1, sum)+apply(fw_81, 2, sum)
#
#order_sum_link<-as.numeric(names(sort(sum.link.weights.sp, decreasing = T)))
#remove<-c(2,3,which(colnames(w_fw[[i]])=="Detritus"))
#order_sum_link<-setdiff(order_sum_link, remove)
#
#sum.link.ext<-ExtinctionOrder(fw_81, order_sum_link, NetworkType = "Trophic")
#ExtinctionPlot(History = sum.link.ext$sims)
#plot(sum.link.ext$sims$NumExt, sum.link.ext$sims$AccSecExt, type = "l")



## loop extinction order by sum of link weights
{
R50sum.link.weight<-c()
extinc.plots<-list()
par(mfrow=c(4,2))
for(i in 1:length(w_fw)){
  fw_i<-w_fw[[i]]
  sp<-colnames(w_fw[[i]])
  a<-dim(fw_i)
  colnames(fw_i)<-c(1:a[1]) # remove node names so extinction function can run
  rownames(fw_i)<-c(1:a[1])

  fw<-graph_from_adjacency_matrix(
    fw_i,
    mode = c("directed"),
    weighted = TRUE,
  )
  sum.link.weights.sp<-apply(fw_i, 1, sum)+apply(fw_i, 2, sum)
  
  order_sum_link<-as.numeric(names(sort(sum.link.weights.sp, decreasing = T)))
  #what species go extinct?
  b<-sort(`names<-`(sum.link.weights.sp, sp), decreasing = T)
  b<-names(b)
  remove2<-c("Autotroph", "Mixotroph", "Detritus")
  b<-setdiff(b, remove2)
  #remove basal nodes from extinction order: autotrophs=2, mixotrophs=3, detritus=15
  remove<-c(2,3,which(colnames(w_fw[[i]])=="Detritus"))
  order_sum_link<-setdiff(order_sum_link, remove)
  
  #degree högsta till lägsta 
  sum.link.extinc<-ExtinctionOrder(fw_i, order_sum_link, NetworkType = "Trophic")
  extinc.plots[i]<-sum.link.extinc
  
  #R50
  dim(fw_i)[1]/2
  R50sum.link.weight[i]<-length(which(sum.link.extinc$sims$TotalExt<  dim(fw_i)[1]/2))
  #ExtinctionPlot(History = extinc.plots) # can work in loop with list 
  year<-c("w1981 sum of linkweights", "w1986 sum of linkweights", "w1991 sum of linkweights", "w1996 sum of linkweights",
          "w2001 sum of linkweights", "w2006 sum of linkweights", "w2011 sum of linkweights")
  plot(extinc.plots[[i]]$NumExt, extinc.plots[[i]]$AccSecExt, type = "l",main = year[i],
       xlab="Number of extinctions", ylab="Accumulated Secondary Extinctions", ylim= c(0,34), xlim = c(0,34))
  abline(v= R50sum.link.weight[[i]])
  legend("left", legend= R50sum.link.weight[[i]], title = "R50=", box.lty = 0, cex = 0.8)
  #legend("right", legend = b[1:R50sum.link.weight[[i]]], box.lty = 0, cex = 0.8)
}
plot(R50sum.link.weight, type = "l", ylab= "R50", xlab = "year",xaxt = "n", main = "R50")
year2<-c("1981", "1986", "1991", "1996", "2001", "2006", "2011")
axis(1, at=1:7, labels= year2[1:7])
}


########## extinction order by Weighted among mod degree
{
R50among.mod.w<-c()
extinc.plots<-list()
par(mfrow=c(4,2))
for(i in 1:length(w_fw)){
  fw_i<-w_fw[[i]]
  sp<-colnames(w_fw[[i]])
  a<-dim(fw_i)
  colnames(fw_i)<-c(1:a[1]) # remove node names so extinction function can run
  rownames(fw_i)<-c(1:a[1])
  
  #among module degree in order
  roles<-w.top.roles[[i]]
  roles[,3]<-cbind(1:a[1])
  sort<-sort(roles[2:3], decreasing = T)
  w.order_among.mod<-as.numeric(unlist(sort[,2]))
  #remove basal nodes from extinction order: autotrophs=2, mixotrophs=3, detritus=
  remove<-c(2,3,which(colnames(w_fw[[i]])=="Detritus"))
  w.order_among.mod<-setdiff(w.order_among.mod, remove)
  
  #what species go extinct?
  b<-sort(roles2[2:3], decreasing = T)
  remove2<-c("Autotroph", "Mixotroph", "Detritus")
  b<-as.character(unlist(b[,2]))
  b<-setdiff(b, remove2)
  
  #extinction högsta till lägsta among module degree
  w.among.mod.extinc<-ExtinctionOrder(fw_i, w.order_among.mod, NetworkType = "Trophic")
  extinc.plots[i]<-w.among.mod.extinc
  #ExtinctionPlot(History = extinc.plots) # did not work in loop with list, something about the format
  
  #R50
  R50among.mod.w[i]<-length(which(w.among.mod.extinc$sims$TotalExt<  dim(fw_i)[1]/2))
  #workaround ploting
  year<-c("w1981 among module deg", "w1986 among module deg", "w1991 among module deg", "w1996 among module deg",
          "w2001 among module deg", "w2006 among module deg", "w2011 among module deg")
  plot(extinc.plots[[i]]$NumExt, extinc.plots[[i]]$AccSecExt, type = "l",main = year[i],
       xlab="Number of extinctions", ylab="Accumulated Secondary Extinctions", ylim = c(0,34), xlim = c(0,34))
  abline(v= R50among.mod.w[[i]])
  legend("left", legend= R50among.mod.w[[i]], title = "R50=", box.lty = 0, cex = 0.8)
  #legend("right", legend = b[1:R50among.mod.w[[i]]], box.lty = 0, cex = 0.8)
}
plot(R50among.mod.w, type = "l", ylab= "R50", xlab = "year",xaxt = "n", main = "R50")
year2<-c("1981", "1986", "1991", "1996", "2001", "2006", "2011")
axis(1, at=1:7, labels= year2[1:7])
}

##### weighted version extinctions order by within module degree#
{
R50within.mod.w<-c()
extinc.plots<-list()
par(mfrow=c(4,2))
for(i in 1:length(w_fw)){
  fw_i<-w_fw[[i]]
  sp<-colnames(w_fw[[i]])
  a<-dim(fw_i)
  colnames(fw_i)<-c(1:a[1]) # remove node names so extinction function can run
  rownames(fw_i)<-c(1:a[1])
  
  #within module degree in order
  roles<-w.top.roles[[i]]
  roles2<-roles
  roles[,3]<-cbind(1:a[1])
  sort<-sort(roles[1:3], decreasing = T)
  w.order_within.mod<-as.numeric(unlist(sort[,3]))
  #remove basal nodes from extinction order: autotrophs=2, mixotrophs=3, detritus=15
  remove<-c(2,3,which(colnames(w_fw[[i]])=="Detritus"))
  w.order_within.mod<-setdiff(w.order_within.mod, remove)
  
  #what species go extinct?
  b<-sort(roles2[1:3], decreasing = T)
  remove2<-c("Autotroph", "Mixotroph", "Detritus")
  b<-as.character(unlist(b[,3]))
  b<-setdiff(b, remove2)
  
  #extinction högsta till lägsta within module degree
  w.within.mod.extinc<-ExtinctionOrder(fw_i, Order = w.order_within.mod, NetworkType = "Trophic")
  extinc.plots[i]<-w.within.mod.extinc
  #ExtinctionPlot(History = extinc.plots) # did not work in loop with list, something about the format
  
  #R50
  R50within.mod.w[i]<-length(which(w.within.mod.extinc$sims$TotalExt<  dim(fw_i)[1]/2))
  #workaround ploting
  year<-c("w1981 within module deg", "w1986 within module deg", "w1991 within module deg", "w1996 within module deg",
          "w2001 within module deg", "w2006 within module deg", "w2011 within module deg")
  plot(extinc.plots[[i]]$NumExt, extinc.plots[[i]]$AccSecExt, type = "l",main = year[i],
       xlab="Number of extinctions", ylab="Accumulated Secondary Extinctions", ylim = c(0,34), xlim = c(0,34))
  abline(v= R50within.mod.w[[i]])
  legend("left", legend= R50within.mod.w[[i]], title = "R50=", box.lty = 0, cex = 0.8)
  #legend("right", legend = b[1:w.order_within.mod[[i]]], box.lty = 0, cex = 0.8)
}
plot(R50within.mod.w, type = "l", ylab= "R50", xlab = "year",xaxt = "n", main = "R50")
year2<-c("1981", "1986", "1991", "1996", "2001", "2006", "2011")
axis(1, at=1:7, labels= year2[1:7])
}


### weighted versions with threshold ####

###############loop extinction order by sum of link weights##############
{
  R50sum.link.weight<-c()
  extinc.plots<-list()
  par(mfrow=c(4,2))
  for(i in 1:length(w_fw)){
    #y<-graph_from_adjacency_matrix(w_fw[[4]], weighted = T, mode = "directed")
    #fw_i<-as.matrix(y,  matrix.type = "adjacency", attr= "weight", sparse= F)
    fw_i<-w_fw[[i]]
    sp<-colnames(w_fw[[i]])
    a<-dim(fw_i)
    colnames(fw_i)<-c(1:a[1]) # remove node names so extinction function can run
    rownames(fw_i)<-c(1:a[1])
    
    fw<-graph_from_adjacency_matrix(
      fw_i,
      mode = c("directed"),
      weighted = TRUE,
    )
    sum.link.weights.sp<-apply(fw_i, 1, sum)+apply(fw_i, 2, sum)
    
    order_sum_link<-as.numeric(names(sort(sum.link.weights.sp, decreasing = T)))
  
    #what species go extinct
    b<-sort(`names<-`(sum.link.weights.sp, sp), decreasing = T)
    b<-names(b)
    remove2<-c("Autotroph", "Mixotroph", "Detritus")
    b<-setdiff(b, remove2)
    
    #remove basal nodes from extinction order: autotrophs=2, mixotrophs=3, detritus=which(colnames(w_fw[[i]])=="Detritus"))
    remove<-c(2,3,which(colnames(w_fw[[i]])=="Detritus"))
    order_sum_link<-setdiff(order_sum_link, remove)
    
    #degree högsta till lägsta 
    #network(w_fw[[1]], ignore.eval=FALSE,names.eval='weight')
    #fw_i<-graph.adjacency(fw_i, mode="directed", weighted=T)
    #fw_i<-as.matrix(fw_i, matrix.type = "adjacency", attr= "weight", sparse= F)
    #fw_i<-as.network(fw_i, matrix.type = "adjacency", attrnames= "weight",names.eval = "weight", ignore.eval=FALSE)
    sum.link.extinc<-ExtinctionOrder(fw_i, order_sum_link, NetworkType = "Trophic", IS=0.3)
    extinc.plots[i]<-sum.link.extinc
    
    #R50
    dim(fw_i)[1]/2
    R50sum.link.weight[i]<-length(which(sum.link.extinc$sims$TotalExt<  dim(fw_i)[1]/2))
    #ExtinctionPlot(History = extinc.plots) # can work in loop with list 
    year<-c("w1981 sum.lw, TH 30%", "w1986 sum.lw, Treshold 30%", "w1991 sum.lw, TH 30%", "w1996 sum.lw, 30%",
            "w2001 sum.lw, TH 30%", "w2006 sum.lw, TH 30%", "w2011 sum.lw, TH 30%")
    plot(extinc.plots[[i]]$NumExt, extinc.plots[[i]]$AccSecExt, type = "l",main = year[i],
         xlab="Number of extinctions", ylab="Accumulated Secondary Extinctions", ylim = c(0,34), xlim = c(0,34))
    abline(v= R50sum.link.weight[[i]])
    legend("left", legend= R50sum.link.weight[[i]], title = "R50=", box.lty = 0, cex = 0.8)
    #legend("right", legend = b[1:10], box.lty = 0, cex = 0.8)
  }
  plot(R50sum.link.weight, type = "l", ylab= "R50", xlab = "year",xaxt = "n", main = "R50")
  year2<-c("1981", "1986", "1991", "1996", "2001", "2006", "2011")
  axis(1, at=1:7, labels= year2[1:7])
}
#as.network(fw_i, matrix.type = "adjacency", attrname="weight",ignore.eval=FALSE, names.eval='weight')

########## extinction order by Weighted among mod degree , TH 30% ###############
{
  R50among.mod.w<-c()
  extinc.plots<-list()
  par(mfrow=c(4,2))
  for(i in 1:length(w_fw)){
    fw_i<-w_fw[[i]]
    sp<-colnames(w_fw[[i]])
    a<-dim(fw_i)
    colnames(fw_i)<-c(1:a[1]) # remove node names so extinction function can run
    rownames(fw_i)<-c(1:a[1])
    
    #among module degree in order
    roles<-w.top.roles[[i]]
    roles2<-roles
    roles[,3]<-cbind(1:a[1])
    sort<-sort(roles[2:3], decreasing = T)
    w.order_among.mod<-as.numeric(unlist(sort[,2]))
    #remove basal nodes from extinction order: autotrophs=2, mixotrophs=3, detritus=15
    remove<-c(2,3,15)
    w.order_among.mod<-setdiff(w.order_among.mod, remove)
    
    #what species go extinct?
    b<-sort(roles2[2:3], decreasing = T)
    remove2<-c("Autotroph", "Mixotroph", "Detritus")
    b<-as.character(unlist(b[,2]))
    b<-setdiff(b, remove2)
    
    #extinction högsta till lägsta among module degree
    w.among.mod.extinc<-ExtinctionOrder(fw_i, w.order_among.mod[1:18], NetworkType = "Trophic", IS=0.3)
    extinc.plots[i]<-w.among.mod.extinc
    #ExtinctionPlot(History = extinc.plots) # did not work in loop with list, something about the format
    
    #R50
    R50among.mod.w[i]<-length(which(w.among.mod.extinc$sims$TotalExt<= dim(fw_i)[1]/2))
    #workaround ploting
    year<-c("w1981 am.mod.deg, TH 30%", "w1986 am.mod.deg, TH 30%", "w1991 am.mod.deg, TH 30%", "w1996 am.mod.deg, TH 30%",
            "w2001 am.mod.deg, TH 30%", "w2006 am.mod.deg, TH 30%", "w2011 am.mod.deg, TH 30%")
    plot(extinc.plots[[i]]$NumExt, extinc.plots[[i]]$AccSecExt, type = "l",main = year[i],
         xlab="Number of extinctions", ylab="Accumulated Secondary Extinctions", ylim = c(0,34), xlim = c(0,34))
    abline(v= R50among.mod.w[[i]])
    legend("left", legend= R50among.mod.w[[i]], title = "R50=", box.lty = 0, cex = 0.8)
    #legend("right", legend = b[1:R50among.mod.w[[i]]], box.lty = 0, cex = 0.8)
  }
  plot(R50among.mod.w, type = "l", ylab= "R50", xlab = "year",xaxt = "n", main = "R50")
  year2<-c("1981", "1986", "1991", "1996", "2001", "2006", "2011")
  axis(1, at=1:7, labels= year2[1:7]) 
}

##### weighted version extinctions order by within module degree#####
{
  R50within.mod.w<-c()
  extinc.plots<-list()
  par(mfrow=c(4,2))
  for(i in 1:length(w_fw)){
    fw_i<-w_fw[[i]]
    sp<-colnames(w_fw[[i]])
    a<-dim(fw_i)
    colnames(fw_i)<-c(1:a[1]) # remove node names so extinction function can run
    rownames(fw_i)<-c(1:a[1])
    
    #within module degree in order
    roles<-w.top.roles[[i]]
    roles2<-roles
    roles[,3]<-cbind(1:a[1])
    sort<-sort(roles[1:3], decreasing = T)
    w.order_within.mod<-as.numeric(unlist(sort[,3]))
    #remove basal nodes from extinction order: autotrophs=2, mixotrophs=3, detritus=which(colnames(w_fw[[i]])=="Detritus")
    remove<-c(2,3,which(colnames(w_fw[[i]])=="Detritus"))
    w.order_within.mod<-setdiff(w.order_within.mod, remove)
    
    #what species go extinct?
    b<-sort(roles2[1:3], decreasing = T)
    remove2<-c("Autotroph", "Mixotroph", "Detritus")
    b<-as.character(unlist(b[,3]))
    b<-setdiff(b, remove2)
    
    #extinction högsta till lägsta within module degree
    w.within.mod.extinc<-ExtinctionOrder(fw_i, w.order_within.mod[1:18], NetworkType = "Trophic", IS=0.3)
    extinc.plots[i]<-w.within.mod.extinc
    #ExtinctionPlot(History = extinc.plots) # did not work in loop with list, something about the format
    
    #R50
    R50within.mod.w[i]<-length(which(w.within.mod.extinc$sims$TotalExt<  dim(fw_i)[1]/2))
    #workaround ploting
    year<-c("w1981 wit.mod.deg, TH 30%", "w1986  wit.mod.deg, TH 30%", "w1991  wit.mod.deg, TH 30%", "w1996  wit.mod.deg, TH 30%",
            "w2001  wit.mod.deg, TH 30%", "w2006  wit.mod.deg, TH 30%", "w2011  wit.mod.deg, TH 30%")
    plot(extinc.plots[[i]]$NumExt, extinc.plots[[i]]$AccSecExt, type = "l",main = year[i],
         xlab="Number of extinctions", ylab="Accumulated Secondary Extinctions", ylim= c(0,34), xlim = c(0,34))
    abline(v= R50within.mod.w[[i]])
    legend("left", legend= R50within.mod.w[[i]], title = "R50=", box.lty = 0, cex = 0.8)
    #legend("right", legend = b[1:R50within.mod.w[[i]]], box.lty = 0, cex = 0.8)
  }
  plot(R50within.mod.w, type = "l", ylab= "R50", xlab = "year",xaxt = "n", main = "R50")
  year2<-c("1981", "1986", "1991", "1996", "2001", "2006", "2011")
  axis(1, at=1:7, labels= year2[1:7]) 
}


#loop for all years, TH 30%
{
  R50sum.link.weight<-list()
  R50index<-list()
  extinc.plots<-list()
  #par(mfrow=c(4,2))
  for(i in 1:length(w_fw)){
    #y<-graph_from_adjacency_matrix(w_fw[[4]], weighted = T, mode = "directed")
    #fw_i<-as.matrix(y,  matrix.type = "adjacency", attr= "weight", sparse= F)
    fw_i<-w_fw[[i]]
    sp<-colnames(w_fw[[i]])
    a<-dim(fw_i)
    colnames(fw_i)<-c(1:a[1]) # remove node names so extinction function can run
    rownames(fw_i)<-c(1:a[1])
    
    fw<-graph_from_adjacency_matrix(
      fw_i,
      mode = c("directed"),
      weighted = TRUE,
    )
    sum.link.weights.sp<-apply(fw_i, 1, sum)+apply(fw_i, 2, sum)
    
    order_sum_link<-as.numeric(names(sort(sum.link.weights.sp, decreasing = T)))
    
     #what species go extinct
     b<-sort(`names<-`(sum.link.weights.sp, sp), decreasing = T)
     b<-names(b)
     remove2<-c("Autotroph", "Mixotroph", "Detritus")
     b<-setdiff(b, remove2)
    
    #remove basal nodes from extinction order: autotrophs=2, mixotrophs=3, detritus=which(colnames(w_fw[[i]])=="Detritus"))
    remove<-c(2,3,which(colnames(w_fw[[i]])=="Detritus"))
    order_sum_link<-setdiff(order_sum_link, remove)
    
    sum.link.extinc<-ExtinctionOrder(fw_i, order_sum_link, NetworkType = "Trophic", IS=0.3) #0.3 leads to removal of all nodes which lose 70 percent of their interaction strength in the Network argument)
    extinc.plots[[i]]<-sum.link.extinc
    R50index[i]<-sum.link.extinc$R50
    #ExtinctionPlot(History=extinc.plots[[i]]$sims, Variable = "AccSecExt")
    #R50
    R50sum.link.weight[[i]]<-length(which(sum.link.extinc$sims$TotalExt<  dim(fw_i)[1]/2))
  }
  
  year<-c("1981 ", "1986", "1991 ", "1996 ", "2001 ", "2006 ", "2011 ")
  par(bg = "grey")
  plot(extinc.plots[[1]]$sims$NumExt, extinc.plots[[1]]$sims$AccSecExt, type = "l",main = "Extinction by sum.link.weigth, loss of 70% incoming flow = extinction",
       xlab="Number of primary extinctions", ylab="Secondary Extinctions", ylim = c(0,20), xlim = c(0,20))
  #abline(v= R50multiTH)
  #legend("topleft", inset = 0.1, legend= R50.WMW, title = "R50=", box.lty = 0)
  #legend("right", legend =  b[1:R50.WMW], title = "Extinction order:", box.lty = 0)
  #plotcol <- c("red","blue","steelblue","black", "cyan3","orange","purple")
  plotcol<-brewer.pal(name = "Spectral", n=length(extinc.plots))
  for (j in 1:length(extinc.plots)){
    lines(extinc.plots[[j]]$sims$NumExt, extinc.plots[[j]]$sims$AccSecExt, col=plotcol[j], lwd=2)
    points(x=extinc.plots[[j]]$sims$NumExt, y=extinc.plots[[j]]$sims$AccSecExt, col=plotcol[j], pch=15, type = "p")
    # abline(v= R50multiTH[j], col=j)
   
  }
  legend("bottomright", legend= R50sum.link.weight, title = "Nr prim.ext. R50", col = plotcol, pch = 15)
  legend("topright", legend = year, title = "Year", col = plotcol, pch = 15)
  legend("topleft", legend= R50index, title = "R50", col = plotcol, pch = 15)
  #legend("topright", legend = TH*100, title = "% left of incoming flow
  #     required to trigger extinction", col = plotcol, pch = 15, box.lty = 0)
  
}

plot(extinc.plots[[1]]$sims$NumExt, extinc.plots[[1]]$sims$AccSecExt, type = "l",main = "Extinction by sum.link.weigth, loss of 70% incoming flow = extinction",
     xlab="Number of primary extinctions", ylab="Secondary Extinctions", ylim = c(0,16), xlim = c(0,16))
lines(extinc.plots[[1]]$sims$NumExt, extinc.plots[[1]]$sims$AccSecExt, col=plotcol[1], lwd=2)
lines(extinc.plots[[7]]$sims$NumExt, extinc.plots[[7]]$sims$AccSecExt, col=plotcol[7], lwd=2)
points(x=extinc.plots[[7]]$sims$NumExt, y=extinc.plots[[7]]$sims$AccSecExt, col=plotcol[7], pch=15, type = "p")
points(x=extinc.plots[[1]]$sims$NumExt, y=extinc.plots[[1]]$sims$AccSecExt, col=plotcol[1], pch=15, type = "p")
legend("bottomright", legend= R50sum.link.weight[c(1,7)], title = "Nr of primary extinctions to reach R50", col = plotcol[c(1,7)], pch = 15)
legend("topright", legend = year[c(1,7)], title = "Year", col = plotcol[c(1,7)], pch = 15)
legend("topleft", legend= R50index[c(1,7)], title = "R50", col = plotcol[c(1,7)], pch = 15)






{
  R50sum.link.weight<-list()
  R50index<-list()
  extinc.plots<-list()
  #par(mfrow=c(4,2))
  for(i in 1:length(w_fw)){
    #y<-graph_from_adjacency_matrix(w_fw[[4]], weighted = T, mode = "directed")
    #fw_i<-as.matrix(y,  matrix.type = "adjacency", attr= "weight", sparse= F)
    fw_i<-w_fw[[i]]
    sp<-colnames(w_fw[[i]])
    a<-dim(fw_i)
    colnames(fw_i)<-c(1:a[1]) # remove node names so extinction function can run
    rownames(fw_i)<-c(1:a[1])
    
    fw<-graph_from_adjacency_matrix(
      fw_i,
      mode = c("directed"),
      weighted = TRUE,
    )
    sum.link.weights.sp<-apply(fw_i, 1, sum)+apply(fw_i, 2, sum)
    
    order_sum_link<-as.numeric(names(sort(sum.link.weights.sp, decreasing = T)))
    
    #what species go extinct
    b<-sort(`names<-`(sum.link.weights.sp, sp), decreasing = T)
    b<-names(b)
    remove2<-c("Autotroph", "Mixotroph", "Detritus")
    b<-setdiff(b, remove2)
    
    #remove basal nodes from extinction order: autotrophs=2, mixotrophs=3, detritus=which(colnames(w_fw[[i]])=="Detritus"))
    remove<-c(2,3,which(colnames(w_fw[[i]])=="Detritus"))
    order_sum_link<-setdiff(order_sum_link, remove)
    
    sum.link.extinc<-ExtinctionOrder(fw_i, order_sum_link, NetworkType = "Trophic", IS=0.1) #0.3 leads to removal of all nodes which lose 70 percent of their interaction strength in the Network argument)
    extinc.plots[[i]]<-sum.link.extinc
    R50index[i]<-sum.link.extinc$R50
    #ExtinctionPlot(History=extinc.plots[[i]]$sims, Variable = "AccSecExt")
    #R50
    R50sum.link.weight[[i]]<-length(which(sum.link.extinc$sims$TotalExt<  dim(fw_i)[1]/2))
  }
  
  year<-c("1981 ", "1986", "1991 ", "1996 ", "2001 ", "2006 ", "2011 ")
  par(bg = "grey")
  plot(extinc.plots[[1]]$sims$NumExt, extinc.plots[[1]]$sims$AccSecExt, type = "l",main = "Extinction by sum.link.weigth, loss of 90% incoming flow = extinction",
       xlab="Number of primary extinctions", ylab="Secondary Extinctions", ylim = c(0,20), xlim = c(0,20))
  #abline(v= R50multiTH)
  #legend("topleft", inset = 0.1, legend= R50.WMW, title = "R50=", box.lty = 0)
  #legend("right", legend =  b[1:R50.WMW], title = "Extinction order:", box.lty = 0)
  #plotcol <- c("red","blue","steelblue","black", "cyan3","orange","purple")
  plotcol<-brewer.pal(name = "Spectral", n=length(extinc.plots))
  for (j in 1:length(extinc.plots)){
    lines(extinc.plots[[j]]$sims$NumExt, extinc.plots[[j]]$sims$AccSecExt, col=plotcol[j], lwd=2)
    points(x=extinc.plots[[j]]$sims$NumExt, y=extinc.plots[[j]]$sims$AccSecExt, col=plotcol[j], pch=15, type = "p")
    # abline(v= R50multiTH[j], col=j)
    
  }
  legend("bottomright", legend= R50sum.link.weight, title = "Nr of prim.ext.R50", col = plotcol, pch = 15)
  legend("topright", legend = year, title = "Year", col = plotcol, pch = 15)
  legend("topleft", legend= R50index, title = "R50", col = plotcol, pch = 15)
  #legend("topright", legend = TH*100, title = "% left of incoming flow
  #     required to trigger extinction", col = plotcol, pch = 15, box.lty = 0)
  
}

plot(extinc.plots[[1]]$sims$NumExt, extinc.plots[[1]]$sims$AccSecExt, type = "l",main = "Extinction by sum.link.weigth, loss of % incoming flow = extinction",
     xlab="Number of primary extinctions", ylab="Secondary Extinctions", ylim = c(0,16), xlim = c(0,16))
lines(extinc.plots[[1]]$sims$NumExt, extinc.plots[[1]]$sims$AccSecExt, col=plotcol[1], lwd=2)
lines(extinc.plots[[7]]$sims$NumExt, extinc.plots[[7]]$sims$AccSecExt, col=plotcol[7], lwd=2)
points(x=extinc.plots[[7]]$sims$NumExt, y=extinc.plots[[7]]$sims$AccSecExt, col=plotcol[7], pch=15, type = "p")
points(x=extinc.plots[[1]]$sims$NumExt, y=extinc.plots[[1]]$sims$AccSecExt, col=plotcol[1], pch=15, type = "p")
legend("bottomright", legend= R50sum.link.weight[c(1,7)], title = "Nr of primary extinctions to reach R50", col = plotcol[c(1,7)], pch = 15)
legend("topright", legend = year[c(1,7)], title = "Year", col = plotcol[c(1,7)], pch = 15)
legend("topleft", legend= R50index[c(1,7)], title = "R50", col = plotcol[c(1,7)], pch = 15)

