#Metaweb loop

setwd("D:/R_Git/GoR-Food-webs/Rscripts")
Metaweb_GOR <- readRDS("../Data/Metaweb_GOR.RDS")
LW_Metaweb_GOR <- readRDS("../Data/LW_Metaweb_GOR.RDS")
Meta.top.roles<-readRDS("../Data/Metaweb_top_roles.RDS")
load("../Data/metawebs_GoR.RData")

#install.packages("NetworkExtinction")
library(NetworkExtinction)
library(network)
library(igraph)
library(RColorBrewer)

metaweb<-as.matrix(Metaweb_GOR, matrix.type = "adjacency", sparse= F)

#MW extinction by degree
{
fw_i<-metaweb
sp<-colnames(metaweb)
a<-dim(fw_i)
colnames(fw_i)<-c(1:a[1])# remove node names so extinction function can run
rownames(fw_i)<-c(1:a[1])

#degree of nodes
fw<-graph_from_adjacency_matrix(
  fw_i,
  mode = c("directed"),
  weighted = T
)
deg.sp<-degree(fw, mode = c("all"),loops = TRUE,normalized = FALSE)

order_degree<-as.numeric(names(sort(deg.sp, decreasing = T)))

#what species go extinct? 
b<-sort(`names<-`(deg.sp, sp), decreasing = T)
b<-names(b)
remove2<-c("Autotroph", "Mixotroph", "Detritus")
b<-setdiff(b, remove2)

#remove basal nodes from extinction order: autotrophs=2, mixotrophs=3, detritus=20
remove<-c(2,3,20)
order_degree<-setdiff(order_degree, remove)

#degree högsta till lägsta 
degree.extinc.mw<-ExtinctionOrder(fw_i, Order = order_degree, NetworkType = "Trophic", IS=0.3)
#R50
dim(fw_i)[1]/2
R50.MW<-length(which(degree.extinc.mw$sims$TotalExt<= dim(fw_i)[1]/2))

plot(degree.extinc.mw$sims$NumExt, degree.extinc.mw$sims$AccSecExt, type = "l",main = "Metaweb unweighted, degree",
     xlab="Number of primary extinctions", ylab="Accumulated Secondary Extinctions", ylim = c(0,34), xlim = c(0,34))
abline(v=R50.MW)
legend("topleft", inset = 0.1, legend= R50.MW, title = "R50=", box.lty = 0)
legend("right", legend=b[1:R50.MW], title="Extinction order", box.lty=0)
}


# metaweb extinction uw, by among mod degree
{
    fw_i<-metaweb
    a<-dim(fw_i)
    colnames(fw_i)<-c(1:a[1]) # remove node names so extinction function can run
    rownames(fw_i)<-c(1:a[1])
    
    #among module degree in order
    roles<-Meta.top.roles[[1]]
    roles2<-roles
    roles[,3]<-cbind(1:a[1])
    sort<-sort(roles[2:3], decreasing = T)
    order_among.mod<-as.numeric(unlist(sort[,2]))
    #remove basal nodes from extinction order: autotrophs=2, mixotrophs=3, detritus=20
    remove<-c(2,3,20)
    order_among.mod<-setdiff(order_among.mod, remove)
    
    #what species go extinct?
    b<-sort(roles2[2:3], decreasing = T)
    remove2<-c("Autotroph", "Mixotroph", "Detritus")
    b<-as.character(unlist(b[,2]))
    b<-setdiff(b, remove2)
    #extinction högsta till lägsta among module degree
    among.mod.extinc<-ExtinctionOrder(fw_i, Order = order_among.mod, NetworkType = "Trophic")
    extinc.plots<-among.mod.extinc
    #ExtinctionPlot(History = extinc.plots) # did not work in loop with list, something about the format
    
    #R50
    dim(fw_i)[1]/2
    R50among.mod.mw<-length(which(among.mod.extinc$sims$TotalExt<= dim(fw_i)[1]/2))
    #R50among.mod[i]<-R50among.mod[i]/dim(fw_i)[1]
    #workaround ploting
    plot(extinc.plots$sims$NumExt, extinc.plots$sims$AccSecExt, type = "l",main = "Metaweb unweighted, among.mod.deg",
         xlab="Number of primary extinctions", ylab="Accumulated Secondary Extinctions", ylim = c(0,34), xlim = c(0,34))
    abline(v=R50among.mod.mw)
    legend("topleft", inset = 0.1, legend= R50among.mod.mw, title = "R50=", box.lty = 0)
    legend("right", legend = b[1:R50among.mod.mw], title = "Extinction order:", box.lty = 0)
  }
 

# Metaweb extinctions, by within module degree
{
  fw_i<-metaweb
  a<-dim(fw_i)
  colnames(fw_i)<-c(1:a[1]) # remove node names so extinction function can run
  rownames(fw_i)<-c(1:a[1])
  
  #among module degree in order
  roles<-Meta.top.roles[[1]]
  roles2<-roles
  roles[,3]<-cbind(1:a[1])
  sort<-sort(roles[1:3], decreasing = T)
  order_within.mod<-as.numeric(unlist(sort[,3]))
  #remove basal nodes from extinction order: autotrophs=1, mixotrophs=2, detritus=20
  remove<-c(2,3,20)
  order_within.mod<-setdiff(order_within.mod, remove)
  b<-sort(roles2[1:3], decreasing = T)
  remove2<-c("Autotroph", "Mixotroph", "Detritus")
  b<-as.character(unlist(b[,3]))
  b<-setdiff(b, remove2)
  #extinction högsta till lägsta among module degree
  within.mod.extinc<-ExtinctionOrder(fw_i, Order = order_within.mod, NetworkType = "Trophic")
  extinc.plots<-within.mod.extinc
  #ExtinctionPlot(History = extinc.plots) # did not work in loop with list, something about the format
  
  #R50
  dim(fw_i)[1]/2
  R50within.mod.mw<-length(which(within.mod.extinc$sims$TotalExt<= dim(fw_i)[1]/2))
  #R50among.mod[i]<-R50among.mod[i]/dim(fw_i)[1]
  #workaround ploting
  plot(extinc.plots$sims$NumExt, extinc.plots$sims$AccSecExt, type = "l",main = "Metaweb unweighted, within.mod.deg",
       xlab="Number of primary extinctions", ylab="Accumulated Secondary Extinctions", ylim = c(0,34), xlim = c(0,34))
  abline(v=R50within.mod.mw)
  legend("topleft", inset = 0.1, legend= R50within.mod.mw, title = "R50=", box.lty = 0)
  legend("right", legend = b[1:R50within.mod.mw], title = "Extinction order:", box.lty = 0)
}





#weighted MW
LW_Metaweb_GOR

# Weighted MW extinction by sum of link weights
{
  fw_i<-as.matrix(LW_Metaweb_GOR, matrix.type = "adjacency", attr= "weight", sparse= F)
  sp<-colnames(metaweb)
  a<-dim(fw_i)
  colnames(fw_i)<-c(1:a[1])# remove node names so extinction function can run
  rownames(fw_i)<-c(1:a[1])
  
  #sum of of links
  fw<-graph_from_adjacency_matrix(
    fw_i,
    mode = c("directed"),
    weighted = T
  )
  sum.link.weights.sp<-apply(fw_i, 1, sum)+apply(fw_i, 2, sum)
  
  order_sum_link<-as.numeric(names(sort(sum.link.weights.sp, decreasing = T)))
  
  #what species go extinct
  b<-sort(`names<-`(sum.link.weights.sp, sp), decreasing = T)
  b<-names(b)
  remove2<-c("Autotroph", "Mixotroph", "Detritus")
  b<-setdiff(b, remove2)
  
  #remove basal nodes from extinction order: autotrophs=2, mixotrophs=3, detritus=20
  remove<-c(2,3,20)
  order_sum_link<-setdiff(order_sum_link, remove)
  #sum links
  sum.link.ext<-ExtinctionOrder(fw_i, order_sum_link[1:24], NetworkType = "Trophic")
  ExtinctionPlot(History = sum.link.ext$sims)
  #R50
  dim(fw_i)[1]/2
  R50.WMW<-length(which(sum.link.ext$sims$TotalExt<= dim(fw_i)[1]/2))
  #ploting
  plot(sum.link.ext$sims$NumExt, sum.link.ext$sims$AccSecExt, type = "l",main = "Weighted metaweb, sum linkweights",
       xlab="Number of primary extinctions", ylab="Accumulated Secondary Extinctions", ylim = c(0,34), xlim = c(0,34))
  abline(v= R50.WMW)
  legend("topleft", inset = 0.1, legend= R50.WMW, title = "R50=", box.lty = 0)
  legend("right", legend = b[1:R50.WMW], title = "Extinction order:", box.lty = 0)
}



# Weighted mw extinction by sum of link weights TH 30%
{
  fw_i<-as.matrix(LW_Metaweb_GOR, matrix.type = "adjacency", attr= "weight", sparse= F)
  sp<-colnames(metaweb)
  a<-dim(fw_i)
  colnames(fw_i)<-c(1:a[1])# remove node names so extinction function can run
  rownames(fw_i)<-c(1:a[1])
  
  #sum of of links
  fw<-graph_from_adjacency_matrix(
    fw_i,
    mode = c("directed"),
    weighted = T
  )
  sum.link.weights.sp<-apply(fw_i, 1, sum)+apply(fw_i, 2, sum)
  
  order_sum_link<-as.numeric(names(sort(sum.link.weights.sp, decreasing = T)))
  
  #what species go extinct
  b<-sort(`names<-`(sum.link.weights.sp, sp), decreasing = T)
  b<-names(b)
  remove2<-c("Autotroph", "Mixotroph", "Detritus")
  b<-setdiff(b, remove2)
  
  #remove basal nodes from extinction order: autotrophs=2, mixotrophs=3, detritus=20
  remove<-c(2,3,20)
  order_sum_link<-setdiff(order_sum_link, remove)
  
  #extinction sum links
  sum.link.ext<-ExtinctionOrder(fw_i, Order=order_sum_link[1:24], NetworkType = "Trophic", IS=0.3)
  ExtinctionPlot(History = sum.link.ext$sims)
  #R50
  dim(fw_i)[1]/2
  R50.WMW<-length(which(sum.link.ext$sims$TotalExt<= dim(fw_i)[1]/2))
  #what species go extinct?
  sp[order_sum_link[1:R50.WMW]]
  #ploting
  plot(sum.link.ext$sims$NumExt, sum.link.ext$sims$AccSecExt, type = "l",main = "Weighted metaweb, extinction by sum.link.weigth, TH 30%",
       xlab="Number of primary extinctions", ylab="Accumulated Secondary Extinctions", ylim = c(0,34), xlim = c(0,34))
  abline(v= R50.WMW)
  legend("topleft", inset = 0.1, legend= R50.WMW, title = "R50=", box.lty = 0)
  legend("right", legend =  b[1:R50.WMW], title = "Extinction order:", box.lty = 0)
}




#weighted metaweb extinction by among module degree
{
  fw_i<-as.matrix(LW_Metaweb_GOR, matrix.type = "adjacency", attr= "weight", sparse= F)
  sp<-colnames(metaweb)
  a<-dim(fw_i)
  colnames(fw_i)<-c(1:a[1]) # remove node names so extinction function can run
  rownames(fw_i)<-c(1:a[1])
  
  #among module degree in order
  roles<-Meta.top.roles[[2]]
  roles2<-roles
  roles[,3]<-cbind(1:a[1])
  sort<-sort(roles[2:3], decreasing = T)
  order_among.mod<-as.numeric(unlist(sort[,2]))
  #remove basal nodes from extinction order: autotrophs=2, mixotrophs=3, detritus=20
  remove<-c(2,3,20)
  order_among.mod<-setdiff(order_among.mod, remove)
  #what species go extinct?
  b<-sort(roles2[1:3], decreasing = T)
  remove2<-c("Autotroph", "Mixotroph", "Detritus")
  b<-as.character(unlist(b[,3]))
  b<-setdiff(b, remove2)
  
  #extinction högsta till lägsta among module degree
  among.mod.extinc<-ExtinctionOrder(fw_i, Order = order_among.mod[1:24], NetworkType = "Trophic")
  extinc.plots<-among.mod.extinc
  #ExtinctionPlot(History = extinc.plots) # did not work in loop with list, something about the format
  
  #R50
  dim(fw_i)[1]/2
  R50among.mod.mw<-length(which(among.mod.extinc$sims$TotalExt<= dim(fw_i)[1]/2))
  #R50among.mod[i]<-R50among.mod[i]/dim(fw_i)[1]
  #workaround ploting
  plot(extinc.plots$sims$NumExt, extinc.plots$sims$AccSecExt, type = "l",main = "Metaweb weighted, among.mod.deg",
       xlab="Number of primary extinctions", ylab="Accumulated Secondary Extinctions", ylim = c(0,34), xlim = c(0,34))
  abline(v=R50among.mod.mw)
  legend("topleft", inset = 0.1, legend= R50among.mod.mw, title = "R50=", box.lty = 0)
  legend("right", legend = b[1:R50among.mod.mw], title = "Extinction order:", box.lty = 0)
}




#weighted metaweb extinction by among module degree TH 30%
{
  fw_i<-as.matrix(LW_Metaweb_GOR, matrix.type = "adjacency", attr= "weight", sparse= F)
  sp<-colnames(metaweb)
  a<-dim(fw_i)
  colnames(fw_i)<-c(1:a[1]) # remove node names so extinction function can run
  rownames(fw_i)<-c(1:a[1])
  
  #among module degree in order
  roles<-Meta.top.roles[[2]]
  roles2<-roles
  roles[,3]<-cbind(1:a[1])
  sort<-sort(roles[2:3], decreasing = T)
  order_among.mod<-as.numeric(unlist(sort[,2]))
  #remove basal nodes from extinction order: autotrophs=2, mixotrophs=3, detritus=20
  remove<-c(2,3,20)
  order_among.mod<-setdiff(order_among.mod, remove)
  
  #what species go extinct?
  b<-sort(roles2[1:3], decreasing = T)
  remove2<-c("Autotroph", "Mixotroph", "Detritus")
  b<-as.character(unlist(b[,3]))
  b<-setdiff(b, remove2)
  #extinction högsta till lägsta among module degree
  #fw_i<-as.network(fw_i, matrix.type = "adjacency", attrnames= "weight",names.eval = "weight", ignore.eval=FALSE)
  among.mod.extinc<-ExtinctionOrder(fw_i, Order = order_among.mod[1:24], NetworkType = "Trophic", IS=0.3)
  extinc.plots<-among.mod.extinc
  #ExtinctionPlot(History = extinc.plots) # did not work in loop with list, something about the format
  
  #R50
  dim(fw_i)[1]/2
  R50among.mod.mw<-length(which(among.mod.extinc$sims$TotalExt<= dim(fw_i)[1]/2))
  #R50among.mod[i]<-R50among.mod[i]/dim(fw_i)[1]
  #workaround ploting
  plot(extinc.plots$sims$NumExt, extinc.plots$sims$AccSecExt, type = "l",main = "Metaweb weighted, among.mod.deg, TH 30%",
       xlab="Number of primary extinctions", ylab="Accumulated Secondary Extinctions", ylim = c(0,34), xlim = c(0,34))
  abline(v=R50among.mod.mw)
  legend("topleft", inset = 0.1, legend= R50among.mod.mw, title = "R50=", box.lty = 0)
  legend("right", legend = b[1:R50among.mod.mw], title = "Extinction order:", box.lty = 0)
}


#weighted metaweb extinction by within module degree
{
  fw_i<-(as.matrix(LW_Metaweb_GOR, matrix.type = "adjacency", attr= "weight", sparse= F))
  sp<-colnames(metaweb)
  a<-dim(fw_i)
  colnames(fw_i)<-c(1:a[1]) # remove node names so extinction function can run
  rownames(fw_i)<-c(1:a[1])
  
  #among module degree in order
  roles<-Meta.top.roles[[2]]
  roles2<-roles
  roles[,3]<-cbind(1:a[1])
  sort<-sort(roles[1:3], decreasing = T)
  order_within.mod<-as.numeric(unlist(sort[,3]))
  #remove basal nodes from extinction order: autotrophs=2, mixotrophs=3, detritus=20
  remove<-c(2,3,20)
  order_within.mod<-setdiff(order_within.mod, remove)
  #what species go extinct?
  b<-sort(roles2[1:3], decreasing = T)
  remove2<-c("Autotroph", "Mixotroph", "Detritus")
  b<-as.character(unlist(b[,3]))
  b<-setdiff(b, remove2)
  
  #extinction högsta till lägsta among module degree
  within.mod.extinc<-ExtinctionOrder(fw_i, Order = order_within.mod, NetworkType = "Trophic")
  extinc.plots<-within.mod.extinc
  #ExtinctionPlot(History = extinc.plots) # did not work in loop with list, something about the format
  
  #R50
  dim(fw_i)[1]/2
  R50within.mod.mw<-length(which(within.mod.extinc$sims$TotalExt<= dim(fw_i)[1]/2))
  #R50among.mod[i]<-R50among.mod[i]/dim(fw_i)[1]
  #workaround ploting
  plot(extinc.plots$sims$NumExt, extinc.plots$sims$AccSecExt, type = "l",main = "Metaweb weighted, within.mod.deg",
       xlab="Number of primary extinctions", ylab="Accumulated Secondary Extinctions", ylim = c(0,34), xlim = c(0,34))
  abline(v=R50within.mod.mw)
  legend("topleft", inset = 0.1, legend= R50within.mod.mw, title = "R50=", box.lty = 0)
  legend("right", legend = b[1:R50within.mod.mw], title = "Extinction order:", box.lty = 0)
}



#weighted metaweb extinction by within module degree TH 30%
{
  #fw_i<-as_adjacency_matrix(LW_Metaweb_GOR, attr = "weight", sparse = T, names = T)
  fw_i<-as.matrix(LW_Metaweb_GOR, matrix.type = "adjacency", attr= "weight", sparse= F)
  sp<-colnames(metaweb)
  a<-dim(fw_i)
  colnames(fw_i)<-c(1:a[1]) # remove node names so extinction function can run
  rownames(fw_i)<-c(1:a[1])
  
  #among module degree in order
  roles<-Meta.top.roles[[2]]
  roles2<-roles
  roles[,3]<-cbind(1:a[1])
  sort<-sort(roles[1:3], decreasing = T)
  order_within.mod<-as.numeric(unlist(sort[,3]))
  
  #remove basal nodes from extinction order: autotrophs=1, mixotrophs=2, detritus=20
  remove<-c(2,3,20)
  order_within.mod<-setdiff(order_within.mod, remove)
  
  #what species go extinct?
  b<-sort(roles2[1:3], decreasing = T)
  remove2<-c("Autotroph", "Mixotroph", "Detritus")
  b<-as.character(unlist(b[,3]))
  b<-setdiff(b, remove2)
  
  #extinction högsta till lägsta among module degree
  #x<-get.data.frame(LW_Metaweb_GOR)
  #x<-as.network(x, matrix.type= "edgelist", attrnamse="weight", ignore.eval=FALSE)
  #as.network(fw_i, matrix.type = "adjacency", attrname="weight",ignore.eval=FALSE, names.eval='weight')
  #as.matrix(LW_Metaweb_GOR, matrix.type = "adjacency", attr= "weight", sparse= F)
  
  within.mod.extinc<-ExtinctionOrder(fw_i, Order = order_within.mod[1:24], NetworkType = "Trophic", IS=0.3)
  extinc.plots<-within.mod.extinc
  #ExtinctionPlot(History = extinc.plots) # did not work in loop with list, something about the format
  
  #R50
  dim(fw_i)[1]/2
  R50within.mod.mw<-length(which(within.mod.extinc$sims$TotalExt<= dim(fw_i)[1]/2))
  #R50among.mod[i]<-R50among.mod[i]/dim(fw_i)[1]
  #workaround plotinga
  plot(extinc.plots$sims$NumExt, extinc.plots$sims$AccSecExt, type = "l",main = "Metaweb weighted, within.mod.deg TH 30%",
       xlab="Number of primary extinctions", ylab="Accumulated Secondary Extinctions", ylim = c(0,34), xlim = c(0,34))
  abline(v=R50within.mod.mw)
  legend("topleft", inset = 0.1, legend= R50within.mod.mw, title = "R50=", box.lty = 0)
  legend("right", legend = b[1:R50within.mod.mw], title = "Extinction order:", box.lty = 0)
}





# Weighted mw extinction by sum of link weights diffrent thresholds, same plot
{
  fw_i<-as.matrix(LW_Metaweb_GOR, matrix.type = "adjacency", attr= "weight", sparse= F)
  sp<-colnames(metaweb)
  a<-dim(fw_i)
  colnames(fw_i)<-c(1:a[1])# remove node names so extinction function can run
  rownames(fw_i)<-c(1:a[1])
  TH<-list() #threshold for extinction
  TH<-c(0,0.01,0.05,0.1,0.2,0.3,0.6,0.99) #threshold for extinction
  extinc.plots<-list()
  R50multiTH<-list()
  R50index<-list()
  #sum of of links
  fw<-graph_from_adjacency_matrix(
    fw_i,
    mode = c("directed"),
    weighted = T
  )
  sum.link.weights.sp<-apply(fw_i, 1, sum)+apply(fw_i, 2, sum)
  
  order_sum_link<-as.numeric(names(sort(sum.link.weights.sp, decreasing = T)))
  
  #what species go extinct
  #b<-sort(`names<-`(sum.link.weights.sp, sp), decreasing = T)
  #b<-names(b)
  #remove2<-c("Autotroph", "Mixotroph", "Detritus")
  #b<-setdiff(b, remove2)
  
  #remove basal nodes from extinction order: autotrophs=2, mixotrophs=3, detritus=20
  remove<-c(2,3,20)
  order_sum_link<-setdiff(order_sum_link, remove)
  
  #extinction sum links
  for(i in 1:length(TH)){
  sum.link.ext<-ExtinctionOrder(fw_i, Order=order_sum_link, NetworkType = "Trophic", IS=TH[[i]])
  extinc.plots[[i]]<-sum.link.ext
  R50index[i]<-sum.link.ext$R50
  #ExtinctionPlot(History = sum.link.ext$sims)
  #R50
  R50multiTH[i]<-length(which(sum.link.ext$sims$TotalExt<= dim(fw_i)[1]/2))
  }
  #what species go extinct?
  #sp[order_sum_link[1:R50.WMW]]
  #ploting
  plot(extinc.plots[[1]]$sims$NumExt, extinc.plots[[1]]$sims$AccSecExt, type = "l",main = "Weighted metaweb, extinction by sum.link.weigth",
       xlab="Number of primary extinctions", ylab="Secondary Extinctions", ylim = c(0,34), xlim = c(0,34))
  #abline(v= R50multiTH)
  #legend("topleft", inset = 0.1, legend= R50.WMW, title = "R50=", box.lty = 0)
  #legend("right", legend =  b[1:R50.WMW], title = "Extinction order:", box.lty = 0)
  plotcol <- c("red","blue","gray","black", "cyan","orange","purple","darkgreen")
  for (j in 1:length(R50multiTH)){
    
  lines(extinc.plots[[j]]$sims$NumExt, extinc.plots[[j]]$sims$AccSecExt, col=plotcol[j], lwd=2)
  #points(R50multiTH[j], 16, col=j, pch=15)
   # abline(v= R50multiTH[j], col=j)
  }
  legend("bottomright", legend= R50multiTH, title = "R50=", col = plotcol, pch = 15)
  #legend("topright", legend = TH*100, title = "Threshold %", col = c(9,8,6,5,4), pch = 15)
  legend("topright", legend = TH*100, title = "% left of incoming flow
         required to trigger extinction", col = plotcol, pch = 15, box.lty = 0)
}
R50index

RandomExtinctions(fw_i)
data(net)

weight.derp<-as.network(fw_i, matrix.type = "adjacency", 
                        ignore.eval = FALSE, names.eval = 'weight')
summary(weight.derp)

derp<-RandomExtinctions(Network= weight.derp, nsim= 100, SimNum =24, NetworkType = "Trophic")
CompareExtinctions(Nullmodel=derp, Hypothesis = sum.link.ext)


#degree vs null model
#derp<-RandomExtinctions(Network= fw_i, nsim= 100, SimNum =, NetworkType = "Trophic")
CompareExtinctions(Nullmodel=derp, Hypothesis = degree.extinc.mw)
CompareExtinctions(Nullmodel=derp, Hypothesis = among.mod.extinc)
CompareExtinctions(Nullmodel=derp, Hypothesis = within.mod.extinc)

set.seed(42)
#null model test
w_null<-RandomExtinctions(Network= fw_i, nsim= 1000, SimNum =31, NetworkType = "Trophic", IS=0.3)
CompareExtinctions(Nullmodel=w_null, Hypothesis = extinc.plots[[8]])

glist<-list()
for(i in 1:length(extinc.plots)){
glist[[i]]<-CompareExtinctions(Nullmodel=w_null, Hypothesis = extinc.plots[[i]])
}
glist

uw_null<-RandomExtinctions(Network= fw_i, nsim= 1000, SimNum =31, NetworkType = "Trophic")
CompareExtinctions(Nullmodel=uw_null, Hypothesis = degree.extinc.mw)

#rewiring test, need to set up so = producers don't worry about rewiring, # intermediate consumders can only rewire to producers
# apex predators can only rewire to intermediate consumers
#dist<-fw_i
#dist[,c(2,3,20)]
#dsit[,c()]
#deg.sp<-degree(LW_Metaweb_GOR, mode = c("in"),loops = TRUE,normalized = FALSE)
#names(sort(deg.sp, decreasing = T))

#rewire.test<-ExtinctionOrder(fw_i, Order=order_sum_link[1:24], NetworkType = "Trophic", IS=0.3,#
#                             Rewiring = function(x){x}, RewiringDist = dist, RewiringProb = 0.5)

#slightly better plot for all thresholds
{  fw_i<-as.matrix(LW_Metaweb_GOR, matrix.type = "adjacency", attr= "weight", sparse= F)
  sp<-colnames(metaweb)
  a<-dim(fw_i)
  colnames(fw_i)<-c(1:a[1])# remove node names so extinction function can run
  rownames(fw_i)<-c(1:a[1])
  TH<-list() #threshold for extinction
  TH<-c(seq(0, 0.9, by = 0.1)) #threshold for extinction
  extinc.plots<-list()
  R50multiTH<-list()
  R50index<-list()
  #sum of of links
  fw<-graph_from_adjacency_matrix(
    fw_i,
    mode = c("directed"),
    weighted = T
  )
  sum.link.weights.sp<-apply(fw_i, 1, sum)+apply(fw_i, 2, sum)
  
  order_sum_link<-as.numeric(names(sort(sum.link.weights.sp, decreasing = T)))
  
  #remove basal nodes from extinction order: autotrophs=2, mixotrophs=3, detritus=20
  remove<-c(2,3,20)
  order_sum_link<-setdiff(order_sum_link, remove)
  
  #extinction sum links
  for(i in 1:length(TH)){
    sum.link.ext<-ExtinctionOrder(fw_i, Order=order_sum_link, NetworkType = "Trophic", IS=TH[[i]]) ##0.3 leads to removal of all nodes which lose 70 percent of their interaction strength in the Network argument)
    extinc.plots[[i]]<-sum.link.ext
    R50index[i]<-sum.link.ext$R50
    #ExtinctionPlot(History = sum.link.ext$sims)
    #R50
    R50multiTH[i]<-length(which(sum.link.ext$sims$TotalExt<= dim(fw_i)[1]/2))
  }
  #what species go extinct?
  #sp[order_sum_link[1:R50.WMW]]
  #ploting
  plot(extinc.plots[[1]]$sims$NumExt, extinc.plots[[1]]$sims$AccSecExt, type = "l",main = "Weighted metaweb, extinction by sum.link.weigth",
       xlab="Number of primary extinctions", ylab="Secondary Extinctions", ylim = c(0,34), xlim = c(0,34))
  #abline(v= R50multiTH)
  #legend("topleft", inset = 0.1, legend= R50.WMW, title = "R50=", box.lty = 0)
  #legend("right", legend =  b[1:R50.WMW], title = "Extinction order:", box.lty = 0)
  #plotcol <- c(topo.colors((11)))
  plotcol<-brewer.pal(name = "Spectral", n=length(extinc.plots))
  
  for (j in 1:length(R50multiTH)){
    
    lines(extinc.plots[[j]]$sims$NumExt, extinc.plots[[j]]$sims$AccSecExt, col=plotcol[j], lwd=2)
    points(extinc.plots[[j]]$sims$NumExt, extinc.plots[[j]]$sims$AccSecExt, col=plotcol[j], lwd=2, type = "p")
    #points(R50multiTH[j], 16, col=j, pch=15)
    # abline(v= R50multiTH[j], col=j)
  }
  #legend("topright", legend = TH*100, title = "Threshold %", col = c(9,8,6,5,4), pch = 15)
  #legend("topright", legend = TH*100, title = "% left of incoming flow
  #      required to trigger extinction", col = plotcol, pch = 15, box.lty = 0)
  legend("topleft", legend= R50index, title = "R50", col = plotcol, pch = 15)
  legend("topright", legend = c("Binary", "90%", "80%", "70%", "60%", "50%", "40%", "30%", "20%","10%"), 
  title = "Threshold,
  % loss of incoming energy
  needed to trigger extinction", col = plotcol, pch = 15, box.lty = 0, inset=0.05)
  legend("bottomright", legend= R50multiTH, title = "Nr prim.ext.R50", col = plotcol, pch = 15)
}

