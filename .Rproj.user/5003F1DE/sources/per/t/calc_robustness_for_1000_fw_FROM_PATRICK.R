
setwd("C:/LocalData/susakort/GitHub/GoR-Food-webs/Data")
setwd("D:/R_Git/GoR-Food-webs/Data")# Patriks setwd 

setwd(".../Data")
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
  R50[j,i]<-degree.extinc.mw$R50 ## R50 value extracted from the extinctionoutput, R50 between 1 and 0
  n.sp.extR50[j,i]<-length(which(degree.extinc.mw$sims$TotalExt<= dim(fw_i)[1]/2)) # number of species extinct at R50
  

  #YOU NEED TO DO THE SAME FOR WEIGHTED NETWORKS
  #YOU CAN USE THE BIOMASS INFO FOR THIS AND THE UNEWEIGHTED fw
  
 #se work in progress below
  
  #INSTEAD OF MATRICES YOU CAN ASLO STORE YOUR OUTPUT IN LISTS
  #THESE LISTS WE CAN ADD TO ARRAYS AND THEN CALCULATE THE MEAN VALUES
  
  }
}

index.median<- apply(index, 1, median)
index.median
round(apply(deg, 1, median), 2)
apply(R50, 1, median)# R50 value extracted from the extinctionoutput, same as extinctionoutput[[i]]$R50
apply(n.sp.extR50,1,median)

extinctionoutput #output is a list, so we have a list of lists again
extinctionoutput[[1]]$sims #this is what I used for the previous plots of extinctions



##### Un-weighted just networks and metrics, no extinctions, will run another loop for that####
index<-matrix(NA, 34, 1000)
index.arr <- list() #list of output for the 1000 food webs per year

for (i in 1:1000){
  
  #i<-1
  sp.pa.year<-1*((t(biomC[,,i])[, match(row.names(un_mw), row.names(biomC[,,1]))])>0) #reorder to use for making networks (1000 per year)
  
  webs<-list()
  for (j in 1:34){
    fw.sel<-which(sp.pa.year[j,]==1)  
    fw<-un_mw[fw.sel, fw.sel]
    
    index[j,i]<-dim(fw)[1] #number of species
    #deg[j,i]<-c(mean(degree(graph.adjacency(fw)))) # mean degree test, will add more later when weights are up and running
    #links
    #generality
    #vulnerability
    #connectance
    
    #INSTEAD OF MATRICES YOU CAN ASLO STORE YOUR OUTPUT IN LISTS
    #THESE LISTS WE CAN ADD TO ARRAYS AND THEN CALCULATE THE MEAN VALUES
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
w.webs.1000[[1]][[1]]# so this would be the first of a 1000 webs for the first year, right?

w.webs.1000[[10]][[34]]
w.webs.1000[[1000]][[34]]

saveRDS(w.webs.1000, file="w.webs.1000.RDS")
w.webs.1000<-readRDS("w.webs.1000.RDS")







#### extinction loop on list of lists or(index.arr= 1000 networks per year) test####

#R50 test run with un-weighted extinction by degree, when fluxes are integrated can change it to/ add extinction by sum of link weights
R50<-matrix(NA, 34, 1000)
#n.sp.extR50<-matrix(NA, 34, 1000)
ext.ind<-list()

for (i in 1:1000) {
  #i<-1
  a<-index.arr[[i]]
  
  ext.output<-list()
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
    deg.sp<-degree(fw1, mode = c("all"),loops = TRUE,normalized = FALSE)
    order_degree<-as.numeric(names(sort(deg.sp, decreasing = T))) #order of nodes by degree
    
    #remove basal nodes from extinction order: autotrophs, mixotrophs, detritus
    order_degree<-setdiff(order_degree, c(which(colnames(fw)=="Autotroph"),
                                          which(colnames(fw)=="Mixotroph"),
                                          which(colnames(fw)=="Detritus")))
    
    #ext. by highest to lowest degree
    degree.extinc.mw<-ExtinctionOrder(fw_i, Order = order_degree, NetworkType = "Trophic",verbose = F)
    ext.output[[j]]<-degree.extinc.mw
    #R50
    R50[j,i]<-degree.extinc.mw$R50 ## R50 value extracted from the ext.output, R50 between 1 and 0
    #n.sp.extR50[i,j]<-length(which(degree.extinc.mw$sims$TotalExt<= dim(fw_i)[1]/2)) # number of species extinct at R50
    
  }
  ext.ind[[i]]<-ext.output
}

ext.ind[[1]][[1]]
ext.ind[[1000]][[34]]
ext.ind[[1]][[1]]$R50
ext.ind[[1000]][[34]]$R50
saveRDS(ext.ind, file="uw.ext.output.1000.rds")
saveRDS(R50,file="uw.R50.output.1000.rds")
uw.ext.output.1000<-readRDS("uw.ext.output.1000.rds")
uw.R50.output.1000<-readRDS("uw.R50.output.1000.rds")
uw.R50.output.1000
apply(R50,1,median)







#### weighted extinctions using w.webs.1000, 1000 webs for each year, sum of link weights to determine extinction order####
w.R50<-matrix(NA, 34, 1000)
#n.sp.extR50<-matrix(NA, 34, 10)
w.ext.ind<-list()

for (i in 1:1000) {
  #i<-1
  a<-w.webs.1000[[i]]
  
  ext.output<-list()
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
    ext.output[[j]]<-sumlink.extinc.mw
    #R50
    w.R50[j,i]<-sumlink.extinc.mw$R50 ## R50 value extracted from the ext.output, R50 between 1 and 0
    #n.sp.extR50[i,j]<-length(which(degree.extinc.mw$sims$TotalExt<= dim(fw_i)[1]/2)) # number of species extinct at R50
    
  }
  w.ext.ind[[i]]<-ext.output
}

w.ext.ind[[1]][[1]]$sims # can be used for the extinction plot that i had earlier
w.ext.ind[[1000]][[34]]
w.R50

saveRDS(w.ext.ind, file = "w.ext.ind.1000.RDS")
saveRDS(w.R50,file="w.R50.output.1000.RDS")

apply(w.R50, 1, median)
apply(R50, 1, median)

hist(R50)
hist(w.R50)

readRDS("w.ext.ind.1000.RDS")
w.R50th30<-readRDS("w.R50.output.1000.RDS")# used loss of 70% as threshold for this one
R50=uw.R50.output.1000
uw.R50<-apply(R50, 1, median)
wei.R50<-apply(w.R50, 1, median)

####ploting R50####
df<-matrix(NA,34,7)
colnames(df)<-c("year","uw_median","w_median","uw_min", "uw_max", "w_min", "w_max")
df[,1]<-c(seq(1981,2014))
df[,2]<-uw.R50
df[,3]<-wei.R50
for (i in 1:34) {
  df[i,4]<-min(R50[i,])
  df[i,5]<-max(R50[i,])
}

for (i in 1:34) {
  df[i,6]<-min(w.R50[i,])
  df[i,7]<-max(w.R50[i,])
}



{plot(df[,1],df[,2], type = "l", col="red", lwd=2, ylim = c(0.1,0.47), xlim = c(1981,2014),
     main = "R50 degree vs R50 sum of link weight, median of 100",
     xlab="Year", ylab="R50 index", xaxt='n')
axis(side = 1, at=seq(1981, 2014, by = 1), las=2)
lines(df[,1],df[,3], col="blue", lwd=2)
}


df<-as.data.frame(df)
g<-ggplot(df, aes(x=year))+
  geom_line(aes(y=uw_median), col="red",linewidth=1.1)+
  geom_ribbon(aes(ymin=uw_min , ymax=uw_max), fill="red", alpha=0.2)+ 
  geom_ribbon(aes(ymin=w_min, ymax=df[,7]), fill="blue",alpha=0.2)+
  geom_line(aes(x=year, y=w_median ), col="blue", linewidth=1.1)+
  xlim(c(1981,2014))+ylim(c(0.05, 0.55))+
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  scale_x_continuous(name = "Year", breaks = 1981:2014)+
  xlab("Year")+ylab("R50 index")


g+ggtitle("R50(binary) vs R50 (loss of 70%= extinction) over time")+theme(plot.title = element_text(hjust = 0.5))


g<-ggplot(df)+
  geom_line(aes(x=year, y=uw_median, colour=factor(uw_median)),size=1.1)+
  geom_ribbon(aes(x=year, ymin=uw_min , ymax=uw_max), fill="red", alpha=0.1)+ 
  geom_ribbon(aes(x=year,ymin=w_min, ymax=df[,7]), fill="blue", alpha=0.1)+
  geom_line(aes(x=year, y=w_median, colour =factor(wR50)),color="blue", size=1.1)+
  xlim(c(1981,2014))+ylim(c(0.05, 0.55))+
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  scale_x_continuous(name = "Year", breaks = 1981:2014)+
  xlab("Year")+ylab("R50 index")

g+labs(colour= "method")

pivot_longer(df, year)
# need to make long format of df for easier plot with legend?




#####code spinets####

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




