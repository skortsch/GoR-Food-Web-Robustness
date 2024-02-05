####GOR robustness nullmodel test####
library(igraph)
library(fluxweb)
library(MASS)
library(NetIndices)
library(fluxweb)
library(NetworkExtinction)
library(network)

source("curveball.R")
source("curv.R")

#time by species matrix
biomMedian<-apply(biomC, c(1,2),median)
time.sp.matrix<-(biomMedian>0)*1
time.sp.matrix<-t(time.sp.matrix)#want to have columns as sp. and year as row

biomMin<-apply(biomC, c(1,2),min)
time.sp.min<-(biomMin>0)*1
time.sp.min<-t(time.sp.min)

biomMax<-apply(biomC, c(1,2),max)
time.sp.max<-(biomMax>0)*1
time.sp.max<-t(time.sp.max)

#load the time by species matrix 
poly<-time.sp.min#time.sp.max, time.sp.min, time.sp.matrix (median)
#poly<-time.sp.min
#poly<-time.sp.max
#load the food web metaweb 
#fw<-read.table("fw.metaweb.baltic.txt")
fw<-get.adjacency(net)
#fw<-un_mw
# need to specify names of basal species, this example (9 species) is from the Barents Sea data set
basal=colnames(poly[TRUE, c("Autotroph","Mixotroph","Detritus")]) # names of basal species, not to be randomized 

model_curveball<-list()
iter<-100 # no of iterations
for (l in 1:iter){ # loop on randomisations
  
  null.model.curveball<-matrix(0, dim(poly)[1], dim(poly)[2]) # empty matrix with rows (number of rows are equal to the number of food webs) to be filled in with swapped elements (i.e. prey items)
  
  cat('iteration ',l,'\n')
  
  poly.ran<-curve_ball(poly) # use the swap alogorithm to create a new dataset (XX rows * XX potential species)
  
  i=1
  while (i <= dim(poly)[1]){ 					# loop on the rows to check that they are valid, i.e. fullfil the constraints
    
    id.s<-which(poly.ran[i,]==1)  # species present in polygon i
    fw.ids<-fw[id.s,id.s]         # reduced foodweb, with only species present
    
    # test for prey
    n.prey=apply(fw.ids, 2, sum) # number of prey for ALL species in the reduced foodweb
    zero.prey.list=names(n.prey[n.prey==0]) # list of species with no prey
    test.4.zero.prey=(sum(zero.prey.list%in%basal==FALSE)<=0) # in Ecography paper: test no more than 1 non-basal species could have no prey
    
    # test for network unicity
    ids.g<-graph.adjacency(fw.ids)	#identify unconnected graphs
    cl<-clusters(ids.g)					    #one cluster means the graph is connected, if more clusters some species are singletons 
    
    if (test.4.zero.prey==TRUE & cl$no==1){ # test that the prey & connected criteria are fullfilled
      cat('polygon',i,'passed \n') # if yes, move to the next polygon
      i=i+1
    }else{ # if no, reshuffle the species matrix and restart at polygon 1
      cat('polygon',i,'failed, zero.prey species = ',sum(zero.prey.list%in%basal==FALSE),'nsg=',cl$no,'\t')
      cat(zero.prey.list[zero.prey.list%in%basal==FALSE],'\n')
      poly.ran<-curve_ball(poly) # use the swap alogorithm to create a new dataset (25 polygons * 233 potential species) 
      i=1 
    }
  }
  model_curveball[[l]]<- as.matrix(poly.ran) # the list of iterations
}

model_curveball_median<-model_curveball
saveRDS(model_curveball_median,file="curv.mod.RDS")
model_curveball_median<-readRDS("curv.mod.RDS")


model_curveball_min<-model_curveball
saveRDS(model_curveball_min, file="curv.mod.min.RDS")

model_curveball_max<-model_curveball
saveRDS(model_curveball_max, file="curv.mod.max.RDS")

model_curveball_min<-readRDS("curv.mod.min.RDS")
model_curveball_max<-readRDS("curv.mod.max.RDS")

#code from MS, sp for 1 year 999 iterations
#fw.1981<-list()
#for (i in 1:iter){
#  id.s<-which(model_curveball[[i]][1,]==1)
#  fw.1981[[i]]<-as.matrix(fw[id.s,id.s])
#}



#multiple years loop
nullmodelwebs<-list()
webs<-list()
for (i in 1:iter){
  a<-model_curveball_max[[i]]#model_curveball_max, model_curveball_min, model_curveball_median
  for (j in 1:34){
  id.s<-which(a[j,]==1)  # species present in YEAR
  webs[[j]]<-as.matrix(fw[id.s,id.s])
  }
  nullmodelwebs[[i]]<-webs # list of list, can use with R50 code loop
}

nullmodelwebs[[1]][[1]]# iterartion xx of year xx
nullmodelwebs[[100]][[34]]


null.min<-nullmodelwebs
null.max<-nullmodelwebs

saveRDS(nullmodelwebs, file = "median.nullnets.RDS")
saveRDS(null.min, file = "min.nullnets.RDS")
saveRDS(null.max, file = "max.nullnets.RDS")
####weighted null####

#load data
#biom.data<-read.table("clipboard", dec=",", header=T)
#biom.dat<-as.data.frame(t(biom.data))
#write.csv(biom.dat, "biom.dat.csv", row.names=F)



#biomass data yearly
biom.data<-apply(biomC, c(1,2),median)
biom.data<-t(biom.data)

biom.data.min<-apply(biomC, c(1,2),min)
biom.data.min<-t(biom.data.min)

biom.data.max<-apply(biomC, c(1,2),max)
biom.data.max<-t(biom.data.max)


#binary randomisations time by species matrix
mods_curve<-readRDS("curv.mod.RDS")
mods_curve<-readRDS("curv.mod.min.RDS")
mods_curve<-readRDS("curv.mod.max.RDS")


weighted_ran_mat<-list()
for (i in 1:length(mods_curve)){
  ran.mat<-mods_curve[[i]]
  for (j in 1:dim(ran.mat)[2]){
    
    id<-which(ran.mat[,j]==1) # row id of species j -
    
    if (length(id)>0){
      non.zero.vec<-biom.data.min[,j][biom.data.min[,j]!=0] #remove the zeros #biom.data.min, biom.data.max, biom.data 
      
      ran.mat[id,j]<-sample(non.zero.vec) #randomise biomass info and put back into the row id's
    }
  }
  weighted_ran_mat[[i]]<-ran.mat
}
saveRDS(weighted_ran_mat, file="weighted_ran_mat.RDS")
saveRDS(weighted_ran_mat, file="weighted_ran_mat.min.RDS")
saveRDS(weighted_ran_mat, file="weighted_ran_mat.max.RDS")
readRDS("weighted_ran_mat.RDS")
weighted_ran_mat[[1]]
weighted_ran_mat[[100]]
#ran.info.1981<-list()
#for (i in 1:length(weighted_ran_mat)){
#  # i<-1
 # ran.info[,4]<-weighted_ran_mat[[i]][1,]
 # ran.info.1981[[i]]<-ran.info
#}
ran.info<-info
rownames(ran.info)
ran.info.all<-list()
wbs<-list()
for (i in 1:length(weighted_ran_mat)){
  b<-weighted_ran_mat[[i]]
  for (j in 1:34) {
    ran.info[,4]<-b[j,]
    wbs[[j]]<-ran.info
  }
  ran.info.all[[i]]<-wbs
}
ran.info.all[[1]][[1]]$meanB


#weight.ran.flux1981<-list()
#for (i in 1:iter){
#  ids<-which(ran.info.1981[[i]]$meanB=="0")
#  fluxes <- fluxing(fw.1981[[i]], ran.info.1981[[i]][-ids,]$meanB, ran.info.1981[[1]][-ids,]$losses, 
#                    ran.info.1981[[i]][-ids,]$efficiencies, ef.level="prey")
#  fluxes <- fluxes*86.4 # conversion from J/sec to kJ/day 
#  # 1 J/sec = 86.4 kJ/day (there are 86400 sec/day)
#  weight.ran.flux1981[[i]]<-fluxes
#}


weight.ran.flux<-list()
y<-list()
for (i in 1:iter){
  mat<-nullmodelwebs[[i]]
  inf<-ran.info.all[[i]]
  for (j in 1:34) {
  ids<-which(inf[[j]]$meanB=="0")
  fluxes <- fluxing(mat[[j]], inf[[j]][-ids,]$meanB, inf[[j]][-ids,]$losses, 
                    inf[[j]][-ids,]$efficiencies, ef.level="prey")
  fluxes <- fluxes*86.4 # conversion from J/sec to kJ/day 
  # 1 J/sec = 86.4 kJ/day (there are 86400 sec/day)
  y[[j]]<-fluxes
  }
  weight.ran.flux[[i]]<-y
}
weight.ran.flux[[1]][[1]]
weight.ran.flux[[100]][[1]]



####Null model R50####

w.null.R50th.bin<-matrix(NA, 34, 100)
w.null.R50th.10<-matrix(NA, 34, 100)
w.null.R50th.20<-matrix(NA, 34, 100)
w.null.R50th.30<-matrix(NA, 34, 100)
w.null.R50th.40<-matrix(NA, 34, 100)
w.null.R50th.50<-matrix(NA, 34, 100)
w.null.R50th.60<-matrix(NA, 34, 100)
w.null.R50th.70<-matrix(NA, 34, 100)
w.null.R50th.80<-matrix(NA, 34, 100)
w.null.R50th.90<-matrix(NA, 34, 100)
#n.sp.extR50<-matrix(NA, 34, 10)
#w.ext.ind.th40<-list()
for (i in 1:100) {
  #i<-1
  a<-weight.ran.flux[[i]]
  
  #ext.output<-list()
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
    bin<-ExtinctionOrder(fw_i, Order = order_sum_link, NetworkType = "Trophic",verbose = F)
    th10<-ExtinctionOrder(fw_i, Order = order_sum_link, NetworkType = "Trophic",verbose = F, IS=0.1)
    th20<-ExtinctionOrder(fw_i, Order = order_sum_link, NetworkType = "Trophic",verbose = F, IS=0.2)
    th30<-ExtinctionOrder(fw_i, Order = order_sum_link, NetworkType = "Trophic",verbose = F, IS=0.3)
    th40<-ExtinctionOrder(fw_i, Order = order_sum_link, NetworkType = "Trophic",verbose = F, IS=0.4)
    th50<-ExtinctionOrder(fw_i, Order = order_sum_link, NetworkType = "Trophic",verbose = F, IS=0.5)
    th60<-ExtinctionOrder(fw_i, Order = order_sum_link, NetworkType = "Trophic",verbose = F, IS=0.6)
    th70<-ExtinctionOrder(fw_i, Order = order_sum_link, NetworkType = "Trophic",verbose = F, IS=0.7)
    th80<-ExtinctionOrder(fw_i, Order = order_sum_link, NetworkType = "Trophic",verbose = F, IS=0.8)
    th90<-ExtinctionOrder(fw_i, Order = order_sum_link, NetworkType = "Trophic",verbose = F, IS=0.9)
    #R50
    w.null.R50th.bin[j,i]<-bin$R50
    w.null.R50th.10[j,i]<-th10$R50
    w.null.R50th.20[j,i]<-th20$R50
    w.null.R50th.30[j,i]<-th30$R50
    w.null.R50th.40[j,i]<-th40$R50
    w.null.R50th.50[j,i]<-th50$R50
    w.null.R50th.60[j,i]<-th60$R50
    w.null.R50th.70[j,i]<-th70$R50
    w.null.R50th.80[j,i]<-th80$R50
    w.null.R50th.90[j,i]<-th90$R50
    # R50 value extracted from the ext.output, R50 between 1 and 0
    #n.sp.extR50[i,j]<-length(which(degree.extinc.mw$sims$TotalExt<= dim(fw_i)[1]/2)) # number of species extinct at R50
  }
  #w.ext.ind.th40[[i]]<-ext.output
}

max.null.th<-list(w.null.R50th.bin,w.null.R50th.10,w.null.R50th.20,w.null.R50th.30,
                     w.null.R50th.40,w.null.R50th.50,w.null.R50th.60,w.null.R50th.70,
                     w.null.R50th.80,w.null.R50th.90)

saveRDS(median.null.th, file = "median.null.th.RDS")
saveRDS(min.null.th, file = "min.null.th.RDS")
saveRDS(max.null.th, file = "max.null.th.RDS")



####uw null metrics####
readRDS("median.nullnets.RDS")
readRDS("min.nullnets.RDS")
readRDS("max.nullnets.RDS")

uw.null.con<-matrix(NA,34,100)
uw.null.sp<-matrix(NA,34,100)
for (i in 1:100){
  
  #i<-1
  uwnullnet<-nullnet[[i]]
    
  for (j in 1:34){
    fw.sel<-which(uwnullnet[j,]==1)  
    fw<-un_mw[fw.sel, fw.sel]
    fw_g<-graph.adjacency(fw, mode = "directed")
    
    uw.null.sp[j,i]<-dim(fw)[1] #number of species
    #connectance 
    uw.null.con[j,i]<-edge_density(fw_g)
    
  }
}
median.null.con<-uw.null.con
min.null.con<-uw.null.con
max.null.con<-uw.null.con

####weigthed con####
w.null.metrics<-list()
for (i in 1:100){
  #i<-1
  a<-weight.ran.flux[[i]]
  met.list<-list()
  
  for (j in 1:34){
    fw<-a[[j]]
    
    met.list[[j]]<-fluxind(fw)
    
  }
  w.null.metrics[[i]]<-met.list
}

w.null.metrics[[100]][[34]]

w.null.con<-matrix(NA,34, 100)
for(i in 1:100){
  #i<-1
  a<-w.null.metrics[[i]]
  
  for( j in 1:34){
    #j<-1
    b<-a[[j]]$qC.w
    w.null.con[j,i]<-b
    
  }
}

median.null.wcon<-w.null.con
min.null.wcon<-w.null.con
max.null.wcon<-w.null.con
saveRDS(median.null.wcon, file = "median.null.wcon.RDS")
saveRDS(min.null.wcon, file = "min.null.wcon.RDS")
saveRDS(max.null.wcon, file = "max.null.wcon.RDS")

apply(w.null.con,1,median)





df.allTH<-matrix(NA,34,11)
df.allTH[,11]<-c(seq(1981,2014))
colnames(df.allTH)<-c("Binary","TH 10%","TH 20%", "TH 30%", "TH 40%", "TH 50%",
                      "TH 60%", "TH 70%", "TH 80%", "TH 90%", "Year")

for (i in 1:10){
  df.allTH[,i]<-apply(median.null.th[[i]],1, median)#median.null.th, min.null.th, max.null.th
}
as.data.frame(df.allTH)
d<-pivot_longer(as.data.frame(df.allTH), cols = 1:10)

ggplot(d, aes(x=Year, y=value, color = variable))+
  geom_line(aes(col=name),linewidth=1.1)+
  xlim(c(1981,2014))+ylim(c(0, 0.5))+
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, size = 25))+
  theme(axis.text.y = element_text(vjust = 0.5, size = 25))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  scale_x_continuous(name = "Year", breaks = scales::pretty_breaks(n = 5))+
  xlab("Year")+ylab("R50")+theme(axis.title=element_text(size=28))+
  theme(axis.title.y = element_text(angle = 0, vjust = 0.5))+
  scale_color_brewer(palette="Spectral")+
  ggtitle("Weighted nullmodel median food webs R50 with Thresholds")+theme(plot.title = element_text(hjust = 0.5, size = 25))+
  labs(color = "Thresholds")+
  theme(legend.text=element_text(size=25),
        legend.title=element_text(size=25))




null.wR50THplots<-list()
for (i in 1:10){
  data<-data.frame(as.vector(median.null.th[[i]]), as.vector(median.null.wcon))
  colnames(data)<-c("R50", "w.Con")
  null.wR50THplots[i]<-list(ggplot(data, aes(w.Con, R50))+
                              geom_point(alpha=1)+
                              stat_smooth(method = "lm",
                                          formula = "y~x",
                                          geom = "smooth",
                                          col="Red")+
                              scale_x_continuous(labels = scales::number_format(accuracy = 0.01), 
                                                 breaks = scales::pretty_breaks(n = 5))+
                              theme(axis.text.x = element_text(angle = 45, vjust = 0.5, size = 15))+
                              theme(axis.text.y = element_text(vjust = 0.5, size = 15))+
                              theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                    panel.background = element_blank(), axis.line = element_line(colour = "black"))+
                              xlab("w.con")+ylab("R50")+
                              theme(axis.title=element_text(size=15)))
  
}
null.wR50THplots[[1]]+ggtitle("null median TH Binary")+null.wR50THplots[[2]]+ggtitle(" null median TH 10%")+null.wR50THplots[[3]]+
  ggtitle("null median TH 20%")+null.wR50THplots[[4]]+ggtitle("null median TH 30%")+null.wR50THplots[[5]]+ggtitle("null median TH 40%")+
  null.wR50THplots[[6]]+ggtitle("null median TH 50%")+null.wR50THplots[[7]]+ggtitle("null median TH 60%")+
  null.wR50THplots[[8]]+ggtitle("null median TH 70%")+null.wR50THplots[[9]]+ggtitle("null median TH 80%")+
  null.wR50THplots[[10]]+ggtitle("null median TH 90%")



null.wR50THplots<-list()
for (i in 1:10){
  data<-data.frame(as.vector(min.null.th[[i]]), as.vector(min.null.wcon))
  colnames(data)<-c("R50", "w.Con")
  null.wR50THplots[i]<-list(ggplot(data, aes(w.Con, R50))+
                          geom_point(alpha=1)+
                          stat_smooth(method = "lm",
                                      formula = "y~x",
                                      geom = "smooth",
                                      col="Red")+
                          scale_x_continuous(labels = scales::number_format(accuracy = 0.01), 
                                             breaks = scales::pretty_breaks(n = 5))+
                          theme(axis.text.x = element_text(angle = 45, vjust = 0.5, size = 15))+
                          theme(axis.text.y = element_text(vjust = 0.5, size = 15))+
                          theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                panel.background = element_blank(), axis.line = element_line(colour = "black"))+
                          xlab("w.con")+ylab("R50")+
                          theme(axis.title=element_text(size=15)))
  
}
null.wR50THplots[[1]]+ggtitle("null min TH Binary")+null.wR50THplots[[2]]+ggtitle("null min TH 10%")+null.wR50THplots[[3]]+
  ggtitle("TH 20%")+null.wR50THplots[[4]]+ggtitle("null min TH 30%")+null.wR50THplots[[5]]+ggtitle("null min TH 40%")+
  null.wR50THplots[[6]]+ggtitle("null min TH 50%")+null.wR50THplots[[7]]+ggtitle("null min TH 60%")+
  null.wR50THplots[[8]]+ggtitle("null min TH 70%")+null.wR50THplots[[9]]+ggtitle("null min TH 80%")+
  null.wR50THplots[[10]]+ggtitle("null min TH 90%")


null.wR50THplots<-list()#min.null.th,min.null.wcon
for (i in 1:10){
  data<-data.frame(as.vector(max.null.th[[i]]), as.vector(max.null.wcon))
  colnames(data)<-c("R50", "w.Con")
  null.wR50THplots[i]<-list(ggplot(data, aes(w.Con, R50))+
                              geom_point(alpha=1)+
                              stat_smooth(method = "lm",
                                          formula = "y~x",
                                          geom = "smooth",
                                          col="Red")+
                              scale_x_continuous(labels = scales::number_format(accuracy = 0.01), 
                                                 breaks = scales::pretty_breaks(n = 5))+
                              theme(axis.text.x = element_text(angle = 45, vjust = 0.5, size = 15))+
                              theme(axis.text.y = element_text(vjust = 0.5, size = 15))+
                              theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                    panel.background = element_blank(), axis.line = element_line(colour = "black"))+
                              xlab("w.con")+ylab("R50")+
                              theme(axis.title=element_text(size=15)))
  
}
null.wR50THplots[[1]]+ggtitle("null max TH Binary")+null.wR50THplots[[2]]+ggtitle("null max TH 10%")+null.wR50THplots[[3]]+
  ggtitle("null max TH 20%")+null.wR50THplots[[4]]+ggtitle("null max TH 30%")+null.wR50THplots[[5]]+ggtitle("null max TH 40%")+
  null.wR50THplots[[6]]+ggtitle("null max TH 50%")+null.wR50THplots[[7]]+ggtitle("null max TH 60%")+
  null.wR50THplots[[8]]+ggtitle("null max TH 70%")+null.wR50THplots[[9]]+ggtitle("null max TH 80%")+
  null.wR50THplots[[10]]+ggtitle("null max TH 90%")





####null ecd plot####
metric.cdf<-ecdf(rand.mod.1981) #mod.vec is your vector with randomized modularity values
emp.cdf<-metric.cdf(obs.mod.1981) # empirical/observed modularity value

plot(metric.cdf)
abline(v=obs.mod.1981)

# want to know probability of getting same as or higher that empirical value 
#trying to test if the observed value is significantly higher than null hypothesis
# significant would usualy be >0.05

1-ecdf(rand.mod.1981)(obs.mod.1981)#p value = 0.383, would mean there is a ~38% chance of getting same or higher mod

#wrong? this would just give me % of mod over empirical value, right???
ecdf(rand.mod.1981)(obs.mod.1981)

# not sure if it works with multimodal dist???
rm2011.cdf<-ecdf(rand.mod.2011) #mod.vec is your vector with randomized modularity values
emp2011.cdf<-rm2011.cdf(obs.mod.2011) # empirical/observed modularity value
emp2011.cdf
plot(rm2011.cdf)
abline(v=obs.mod.2011)


plot(ecdf(rand.mod.1981))
points(x=obs.mod.1981, y=(ecdf(rand.mod.1981)(obs.mod.1981)), col="red", pch=16)




max.null.wcon
min.null.wcon
median.null.wcon

lol<-list()
for (i in 1:34){
median.null.wcon[i,]#vector with randomized w.con values for year i
c<-apply(w.con,1,median)# median value for emp webs, all years
c[i] #empirical/observed median R50 value for 1981
lol[[i]]<-ecdf(median.null.wcon[i,])(c[i])#
}
unlist(lol)


lol<-list()
for (i in 1:34){
  min.null.wcon[i,]#vector with randomized w.con values for year i
  c<-apply(w.con,1,min)# median value for emp webs, all years
  c[i] #empirical/observed median value for year i
  lol[[i]]<-ecdf(min.null.wcon[i,])(c[i])#
}
unlist(lol)

lol<-list()
for (i in 1:34){
  max.null.wcon[i,]#vector with randomized w.con values for year i
  c<-apply(w.con,1,max)# median value for emp webs, all years
  c[i] #empirical/observed median value for year i
  lol[[i]]<-ecdf(max.null.wcon[i,])(c[i])#
}
unlist(lol)

years<-seq(1981,2014)
for (i in 1:34) {
  c<-apply(w.con,1,median)# median value for emp webs, all years
  plot(ecdf(median.null.wcon[i,]), main=years[i], xlab="null median W.con", xlim=c(0.055,0.18))
  abline(v=c[i])
  #points(x=c[i], col="red", pch=16, cex= 1.5)
}

for (i in 1:34) {
  c<-apply(w.con,1,min)# median value for emp webs, all years
  plot(ecdf(min.null.wcon[i,]), main=years[i], xlab="null min W.con")
  abline(v=c[i])
  #points(x=c[i], col="red", pch=16, cex= 1.5)
}

for (i in 1:34){
  c<-apply(w.con,1,max)# median value for emp webs, all years, c[i] for specific years
  plot(ecdf(max.null.wcon[1,]), main=years[i], xlab="null max W.con")
  abline(v=c[i])
  points(x=c[i], y=(ecdf(max.null.wcon[i])(c[i])), col="red", pch=16, cex= 1.5)
}


metric.cdf<-ecdf(median.null.wcon[1,]) #mod.vec is your vector with randomized  values
emp.cdf<-metric.cdf(apply(w.con,1,median))#) # empirical/observed value

append(metric.cdf, emp.cdf, after = length(metric.cdf))

plot(metric.cdf,emp.cdf, xlim=c(0.01,0.20))
abline(v=emp.cdf[1])


