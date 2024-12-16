# random deletion R50
library(igraph)
library(MASS)
library(NetIndices)
library(fluxweb)
library(NetworkExtinction)
library(network)
library(dplyr)
library(patchwork)
library(tidyverse)


webs<-readRDS("index.arr.RDS")

uw.R50rand.th10<-matrix(NA, 34, 1000)
for (i in 1:1000) {
  #i<-1
  a<-webs[[1]]
  
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
    sp<-c(1:dim(fw)[1])
    order_rand<-sample(sp, , replace=F)# random order
    
    #remove basal nodes from random extinction order: autotrophs, mixotrophs, detritus
    order_rand<-setdiff(order_rand, c(which(colnames(fw)=="Autotroph"),
                                      which(colnames(fw)=="Mixotroph"),
                                      which(colnames(fw)=="Detritus")))
    
    #ext. by rand order
    rand.extinc.mw<-ExtinctionOrder(fw_i, Order = order_rand, NetworkType = "Trophic",
                                    verbose = F)
    # IS is the threshold for extinctions. i.e 0.3 leads to removal of all nodes which lose 70percent of their interaction strength
    
    #R50
    uw.R50rand[j,i]<-rand.extinc.mw$R50 ## R50 value extracted from the ext.output, R50 between 1 and 0
    
  }
}
uw.R50rand#no th


uw.R50rand.th10<-matrix(NA, 34, 1000)
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
    sp<-c(1:dim(fw)[1])
    order_rand<-sample(sp, replace=F)# random order
    
    #remove basal nodes from random extinction order: autotrophs, mixotrophs, detritus
    order_rand<-setdiff(order_rand, c(which(colnames(fw)=="Autotroph"),
                                      which(colnames(fw)=="Mixotroph"),
                                      which(colnames(fw)=="Detritus")))
    
    #ext. by rand order
    rand.extinc.mw<-ExtinctionOrder(fw_i, Order = order_rand, NetworkType = "Trophic",
                                    verbose = F, IS=0.1)
    # IS is the threshold for extinctions. i.e 0.3 leads to removal of all nodes which lose 70percent of their interaction strength
    
    #R50
    uw.R50rand.th10[j,i]<-rand.extinc.mw$R50 ## R50 value extracted from the ext.output, R50 between 1 and 0
    
  }
}
uw.R50rand.th10


uw.R50rand.th20<-matrix(NA, 34, 1000)
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
    sp<-c(1:dim(fw)[1])
    order_rand<-sample(sp, replace=F)# random order
    
    #remove basal nodes from random extinction order: autotrophs, mixotrophs, detritus
    order_rand<-setdiff(order_rand, c(which(colnames(fw)=="Autotroph"),
                                      which(colnames(fw)=="Mixotroph"),
                                      which(colnames(fw)=="Detritus")))
    
    #ext. by rand order
    rand.extinc.mw<-ExtinctionOrder(fw_i, Order = order_rand, NetworkType = "Trophic",
                                    verbose = F, IS=0.2)
    # IS is the threshold for extinctions. i.e 0.3 leads to removal of all nodes which lose 70percent of their interaction strength
    
    #R50
    uw.R50rand.th20[j,i]<-rand.extinc.mw$R50 ## R50 value extracted from the ext.output, R50 between 1 and 0
    
  }
}
uw.R50rand.th20


uw.R50rand.th30<-matrix(NA, 34, 1000)
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
    sp<-c(1:dim(fw)[1])
    order_rand<-sample(sp, replace=F)# random order
    
    #remove basal nodes from random extinction order: autotrophs, mixotrophs, detritus
    order_rand<-setdiff(order_rand, c(which(colnames(fw)=="Autotroph"),
                                      which(colnames(fw)=="Mixotroph"),
                                      which(colnames(fw)=="Detritus")))
    
    #ext. by rand order
    rand.extinc.mw<-ExtinctionOrder(fw_i, Order = order_rand, NetworkType = "Trophic",
                                    verbose = F, IS=0.3)
    # IS is the threshold for extinctions. i.e 0.3 leads to removal of all nodes which lose 70percent of their interaction strength
    
    #R50
    uw.R50rand.th30[j,i]<-rand.extinc.mw$R50 ## R50 value extracted from the ext.output, R50 between 1 and 0
    
  }
}
uw.R50rand.th30


uw.R50rand.th40<-matrix(NA, 34, 1000)
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
    sp<-c(1:dim(fw)[1])
    order_rand<-sample(sp, replace=F)# random order
    
    #remove basal nodes from random extinction order: autotrophs, mixotrophs, detritus
    order_rand<-setdiff(order_rand, c(which(colnames(fw)=="Autotroph"),
                                      which(colnames(fw)=="Mixotroph"),
                                      which(colnames(fw)=="Detritus")))
    
    #ext. by rand order
    rand.extinc.mw<-ExtinctionOrder(fw_i, Order = order_rand, NetworkType = "Trophic",
                                    verbose = F, IS=0.4)
    # IS is the threshold for extinctions. i.e 0.3 leads to removal of all nodes which lose 70percent of their interaction strength
    
    #R50
    uw.R50rand.th40[j,i]<-rand.extinc.mw$R50 ## R50 value extracted from the ext.output, R50 between 1 and 0
    
  }
}
uw.R50rand.th40


uw.R50rand.th50<-matrix(NA, 34, 1000)
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
    sp<-c(1:dim(fw)[1])
    order_rand<-sample(sp, replace=F)# random order
    
    #remove basal nodes from random extinction order: autotrophs, mixotrophs, detritus
    order_rand<-setdiff(order_rand, c(which(colnames(fw)=="Autotroph"),
                                      which(colnames(fw)=="Mixotroph"),
                                      which(colnames(fw)=="Detritus")))
    
    #ext. by rand order
    rand.extinc.mw<-ExtinctionOrder(fw_i, Order = order_rand, NetworkType = "Trophic",
                                    verbose = F, IS=0.5)
    # IS is the threshold for extinctions. i.e 0.3 leads to removal of all nodes which lose 70percent of their interaction strength
    
    #R50
    uw.R50rand.th50[j,i]<-rand.extinc.mw$R50 ## R50 value extracted from the ext.output, R50 between 1 and 0
    
  }
}
uw.R50rand.th50


uw.R50rand.th60<-matrix(NA, 34, 1000)
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
    sp<-c(1:dim(fw)[1])
    order_rand<-sample(sp, replace=F)# random order
    
    #remove basal nodes from random extinction order: autotrophs, mixotrophs, detritus
    order_rand<-setdiff(order_rand, c(which(colnames(fw)=="Autotroph"),
                                      which(colnames(fw)=="Mixotroph"),
                                      which(colnames(fw)=="Detritus")))
    
    #ext. by rand order
    rand.extinc.mw<-ExtinctionOrder(fw_i, Order = order_rand, NetworkType = "Trophic",
                                    verbose = F, IS=0.6)
    # IS is the threshold for extinctions. i.e 0.3 leads to removal of all nodes which lose 70percent of their interaction strength
    
    #R50
    uw.R50rand.th60[j,i]<-rand.extinc.mw$R50 ## R50 value extracted from the ext.output, R50 between 1 and 0
    
  }
}
uw.R50rand.th60


uw.R50rand.th70<-matrix(NA, 34, 1000)
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
    sp<-c(1:dim(fw)[1])
    order_rand<-sample(sp, replace=F)# random order
    
    #remove basal nodes from random extinction order: autotrophs, mixotrophs, detritus
    order_rand<-setdiff(order_rand, c(which(colnames(fw)=="Autotroph"),
                                      which(colnames(fw)=="Mixotroph"),
                                      which(colnames(fw)=="Detritus")))
    
    #ext. by rand order
    rand.extinc.mw<-ExtinctionOrder(fw_i, Order = order_rand, NetworkType = "Trophic",
                                    verbose = F, IS=0.7)
    # IS is the threshold for extinctions. i.e 0.3 leads to removal of all nodes which lose 70percent of their interaction strength
    
    #R50
    uw.R50rand.th70[j,i]<-rand.extinc.mw$R50 ## R50 value extracted from the ext.output, R50 between 1 and 0
    
  }
}
uw.R50rand.th70


uw.R50rand.th80<-matrix(NA, 34, 1000)
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
    sp<-c(1:dim(fw)[1])
    order_rand<-sample(sp, replace=F)# random order
    
    #remove basal nodes from random extinction order: autotrophs, mixotrophs, detritus
    order_rand<-setdiff(order_rand, c(which(colnames(fw)=="Autotroph"),
                                      which(colnames(fw)=="Mixotroph"),
                                      which(colnames(fw)=="Detritus")))
    
    #ext. by rand order
    rand.extinc.mw<-ExtinctionOrder(fw_i, Order = order_rand, NetworkType = "Trophic",
                                    verbose = F, IS=0.8)
    # IS is the threshold for extinctions. i.e 0.3 leads to removal of all nodes which lose 70% of their interaction strength
    
    #R50
    uw.R50rand.th80[j,i]<-rand.extinc.mw$R50 ## R50 value extracted from the ext.output, R50 between 1 and 0
    
  }
}
uw.R50rand.th80


uw.R50rand.th90<-matrix(NA, 34, 1000)
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
    sp<-c(1:dim(fw)[1])
    order_rand<-sample(sp,replace=F)# random order
    
    #remove basal nodes from random extinction order: autotrophs, mixotrophs, detritus
    order_rand<-setdiff(order_rand, c(which(colnames(fw)=="Autotroph"),
                                      which(colnames(fw)=="Mixotroph"),
                                      which(colnames(fw)=="Detritus")))
    
    #ext. by rand order
    rand.extinc.mw<-ExtinctionOrder(fw_i, Order = order_rand, NetworkType = "Trophic",
                                    verbose = F, IS=0.9)
    # IS is the threshold for extinctions. i.e 0.3 leads to removal of all nodes which lose 70percent of their interaction strength
    
    #R50
    uw.R50rand.th90[j,i]<-rand.extinc.mw$R50 ## R50 value extracted from the ext.output, R50 between 1 and 0
    
  }
}
uw.R50rand.th90




#all in one loop so the rand deletion sequence is same for the ths, to be comparable to each other for the ths
uw.R50rand<-matrix(NA, 34, 1000)
uw.R50rand.th10<-matrix(NA, 34, 1000)
uw.R50rand.th20<-matrix(NA, 34, 1000)
uw.R50rand.th30<-matrix(NA, 34, 1000)
uw.R50rand.th40<-matrix(NA, 34, 1000)
uw.R50rand.th50<-matrix(NA, 34, 1000)
uw.R50rand.th60<-matrix(NA, 34, 1000)
uw.R50rand.th70<-matrix(NA, 34, 1000)
uw.R50rand.th80<-matrix(NA, 34, 1000)
uw.R50rand.th90<-matrix(NA, 34, 1000)

set.seed(33985)
for (i in 1:1000) {
  a<-webs[[i]]
  
  for (j in 1:34) {
    fw<-a[[j]]
    
    fw_i<-fw
    colnames(fw_i)<-c(1:dim(fw)[1])# remove node names so extinction function can run
    rownames(fw_i)<-c(1:dim(fw)[1])# remove node names so extinction function can run
    
    sp<-c(1:dim(fw)[1])
    order_rand<-sample(sp, replace=F)# random order
    
    #remove basal nodes from random extinction order: autotrophs, mixotrophs, detritus
    order_rand<-setdiff(order_rand, c(which(colnames(fw)=="Autotroph"),
                                      which(colnames(fw)=="Mixotroph"),
                                      which(colnames(fw)=="Detritus")))
    
    #ext. by rand order, same order used for all ths so they are comparable
    rand.extinc.mw<-ExtinctionOrder(fw_i, Order = order_rand, NetworkType = "Trophic",
                                    verbose = F)
    uw.R50rand[j,i]<-rand.extinc.mw$R50 
    rand.extinc.mw<-ExtinctionOrder(fw_i, Order = order_rand, NetworkType = "Trophic",
                                    verbose = F, IS=0.1)
    uw.R50rand.th10[j,i]<-rand.extinc.mw$R50 
    rand.extinc.mw<-ExtinctionOrder(fw_i, Order = order_rand, NetworkType = "Trophic",
                                    verbose = F, IS=0.2)
    uw.R50rand.th20[j,i]<-rand.extinc.mw$R50 
    rand.extinc.mw<-ExtinctionOrder(fw_i, Order = order_rand, NetworkType = "Trophic",
                                    verbose = F, IS=0.3)
    uw.R50rand.th30[j,i]<-rand.extinc.mw$R50 
    rand.extinc.mw<-ExtinctionOrder(fw_i, Order = order_rand, NetworkType = "Trophic",
                                    verbose = F, IS=0.4)
    uw.R50rand.th40[j,i]<-rand.extinc.mw$R50 
    rand.extinc.mw<-ExtinctionOrder(fw_i, Order = order_rand, NetworkType = "Trophic",
                                    verbose = F, IS=0.5)
    uw.R50rand.th50[j,i]<-rand.extinc.mw$R50 
    rand.extinc.mw<-ExtinctionOrder(fw_i, Order = order_rand, NetworkType = "Trophic",
                                    verbose = F, IS=0.6)
    uw.R50rand.th60[j,i]<-rand.extinc.mw$R50 
    rand.extinc.mw<-ExtinctionOrder(fw_i, Order = order_rand, NetworkType = "Trophic",
                                    verbose = F, IS=0.7)
    uw.R50rand.th70[j,i]<-rand.extinc.mw$R50 
    rand.extinc.mw<-ExtinctionOrder(fw_i, Order = order_rand, NetworkType = "Trophic",
                                    verbose = F, IS=0.8)
    uw.R50rand.th80[j,i]<-rand.extinc.mw$R50 
    rand.extinc.mw<-ExtinctionOrder(fw_i, Order = order_rand, NetworkType = "Trophic",
                                    verbose = T, IS=0.9) #verbose=T here so you can se that the loop is running, it just prints progressbar
    uw.R50rand.th90[j,i]<-rand.extinc.mw$R50 
  }
}

uw.rand_del_th_list<-list(uw.R50rand,
                        uw.R50rand.th10,
                        uw.R50rand.th20,
                        uw.R50rand.th30,
                        uw.R50rand.th40,
                        uw.R50rand.th50,
                        uw.R50rand.th60,
                        uw.R50rand.th70,
                        uw.R50rand.th80,
                        uw.R50rand.th90)
saveRDS(uw.rand_del_th_list,"rand_del_all.RDS")
uw.rand_del_th_list

m<-matrix(NA,1000,34)
colnames(m)<-c(1981:2014)
row.names(m)<-c(1:1000)
mlist<-list()
for (j in 1:10){
  a<-uw.rand_del_th_list[[j]]
  for (i in 1:34) {
    m[,i]<-a[i,]
  }
  mlist[[j]]<-m
}

THlist<-list("binary","10%", "20%","30%","40%","50%", "60%", "70%", "80%", "90%")
df.R50list<-list()
for (i in 1:10) {
  m2<-as.data.frame(as.table(mlist[[i]]))
  durp<-as_tibble(m2)
  durp$Var2<-as.numeric(as.character(durp$Var2))
  colnames(durp)<-c("iter", "year", THlist[[i]][1])
  df.R50list[[i]]<-durp
}

Con<-readRDS("Con_34000.RDS")
m3<-matrix(NA,1000,34)
colnames(m3)<-c(1981:2014)
row.names(m3)<-c(1:1000)
for (i in 1:34) {
  m3[,i]<-Con[i,]
}
m3
m4<-as.data.frame(as.table(m3))
durp2<-as_tibble(m4)
durp2$Var2<-as.numeric(as.character(durp2$Var2))
colnames(durp2)<-c("iter", "year", "Con")

urp<-full_join(durp2, df.R50list[[1]])
for (i in 1:10) {
  urp<-left_join(urp, df.R50list[[i]])
}
urp
urp3<-pivot_longer(urp, 4:13)
colnames(urp3)[4]<-c("threshold")
colnames(urp3)[5]<-c("R50")
urp3

ggplot(urp3, aes(x=Con, y=R50, colour=year))+
  geom_point()+
  stat_smooth(method = "glm",
              formula = "y~x",
              geom = "smooth",
              col="Red")+
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, size = 18))+
  theme(axis.text.y = element_text(angle = 0, vjust = 0.5, size = 18))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  scale_x_continuous(labels = scales::number_format(accuracy = 0.01))+
  scale_color_gradient(low="yellow", high="Blue")+
  theme(plot.title = element_text(hjust = 0))+xlab("Connectance")+ylab("R50")+
  theme(strip.text.x = element_text(size = 15))+
  theme(axis.title=element_text(size=28))+theme(legend.text=element_text(size=25),
                                                legend.title=element_text(size=25))+
  facet_wrap(~factor(threshold, levels=c("binary","10%", "20%","30%","40%",
                                         "50%", "60%", "70%", "80%", "90%")),ncol=5)+
  theme(strip.background = element_blank())+theme(panel.border=element_rect(fill=NA))+
  theme(strip.text.x = element_text(size = 18))+
  theme(panel.spacing.x = unit(1, "lines"))+
  ggtitle("Random deletions unweighted webs")+
  theme(plot.title = element_text(size=25))

df.r50<-matrix(NA,34,7)
colnames(df.r50)<-c("Median","upperCI95","lowerCI95","upperCI50","lowerCI50","Year","threshold")
df.r50[,6]<-c(1981:2014)
dflist<-list()
THlist<-list("binary","10%", "20%","30%","40%","50%", "60%", "70%", "80%", "90%")
for (j in 1:10){
  a<-uw.rand_del_th_list[[j]]
  df.r50[,7]<-THlist[[j]][1]
  #a[1,]
  for (i in 1:34){
    df.r50[i,1]<-median(a[i,])
    df.r50[i,2]<-quantile(a[i,], probs = 0.975)
    df.r50[i,3]<-quantile(a[i,], probs = 0.05/2)
    df.r50[i,4]<-quantile(a[i,], probs = 0.75)
    df.r50[i,5]<-quantile(a[i,], probs = 0.25/2)
    dflist[[j]]<-df.r50
  }
}
dflist[[1]]
df.r50.time<-as.data.frame(dflist[[1]])
for (i in 1:10) {
  df.r50.time<-full_join(df.r50.time, as.data.frame(dflist[[i]]))
}
df.r50.time
#df.r50.time<-pivot_longer(df.r50.time, 1:5)

df.r50.time$Median<-as.numeric(as.character(df.r50.time$Median))
df.r50.time$Year<-as.numeric(as.character(df.r50.time$Year))
df.r50.time$upperCI95<-as.numeric(as.character(df.r50.time$upperCI95))
df.r50.time$lowerCI95<-as.numeric(as.character(df.r50.time$lowerCI95))
df.r50.time$lowerCI50<-as.numeric(as.character(df.r50.time$lowerCI50))
df.r50.time$upperCI50<-as.numeric(as.character(df.r50.time$upperCI50))

df.r50.time<-pivot_longer(df.r50.time, 6)
df.r50.time$value<-as.numeric(as.character(df.r50.time$value))
#filter(df.r50.time, name == "Median")

ggplot(df.r50.time, aes(x=value, y=Median, colour=threshold, fill=threshold))+
  geom_line(linewidth=1.1)+
  theme(axis.text.x = element_text(angle = 45,vjust = 0.5, size = 16))+
  geom_ribbon(aes(ymin=upperCI95, ymax=lowerCI95),alpha=0.3,colour = NA)+
  geom_ribbon(aes(ymin=upperCI50, ymax=lowerCI50),alpha=0.5,colour = NA)+
  theme(axis.text.y = element_text(angle = 0, vjust = 0.5, size = 16))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  scale_x_continuous(name = "Year", breaks = scales::pretty_breaks(n = 5))+
  xlab("Year")+ylab("R50")+theme(axis.title=element_text(size=25))+
  theme(legend.text=element_text(size=25),legend.title=element_text(size=25))+
  theme(legend.title = element_text("Thresholds"))+theme(legend.position = "none")+
  facet_wrap(~factor(threshold, levels=c("binary","10%", "20%","30%","40%",
                                         "50%", "60%", "70%", "80%", "90%")),ncol=5)+
  theme(strip.text.x = element_text(size = 16))+
  theme(strip.background = element_blank())+theme(panel.border=element_rect(fill=NA))+
  scale_color_brewer(palette = "Spectral", direction =1)+scale_fill_brewer(palette = "Spectral", direction =1)+
  ggtitle("Random deletions unweighted webs")+
  theme(plot.title = element_text(size=25))
# scale_fill_manual(values=c("green","orange","blue","black","red","yellow","grey","orange","dark green","dark blue"))



####weighted rand del####
w.R50rand<-matrix(NA, 34, 1000)
w.R50rand.th10<-matrix(NA, 34, 1000)
w.R50rand.th20<-matrix(NA, 34, 1000)
w.R50rand.th30<-matrix(NA, 34, 1000)
w.R50rand.th40<-matrix(NA, 34, 1000)
w.R50rand.th50<-matrix(NA, 34, 1000)
w.R50rand.th60<-matrix(NA, 34, 1000)
w.R50rand.th70<-matrix(NA, 34, 1000)
w.R50rand.th80<-matrix(NA, 34, 1000)
w.R50rand.th90<-matrix(NA, 34, 1000)

set.seed(33985)
for (i in 1:1000) {
  a<-w.webs.1000[[i]]
  
  for (j in 1:34) {
    fw<-a[[j]]
    
    fw_i<-fw
    colnames(fw_i)<-c(1:dim(fw)[1])# remove node names so extinction function can run
    rownames(fw_i)<-c(1:dim(fw)[1])# remove node names so extinction function can run
    
    sp<-c(1:dim(fw)[1]) #sp in vector to randomize
    order_rand<-sample(sp, replace=F)# random order
    
    #remove basal nodes from random extinction order: autotrophs, mixotrophs, detritus
    order_rand<-setdiff(order_rand, c(which(colnames(fw)=="Autotroph"),
                                      which(colnames(fw)=="Mixotroph"),
                                      which(colnames(fw)=="Detritus")))
    
    #ext. by rand order, same order used for all ths so they are comparable
    rand.extinc.mw<-ExtinctionOrder(fw_i, Order = order_rand, NetworkType = "Trophic",
                                    verbose = F)
    w.R50rand[j,i]<-rand.extinc.mw$R50
    rand.extinc.mw<-ExtinctionOrder(fw_i, Order = order_rand, NetworkType = "Trophic",
                                    verbose = F, IS=0.1)
    w.R50rand.th10[j,i]<-rand.extinc.mw$R50
    rand.extinc.mw<-ExtinctionOrder(fw_i, Order = order_rand, NetworkType = "Trophic",
                                    verbose = F, IS=0.2)
    w.R50rand.th20[j,i]<-rand.extinc.mw$R50
    rand.extinc.mw<-ExtinctionOrder(fw_i, Order = order_rand, NetworkType = "Trophic",
                                    verbose = F, IS=0.3)
    w.R50rand.th30[j,i]<-rand.extinc.mw$R50
    rand.extinc.mw<-ExtinctionOrder(fw_i, Order = order_rand, NetworkType = "Trophic",
                                    verbose = F, IS=0.4)
    w.R50rand.th40[j,i]<-rand.extinc.mw$R50
    rand.extinc.mw<-ExtinctionOrder(fw_i, Order = order_rand, NetworkType = "Trophic",
                                    verbose = F, IS=0.5)
    w.R50rand.th50[j,i]<-rand.extinc.mw$R50
    rand.extinc.mw<-ExtinctionOrder(fw_i, Order = order_rand, NetworkType = "Trophic",
                                    verbose = F, IS=0.6)
    w.R50rand.th60[j,i]<-rand.extinc.mw$R50
    rand.extinc.mw<-ExtinctionOrder(fw_i, Order = order_rand, NetworkType = "Trophic",
                                    verbose = F, IS=0.7)
    w.R50rand.th70[j,i]<-rand.extinc.mw$R50 
    rand.extinc.mw<-ExtinctionOrder(fw_i, Order = order_rand, NetworkType = "Trophic",
                                    verbose = F, IS=0.8)
    w.R50rand.th80[j,i]<-rand.extinc.mw$R50 
    rand.extinc.mw<-ExtinctionOrder(fw_i, Order = order_rand, NetworkType = "Trophic",
                                    verbose = T, IS=0.9) #verbose=T here so you can se that the loop is running, it just prints progressbar
    w.R50rand.th90[j,i]<-rand.extinc.mw$R50 
  }
}

w.rand_del_th_list<-list(w.R50rand,
                       w.R50rand.th10,
                       w.R50rand.th20,
                       w.R50rand.th30,
                       w.R50rand.th40,
                       w.R50rand.th50,
                       w.R50rand.th60,
                       w.R50rand.th70,
                       w.R50rand.th80,
                       w.R50rand.th90)

#saveRDS(w.rand_del_th_list,"w_rand_del_all.RDS")


m<-matrix(NA,1000,34)
colnames(m)<-c(1981:2014)
row.names(m)<-c(1:1000)
mlist<-list()
for (j in 1:10){
  a<-w.rand_del_th_list[[1]]
  for (i in 1:34) {
    m[,i]<-a[i,]
  }
  mlist[[j]]<-m
}

THlist<-list("binary","10%", "20%","30%","40%","50%", "60%", "70%", "80%", "90%")
df.R50list<-list()
for (i in 1:10) {
  m2<-as.data.frame(as.table(mlist[[i]]))
  durp<-as_tibble(m2)
  durp$Var2<-as.numeric(as.character(durp$Var2))
  colnames(durp)<-c("iter", "year", THlist[[i]][1])
  df.R50list[[i]]<-durp
}

w.con<-readRDS("W.con_34000.RDS")
m3<-matrix(NA,1000,34)
colnames(m3)<-c(1981:2014)
row.names(m3)<-c(1:1000)
for (i in 1:34) {
  m3[,i]<-w.con[i,]
}
m3
m4<-as.data.frame(as.table(m3))
durp2<-as_tibble(m4)
durp2$Var2<-as.numeric(as.character(durp2$Var2))
colnames(durp2)<-c("iter", "year", "w.con")

urp<-full_join(durp2, df.R50list[[1]])
for (i in 1:10) {
  urp<-left_join(urp, df.R50list[[i]])
}
urp
urp2<-pivot_longer(urp, 4:13)
colnames(urp2)[4]<-c("threshold")
colnames(urp2)[5]<-c("R50")
urp2

ggplot(urp2, aes(x=w.con, y=R50, colour=year))+
  geom_point()+
  stat_smooth(method = "glm",
              formula = "y~x",
              geom = "smooth",
              col="Red")+
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, size = 18))+
  theme(axis.text.y = element_text(angle = 0, vjust = 0.5, size = 18))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  scale_x_continuous(labels = scales::number_format(accuracy = 0.01))+
  scale_color_gradient(low="yellow", high="Blue")+
  theme(plot.title = element_text(hjust = 0))+xlab("Weighted Connectance")+ylab("R50")+
  theme(strip.text.x = element_text(size = 15))+
  theme(axis.title=element_text(size=28))+theme(legend.text=element_text(size=25),
                                                legend.title=element_text(size=25))+
  facet_wrap(~factor(threshold, levels=c("binary","10%", "20%","30%","40%",
                                         "50%", "60%", "70%", "80%", "90%")),ncol=5)+
  theme(strip.background = element_blank())+theme(panel.border=element_rect(fill=NA))+
  theme(strip.text.x = element_text(size = 18))+
  ggtitle("Random deletions weighted webs")+
  theme(plot.title = element_text(size=25))
#ncol=5, scale_color_gradient(low="yellow", high="Blue")+scale_color_viridis(name = "Year")+


df.r50<-matrix(NA,34,7)
colnames(df.r50)<-c("Median","upperCI95","lowerCI95","upperCI50","lowerCI50","Year","threshold")
df.r50[,6]<-c(1981:2014)
dflist<-list()
THlist<-list("binary","10%", "20%","30%","40%","50%", "60%", "70%", "80%", "90%")
for (j in 1:10){
  a<-w.rand_del_th_list[[j]]
  df.r50[,7]<-THlist[[j]][1]
  #a[1,]
  for (i in 1:34){
    df.r50[i,1]<-median(a[i,])
    df.r50[i,2]<-quantile(a[i,], probs = 0.975)
    df.r50[i,3]<-quantile(a[i,], probs = 0.05/2)
    df.r50[i,4]<-quantile(a[i,], probs = 0.75)
    df.r50[i,5]<-quantile(a[i,], probs = 0.25/2)
    dflist[[j]]<-df.r50
  }
}
dflist[[1]]
df.r50.time<-as.data.frame(dflist[[1]])
for (i in 1:10) {
  df.r50.time<-full_join(df.r50.time, as.data.frame(dflist[[i]]))
}
df.r50.time
#df.r50.time<-pivot_longer(df.r50.time, 1:5)



df.r50.time$Median<-as.numeric(as.character(df.r50.time$Median))
df.r50.time$Year<-as.numeric(as.character(df.r50.time$Year))
df.r50.time$upperCI95<-as.numeric(as.character(df.r50.time$upperCI95))
df.r50.time$lowerCI95<-as.numeric(as.character(df.r50.time$lowerCI95))
df.r50.time$lowerCI50<-as.numeric(as.character(df.r50.time$lowerCI50))
df.r50.time$upperCI50<-as.numeric(as.character(df.r50.time$upperCI50))


df.r50.time<-pivot_longer(df.r50.time, 6)
df.r50.time$value<-as.numeric(as.character(df.r50.time$value))

ggplot(df.r50.time, aes(x=value, y=Median, colour=threshold, fill=threshold))+
  geom_line(linewidth=1.1)+
  theme(axis.text.x = element_text(angle = 45,vjust = 0.5, size = 16))+
  geom_ribbon(aes(ymin=upperCI95, ymax=lowerCI95),alpha=0.3,colour = NA)+
  geom_ribbon(aes(ymin=upperCI50, ymax=lowerCI50),alpha=0.5,colour = NA)+
  theme(axis.text.y = element_text(angle = 0, vjust = 0.5, size = 16))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  scale_x_continuous(name = "Year", breaks = scales::pretty_breaks(n = 5))+
  xlab("Year")+ylab("R50")+theme(axis.title=element_text(size=25))+
  theme(legend.text=element_text(size=25),legend.title=element_text(size=25))+
  theme(legend.title = element_text("Thresholds"))+theme(legend.position = "none")+
  facet_wrap(~factor(threshold, levels=c("binary","10%", "20%","30%","40%",
                                         "50%", "60%", "70%", "80%", "90%")),ncol=5)+
  theme(strip.text.x = element_text(size = 16))+
  theme(strip.background = element_blank())+theme(panel.border=element_rect(fill=NA))+
  scale_color_brewer(palette = "Spectral", direction =1)+scale_fill_brewer(palette = "Spectral", direction =1)+
  ggtitle("Random deletions weighted webs")+
  theme(plot.title = element_text(size=25))
# scale_fill_manual(values=c("green","orange","blue","black","red","yellow","grey","orange","dark green","dark blue"))



### Rand deletions A/C###
