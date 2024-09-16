####### LOW TO HIGH DELETION SEQUENCES #######

#Script for reviewers of the ms "Food web robustness depends on network type and threshold for extinction"

####READ ME####
#This file contains the script to produce the results from the SI of the manuscript for the low to high deletion sequence
#This script assumes you have run the loop for connectance, weighted connectance and relative ascendency in the "Script for reviewers.R"

#### packages ####

install.packages(c("igraph", "ggplot2","MASS", "NetIndecies", "fluxweb", "NetworkExtinctions",
                   "network", "dplyr", "patchwork", "scales", "tidyverse", "fmsb"))
library(ggplot2)
library(igraph)
library(MASS)
library(NetIndices)
library(fluxweb)
library(NetworkExtinction)
library(network)
library(dplyr)
library(patchwork)
library(scales)
library(tidyverse)
library(fmsb)

#### Data ####
#Food webs used in the simulations stored in the files "uw.webs.1000.RDS" and "w.webs.1000.RDS"

#read data
uw.webs.1000<-readRDS("uw.webs.1000.RDS")#the unweighted webs, format is list containing interaction matrices

w.webs.1000<-readRDS("w.webs.1000.RDS")#the weighted webs, format is list containing interaction matrices

#The lists have a length if 1000, each entry has 34 webs inside, one for each year in the time series. So 34000 food webs total
#uw.webs.1000[[1]][[1]] is the first of 1000 un-weighted web for the year 1980, uw.webs.1000[[1]][[2]] is the first un-weighted web for 1981 and so on.
#same for the weighted webs


uw.R50low_high<-matrix(NA, 34, 1000)
uw.R50low_high.th10<-matrix(NA, 34, 1000)
uw.R50low_high.th20<-matrix(NA, 34, 1000)
uw.R50low_high.th30<-matrix(NA, 34, 1000)
uw.R50low_high.th40<-matrix(NA, 34, 1000)
uw.R50low_high.th50<-matrix(NA, 34, 1000)
uw.R50low_high.th60<-matrix(NA, 34, 1000)
uw.R50low_high.th70<-matrix(NA, 34, 1000)
uw.R50low_high.th80<-matrix(NA, 34, 1000)
uw.R50low_high.th90<-matrix(NA, 34, 1000)


#This loop goes through all 34000 webs, and simulating deletions at different thresholds
#NOTE: Has a very long run-time.  
#See file "uw.low_high_R50.RDS" for R50 results to compare output.

#IMPORTANT NOTE: Will give several warnings at the end! This does not mean there was an error. The loop works as intended.
#For example, the following Warning messages show up:
#Warning messages:
#: In network::as.matrix.network.adjacency(Temp, attrname = "weight"):
#There is no edge attribute named weight" 
#NOTE:This one comes up when the webs are disassembled to the point that 
#only the basal nodes are left, the basal nodes are not connected, therefore there are no links that have link weights

#": In ExtinctionOrder(fw_i, Order = order, NetworkType = "Trophic",  :
#Your network became completely unconnected before all primary extinctions were simulated. This happened at extinction step X out of Y"
#NOTE: Self explanatory. This happens when extinctions trigger cascades and disassembled the network completely.

#": In ExtinctionOrder(fw_i, Order = order, NetworkType = "Trophic",  :
#Primary extinctions of X, Y skipped due to their prior extinction as secondary extinctions."
#NOTE: Self explanatory. The extinction already happend, so was skipped in the order.


#### Unweighted webs ####
webs<-uw.webs.1000
for (i in 1:1000){
  a<-webs[[i]]
  
  for (j in 1:34) {
    fw<-a[[j]]
    
    fw_i<-fw
    colnames(fw_i)<-c(1:dim(fw)[1])# remove node names so extinction function can run
    rownames(fw_i)<-c(1:dim(fw)[1])# remove node names so extinction function can run
    
  
    deg.sp<-apply(fw_i, 1, sum)+apply(fw_i, 2, sum)
    order_low_high<-as.numeric(names(sort(deg.sp, decreasing = F)))
    
    #remove basal nodes from random extinction order: autotrophs, mixotrophs, detritus
    order_low_high<-setdiff(order_low_high, c(which(colnames(fw)=="Autotroph"),
                                              which(colnames(fw)=="Mixotroph"),
                                              which(colnames(fw)=="Detritus")))
    
    #ext. by low to high degree order
    low_high_ext<-ExtinctionOrder(fw_i, Order = order_low_high, NetworkType = "Trophic",
                                  verbose = F)
    uw.R50low_high[j,i]<-low_high_ext$R50 
    low_high_ext<-ExtinctionOrder(fw_i, Order = order_low_high, NetworkType = "Trophic",
                                  verbose = F, IS=0.1)
    uw.R50low_high.th90[j,i]<-low_high_ext$R50 
    low_high_ext<-ExtinctionOrder(fw_i, Order = order_low_high, NetworkType = "Trophic",
                                  verbose = F, IS=0.2)
    uw.R50low_high.th80[j,i]<-low_high_ext$R50 
    low_high_ext<-ExtinctionOrder(fw_i, Order = order_low_high, NetworkType = "Trophic",
                                  verbose = F, IS=0.3)
    uw.R50low_high.th70[j,i]<-low_high_ext$R50 
    low_high_ext<-ExtinctionOrder(fw_i, Order = order_low_high, NetworkType = "Trophic",
                                  verbose = F, IS=0.4)
    uw.R50low_high.th60[j,i]<-low_high_ext$R50 
    low_high_ext<-ExtinctionOrder(fw_i, Order = order_low_high, NetworkType = "Trophic",
                                  verbose = F, IS=0.5)
    uw.R50low_high.th50[j,i]<-low_high_ext$R50 
    low_high_ext<-ExtinctionOrder(fw_i, Order = order_low_high, NetworkType = "Trophic",
                                  verbose = F, IS=0.6)
    uw.R50low_high.th40[j,i]<-low_high_ext$R50 
    low_high_ext<-ExtinctionOrder(fw_i, Order = order_low_high, NetworkType = "Trophic",
                                  verbose = F, IS=0.7)
    uw.R50low_high.th30[j,i]<-low_high_ext$R50 
    low_high_ext<-ExtinctionOrder(fw_i, Order = order_low_high, NetworkType = "Trophic",
                                  verbose = F, IS=0.8)
    uw.R50low_high.th20[j,i]<-low_high_ext$R50 
    low_high_ext<-ExtinctionOrder(fw_i, Order = order_low_high, NetworkType = "Trophic",
                                  verbose = T, IS=0.9) #verbose=T here so you can se that the loop is running, it just prints progressbar
    uw.R50low_high.th10[j,i]<-low_high_ext$R50 
  }
}

uw.low_high_del_th_list<-list(uw.R50low_high,
                              uw.R50low_high.th90,
                              uw.R50low_high.th80,
                              uw.R50low_high.th70,
                              uw.R50low_high.th60,
                              uw.R50low_high.th50,
                              uw.R50low_high.th40,
                              uw.R50low_high.th30,
                              uw.R50low_high.th20,
                              uw.R50low_high.th10)

#Make data frame
#the following is a series of loops that produce and join data frames into one data frame for use in plots and analysis
#assumes you have calculated connectance, 
#weighted connectance and relative ascendency for the webs in the script "Main_results.R"
m<-matrix(NA,1000,34)
colnames(m)<-c(1981:2014)
row.names(m)<-c(1:1000)
mlist<-list()
for (j in 1:10){
  a<-uw.low_high_del_th_list[[j]]
  for (i in 1:34) {
    m[,i]<-a[i,]
  }
  mlist[[j]]<-m
}

THlist<-list("Links/no links","90%","80%","70%","60%","50%", "40%", "30%", "20%", "10%")
df.R50list<-list()
for (i in 1:10) {
  m2<-as.data.frame(as.table(mlist[[i]]))
  d<-as_tibble(m2)
  d$Var2<-as.numeric(as.character(d$Var2))
  colnames(d)<-c("iter", "year", THlist[[i]][1])
  df.R50list[[i]]<-d
}

m3<-matrix(NA,1000,34)
colnames(m3)<-c(1981:2014)
row.names(m3)<-c(1:1000)
for (i in 1:34) {
  m3[,i]<-w.con[i,]
}

m4<-as.data.frame(as.table(m3))
d2<-as_tibble(m4)
d2$Var2<-as.numeric(as.character(d2$Var2))
colnames(d2)<-c("iter", "year", "w.con")

u<-full_join(d2, df.R50list[[1]])
for (i in 1:10) {
  u<-left_join(u, df.R50list[[i]])
}
u

m3<-matrix(NA,1000,34)
colnames(m3)<-c(1981:2014)
row.names(m3)<-c(1:1000)
for (i in 1:34) {
  m3[,i]<-rel.asc[i,]
}
m3
m4<-as.data.frame(as.table(m3))
AC<-as_tibble(m4)
AC$Var2<-as.numeric(as.character(AC$Var2))
colnames(AC)<-c("iter", "year", "A/C")

uw.DF.low_high<-right_join(u, AC)

uw.DF.low_high.plots<-pivot_longer(uw.DF.low_high, 4:13)# make longformat df for use in ggplot
colnames(uw.DF.low_high.plots)[5]<-c("threshold")
colnames(uw.DF.low_high.plots)[6]<-c("R50")

uw.DF.low_high# data frame used for the correlation analysis
uw.DF.low_high.plots#long format data frame used for ggplot


#### weighted webs ####
#This loop goes through all 34000 webs, and simulating deletions at different thresholds
#NOTE: Has a very long run-time.  
#See file "w.low_high_R50.RDS" for R50 results to compare output.


#IMPORTANT NOTE: Will give several warnings at the end! This does not mean there was an error.
#For example, the following Warning messages show up:

#"1: In network::as.matrix.network.adjacency(Temp, attrname = "weight"):
#There is no edge attribute named weight"
#NOTE:This one comes up when the weighted webs are disassembled to the point that 
#only the basal nodes are left, the basal nodes are not connected, therefore there are no links that have link weights

#"2: In ExtinctionOrder(fw_i, Order = order, NetworkType = "Trophic",  :
#Your network became completely unconnected before all primary extinctions were simulated. This happened at extinction step X out of Y"
#NOTE: Self explanatory. This happens when extinctions trigger cascades and disassembled the network completely.

#"3: In ExtinctionOrder(fw_i, Order = order, NetworkType = "Trophic",  :
#Primary extinctions of X, Y skipped due to their prior extinction as secondary extinctions."
#NOTE: Self explanatory. The extinction already happend, so was skipped in the order.

#weighted webs
webs<-w.webs.1000
for (i in 1:1000){
  a<-webs[[i]]
  
  for (j in 1:34) {
    fw<-a[[j]]
    
    fw_i<-fw
    colnames(fw_i)<-c(1:dim(fw)[1])# remove node names so extinction function can run
    rownames(fw_i)<-c(1:dim(fw)[1])# remove node names so extinction function can run
    

    deg.sp<-apply(fw_i, 1, sum)+apply(fw_i, 2, sum)
    order_low_high<-as.numeric(names(sort(deg.sp, decreasing = F)))
    
    #remove basal nodes from random extinction order: autotrophs, mixotrophs, detritus
    order_low_high<-setdiff(order_low_high, c(which(colnames(fw)=="Autotroph"),
                                              which(colnames(fw)=="Mixotroph"),
                                              which(colnames(fw)=="Detritus")))
    
    #ext. by low to high degree order
    low_high_ext<-ExtinctionOrder(fw_i, Order = order_low_high, NetworkType = "Trophic",
                                  verbose = F)
    w.R50low_high[j,i]<-low_high_ext$R50 
    low_high_ext<-ExtinctionOrder(fw_i, Order = order_low_high, NetworkType = "Trophic",
                                  verbose = F, IS=0.1)
    w.R50low_high.th90[j,i]<-low_high_ext$R50 
    low_high_ext<-ExtinctionOrder(fw_i, Order = order_low_high, NetworkType = "Trophic",
                                  verbose = F, IS=0.2)
    w.R50low_high.th80[j,i]<-low_high_ext$R50 
    low_high_ext<-ExtinctionOrder(fw_i, Order = order_low_high, NetworkType = "Trophic",
                                  verbose = F, IS=0.3)
    w.R50low_high.th70[j,i]<-low_high_ext$R50 
    low_high_ext<-ExtinctionOrder(fw_i, Order = order_low_high, NetworkType = "Trophic",
                                  verbose = F, IS=0.4)
    w.R50low_high.th60[j,i]<-low_high_ext$R50 
    low_high_ext<-ExtinctionOrder(fw_i, Order = order_low_high, NetworkType = "Trophic",
                                  verbose = F, IS=0.5)
    w.R50low_high.th50[j,i]<-low_high_ext$R50 
    low_high_ext<-ExtinctionOrder(fw_i, Order = order_low_high, NetworkType = "Trophic",
                                  verbose = F, IS=0.6)
    w.R50low_high.th40[j,i]<-low_high_ext$R50 
    low_high_ext<-ExtinctionOrder(fw_i, Order = order_low_high, NetworkType = "Trophic",
                                  verbose = F, IS=0.7)
    w.R50low_high.th30[j,i]<-low_high_ext$R50 
    low_high_ext<-ExtinctionOrder(fw_i, Order = order_low_high, NetworkType = "Trophic",
                                  verbose = F, IS=0.8)
    w.R50low_high.th20[j,i]<-low_high_ext$R50 
    low_high_ext<-ExtinctionOrder(fw_i, Order = order_low_high, NetworkType = "Trophic",
                                  verbose = T, IS=0.9) #verbose=T here so you can se that the loop is running, it just prints progressbar
    w.R50low_high.th10[j,i]<-low_high_ext$R50 
  }
}

w.low_high_del_th_list<-list(w.R50low_high,
                              w.R50low_high.th90,
                              w.R50low_high.th80,
                              w.R50low_high.th70,
                              w.R50low_high.th60,
                              w.R50low_high.th50,
                              w.R50low_high.th40,
                              w.R50low_high.th30,
                              w.R50low_high.th20,
                              w.R50low_high.th10)


#Make data frame
#the following is a series of loops that produce and join data frames into one data frame for use in plots and analysis
#assumes you have calculated connectance, 
#weighted connectance and relative ascendency for the webs in the script "Main_results.R"
m<-matrix(NA,1000,34)
colnames(m)<-c(1981:2014)
row.names(m)<-c(1:1000)
mlist<-list()
for (j in 1:10){
  a<-w.low_high_del_th_list[[j]]
  for (i in 1:34) {
    m[,i]<-a[i,]
  }
  mlist[[j]]<-m
}

THlist<-list("Links/no links","90%","80%","70%","60%","50%", "40%", "30%", "20%", "10%")
df.R50list<-list()
for (i in 1:10) {
  m2<-as.data.frame(as.table(mlist[[i]]))
  d<-as_tibble(m2)
  d$Var2<-as.numeric(as.character(d$Var2))
  colnames(d)<-c("iter", "year", THlist[[i]][1])
  df.R50list[[i]]<-d
}

m3<-matrix(NA,1000,34)
colnames(m3)<-c(1981:2014)
row.names(m3)<-c(1:1000)
for (i in 1:34) {
  m3[,i]<-w.con[i,]
}

m4<-as.data.frame(as.table(m3))
d2<-as_tibble(m4)
d2$Var2<-as.numeric(as.character(d2$Var2))
colnames(d2)<-c("iter", "year", "w.con")

u<-full_join(d2, df.R50list[[1]])
for (i in 1:10) {
  u<-left_join(u, df.R50list[[i]])
}
u

m3<-matrix(NA,1000,34)
colnames(m3)<-c(1981:2014)
row.names(m3)<-c(1:1000)
for (i in 1:34) {
  m3[,i]<-w.rel.asc[i,]
}
m3
m4<-as.data.frame(as.table(m3))
AC<-as_tibble(m4)
AC$Var2<-as.numeric(as.character(AC$Var2))
colnames(AC)<-c("iter", "year", "A/C")

w.DF.low_high<-right_join(u, AC)

w.DF.low_high.plots<-pivot_longer(w.DF.low_high, 4:13)# make longformat df for use in ggplot
colnames(w.DF.low_high.plots)[5]<-c("threshold")
colnames(w.DF.low_high.plots)[6]<-c("R50")

w.DF.low_high# dataframe used for the correlation analysis
w.DF.low_high.plots#longformat dataframe used for ggplot


#### SI fig. 3 plot ####
df.allTH<-matrix(NA,34,12)
df.allTH[,11]<-c(seq(1981,2014))
df.allTH[,12]<-c("a")
colnames(df.allTH)<-c("Links/no links","90%","80%", "70%", "60%",
                      "50%","40%", "30%", "20%", "10%", "Year", "webs")
for (i in 1:10){
  df.allTH[,i]<-apply(uw.low_high_del_th_list[[i]],1, median)
}
as.data.frame(df.allTH)
d.uw.lh<-pivot_longer(as.data.frame(df.allTH), cols = 1:10)

df.allTH<-matrix(NA,34,12)
df.allTH[,11]<-c(seq(1981,2014))
df.allTH[,12]<-c("b")
colnames(df.allTH)<-c("Links/no links","90%","80%", "70%", "60%",
                      "50%","40%", "30%", "20%", "10%", "Year", "webs")
for (i in 1:10){
  df.allTH[,i]<-apply(w.low_high_del_th_list[[i]],1, median)
}
as.data.frame(df.allTH)
d.w.lh<-pivot_longer(as.data.frame(df.allTH), cols = 1:10)

#combine
d.comb.lh<-full_join(d.uw.lh, d.w.lh)
d.comb.lh$Year<-as.numeric(as.character(d.comb.lh$Year))
d.comb.lh$value<-as.numeric(as.character(d.comb.lh$value))
d.comb.lh
d.comb.lh$name<-factor(d.comb.lh$name, levels = c("Links/no links","90%","80%", "70%", "60%",
                                                  "50%","40%", "30%", "20%", "10%"))

si.fig3<-ggplot(d.comb.lh, aes(x=Year, y=value, color = name))+
  geom_line(aes(group=name),linewidth=1.1)+
  theme(axis.text.x = element_text(angle = -90, vjust = 0.5, size = 25))+
  theme(axis.text.y = element_text(angle = 0,vjust = 0.5, size = 25))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  xlab("Year")+ylab("R50")+theme(axis.title=element_text(size=28), axis.title.x = element_text(angle = 0),
                                 axis.title.y.left = element_text(angle = 90))+
  ylim(0,0.5)+
  theme(axis.title.y = element_text(angle = 0, vjust = 0.5))+
  scale_color_brewer(palette="Spectral")+
  facet_wrap(~factor(webs),ncol=2)+labs(color = "Thresholds")+
  theme(strip.background = element_blank())+theme(panel.border=element_rect(fill=NA))+
  theme(strip.text.x = element_text(size = 18,face="bold", hjust = -0.01))+
  theme(legend.text=element_text(size=25))+theme(legend.title=element_text(size=25))+
  coord_cartesian(clip = "off")

si.fig3


#### SI fig 5 plot ####

si.fig5a<-ggplot(uw.DF.low_high.plots, aes(x=con, y=R50, colour=year))+
  geom_point()+
  stat_smooth(method = "glm",
              formula = "y~x",
              geom = "smooth",
              col="Red")+
  theme(axis.text.x = element_text(angle = -90, vjust = 0.5, size = 17))+
  theme(axis.text.y = element_text(angle = 0, vjust = 0.5, size = 17))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  scale_x_continuous(labels = scales::number_format(accuracy = 0.01))+
  ylim(0, 0.5)+
  scale_color_gradient(low="yellow", high="Blue")+
  theme(plot.title = element_text(hjust = 0))+xlab("Connectance")+ylab("R50")+
  theme(strip.text.x = element_text(size = 15))+
  theme(axis.title=element_text(size=22))+theme(legend.text=element_text(size=22),
                                                legend.title=element_text(size=22))+
  facet_wrap(~factor(threshold, levels=c("Links/no links","90%", "80%","70%","60%",
                                         "50%", "40%", "30%", "20%", "10%")),ncol=5)+
  theme(strip.background = element_blank())+theme(panel.border=element_rect(fill=NA))+
  theme(strip.text.x = element_text(size = 12))


si.fig5b<-ggplot(w.DF.low_high.plots, aes(x=w.con, y=R50, colour=year))+
  geom_point()+
  stat_smooth(method = "glm",
              formula = "y~x",
              geom = "smooth",
              col="Red")+
  theme(axis.text.x = element_text(angle = -90, vjust = 0.5, size = 17))+
  theme(axis.text.y = element_text(angle = 0, vjust = 0.5, size = 17))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  scale_x_continuous(labels = scales::number_format(accuracy = 0.01))+
  ylim(0, 0.5)+
  scale_color_gradient(low="yellow", high="Blue")+
  theme(plot.title = element_text(hjust = 0))+xlab("Weighted connectance")+ylab("R50")+
  theme(strip.text.x = element_text(size = 15))+
  theme(axis.title=element_text(size=22))+theme(legend.text=element_text(size=22),
                                                legend.title=element_text(size=22))+
  facet_wrap(~factor(threshold, levels=c("Links/no links","90%", "80%","70%","60%",
                                         "50%", "40%", "30%", "20%", "10%")),ncol=5)+
  theme(strip.background = element_blank())+theme(panel.border=element_rect(fill=NA))+
  theme(strip.text.x = element_text(size = 12))

#combined
si.fig5a+theme(legend.position = "none")+labs(tag = "a")+theme(plot.tag = element_text(size = 18,face = "bold"))+
  si.fig5b+theme(axis.title.y=element_blank())+theme(axis.text.y=element_blank())+labs(tag = "b")+
  theme(plot.tag = element_text(size = 18,face = "bold"))




#### SI fig 6 plot ####

si.fig6a<-ggplot(uw.DF.low_high.plots, aes(x=`A/C`, y=R50, colour=year))+
  geom_point()+
  stat_smooth(method = "glm",
              formula = "y~x",
              geom = "smooth",
              col="Red")+
  theme(axis.text.x = element_text(angle = -90, vjust = 0.5, size = 17))+
  theme(axis.text.y = element_text(angle = 0, vjust = 0.5, size = 17))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  scale_x_continuous(labels = scales::number_format(accuracy = 0.01))+
  ylim(0, 0.5)+
  scale_color_gradient(low="yellow", high="Blue")+
  theme(plot.title = element_text(hjust = 0))+xlab("Relative ascendency")+ylab("R50")+
  theme(strip.text.x = element_text(size = 15))+
  theme(axis.title=element_text(size=22))+theme(legend.text=element_text(size=22),
                                                legend.title=element_text(size=22))+
  facet_wrap(~factor(threshold, levels=c("Links/no links","90%", "80%","70%","60%",
                                         "50%", "40%", "30%", "20%", "10%")),ncol=5)+
  theme(strip.background = element_blank())+theme(panel.border=element_rect(fill=NA))+
  theme(strip.text.x = element_text(size = 12))

si.fig6b<-ggplot(w.DF.low_high.plots, aes(x=`A/C`, y=R50, colour=year))+
  geom_point()+
  stat_smooth(method = "glm",
              formula = "y~x",
              geom = "smooth",
              col="Red")+
  theme(axis.text.x = element_text(angle = -90, vjust = 0.5, size = 17))+
  theme(axis.text.y = element_text(angle = 0, vjust = 0.5, size = 17))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  scale_x_continuous(labels = scales::number_format(accuracy = 0.01))+
  ylim(0, 0.5)+
  scale_color_gradient(low="yellow", high="Blue")+
  theme(plot.title = element_text(hjust = 0))+xlab("Relative accendancy")+ylab("R50")+
  theme(strip.text.x = element_text(size = 15))+
  theme(axis.title=element_text(size=22))+theme(legend.text=element_text(size=22),
                                                legend.title=element_text(size=22))+
  facet_wrap(~factor(threshold, levels=c("Links/no links","90%", "80%","70%","60%",
                                         "50%", "40%", "30%", "20%", "10%")),ncol=5)+
  theme(strip.background = element_blank())+theme(panel.border=element_rect(fill=NA))+
  theme(strip.text.x = element_text(size = 12))

#combined
si.fig6a+theme(legend.position = "none")+labs(tag = "a")+theme(plot.tag = element_text(size = 18,face = "bold"))+
  si.fig6b+theme(axis.title.y=element_blank())+theme(axis.text.y=element_blank())+labs(tag = "b")+
  theme(plot.tag = element_text(size = 18,face = "bold"))
