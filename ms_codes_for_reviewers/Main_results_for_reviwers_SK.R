#Script for reviewers of the ms "Food web robustness depends on network type and threshold for extinction"

#This script contains code to replicate the main results from the manuscript.
#The script is divided into sections for the different steps.
#Some notes and annotations are included in each section.

#R-version was used: R 4.2.2, it has not been tested if it runs on other versions.

#For the code for the other deletion sequences presented in the SI of the manuscript
#see the other two scripts, "Low_to_high_sequence.R" and "random_sequence.R"

#### packages ####

install.packages(c("igraph", "ggplot2","MASS", "NetIndecies", "fluxweb", "NetworkExtinctions",
                   "network", "dplyr", "patchwork", "scales", "tidyverse", "fmsb"))
#load
library(ggplot2) #!!!! DO YOU NEED ggplot  when you have tidyverse???
library(igraph) #network calculations and visualization
library(MASS)
library(NetIndices) #calculation of information-theoretic network (AMI) indices
library(fluxweb) #bioenergetic modeling
library(NetworkExtinction) #secondary extinctions using thresholds
library(network)
library(dplyr) #!!!! DO YOU NEED dplyr when you have tidyverse???
library(patchwork) #plotting
library(scales)
library(tidyverse) #for data wrangling
library(fmsb)

#### Data ####
#Foodwebs used in the simulations is stored in the files "uw.webs.1000.RDS" and "w.webs.1000.RDS"

#read data
uw.webs.1000<-readRDS("uw.webs.1000.RDS")#the unweighted webs, format is list containing interaction matrices

w.webs.1000<-readRDS("w.webs.1000.RDS")#the weighted webs, format is list containing interaction matrices

#The lists have a length if 1000, each entry has 34 webs inside, one for each year in the time series. So 34000 food webs total
#uw.webs.1000[[1]][[1]] is the first of 1000 un-weighted web for the year 1980, uw.webs.1000[[1]][[2]] is the first un-weighted web for 1981 and so on.
#same for the weighted webs


######### NETWORK METRICS ############

####Connectance calculation in the unweighted networks####
#Connectance is The fraction of all possible realized links in a network.

Con<-matrix(NA, 34, 1000)#set up matrix to store the output
#loop to calculate connectance for all 34000 webs.
for (i in 1:1000){
  a<-uw.webs.1000[[i]]
  
  for (j in 1:34){
    fw<-a[[j]]
    fw_g<-graph.adjacency(fw, mode = "directed")
    Con[j,i]<-edge_density(fw_g)#connectance edge_density() is an igraph function to get the connectance of a web
    
  }
}
Con

####Custom code for function to calculate link weighted connectance####
# fluxind(), custom function, tutorial: https://rfrelat.github.io/BalticFoodWeb.html
# This custom function is used in a loop in the next section, run once before moving on so it is in the .GlobalEnv
fluxind<- function(fluxes){
  res<-c()
  # The flux matrix
  W.net<-as.matrix(fluxes) #fluxmatrix from fluxweb
  
  ### Taxon-specific Shannon indices of inflows
  # sum of k species inflows --> colsums
  sum.in<-apply(W.net, 2, sum)
  
  # Diversity of k species inflows
  H.in.mat<- t(t(W.net)/sum.in)*t(log(t(W.net)/sum.in)) #columns divided by the total col sum
  H.in.mat[!is.finite(H.in.mat)] <- 0 #converts NaN to 0's
  H.in<- apply(H.in.mat, 2, sum)*-1
  
  ### Taxon-specific Shannon indices of outflows
  # sum of k speies outflows --> rowsums
  sum.out<-apply(W.net, 1, sum)
  
  # Diversity of k species outflows
  H.out.mat<- (W.net/sum.out)*log(W.net/sum.out) #rows divided by the total row sum
  H.out.mat[!is.finite(H.out.mat)] <- 0 #converts NaN to 0's
  H.out<- apply(H.out.mat, 1, sum)*-1
  
  # Effective number of prey or resources = N(R,k) 
  # The reciprocal of H(R,k) --> N (R,k) is the equivalent number of prey for species k
  # N.res<-exp(H.in)
  N.res<-ifelse(sum.in==0, H.in, exp(H.in))
  
  # Effective number of predators or consumers = N(C,k) 
  # The reciprocal of H(C,k) --> N (C,k) is the equivalent number of predators for species k
  # N.con<-exp(H.out)
  N.con<-ifelse(sum.out==0, H.out, exp(H.out))
  
  ### Quantitative Weighted and unweighted link density and weighted connectance
  no.species<-ncol(W.net)
  # The unweighted link density (LDw) is:
  # LD.uw<- sum(N.res/no.species) + sum(N.con/no.species)/2
  #I think parenthesis are missing in above formula
  #I decided to follow the equation in the manuscript
  res$qLD.uw<- 1/(2*no.species)* (sum(N.res) + sum(N.con))
  
  #unweighted connectance
  res$qC.uw<-res$qLD.uw/no.species
  
  # The weighted link density (LDw) is:
  # In the weighted version the effective number of predators for species i is weighted by i's 
  # contribution to the total outflow
  # the same is the case for the inflows
  
  tot.mat<- sum(W.net)
  # LD.w <- (sum((sum.in/tot.mat)*N.res) + sum((sum.out/tot.mat)*N.con))/2
  # equivalent to next formula, but next one is closer to manuscript
  res$qLD.w <- 1/(2*tot.mat)*(sum(sum.in*N.res) + sum(sum.out*N.con))
  
  #Weighted connectance
  res$qC.w<- res$qLD.w/no.species
  
  #total number of qualitative links
  res$qL <- no.species * res$qLD.uw
  
  # positional.index
  pos.ind<- sum.in*N.res/(sum.in*N.res+sum.out*N.con) #postional index
  basal.sp<-pos.ind[pos.ind==0] #basal species = 0
  top.sp<- pos.ind[pos.ind==1] #top predator =1
  #defintion according to Bersier et al. 2002 top species = [0.99, 1]
  #int.sp<- pos.ind[pos.ind>0&pos.ind<0.99]# intermediate are all species that are not basal nor top
  
  con.sp<-length(pos.ind)-length(basal.sp)# all consumer taxa except basal
  # unweighted quantitative Generality
  res$qG.uw<- sum(N.res)/con.sp
  # weighted quantitative Generality
  res$qG.w<-sum(sum.in*N.res/sum(W.net))
  
  res.sp<- length(pos.ind)-length(top.sp)
  #unweighted quantitative Vulnerability
  res$qV.uw<- sum(N.con)/res.sp
  # weighted quantitative Vulnerability
  res$qV.w<-sum(sum.out*N.con/sum(W.net))
  
  #Standard deviation of stanardized unweighted quantitative Generality
  s<-length(pos.ind)
  res$qGsd.uw<-sd(s*N.res/sum(N.res)) #the mean  of (s*N.res/sum(N.res)  is = 1
  
  #Standard deviation of stanardized weighted quantitative Generality
  res$qGsd.w<-sd(s*sum.in*N.res/sum(N.res*sum.in))
  
  #Standard deviation of stanardized unweighted quantitative Vulnerability
  res$qVsd.uw<-sd(s*N.con/sum(N.con)) #the mean of  s*N.con/sum(N.con) is = 1
  
  #Standard deviation of stanardized weighted quantitative Vulnerability
  res$qVsd.w<-sd(s*sum.out*N.con/sum(N.con*sum.out))
  
  return(res)
}

####Weighted connectance list####
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


w.con<-matrix(NA, 34, 1000)#to store output in
for(i in 1:1000){
  #i<-1
  a<-w.metrics[[i]]
  
  for( j in 1:34){
    #j<-1
    b<-a[[j]]$qC.w
    w.con[j,i]<-b
  }
}

w.con


####Relative ascendency calculations####

#calculations with unweighted webs
uw.rel.asc<-matrix(NA,34,1000)#seting up output matrix
#loop used to get relative ascendency for all 34000 webs.
for (i in 1:1000) {
  a<-uw.webs.1000[[i]]
  for (j in 1:34) {
    b<-a[[j]]
    asc<-AscInd(b) # AscInd() is a function in NetIndices package used to calculate relative ascendency
    uw.rel.asc[j,i]<-asc[1,4]
  }
}
uw.rel.asc#output in matrix format, unweighted webs

#repeat with weighted webs
w.rel.asc<-matrix(NA,34,1000)
#loop used to get relative ascendency for all 34000 webs.
for (i in 1:1000) {
  a<-w.webs.1000[[i]]
  for (j in 1:34) {
    b<-a[[j]]
    asc<-AscInd(b) # AscInd() is a function in NetIndices package used to calculate relative ascendency
    w.rel.asc[j,i]<-asc[1,4]
  }
}
w.rel.asc #output in matrix format, weighted webs
######################################################################################################


########## FIGURE 2 ################

#Connectance
df.con<-matrix(NA,34,6)
colnames(df.con)<-c("Median","upperCI95","lowerCI95","upperCI50","lowerCI50","Year")
df.con[,6]<-c(1981:2014)
df.con[,1]<-apply(Con,1,median)
# calculate the CI
for (i in 1:34){
  df.con[i,2]<-quantile(Con[i,], probs = 0.975)
  df.con[i,3]<-quantile(Con[i,], probs = 0.05/2)
  df.con[i,4]<-quantile(Con[i,], probs = 0.75)
  df.con[i,5]<-quantile(Con[i,], probs = 0.25/2)
}
df.con<-as.data.frame(df.con)
g.con<-ggplot(df.con, aes(x=Year))+
  geom_line(aes(y=Median), col="red", linewidth=1.1)+
  geom_ribbon(aes(ymin=upperCI95, ymax=lowerCI95), fill="red",alpha=0.2)+
  geom_ribbon(aes(ymin=upperCI50, ymax=lowerCI50), fill="red",alpha=0.3)+
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, size = 18))+
  theme(axis.text.y = element_text(angle = 45, vjust = 0.5, size = 18))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  scale_x_continuous(name = "Year", breaks = scales::pretty_breaks(n = 5))+
  ylab("Connectance")+theme(axis.title=element_text(size=18))+
  scale_y_continuous(labels = label_number(accuracy = 0.01),
                     breaks = scales::pretty_breaks(n = 5))

#weighted connectance
df.wcon<-matrix(NA,34,6)
colnames(df.wcon)<-c("Median","upperCI95","lowerCI95","upperCI50","lowerCI50","Year")
df.wcon[,6]<-c(1981:2014)
df.wcon[,1]<-apply(w.con,1,median)
for (i in 1:34){
  df.wcon[i,2]<-quantile(w.con[i,], probs = 0.975)
  df.wcon[i,3]<-quantile(w.con[i,], probs = 0.05/2)
  df.wcon[i,4]<-quantile(w.con[i,], probs = 0.75)
  df.wcon[i,5]<-quantile(w.con[i,], probs = 0.25/2)
}
df.wcon<-as.data.frame(df.wcon)

g.wcon<-ggplot(df.wcon, aes(x=Year))+
  geom_line(aes(y=Median), col="orange", linewidth=1.1)+
  geom_ribbon(aes(ymin=upperCI95, ymax=lowerCI95), fill="orange",alpha=0.2)+
  geom_ribbon(aes(ymin=upperCI50, ymax=lowerCI50), fill="orange",alpha=0.3)+
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, size = 18))+
  theme(axis.text.y = element_text(angle = 45, vjust = 0.5, size = 18))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  scale_x_continuous(name = "Year", breaks = scales::pretty_breaks(n = 5))+
  scale_y_continuous(labels = label_number(accuracy = 0.01),
                     breaks = scales::pretty_breaks(n = 4))+
  xlab("Year")+ylab("Weighted Connectance")+theme(axis.title=element_text(size=18))

#Relatice ascendency in unweighted webs
df<-matrix(NA,34,6)
colnames(df)<-c("Year","Median","upCI95","lowCI95", "upCI50", "lowCI50")
df[,1]<-c(seq(1981,2014))
df[,2]<-apply(uw.rel.asc,1,median)
for (i in 1:34){
  df[i,3]<-quantile(uw.rel.asc[i,], probs = 0.975)
  df[i,4]<-quantile(uw.rel.asc[i,], probs = 0.05/2)
  df[i,5]<-quantile(uw.rel.asc[i,], probs = 0.75)
  df[i,6]<-quantile(uw.rel.asc[i,], probs = 0.25/2)
}
df<-as.data.frame(df)

g.uw.rel.asc<-ggplot(df, aes(x=Year))+
  geom_line(aes(y=Median), col="blue", linewidth=1.1)+
  geom_ribbon(aes(ymin=upCI95, ymax=lowCI95), fill="blue",alpha=0.2)+
  geom_ribbon(aes(ymin=upCI50, ymax=lowCI50), fill="blue",alpha=0.3)+
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, size = 18))+
  theme(axis.text.y = element_text(angle = 45, vjust = 0.5, size = 18))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  scale_x_continuous(name = "Year", breaks = scales::pretty_breaks(n = 5))+
  scale_y_continuous(breaks = scales::pretty_breaks(n = 5),
                     labels = label_number(accuracy = 0.01))+
  xlab("Year")+ylab("A/C unweighted webs")+theme(axis.title=element_text(size=18))

# relative ascendecny weoihgted webs
df<-matrix(NA,34,6)
colnames(df)<-c("Year","Median","upCI95","lowCI95", "upCI50", "lowCI50")
df[,1]<-c(seq(1981,2014))
df[,2]<-apply(w.rel.asc,1,median)

for (i in 1:34){
  df[i,3]<-quantile(w.rel.asc[i,], probs = 0.975)
  df[i,4]<-quantile(w.rel.asc[i,], probs = 0.05/2)
  df[i,5]<-quantile(w.rel.asc[i,], probs = 0.75)
  df[i,6]<-quantile(w.rel.asc[i,], probs = 0.25/2)
}

df<-as.data.frame(df)
g.rel.asc<-ggplot(df, aes(x=Year))+
  geom_line(aes(y=Median), col="cyan", linewidth=1.1)+
  geom_ribbon(aes(ymin=upCI95, ymax=lowCI95), fill="cyan",alpha=0.2)+
  geom_ribbon(aes(ymin=upCI50, ymax=lowCI50), fill="cyan",alpha=0.3)+
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, size = 18))+
  theme(axis.text.y = element_text(angle = 45, vjust = 0.5, size = 18))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  scale_x_continuous(name = "Year", breaks = scales::pretty_breaks(n = 5))+
  scale_y_continuous(breaks = scales::pretty_breaks(n = 5),
                     labels = label_number(accuracy = 0.01))+
  xlab("Year")+ylab("A/C weighted webs")+theme(axis.title=element_text(size=18))

#all in one plot
fig.2<-g.con+labs(tag = "a")+theme(plot.tag = element_text(hjust = 1, size = 18,face = "bold"))+
  theme(axis.text.x = element_text(angle = -90, vjust = 0.5, size = 18))+
  theme(axis.text.y = element_text(angle = 0, vjust = 0.5, size = 18))+
  g.wcon+labs(tag = "b")+theme(plot.tag = element_text(hjust = 0.5, size = 18,face = "bold"))+
  theme(axis.text.x = element_text(angle = -90, vjust = 0.5, size = 18))+
  theme(axis.text.y = element_text(angle = 0, vjust = 0.5, size = 18))+
  g.uw.rel.asc+labs(tag = "c")+theme(plot.tag = element_text(hjust = 0.5, size = 18,face = "bold"))+
  theme(axis.text.x = element_text(angle = -90, vjust = 0.5, size = 18))+
  theme(axis.text.y = element_text(angle = 0, vjust = 0.5, size = 18))+
  g.rel.asc+labs(tag = "d")+theme(plot.tag = element_text(hjust = 0.5, size = 18,face = "bold"))+ theme(axis.text.x = element_text(angle = -90, vjust = 0.5, size = 18))+
  theme(axis.text.y = element_text(angle = 0, vjust = 0.5, size = 18))+
  plot_layout(ncol = 2)

#print figure
fig.2



####Species deletion simulation and R50 calculation####

#set up matrixes to store outputs
uw.R50<-matrix(NA, 34, 100)
uw.R50.th10<-matrix(NA, 34, 100)
uw.R50.th20<-matrix(NA, 34, 100)
uw.R50.th30<-matrix(NA, 34, 100)
uw.R50.th40<-matrix(NA, 34, 100)
uw.R50.th50<-matrix(NA, 34, 100)
uw.R50.th60<-matrix(NA, 34, 100)
uw.R50.th70<-matrix(NA, 34, 100)
uw.R50.th80<-matrix(NA, 34, 100)
uw.R50.th90<-matrix(NA, 34, 100)

#The following loop goes through all 34000 webs, and simulating deletions at different thresholds
#NOTE: Has a very long run-time. Recommended to leave running over night.
#See file "uw.R50results.RDS" for R50 results used in the manuscript to compare output.

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

webs<-uw.webs.100
for (i in 1:100){
  a<-webs[[i]] 
  
  for (j in 1:34) {
    fw<-a[[j]]
    
    fw_i<-fw
    colnames(fw_i)<-c(1:dim(fw)[1])# changes node names to numbers so extinction function can run, could not run with species names
    rownames(fw_i)<-c(1:dim(fw)[1])# changes node names to numbers so extinction function can run, could not run with species names
    
    deg.sp<-apply(fw_i, 1, sum)+apply(fw_i, 2, sum)#summing species in and out going links to get degree of the species
    order_high_low<-as.numeric(names(sort(deg.sp, decreasing = T)))#extion order set to from high to low degree
    
    #remove basal nodes from extinction order: autotrophs, mixotrophs, detritus
    order<-setdiff(order_high_low, c(which(colnames(fw)=="Autotroph"),
                                              which(colnames(fw)=="Mixotroph"),
                                              which(colnames(fw)=="Detritus")))
    
   
    #ExtinctionOrder() function from NetworkExtinction package. Simulates extinctions in order given and at set thresholds.
    
    #no threshold
    ext<-ExtinctionOrder(fw_i, Order = order, NetworkType = "Trophic",
                                  verbose = F)
    uw.R50[j,i]<-ext$R50#extracting R50 from output
    #90% threshold, IS=0.1 means that a species is removed when 10% of its energy intake is left
    ext<-ExtinctionOrder(fw_i, Order = order, NetworkType = "Trophic",
                                  verbose = F, IS=0.1)
    uw.R50.th90[j,i]<-ext$R50 
    #80% threshold, IS=0.2 means that a species is removed when 20% of its energy intake is left
    ext<-ExtinctionOrder(fw_i, Order = order, NetworkType = "Trophic",
                                  verbose = F, IS=0.2)
    uw.R50.th80[j,i]<-ext$R50
    #70% threshold, IS=0.3 means that a species is removed when 30% of its energy intake is left
    ext<-ExtinctionOrder(fw_i, Order = order, NetworkType = "Trophic",
                                  verbose = F, IS=0.3)
    uw.R50.th70[j,i]<-ext$R50
    #60% threshold, IS=0.4 means that a species is removed when 40% of its energy intake is left
    ext<-ExtinctionOrder(fw_i, Order = order, NetworkType = "Trophic",
                                  verbose = F, IS=0.4)
    uw.R50.th60[j,i]<-ext$R50 
    #50% threshold, IS=0.5 means that a species is removed when 50% of its energy intake is left
    ext<-ExtinctionOrder(fw_i, Order = order, NetworkType = "Trophic",
                                  verbose = F, IS=0.5)
    #threshold, IS=0.6 means that a species is removed when 40% of its energy intake is left
    uw.R50.th50[j,i]<-ext$R50 
    ext<-ExtinctionOrder(fw_i, Order = order, NetworkType = "Trophic",
                                  verbose = F, IS=0.6)
    #threshold, IS=0.7 means that a species is removed when 60% of its energy intake is left
    uw.R50.th40[j,i]<-ext$R50 
    ext<-ExtinctionOrder(fw_i, Order = order, NetworkType = "Trophic",
                                  verbose = F, IS=0.7)
    uw.R50.th30[j,i]<-ext$R50 
    
    #80% threshold, IS=0.8 means that a species is removed when 20% of its energy intake is left
    ext<-ExtinctionOrder(fw_i, Order = order, NetworkType = "Trophic",
                                  verbose = F, IS=0.8)
    uw.R50.th20[j,i]<-ext$R50
    #10%, IS=0.9 means a species is removed when 90% of of its energy intake is left
    ext<-ExtinctionOrder(fw_i, Order = order, NetworkType = "Trophic",
                                  verbose = T, IS=0.9) #verbose=T here so you can se that the loop is running, it just prints progress bar
    uw.R50.th10[j,i]<-ext$R50 
  }
}
#output with all the R50 values. 
#See section "Data frames" for code to make a easy to use data frame from the output
uwR50.list<-list(uw.R50, uw.R50.th90, uw.R50.th80, uw.R50.th70, uw.R50.th60,
                      uw.R50.th50,uw.R50.th40, uw.R50.th30, uw.R50.th20, uw.R50.th10)




#### weighted webs R50 simulations ####
#The following loop goes through all 34000 webs, and simulating deletions at different thresholds
#NOTE: Has a very long run-time. Recommended to leave running over night.

#IMPORTANT NOTE: Will give several warnings at the end! This does not mean there was an error.
#For example, the following Warning messages show up:

#"1: In network::as.matrix.network.adjacency(Temp, attrname = "weight"):
#There is no edge attribute named weight"
#NOTE:This one comes up when the webs are disassembled to the point that 
#only the basal nodes are left, the basal nodes are not connected, therefore there are no links that have link weights

#"2: In ExtinctionOrder(fw_i, Order = order, NetworkType = "Trophic",  :
#Your network became completely unconnected before all primary extinctions were simulated. This happened at extinction step X out of Y"
#NOTE: Self explanatory. This happens when extinctions trigger cascades and disassembled the network completely.

#"3: In ExtinctionOrder(fw_i, Order = order, NetworkType = "Trophic",  :
#Primary extinctions of X, Y skipped due to their prior extinction as secondary extinctions."
#NOTE: Self explanatory. The extinction already happend, so was skipped in the order.

#See file "w.R50list.RDS" for R50 results used in the manuscript for comparison.
w.R50<-matrix(NA, 34, 1000)
w.R50.th10<-matrix(NA, 34, 1000)
w.R50.th20<-matrix(NA, 34, 1000)
w.R50.th30<-matrix(NA, 34, 1000)
w.R50.th40<-matrix(NA, 34, 1000)
w.R50.th50<-matrix(NA, 34, 1000)
w.R50.th60<-matrix(NA, 34, 1000)
w.R50.th70<-matrix(NA, 34, 1000)
w.R50.th80<-matrix(NA, 34, 1000)
w.R50.th90<-matrix(NA, 34, 1000)

webs<-w.webs.1000
for (i in 1:1000){
  a<-webs[[i]]
  
  for (j in 1:34) {
    fw<-a[[i]]
    
    fw_i<-fw
    colnames(fw_i)<-c(1:dim(fw)[1])# changes node names to numbers so extinction function can run, cold not run with species names
    rownames(fw_i)<-c(1:dim(fw)[1])# changes node names to numbers so extinction function can run, cold not run with species names
    
    sum_of_links.sp<-apply(fw_i, 1, sum)+apply(fw_i, 2, sum)#summing species in and out going links to get sum of link weights of the species
    order_high_low<-as.numeric(names(sort(sum_of_links.sp, decreasing = T)))#extincion order set to from high to low sum of link weights
    
    #remove basal nodes from extinction order: autotrophs, mixotrophs, detrits
    order<-setdiff(order_high_low, c(which(colnames(fw)=="Autotroph"),
                                     which(colnames(fw)=="Mixotroph"),
                                     which(colnames(fw)=="Detrits")))
    
    
    #ExtinctionOrder() function from NetworkExtinction package. Simulates extinctions in order given and at set thresholds.
    
    #no threshold
    ext<-ExtinctionOrder(fw_i, Order = order, NetworkType = "Trophic",
                         verbose = F)#NetworkType=Trophic, means secondary extinctions only occur for any consumer, but not producers. 
    w.R50[j,i]<-ext$R50#extracting R50 from outpt
    #90% threshold, IS=0.1 means that a species is removed when 10% of its energy intake is left
    ext<-ExtinctionOrder(fw_i, Order = order, NetworkType = "Trophic",
                         verbose = F, IS=0.1)
    w.R50.th90[j,i]<-ext$R50 
    #80% threshold, IS=0.2 means that a species is removed when 20% of its energy intake is left
    ext<-ExtinctionOrder(fw_i, Order = order, NetworkType = "Trophic",
                         verbose = F, IS=0.2)
    w.R50.th80[j,i]<-ext$R50
    #70% threshold, IS=0.3 means that a species is removed when 30% of its energy intake is left
    ext<-ExtinctionOrder(fw_i, Order = order, NetworkType = "Trophic",
                         verbose = F, IS=0.3)
    w.R50.th70[j,i]<-ext$R50
    #60% threshold, IS=0.4 means that a species is removed when 40% of its energy intake is left
    ext<-ExtinctionOrder(fw_i, Order = order, NetworkType = "Trophic",
                         verbose = F, IS=0.4)
    w.R50.th60[j,i]<-ext$R50 
    #50% threshold, IS=0.5 means that a species is removed when 50% of its energy intake is left
    ext<-ExtinctionOrder(fw_i, Order = order, NetworkType = "Trophic",
                         verbose = F, IS=0.5)
    #threshold, IS=0.6 means that a species is removed when 40% of its energy intake is left
    w.R50.th50[j,i]<-ext$R50 
    ext<-ExtinctionOrder(fw_i, Order = order, NetworkType = "Trophic",
                         verbose = F, IS=0.6)
    #threshold, IS=0.7 means that a species is removed when 60% of its energy intake is left
    w.R50.th40[j,i]<-ext$R50 
    ext<-ExtinctionOrder(fw_i, Order = order, NetworkType = "Trophic",
                         verbose = F, IS=0.7)
    w.R50.th30[j,i]<-ext$R50 
    #80% threshold, IS=0.8 means that a species is removed when 20% of its energy intake is left
    ext<-ExtinctionOrder(fw_i, Order = order, NetworkType = "Trophic",
                         verbose = F, IS=0.8)
    w.R50.th20[j,i]<-ext$R50
    #10%, IS=0.9 means a species is removed when 90% of of its energy intake is left
    ext<-ExtinctionOrder(fw_i, Order = order, NetworkType = "Trophic",
                         verbose = T, IS=0.9) #verbose=T here so yo can see that the loop is running, it just prints progress bar
    w.R50.th10[j,i]<-ext$R50 
  }
}
#output with all the R50 values. 
#See section "Data frames" for code to make a easy to use data frame from the output
wR50.list<-list(w.R50, w.R50.th90, w.R50.th80, w.R50.th70, w.R50.th60,
                 w.R50.th50,w.R50.th40, w.R50.th30, w.R50.th20, w.R50.th10)

#compare to "w.R50results.RDS"



######## FIGURE 3 #######

#Code to rearrange the unweighted and weighted R50 outputs in a dataframe to plot fig.2 in the manuscript
#if you ran the whole extinction simulation for the unweighted and weighted webs, run the code as is.

#if you did not run the whole extinction simulation for the unweighted webs, first run "uwR50.list<-readRDS("uw.R50results.RDS")"
#if you did not run the whole extinction simulation for the weighted webs, first run "wR50.list<-readRDS("w.R50results.RDS")"
df.allTH<-matrix(NA,34,12)
df.allTH[,11]<-c(seq(1981,2014))
df.allTH[,12]<-c("a")
colnames(df.allTH)<-c("Links/no links","90%","80%", "70%", "60%",
                      "50%","40%", "30%", "20%", "10%", "Year", "webs")
for (i in 1:10){
  df.allTH[,i]<-apply(uwR50.list[[i]],1, median, na.rm=T)
}
as.data.frame(df.allTH)
d.uw<-pivot_longer(as.data.frame(df.allTH), cols = 1:10)

df.allTH<-matrix(NA,34,12)
df.allTH[,11]<-c(seq(1981,2014))
df.allTH[,12]<-c("b")
colnames(df.allTH)<-c("Links/no links","90%","80%", "70%", "60%",
                      "50%","40%", "30%", "20%", "10%", "Year", "webs")
for (i in 1:10){
  df.allTH[,i]<-apply(wR50.list[[i]],1, median, na.rm=T)
}
as.data.frame(df.allTH)
d.w<-pivot_longer(as.data.frame(df.allTH), cols = 1:10)
#combineinto one 
d.comb<-full_join(d.uw, d.w)
d.comb$Year<-as.numeric(as.character(d.comb$Year))
d.comb$value<-as.numeric(as.character(d.comb$value))
d.comb
d.comb$name<-factor(d.comb$name, levels = c("Links/no links","90%","80%", "70%", "60%",
                                            "50%","40%", "30%", "20%", "10%"))

fig.3<-ggplot(d.comb, aes(x=Year, y=value, color = name))+
  geom_line(aes(group=name),linewidth=1.1)+
  theme(axis.text.x = element_text(angle = -90, vjust = 0.5, size = 25))+
  theme(axis.text.y = element_text(angle = 0,vjust = 0.5, size = 25))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  xlab("Year")+ylab("R50")+theme(axis.title=element_text(size=28), axis.title.x = element_text(angle = 0),
                                 axis.title.y.left = element_text(angle = 90))+
  scale_y_continuous(breaks = scales::pretty_breaks(n = 5))+
  theme(axis.title.y = element_text(angle = 0, vjust = 0.5))+
  scale_color_brewer(palette="Spectral")+
  facet_wrap(~factor(webs),ncol=2)+labs(color = "Thresholds")+
  theme(strip.background = element_blank())+theme(panel.border=element_rect(fill=NA))+
  theme(strip.text.x = element_text(size = 18,face="bold", hjust = -0.01))+
  theme(legend.text=element_text(size=25))+theme(legend.title=element_text(size=25))+
  coord_cartesian(clip = "off")
fig.3



####Data frames####
#Making data-frame for ease of use in plotting and analysis#
#the following is a series of loops that produce and join data frames into one data frame for use in plots and analysis

#if you ran the whole extinction simulation for the unweighted webs, run the code as is.
#if you did not run the whole extinction simulation for the weighted webs, first run "uwR50.list<-readRDS("uw.R50results.RDS")"
#see file named "uw_results_df.RDS" for comparison or if you want to skip this part 
m<-matrix(NA,1000,34)
colnames(m)<-c(1981:2014)
row.names(m)<-c(1:1000)
mlist<-list()
for (j in 1:10){
  a<-uwR50.list[[j]]
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
  m3[,i]<-Con[i,]
}

m4<-as.data.frame(as.table(m3))
d2<-as_tibble(m4)
d2$Var2<-as.numeric(as.character(d2$Var2))
colnames(d2)<-c("iter", "year", "con")

u<-full_join(d2, df.R50list[[1]])
for (i in 1:10) {
  u<-left_join(u, df.R50list[[i]])
}
u

m3<-matrix(NA,1000,34)
colnames(m3)<-c(1981:2014)
row.names(m3)<-c(1:1000)
for (i in 1:34) {
  m3[,i]<-uw.rel.asc[i,]
}
m3
m4<-as.data.frame(as.table(m3))
AC<-as_tibble(m4)
AC$Var2<-as.numeric(as.character(AC$Var2))
colnames(AC)<-c("iter", "year", "A/C")

uw.DF<-right_join(u, AC)

uw.DF.plots<-pivot_longer(uw.DF, 4:13)# make longformat df for use in ggplot
colnames(uw.DF.plots)[5]<-c("threshold")
colnames(uw.DF.plots)[6]<-c("R50")

uw.DF#dataframe used for the correlation analysis for the unweighted webs
uw.DF.plots# dataframe used for the plots for the unweighted webs

#same for the weighted version
#if you ran the whole extinction simulation for the weighted webs run the code as is.
#if you did not run the whole extinction simulation for the weighted webs, first run "wR50.list<-readRDS("w.R50results.RDS")"
#see file named "w_results_df.RDS" for comparison or if you want to skip this part 
m<-matrix(NA,1000,34)
colnames(m)<-c(1981:2014)
row.names(m)<-c(1:1000)
mlist<-list()
for (j in 1:10){
  a<-wR50.list[[j]]
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

w.DF<-right_join(u, AC)

w.DF.plots<-pivot_longer(w.DF, 4:13)# make longformat df for use in ggplot
colnames(w.DF.plots)[5]<-c("threshold")
colnames(w.DF.plots)[6]<-c("R50")

w.DF# dataframe used for the correlation analysis
w.DF.plots#longformat dataframe used for ggplot


####Correlation of R50 and the metrics####
##qqplot to see if R50 data is normally distributed
ggplot(uw.DF.plots, aes(sample=R50))+
  stat_qq()+
  stat_qq_line()+
  facet_wrap(~factor(threshold,levels=c("Links/no links","10%", "20%","30%","40%",
                                        "50%", "60%", "70%", "80%", "90%")),ncol=5)+
  theme(strip.background = element_blank())+theme(panel.border=element_rect(fill=NA))+
  ggtitle(element_text("unweighted webs R50 qqplot at diffrent THs"))

#R50 and connectance
#Using glm with link=identity for the response variable with normal distribution (or close to)
#and glm with link= log  for the more skewed ones.

# the no threshold scenario 
mod<-glm(uw.DF$`Links/no links`~uw.DF$con, family = gaussian(link = "log")) 
summary(mod)
NagelkerkeR2(mod)

# the 90% threshold scenario
mod<-glm(uw.DF$`90%`~uw.DF$con, family = gaussian(link = "identity")) 
summary(mod)
NagelkerkeR2(mod)

# the 80% threshold scenario
mod<-glm(uw.DF$`80%`~uw.DF$con, family = gaussian(link = "identity")) 
summary(mod)
NagelkerkeR2(mod)

# the 70% threshold scenario
mod<-glm(uw.DF$`70%`~uw.DF$con, family = gaussian(link = "identity")) 
summary(mod)
NagelkerkeR2(mod)

# the 60% threshold scenario
mod<-glm(uw.DF$`60%`~uw.DF$con, family = gaussian(link = "log"))
summary(mod)
NagelkerkeR2(mod)

# the 50% threshold scenario
mod<-glm(uw.DF$`50%`~uw.DF$con, family = gaussian(link = "log")) 
summary(mod)
NagelkerkeR2(mod)

# the 40% threshold scenario
mod<-glm(uw.DF$`40%`~uw.DF$con, family = gaussian(link = "log"))
summary(mod)
NagelkerkeR2(mod)

# the 30% threshold scenario
mod<-glm(uw.DF$`30%`~uw.DF$con, family = gaussian(link = "log")) 
summary(mod)
NagelkerkeR2(mod)

# the 20% threshold scenario
mod<-glm(uw.DF$`20%`~uw.DF$con, family = gaussian(link = "log")) 
summary(mod)
NagelkerkeR2(mod)

# the 10% threshold scenario
mod<-glm(uw.DF$`10%`~uw.DF$con, family = gaussian(link = "log")) 
summary(mod)
NagelkerkeR2(mod)



#R50 and relative ascendency in unweighted webs
mod<-glm(uw.DF$`Links/no links`~uw.DF$`A/C`, family = gaussian(link = "log"))
summary(mod)
NagelkerkeR2(mod)
mod<-glm(uw.DF$`90%`~uw.DF$`A/C`, family = gaussian(link = "identity"))
summary(mod)
NagelkerkeR2(mod)
mod<-glm(uw.DF$`80%`~uw.DF$`A/C`, family = gaussian(link = "identity"))
summary(mod)
NagelkerkeR2(mod)
mod<-glm(uw.DF$`70%`~uw.DF$`A/C`, family = gaussian(link = "identity"))
summary(mod)
NagelkerkeR2(mod)
mod<-glm(uw.DF$`60%`~uw.DF$`A/C`, family = gaussian(link = "log"))
summary(mod)
NagelkerkeR2(mod)
mod<-glm(uw.DF$`50%`~uw.DF$`A/C`, family = gaussian(link = "log"))
summary(mod)
NagelkerkeR2(mod)
mod<-glm(uw.DF$`40%`~uw.DF$`A/C`, family = gaussian(link = "log"))
summary(mod)
NagelkerkeR2(mod)
mod<-glm(uw.DF$`30%`~uw.DF$`A/C`, family = gaussian(link = "log"))
summary(mod)
NagelkerkeR2(mod)
mod<-glm(uw.DF$`20%`~uw.DF$`A/C`, family = gaussian(link = "log"))
summary(mod)
NagelkerkeR2(mod)
mod<-glm(uw.DF$`10%`~uw.DF$`A/C`, family = gaussian(link = "log"))
summary(mod)
NagelkerkeR2(mod)


#Weighted webs
##qqplot to see if R50 data is normally distributed
ggplot(w.DF.plots, aes(sample=R50))+
  stat_qq()+
  stat_qq_line()+
  facet_wrap(~factor(threshold,levels=c("Links/no links","10%", "20%","30%","40%",
                                        "50%", "60%", "70%", "80%", "90%")),ncol=5)+
  theme(strip.background = element_blank())+theme(panel.border=element_rect(fill=NA))+
  ggtitle(element_text("weighted webs R50 qqplot at diffrent THs"))

#R50 and weighted connectance
#Using glm with link=identity for the response variable with normal distribution (or close to)
#and glm with link= log  for the more skewed ones.
mod<-glm(w.DF$`Links/no links`~w.DF$w.con, family = gaussian(link = "log"))
summary(mod)
NagelkerkeR2(mod)
mod<-glm(w.DF$`90%`~w.DF$w.con, family = gaussian(link = "log"))
summary(mod)
NagelkerkeR2(mod)
mod<-glm(w.DF$`80%`~w.DF$w.con, family = gaussian(link = "log"))
summary(mod)
NagelkerkeR2(mod)
mod<-glm(w.DF$`70%`~w.DF$w.con, family = gaussian(link = "log"))
summary(mod)
NagelkerkeR2(mod)
mod<-glm(w.DF$`60%`~w.DF$w.con, family = gaussian(link = "log"))
summary(mod)
NagelkerkeR2(mod)
mod<-glm(w.DF$`50%`~w.DF$w.con, family = gaussian(link = "log"))
summary(mod)
NagelkerkeR2(mod)
mod<-glm(w.DF$`40%`~w.DF$w.con, family = gaussian(link = "log"))
summary(mod)
NagelkerkeR2(mod)
mod<-glm(w.DF$`30%`~w.DF$w.con, family = gaussian(link = "log"))
summary(mod)
NagelkerkeR2(mod)
mod<-glm(w.DF$`20%`~w.DF$w.con, family = gaussian(link = "log"))
summary(mod)
NagelkerkeR2(mod)
mod<-glm(w.DF$`10%`~w.DF$w.con, family = gaussian(link = "log"))
summary(mod)
NagelkerkeR2(mod)

#R50 and relative ascendency in weighted webs
#Using glm with link=identity for the response variable with normal distribution (or close to)
#and glm with link= log  for the more skewed ones.
mod<-glm(w.DF$`Links/no links`~w.DF$`A/C`, family = gaussian(link = "log"))
summary(mod)
NagelkerkeR2(mod)
mod<-glm(w.DF$`90%`~w.DF$`A/C`, family = gaussian(link = "log"))
summary(mod)
NagelkerkeR2(mod)
mod<-glm(w.DF$`80%`~w.DF$`A/C`, family = gaussian(link = "log"))
summary(mod)
NagelkerkeR2(mod)
mod<-glm(w.DF$`70%`~w.DF$`A/C`, family = gaussian(link = "log"))
summary(mod)
NagelkerkeR2(mod)
mod<-glm(w.DF$`60%`~w.DF$`A/C`, family = gaussian(link = "log"))
summary(mod)
NagelkerkeR2(mod)
mod<-glm(w.DF$`50%`~w.DF$`A/C`, family = gaussian(link = "log"))
summary(mod)
NagelkerkeR2(mod)
mod<-glm(w.DF$`40%`~w.DF$`A/C`, family = gaussian(link = "log"))
summary(mod)
NagelkerkeR2(mod)
mod<-glm(w.DF$`30%`~w.DF$`A/C`, family = gaussian(link = "log"))
summary(mod)
NagelkerkeR2(mod)
mod<-glm(w.DF$`20%`~w.DF$`A/C`, family = gaussian(link = "log"))
summary(mod)
NagelkerkeR2(mod)
mod<-glm(w.DF$`10%`~w.DF$`A/C`, family = gaussian(link = "log"))
summary(mod)
NagelkerkeR2(mod)


########### FIGURE 4 & 5 ################
#Connectance
fig4a<-ggplot(uw.DF.plots, aes(x=con, y=R50, colour=year))+
  geom_point()+
  stat_smooth(method = "glm",
              formula = "y~x",
              geom = "smooth",
              col="Red")+
  theme(axis.text.x = element_text(angle = -90, vjust = 0.5, size = 15))+
  theme(axis.text.y = element_text(angle = 0, vjust = 0.5, size = 15))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  scale_x_continuous(labels = scales::number_format(accuracy = 0.01))+
  scale_color_gradient(low="yellow", high="Blue",
                       limits = c(1980, 2014), 
                       breaks = seq(1980, 2014, by = 10))+
  theme(plot.title = element_text(hjust = 0))+xlab("Connectance")+ylab("R50")+
  theme(strip.text.x = element_text(size = 15))+
  theme(axis.title=element_text(size=22))+theme(legend.text=element_text(size=22),
                                                legend.title=element_text(size=22))+
  facet_wrap(~factor(threshold, levels=c("Links/no links","90%", "80%","70%","60%",
                                         "50%", "40%", "30%", "20%", "10%")),ncol=5)+
  theme(strip.background = element_blank())+theme(panel.border=element_rect(fill=NA))+
  theme(strip.text.x = element_text(size = 12))

#weighted connectance 
fig4b<-ggplot(w.DF.plots, aes(x=w.con, y=R50, colour=year))+
  geom_point()+
  stat_smooth(method = "glm",
              formula = "y~x",
              geom = "smooth",
              col="Red")+
  theme(axis.text.x = element_text(angle = -90, vjust = 0.5, size = 15))+
  theme(axis.text.y = element_text(angle = 0, vjust = 0.5, size = 15))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  scale_x_continuous(labels = scales::number_format(accuracy = 0.01))+
  scale_color_gradient(low="yellow", high="Blue",
                       limits = c(1980, 2014), 
                       breaks = seq(1980, 2014, by = 10))+
  theme(plot.title = element_text(hjust = 0))+xlab("Weighted Connectance")+ylab("R50")+
  theme(strip.text.x = element_text(size = 15))+
  theme(axis.title=element_text(size=22))+theme(legend.text=element_text(size=22),
                                                legend.title=element_text(size=22))+
  facet_wrap(~factor(threshold, levels=c("Links/no links","90%", "80%","70%","60%",
                                         "50%", "40%", "30%", "20%", "10%")),ncol=5)+
  theme(strip.background = element_blank())+theme(panel.border=element_rect(fill=NA))+
  theme(strip.text.x = element_text(size = 12))

#combine into fig.4
fig4a+theme(legend.position = "none")+labs(tag = "a")+theme(plot.tag = element_text(size = 20,face = "bold"))+ylim(0,0.5)+
  fig4b+theme(axis.title.y=element_blank())+theme(axis.text.y=element_blank())+labs(tag = "b")+
  theme(plot.tag = element_text(size = 20,face = "bold"))+ylim(0,0.5)
  


#relative ascendency unweighted webs
fig5a<-ggplot(uw.DF.plots, aes(x=`A/C`, y=R50, colour=year))+
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
  scale_color_gradient(low="yellow", high="Blue",
                       limits = c(1980, 2014), 
                       breaks = seq(1980, 2014, by = 10))+
  theme(plot.title = element_text(hjust = 0))+xlab("Relative ascendency")+ylab("R50")+
  theme(strip.text.x = element_text(size = 15))+
  theme(axis.title=element_text(size=22))+theme(legend.text=element_text(size=22),
                                                legend.title=element_text(size=22))+
  facet_wrap(~factor(threshold, levels=c("Links/no links","90%", "80%","70%","60%",
                                         "50%", "40%", "30%", "20%", "10%")),ncol=5)+
  theme(strip.background = element_blank())+theme(panel.border=element_rect(fill=NA))+
  theme(strip.text.x = element_text(size = 12))


#relative ascendency weighted webs
fig5b<-ggplot(w.DF.plots, aes(x=`A/C`, y=R50, colour=year))+
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
  scale_color_gradient(low="yellow", high="Blue",
                       limits = c(1980, 2014), 
                       breaks = seq(1980, 2014, by = 10))+
  theme(plot.title = element_text(hjust = 0))+xlab("Relative ascendency")+ylab("R50")+
  theme(strip.text.x = element_text(size = 15))+
  theme(axis.title=element_text(size=22))+theme(legend.text=element_text(size=22),
                                                legend.title=element_text(size=22))+
  facet_wrap(~factor(threshold, levels=c("Links/no links","90%", "80%","70%","60%",
                                         "50%", "40%", "30%", "20%", "10%")),ncol=5)+
  theme(strip.background = element_blank())+theme(panel.border=element_rect(fill=NA))+
  theme(strip.text.x = element_text(size = 12))

#Combined into fig.5
fig5a+theme(legend.position = "none")+labs(tag = "a")+theme(plot.tag = element_text(size = 18,face = "bold"))+ylim(0,0.5)+
  fig5b+theme(axis.title.y=element_blank())+theme(axis.text.y=element_blank())+labs(tag = "b")+
  theme(plot.tag = element_text(size = 18,face = "bold"))+ylim(0,0.5)



#### SI fig 11 & 12 ####

# SI fig 11
df.siplot<-matrix(NA,34000,2)
df.siplot[,1]<-as.vector(w.rel.asc)
df.siplot[,2]<-as.vector(w.con)

df.si<-as.data.frame(df.siplot)

si.plot<-ggplot(data=df.si, aes(x=V2, y=V1))+
  geom_point(shape=1, alpha=1)+
  stat_smooth(method = "glm",
              formula = "y~x",
              geom = "smooth",
              col="Red")+
  scale_x_continuous(labels = scales::number_format(accuracy = 0.01))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  xlab("Weighted connectance")+ylab("A/C weighted webs")+
  theme(axis.text.x = element_text(size = 17))+
  theme(axis.text.y = element_text(size = 17))+
  theme(axis.title=element_text(size=22))

summary(glm(as.vector(w.rel.asc)~as.vector(w.con)))
NagelkerkeR2(glm(as.vector(w.rel.asc)~as.vector(w.con)))


# SI fig 12 
df.siplot<-matrix(NA,34000,2)
df.siplot[,1]<-as.vector(uw.rel.asc)
df.siplot[,2]<-as.vector(Con)

df.si<-as.data.frame(df.siplot)

si.plot2<-ggplot(data=df.si, aes(x=V2, y=V1))+
  geom_point(shape=1, alpha=1)+
  stat_smooth(method = "glm",
              formula = "y~x",
              geom = "smooth",
              col="Red")+
  scale_x_continuous(labels = scales::number_format(accuracy = 0.01))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  xlab("Connectance")+ylab("A/C unweighted webs")+
  theme(axis.text.x = element_text(size = 17))+
  theme(axis.text.y = element_text(size = 17))+
  theme(axis.title=element_text(size=22))

summary(glm(as.vector(uw.rel.asc)~as.vector(Con)))
NagelkerkeR2(glm(as.vector(uw.rel.asc)~as.vector(Con)))
