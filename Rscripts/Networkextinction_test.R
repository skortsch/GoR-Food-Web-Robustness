### Triggering extinctions


#Patrick's wd 
#setwd("D:/R_Git/GoR-Food-webs/Rscripts")

#susanne's wd
setwd("C:/LocalData/susakort/GitHub/GoR-Food-webs/Rscripts")

w_fw<-readRDS("../Data/weighted_foodwebs.rds")

#load("D:/R_Git/GoR-Food-webs/Data/uw_fw_GOR.RData")
#load("D:/R_Git/GoR-Food-webs/Data/w_fw_GOR.RData")

#load GoR data
load("../Data/uw_fw_GOR.RData")
load("../Data/w_fw_GOR.RData")
load("../Data/metawebs_GoR.RData")
              
library(NetworkExtinction)
library(network)
library(igraph)

data("net")
SimulateExtinctions(Network = net, Method = "Mostconnected")#exemple, fungerar

#annat exempel network
data("chilean_intertidal")
chilean_intertidal
Mostconnected <- SimulateExtinctions(Network = chilean_intertidal,
                                     Method = "Mostconnected",
                                     NetworkType = "Trophic", # default argument
                                     IS = 0, # default argument
                                     Rewiring = FALSE # default argument
)

Mostconnected$sims
as.matrix(chilean_intertidal)

#w_mw
load("../Data/metawebs_GoR.RData")
w_mw<-as.matrix(w_mw)
rownames(w_mw)<-c(1:dim(w_mw)[1])
colnames(w_mw)<-c(1:dim(w_mw)[1])
w_mw<-as.network(w_mw)
as.network.matrix(w_mw,matrix.type="adjacency")


#unweighted metaweb GoR
un_mw<-as.matrix(un_mw)
rownames(un_mw)<-c(1:dim(un_mw)[1])
colnames(un_mw)<-c(1:dim(un_mw)[1])
un_mw<-as.network(un_mw)
un_mw<-as.network.matrix(un_mw,matrix.type="adjacency")
SimulateExtinctions(un_mw, NetworkType = "Trophic", Method = "Mostconnected")# fungerar inte -->crash, fatal error

Network<-un_mw
Network<-w_mw
#Network<-net
Network <- .DataInit(x = Network)
#if (!NetworkType %in% c("Trophic", "Mutualistic")) {
#  stop("Please specify NetworkType as either 'Trophic' or 'Mutualistic'")
edgelist <- network::as.matrix.network.edgelist(Network, matrix.type = "edgelist")
Grado <- NULL
Conected <- data.frame(ID = 1:network::network.size(Network), Grado = sna::degree(edgelist, c("total")))
Conected <- dplyr::arrange(Conected, desc(Grado))$ID
DF <- ExtinctionOrder(Network = Network, Order = Conected, NetworkType = "Trophic")
#weighted network
DF <- ExtinctionOrder(Network = Network, Order = Conected, NetworkType = "Trophic", IS = 0.5)


#nyt format?
test81<-as.network(fw.1981)
test81# är nu class: network, ser nu ut som "net" eller "Chilean_intertidal"

SimulateExtinctions(test81, Method = "Mostconnected")# fungerar inte -->crash, fatal error

#inga namn? VErkar vara ett av problemen
fw.1981
colnames(fw.1981)<-c(1:25)
rownames(fw.1981)<-c(1:25)
fw.1981

SimulateExtinctions(fw.1981, Method = "Mostconnected")# error in data.frame(ID = 1:network::network.size(Temp), Grado = sna::degree(edgelist,  : 
                                                      #arguments imply differing number of rows: 2, 0

#men det här går
ExtinctionOrder(fw.1981, Order = c(2,3,4,5,15), NetworkType = "Trophic")#fungerar, så går att använda om vi har en extinction order vi vill använda, tex top roles eller link sum
ExtinctionOrder(fw.1981, Order = c(1:25), NetworkType = "Trophic") #fungerar

RandomExtinctions(fw.1981, NetworkType = "Trophic")#fungerar


wfw.1981
colnames(wfw.1981)<-c(1:25)
rownames(wfw.1981)<-c(1:25)
wfw.1981

ExtinctionOrder(wfw.1981, Order = c(2,3,4,5,15), NetworkType = "Trophic") #namn är problemet, fungerar nu med weighted också
ExtinctionOrder(wfw.1981, Order = c(1:25), NetworkType = "Trophic")

#Fortfarande något problem???
SimulateExtinctions(wfw.1981, Method = "Mostconnected")# error in data.frame(ID = 1:network::network.size(Temp), Grado = sna::degree(edgelist,  : 
                                                      #arguments imply differing number of rows: 2, 0
