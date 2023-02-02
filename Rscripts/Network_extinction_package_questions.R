###Code for testing

library(NetworkExtinction)
library(network)
library(igraph)

load("../Data/metawebs_GoR.RData")

#UNWEIGHTED NETWORKS: question 1

#unweighted metaweb GoR
un_mw<-as.matrix(un_mw)
rownames(un_mw)<-c(1:dim(un_mw)[1])
colnames(un_mw)<-c(1:dim(un_mw)[1])
un_mw<-as.network(un_mw)
#un_mw<-as.network.matrix(un_mw,matrix.type="adjacency")

SimulateExtinctions(un_mw, NetworkType = "Trophic", Method = "Mostconnected")# fungerar inte -->crash, fatal error

#now run the code like this
Network<-un_mw
#Network<-net
#Network <- .DataInit(x = Network)
#if (!NetworkType %in% c("Trophic", "Mutualistic")) {
#  stop("Please specify NetworkType as either 'Trophic' or 'Mutualistic'")
edgelist <- network::as.matrix.network.edgelist(Network, matrix.type = "edgelist")
Grado <- NULL
Conected <- data.frame(ID = 1:network::network.size(Network), Grado = sna::degree(edgelist, c("total")))
Conected <- dplyr::arrange(Conected, desc(Grado))$ID
DF <- ExtinctionOrder(Network = Network, Order = Conected, NetworkType = "Trophic")


#Question1: why does the function not work? 
#I can make the code run when I use part of the function

###WEIGHTED NETWORKS: question 2

#w_mw
load("../Data/metawebs_GoR.RData")
w_mw<-as.matrix(w_mw)
rownames(w_mw)<-c(1:dim(w_mw)[1])
colnames(w_mw)<-c(1:dim(w_mw)[1])
w_mw<-as.network(w_mw)
#as.network.matrix(w_mw,matrix.type="adjacency")

#Question 2: where are the edge attributes? Our matrix contains edge values

Network<-w_mw
#Network<-net
#Network <- .DataInit(x = Network)
#if (!NetworkType %in% c("Trophic", "Mutualistic")) {
#  stop("Please specify NetworkType as either 'Trophic' or 'Mutualistic'")
edgelist <- network::as.matrix.network.edgelist(Network, matrix.type = "edgelist")
Grado <- NULL
Conected <- data.frame(ID = 1:network::network.size(Network), Grado = sna::degree(edgelist, c("total")))
Conected <- dplyr::arrange(Conected, desc(Grado))$ID
#weighted network
DF <- ExtinctionOrder(Network = Network, Order = Conected, NetworkType = "Trophic", IS = 0.5)
