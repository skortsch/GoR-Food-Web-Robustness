setwd("D:/R_Git/GoR-Food-webs/Rscripts")

w_fw<-readRDS("../Data/weighted_foodwebs.rds")

load("D:/R_Git/GoR-Food-webs/Data/uw_fw_GOR.RData")
load("D:/R_Git/GoR-Food-webs/Data/w_fw_GOR.RData")
              
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

SimulateExtinctions(fw.1981, Method = "Mostconnected")# fungerar inte -->crash, fatal error

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
