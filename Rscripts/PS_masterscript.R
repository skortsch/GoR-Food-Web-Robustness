####Load packagees#####
library(igraph)# Package containing tools to analyse food webs
library(MASS)
library(NetIndices)
library(multiweb)
library(fluxweb)
library(infomapecology)

#####Load data####
load("BalticFW.Rdata")
netmatrix
net
Metaweb.el<-net
Metaweb.el

colFG<- c("orange", "darkgrey", "blue", "green", "cyan")
info$colfg <- colFG[as.numeric(info$fg)]

sp.ab<-read.csv("sp_abrev.csv",  sep=";", check.names=FALSE)

plotfw(Metaweb.el, vertex.size =10, vertex.frame.color = "gray50", edge.arrow.size=0.2,
       edge.width=0.15, main="Metaweb", vertex.color=info$colfg,
       labcex=0.65, vertex.label.dist= 1, vertex.label=sp.ab[,2])
legend(1.05, y=0.95, legend = sp.ab[,2], cex = 0.6, bty = "n")
legend(1.2, y=0.95, legend = sp.ab[,1], cex = 0.6, bty = "n")

plotfw(W_Mweb.el, vertex.size =10, vertex.frame.color = "gray50", edge.arrow.size=0.2,
       edge.width=wid_Meta, main="Metaweb", vertex.color=info$colfg,
       labcex=0.65, vertex.label.dist= 1, vertex.label=sp.ab[,2])
legend(1.05, y=0.95, legend = sp.ab[,2], cex = 0.6, bty = "n")
legend(1.2, y=0.95, legend = sp.ab[,1], cex = 0.6, bty = "n")

### infomap Metaweb #####
ll_MW<-get.data.frame(Metaweb.el)
ll_MW# ger en edgelist men utan Weight kolumn, funkar inte i create_monolayer_object utan.
ll_MW<-cbind(ll_MW, Weight = 1)#lägg till Weight column med värde 1
ll_MW# funkar nu i create_monolayer_object
MW_obj<-create_monolayer_object(ll_MW, directed = T, bipartite = F)

MW_infMap<-run_infomap_monolayer(MW_obj, two_level = T, seed=123 ,flow_model = "rawdir", trials = 100, silent = T,
                                 ... = '--markov-time 1.05')

#view(MW_infMap$modules) # av nån naledning är ordningen på nodes annorlunda än i matrix -> blir fel med assign col
# tex är Myo & Clu i en egen modul och inte Gan & Spr

#Sortera modul-output enlingt ordningen arterna kommer i matrisen
IM.mod<-MW_infMap$modules
Sp.list<-c(row.names(netmatrix))
Orderd.IM.mod<-IM.mod[match(Sp.list, IM.mod$node_name),]
mod.aff_MW<-Orderd.IM.mod$module_level1
mod.col_MW<-c("green", "cyan")
mod.col_MW<-c("cyan","green")

#bättre plot
plotfw(Metaweb.el, vertex.size =10, vertex.frame.color = "gray50", edge.arrow.size=0.2,
     edge.width=0.15, main="Metaweb modularity", vertex.color=mod.col_MW[mod.aff_MW],
     labcex=0.65, vertex.label.dist= 1, vertex.label=sp.ab[,2], vertex.label.font=2)
legend(1.05, y=1.3, legend = sp.ab[,2], cex = 0.65, bty = "n", text.font = 2)
legend(1.2, y=1.3, legend = sp.ab[,1], cex = 0.65, bty = "n", text.font = 2)

#Weighted Metaweb
biomass <- info$meanB
netmatrix <- get.adjacency(net, sparse=F)
fluxes <- fluxing(netmatrix, biomass, info$losses, info$efficiencies, ef.level="prey")
fluxes <- fluxes*86.4
W_Mweb.el<- graph_from_adjacency_matrix(fluxes, weighted=TRUE)


## infomap Weighted
ll_WMW<-get.data.frame(W_Mweb.el) # gör link list med weights
ll_WMW # weighted linklist i formatet from, to, weight
WMW_obj<-create_monolayer_object(x=ll_WMW, directed = T, bipartite = F)
WMW_obj # är nu en monolayer objekt som går att använda

WMW_infmap<-run_infomap_monolayer(WMW_obj, two_level = T, seed=123 ,flow_model = "rawdir", trials = 100, silent = T,
                                  ... = '--markov-time 1.05')

IM.mod<-WMW_infmap$modules
Sp.list<-c(row.names(netmatrix))
Orderd.IM.mod<-IM.mod[match(Sp.list, IM.mod$node_name),]
mod.aff_WMW<-Orderd.IM.mod$module_level1
mod.col_WMW<-c("green", "cyan", "magenta")
mod.col_WMW<-c("cyan","magenta","green")

# från BFW script# Set visual parameters
Escale <- 30 #multiplying coefficient
Emin <- 0.1 #minimum width 
#Calculate the width of the arrows
wid_Meta <- Emin+(sqrt(E(W_Mweb.el)$weight)/max(sqrt(E(W_Mweb.el)$weight))*Escale)

plotfw(W_Mweb.el, vertex.size =10, vertex.frame.color = "gray50", edge.arrow.size=0.2,
       edge.width=wid_Meta, main="Metaweb weighted modularity", vertex.color=mod.col_WMW[mod.aff_WMW],
       labcex=0.65, vertex.label.dist= 1, vertex.label=sp.ab[,2], vertex.label.font=2)
legend(1.05, y=1.3, legend = sp.ab[,2], cex = 0.65, bty = "n", text.font = 2)
legend(1.2, y=1.3, legend = sp.ab[,1], cex = 0.65, bty = "n", text.font = 2)


##### Enskilda årtal
{
uwfw<-readRDS("unweighted_foodwebs.rds")
fw.1981<-uwfw[[1]]
fw.1986<-uwfw[[2]]
fw.1991<-uwfw[[3]]
fw.1996<-uwfw[[4]]
fw.2001<-uwfw[[5]]
fw.2006<-uwfw[[6]]
fw.2011<-uwfw[[7]]
}

{
fw.1981.el<-graph.adjacency(fw.1981)
fw.1986.el<-graph.adjacency(fw.1986)
fw.1991.el<-graph.adjacency(fw.1991)
fw.1996.el<-graph.adjacency(fw.1996)
fw.2001.el<-graph.adjacency(fw.2001)
fw.2006.el<-graph.adjacency(fw.2006)
fw.2011.el<-graph.adjacency(fw.2011)
}

{
  wfw<-readRDS("weighted_foodwebs.rds")
  wfw.1981<-wfw[[1]]
  wfw.1986<-wfw[[2]]
  wfw.1991<-wfw[[3]]
  wfw.1996<-wfw[[4]]
  wfw.2001<-wfw[[5]]
  wfw.2006<-wfw[[6]]
  wfw.2011<-wfw[[7]]
}
#edgelist
{
  wfw.1986.el<-graph.adjacency(wfw.1981, weighted = TRUE)
  wfw.1986.el<-graph.adjacency(wfw.1986, weighted = TRUE)
  wfw.1991.el<-graph.adjacency(wfw.1991, weighted = TRUE)
  wfw.1996.el<-graph.adjacency(wfw.1996, weighted = TRUE)
  wfw.2001.el<-graph.adjacency(wfw.2001, weighted = TRUE)
  wfw.2006.el<-graph.adjacency(wfw.2006, weighted = TRUE)
  wfw.2011.el<-graph.adjacency(wfw.2011, weighted = TRUE)
}
#1981
ll_1981_uw<-get.data.frame(fw.1981.el)
ll_1981_uw<-cbind(ll_1981_uw, Weight = 1)#lägg till Weight column med värde 1
uw_nw_obj_1981<-create_monolayer_object(ll_1981_uw, directed = T, bipartite = F)
IM_uw_1981<-run_infomap_monolayer(uw_nw_obj_1981, two_level = T, flow_model = "rawdir", trials = 100, silent = T,
                                  ... = '--markov-time 1.05')

IM.mod<-IM_uw_1981$modules
Sp.list<-c(row.names(fw.1981))
Orderd.IM.mod<-IM.mod[match(Sp.list, IM.mod$node_name),]
mod.aff_UW1981<-Orderd.IM.mod$module_level1
mod.col_UW1981<-c("green", "cyan")
mod.col_UW1981<-c("cyan","green")
plotfw(fw.1981.el, col = 1, labcex = 0.65, vertex.label.dist= 1, edge.arrow.size=0.2, vertex.size =10, main="1981", 
       vertex.color=mod.col_UW1981[mod.aff_UW1981], vertex.label=sp.ab[,3])
legend(1.05, y=0.95, legend = sp.ab[,2], cex = 0.6, bty = "n")
legend(1.2, y=0.95, legend = sp.ab[,1], cex = 0.6, bty = "n")
modularity(fw.1981.el,mod.aff_UW1981)

#1981 weight
ll_1981w<-get.data.frame(wfw.1981.el)
w_1981<-create_monolayer_object(ll_1981w, directed = T, bipartite = F)
IM_w_1981<-run_infomap_monolayer(w_1981, two_level = T, flow_model = "rawdir", trials = 100, silent = T,
                                 ... = '--markov-time 1.05')


IM.mod<-IM_w_1981$modules
Sp.list<-c(row.names(wfw.1981))
Orderd.IM.mod<-IM.mod[match(Sp.list, IM.mod$node_name),]
mod.aff_W1981<-Orderd.IM.mod$module_level1
mod.col_W1981<-c("green", "cyan", "magenta")

wid_1981 <- Emin+(sqrt(E(wfw.1981.el)$weight)/max(sqrt(E(wfw.1981.el)$weight))*Escale)

plotfw(wfw.1981.el, col = 1, labcex = 0.65, vertex.label.dist= 1, edge.arrow.size=0.2, vertex.size =10, 
       edge.width=wid_1981, main="1981 weighted", vertex.color=mod.col_W1981[mod.aff_W1981],
       vertex.label=sp.ab[,3])
legend(1.05, y=0.95, legend = sp.ab[,2], cex = 0.6, bty = "n")
legend(1.2, y=0.95, legend = sp.ab[,1], cex = 0.6, bty = "n")

# 1986 uw
ll_1986_uw<-get.data.frame(fw.1986.el)
ll_1986_uw<-cbind(ll_1986_uw, Weight = 1)#lägg till Weight column med värde 1
uw_nw_obj_1986<-create_monolayer_object(ll_1986_uw, directed = T, bipartite = F)
IM_uw_1986<-run_infomap_monolayer(uw_nw_obj_1986, two_level = T, flow_model = "rawdir", seed = 123, trials = 100, silent = T,
                                  ... = '--markov-time 1.05')

IM.mod<-IM_uw_1986$modules
Sp.list<-c(row.names(fw.1986))
Orderd.IM.mod<-IM.mod[match(Sp.list, IM.mod$node_name),]
mod.aff_UW1986<-Orderd.IM.mod$module_level1
mod.col_UW1986<-c("green", "cyan")
mod.col_UW1986<-c("cyan", "green")
plotfw(fw.1986.el, col = 1, labcex = 0.65, vertex.label.dist= 1, edge.arrow.size=0.2, vertex.size =10, main="1986", 
       vertex.color=mod.col_UW1986[mod.aff_UW1986],vertex.label=sp.ab[,4])
legend(1.05, y=0.95, legend = sp.ab[,2], cex = 0.6, bty = "n")
legend(1.2, y=0.95, legend = sp.ab[,1], cex = 0.6, bty = "n")

# 1986 weighted
ll_1986w<-get.data.frame(wfw.1986.el)
w_1986<-create_monolayer_object(ll_1986w, directed = T, bipartite = F)
IM_w_1986<-run_infomap_monolayer(w_1986, two_level = T, flow_model = "rawdir", trials = 100, silent = T,
                                 ... = '--markov-time 1.05')

IM.mod<-IM_w_1986$modules
Sp.list<-c(row.names(wfw.1986))
Orderd.IM.mod<-IM.mod[match(Sp.list, IM.mod$node_name),]
mod.aff_W1986<-Orderd.IM.mod$module_level1
mod.col_W1986<-c("green", "cyan", "magenta")

wid_1986 <- Emin+(sqrt(E(wfw.1986.el)$weight)/max(sqrt(E(wfw.1986.el)$weight))*Escale)

plotfw(wfw.1986.el, col = 1, labcex = 0.65, vertex.label.dist= 1, edge.arrow.size=0.2, vertex.size =10, 
       edge.width=wid_1986, main="1986 weighted", vertex.color=mod.col_W1986[mod.aff_W1986]
       ,vertex.label=sp.ab[,4])
legend(1.05, y=0.95, legend = sp.ab[,2], cex = 0.6, bty = "n")
legend(1.2, y=0.95, legend = sp.ab[,1], cex = 0.6, bty = "n")

#1991 uw
ll_1991_uw<-get.data.frame(fw.1991.el)
ll_1991_uw<-cbind(ll_1991_uw, Weight = 1)#lägg till Weight column med värde 1
uw_nw_obj_1991<-create_monolayer_object(ll_1991_uw, directed = T, bipartite = F)
IM_uw_1991<-run_infomap_monolayer(uw_nw_obj_1991, two_level = T, flow_model = "rawdir", seed = 123, trials = 100, silent = T,
                                  ... = '--markov-time 1.05')

IM.mod<-IM_uw_1991$modules
Sp.list<-c(row.names(fw.1991))
Orderd.IM.mod<-IM.mod[match(Sp.list, IM.mod$node_name),]
mod.aff_UW1991<-Orderd.IM.mod$module_level1
mod.col_UW1991<-c("cyan","green")

plotfw(fw.1991.el, col = 1, labcex = 0.65, vertex.label.dist= 1, edge.arrow.size=0.2, vertex.size =10, main="1991", 
       vertex.color=mod.col_UW1991[mod.aff_UW1991],vertex.label=sp.ab[,5])
legend(1.05, y=0.95, legend = sp.ab[,2], cex = 0.6, bty = "n")
legend(1.2, y=0.95, legend = sp.ab[,1], cex = 0.6, bty = "n")

#1991 weight
ll_1991w<-get.data.frame(wfw.1991.el)
w_1991<-create_monolayer_object(ll_1991w, directed = T, bipartite = F)
IM_w_1991<-run_infomap_monolayer(w_1991, two_level = T, flow_model = "rawdir", seed = 123, trials = 100, silent = T,
                                 ... = '--markov-time 1.05')

IM.mod<-IM_w_1991$modules
Sp.list<-c(row.names(wfw.1991))
Orderd.IM.mod<-IM.mod[match(Sp.list, IM.mod$node_name),]
mod.aff_W1991<-Orderd.IM.mod$module_level1
mod.col_W1991<-c("green", "cyan", "magenta", "blue")

wid_1991 <- Emin+(sqrt(E(wfw.1991.el)$weight)/max(sqrt(E(wfw.1991.el)$weight))*Escale)

plotfw(wfw.1991.el, col = 1, labcex = 0.65, vertex.label.dist= 1, edge.arrow.size=0.2, vertex.size =10, 
       edge.width=wid_1991, main="1991 weighted", vertex.color=mod.col_W1991[mod.aff_W1991]
       ,vertex.label=sp.ab[,5])
legend(1.05, y=0.95, legend = sp.ab[,2], cex = 0.6, bty = "n")
legend(1.2, y=0.95, legend = sp.ab[,1], cex = 0.6, bty = "n")

# 1996 uw
ll_1996_uw<-get.data.frame(fw.1996.el)
ll_1996_uw<-cbind(ll_1996_uw, Weight = 1)#lägg till Weight column med värde 1
uw_nw_obj_1996<-create_monolayer_object(ll_1996_uw, directed = T, bipartite = F)
IM_uw_1996<-run_infomap_monolayer(uw_nw_obj_1996, two_level = T, flow_model = "rawdir", seed = 123, trials = 100, silent = T,
                                  ... = '--markov-time 1.05')


IM.mod<-IM_uw_1996$modules
Sp.list<-c(row.names(fw.1996))
Orderd.IM.mod<-IM.mod[match(Sp.list, IM.mod$node_name),]
mod.aff_UW1996<-Orderd.IM.mod$module_level1
mod.col_UW1996<-c("green", "cyan", "magenta")
mod.col_UW1996<-c("cyan", "green", "magenta")
plotfw(fw.1996.el, col = 1, labcex = 0.65, vertex.label.dist= 1, edge.arrow.size=0.2, vertex.size =10, main="1996", 
       vertex.color=mod.col_UW1996[mod.aff_UW1996],vertex.label=sp.ab[,6])
legend(1.05, y=0.95, legend = sp.ab[,2], cex = 0.6, bty = "n")
legend(1.2, y=0.95, legend = sp.ab[,1], cex = 0.6, bty = "n")

# 1996 weight
ll_1996w<-get.data.frame(wfw.1996.el)
w_1996<-create_monolayer_object(ll_1996w, directed = T, bipartite = F)
IM_w_1996<-run_infomap_monolayer(w_1996, two_level = T, flow_model = "rawdir", seed = 123, trials = 100, silent = T,
                                 ... = '--markov-time 1.05')

IM.mod<-IM_w_1996$modules
Sp.list<-c(row.names(wfw.1996))
Orderd.IM.mod<-IM.mod[match(Sp.list, IM.mod$node_name),]
mod.aff_W1996<-Orderd.IM.mod$module_level1
mod.col_W1996<-c("cyan","green")
mod.col_W1996<-c("green", "cyan", "magenta")
wid_1996 <- Emin+(sqrt(E(wfw.1996.el)$weight)/max(sqrt(E(wfw.1996.el)$weight))*Escale)

plotfw(wfw.1996.el, col = 1, labcex = 0.65, vertex.label.dist= 1, edge.arrow.size=0.2, vertex.size =10, 
       edge.width=wid_1996, main="1996 weighted", vertex.color=mod.col_W1996[mod.aff_W1996]
       ,vertex.label=sp.ab[,6])
legend(1.05, y=0.95, legend = sp.ab[,2], cex = 0.6, bty = "n")
legend(1.2, y=0.95, legend = sp.ab[,1], cex = 0.6, bty = "n")


# 2001 uw
ll_2001_uw<-get.data.frame(fw.2001.el)
ll_2001_uw<-cbind(ll_2001_uw, Weight = 1)#lägg till Weight column med värde 1
uw_nw_obj_2001<-create_monolayer_object(ll_2001_uw, directed = T, bipartite = F)
IM_uw_2001<-run_infomap_monolayer(uw_nw_obj_2001, two_level = T, flow_model = "rawdir", seed = 123, trials = 100, silent = T,
                                  ... = '--markov-time 1.05')

IM.mod<-IM_uw_2001$modules
Sp.list<-c(row.names(fw.2001))
Orderd.IM.mod<-IM.mod[match(Sp.list, IM.mod$node_name),]
mod.aff_UW2001<-Orderd.IM.mod$module_level1
mod.col_UW2001<-c("green", "cyan", "magenta")
mod.col_UW2001<-c("cyan", "green", "magenta")
plotfw(fw.2001.el, col = 1, labcex = 0.65, vertex.label.dist= 1, edge.arrow.size=0.2, vertex.size =10, main="2001", 
       vertex.color=mod.col_UW2001[mod.aff_UW2001],vertex.label=sp.ab[,7])
legend(1.05, y=0.95, legend = sp.ab[,2], cex = 0.6, bty = "n")
legend(1.2, y=0.95, legend = sp.ab[,1], cex = 0.6, bty = "n")

# 2001 weight
ll_2001w<-get.data.frame(wfw.2001.el)
w_2001<-create_monolayer_object(ll_2001w, directed = T, bipartite = F)
IM_w_2001<-run_infomap_monolayer(w_2001, two_level = T, flow_model = "rawdir", seed = 123, trials = 100, silent = T,
                                 ... = '--markov-time 1.05')

IM.mod<-IM_w_2001$modules
Sp.list<-c(row.names(wfw.2001))
Orderd.IM.mod<-IM.mod[match(Sp.list, IM.mod$node_name),]
mod.aff_W2001<-Orderd.IM.mod$module_level1
mod.col_W2001<-c("green", "cyan", "magenta", "blue")

wid_2001<- Emin+(sqrt(E(wfw.2001.el)$weight)/max(sqrt(E(wfw.2001.el)$weight))*Escale)

plotfw(wfw.2001.el, col = 1, labcex = 0.65, vertex.label.dist= 1, edge.arrow.size=0.2, vertex.size =10, 
       edge.width=wid_2001, main="2001 weighted", vertex.color=mod.col_W2001[mod.aff_W2001]
       ,vertex.label=sp.ab[,7])
legend(1.05, y=0.95, legend = sp.ab[,2], cex = 0.6, bty = "n")
legend(1.2, y=0.95, legend = sp.ab[,1], cex = 0.6, bty = "n")


# 2006 uw
ll_2006_uw<-get.data.frame(fw.2006.el)
ll_2006_uw<-cbind(ll_2006_uw, Weight = 1)#lägg till Weight column med värde 1
uw_nw_obj_2006<-create_monolayer_object(ll_2006_uw, directed = T, bipartite = F)
IM_uw_2006<-run_infomap_monolayer(uw_nw_obj_2006, two_level = T, flow_model = "rawdir", trials = 100, silent = T,
                                  ... = '--markov-time 1.05')

IM.mod<-IM_uw_2006$modules
Sp.list<-c(row.names(fw.2006))
Orderd.IM.mod<-IM.mod[match(Sp.list, IM.mod$node_name),]
mod.aff_UW2006<-Orderd.IM.mod$module_level1
mod.col_UW2006<-c("cyan","green")

plotfw(fw.2006.el, col = 1, labcex = 0.65, vertex.label.dist= 1, edge.arrow.size=0.2, vertex.size =10, main="2006", 
       vertex.color=mod.col_UW2006[mod.aff_UW2006],vertex.label=sp.ab[,8])
legend(1.05, y=0.95, legend = sp.ab[,2], cex = 0.6, bty = "n")
legend(1.2, y=0.95, legend = sp.ab[,1], cex = 0.6, bty = "n")

# 2006 weight
ll_2006w<-get.data.frame(wfw.2006.el)
w_2006<-create_monolayer_object(ll_2006w, directed = T, bipartite = F)
IM_w_2006<-run_infomap_monolayer(w_2006, two_level = T, flow_model = "rawdir", seed = 123, trials = 100, silent = T,
                                 ... = '--markov-time 1.05')

IM.mod<-IM_w_2006$modules
Sp.list<-c(row.names(wfw.2006))
Orderd.IM.mod<-IM.mod[match(Sp.list, IM.mod$node_name),]
mod.aff_W2006<-Orderd.IM.mod$module_level1
mod.col_W2006<-c("green", "cyan", "magenta","blue")
wid_2006<- Emin+(sqrt(E(wfw.2006.el)$weight)/max(sqrt(E(wfw.2006.el)$weight))*Escale)

plotfw(wfw.2006.el, col = 1, labcex = 0.65, vertex.label.dist= 1, edge.arrow.size=0.2, vertex.size =10, 
       edge.width=wid_2006, main="2006 weighted", vertex.color=mod.col_W2006[mod.aff_W2006]
       ,vertex.label=sp.ab[,8])
legend(1.05, y=0.95, legend = sp.ab[,2], cex = 0.6, bty = "n")
legend(1.2, y=0.95, legend = sp.ab[,1], cex = 0.6, bty = "n")


# 2011 uw
ll_2011_uw<-get.data.frame(fw.2011.el)
ll_2011_uw<-cbind(ll_2011_uw, Weight = 1)#lägg till Weight column med värde 1
uw_nw_obj_2011<-create_monolayer_object(ll_2011_uw, directed = T, bipartite = F)
IM_uw_2011<-run_infomap_monolayer(uw_nw_obj_2011, two_level = T, flow_model = "rawdir", seed = 123, trials = 100, silent = T,
                                  ... = '--markov-time 1.05')

IM.mod<-IM_uw_2011$modules
Sp.list<-c(row.names(fw.2011))
Orderd.IM.mod<-IM.mod[match(Sp.list, IM.mod$node_name),]
mod.aff_UW2011<-Orderd.IM.mod$module_level1
mod.col_UW2011<-c("cyan","green")

plotfw(fw.2011.el, labcex = 0.65, vertex.label.dist= 1, edge.arrow.size=0.2, vertex.size =10, main="2011", 
       vertex.color=mod.col_UW2011[mod.aff_UW2011],vertex.label=sp.ab[,9])
legend(1.05, y=0.95, legend = sp.ab[,2], cex = 0.6, bty = "n")
legend(1.2, y=0.95, legend = sp.ab[,1], cex = 0.6, bty = "n")

#2011 weighted
ll_2011w<-get.data.frame(wfw.2011.el)
w_2011<-create_monolayer_object(ll_2011w, directed = T, bipartite = F)
IM_w_2011<-run_infomap_monolayer(w_2011, two_level = T, seed = 123, flow_model = "rawdir", trials = 100, silent = T,
                                 ... = '--markov-time 1.05')

IM.mod<-IM_w_2011$modules
Sp.list<-c(row.names(wfw.2011))
Orderd.IM.mod<-IM.mod[match(Sp.list, IM.mod$node_name),]
mod.aff_W2011<-Orderd.IM.mod$module_level1
mod.col_W2011<-c("green", "cyan", "magenta", "blue", "yellow", "red")

wid_2011<- Emin+(sqrt(E(wfw.2011.el)$weight)/max(sqrt(E(wfw.2011.el)$weight))*Escale)

plotfw(wfw.2011.el, labcex = 0.65, vertex.label.dist= 1, edge.arrow.size=0.2, vertex.size =10, 
       edge.width=wid_2011, main="2011 weighted", vertex.color=mod.col_W2011[mod.aff_W2011],
       vertex.label=sp.ab[,9])
legend(1.05, y=0.95, legend = sp.ab[,2], cex = 0.6, bty = "n")
legend(1.2, y=0.95, legend = sp.ab[,1], cex = 0.6, bty = "n")


# ... = '--markov-time 1.05'

########topological roles###########

Top_rol_Meta<-Troles(Metaweb.el, mod.aff_MW)
Top_rol_Meta<-round(Top_rol_Meta, 2)
Top_rol_Meta[3]<-cbind(row.names(netmatrix))
plot(Top_rol_Meta[,2],Top_rol_Meta[,1], xlab = "PC (among_mod_conn)", ylab = "z-score (within_mod_degree)", xlim = c(-0.2,0.8),
     ylim= c(-3,3), main = "Metaweb_uw", col= info$colfg, cex=1.5, pch=16)
text(Top_rol_Meta[,2], Top_rol_Meta[,1]-0.05, labels = sp.ab[,2], cex = 0.5)
text(x=0.75, y=2.7, lab= "Netwrok connector", cex=0.8)
text(x=0.75, y=1.9, lab= "connectors", cex=0.8)
text(x=-0.1, y=2.7, lab= "module connectors", cex=0.8)
text(x=-0.1, y=1.9, lab= "peripherals", cex=0.8)
abline(v=0.625, h=2.5)

Top_rol_Meta_W<-Troles(Metaweb.el, mod.aff_WMW)
Top_rol_Meta_W<-round(Top_rol_Meta_W, 2)
Top_rol_Meta_W[,3]<-cbind(row.names(netmatrix))
plot(Top_rol_Meta_W[,2],Top_rol_Meta_W[,1], xlab = "PC (among_mod_conn)", ylab = "z-score (within_mod_degree)",xlim = c(-0.2,0.8), 
     ylim= c(-3,3), main = "Metaweb_W", col= info$colfg, cex=1.5, pch=16)
text(Top_rol_Meta_W[,2],Top_rol_Meta_W[,1]-0.05, labels = sp.ab[,2], cex = 0.5)
text(x=0.75, y=2.7, lab= "Netwrok connector", cex=0.8)
text(x=0.75, y=1.9, lab= "module connectors", cex=0.8)
text(x=-0.1, y=2.7, lab= "module hub", cex=0.8)
text(x=-0.1, y=1.9, lab= "peripherals", cex=0.8)
abline(v=0.625, h=2.5)

{
  sp.fg.1981<-read.csv("fg.1981.csv", sep=";", check.names = F, stringsAsFactors = T)
  sp.fg.1986<-read.csv("fg.1986.csv", sep=";", check.names = F, stringsAsFactors = T)
  sp.fg.1991<-read.csv("fg.1991.csv", sep=";", check.names = F, stringsAsFactors = T)
  sp.fg.1996<-read.csv("fg.1996.csv", sep=";", check.names = F, stringsAsFactors = T)
  sp.fg.2001<-read.csv("fg.2001.csv", sep=";", check.names = F, stringsAsFactors = T)
  sp.fg.2006<-read.csv("fg.2006.csv", sep=";", check.names = F, stringsAsFactors = T)
  sp.fg.2011<-read.csv("fg.2011.csv", sep=";", check.names = F, stringsAsFactors = T)
}

Top_rol_1981<-Troles(fw.1981.el, mod.aff_UW1981)
Top_rol_1981<-round(Top_rol_1981, 2)
Top_rol_1981[,3]<-cbind(row.names(fw.1981))
plot(Top_rol_1981[,2],Top_rol_1981[,1], xlab = "PC (among_mod_conn)", ylab = "z-score (within_mod_degree)",xlim = c(-0.2,0.8), 
     ylim= c(-3,3), main = "1981_uw", col= colFG[as.numeric(sp.fg.1981[,2])], cex=1.5, pch=16)
text(Top_rol_1981[,2],Top_rol_1981[,1]-0.05, labels = sp.ab[,3], cex = 0.5)
text(x=0.75, y=2.7, lab= "Hub netwrok connector", cex=0.8)
text(x=0.75, y=2.2, lab= "module connectors", cex=0.8)
text(x=-0.1, y=2.7, lab= "module hub", cex=0.8)
text(x=-0.1, y=2.2, lab= "peripherals", cex=0.8)
abline(v=0.625, h=2.5)

Top_rol_1981_w<-Troles(fw.1981.el, mod.aff_W1981)
Top_rol_1981_w<-round(Top_rol_1981_w, 2)
Top_rol_1981_w[,3]<-cbind(row.names(fw.1981))
plot(Top_rol_1981_w[,2],Top_rol_1981_w[,1], xlab = "PC (among_mod_conn)", ylab = "z-score (within_mod_degree)",xlim = c(-0.2,0.8), 
     ylim= c(-3,3), main = "1981_w", col= colFG[as.numeric(sp.fg.1981[,2])], cex=1.5, pch=16)
text(Top_rol_1981_w[,2],Top_rol_1981_w[,1]-0.05, labels = sp.ab[,3], cex = 0.5)
text(x=0.75, y=2.7, lab= "Hub netwrok connector", cex=0.8)
text(x=0.75, y=2.2, lab= "module connectors", cex=0.8)
text(x=-0.1, y=2.7, lab= "module hub", cex=0.8)
text(x=-0.1, y=2.2, lab= "peripherals", cex=0.8)
abline(v=0.625, h=2.5)


Top_rol_1986<-Troles(fw.1986.el, mod.aff_UW1986)
Top_rol_1986<-round(Top_rol_1986, 2)
Top_rol_1986[,3]<-cbind(row.names(fw.1986))
plot(Top_rol_1986[,2],Top_rol_1986[,1], xlab = "PC (among_mod_conn)", ylab = "z-score (within_mod_degree)",xlim = c(-0.2,0.8), 
     ylim= c(-3,3), main = "1986_uw", col=colFG[as.numeric(sp.fg.1986[,2])], cex=1.5, pch=16)
text(Top_rol_1986[,2],Top_rol_1986[,1]-0.05, labels = sp.ab[,4], cex = 0.5)
text(x=0.75, y=2.7, lab= "Hub netwrok connector", cex=0.8)
text(x=0.75, y=2.2, lab= "module connectors", cex=0.8)
text(x=-0.1, y=2.7, lab= "module hub", cex=0.8)
text(x=-0.1, y=2.2, lab= "peripherals", cex=0.8)
abline(v=0.625, h=2.5)

Top_rol_1986_w<-Troles(fw.1986.el, mod.aff_W1986)
Top_rol_1986_w<-round(Top_rol_1986_w, 2)
Top_rol_1986_w[,3]<-cbind(row.names(fw.1986))
plot(Top_rol_1986_w[,2],Top_rol_1986_w[,1], xlab = "PC (among_mod_conn)", ylab = "z-score (within_mod_degree)",xlim = c(-0.2,0.8), 
     ylim= c(-3,3), main = "1986_w", col=colFG[as.numeric(sp.fg.1986[,2])], cex=1.5, pch=16)
text(Top_rol_1986_w[,2],Top_rol_1986_w[,1]-0.05, labels = sp.ab[,4], cex = 0.5)
text(x=0.75, y=2.7, lab= "Hub netwrok connector", cex=0.8)
text(x=0.75, y=2.2, lab= "module connectors", cex=0.8)
text(x=-0.1, y=2.7, lab= "module hub", cex=0.8)
text(x=-0.1, y=2.2, lab= "peripherals", cex=0.8)
abline(v=0.625, h=2.5)


Top_rol_1991<-Troles(fw.1991.el, mod.aff_UW1991)
Top_rol_1991<-round(Top_rol_1991, 2)
Top_rol_1991[,3]<-cbind(row.names(fw.1991))
plot(Top_rol_1991[,2],Top_rol_1991[,1], xlab = "PC (among_mod_conn)", ylab = "z-score (within_mod_degree)",xlim = c(-0.2,0.8), 
     ylim= c(-3,3), main = "1991_uw", col=colFG[as.numeric(sp.fg.1991[,2])], cex=1.5, pch=16)
text(Top_rol_1991[,2],Top_rol_1991[,1]-0.05, labels = sp.ab[,5], cex = 0.5)
text(x=0.75, y=2.7, lab= "Hub netwrok connector", cex=0.8)
text(x=0.75, y=2.2, lab= "module connectors", cex=0.8)
text(x=-0.1, y=2.7, lab= "module hub", cex=0.8)
text(x=-0.1, y=2.2, lab= "peripherals", cex=0.8)
abline(v=0.625, h=2.5)

Top_rol_1991_w<-Troles(fw.1991.el, mod.aff_W1991)
Top_rol_1991_w<-round(Top_rol_1991_w, 2)
Top_rol_1991_w[,3]<-cbind(row.names(fw.1991))
plot(Top_rol_1991_w[,2],Top_rol_1991_w[,1], xlab = "PC (among_mod_conn)", ylab = "z-score (within_mod_degree)",xlim = c(-0.2,0.8), 
     ylim= c(-3,3), main = "1991_w", col=colFG[as.numeric(sp.fg.1991[,2])], cex=1.5, pch=16)
text(Top_rol_1991_w[,2],Top_rol_1991_w[,1]-0.05, labels = sp.ab[,5], cex = 0.5)
text(x=0.75, y=2.7, lab= "Hub netwrok connector", cex=0.8)
text(x=0.75, y=2.2, lab= "module connectors", cex=0.8)
text(x=-0.1, y=2.7, lab= "module hub", cex=0.8)
text(x=-0.1, y=2.2, lab= "peripherals", cex=0.8)
abline(v=0.625, h=2.5)


Top_rol_1996<-Troles(fw.1996.el, mod.aff_UW1996)
Top_rol_1996<-round(Top_rol_1996, 2)
Top_rol_1996[,3]<-cbind(row.names(fw.1996))
plot(Top_rol_1996[,2],Top_rol_1996[,1], xlab = "PC (among_mod_conn)", ylab = "z-score (within_mod_degree)",xlim = c(-0.2,0.8), 
     ylim= c(-3,3), main = "1996_uw", col=colFG[as.numeric(sp.fg.1996[,2])], cex=1.5, pch=16)
text(Top_rol_1996[,2],Top_rol_1996[,1]-0.05, labels = sp.ab[,6], cex = 0.5)
text(x=0.75, y=2.7, lab= "Hub netwrok connector", cex=0.8)
text(x=0.75, y=2.2, lab= "module connectors", cex=0.8)
text(x=-0.1, y=2.7, lab= "module hub", cex=0.8)
text(x=-0.1, y=2.2, lab= "peripherals", cex=0.8)
abline(v=0.625, h=2.5)

Top_rol_1996_w<-Troles(fw.1996.el, mod.aff_W1996)
Top_rol_1996_w<-round(Top_rol_1996_w, 2)
Top_rol_1996_w[,3]<-cbind(row.names(fw.1996))
plot(Top_rol_1996_w[,2],Top_rol_1996_w[,1], xlab = "PC (among_mod_conn)", ylab = "z-score (within_mod_degree)",xlim = c(-0.2,0.8), 
     ylim= c(-3,3), main = "1996_w", col=colFG[as.numeric(sp.fg.1996[,2])], cex=1.5, pch=16)
text(Top_rol_1996_w[,2],Top_rol_1996_w[,1]-0.05, labels = sp.ab[,6], cex = 0.5)
text(x=0.75, y=2.7, lab= "Hub netwrok connector", cex=0.8)
text(x=0.75, y=2.2, lab= "module connectors", cex=0.8)
text(x=-0.1, y=2.7, lab= "module hub", cex=0.8)
text(x=-0.1, y=2.2, lab= "peripherals", cex=0.8)
abline(v=0.625, h=2.5)



Top_rol_2001<-Troles(fw.2001.el, mod.aff_UW2001)
Top_rol_2001<-round(Top_rol_2001, 2)
Top_rol_2001[,3]<-cbind(row.names(fw.2001))
plot(Top_rol_2001[,2],Top_rol_2001[,1], xlab = "PC (among_mod_conn)", ylab = "z-score (within_mod_degree)",xlim = c(-0.2,0.8), 
     ylim= c(-3,3), main = "2001_uw", col=colFG[as.numeric(sp.fg.2001[,2])], cex=1.5, pch=16)
text(Top_rol_2001[,2],Top_rol_2001[,1]-0.05, labels = sp.ab[,7], cex = 0.5)
text(x=0.75, y=2.7, lab= "Hub netwrok connector", cex=0.8)
text(x=0.75, y=2.2, lab= "module connectors", cex=0.8)
text(x=-0.1, y=2.7, lab= "module hub", cex=0.8)
text(x=-0.1, y=2.2, lab= "peripherals", cex=0.8)
abline(v=0.625, h=2.5)

Top_rol_2001_w<-Troles(fw.2001.el, mod.aff_W2001)
Top_rol_2001_w<-round(Top_rol_2001_w, 2)
Top_rol_2001_w[,3]<-cbind(row.names(fw.2001))
plot(Top_rol_2001_w[,2],Top_rol_2001_w[,1], xlab = "PC (among_mod_conn)", ylab = "z-score (within_mod_degree)",xlim = c(-0.2,0.8), 
     ylim= c(-3,3), main = "2001_w", col=colFG[as.numeric(sp.fg.2001[,2])], cex=1.5, pch=16)
text(Top_rol_2001_w[,2],Top_rol_2001_w[,1]-0.05, labels = sp.ab[,7], cex = 0.5)
text(x=0.75, y=2.7, lab= "Hub netwrok connector", cex=0.8)
text(x=0.75, y=2.2, lab= "module connectors", cex=0.8)
text(x=-0.1, y=2.7, lab= "module hub", cex=0.8)
text(x=-0.1, y=2.2, lab= "peripherals", cex=0.8)
abline(v=0.625, h=2.5)



Top_rol_2006<-Troles(fw.2006.el, mod.aff_UW2006)
Top_rol_2006<-round(Top_rol_2006, 2)
Top_rol_2006[,3]<-cbind(row.names(fw.2006))
plot(Top_rol_2006[,2],Top_rol_2006[,1], xlab = "PC (among_mod_conn)", ylab = "z-score (within_mod_degree)",xlim = c(-0.2,0.8), 
     ylim= c(-3,3), main = "2006_uw", col=colFG[as.numeric(sp.fg.2006[,2])], cex=1.5, pch=16)
text(Top_rol_2006[,2],Top_rol_2006[,1]-0.05, labels = sp.ab[,8], cex = 0.5)
text(x=0.75, y=2.7, lab= "Hub netwrok connector", cex=0.8)
text(x=0.75, y=2.2, lab= "module connectors", cex=0.8)
text(x=-0.1, y=2.7, lab= "module hub", cex=0.8)
text(x=-0.1, y=2.2, lab= "peripherals", cex=0.8)
abline(v=0.625, h=2.5)

Top_rol_2006_w<-Troles(fw.2006.el, mod.aff_W2006)
Top_rol_2006_w<-round(Top_rol_2006_w, 2)
Top_rol_2006_w[,3]<-cbind(row.names(fw.2006))
plot(Top_rol_2006_w[,2],Top_rol_2006_w[,1], xlab = "PC (among_mod_conn)", ylab = "z-score (within_mod_degree)",xlim = c(-0.2,0.8), 
     ylim= c(-3,3), main = "2006_w", col=colFG[as.numeric(sp.fg.2006[,2])], cex=1.5, pch=16)
text(Top_rol_2006_w[,2],Top_rol_2006_w[,1]-0.05, labels = sp.ab[,8], cex = 0.5)
text(x=0.75, y=2.7, lab= "Hub netwrok connector", cex=0.8)
text(x=0.75, y=2.2, lab= "module connectors", cex=0.8)
text(x=-0.1, y=2.7, lab= "module hub", cex=0.8)
text(x=-0.1, y=2.2, lab= "peripherals", cex=0.8)
abline(v=0.625, h=2.5)



Top_rol_2011<-Troles(fw.2011.el, mod.aff_UW2011)
Top_rol_2011<-round(Top_rol_2011, 2)
Top_rol_2011[,3]<-cbind(row.names(fw.2011))
plot(Top_rol_2011[,2],Top_rol_2011[,1], xlab = "PC (among_mod_conn)", ylab = "z-score (within_mod_degree)",xlim = c(-0.2,0.8), 
     ylim= c(-3,3), main = "2011_uw", col=colFG[as.numeric(sp.fg.2011[,2])], cex=1.5, pch=16)
text(Top_rol_2011[,2],Top_rol_2011[,1]-0.05, labels = sp.ab[,9], cex = 0.5)
text(x=0.75, y=2.7, lab= "Hub netwrok connector", cex=0.8)
text(x=0.75, y=2.2, lab= "module connectors", cex=0.8)
text(x=-0.1, y=2.7, lab= "module hub", cex=0.8)
text(x=-0.1, y=2.2, lab= "peripherals", cex=0.8)
abline(v=0.625, h=2.5)

Top_rol_2011_w<-Troles(fw.2011.el, mod.aff_W2011)
Top_rol_2011_w<-round(Top_rol_2011_w, 2)
Top_rol_2011_w[,3]<-cbind(row.names(fw.2011))
plot(Top_rol_2011_w[,2],Top_rol_2011_w[,1], xlab = "PC (among_mod_conn)", ylab = "z-score (within_mod_degree)",xlim = c(-0.2,0.8), 
     ylim= c(-3,3), main = "2011_w", col=colFG[as.numeric(sp.fg.2011[,2])], cex=1.5, pch=16)
text(Top_rol_2011_w[,2],Top_rol_2011_w[,1]-0.05, labels = sp.ab[,9], cex = 0.5)
text(x=0.75, y=2.7, lab= "Hub netwrok connector", cex=0.8)
text(x=0.75, y=2.2, lab= "module connectors", cex=0.8)
text(x=-0.1, y=2.7, lab= "module hub", cex=0.8)
text(x=-0.1, y=2.2, lab= "peripherals", cex=0.8)
abline(v=0.625, h=2.5)

#########Plot av top_roles för varje FG#########
#zooplankton
Top_rol_1981[,4]<-c(sp.fg.1981[2], use.names = F)
Top_rol_1981[,5]<-c(1981)
Trol_uw1981<-filter(Top_rol_1981, grepl("Zooplankton", V4))
Top_rol_1986[,4]<-c(sp.fg.1986[2], use.names = F)
Top_rol_1986[,5]<-c(1986)
Trol_uw1986<-filter(Top_rol_1986, grepl("Zooplankton", V4))
Top_rol_1991[,4]<-c(sp.fg.1991[2], use.names = F)
Top_rol_1991[,5]<-c(1991)
Trol_uw1991<-filter(Top_rol_1991, grepl("Zooplankton", V4))
Top_rol_1996[,4]<-c(sp.fg.1996[2], use.names = F)
Top_rol_1996[,5]<-c(1996)
Trol_uw1996<-filter(Top_rol_1996, grepl("Zooplankton", V4))
Top_rol_2001[,4]<-c(sp.fg.2001[2], use.names = F)
Top_rol_2001[,5]<-c(2001)
Trol_uw2001<-filter(Top_rol_2001, grepl("Zooplankton", V4))
Top_rol_2006[,4]<-c(sp.fg.2006[2], use.names = F)
Top_rol_2006[,5]<-c(2006)
Trol_uw2006<-filter(Top_rol_2006, grepl("Zooplankton", V4))
Top_rol_2011[,4]<-c(sp.fg.2011[2], use.names = F)
Top_rol_2011[,5]<-c(2011)
Trol_uw2011<-filter(Top_rol_2011, grepl("Zooplankton", V4))
Trol_zoopl_uw<-vctrs::vec_c(Trol_uw1981, Trol_uw1986, Trol_uw1991, Trol_uw1996, Trol_uw2001, Trol_uw2006, Trol_uw2011)
#plot(Trol_zoopl_uw[,2],Trol_zoopl_uw[,1], xlab = "PC (among_mod_conn)", ylab = "z-score (within_mod_degree)",xlim = c(-0.2,0.8), 
#     ylim= c(-3,3), main = "Zooplankton over time uw", col=factor(Trol_zoopl_uw$V5), cex=1.5, pch=16)
#text(Trol_zoopl_uw[,2],Trol_zoopl_uw[,1]-0.05, labels = Trol_zoopl_uw[,3], cex = 0.5)
#text(x=0.75, y=2.7, lab= "Hub netwrok connector", cex=0.8)
#text(x=0.75, y=2.2, lab= "module connectors", cex=0.8)
#text(x=-0.1, y=2.7, lab= "module hub", cex=0.8)
#text(x=-0.1, y=2.2, lab= "peripherals", cex=0.8)
#abline(v=0.6, h=2.5)
#legend("bottomright", legend =c("1981","1986","1991", "1996", "2001", "2006", "2011"), fill = c(1,2,3,4,5,6,7))
#Phytoplankton
Top_rol_1981[,4]<-c(sp.fg.1981[2], use.names = F)
Top_rol_1981[,5]<-c(1981)
Trol_uw1981<-filter(Top_rol_1981, grepl("Phytoplankton", V4))
Top_rol_1986[,4]<-c(sp.fg.1986[2], use.names = F)
Top_rol_1986[,5]<-c(1986)
Trol_uw1986<-filter(Top_rol_1986, grepl("Phytoplankton", V4))
Top_rol_1991[,4]<-c(sp.fg.1991[2], use.names = F)
Top_rol_1991[,5]<-c(1991)
Trol_uw1991<-filter(Top_rol_1991, grepl("Phytoplankton", V4))
Top_rol_1996[,4]<-c(sp.fg.1996[2], use.names = F)
Top_rol_1996[,5]<-c(1996)
Trol_uw1996<-filter(Top_rol_1996, grepl("Phytoplankton", V4))
Top_rol_2001[,4]<-c(sp.fg.2001[2], use.names = F)
Top_rol_2001[,5]<-c(2001)
Trol_uw2001<-filter(Top_rol_2001, grepl("Phytoplankton", V4))
Top_rol_2006[,4]<-c(sp.fg.2006[2], use.names = F)
Top_rol_2006[,5]<-c(2006)
Trol_uw2006<-filter(Top_rol_2006, grepl("Phytoplankton", V4))
Top_rol_2011[,4]<-c(sp.fg.2011[2], use.names = F)
Top_rol_2011[,5]<-c(2011)
Trol_uw2011<-filter(Top_rol_2011, grepl("Phytoplankton", V4))
Trol_phytopl_uw<-vctrs::vec_c(Trol_uw1981, Trol_uw1986, Trol_uw1991, Trol_uw1996, Trol_uw2001, Trol_uw2006, Trol_uw2011)
Trol_pla_uw<-combine(Trol_zoopl_uw, Trol_phytopl_uw)


#par(xpd = T, mar = par()$mar + c(0,0,0,10))
#plot(Trol_zoopl_uw[,5], Trol_zoopl_uw[,1], ylab = "z-score (within_mod_degree)", xlab= "year",
#     main = "Phyto- & zooplankton within module degree over time unweighted", xlim=c(1980, 2012),xaxt="n",
#     pch=as.numeric(as.factor(Trol_zoopl_uw$V3)), col=as.factor(Trol_zoopl_uw$V3),ylim = c(-2.7,2))
#axis(1, at=seq(1981, 2011, by=5))
#legend(2013, 0.6, legend = levels(as.factor(Trol_zoopl_uw$V3)), pch=c(1:14), cex = 0.8, col=c(1:14))

#plot(Trol_zoopl_uw[,5], Trol_zoopl_uw[,2], ylab = "PC (among_mod_conn)", xlab= "year",
#     main = "Phyto- & zooplankton among module degree over time unweighted", xlim=c(1980, 2012), xaxt="n",
#     pch=as.numeric(as.factor(Trol_zoopl_uw$V3)), col=as.factor(Trol_zoopl_uw$V3),ylim = c(-0.1,0.7))
#axis(1, at=seq(1981, 2011, by=5))
#legend(2013, 0.6, legend = levels(as.factor(Trol_zoopl_uw$V3)), pch=c(1:14), cex = 0.8, col=c(1:14))


Plauw1<-ggplot(data=Trol_pla_uw, aes(Trol_pla_uw[,5],Trol_pla_uw[,1],colour=Trol_pla_uw[,3]))+
  geom_line(size=1.1)+geom_point(size=2)+labs(title = "Phyto- & zooplankton within-module degree over time unweighted",tag = "a)")+xlab("Year")+ylab("z-score")+
  labs(color="Species")+ylim(-2.65,2.5)+xlim(1980,2012)+scale_x_continuous(breaks= c(1981,1986,1991,1996,2001,2006,2011))+
  theme(plot.title = element_text(size = 17, face = "bold"))+
  theme(text = element_text(size = 17, face="bold"))+
  scale_color_manual(values = c("red","blue","green","cyan","black","orange","firebrick2","darkcyan","purple","coral","springgreen4","slategray","brown","magenta"))+
  theme_minimal()
  
Plauw2<-ggplot(data=Trol_pla_uw, aes(Trol_pla_uw[,5],Trol_pla_uw[,2],colour=Trol_pla_uw[,3]))+
  geom_line(size=1.1)+geom_point(size=2)+labs(title = "Phyto- & zooplankton among-module degree over time unweighted", tag = "b)")+xlab("Year")+ylab("PC")+
  labs(color="Species")+ylim(-0.05,0.65)+xlim(1980,2012)+scale_x_continuous(breaks= c(1981,1986,1991,1996,2001,2006,2011))+
  theme(plot.title = element_text(size = 17, face = "bold"))+
  theme(text = element_text(size = 17, face="bold"))+
  scale_color_manual(values = c("red","blue","green","cyan","black","orange","firebrick2","darkcyan","purple","coral","springgreen4","slategray","brown","magenta"))+
  theme_minimal()


#Phytoplankton
#Top_rol_1981[,4]<-c(sp.fg.1981[2], use.names = F)
#Top_rol_1981[,5]<-c(1981)
#Trol_uw1981<-filter(Top_rol_1981, grepl("Phytoplankton", V4))
#Top_rol_1986[,4]<-c(sp.fg.1986[2], use.names = F)
#Top_rol_1986[,5]<-c(1986)
#Trol_uw1986<-filter(Top_rol_1986, grepl("Phytoplankton", V4))
#Top_rol_1991[,4]<-c(sp.fg.1991[2], use.names = F)
#Top_rol_1991[,5]<-c(1991)
#Trol_uw1991<-filter(Top_rol_1991, grepl("Phytoplankton", V4))
#Top_rol_1996[,4]<-c(sp.fg.1996[2], use.names = F)
#Top_rol_1996[,5]<-c(1996)
#Trol_uw1996<-filter(Top_rol_1996, grepl("Phytoplankton", V4))
#Top_rol_2001[,4]<-c(sp.fg.2001[2], use.names = F)
#Top_rol_2001[,5]<-c(2001)
#Trol_uw2001<-filter(Top_rol_2001, grepl("Phytoplankton", V4))
#Top_rol_2006[,4]<-c(sp.fg.2006[2], use.names = F)
#Top_rol_2006[,5]<-c(2006)
#Trol_uw2006<-filter(Top_rol_2006, grepl("Phytoplankton", V4))
#Top_rol_2011[,4]<-c(sp.fg.2011[2], use.names = F)
#Top_rol_2011[,5]<-c(2011)
#Trol_uw2011<-filter(Top_rol_2011, grepl("Phytoplankton", V4))
#Trol_zoopl_uw<-vctrs::vec_c(Trol_uw1981, Trol_uw1986, Trol_uw1991, Trol_uw1996, Trol_uw2001, Trol_uw2006, Trol_uw2011)
#plot(Trol_zoopl_uw[,2],Trol_zoopl_uw[,1], xlab = "PC (among_mod_conn)", ylab = "z-score (within_mod_degree)",xlim = c(-0.2,0.8), 
#     ylim= c(-3,3), main = "Phytoplankton over time uw", col=factor(Trol_zoopl_uw$V5), cex=1.5, pch=16)
#text(Trol_zoopl_uw[,2],Trol_zoopl_uw[,1]-0.05, labels = Trol_zoopl_uw[,3], cex = 0.5)
#text(x=0.75, y=2.7, lab= "Hub netwrok connector", cex=0.8)
#text(x=0.75, y=2.2, lab= "module connectors", cex=0.8)
#text(x=-0.1, y=2.7, lab= "module hub", cex=0.8)
#text(x=-0.1, y=2.2, lab= "peripherals", cex=0.8)
#abline(v=0.6, h=2.5)
#legend("bottomright", legend =c("1981","1986","1991", "1996", "2001", "2006", "2011"), fill = c(1,2,3,4,5,6,7))

#plot(Trol_zoopl_uw[,5], Trol_zoopl_uw[,1], ylab = "z-score (within_mod_degree)", xlab= "year",
#     main = "Phytoplankton within module degree over time uweighted", xlim=c(1980, 2012),xaxt="n",
#     pch=as.numeric(as.factor(Trol_zoopl_uw$V3)), col=as.factor(Trol_zoopl_w$V3),ylim = c(-2.7,2))
#axis(1, at=seq(1981, 2011, by=5))
#legend(2013, 0.6, legend = levels(as.factor(Trol_zoopl_w$V3)), pch=c(1:14), cex = 0.8, col=c(1:14))

#plot(Trol_zoopl_uw[,5], Trol_zoopl_uw[,2], ylab = "PC (among_mod_conn)", xlab= "year",
#     main = "Phytoplankton among module degree over time uweighted", xlim=c(1980, 2012), xaxt="n",
#     pch=as.numeric(as.factor(Trol_zoopl_uw$V3)), col=as.factor(Trol_zoopl_w$V3),ylim = c(-0.1,0.6))
#axis(1, at=seq(1981, 2011, by=5))
#legend(2013, 0.6, legend = levels(as.factor(Trol_zoopl_w$V3)), pch=c(1:14), cex = 0.8, col=c(1:14))


#Benthos
Top_rol_1981[,4]<-c(sp.fg.1981[2], use.names = F)
Top_rol_1981[,5]<-c(1981)
Trol_uw1981<-filter(Top_rol_1981, grepl("Benthos", V4))
Top_rol_1986[,4]<-c(sp.fg.1986[2], use.names = F)
Top_rol_1986[,5]<-c(1986)
Trol_uw1986<-filter(Top_rol_1986, grepl("Benthos", V4))
Top_rol_1991[,4]<-c(sp.fg.1991[2], use.names = F)
Top_rol_1991[,5]<-c(1991)
Trol_uw1991<-filter(Top_rol_1991, grepl("Benthos", V4))
Top_rol_1996[,4]<-c(sp.fg.1996[2], use.names = F)
Top_rol_1996[,5]<-c(1996)
Trol_uw1996<-filter(Top_rol_1996, grepl("Benthos", V4))
Top_rol_2001[,4]<-c(sp.fg.2001[2], use.names = F)
Top_rol_2001[,5]<-c(2001)
Trol_uw2001<-filter(Top_rol_2001, grepl("Benthos", V4))
Top_rol_2006[,4]<-c(sp.fg.2006[2], use.names = F)
Top_rol_2006[,5]<-c(2006)
Trol_uw2006<-filter(Top_rol_2006, grepl("Benthos", V4))
Top_rol_2011[,4]<-c(sp.fg.2011[2], use.names = F)
Top_rol_2011[,5]<-c(2011)
Trol_uw2011<-filter(Top_rol_2011, grepl("Benthos", V4))
Trol_benth_uw<-vctrs::vec_c(Trol_uw1981, Trol_uw1986, Trol_uw1991, Trol_uw1996, Trol_uw2001, Trol_uw2006, Trol_uw2011)
#plot(Trol_zoopl_uw[,2],Trol_zoopl_uw[,1], xlab = "PC (among_mod_conn)", ylab = "z-score (within_mod_degree)",xlim = c(-0.2,0.8), 
#     ylim= c(-3,3), main = "Benthos over time uw", col=factor(Trol_zoopl_uw$V5), cex=1.5, pch=16)
#text(Trol_zoopl_uw[,2],Trol_zoopl_uw[,1]-0.05, labels = Trol_zoopl_uw[,3], cex = 0.5)
#text(x=0.75, y=2.7, lab= "Hub netwrok connector", cex=0.8)
#text(x=0.75, y=2.2, lab= "module connectors", cex=0.8)
#text(x=-0.1, y=2.7, lab= "module hub", cex=0.8)
#text(x=-0.1, y=2.2, lab= "peripherals", cex=0.8)
#abline(v=0.6, h=2.5)
#legend("bottomright", legend =c("1981","1986","1991", "1996", "2001", "2006", "2011"), fill = c(1,2,3,4,5,6,7))


#plot(Trol_zoopl_uw[,5], Trol_zoopl_uw[,1], ylab = "z-score (within_mod_degree)", xlab= "year",
#     main = "Benthos within module degree over time unweighted", xlim=c(1980, 2012),xaxt="n",
#     pch=as.numeric(as.factor(Trol_zoopl_uw$V3)), col= as.factor(Trol_zoopl_uw$V3),ylim = c(-2.7,2))
#axis(1, at=seq(1981, 2011, by=5))
#legend(2013, 0.6, legend = levels(as.factor(Trol_zoopl_uw$V3)), pch=c(1:54), cex = 0.8, col=c(1:54))

#plot(Trol_zoopl_uw[,5], Trol_zoopl_uw[,2], ylab = "PC (among_mod_conn)", xlab= "year",
#     main = "Benthos among module degree over time unweighted", xlim=c(1980, 2012), xaxt="n",
#     pch=as.numeric(as.factor(Trol_zoopl_uw$V3)), col=as.factor(Trol_zoopl_uw$V3),ylim = c(-0.1,0.7))
#axis(1, at=seq(1981, 2011, by=5))
#legend(2013, 0.6, legend = levels(as.factor(Trol_zoopl_uw$V3)), pch=c(1:14), cex = 0.8, col=c(1:14))


Beuw1<-ggplot(data=Trol_benth_uw, aes(Trol_benth_uw[,5],Trol_benth_uw[,1],colour=Trol_benth_uw[,3]))+
  geom_line(size=1.1)+geom_point(size=2)+labs(title = "Benthos within-module degree over time unweighted", tag = "a)")+xlab("Year")+ylab("z-score")+
  labs(color="Species")+ylim(-2.65,2.5)+xlim(1980,2012)+scale_x_continuous(breaks= c(1981,1986,1991,1996,2001,2006,2011))+
  theme(plot.title = element_text(size = 17, face = "bold"))+
  theme(text = element_text(size = 17, face="bold"))+
  #scale_color_manual(values = c("red","blue","green","cyan","black","yellow","pink","aquamarine","purple","coral","lightgreen","slategray","brown","magenta"))
  scale_color_manual(values = c("red","blue","green","cyan","black","orange","firebrick2","darkcyan","purple","coral","springgreen4","slategray","brown","magenta"))+
  theme_minimal()

Beuw2<-ggplot(data=Trol_benth_uw, aes(Trol_benth_uw[,5],Trol_benth_uw[,2],colour=Trol_benth_uw[,3]))+
  geom_line(size=1.1)+geom_point(size=2)+labs(title = "Benthos among-module degree over time unweighted", tags= "b)")+xlab("Year")+ylab("PC")+
  labs(color="Species")+ylim(-0.05,0.65)+xlim(1980,2012)+scale_x_continuous(breaks= c(1981,1986,1991,1996,2001,2006,2011))+
  theme(plot.title = element_text(size = 17, face = "bold"))+
  theme(text = element_text(size = 17, face="bold"))+
  scale_color_manual(values = c("red","blue","green","cyan","black","orange","firebrick2","darkcyan","purple","coral","springgreen4","slategray","brown","magenta"))+
  theme_minimal()


#Fish
Top_rol_1981[,4]<-c(sp.fg.1981[2], use.names = F)
Top_rol_1981[,5]<-c(1981)
Trol_uw1981<-filter(Top_rol_1981, grepl("Fish", V4))
Top_rol_1986[,4]<-c(sp.fg.1986[2], use.names = F)
Top_rol_1986[,5]<-c(1986)
Trol_uw1986<-filter(Top_rol_1986, grepl("Fish", V4))
Top_rol_1991[,4]<-c(sp.fg.1991[2], use.names = F)
Top_rol_1991[,5]<-c(1991)
Trol_uw1991<-filter(Top_rol_1991, grepl("Fish", V4))
Top_rol_1996[,4]<-c(sp.fg.1996[2], use.names = F)
Top_rol_1996[,5]<-c(1996)
Trol_uw1996<-filter(Top_rol_1996, grepl("Fish", V4))
Top_rol_2001[,4]<-c(sp.fg.2001[2], use.names = F)
Top_rol_2001[,5]<-c(2001)
Trol_uw2001<-filter(Top_rol_2001, grepl("Fish", V4))
Top_rol_2006[,4]<-c(sp.fg.2006[2], use.names = F)
Top_rol_2006[,5]<-c(2006)
Trol_uw2006<-filter(Top_rol_2006, grepl("Fish", V4))
Top_rol_2011[,4]<-c(sp.fg.2011[2], use.names = F)
Top_rol_2011[,5]<-c(2011)
Trol_uw2011<-filter(Top_rol_2011, grepl("Fish", V4))
Trol_fish_uw<-vctrs::vec_c(Trol_uw1981, Trol_uw1986, Trol_uw1991, Trol_uw1996, Trol_uw2001, Trol_uw2006, Trol_uw2011)
#plot(Trol_zoopl_uw[,2],Trol_zoopl_uw[,1], xlab = "PC (among_mod_conn)", ylab = "z-score (within_mod_degree)",xlim = c(-0.2,0.8), 
#     ylim= c(-3,3), main = "Fish over time uw", col=factor(Trol_zoopl_uw$V5), cex=1.5, pch=16)
#text(Trol_zoopl_uw[,2],Trol_zoopl_uw[,1]-0.05, labels = Trol_zoopl_uw[,3], cex = 0.5)
#text(x=0.75, y=2.7, lab= "Hub netwrok connector", cex=0.8)
#text(x=0.75, y=2.2, lab= "module connectors", cex=0.8)
#text(x=-0.1, y=2.7, lab= "module hub", cex=0.8)
#text(x=-0.1, y=2.2, lab= "peripherals", cex=0.8)
#abline(v=0.6, h=2.5)
#legend("bottomright", legend =c("1981","1986","1991", "1996", "2001", "2006", "2011"), fill = c(1,2,3,4,5,6,7))

#plot(Trol_zoopl_uw[,5], Trol_zoopl_uw[,1], ylab = "z-score (within_mod_degree)", xlab= "year",
#     main = "Fish within module degree over time unweighted", xlim=c(1980, 2012),xaxt="n",
#     pch=as.numeric(as.factor(Trol_zoopl_uw$V3)), col=as.factor(Trol_zoopl_uw$V3),ylim = c(-2.7,2))
#axis(1, at=seq(1981, 2011, by=5))
#legend(2013, 0.6, legend = levels(as.factor(Trol_zoopl_w$V3)), pch=c(1:14), cex = 0.8, col=c(1:14))

#plot(Trol_zoopl_uw[,5], Trol_zoopl_uw[,2], ylab = "PC (among_mod_conn)", xlab= "year",
#     main = "Fish among module degree over time unweighted", xlim=c(1980, 2012), xaxt="n",
#     pch=as.numeric(as.factor(Trol_zoopl_uw$V3)), col=as.factor(Trol_zoopl_uw$V3),ylim = c(-0.1,0.7))
#axis(1, at=seq(1981, 2011, by=5))
#legend(2013, 0.6, legend = levels(as.factor(Trol_zoopl_uw$V3)), pch=c(1:14), cex = 0.8, col=c(1:14))


Fuw1<-ggplot(data=Trol_fish_uw, aes(Trol_fish_uw[,5],Trol_fish_uw[,1],colour=Trol_fish_uw[,3]))+
  geom_line(size=1.1)+geom_point(size=2)+labs(title = "Fish within-module degree over time unweighted", tags= "a)")+xlab("Year")+ylab("z-score")+
  labs(color="Species")+ylim(-2.65,2.5)+xlim(1980,2012)+scale_x_continuous(breaks= c(1981,1986,1991,1996,2001,2006,2011))+
  theme(plot.title = element_text(size = 17, face = "bold"))+
  theme(text = element_text(size = 17, face="bold"))+
  scale_color_manual(values = c("red","blue","green","cyan","black","orange","firebrick2","darkcyan","purple","coral","springgreen4","slategray","brown","magenta"))+
  theme_minimal()

Fuw2<-ggplot(data=Trol_fish_uw, aes(Trol_fish_uw[,5],Trol_fish_uw[,2],colour=Trol_fish_uw[,3]))+
  geom_line(size=1.1)+geom_point(size=2)+labs(title = "Fish among-module degree over time unweighted", tag = "b)")+xlab("Year")+ylab("PC")+
  labs(color="Species")+ylim(-0.05,0.65)+xlim(1980,2012)+scale_x_continuous(breaks= c(1981,1986,1991,1996,2001,2006,2011))+
  theme(plot.title = element_text(size = 17, face = "bold"))+
  theme(text = element_text(size = 17, face="bold"))+
  scale_color_manual(values = c("red","blue","green","cyan","black","orange","firebrick2","darkcyan","purple","coral","springgreen4","slategray","brown","magenta"))+
  theme_minimal()

#weighted zooplankton
Top_rol_1981_w[,4]<-c(sp.fg.1981[2], use.names = F)
Top_rol_1981_w[,5]<-c(1981)
Trol_w1981<-filter(Top_rol_1981_w, grepl("Zooplankton", V4))
Top_rol_1986_w[,4]<-c(sp.fg.1986[2], use.names = F)
Top_rol_1986_w[,5]<-c(1986)
Trol_w1986<-filter(Top_rol_1986_w, grepl("Zooplankton", V4))
Top_rol_1991_w[,4]<-c(sp.fg.1991[2], use.names = F)
Top_rol_1991_w[,5]<-c(1991)
Trol_w1991<-filter(Top_rol_1991_w, grepl("Zooplankton", V4))
Top_rol_1996_w[,4]<-c(sp.fg.1996[2], use.names = F)
Top_rol_1996_w[,5]<-c(1996)
Trol_w1996<-filter(Top_rol_1996_w, grepl("Zooplankton", V4))
Top_rol_2001_w[,4]<-c(sp.fg.2001[2], use.names = F)
Top_rol_2001_w[,5]<-c(2001)
Trol_w2001<-filter(Top_rol_2001_w, grepl("Zooplankton", V4))
Top_rol_2006_w[,4]<-c(sp.fg.2006[2], use.names = F)
Top_rol_2006_w[,5]<-c(2006)
Trol_w2006<-filter(Top_rol_2006_w, grepl("Zooplankton", V4))
Top_rol_2011_w[,4]<-c(sp.fg.2011[2], use.names = F)
Top_rol_2011_w[,5]<-c(2011)
Trol_w2011<-filter(Top_rol_2011_w, grepl("Zooplankton", V4))
Trol_zoopl_w<-vctrs::vec_c(Trol_w1981, Trol_w1986, Trol_w1991, Trol_w1996, Trol_w2001, Trol_w2006, Trol_w2011)
#plot(Trol_zoopl_w[,2],Trol_zoopl_w[,1], xlab = "PC (among_mod_conn)", ylab = "z-score (within_mod_degree)",xlim = c(-0.2,0.8), 
#     ylim= c(-3,3), main = "Zooplankton over time w", col=factor(Trol_zoopl_w$V5), cex=1.5, pch=16)
#text(Trol_zoopl_w[,2],Trol_zoopl_w[,1]-0.05, labels = Trol_zoopl_w[,3], cex = 0.5)
#text(x=0.75, y=2.7, lab= "Hub netwrok connector", cex=0.8)
#text(x=0.75, y=2.2, lab= "module connectors", cex=0.8)
#text(x=-0.1, y=2.7, lab= "module hub", cex=0.8)
#text(x=-0.1, y=2.2, lab= "peripherals", cex=0.8)
#abline(v=0.6, h=2.5)
#legend("bottomright", legend =c("1981","1986","1991", "1996", "2001", "2006", "2011"), fill = c(1,2,3,4,5,6,7))

#phytoplankton w
Top_rol_1981_w[,4]<-c(sp.fg.1981[2], use.names = F)
Top_rol_1981_w[,5]<-c(1981)
Trol_w1981<-filter(Top_rol_1981_w, grepl("Phytoplankton", V4))
Top_rol_1986_w[,4]<-c(sp.fg.1986[2], use.names = F)
Top_rol_1986_w[,5]<-c(1986)
Trol_w1986<-filter(Top_rol_1986_w, grepl("Phytoplankton", V4))
Top_rol_1991_w[,4]<-c(sp.fg.1991[2], use.names = F)
Top_rol_1991_w[,5]<-c(1991)
Trol_w1991<-filter(Top_rol_1991_w, grepl("Phytoplankton", V4))
Top_rol_1996_w[,4]<-c(sp.fg.1996[2], use.names = F)
Top_rol_1996_w[,5]<-c(1996)
Trol_w1996<-filter(Top_rol_1996_w, grepl("Phytoplankton", V4))
Top_rol_2001_w[,4]<-c(sp.fg.2001[2], use.names = F)
Top_rol_2001_w[,5]<-c(2001)
Trol_w2001<-filter(Top_rol_2001_w, grepl("Phytoplankton", V4))
Top_rol_2006_w[,4]<-c(sp.fg.2006[2], use.names = F)
Top_rol_2006_w[,5]<-c(2006)
Trol_w2006<-filter(Top_rol_2006_w, grepl("Phytoplankton", V4))
Top_rol_2011_w[,4]<-c(sp.fg.2011[2], use.names = F)
Top_rol_2011_w[,5]<-c(2011)
Trol_w2011<-filter(Top_rol_2011_w, grepl("Phytoplankton", V4))
Trol_phytopl_w<-vctrs::vec_c(Trol_w1981, Trol_w1986, Trol_w1991, Trol_w1996, Trol_w2001, Trol_w2006, Trol_w2011)
Trol_pla_w<-combine(Trol_zoopl_w, Trol_phytopl_w)

#plot(Trol_zoopl_w[,5], Trol_zoopl_w[,1], ylab = "z-score (within_mod_degree)", xlab= "year",
#    main = "Phyto- & Zooolankton within module degree over time weighted ", xlim=c(1980, 2012),xaxt="n",
#     pch=as.numeric(as.factor(Trol_zoopl_w$V3)), col=as.factor(Trol_zoopl_w$V3), ylim = c(-2.7,2))
#axis(1, at=seq(1981, 2011, by=5))
#legend(2013, 0.6, legend = levels(as.factor(Trol_zoopl_w$V3)), pch=c(1:14), cex = 0.8, col=c(1:14))

#plot(Trol_zoopl_w[,5], Trol_zoopl_w[,2], ylab = "PC (among_mod_conn)", xlab= "year",
#     main = "Phyto- & Zooplankton among module degree over time weighted", xlim=c(1980, 2012), xaxt="n",
#     pch=as.numeric(as.factor(Trol_zoopl_w$V3)), col=as.factor(Trol_zoopl_w$V3), ylim = c(-0.1,0.7))
#axis(1, at=seq(1981, 2011, by=5))
#legend(2013, 0.6, legend = levels(as.factor(Trol_zoopl_w$V3)), pch=c(1:14), cex = 0.8, col=c(1:14))


Plaw1<-ggplot(data=Trol_pla_w, aes(Trol_pla_w[,5],Trol_pla_w[,1],colour=Trol_pla_w[,3]))+
  geom_line(size=1.1)+geom_point(size=2)+labs(title = "Phyto- & Zooolankton within-module degree over time weighted", tag = "c)")+xlab("Year")+ylab("z-score")+
  labs(color="Species")+ylim(-2.65,2.5)+xlim(1980,2012)+scale_x_continuous(breaks= c(1981,1986,1991,1996,2001,2006,2011))+
  theme(plot.title = element_text(size = 17, face = "bold"))+
  theme(text = element_text(size = 17, face="bold"))+
  scale_color_manual(values = c("red","blue","green","cyan","black","orange","firebrick2","darkcyan","purple","coral","springgreen4","slategray","brown","magenta"))+
  theme_minimal()

Plaw2<-ggplot(data=Trol_pla_w, aes(Trol_pla_w[,5],Trol_pla_w[,2],colour=Trol_pla_w[,3]))+
  geom_line(size=1.1)+geom_point(size=2)+labs(title = "Phyto- & Zooolankton among-module degree over time weighted", tag = "d)")+xlab("Year")+ylab("PC")+
  labs(color="Species")+ylim(-0.05,0.65)+xlim(1980,2012)+scale_x_continuous(breaks= c(1981,1986,1991,1996,2001,2006,2011))+
  theme(plot.title = element_text(size = 17, face = "bold"))+
  theme(text = element_text(size = 17, face="bold"))+
  scale_color_manual(values = c("red","blue","green","cyan","black","orange","firebrick2","darkcyan","purple","coral","springgreen4","slategray","brown","magenta"))+
  theme_minimal()


#phytoplankton w
#Top_rol_1981_w[,4]<-c(sp.fg.1981[2], use.names = F)
#Top_rol_1981_w[,5]<-c(1981)
#Trol_w1981<-filter(Top_rol_1981_w, grepl("Phytoplankton", V4))
#Top_rol_1986_w[,4]<-c(sp.fg.1986[2], use.names = F)
#Top_rol_1986_w[,5]<-c(1986)
#Trol_w1986<-filter(Top_rol_1986_w, grepl("Phytoplankton", V4))
#Top_rol_1991_w[,4]<-c(sp.fg.1991[2], use.names = F)
#Top_rol_1991_w[,5]<-c(1991)
#Trol_w1991<-filter(Top_rol_1991_w, grepl("Phytoplankton", V4))
#Top_rol_1996_w[,4]<-c(sp.fg.1996[2], use.names = F)
#Top_rol_1996_w[,5]<-c(1996)
#Trol_w1996<-filter(Top_rol_1996_w, grepl("Phytoplankton", V4))
#Top_rol_2001_w[,4]<-c(sp.fg.2001[2], use.names = F)
#Top_rol_2001_w[,5]<-c(2001)
#Trol_w2001<-filter(Top_rol_2001_w, grepl("Phytoplankton", V4))
#Top_rol_2006_w[,4]<-c(sp.fg.2006[2], use.names = F)
#Top_rol_2006_w[,5]<-c(2006)
#Trol_w2006<-filter(Top_rol_2006_w, grepl("Phytoplankton", V4))
#Top_rol_2011_w[,4]<-c(sp.fg.2011[2], use.names = F)
#Top_rol_2011_w[,5]<-c(2011)
#Trol_w2011<-filter(Top_rol_2011_w, grepl("Phytoplankton", V4))
#Trol_zoopl_w<-vctrs::vec_c(Trol_w1981, Trol_w1986, Trol_w1991, Trol_w1996, Trol_w2001, Trol_w2006, Trol_w2011)
#plot(Trol_zoopl_w[,2],Trol_zoopl_w[,1], xlab = "PC (among_mod_conn)", ylab = "z-score (within_mod_degree)",xlim = c(-0.2,0.8), 
 #    ylim= c(-3,3), main = "Phytoplankton over time w", col=factor(Trol_zoopl_w$V5), cex=1.5, pch=16)
#text(Trol_zoopl_w[,2],Trol_zoopl_w[,1]-0.05, labels = Trol_zoopl_w[,3], cex = 0.5)
#text(x=0.75, y=2.7, lab= "Hub netwrok connector", cex=0.8)
#text(x=0.75, y=2.2, lab= "module connectors", cex=0.8)
#text(x=-0.1, y=2.7, lab= "module hub", cex=0.8)
#text(x=-0.1, y=2.2, lab= "peripherals", cex=0.8)
#abline(v=0.6, h=2.5)
#legend("bottomright", legend =c("1981","1986","1991", "1996", "2001", "2006", "2011"), fill = c(1,2,3,4,5,6,7))


#plot(Trol_zoopl_w[,5], Trol_zoopl_w[,1], ylab = "z-score (within_mod_degree)", xlab= "year",
#     main = "Phytolankton within module degree over time weighted ", xlim=c(1980, 2012),xaxt="n",
#     pch=as.numeric(as.factor(Trol_zoopl_w$V3)), col=as.factor(Trol_zoopl_w$V3),ylim = c(-2.7,2))
#axis(1, at=seq(1981, 2011, by=5))
#legend(2013, 0.6, legend = levels(as.factor(Trol_zoopl_w$V3)), pch=c(1:14), cex = 0.8, col=c(1:14))

#plot(Trol_zoopl_w[,5], Trol_zoopl_w[,2], ylab = "PC (among_mod_conn)", xlab= "year",
#     main = "Phytoplankton among module degree over time weighted", xlim=c(1980, 2012), ylim = c(-0.1,0.6), xaxt="n", 
#     pch=as.numeric(as.factor(Trol_zoopl_w$V3)), col=as.factor(Trol_zoopl_w$V3))
#axis(1, at=seq(1981, 2011, by=5))
#legend(2013, 0.1, legend = levels(as.factor(Trol_zoopl_w$V3)), pch=c(1:14), cex = 0.8, col=c(1:14))


#Benthos W
Top_rol_1981_w[,4]<-c(sp.fg.1981[2], use.names = F)
Top_rol_1981_w[,5]<-c(1981)
Trol_w1981<-filter(Top_rol_1981_w, grepl("Benthos", V4))
Top_rol_1986_w[,4]<-c(sp.fg.1986[2], use.names = F)
Top_rol_1986_w[,5]<-c(1986)
Trol_w1986<-filter(Top_rol_1986_w, grepl("Benthos", V4))
Top_rol_1991_w[,4]<-c(sp.fg.1991[2], use.names = F)
Top_rol_1991_w[,5]<-c(1991)
Trol_w1991<-filter(Top_rol_1991_w, grepl("Benthos", V4))
Top_rol_1996_w[,4]<-c(sp.fg.1996[2], use.names = F)
Top_rol_1996_w[,5]<-c(1996)
Trol_w1996<-filter(Top_rol_1996_w, grepl("Benthos", V4))
Top_rol_2001_w[,4]<-c(sp.fg.2001[2], use.names = F)
Top_rol_2001_w[,5]<-c(2001)
Trol_w2001<-filter(Top_rol_2001_w, grepl("Benthos", V4))
Top_rol_2006_w[,4]<-c(sp.fg.2006[2], use.names = F)
Top_rol_2006_w[,5]<-c(2006)
Trol_w2006<-filter(Top_rol_2006_w, grepl("Benthos", V4))
Top_rol_2011_w[,4]<-c(sp.fg.2011[2], use.names = F)
Top_rol_2011_w[,5]<-c(2011)
Trol_w2011<-filter(Top_rol_2011_w, grepl("Benthos", V4))
Trol_benth_w<-vctrs::vec_c(Trol_w1981, Trol_w1986, Trol_w1991, Trol_w1996, Trol_w2001, Trol_w2006, Trol_w2011)
#plot(Trol_zoopl_w[,2],Trol_zoopl_w[,1], xlab = "PC (among_mod_conn)", ylab = "z-score (within_mod_degree)",xlim = c(-0.2,0.8), 
 #    ylim= c(-3,3), main = "Benthos over time w", col=factor(Trol_zoopl_w$V5), cex=1.5, pch=16)
#text(Trol_zoopl_w[,2],Trol_zoopl_w[,1]-0.05, labels = Trol_zoopl_w[,3], cex = 0.5)
#text(x=0.75, y=2.7, lab= "Hub netwrok connector", cex=0.8)
#text(x=0.75, y=2.2, lab= "module connectors", cex=0.8)
#text(x=-0.1, y=2.7, lab= "module hub", cex=0.8)
#text(x=-0.1, y=2.2, lab= "peripherals", cex=0.8)
#abline(v=0.6, h=2.5)
#legend("bottomright", legend =c("1981","1986","1991", "1996", "2001", "2006", "2011"), fill = c(1,2,3,4,5,6,7))

#nya plots top roles x och y skillt över tid 

#plot(Trol_zoopl_w[,5], Trol_zoopl_w[,1], ylab = "z-score (within_mod_degree)", xlab= "year",
#     main = "Benthos within module degree over time weighted", xlim=c(1980, 2012),xaxt="n",
#     pch=as.numeric(as.factor(Trol_zoopl_w$V3)), col=as.factor(Trol_zoopl_w$V3), ylim = c(-2.7,2))
#axis(1, at=seq(1981, 2011, by=5))
#legend(2013, 0.6, legend = levels(as.factor(Trol_zoopl_w$V3)), pch=c(1:14), cex = 0.8, col=c(1:14))

#plot(Trol_zoopl_w[,5], Trol_zoopl_w[,2], ylab = "PC (among_mod_conn)", xlab= "year",
#     main = "Benthos among module degree over time weighted", xlim=c(1980, 2012), xaxt="n",
#     pch=as.numeric(as.factor(Trol_zoopl_w$V3)), col=as.factor(Trol_zoopl_w$V3), ylim = c(-0.1,0.7))
#axis(1, at=seq(1981, 2011, by=5))
#legend(2013, 0.6, legend = levels(as.factor(Trol_zoopl_w$V3)), pch=c(1:14), cex = 0.8, col=c(1:14))


Bew1<-ggplot(data=Trol_benth_w, aes(Trol_benth_w[,5],Trol_benth_w[,1],colour=Trol_benth_w[,3]))+
  geom_line(size=1.1)+geom_point(size=2)+labs(title = "Benthos within-module degree over time weighted", tag = "c)")+xlab("Year")+ylab("z-score")+
  labs(color="Species")+ylim(-2.65,2.5)+xlim(1980,2012)+scale_x_continuous(breaks= c(1981,1986,1991,1996,2001,2006,2011))+
  theme(plot.title = element_text(size = 17, face = "bold"))+
  theme(text = element_text(size = 17, face="bold"))+
  scale_color_manual(values = c("red","blue","green","cyan","black","orange","firebrick2","darkcyan","purple","coral","springgreen4","slategray","brown","magenta"))+
  theme_minimal()

Bew2<-ggplot(data=Trol_benth_w, aes(Trol_benth_w[,5],Trol_benth_w[,2],colour=Trol_benth_w[,3]))+
  geom_line(size=1.1)+geom_point(size=2)+labs(title = "Benthos among-module degree over time weighted", tag = "d)")+xlab("Year")+ylab("PC")+
  labs(color="Species")+ylim(-0.05,0.65)+xlim(1980,2012)+scale_x_continuous(breaks= c(1981,1986,1991,1996,2001,2006,2011))+
  theme(plot.title = element_text(size = 17, face = "bold"))+
  theme(text = element_text(size = 17, face="bold"))+
  scale_color_manual(values = c("red","blue","green","cyan","black","orange","firebrick2","darkcyan","purple","coral","springgreen4","slategray","brown","magenta"))+
  theme_minimal()


#Fish W
Top_rol_1981_w[,4]<-c(sp.fg.1981[2], use.names = F)
Top_rol_1981_w[,5]<-c(1981)
Trol_w1981<-filter(Top_rol_1981_w, grepl("Fish", V4))
Top_rol_1986_w[,4]<-c(sp.fg.1986[2], use.names = F)
Top_rol_1986_w[,5]<-c(1986)
Trol_w1986<-filter(Top_rol_1986_w, grepl("Fish", V4))
Top_rol_1991_w[,4]<-c(sp.fg.1991[2], use.names = F)
Top_rol_1991_w[,5]<-c(1991)
Trol_w1991<-filter(Top_rol_1991_w, grepl("Fish", V4))
Top_rol_1996_w[,4]<-c(sp.fg.1996[2], use.names = F)
Top_rol_1996_w[,5]<-c(1996)
Trol_w1996<-filter(Top_rol_1996_w, grepl("Fish", V4))
Top_rol_2001_w[,4]<-c(sp.fg.2001[2], use.names = F)
Top_rol_2001_w[,5]<-c(2001)
Trol_w2001<-filter(Top_rol_2001_w, grepl("Fish", V4))
Top_rol_2006_w[,4]<-c(sp.fg.2006[2], use.names = F)
Top_rol_2006_w[,5]<-c(2006)
Trol_w2006<-filter(Top_rol_2006_w, grepl("Fish", V4))
Top_rol_2011_w[,4]<-c(sp.fg.2011[2], use.names = F)
Top_rol_2011_w[,5]<-c(2011)
Trol_w2011<-filter(Top_rol_2011_w, grepl("Fish", V4))
Trol_fish_w<-vctrs::vec_c(Trol_w1981, Trol_w1986, Trol_w1991, Trol_w1996, Trol_w2001, Trol_w2006, Trol_w2011)
#plot(Trol_zoopl_w[,2],Trol_zoopl_w[,1],xlab = "PC (among_mod_conn)", ylab = "z-score (within_mod_degree)",xlim = c(-0.2,0.8), 
    # ylim= c(-3,3), main = "Fish over time w", col=factor(Trol_zoopl_w$V5), cex=1.5, pch=16)
 # text(Trol_zoopl_w[,2],Trol_zoopl_w[,1]-0.05, labels = Trol_zoopl_w[,3], cex = 0.5)
#text(x=0.75, y=2.7, lab= "Hub netwrok connector", cex=0.8)
#text(x=0.75, y=2.2, lab= "module connectors", cex=0.8)
#text(x=-0.1, y=2.7, lab= "module hub", cex=0.8)
#text(x=-0.1, y=2.2, lab= "peripherals", cex=0.8)
#abline(v=0.6, h=2.5)
#legend("bottomright", legend =c("1981","1986","1991", "1996", "2001", "2006", "2011"), fill= c(1,2,3,4,5,6,7))

#filter(Trol_zoopl_w,among_module_conn>0.6)

#nya plots top roles x och y skillt över tid 

#plot(Trol_zoopl_w[,5], Trol_zoopl_w[,1], ylab = "z-score (within_mod_degree)", xlab= "year",
#     main = "Fish within module degree over time weighted", xlim=c(1980, 2012),xaxt="n",
#     pch=as.numeric(as.factor(Trol_zoopl_w$V3)),col=as.factor(Trol_zoopl_w$V3), ylim = c(-2.7,2))
#axis(1, at=seq(1981, 2011, by=5))
#text(Trol_zoopl_w[,5]+1, Trol_zoopl_w[,1]-0.01, labels = Trol_zoopl_w[,3], cex = 0.6)
#legend(2013, 1.8, legend = levels(as.factor(Trol_zoopl_w$V3)), pch=c(1:14), cex = 0.8, col=c(1:14))

#plot(Trol_zoopl_w[,5], Trol_zoopl_w[,2], ylab = "PC (among_mod_conn)", xlab= "year",
#     main = "Fish among module degree over time weighted", xlim=c(1980, 2012), xaxt="n",
#     pch=as.numeric(as.factor(Trol_zoopl_w$V3)), col=as.factor(Trol_zoopl_w$V3), ylim = c(-0.1,0.7))
#axis(1, at=seq(1981, 2011, by=5))
#text(Trol_zoopl_w[,5]+1.2, Trol_zoopl_w[,2], labels = Trol_zoopl_w[,3], cex = 0.6)
#legend(2013, 0.6, legend = levels(as.factor(Trol_zoopl_w$V3)), pch=c(1:14), cex = 0.8, col=c(1:14))


Fw1<-ggplot(data=Trol_fish_w, aes(Trol_fish_w[,5],Trol_fish_w[,1],colour=Trol_fish_w[,3]))+
  geom_line(size=0.1)+geom_point(size=3)+labs(title = "Fish within-module degree over time weighted", tag = "c)")+xlab("Year")+ylab("z-score")+
  labs(color="Species")+ylim(-2.65,2.5)+xlim(1980,2012)+scale_x_continuous(breaks= c(1981,1986,1991,1996,2001,2006,2011))+
  theme(plot.title = element_text(size = 17, face = "bold"))+
  theme(text = element_text(size = 17, face="bold"))+
  scale_color_manual(values = c("red","blue","green","cyan","black","orange","firebrick2","darkcyan","purple","coral","springgreen4","slategray","brown","magenta"))+
  theme_minimal()

Fw2<-ggplot(data=Trol_fish_w, aes(Trol_fish_w[,5],Trol_fish_w[,2],colour=Trol_fish_w[,3]))+
  geom_line(size=0.1)+geom_point(size=3)+xlab("Year")+ylab("PC")+
  labs(color="Species")+labs(title = "Weighted Fish among-module degree over time", tag = "d)")+ylim(-0.05,0.65)+xlim(1980,2012)+scale_x_continuous(breaks= c(1981,1986,1991,1996,2001,2006,2011))+
  theme(plot.title = element_text(size = 17, face = "bold"))+
  theme(text = element_text(size = 17, face="bold"))+
  scale_color_manual(values = c("red","blue","green","cyan","black","orange","firebrick2","darkcyan","purple","coral","springgreen4","slategray","brown","magenta"))+
  theme_minimal()


############## nya plots, skippar linje vid år arten inte finns med
write.csv2(Trol_fish_w, "fishroles.csv", row.names = F)# adda manulet till året om arten inte fanns 
write.csv2(Trol_benth_w, "bentroles.csv", row.names = F)
write.csv2(Trol_pla_w, "planktroles.csv", row.names = F)

fishroles.csv<-read.csv2("fishroles.csv")
new.w.fish1<-ggplot(data=fishroles.csv, aes(fishroles.csv[,5],fishroles.csv[,1],colour=fishroles.csv[,3]))+
  geom_line(size=0.1)+geom_point(size=3,alpha=0.6)+labs(title = "Weighted Fish within-module degree over time", tag = "c)")+xlab("Year")+ylab("z-score")+
  labs(color="Species")+ylim(-2.65,2.5)+xlim(1980,2012)+scale_x_continuous(breaks= c(1981,1986,1991,1996,2001,2006,2011))+
  theme(plot.title = element_text(size = 18, face = "bold"))+
  theme(text = element_text(size = 18, face="bold"))+
  scale_color_manual(values = c("red","blue","green","cyan","black","orange","firebrick2","darkcyan","purple","coral","springgreen4","slategray","brown","magenta"))+
  theme_minimal()
new.w.fish2<-ggplot(data=fishroles.csv, aes(fishroles.csv[,5],fishroles.csv[,2],colour=fishroles.csv[,3]))+
  geom_line(size=0.1)+geom_point(size=3,alpha=0.6)+labs(title = "Weighted Fish among-module degree over time", tag = "d)")+xlab("Year")+ylab("PC-score")+
  labs(color="Species")+ylim(-0.05,0.65)+xlim(1980,2012)+scale_x_continuous(breaks= c(1981,1986,1991,1996,2001,2006,2011))+
  theme(plot.title = element_text(size = 18, face = "bold"))+
  theme(text = element_text(size = 18, face="bold"))+
  scale_color_manual(values = c("red","blue","green","cyan","black","orange","firebrick2","darkcyan","purple","coral","springgreen4","slategray","brown","magenta"))+
  theme_minimal()


bentroles.csv<-read.csv2("bentroles.csv")
new.w.bent1<-ggplot(data=bentroles.csv, aes(bentroles.csv[,5],bentroles.csv[,1],colour=bentroles.csv[,3]))+
  geom_line(size=0.5)+geom_point(size=3,alpha=0.6)+labs(title = "Weighted benthos within-module degree over time weighted", tag = "c)")+xlab("Year")+ylab("z-score")+
  labs(color="Species")+ylim(-2.65,2.5)+xlim(1980,2012)+scale_x_continuous(breaks= c(1981,1986,1991,1996,2001,2006,2011))+
  theme(plot.title = element_text(size = 18, face = "bold"))+
  theme(text = element_text(size = 18, face="bold"))+
  scale_color_manual(values = c("red","blue","green","cyan","black","orange","firebrick2","darkcyan","purple","coral","springgreen4","slategray","brown","magenta"))+
  theme_minimal()
new.w.bent2<-ggplot(data=bentroles.csv, aes(bentroles.csv[,5],bentroles.csv[,2],colour=bentroles.csv[,3]))+
  geom_line(size=0.5)+geom_point(size=3,alpha=0.6)+labs(title = "Weighted benthos among-module degree over time", tag = "d)")+xlab("Year")+ylab("PC-score")+
  labs(color="Species")+ylim(-0.05,0.65)+xlim(1980,2012)+scale_x_continuous(breaks= c(1981,1986,1991,1996,2001,2006,2011))+
  theme(plot.title = element_text(size = 18, face = "bold"))+
  theme(text = element_text(size = 18, face="bold"))+
  scale_color_manual(values = c("red","blue","green","cyan","black","orange","firebrick2","darkcyan","purple","coral","springgreen4","slategray","brown","magenta"))+
  theme_minimal()


planktroles.csv<-read.csv2("planktroles.csv")
new.w.plank1<-ggplot(data=planktroles.csv, aes(planktroles.csv[,5],planktroles.csv[,1],colour=planktroles.csv[,3]))+
  geom_line(size=0.5)+geom_point(size=3,alpha=0.6)+labs(title = "Weighted Plankton within-module degree over time weighted", tag = "c)")+xlab("Year")+ylab("z-score")+
  labs(color="Species")+ylim(-2.65,2.5)+xlim(1980,2012)+scale_x_continuous(breaks= c(1981,1986,1991,1996,2001,2006,2011))+
  theme(plot.title = element_text(size = 18, face = "bold"))+
  theme(text = element_text(size = 18, face="bold"))+
  scale_color_manual(values = c("red","blue","green","cyan","black","orange","firebrick2","darkcyan","purple","coral","springgreen4","slategray","brown","magenta"))+
  theme_minimal()
new.w.plank2<-ggplot(data=planktroles.csv, aes(planktroles.csv[,5],planktroles.csv[,2],colour=planktroles.csv[,3]))+
  geom_line(size=0.5)+geom_point(size=3,alpha=0.6)+labs(title = "Weighted Plankton among-module degree over time", tag = "d)")+xlab("Year")+ylab("PC-score")+
  labs(color="Species")+ylim(-0.05,0.65)+xlim(1980,2012)+scale_x_continuous(breaks= c(1981,1986,1991,1996,2001,2006,2011))+
  theme(plot.title = element_text(size = 18, face = "bold"))+
  theme(text = element_text(size = 18, face="bold"))+
  scale_color_manual(values = c("red","blue","green","cyan","black","orange","firebrick2","darkcyan","purple","coral","springgreen4","slategray","brown","magenta"))+
  theme_minimal()

write.csv2(Trol_fish_uw, "fishroles_uw.csv", row.names = F)# adda manuelt till året som fisken inte fanns 
write.csv2(Trol_benth_uw, "bentroles_uw.csv", row.names = F)
write.csv2(Trol_pla_uw, "planktroles_uw.csv", row.names = F)

uwfishroles.csv<-read.csv2("fishroles_uw.csv")
uwbentroles.csv<-read.csv2("bentroles_uw.csv")
uwplanktroles.csv<-read.csv2("planktroles_uw.csv")

new.uw.bent1<-ggplot(data=uwbentroles.csv, aes(uwbentroles.csv[,5],uwbentroles.csv[,1],colour=uwbentroles.csv[,3]))+
  geom_line(size=0.5)+geom_point(size=3,alpha=0.6)+labs(title = "Benthos within-module degree over time weighted", tag = "a)")+xlab("Year")+ylab("z-score")+
  labs(color="Species")+ylim(-2.65,2.5)+xlim(1980,2012)+scale_x_continuous(breaks= c(1981,1986,1991,1996,2001,2006,2011))+
  theme(plot.title = element_text(size = 18, face = "bold"))+
  theme(text = element_text(size = 18, face="bold"))+
  scale_color_manual(values = c("red","blue","green","cyan","black","orange","firebrick2","darkcyan","purple","coral","springgreen4","slategray","brown","magenta"))+
  theme_minimal()
new.uw.bent2<-ggplot(data=uwbentroles.csv, aes(uwbentroles.csv[,5],uwbentroles.csv[,2],colour=uwbentroles.csv[,3]))+
  geom_line(size=0.5)+geom_point(size=3,alpha=0.6)+labs(title = "Benthos among-module degree over time", tag = "b)")+xlab("Year")+ylab("PC-score")+
  labs(color="Species")+ylim(-0.05,0.65)+xlim(1980,2012)+scale_x_continuous(breaks= c(1981,1986,1991,1996,2001,2006,2011))+
  theme(plot.title = element_text(size = 18, face = "bold"))+
  theme(text = element_text(size = 18, face="bold"))+
  scale_color_manual(values = c("red","blue","green","cyan","black","orange","firebrick2","darkcyan","purple","coral","springgreen4","slategray","brown","magenta"))+
  theme_minimal()

new.uw.plank1<-ggplot(data=uwplanktroles.csv, aes(uwplanktroles.csv[,5],uwplanktroles.csv[,1],colour=uwplanktroles.csv[,3]))+
  geom_line(size=0.5)+geom_point(size=3,alpha=0.6)+labs(title = "Plankton within-module degree over time", tag = "a)")+xlab("Year")+ylab("z-score")+
  labs(color="Species")+ylim(-2.65,2.5)+xlim(1980,2012)+scale_x_continuous(breaks= c(1981,1986,1991,1996,2001,2006,2011))+
  theme(plot.title = element_text(size = 18, face = "bold"))+
  theme(text = element_text(size = 18, face="bold"))+
  scale_color_manual(values = c("red","blue","green","cyan","black","orange","firebrick2","darkcyan","purple","coral","springgreen4","slategray","brown","magenta"))+
  theme_minimal()
new.uw.plank2<-ggplot(data=uwplanktroles.csv, aes(uwplanktroles.csv[,5],uwplanktroles.csv[,2],colour=uwplanktroles.csv[,3]))+
  geom_line(size=0.5)+geom_point(size=3,alpha=0.6)+labs(title = "Plankton among-module degree over time", tag = "b)")+xlab("Year")+ylab("PC-score")+
  labs(color="Species")+ylim(-0.05,0.65)+xlim(1980,2012)+scale_x_continuous(breaks= c(1981,1986,1991,1996,2001,2006,2011))+
  theme(plot.title = element_text(size = 18, face = "bold"))+
  theme(text = element_text(size = 18, face="bold"))+
  scale_color_manual(values = c("red","blue","green","cyan","black","orange","firebrick2","darkcyan","purple","coral","springgreen4","slategray","brown","magenta"))+
  theme_minimal()

new.uw.fish1<-ggplot(data=uwfishroles.csv, aes(uwfishroles.csv[,5],uwfishroles.csv[,1],colour=uwfishroles.csv[,3]))+
  geom_line(size=0.1)+geom_point(size=3,alpha=0.6)+labs(title = "Fish within-module degree over time", tag = "a)")+xlab("Year")+ylab("z-score")+
  labs(color="Species")+ylim(-2.65,2.5)+xlim(1980,2012)+scale_x_continuous(breaks= c(1981,1986,1991,1996,2001,2006,2011))+
  theme(plot.title = element_text(size = 18, face = "bold"))+
  theme(text = element_text(size = 18, face="bold"))+
  scale_color_manual(values = c("red","blue","green","cyan","black","orange","firebrick2","darkcyan","purple","coral","springgreen4","slategray","brown","magenta"))+
  theme_minimal()
new.uw.fish2<-ggplot(data=uwfishroles.csv, aes(uwfishroles.csv[,5],uwfishroles.csv[,2],colour=uwfishroles.csv[,3]))+
  geom_line(size=0.1)+geom_point(size=3,alpha=0.6)+labs(title = "Fish among-module degree over time", tag = "b)")+xlab("Year")+ylab("PC-score")+
  labs(color="Species")+ylim(-0.05,0.65)+xlim(1980,2012)+scale_x_continuous(breaks= c(1981,1986,1991,1996,2001,2006,2011))+
  theme(plot.title = element_text(size = 18, face = "bold"))+
  theme(text = element_text(size = 18, face="bold"))+
  scale_color_manual(values = c("red","blue","green","cyan","black","orange","firebrick2","darkcyan","purple","coral","springgreen4","slategray","brown","magenta"))+
  theme_minimal()

new.uw.fish1+theme(legend.position = "none")+new.uw.fish2+theme(legend.position = "none")+new.w.fish1+theme(legend.position = "none")+new.w.fish2

#new.uw.fish1+theme(legend.position = "none")+new.w.fish1+theme(legend.position = "none")+new.uw.fish2+theme(legend.position = "none")+new.w.fish2+
#  theme(legend.text= element_text(face = "bold", size = 12))+theme(legend.title = element_text(face = "bold"))

new.uw.fish1+theme(legend.position = "none")+new.uw.fish2+theme(legend.position = "none")+new.w.fish1+theme(legend.position = "none")+new.w.fish2+
  theme(legend.text= element_text(face = "bold", size = 12))+theme(legend.title = element_text(face = "bold"))

new.uw.plank1+theme(legend.position = "none")+
  new.uw.plank2+theme(legend.position = "none")+new.w.plank1+theme(legend.position = "none")+
  new.w.plank2+theme(legend.text= element_text(face = "bold", size = 12))+theme(legend.title = element_text(face = "bold"))

new.uw.bent1+theme(legend.position = "none")+new.uw.bent2+theme(legend.position = "none")+
  new.w.bent1+theme(legend.position = "none")+new.w.bent2+theme(legend.text= element_text(face = "bold", size = 12))+theme(legend.title = element_text(face = "bold"))


#ggplot(data=Trol_fish_w)+
#  geom_line(size=1, aes(Trol_fish_w[,5],Trol_fish_w[,2],colour=Trol_fish_w[,3]))+
#  geom_point(aes(Trol_fish_w[,5],Trol_fish_w[,2],shape=Trol_fish_w[,3], colour=Trol_fish_w[,3]), size=4)+
#  scale_shape_manual(values=seq(0,14))+xlab("Year")+ylab("PC")+labs(color="Species")+labs(title = "Fish among module degree over time weighted", tag = "d)")+ylim(-0.05,0.65)+xlim(1980,2012)+scale_x_continuous(breaks= c(1981,1986,1991,1996,2001,2006,2011))+
#  theme(plot.title = element_text(size = 19, face = "bold"))+
#  theme(text = element_text(size = 17, face="bold"))+
#  scale_color_manual(values = c("red","blue","green","cyan","black","yellow","pink","aquamarine","purple","coral","lightgreen","slategray","brown","magenta"))



#+scale_color_manual(values = c("red","blue","green","cyan","black","yellow","pink","aquamarine","purple","coral","lightgreen","slategray","brown","magenta"))

library(patchwork)

Plauw1+theme(legend.position = "none")+Plaw1+theme(legend.position = "none")+Plauw2+theme(legend.position = "none")+Plaw2

Beuw1+theme(legend.position = "none")+Bew1+theme(legend.position = "none")+Beuw2+theme(legend.position = "none")+Bew2

Fuw1+theme(legend.position = "none")+Fw1+theme(legend.position = "none")+Fuw2+theme(legend.position = "none")+Fw2


#########Barplot of the number of species within each module###########
alltest2<-read.csv("alltest2.csv", sep = ";", header = T)
Alltest3<-read.csv("alltest3.csv", sep = ";", header = T)
correct_metamoduw<-alltest2[,3]
correct_metamodw<-Alltest3[,3]
meta.sp.mod<-table(info$fg, correct_metamoduw)
barplot(meta.sp.mod, col = colFG, xlab="", ylab = "no. of species", ylim = c(0,34), main = "Metaweb functional groups per module")
meta.sp.mod.w<-table(info$fg, correct_metamodw)
barplot(meta.sp.mod.w, col = colFG, xlab="", ylab = "no. of species", ylim = c(0,34), main = "Weighted Metaweb functional groups per module")


sp.mod.1981<-table(sp.fg.1981[,2], mod.aff_UW1981)
barplot(sp.mod.1981, col = colFG, xlab="", ylab = "no. of species", ylim = c(0,34), main = "1981 functional groups per module")
sp.mod.1981w<-table(sp.fg.1981[,2], mod.aff_W1981)
barplot(sp.mod.1981w, col = colFG, xlab="", ylab = "no. of species", ylim = c(0,34), main = "Weighted 1981 functional groups per module")


sp.mod.1986<-table(sp.fg.1986[,2], mod.aff_UW1986)
barplot(sp.mod.1986, col = colFG, xlab="", ylab = "no. of species", ylim = c(0,34), main = "1986 functional groups per module")
sp.mod.1986w<-table(sp.fg.1986[,2], mod.aff_W1986)
barplot(sp.mod.1986w, col = colFG, xlab="", ylab = "no. of species", ylim = c(0,34), main = "Weighted 1986 functional groups per module")


sp.mod.1991<-table(sp.fg.1991[,2], mod.aff_UW1991)
barplot(sp.mod.1991, col = colFG, xlab="", ylab = "no. of species", ylim = c(0,34), main = "1991 functional groups per module")
sp.mod.1991w<-table(sp.fg.1991[,2], mod.aff_W1991)
barplot(sp.mod.1991w, col = colFG, xlab="", ylab = "no. of species", ylim = c(0,34), main = "Weighted 1991 functional groups per module")


sp.mod.1996<-table(sp.fg.1996[,2], mod.aff_UW1996)
barplot(sp.mod.1996, col = colFG, xlab="", ylab = "no. of species", ylim = c(0,34), main = "1996 functional groups per module")
sp.mod.1996w<-table(sp.fg.1996[,2], mod.aff_W1996)
barplot(sp.mod.1996w, col = colFG, xlab="", ylab = "no. of species", ylim = c(0,34), main = "Weighted 1996 functional groups per module")


sp.mod.2001<-table(sp.fg.2001[,2], mod.aff_UW2001)
barplot(sp.mod.2001, col = colFG, xlab="", ylab = "no. of species", ylim = c(0,34), main = "2001 functional groups per module")
sp.mod.2001w<-table(sp.fg.2001[,2], mod.aff_W2001)
barplot(sp.mod.2001w, col = colFG, xlab="", ylab = "no. of species", ylim = c(0,34), main = "Weighted 2001 functional groups per module")


sp.mod.2006<-table(sp.fg.2006[,2], mod.aff_UW2006)
barplot(sp.mod.2006, col = colFG, xlab="", ylab = "no. of species", ylim = c(0,34), main = "2006 functional groups per module")
sp.mod.2006w<-table(sp.fg.2006[,2], mod.aff_W2006)
barplot(sp.mod.2006w, col = colFG, xlab="", ylab = "no. of species", ylim = c(0,34), main = "Weighted 2006 functional groups per module")


sp.mod.2011<-table(sp.fg.2011[,2], mod.aff_UW2011)
barplot(sp.mod.2011, col = colFG, xlab="", ylab = "no. of species", ylim = c(0,34), main = "2011 functional groups per module")
sp.mod.2011w<-table(sp.fg.2011[,2], mod.aff_W2011)
barplot(sp.mod.2011w, col = colFG, xlab="", ylab = "no. of species", ylim = c(0,34), main = "Weighted 2011 functional groups per module")

der<-table(sp.fg.1981[,2], rep(c(1981),each=length(sp.fg.1981[,2])))
barplot(der,col = colFG)
der2<-table(sp.fg.1996[,2], rep(c(1),each=length(sp.fg.1996[,2])))
barplot(der2, der2,col = colFG)



#Alluvial layer diagram av modules over time
#??????????
#devtools::install_github("mbojan/alluvial")
library(alluvial)

Alluvtest<-read.csv("Alluvial test.csv", sep=";", stringsAsFactors = F, header = T)
Alluvtest[Alluvtest==0]<-"Absent"
Alluvtest[Alluvtest==1]<-"Present"
#alluvial(Alluvtest[,1:10],freq =1 , col = info$colfg, axis_labels = c("Taxa","Funktion group","Metaweb","1981","1986","1991","1996","2001","2006",2011))
alluvial(Alluvtest[,2:10],freq = Alluvtest[,11], cex = 0.5,col = info$colfg, mar = c(4,3,3,3), gap.width = 0.3,
         alpha = 0.7, cw=0.15, xw= 0.1, ann = T, cex.axis = 0.5)

alltest2<-read.csv("alltest2.csv", sep = ";", header = T)
#alltest2[alltest2==0]<-"not present"
alltest2[alltest2==1]<-"Module 1"
alltest2[alltest2==2]<-"Module 2"
alluvial(alltest2[3:10], freq = 1, cex = 0.8, col = info$colfg, title = "Unweighted module membership over time", mar=c(4,3,3,3), gap.width = 0.3,
         alpha = 0.7, cw=0.19, xw= 0.1, ann = T, cex.axis = 0.8,
         axis_labels = c("Metaweb","1981","1986","1991","1996","2001","2006",2011))


alluvial(alltest2[2:10], freq = 1, cex = 0.5, col = c(1:34), title = "Unweighted module membership over time", mar = c(6,3,3,9), gap.width = 0.3,
         alpha = 0.7, cw=0.19, xw= 0.1, ann = T, cex.axis = 0.5,
         axis_labels = c("Funktion group","Metaweb","1981","1986","1991","1996","2001","2006",2011))
legend(9.3,0.759,legend = sp.ab[,2], cex = 0.5, fill = c(1:34), bty="n")


#weighted networks
Alluvtest3<-read.csv("alltest3.csv", sep=";", header = T)
#Alluvtest3[Alluvtest3==0]<-"not present"
Alluvtest3[Alluvtest3==1]<-"Module 1"
Alluvtest3[Alluvtest3==2]<-"Module 2"
Alluvtest3[Alluvtest3==3]<-"Module 3"
alluvial(Alluvtest3[3:10], freq = 1, cex = 0.8, col = info$colfg, title = "Weighted module membership over time", mar = c(4,3,3,3), gap.width = 0.3,
         alpha = 0.7, cw=0.19, xw= 0.1, ann = T, cex.axis = 0.8,
         axis_labels = c("Metaweb","1981","1986","1991","1996","2001","2006",2011))


alluvial(Alluvtest3[4:10], freq = 1, cex = 0.5, col = c(1:34), title = "Weighted module membership over time", mar = c(6,3,3,9), gap.width = 0.3,
         alpha = 0.7, cw=0.19, xw= 0.1, ann = T, cex.axis = 0.5)
legend(9.3,0.759,legend = sp.ab[,2], cex = 0.5, fill = c(1:34), bty="n")

# this one works!!!
alltest2<-read.csv("alltest2.csv", sep = ";", header = T)
#alltest2[alltest2==0]<-"not present"
alltest2[alltest2==1]<-"Module 2"
alltest2[alltest2==2]<-"Module 1"
fixed.all(alltest2[3:10], freq = 1, cex = 0.8, col = info$colfg, gap.width = 0.3,
         alpha = 0.7, cw=0.19, xw= 0.1, cex.axis = 0.8,
         axis_labels = c("Metaweb","1981","1986","1991","1996","2001","2006",2011))

mtext("Unweighted module membership over time",3, line=2.5, font=2, cex = 2)

Alluvtest3<-read.csv("alltest3.csv", sep=";", header = T)
#Alluvtest3[Alluvtest3==0]<-"not present"
Alluvtest3[Alluvtest3==1]<-"Module 1"
Alluvtest3[Alluvtest3==2]<-"Module 2"
Alluvtest3[Alluvtest3==3]<-"Module 3"

fixed.all(Alluvtest3[3:10], freq = 1, cex = 0.8, col = info$colfg, gap.width = 0.3,
         alpha = 0.7, cw=0.19, xw= 0.1, cex.axis = 0.8,
         axis_labels = c("Metaweb","1981","1986","1991","1996","2001","2006",2011))

mtext("Weighted module membership over time",3, line=2.5, font=2, cex = 2)


###
#diagram instead of table with struckture
met.table<-read.csv("Stand_met.csv", sep=";", check.names = F,header = TRUE)

sp<-ggplot(data=met.table, aes(met.table[,1],met.table[,2]))+
  geom_line(size=1.1)+geom_point(size=2)+labs(title = "No. species", tag = "a)")+xlab("Year")+ylab("No. species")+
  ylim(21,33)+xlim(1980,2012)+scale_x_continuous(breaks= c(1981,1986,1991,1996,2001,2006,2011))+
  theme(plot.title = element_text(size = 15, face = "bold"))+
  theme(text = element_text(size = 15, face="bold"))+
  theme(axis.text.x = element_text(angle = 90))

link<-ggplot(data=met.table, aes(met.table[,1],met.table[,3]))+
  geom_line(size=1.1)+geom_point(size=2)+labs(title = "Links", tag = "b)")+xlab("Year")+ylab("No. links")+
  ylim(110,179)+xlim(1980,2012)+scale_x_continuous(breaks= c(1981,1986,1991,1996,2001,2006,2011))+
  theme(plot.title = element_text(size = 15, face = "bold"))+
  theme(text = element_text(size = 15, face="bold"))+
  theme(axis.text.x = element_text(angle = 90))

con<-ggplot(data=met.table, aes(met.table[,1],met.table[,4]*100))+
  geom_line(size=1.1)+geom_point(size=2)+labs(title = "Connectance", tag = "c)")+xlab("Year")+ylab("Connectance %")+
  ylim(17,20)+xlim(1980,2012)+scale_x_continuous(breaks= c(1981,1986,1991,1996,2001,2006,2011))+
  theme(plot.title = element_text(size = 15, face = "bold"))+
  theme(text = element_text(size = 15, face="bold"))+
  theme(axis.text.x = element_text(angle = 90))

meandeg<-ggplot(data=met.table, aes(met.table[,1],met.table[,5]))+
  geom_line(size=1.1)+geom_point(size=2)+labs(title = "Mean degree", tag = "d)")+xlab("Year")+ylab("Mean degree")+
  ylim(9,12)+xlim(1980,2012)+scale_x_continuous(breaks= c(1981,1986,1991,1996,2001,2006,2011))+
  theme(plot.title = element_text(size = 15, face = "bold"))+
  theme(text = element_text(size = 15, face="bold"))+
  theme(axis.text.x = element_text(angle = 90))

meangen<-ggplot(data=met.table, aes(met.table[,1],met.table[,8]))+
  geom_line(size=1.1)+geom_point(size=2)+labs(title = "Mean generality", tag = "e)")+xlab("Year")+ylab("Mean generality")+
  ylim(4.5,6.5)+xlim(1980,2012)+scale_x_continuous(breaks= c(1981,1986,1991,1996,2001,2006,2011))+
  theme(plot.title = element_text(size = 15, face = "bold"))+
  theme(text = element_text(size = 15, face="bold"))+
  theme(axis.text.x = element_text(angle = 90))

meanvuln<-ggplot(data=met.table, aes(met.table[,1],met.table[,9]))+
  geom_line(size=1.1)+geom_point(size=2)+labs(title = "Mean vulnerability", tag = "f)")+xlab("Year")+ylab("Mean vulnerability")+
  ylim(4.5,6.5)+xlim(1980,2012)+scale_x_continuous(breaks= c(1981,1986,1991,1996,2001,2006,2011))+
  theme(plot.title = element_text(size = 15, face = "bold"))+
  theme(text = element_text(size = 15, face="bold"))+
  theme(axis.text.x = element_text(angle = 90))


wqg<-ggplot(data=met.table, aes(met.table[,1],met.table[,10]))+
  geom_line(size=1.1)+geom_point(size=2)+labs(title = "Weighted generality", tag = "h)")+xlab("Year")+ylab("Weighted generality")+
  ylim(1,1.5)+xlim(1980,2012)+scale_x_continuous(breaks= c(1981,1986,1991,1996,2001,2006,2011))+
  theme(plot.title = element_text(size = 15, face = "bold"))+
  theme(text = element_text(size = 15, face="bold"))+
  theme(axis.text.x = element_text(angle = 90))

wqv<-ggplot(data=met.table, aes(met.table[,1],met.table[,11]))+
  geom_line(size=1.1)+geom_point(size=2)+labs(title = "Weighted vulnerability", tag = "i)")+xlab("Year")+ylab("Weighted vulnerability")+
  ylim(1.1,5.5)+xlim(1980,2012)+scale_x_continuous(breaks= c(1981,1986,1991,1996,2001,2006,2011))+
  theme(plot.title = element_text(size = 15, face = "bold"))+
  theme(text = element_text(size = 15, face="bold"))+
  theme(axis.text.x = element_text(angle = 90))

modplot<-ggplot(data=met.table, aes(met.table[,1],met.table[,12]))+
  geom_line(size=1.1)+geom_point(size=2)+labs(title = "Modularity", tag = "g)")+xlab("Year")+ylab("Modularity")+
  ylim(0.17,0.24)+xlim(1980,2012)+scale_x_continuous(breaks= c(1981,1986,1991,1996,2001,2006,2011))+
  theme(plot.title = element_text(size = 15, face = "bold"))+
  theme(text = element_text(size = 15, face="bold"))+
  theme(axis.text.x = element_text(angle = 90))

wmodplot<-ggplot(data=met.table, aes(met.table[,1],met.table[,13]))+
  geom_line(size=1.1)+geom_point(size=2)+labs(title = "Weighted modularity", tag = "j)")+xlab("Year")+ylab("Modularity")+
  ylim(0,0.48)+xlim(1980,2012)+scale_x_continuous(breaks= c(1981,1986,1991,1996,2001,2006,2011))+
  theme(plot.title = element_text(size = 15, face = "bold"))+
  theme(text = element_text(size = 15, face="bold"))+
  theme(axis.text.x = element_text(angle = 90))
  
sp+link+con+meandeg+meangen+meanvuln+modplot+wqg+wqv+wmodplot


modplot<-ggplot(data=met.table, aes(met.table[,1],met.table[,12]))+
  geom_line(size=1.5)+geom_point(size=3)+labs(title = "Modularity")+xlab("Year")+ylab("Modularity")+
  ylim(0.0,0.5)+xlim(1980,2012)+scale_x_continuous(breaks= c(1981,1986,1991,1996,2001,2006,2011))+
  theme(plot.title = element_text(size = 15, face = "bold"))+
  theme(text = element_text(size = 15, face="bold"))+
  theme(axis.text.x = element_text(angle = 90))

wmodplot<-ggplot(data=met.table, aes(met.table[,1],met.table[,13]))+
  geom_line(size=1.5)+geom_point(size=3)+labs(title = "Weighted modularity")+xlab("Year")+ylab("Modularity")+
  ylim(0,0.5)+xlim(1980,2012)+scale_x_continuous(breaks= c(1981,1986,1991,1996,2001,2006,2011))+
  theme(plot.title = element_text(size = 15, face = "bold"))+
  theme(text = element_text(size = 15, face="bold"))+
  theme(axis.text.x = element_text(angle = 90))
modplot+wmodplot

#######Null model test#######

#Metaweb nodes=34, links=207, modularity= 0.1971108
#modularity_MW<-modularity(Metaweb.el,mod.aff_MW)

#rand_webs_mod_list=vector(length=1000)# empty vector
#set.seed(1)
#for (i in 1:1000){
#rand_graph=erdos.renyi.game(34,207,type="gnm", directed = T, loops = F) # random graph
#llrand<-get.data.frame(rand_graph)
#llrand<-cbind(llrand, Weight = 1)#lägg till Weight column med värde 1, behövs för create_monolayer_object funktionen
                                  #funkar inte av nån anledning nu, creat_monolayer_object funktionen ger error
#Can't join on `x$from` x `y$from` because of incompatible types. #i `x$from` is of type <double>>. #i `y$from` is of type <character>>.
#llrand$from<-as.character(llrand$from)#node number till char för att motsvara "name"
#llrand$to<-as.character(llrand$to)#node number till char för att motsvara "namne"
#rand_mon_obj<-create_monolayer_object(llrand, directed = T, bipartite = F) #fungerar nu
#rand_infMap<-run_infomap_monolayer(rand_mon_obj, two_level = T, flow_model = "rawdir", seed = 123, trials = 100, silent = T,
#                                   ... = '--markov-time 1.05')
#rand.mod<-rand_infMap$modules$module_level1
#rand_webs_mod_list[i]=modularity(rand_graph,rand.mod) #store the modularity of random graph as the element of the vector
#}
#rand_webs_mod_list

#plot(rand_webs_mod_list)

## Null model curveball##

source("curveball.R") 
source("curv.R")

#load the time by species matrix 
poly<-read.table("clipboard", header = T) #you can change the name of course but here this is called poly
poly[7,22]=0
#poly2<-read.table("bok1.csv", header = T, sep = ";") # Neogobus completly removed. Use this one??
#poly<-read.csv("bok1.csv", header = T, sep = ";")
#poly<-read.table("bok2.csv", header = T, check.names = F, sep = ";") # neogobius changed to 0 for 2011, then set as not to be randomized
#rownames(poly)<-c(1981, 1986, 1991, 1996, 2001, 2006, 2011)
#load the food web metaweb 
#fw<-read.table("fw.metaweb.baltic.txt")
fw<-get.adjacency(net)
#poly<-poly[,-22]
# need to specify names of basal species, this example (9 species) is from the Barents Sea data set
basal=colnames(poly[TRUE,c("Autotroph","Mixotroph","Detritus","Neogobius melanostomus")]) # names of basal species, not to be randomized 
                                          # Neogobius as 0 and not to be randomized
model_curveball<-list()
iter<-999  # no of iterations
for (l in 1:iter){ # loop on randomisations
  
  null.model.curveball<-matrix(0, dim(poly)[1], dim(poly)[2]) 			# empty matrix with rows (number of rows are equal to the number of food webs) to be filled in with swapped elements (i.e. prey items)
  
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

#saveRDS(model_curveball, "model_curveball.rds") # old curveball
#model_curveball_long<-model_curveball

fw.1981<-list()
for (i in 1:iter){
  id.s<-which(model_curveball[[i]][1,]==1)  # species present in YEAR i
  fw.1981[[i]]<-as.matrix(fw[id.s,id.s])
}



rand.mod.1981<-c()
for (i in 1:iter){
  adj<-graph.adjacency(fw.1981[[i]])
  el<-get.data.frame(adj)
  el<-cbind(el, Weight = 1)
  ran.mon.obj<-create_monolayer_object(el, directed = T, bipartite = F)
  rand.infMap<-run_infomap_monolayer(ran.mon.obj, two_level = T, flow_model = "rawdir", seed = 123, trials = 100, silent = T,
                                    ... = '--markov-time 1.05')
  rand.mod<-rand.infMap$modules$module_level1
  rand.mod.1981[i]<-modularity(adj,rand.mod)
}


fw.1986<-list()
for (i in 1:iter){
  id.s<-which(model_curveball[[i]][2,]==1)  # species present in YEAR i
  fw.1986[[i]]<-as.matrix(fw[id.s,id.s])
}

rand.mod.1986<-c()
for (i in 1:iter){
  adj<-graph.adjacency(fw.1986[[i]])
  el<-get.data.frame(adj)
  el<-cbind(el, Weight = 1)
  ran.mon.obj<-create_monolayer_object(el, directed = T, bipartite = F)
  rand.infMap<-run_infomap_monolayer(ran.mon.obj, two_level = T, flow_model = "rawdir", seed = 123, trials = 100, silent = T,
                                     ... = '--markov-time 1.05')
  rand.mod<-rand.infMap$modules$module_level1
  rand.mod.1986[i]<-modularity(adj,rand.mod)
}


fw.1991<-list()
for (i in 1:iter){
  id.s<-which(model_curveball[[i]][3,]==1)  # species present in YEAR i
  fw.1991[[i]]<-as.matrix(fw[id.s,id.s])
}

rand.mod.1991<-c()
for (i in 1:iter){
  adj<-graph.adjacency(fw.1991[[i]])
  el<-get.data.frame(adj)
  el<-cbind(el, Weight = 1)
  ran.mon.obj<-create_monolayer_object(el, directed = T, bipartite = F)
  rand.infMap<-run_infomap_monolayer(ran.mon.obj, two_level = T, flow_model = "rawdir", seed = 123, trials = 100, silent = T,
                                     ... = '--markov-time 1.05')
  rand.mod<-rand.infMap$modules$module_level1
  rand.mod.1991[i]<-modularity(adj,rand.mod)
}


fw.1996<-list()
for (i in 1:iter){
  id.s<-which(model_curveball[[i]][4,]==1)  # species present in YEAR i
  fw.1996[[i]]<-as.matrix(fw[id.s,id.s])
}

rand.mod.1996<-c()
for (i in 1:iter){
  adj<-graph.adjacency(fw.1996[[i]])
  el<-get.data.frame(adj)
  el<-cbind(el, Weight = 1)
  ran.mon.obj<-create_monolayer_object(el, directed = T, bipartite = F)
  rand.infMap<-run_infomap_monolayer(ran.mon.obj, two_level = T, flow_model = "rawdir", seed = 123, trials = 100, silent = T,
                                     ... = '--markov-time 1.05')
  rand.mod<-rand.infMap$modules$module_level1
  rand.mod.1996[i]<-modularity(adj,rand.mod)
}

fw.2001<-list()
for (i in 1:iter){
  id.s<-which(model_curveball[[i]][5,]==1)  # species present in YEAR i
  fw.2001[[i]]<-as.matrix(fw[id.s,id.s])
}

rand.mod.2001<-c()
for (i in 1:iter){
  adj<-graph.adjacency(fw.2001[[i]])
  el<-get.data.frame(adj)
  el<-cbind(el, Weight = 1)
  ran.mon.obj<-create_monolayer_object(el, directed = T, bipartite = F)
  rand.infMap<-run_infomap_monolayer(ran.mon.obj, two_level = T, flow_model = "rawdir", seed = 123, trials = 100, silent = T,
                                     ... = '--markov-time 1.05')
  rand.mod<-rand.infMap$modules$module_level1
  rand.mod.2001[i]<-modularity(adj,rand.mod)
}

fw.2006<-list()
for (i in 1:iter){
  id.s<-which(model_curveball[[i]][6,]==1)  # species present in YEAR i
  fw.2006[[i]]<-as.matrix(fw[id.s,id.s])
}

rand.mod.2006<-c()
for (i in 1:iter){
  adj<-graph.adjacency(fw.2006[[i]])
  el<-get.data.frame(adj)
  el<-cbind(el, Weight = 1)
  ran.mon.obj<-create_monolayer_object(el, directed = T, bipartite = F)
  rand.infMap<-run_infomap_monolayer(ran.mon.obj, two_level = T, flow_model = "rawdir", seed = 123, trials = 100, silent = T,
                                     ... = '--markov-time 1.05')
  rand.mod<-rand.infMap$modules$module_level1
  rand.mod.2006[i]<-modularity(adj,rand.mod)
}

fw.2011<-list()
for (i in 1:iter){
  id.s<-which(model_curveball[[i]][7,]==1)  # species present in YEAR i
  fw.2011[[i]]<-as.matrix(fw[id.s,id.s])
}

rand.mod.2011<-c()
for (i in 1:iter){
  adj<-graph.adjacency(fw.2011[[i]])
  el<-get.data.frame(adj)
  el<-cbind(el, Weight = 1)
  ran.mon.obj<-create_monolayer_object(el, directed = T, bipartite = F)
  rand.infMap<-run_infomap_monolayer(ran.mon.obj, two_level = T, flow_model = "rawdir", seed = 123, trials = 100, silent = T,
                                     ... = '--markov-time 1.05')
  rand.mod<-rand.infMap$modules$module_level1
  rand.mod.2011[i]<-modularity(adj,rand.mod)
}


hist(rand.mod.1981)
hist(rand.mod.1986)
hist(rand.mod.1991)
hist(rand.mod.1996)
hist(rand.mod.2001)
hist(rand.mod.2011)


obs.mod.1981<-modularity(fw.1981.el, IM_uw_1981$modules$module_level1)
obs.mod.1986<-modularity(fw.1986.el, IM_uw_1986$modules$module_level1)
obs.mod.1991<-modularity(fw.1991.el, IM_uw_1991$modules$module_level1)
obs.mod.1996<-modularity(fw.1996.el, IM_uw_1996$modules$module_level1)
obs.mod.2001<-modularity(fw.2001.el, IM_uw_2001$modules$module_level1)
obs.mod.2006<-modularity(fw.2006.el, IM_uw_2006$modules$module_level1)
obs.mod.2011<-modularity(fw.2011.el, IM_uw_2011$modules$module_level1)




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
ecdf(rand.mod.1986)(obs.mod.1986)
ecdf(rand.mod.1991)(obs.mod.1991)
ecdf(rand.mod.1996)(obs.mod.1996)
ecdf(rand.mod.2001)(obs.mod.2001)
ecdf(rand.mod.2006)(obs.mod.2006)
ecdf(rand.mod.2011)(obs.mod.2011)

# not sure if it works with multimodal dist???
rm2011.cdf<-ecdf(rand.mod.2011) #mod.vec is your vector with randomized modularity values
emp2011.cdf<-rm2011.cdf(obs.mod.2011) # empirical/observed modularity value
emp2011.cdf
plot(rm2011.cdf)
abline(v=obs.mod.2011)


plot(ecdf(rand.mod.1981))
points(x=obs.mod.1981, y=(ecdf(rand.mod.1981)(obs.mod.1981)), col="red", pch=16)

plot(ecdf(rand.mod.1986))
points(x=obs.mod.1986, y=ecdf(rand.mod.1986)(obs.mod.1986), col="red", pch=16)

plot(ecdf(rand.mod.1991))
points(x=obs.mod.1991, y=ecdf(rand.mod.1991)(obs.mod.1991), col="red", pch=16)

plot(ecdf(rand.mod.1996))
points(x=obs.mod.1996, y=ecdf(rand.mod.1996)(obs.mod.1996), col="red", pch=16)

plot(ecdf(rand.mod.2001))
points(x=obs.mod.2001, y=ecdf(rand.mod.2001)(obs.mod.2001), col="red", pch=16)

plot(ecdf(rand.mod.2006))
points(x=obs.mod.2006, y=ecdf(rand.mod.2006)(obs.mod.2006), col="red", pch=16)

plot(ecdf(rand.mod.2011))
points(x=obs.mod.2011, y=ecdf(rand.mod.2011)(obs.mod.2011), col="red", pch=16)





### weightwed null model

#biomass data yearly
biom.data<-read.csv("biom.dat.csv")
rownames(biom.data)<-c(1981, 1986, 1991, 1996, 2001, 2006, 2011)

#binary randomisations time by species matrix
#mods_curve<-readRDS("model_curveball.rds")
mods_curve<-model_curveball
#remove neogobius
biom.data<-biom.data[,-22]

weighted_ran_mat<-list()
for (i in 1:length(mods_curve)){
  ran.mat<-mods_curve[[i]][,-22] #removes neogobius
  for (j in 1:dim(ran.mat)[2]){
    
    id<-which(ran.mat[,j]==1) # row id of species j 
    
    if (length(id)>0){
      non.zero.vec<-biom.data[,j][biom.data[,j]!=0] #remove the zeros
      
      ran.mat[id,j]<-sample(non.zero.vec) #randomise biomass info and put back into the row id's
    }
  }
  weighted_ran_mat[[i]]<-ran.mat
}



#str(info)

ran.info<-info[-22,]
rownames(ran.info)

ran.info.1981<-list()
for (i in 1:length(weighted_ran_mat)){
 # i<-1
  ran.info[,4]<-weighted_ran_mat[[i]][1,]
  ran.info.1981[[i]]<-ran.info
}

#ids<-which(ran.info.1981[[1]]$meanB=="0")
#fluxes <- fluxing(fw.1981[[1]], ran.info.1981[[1]][-ids,]$meanB, ran.info.1981[[1]][-ids,]$losses, 
#                  ran.info.1981[[1]][-ids,]$efficiencies, ef.level="prey")
#rownames(fw.1981[[1]])
#ids<-which(ran.info.1981[[1]]$meanB=="0")
#rownames(ran.info.1981[[1]][-ids,])

#fw.1981<-list()
#for (i in 1:iter){
#  id.s<-which(model_curveball[[i]][1,]==1)  # species present in YEAR i
#  fw.1981[[i]]<-as.matrix(fw[id.s,id.s])
#}


#make loop
weight.ran.flux1981<-list()
for (i in 1:iter){
  ids<-which(ran.info.1981[[i]]$meanB=="0")
  fluxes <- fluxing(fw.1981[[i]], ran.info.1981[[i]][-ids,]$meanB, ran.info.1981[[1]][-ids,]$losses, 
                    ran.info.1981[[i]][-ids,]$efficiencies, ef.level="prey")
  fluxes <- fluxes*86.4 # conversion from J/sec to kJ/day 
  # 1 J/sec = 86.4 kJ/day (there are 86400 sec/day)
  weight.ran.flux1981[[i]]<-fluxes
}

w.rand.mod.1981<-c()
for (i in 1:iter){
  adj<-graph_from_adjacency_matrix(weight.ran.flux1981[[i]], weighted=TRUE)
  el<-get.data.frame(adj)
  ran.mon.obj<-create_monolayer_object(el, directed = T, bipartite = F)
  rand.infMap<-run_infomap_monolayer(ran.mon.obj, two_level = T, flow_model = "rawdir", seed = 123, trials = 100, silent = T,
                                     ... = '--markov-time 1.05')
  rand.mod<-rand.infMap$modules$module_level1
  w.rand.mod.1981[i]<-modularity(adj,rand.mod, weights =el$weight)
}

w.obs.mod1981<-modularity(wfw.1981.el, IM_w_1981$modules$module_level1, weights = get.data.frame(wfw.1981.el)$weight)
ecdf(w.rand.mod.1981)(w.obs.mod1981)
plot(ecdf(w.rand.mod.1981))


ran.info.1986<-list()
for (i in 1:length(weighted_ran_mat)){
  ran.info[,4]<-weighted_ran_mat[[i]][2,]
  ran.info.1986[[i]]<-ran.info
}

weight.ran.flux1986<-list()
for (i in 1:iter){
  ids<-which(ran.info.1986[[i]]$meanB=="0")
  fluxes <- fluxing(fw.1986[[i]], ran.info.1986[[i]][-ids,]$meanB, ran.info.1986[[1]][-ids,]$losses, 
                    ran.info.1986[[i]][-ids,]$efficiencies, ef.level="prey")
  fluxes <- fluxes*86.4 # conversion from J/sec to kJ/day 
  # 1 J/sec = 86.4 kJ/day (there are 86400 sec/day)
  weight.ran.flux1986[[i]]<-fluxes
}

w.rand.mod.1986<-c()
for (i in 1:iter){
  adj<-graph_from_adjacency_matrix(weight.ran.flux1986[[i]], weighted=TRUE)
  el<-get.data.frame(adj)
  ran.mon.obj<-create_monolayer_object(el, directed = T, bipartite = F)
  rand.infMap<-run_infomap_monolayer(ran.mon.obj, two_level = T, flow_model = "rawdir", seed = 123, trials = 100, silent = T,
                                     ... = '--markov-time 1.05')
  rand.mod<-rand.infMap$modules$module_level1
  w.rand.mod.1986[i]<-modularity(adj,rand.mod, weights =el$weight)
}

w.obs.mod1986<-modularity(wfw.1986.el, IM_w_1986$modules$module_level1, weights = get.data.frame(wfw.1986.el)$weight)
ecdf(w.rand.mod.1986)(w.obs.mod1986)
plot(ecdf(w.rand.mod.1986))


ran.info.1991<-list()
for (i in 1:length(weighted_ran_mat)){
  ran.info[,4]<-weighted_ran_mat[[i]][3,]
  ran.info.1991[[i]]<-ran.info
}

weight.ran.flux1991<-list()
for (i in 1:iter){
  ids<-which(ran.info.1991[[i]]$meanB=="0")
  fluxes <- fluxing(fw.1991[[i]], ran.info.1991[[i]][-ids,]$meanB, ran.info.1991[[1]][-ids,]$losses, 
                    ran.info.1991[[i]][-ids,]$efficiencies, ef.level="prey")
  fluxes <- fluxes*86.4 # conversion from J/sec to kJ/day 
  # 1 J/sec = 86.4 kJ/day (there are 86400 sec/day)
  weight.ran.flux1991[[i]]<-fluxes
}

w.rand.mod.1991<-c()
for (i in 1:iter){
  adj<-graph_from_adjacency_matrix(weight.ran.flux1991[[i]], weighted=TRUE)
  el<-get.data.frame(adj)
  ran.mon.obj<-create_monolayer_object(el, directed = T, bipartite = F)
  rand.infMap<-run_infomap_monolayer(ran.mon.obj, two_level = T, flow_model = "rawdir", seed = 123, trials = 100, silent = T,
                                     ... = '--markov-time 1.05')
  rand.mod<-rand.infMap$modules$module_level1
  w.rand.mod.1991[i]<-modularity(adj,rand.mod, weights =el$weight)
}

w.obs.mod1991<-modularity(wfw.1991.el, IM_w_1991$modules$module_level1, weights = get.data.frame(wfw.1991.el)$weight)
ecdf(w.rand.mod.1991)(w.obs.mod1991)
plot(ecdf(w.rand.mod.1991))


ran.info.1996<-list()
for (i in 1:length(weighted_ran_mat)){
  ran.info[,4]<-weighted_ran_mat[[i]][4,]
  ran.info.1996[[i]]<-ran.info
}

weight.ran.flux1996<-list()
for (i in 1:iter){
  ids<-which(ran.info.1996[[i]]$meanB=="0")
  fluxes <- fluxing(fw.1996[[i]], ran.info.1996[[i]][-ids,]$meanB, ran.info.1996[[1]][-ids,]$losses, 
                    ran.info.1996[[i]][-ids,]$efficiencies, ef.level="prey")
  fluxes <- fluxes*86.4 # conversion from J/sec to kJ/day 
  # 1 J/sec = 86.4 kJ/day (there are 86400 sec/day)
  weight.ran.flux1996[[i]]<-fluxes
}

w.rand.mod.1996<-c()
for (i in 1:iter){
  adj<-graph_from_adjacency_matrix(weight.ran.flux1996[[i]], weighted=TRUE)
  el<-get.data.frame(adj)
  ran.mon.obj<-create_monolayer_object(el, directed = T, bipartite = F)
  rand.infMap<-run_infomap_monolayer(ran.mon.obj, two_level = T, flow_model = "rawdir", seed = 123, trials = 100, silent = T,
                                     ... = '--markov-time 1.05')
  rand.mod<-rand.infMap$modules$module_level1
  w.rand.mod.1996[i]<-modularity(adj,rand.mod, weights =el$weight)
}

w.obs.mod1996<-modularity(wfw.1996.el, IM_w_1996$modules$module_level1, weights = get.data.frame(wfw.1996.el)$weight)
ecdf(w.rand.mod.1996)(w.obs.mod1996)
plot(ecdf(w.rand.mod.1996))


ran.info.2001<-list()
for (i in 1:length(weighted_ran_mat)){
  ran.info[,4]<-weighted_ran_mat[[i]][5,]
  ran.info.2001[[i]]<-ran.info
}

weight.ran.flux2001<-list()
for (i in 1:iter){
  ids<-which(ran.info.2001[[i]]$meanB=="0")
  fluxes <- fluxing(fw.2001[[i]], ran.info.2001[[i]][-ids,]$meanB, ran.info.2001[[1]][-ids,]$losses, 
                    ran.info.2001[[i]][-ids,]$efficiencies, ef.level="prey")
  fluxes <- fluxes*86.4 # conversion from J/sec to kJ/day 
  # 1 J/sec = 86.4 kJ/day (there are 86400 sec/day)
  weight.ran.flux2001[[i]]<-fluxes
}

w.rand.mod.2001<-c()
for (i in 1:iter){
  adj<-graph_from_adjacency_matrix(weight.ran.flux2001[[i]], weighted=TRUE)
  el<-get.data.frame(adj)
  ran.mon.obj<-create_monolayer_object(el, directed = T, bipartite = F)
  rand.infMap<-run_infomap_monolayer(ran.mon.obj, two_level = T, flow_model = "rawdir", seed = 123, trials = 100, silent = T,
                                     ... = '--markov-time 1.05')
  rand.mod<-rand.infMap$modules$module_level1
  w.rand.mod.2001[i]<-modularity(adj,rand.mod, weights =el$weight)
}

w.obs.mod2001<-modularity(wfw.2001.el, IM_w_2001$modules$module_level1, weights = get.data.frame(wfw.2001.el)$weight)
ecdf(w.rand.mod.2001)(w.obs.mod2001)
plot(ecdf(w.rand.mod.2001))


ran.info.2006<-list()
for (i in 1:length(weighted_ran_mat)){
  ran.info[,4]<-weighted_ran_mat[[i]][6,]
  ran.info.2006[[i]]<-ran.info
}

weight.ran.flux2006<-list()
for (i in 1:iter){
  ids<-which(ran.info.2006[[i]]$meanB=="0")
  fluxes <- fluxing(fw.2006[[i]], ran.info.2006[[i]][-ids,]$meanB, ran.info.2006[[1]][-ids,]$losses, 
                    ran.info.2006[[i]][-ids,]$efficiencies, ef.level="prey")
  fluxes <- fluxes*86.4 # conversion from J/sec to kJ/day 
  # 1 J/sec = 86.4 kJ/day (there are 86400 sec/day)
  weight.ran.flux2006[[i]]<-fluxes
}

w.rand.mod.2006<-c()
for (i in 1:iter){
  adj<-graph_from_adjacency_matrix(weight.ran.flux2006[[i]], weighted=TRUE)
  el<-get.data.frame(adj)
  ran.mon.obj<-create_monolayer_object(el, directed = T, bipartite = F)
  rand.infMap<-run_infomap_monolayer(ran.mon.obj, two_level = T, flow_model = "rawdir", seed = 123, trials = 100, silent = T,
                                     ... = '--markov-time 1.05')
  rand.mod<-rand.infMap$modules$module_level1
  w.rand.mod.2006[i]<-modularity(adj,rand.mod, weights =el$weight)
}

w.obs.mod2006<-modularity(wfw.2006.el, IM_w_2006$modules$module_level1, weights = get.data.frame(wfw.2006.el)$weight)
ecdf(w.rand.mod.2006)(w.obs.mod2006)
plot(ecdf(w.rand.mod.2006))


ran.info.2011<-list()
for (i in 1:length(weighted_ran_mat)){
  ran.info[,4]<-weighted_ran_mat[[i]][7,]
  ran.info.2011[[i]]<-ran.info
}

weight.ran.flux2011<-list()
for (i in 1:iter){
  ids<-which(ran.info.2011[[i]]$meanB=="0")
  fluxes <- fluxing(fw.2011[[i]], ran.info.2011[[i]][-ids,]$meanB, ran.info.2011[[1]][-ids,]$losses, 
                    ran.info.2011[[i]][-ids,]$efficiencies, ef.level="prey")
  fluxes <- fluxes*86.4 # conversion from J/sec to kJ/day 
  # 1 J/sec = 86.4 kJ/day (there are 86400 sec/day)
  weight.ran.flux2011[[i]]<-fluxes
}

w.rand.mod.2011<-c()
for (i in 1:iter){
  adj<-graph_from_adjacency_matrix(weight.ran.flux2011[[i]], weighted=TRUE)
  el<-get.data.frame(adj)
  ran.mon.obj<-create_monolayer_object(el, directed = T, bipartite = F)
  rand.infMap<-run_infomap_monolayer(ran.mon.obj, two_level = T, flow_model = "rawdir", seed = 123, trials = 100, silent = T,
                                     ... = '--markov-time 1.05')
  rand.mod<-rand.infMap$modules$module_level1
  w.rand.mod.2011[i]<-modularity(adj,rand.mod, weights =el$weight)
}

w.obs.mod2011<-modularity(wfw.2011.el, IM_w_2011$modules$module_level1, weights = get.data.frame(wfw.2011.el)$weight)
ecdf(w.rand.mod.2011)(w.obs.mod2011)
plot(ecdf(w.rand.mod.2011))

w.obs.mod1981
w.obs.mod1986
w.obs.mod1991
w.obs.mod1996
w.obs.mod2001
w.obs.mod2006
w.obs.mod2011

hist(w.rand.mod.1981)
hist(w.rand.mod.1986)
hist(w.rand.mod.1991)
hist(w.rand.mod.1996)
hist(w.rand.mod.2001)
hist(w.rand.mod.2006)
hist(w.rand.mod.2011)

ecdf(w.rand.mod.1981)(w.obs.mod1981)
ecdf(w.rand.mod.1986)(w.obs.mod1986)
ecdf(w.rand.mod.1991)(w.obs.mod1991)
ecdf(w.rand.mod.1996)(w.obs.mod1996)
ecdf(w.rand.mod.2001)(w.obs.mod2001)
ecdf(w.rand.mod.2006)(w.obs.mod2006)
ecdf(w.rand.mod.2011)(w.obs.mod2011)


ecdf(rand.mod.1981)(obs.mod.1981)
ecdf(rand.mod.1986)(obs.mod.1986)
ecdf(rand.mod.1991)(obs.mod.1991)
ecdf(rand.mod.1996)(obs.mod.1996)
ecdf(rand.mod.2001)(obs.mod.2001)
ecdf(rand.mod.2006)(obs.mod.2006)
ecdf(rand.mod.2011)(obs.mod.2011)


par( mfrow= c(2,7) )
plot(ecdf(rand.mod.1981), main="1981", xlab="Modularity", xlim=c(-0.1,0.5))
points(x=obs.mod.1981, y=(ecdf(rand.mod.1981)(obs.mod.1981)), col="red", pch=16, cex= 1.5)
plot(ecdf(rand.mod.1986), main="1986", xlab="Modularity",xlim=c(-0.1,0.5))
points(x=obs.mod.1986, y=ecdf(rand.mod.1986)(obs.mod.1986), col="red", pch=16, cex= 1.5)
plot(ecdf(rand.mod.1991), main="1991", xlab="Modularity", xlim=c(-0.1,0.5))
points(x=obs.mod.1991, y=ecdf(rand.mod.1991)(obs.mod.1991), col="red", pch=16, cex= 1.5)
plot(ecdf(rand.mod.1996), main="1996", xlab="Modularity", xlim=c(-0.1,0.5))
points(x=obs.mod.1996, y=ecdf(rand.mod.1996)(obs.mod.1996), col="red", pch=16, cex= 1.5)
plot(ecdf(rand.mod.2001), main="2001", xlab="Modularity", xlim=c(-0.1,0.5))
points(x=obs.mod.2001, y=ecdf(rand.mod.2001)(obs.mod.2001), col="red", pch=16, cex= 1.5)
plot(ecdf(rand.mod.2006), main="2006", xlab="Modularity",xlim=c(-0.1,0.5))
points(x=obs.mod.2006, y=ecdf(rand.mod.2006)(obs.mod.2006), col="red", pch=16, cex= 1.5)
plot(ecdf(rand.mod.2011), main="2011", xlab="Modularity",xlim=c(-0.1,0.5))
points(x=obs.mod.2011, y=ecdf(rand.mod.2011)(obs.mod.2011), col="red", pch=16, cex= 1.5)

plot(ecdf(w.rand.mod.1981), do.points=T, main="Weighted 1981", xlab="Modularity",xlim=c(-0.1,0.5))
points(x=w.obs.mod1981, y=(ecdf(w.rand.mod.1981)(w.obs.mod1981)), col="red", pch=16, cex= 1.5)
plot(ecdf(w.rand.mod.1986), do.points=T, main="Weighted 1986", xlab="Modularity",xlim=c(-0.1,0.5))
points(x=w.obs.mod1986, y=ecdf(w.rand.mod.1986)(w.obs.mod1986), col="red", pch=16, cex= 1.5)
plot(ecdf(w.rand.mod.1991), do.points=T, main="Weighted 1991", xlab="Modularity",xlim=c(-0.1,0.5))
points(x=w.obs.mod1991, y=ecdf(w.rand.mod.1991)(w.obs.mod1991), col="red", pch=16, cex= 1.5)
plot(ecdf(w.rand.mod.1996), do.points=T, main="Weighted 1996", xlab="Modularity",xlim=c(-0.1,0.5))
points(x=w.obs.mod1996, y=ecdf(w.rand.mod.1996)(w.obs.mod1996), col="red", pch=16, cex= 1.5)
plot(ecdf(w.rand.mod.2001), do.points=T, main="Weighted 2001", xlab="Modularity",xlim=c(-0.1,0.5))
points(x=w.obs.mod2001, y=ecdf(w.rand.mod.2001)(w.obs.mod2001), col="red", pch=16, cex= 1.5)
plot(ecdf(w.rand.mod.2006), do.points=T, main="Weighted 2006", xlab="Modularity",xlim=c(-0.1,0.5))
points(x=w.obs.mod2006, y=ecdf(w.rand.mod.2006)(w.obs.mod2006), col="red", pch=16, cex= 1.5)
plot(ecdf(w.rand.mod.2011), do.points=T, main="Weighted 2011", xlab="Modularity",xlim=c(-0.1,0.5))
points(x=w.obs.mod2011, y=ecdf(w.rand.mod.2011)(w.obs.mod2011), col="red", pch=16, cex= 1.5)
