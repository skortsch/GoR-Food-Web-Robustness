#Metaweb of Bothnian sea

setwd("D:/R_Git/GoR-Food-webs/Bothnian_Sea")
library(igraph)		# Network analysis
library(NetIndices)	# Trophic level and omnivory 

source("TL function.R") # Load Owen Petchey code on TL and omnivory
source("SWTL Omnivory.R")
source("Paths function.R")

#load
pw<-read.table("pw.txt") #pairwise interaction list
TS.in.regions<-read.csv("Species_in_regions.csv", sep=";",check.names=FALSE) #species

poly.pa<-t(TS.in.regions[,-c(1,2)])
colnames(poly.pa)<-TS.in.regions[,1]

fw<-read.table("fw.metaweb.baltic.txt", h=T,  skipNul = TRUE)

#Bothnian Sea fw matrix
BS.ids<-which(poly.pa[3,]=="1")
BS.fw<-fw[BS.ids, BS.ids]
write.table(BS.fw,"BS_fw_mat.txt")
#Bothnian Bay fw matrix
BB.ids<-which(poly.pa[4,]=="1")
BB.fw<-fw[BB.ids, BB.ids]
write.table(BB.fw,"BS_fw_mat.txt")


#Metaweb of Bothnian sea
#Meta.BS.ids<-which(poly.pa[3:4,]=="1")
#Meta.BS.fw<-fw[Meta.BS.ids, Meta.BS.ids]


BS.edge.list<-graph.adjacency(as.matrix(BS.fw))
BB.edge.list<-graph.adjacency(as.matrix(BB.fw))

rbind(as_data_frame(BS.edge.list), as_data_frame(BB.edge.list))

V(BS.edge.list)$name
V(BB.edge.list)$name

attrs <- rbind(as_data_frame(BS.edge.list, "vertices"), as_data_frame(BB.edge.list, "vertices")) %>% unique()
el <- rbind(as_data_frame(BS.edge.list), as_data_frame(BB.edge.list))
new_g <- graph_from_data_frame(el, directed = T, vertices = attrs)

#Meta.GoB.fw<-rbind(as_data_frame(BS.edge.list), as_data_frame(BB.edge.list))
GoB_meta<-new_g

g<-graph.union(BS.edge.list, BB.edge.list)
#test<-graph.edgelist(as.matrix(Meta.GoB.fw))
#test2<-as.matrix(test, matrix.type = "adjacency", sparse=F)

#test2[test2 > 1] <- 1
#Meta.GoB.list<-graph.adjacency(test2)
t<-as.matrix(g, matrix.type = "adjacency", sparse=F)

TL.MetaGob<-TLs(t)[,1]
TL.BB<-TLs(as.matrix(BB.fw))[,1] 

layout.matrix.GoB<-matrix(nrow=length(V(g)), ncol=2)
layout.matrix.GoB[,1]<-runif(length(V(g))) # randomly assign along x-axis
layout.matrix.GoB[,2]<-TL.MetaGob # y-axis value based on trophic level
backup<-layout.matrix.GoB


FG<-read.table("FG.txt") # saknar vissa arter, därför blir FG.GoB fel!!!!

FG.GoB<- FG[which(FG[,1]%in%rownames(t)),] # något går fel, FG.gob blir 35 istället för 36 som är dim av matrix
groups.GoB<- FG.GoB[,3][order(match(FG.GoB[,1],colnames(t)))] # 1=non-organic,2=phytopl, 3=zoopl, 4=benthos, 5=fish, 6=birds, 7=mammals

#manual col, måste fixa senare
a<-matrix(ncol = 2, nrow = 36)
a[,1]<-as.matrix(rownames(t))
a[,2]<-c(2,2,1,2,2,2,2,2,2,3,3,3,3,3,3,3,3,4,4,4,4,4,5,5,5,5,5,5,5,5,5,5,3,3,5,5)
b<-c(2,2,1,2,2,2,2,2,2,3,3,3,3,3,3,3,3,4,4,4,4,4,5,5,5,5,5,5,5,5,5,5,3,3,5,5)

colbar.FG<- c("darkgrey","lightgreen","cyan","orange","blue","magenta","lightpink")
plot(g, layout=layout.matrix.GoB, vertex.size = 1/5*degree(g)+4, vertex.frame.color = "gray50", edge.arrow.size=0.2,
     vertex.color=colbar.FG[b],vertex.label = colnames(as.matrix(g)), vertex.label.cex=0.55, vertex.size=0.1,
     edge.color="grey40", edge.width=0.15, vertex.label.dist=0.7, vertex.label.color= "black")

tt<-layout_on_grid(g)
tt[,2]<-layout.matrix.GoB[,2]
tt[14,1]<-0
tt[6,1]<-4.5
tt[35,1]<-5
plot(g,layout=tt, vertex.size = 1/5*degree(g)+4, vertex.frame.color = "gray50", edge.arrow.size=0.2,
     vertex.color=colbar.FG[groups.GoB],vertex.label = colnames(as.matrix(g)), vertex.label.cex=0.55, vertex.size=0.1,
     edge.color="grey40", edge.width=0.15, vertex.label.dist=-0.5, vertex.label.color= "black")

plot(g,layout=tt, vertex.size = 2/5*degree(g)+4, vertex.frame.color = "gray50", edge.arrow.size=0.2,
     vertex.color=colbar.FG[b],vertex.label = NA, vertex.label.cex=0.55, vertex.size=0.1,
     edge.color="grey40", edge.width=0.15, vertex.label.dist=-0.5, vertex.label.color= "black")


#plotfw(g,vertex.size = 1/5*degree(g)+4, vertex.frame.color = "gray50", edge.arrow.size=0.2,
#       vertex.color=colbar.FG[b],vertex.label = colnames(as.matrix(g)), labcex=0.4, vertex.size=0.1,
#       edge.color="grey40", edge.width=0.15, vertex.label.dist=-0.5, vertex.label.color= "black")



#alla ater med
#load input files
pw<-read.table("pw.txt") #pairwise interaction list
TS.in.regions<-read.csv("Species_in_regions.csv", sep=";",check.names=FALSE) #species

poly.pa<-t(TS.in.regions[,-c(1,2)])
colnames(poly.pa)<-TS.in.regions[,1]

#make meta food web

TS<-TS.in.regions[,1]

fw<- matrix(0, length(TS), length(TS))
row.names(fw)<-TS
colnames(fw)<-TS

# Convert from pairwise list to matrix
for (i in 1:dim(fw)[2]){
  ID<-0
  for (j in 1:dim(pw)[1]){
    if(pw$PREDATOR[j]==as.vector(TS[i])) {ID<-c(ID,which(row.names(fw)==pw$PREY[j]))}
  }
  fw[ID,i]<-1
}

write.table(fw, "fw.txt")


#make local food webs

#Southern Baltic
#sp.SB<-which(poly.pa[1,]==1) 				# sp.p a
#fw.SB<-fw[sp.SB,sp.SB]

#Baltic Proper
#sp.BP<-which(poly.pa[2,]==1) 				# sp.p a
#fw.BP<-fw[sp.BP,sp.BP]

#Bothnian Sea
sp.BS<-which(poly.pa[3,]==1) 				# sp.p a
fw.BS<-fw[sp.BS,sp.BS]

#Bothnian Bay
sp.BB<-which(poly.pa[4,]==1) 				# sp.p a
fw.BB<-fw[sp.BB,sp.BB]



BS.el<-graph.adjacency(as.matrix(fw.BS))
BB.el<-graph.adjacency(as.matrix(fw.BB))

rbind(as_data_frame(BS.el), as_data_frame(BB.el))

V(BS.el)$name
V(BB.el)$name

attrs <- rbind(as_data_frame(BB.el, "vertices"), as_data_frame(BB.el, "vertices")) %>% unique()
el <- rbind(as_data_frame(BS.el), as_data_frame(BB.el))
new_g <- graph_from_data_frame(el, directed = T)

#Meta.GoB.fw<-rbind(as_data_frame(BS.edge.list), as_data_frame(BB.edge.list))
GoB_meta<-new_g

#g<-graph.union(BS.edge.list, BB.edge.list)
#test<-graph.edgelist(as.matrix(Meta.GoB.fw))
#test2<-as.matrix(test, matrix.type = "adjacency", sparse=F)

#test2[test2 > 1] <- 1
#Meta.GoB.list<-graph.adjacency(test2)
t<-as.matrix(GoB_meta, matrix.type = "adjacency", sparse=F)

as.matrix(t)[,1]
TL.MetaGob<-TLs(as.matrix(t))[,1]

layout.matrix.GoB<-matrix(nrow=length(V(GoB_meta)), ncol=2)
layout.matrix.GoB[,1]<-runif(length(V(GoB_meta))) # randomly assign along x-axis
layout.matrix.GoB[,2]<-TL.MetaGob # y-axis value based on trophic level

FG<-read.table("FG.txt") # saknar vissa arter, därför blir FG.GoB fel!!!!


FG.GoB<-FG[which(FG[,1]%in%colnames(t)),] # något går fel, FG.gob blir 31 istället för 35 som är dim av matrix
#sp.names do not all match

FG2<-read.table("FG2.txt") # better but some still missing
FG.GoB<-FG2[which(FG2[,1]%in%colnames(t)),] # något går fel, FG.gob blir 31 istället för 35 som är dim av matrix
groups.GoB<- FG.GoB[,3][order(match(FG.GoB[,1],colnames(t)))] # 1=non-organic,2=phytopl, 3=zoopl, 4=benthos, 5=fish, 6=birds, 7=mammals

plot(GoB_meta, layout= layout_on_grid(GoB_meta) ,vertex.size = 7, vertex.frame.color = "gray50", edge.arrow.size=0.2,
     vertex.color=colbar.FG[groups.GoB],vertex.label = colnames(as.matrix(GoB_meta)), vertex.label.cex=0.55, vertex.size=0.1,
     edge.color="grey40", edge.width=0.15, vertex.label.dist=0.7, vertex.label.color= "black")



