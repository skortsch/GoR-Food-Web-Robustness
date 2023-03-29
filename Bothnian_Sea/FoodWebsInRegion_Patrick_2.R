#### Food webs along gradients in the Batic Sea ####

setwd("E:/BONUS BLUEWEBS project/R scripts food webs")

# Load packages
library(igraph)		# Network analysis
library(NetIndices)	# Trophic level and omnivory 

source("TL function.R") # Load Owen Petchey code on TL and omnivory
source("SWTL Omnivory.R")
source("Paths function.R")


#load input files
#load input files
pw<-read.table("pw.txt") #pairwise interaction list
TS.in.regions<-read.csv("Species_in_regions.csv", sep=";",check.names=FALSE) #species

poly.pa<-t(TS.in.regions[,-c(1,2)])
colnames(poly.pa)<-TS.in.regions[,1]

#make meta food web
#fw<- matrix(0, length(TS), length(TS))
#row.names(fw)<-TS
#colnames(fw)<-TS

# Convert from pairwise list to matrix
#for (i in 1:dim(fw)[2]){
#  ID<-0
#  for (j in 1:dim(pw)[1]){
#    if(pw$PREDATOR[j]==as.vector(TS[i])) {ID<-c(ID,which(row.names(fw)==pw$PREY[j]))}
 # }
#  fw[ID,i]<-1
#}

#write.table(fw, "fw.metaweb.baltic.txt")

fw<-read.table("fw.metaweb.baltic.txt", h=T,  skipNul = TRUE)

#Southern Baltic fw matrix
SB.ids<-which(poly.pa[1,]=="1")
SB.fw<-fw[SB.ids, SB.ids]
write.table(SB.fw,"SB_fw_mat.txt")
#Baltic Proper fw matrix
BP.ids<-which(poly.pa[2,]=="1")
BP.fw<-fw[BP.ids, BP.ids]
write.table(BP.fw,"BP_fw_mat.txt")
#Bothnian Sea fw matrix
BS.ids<-which(poly.pa[3,]=="1")
BS.fw<-fw[BS.ids, BS.ids]
write.table(BS.fw,"BS_fw_mat.txt")
#Bothnian Bay fw matrix
BB.ids<-which(poly.pa[4,]=="1")
BB.fw<-fw[BB.ids, BB.ids]
write.table(BB.fw,"BS_fw_mat.txt")


#make food webs

metric.names<-c("S","L", "LD", "C","Clust","Mod","PercOmni", "Omni", "PredPrey", "Can", "Bas", "Top", "Int", "GenSD", "VulSD", "meanSWTL", "meanLWTL", "meanPath","longPath") # "shortPath" 
fwind<-matrix(NA,dim(poly.pa)[1],length(metric.names))			# matrix to be used for the loop
colnames(fwind)<-c(metric.names)
rownames(fwind)<-rownames(poly.pa)


### The loop to calculate the metrics and fill them into the matrix "fwind.p" given above

for (i in 1:dim(poly.pa)[1]) {
  print(i)
  
  sp.p<-which(poly.pa[i,]==1) 				# sp.p a
  fw0.p<-fw[sp.p,sp.p]					#s remove singeltons from polygon
  
  id1<-which(apply(fw0.p,1,sum)==0)
  id2<-which(apply(fw0.p,2,sum)==0)
  idd1<- which(1*(id1%in%id2)==1)
  if (length(idd1)==0) fw1.p=fw0.p else fw1.p=fw0.p[-id1[idd1],-id1[idd1]]			
  
  fw1.pi<-graph.adjacency(fw1.p)
  
  # Fundamental complexity measures
  fwind[i,1]<- dim(fw1.p)[1]			 	# Number of nodes (trophospecies)
  fwind[i,2]<- sum(fw1.p)    				# Number of links (feeding relationships)
  fwind[i,3]<- sum(fw1.p)/dim(fw1.p)[1]		# Link density
  fwind[i,4]<-fwind[i,2]/fwind[i,1]^2      		# Connectance (realized links relative to maximum possible n links)

  # Community structureal measures
  Clust<- transitivity(fw1.pi,vids=NULL)
  fwind[i, 5]<- Clust
  spingc <- spinglass.community(fw1.pi) 		# Compartmentalisation estimation based on simulated annealing (similar procedure to Guimera et al. paper, 2010)
  spingc$modularity 	
  fwind[i, 6]<- spingc$modularity			# cannot work with unconnected graph, code not finalized!
  
  # Type of taxa
  fwind[i,7]<-Fraction.omnivores(fw0.p) 
  fwind[i,8]<-Level.omnivory(fw0.p)
  
  d.in<-degree(fw1.pi,mode="in")   					          	# in-degree of the trophospecies, indicates how many prey a taxon has
  d.out<-degree(fw1.pi,mode="out") 
  Prey<- (fwind[i,1]-sum(1*(d.out==0)))/fwind[i,1] 			# Fraction of taxa with at least one predator in the web
  Pred<- (fwind[i,1]-sum(1*(d.in==0)))/fwind[i,1] 			# Fraction of taxa with at least one resource in the web
  
  fwind[i, 9]<- Pred/Prey							                # Pred/Prey --> Vulnerability
  fwind[i, 10]<- sum(diag(fw1.p))/dim(fw1.p)[1]				# Fraction of cannibals
  fwind[i,11] <-sum(1*(d.in==0))/fwind[i,1]					  # proportion basal trophospecies 
  fwind[i,12]	<-sum(1*(d.out==0))/fwind[i,1]				  # proportion top trophospecies
  fwind[i,13]	<-1-fwind[i,11]-fwind[i,12]		
  fwind[i, 14]<- sd(colSums(fw1.p)/(sum(fw1.p)/dim(fw1.p)[1]))	# Standard deviation of mean generality
  fwind[i, 15]<- sd(rowSums(fw1.p)/(sum(fw1.p)/dim(fw1.p)[1]))	# Standard deviation of mean generality
  d.all<-degree(fw1.pi,mode="all") 	
  #fwind[i, 16]<- sd(d.all)		# Link density
  
  # Chain charcteristics							
  TrophInd(fw0.p)
  #fwind[i,16]<-mean(TrophInd(fw0.p)[,1])
  TL<-TLs(fw0.p) 						# Trophic level
  #TL[1]<-0							# Set detritus to 0	
  fwind[i,16]<-mean(TL[,1])
  fwind[i,17]<-mean(TL[,2])
  
  
  #fwind[i,18]<- average.path.length(fw1.pi)
  
  fwind[i,18]<- average.path.length(fw1.pi)
  paths<- shortest.paths(fw1.pi)		 				# calculates the length of the shortest path from or to vertices in the network
  #Longest.path<-Paths(fw0.p)[,2]
  fwind[i,19]<- mean(Paths(fw0.p)[,2])
  
}

fwind # Find metrics in the matrix fwind.p
write.table(fwind, "fw.emp.txt")

fwind.sel<-fwind[,c(1:16, 18)] #select metrics
fwind<-write.table(fwind.sel, "fw.output.txt")

#Plot metrics in space
# Code on the way


#draw food webs
#subsample adjacency matrices for the food webs
#Southern Baltic fw matrix
SB.ids<-which(poly.pa[1,]=="1")
SB.fw<-fw[SB.ids, SB.ids]
write.table(SB.fw,"SB_fw_mat.txt")
#Baltic Proper fw matrix
BP.ids<-which(poly.pa[2,]=="1")
BP.fw<-fw[BP.ids, BP.ids]
write.table(BP.fw,"BP_fw_mat.txt")
#Bothnian Sea fw matrix
BS.ids<-which(poly.pa[3,]=="1")
BS.fw<-fw[BS.ids, BS.ids]
write.table(BS.fw,"BS_fw_mat.txt")
#Bothnian Bay fw matrix
BB.ids<-which(poly.pa[4,]=="1")
BB.fw<-fw[BB.ids, BB.ids]
write.table(BB.fw,"BS_fw_mat.txt")


#make igraph edge.lists
SB.edge.list<-graph.adjacency(as.matrix(SB.fw))
BP.edge.list<-graph.adjacency(as.matrix(BP.fw))
BS.edge.list<-graph.adjacency(as.matrix(BS.fw))
BB.edge.list<-graph.adjacency(as.matrix(BB.fw))

### Food web Network Figures ###
# Centrality measures for all 4 regions
Cent.metrics.SB<- data.frame(
  deg.SB=degree(SB.edge.list),
  bet.SB=betweenness(SB.edge.list) ,
  eig.SB=evcent(SB.edge.list)$vector ,
  clo.SB=closeness(SB.edge.list)
)
Cent.metrics.BP<- data.frame(
  deg.BP=degree(BP.edge.list),
  bet.BP=betweenness(BP.edge.list) ,
  eig.BP=evcent(BP.edge.list)$vector ,
  clo.BP=closeness(BP.edge.list)
)
Cent.metrics.BS<- data.frame(
  deg.BS=degree(BS.edge.list),
  bet.BS=betweenness(BS.edge.list) ,
  eig.BS=evcent(BS.edge.list)$vector ,
  clo.BS=closeness(BS.edge.list)
)
Cent.metrics.BB<- data.frame(
  deg.BB=degree(BB.edge.list),
  bet.BB=betweenness(BB.edge.list) ,
  eig.BB=evcent(BB.edge.list)$vector ,
  clo.BB=closeness(BB.edge.list)
)

# Which taxa belong to which functional group
FG<-read.table("FG.txt")
FG.SB<- FG[which(FG[,1]%in%rownames(SB.fw)),]
FG.BP<- FG[which(FG[,1]%in%rownames(BP.fw)),] 
FG.BS<- FG[which(FG[,1]%in%rownames(BS.fw)),] 
FG.BB<- FG[which(FG[,1]%in%rownames(BB.fw)),] 

# Trophic level as a plotting parameter.
as.matrix(SB.fw)
TL.SB<-TLs(as.matrix(SB.fw))[,1] 
TL.BP<-TLs(as.matrix(BP.fw))[,1] 
TL.BS<-TLs(as.matrix(BS.fw))[,1] 
TL.BB<-TLs(as.matrix(BB.fw))[,1] 

# plot food webs according to TL
# creates a two-column matrix identifying the x and y values for each node.
layout.matrix.SB<-matrix(nrow=length(V(SB.edge.list)), ncol=2)
layout.matrix.SB[,1]<-runif(length(V(SB.edge.list))) # randomly assign along x-axis
layout.matrix.SB[,2]<-TL.SB # y-axis value based on trophic level

layout.matrix.BP<-matrix(nrow=length(V(BP.edge.list)), ncol=2)
layout.matrix.BP[,1]<-runif(length(V(BP.edge.list))) # randomly assign along x-axis
layout.matrix.BP[,2]<-TL.BP # y-axis value based on trophic level

layout.matrix.BS<-matrix(nrow=length(V(BS.edge.list)), ncol=2)
layout.matrix.BS[,1]<-runif(length(V(BS.edge.list))) # randomly assign along x-axis
layout.matrix.BS[,2]<-TL.BS # y-axis value based on trophic level

layout.matrix.BB<-matrix(nrow=length(V(BB.edge.list)), ncol=2)
layout.matrix.BB[,1]<-runif(length(V(BB.edge.list))) # randomly assign along x-axis
layout.matrix.BB[,2]<-TL.BB # y-axis value based on trophic level

match()
groups.SB<- FG.SB[,3][order(match(FG.SB[,1],colnames(as.matrix(SB.fw))))] # 1=non-living,2=phytopl, 3=zoopl, 4=benthos, 5=fish, 6=birds, 7=mammals
groups.BP<- FG.BP[,3][order(match(FG.BP[,1],colnames(as.matrix(BP.fw))))] # 1=non-organic,2=phytopl, 3=zoopl, 4=benthos, 5=fish, 6=birds, 7=mammals
groups.BS<- FG.BS[,3][order(match(FG.BS[,1],colnames(as.matrix(BS.fw))))] # 1=non-organic,2=phytopl, 3=zoopl, 4=benthos, 5=fish, 6=birds, 7=mammals
groups.BB<- FG.BB[,3][order(match(FG.BB[,1],colnames(as.matrix(BB.fw))))] # 1=non-organic,2=phytopl, 3=zoopl, 4=benthos, 5=fish, 6=birds, 7=mammals

colbar.FG<- c("darkgrey","lightgreen","cyan","orange","blue","magenta","lightpink")

par(mai=c(0,0,0,0))
par(mar=c(.5,.5,.5,.5),mfrow=c(1,1))
# Plot
plot(SB.edge.list, layout=layout.matrix.SB, vertex.size = 1/5*degree(SB.edge.list)+4, vertex.frame.color = "gray50", edge.arrow.size=F,  vertex.color=colbar.FG[groups.SB],vertex.label = NA, edge.color="grey40", edge.width=0.15)
plot(BP.edge.list, layout=layout.matrix.BP, vertex.size = 1/5*degree(SB.edge.list)+4, vertex.frame.color = "gray50", edge.arrow.size=F,  vertex.color=colbar.FG[groups.BP],vertex.label = NA, edge.color="grey40", edge.width=0.15)
plot(BS.edge.list, layout=layout.matrix.BS, vertex.size = 1/5*degree(BS.edge.list)+4, vertex.frame.color = "gray50", edge.arrow.size=F,  vertex.color=colbar.FG[groups.BS],vertex.label = NA, edge.color="grey40", edge.width=0.15)
plot(BB.edge.list, layout=layout.matrix.BB, vertex.size = 1/5*degree(BB.edge.list)+4, vertex.frame.color = "gray50", edge.arrow.size=F,  vertex.color=colbar.FG[groups.BB],vertex.label = NA, edge.color="grey40", edge.width=0.15)






