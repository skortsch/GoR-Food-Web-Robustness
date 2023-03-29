#### Food webs along gardients in the Batic Sea ####

#set working directory
setwd("D:/BONUS BLUEWEBS project/R scripts food webs")

# Load packages
library(igraph)		# Network analysis
library(NetIndices)	# Trophic level and omnivory 

source("TL function.R") # 
source("SWTL Omnivory.R")
source("Paths function.R")

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
sp.SB<-which(poly.pa[1,]==1) 				# sp.p a
fw.SB<-fw[sp.SB,sp.SB]

#Baltic Proper
sp.BP<-which(poly.pa[2,]==1) 				# sp.p a
fw.BP<-fw[sp.BP,sp.BP]

#Bothnian Sea
sp.BS<-which(poly.pa[2,]==1) 				# sp.p a
fw.BS<-fw[sp.BS,sp.BS]

#Bothnian Bay
sp.BB<-which(poly.pa[2,]==1) 				# sp.p a
fw.BB<-fw[sp.BB,sp.BB]


#food web metrics

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
  #spingc <- spinglass.community(fw1.pi) 		# Compartmentalisation estimation based on simulated annealing (similar procedure to Guimera et al. paper, 2010)
  #spingc$modularity 	
  #fwind[i, 6]<- spingc$modularity			# cannot work with unconnected graph, code not finalized!
  
  # Type of taxa
  #wind[i,7]<-Fraction.omnivores(fw0.p) 
  #fwind[i,8]<-Level.omnivory(fw0.p)
  
  d.in<-degree(fw1.pi,mode="in")   						# in-degree of the trophospecies, indicates how many prey a taxon has
  d.out<-degree(fw1.pi,mode="out") 
  Prey<- (fwind[i,1]-sum(1*(d.out==0)))/fwind[i,1] 			# Fraction of taxa with at least one predator in the web
  Pred<- (fwind[i,1]-sum(1*(d.in==0)))/fwind[i,1] 			# Fraction of taxa with at least one resource in the web
  
  fwind[i, 9]<- Pred/Prey							# Pred/Prey --> Vulnerability
  fwind[i, 10]<- sum(diag(fw1.p))/dim(fw1.p)[1]				# Fraction of cannibals
  fwind[i,11] <-sum(1*(d.in==0))/fwind[i,1]					# proportion basal trophospecies 
  fwind[i,12]	<-sum(1*(d.out==0))/fwind[i,1]				# proportion top trophospecies
  fwind[i,13]	<-1-fwind[i,11]-fwind[i,12]		
  fwind[i, 14]<- sd(colSums(fw1.p)/(sum(fw1.p)/dim(fw1.p)[1]))	# Standard deviation of mean generality
  fwind[i, 15]<- sd(rowSums(fw1.p)/(sum(fw1.p)/dim(fw1.p)[1]))	# Standard deviation of mean generality
  d.all<-degree(fw1.pi,mode="all") 	
  #fwind[i, 16]<- sd(d.all)		# Link density
  
  # Chain charcteristics							
  TrophInd(fw0.p)
  fwind[i,16]<-mean(TrophInd(fw0.p)[,1])
  #TL<-TLs(fw0.p) 						# Trophic level
  #TL[1]<-0							# Set detritus to 0	
  #fwind[i,16]<-mean(TL[,1])
  #fwind[i,17]<-mean(TL[,2])
  
  
  fwind[i,18]<- average.path.length(fw1.pi)
  
  paths<- shortest.paths(fw1.pi)		 	# calculates the length of the shortest path from or to vertices in the network
  #fwind[i,19]<- mean(Paths(fw0.p)[,1])
  #fwind[i,19]<- mean(Paths(fw0.p)[,2])
  
}

fwind # Find metrics in the matrix fwind.p


