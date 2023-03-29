TLs<-function(web){

library(NetIndices)
library(igraph)  

basal <- rownames(subset(web, apply(web, 2, sum)==0) & apply(web, 1, sum)!= 0)
edge.list_web <- graph.adjacency(web, mode = "directed");
paths_prey <- shortest.paths(graph = edge.list_web, v= V(edge.list_web),
to = V(edge.list_web)[basal], mode = "in", weights = NULL, algorithm = "unweighted")
paths_prey[is.infinite(paths_prey)] <- NA
shortest_paths <- as.matrix(apply(paths_prey, 1, min, na.rm=TRUE))
longest_paths <- as.matrix(apply(paths_prey, 1, max, na.rm=TRUE))
in_deg <- apply(web, 2, sum) 				## Take each columns, sum  the rows
out_deg <- apply(web, 1, sum) 			## Take each rows, sum the columns


# Shortest TL
sTL <- 1 + shortest_paths  # Commonly, detritus have a TL value of 1. (Shortest path to basal = 0)
# Longest TL
lTL <- 1 + longest_paths

S<-dim(web)[1]

# Creating the matrix  
short_TL_matrix <- matrix(NA, S, S)
long_TL_matrix <- matrix(NA, S, S)
prey_ave <- rep(NA, S) # averaged prey
chain_ave<-rep(NA, S)
#
for(j in 1:S){
  for(i in 1:S){
    lij <- web[i, j] 					# is the interaction matrix
    prey_ave[j] <- 1 / in_deg[j] 			# if Bas species in, no in-degrees ; re-attribute a value of 0
    short_TL_matrix[i,j] <- lij * sTL[i]   	# Shortest path in matrix
	long_TL_matrix[i,j] <- lij * lTL[i]
  }  
}

prey_ave[which(prey_ave == Inf)] <- 0

sumShortTL <- as.matrix(apply(short_TL_matrix, 2, sum)) # sum all shortest path
sumLongTL<- as.matrix(apply(long_TL_matrix, 2, sum)) # sum all shortest path

# Short-weighted TL
# weigth by the number of prey

SWTL <- matrix(data = NA, nrow = S, ncol = 1)
for(i in 1:S){
  SWTL[i] <- 1 + (prey_ave[i] * sumShortTL[i])  
}

LWTL <- matrix(data = NA, nrow = S, ncol = 1)
for(i in 1:S){
  LWTL[i] <- 1 + (prey_ave[i] * sumLongTL[i])  
}


TL_NI<-TrophInd(web)
TL_NetInd<- TL_NI[,1]

TLS<-cbind(SWTL, LWTL, TL_NetInd)
colnames(TLS)<-c("SWTL", "LWTL", "NetInd_TL")
rownames(TLS)<-rownames(web)

TLS

}





