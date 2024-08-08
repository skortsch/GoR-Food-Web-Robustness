####understanding weighted connectance and AMI and relative ascendency####
library(igraph)
library(NetIndices)

# 5x5 test-matrix
mat<-matrix(ncol = 5, nrow = 5)
mat[,1]<-0
mat[,5]<-c(0,1,1,1,0)
mat[,2]<-c(1,0,0,0,0)
mat[,4]<-c(0,1,1,0,0)
mat[,3]<-c(1,1,0,0,0)
fw<-graph_from_adjacency_matrix(mat, mode = "directed")
plot(fw,edge.width = E(fw)$weight*1.2, edge.arrow.size=1,layout=layout.circle(fw))


#weighted test matrix
w.mat<-mat
w.mat[1,]<-c(0,10,7,0,0)
w.mat[2,]<-c(0,0,6,3,3)
w.mat[3,]<-c(0,0,0,1,1)
w.mat
wfw<-graph_from_adjacency_matrix(w.mat, mode = "directed", weighted = T)
plot(wfw, edge.width = E(wfw)$weight*1.2, edge.arrow.size=1,layout=layout.circle(fw))

#fluxind function from Susanne:
### Quantitaive unweighted and weighted link density and directed connectance ###
# written by Susanne Kortsch 27.08.2019

# Method for calculating weighted directed connectance (Cw) after Bersier et al. 2002 (and Banasek-Richter et al.(2009))
# The method is based on information-theoretical indices, ex. the Shannon measure of entropy (or uncertainity).
# For a given number x of events, H reaches its maximum when all events occur in equal proportion (H=log x)

# A food web with S species can be represented by an S-by-S quantitative food-web matrix b = [bij]; 
# the value of the element bij means the amount of biomass passing from taxon i to taxon j per unit area and time; 
# for taxon k, we can measure the diversity of the biomass (H) coming from its resources (R,k) 
# and of that going to its consumers (C,k), i.e. the taxon-specific Shannon indices of inflows and outflows.

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
  
  res$H.in<-H.in # added for easy output by P.S, not part of original function
  res$H.out<-H.out # added for easy output by P.S
  res$N.con<-N.con # added for easy output by P.S
  res$N.res<-N.res # added for easy output by P.S
  return(res)
}

####dissecting qLDw to understand weighted con####
fluxind(w.mat)$H.in #shannon in div for species 1-5
fluxind(w.mat)$H.out #shannon out
fluxind(w.mat)$N.con #consumers
fluxind(w.mat)$N.res #resourses
fluxind(w.mat)$qLD.w

fluxind(mat)$H.in #shannon in div
fluxind(mat)$H.out # shannon out
fluxind(mat)$N.con #effecitve consumers
fluxind(mat)$N.res #effective resources
fluxind(mat)$qLD.w

fluxind(mat)$qC.w #weighted connectance output, should be the "weighted weighted" connectance right
fluxind(mat)$qC.uw

fluxind(w.mat)$qC.w# #weighted connectance output, weighted weighted connectance right. lower w.con because more skewd
fluxind(w.mat)$qC.uw

####testing different matrices for w.con####
w.mat2<-w.mat
w.mat2[w.mat2>=10]<-24 #increase skew
wfw2<-graph_from_adjacency_matrix(w.mat2, mode = "directed", weighted = T)
plot(wfw2, edge.width = E(wfw2)$weight, edge.arrow.size=0.5, layout=layout.circle(fw))
fluxind(w.mat2)$qC.w#lower when more skewed

hist(w.mat2, breaks=100)

w.mat3<-w.mat2
w.mat3[w.mat3>=10]<-50
wfw3<-graph_from_adjacency_matrix(w.mat3, mode = "directed", weighted = T)
plot(wfw3, edge.width = E(wfw3)$weight, edge.arrow.size=0.5, layout=layout.circle(fw))
fluxind(w.mat3)$qC.w #lower when more skewed

hist(w.mat3, breaks=100)

####AMI and A/C####
UncInd() #function to to get AMI, output is a list, use $AMI to get value
AscInd() #function to get relative ascendancy, output is a table, row 1, col 4 [1,4] is the A/C


mat #unweighted version
UncInd(mat)$AMI
UncInd(mat)$AMI/UncInd(mat)$HR #AMI  (Normalized as in Cannin et al. with max uncertainty)
AscInd(mat)[1,4]# A/C
fluxind(mat)$qC.w


w.mat
wfw<-graph_from_adjacency_matrix(w.mat, mode = "directed", weighted = T)
plot(wfw, edge.width = E(wfw)$weight*1.2, edge.arrow.size=1,layout=layout.circle(fw))
UncInd(w.mat)$AMI
UncInd(w.mat)$AMI/UncInd(w.mat)$HR #AMI (Normalized as in Canning et al. with max uncertainty)
AscInd(w.mat)[1,4]#A/C
fluxind(w.mat)$qC.w

#AMI goes up a bit, A/C also goes up from last web
w.mat2<-w.mat
wfw2<-graph_from_adjacency_matrix(w.mat2, mode = "directed", weighted = T)
plot(wfw2, edge.width = E(wfw2)$weight*1.2, edge.arrow.size=1,layout=layout.circle(fw))
w.mat2[w.mat2>=10]<-24
UncInd(w.mat2)$AMI
UncInd(w.mat2)$AMI/UncInd(w.mat2)$HR #
AscInd(w.mat2)[1,4]#A/C
fluxind(w.mat2)$qC.w

#AMI goes down a bit when skew increases, A/C goes up a bit from last web
w.mat3<-w.mat2
w.mat3[w.mat3>=10]<-50
wfw3<-graph_from_adjacency_matrix(w.mat3, mode = "directed", weighted = T)
plot(wfw3, edge.width = E(wfw3)$weight*1.2, edge.arrow.size=1,layout=layout.circle(fw))
UncInd(w.mat3)$AMI
UncInd(w.mat3)$AMI/UncInd(w.mat3)$HR #
AscInd(w.mat3)[1,4]#A/C
fluxind(w.mat3)$qC.w


#AMI goes down when more flows are stronger, A/C also drops from last web
w.mat4<-w.mat3
w.mat4[1,3]<-25
wfw4<-graph_from_adjacency_matrix(w.mat4, mode = "directed", weighted = T)
plot(wfw4, edge.width = E(wfw4)$weight, edge.arrow.size=1, layout=layout.circle(fw))
UncInd(w.mat4)$AMI
UncInd(w.mat4)$AMI/UncInd(w.mat4)$HR #
AscInd(w.mat4)[1,4]#A/C
fluxind(w.mat4)$qC.w



#similar matrix to first w.mat, AMI is slightly higher, A/C slightly higher, w.con slightly/ barley lower
w.mat5<-matrix(c(0,12,9,0,0,
                 0,0,5,4,3,
                 0,0,0,2,1,
                 0,0,0,0,1,
                 0,0,0,0,0),5,5,byrow = T)
wfw5<-graph_from_adjacency_matrix(w.mat5, mode = "directed", weighted = T)
plot(wfw5, edge.width = E(wfw5)$weight, edge.arrow.size=1, layout=layout.circle(fw))

UncInd(w.mat5)$AMI
UncInd(w.mat5)$AMI/UncInd(w.mat5)$HR #
AscInd(w.mat5)[1,4]#A/C
fluxind(w.mat5)$qC.w


#AMI is high but can still go higher, but A/C is maxed... (the normalized AMI is maxed out at 1)
#The "web" is now constrained to a chain, so no uncertainty where energy flows, just how much that flows?
w.mat6<-matrix(c(0,12,0,0,0,
                 0,0,6,0,0,
                 0,0,0,3,0,
                 0,0,0,0,1.5,
                 0,0,0,0,0),5,5,byrow = T)
wfw6<-graph_from_adjacency_matrix(w.mat6, mode = "directed", weighted = T)
plot(wfw6, edge.width = E(wfw6)$weight, edge.arrow.size=1, layout=layout.circle(fw))

UncInd(w.mat6)$AMI
UncInd(w.mat6)$AMI/UncInd(w.mat6)$HR #
AscInd(w.mat6)[1,4]#A/C
fluxind(w.mat6)$qC.w

#gives MAX AMI and A/C
#The "web" is now constrained to a chain, there is no uncertainty in the flows, come from one goes to next, and all flows are the same
w.mat7<-matrix(c(0,12,0,0,0,
                 0,0,12,0,0,
                 0,0,0,12,0,
                 0,0,0,0,12,
                 0,0,0,0,0),5,5,byrow = T)
wfw7<-graph_from_adjacency_matrix(w.mat7, mode = "directed", weighted = T)
plot(wfw7, edge.width = E(wfw7)$weight, edge.arrow.size=1, layout=layout.circle(fw))
UncInd(w.mat7)$AMI
UncInd(w.mat7)$AMI/UncInd(w.mat7)$HR #
AscInd(w.mat7)[1,4]#A/C
fluxind(w.mat7)$qC.w


#lowest AMI and A/C that i can get without loops included
w.mat8<-matrix(c(0,12,12,12,12,
                 12,0,12,12,12,
                 12,12,0,12,12,
                 12,12,12,0,12,
                 12,12,12,12,0),5,5,byrow = T)
wfw8<-graph_from_adjacency_matrix(w.mat8, mode = "directed", weighted = T)
plot(wfw8, edge.width = E(wfw8)$weight, edge.arrow.size=1, layout=layout.circle(fw))
UncInd(w.mat8)$AMI
UncInd(w.mat8)$AMI/UncInd(w.mat8)$HR #
AscInd(w.mat8)[1,4]#A/C
fluxind(w.mat8)$qC.w


#minimal AMI and A/C, loops includded
#there is an equal probability of flow moving from one taxon to any other taxon and the uncertainty in the direction of flow is greatest?
w.mat9<-matrix(c(12,12,12,12,12,
                 12,12,12,12,12,
                 12,12,12,12,12,
                 12,12,12,12,12,
                 12,12,12,12,12),5,5,byrow = T)
wfw9<-graph_from_adjacency_matrix(w.mat9, mode = "directed", weighted = T)
plot(wfw9, edge.width = E(wfw9)$weight, edge.arrow.size=1, layout=layout.circle(fw))
UncInd(w.mat9)$AMI
UncInd(w.mat9)$AMI/UncInd(w.mat9)$HR #
AscInd(w.mat9)[1,4]# A/C
fluxind(w.mat9)$qC.w


#very low AMI and A/C
#probability of flow moving from one taxon to any other taxon is no longer the same for all, since some flows are weaker?
#but it is still relatively chaotic/uncertain? 
w.mat10<-matrix(seq(25,1),5,5,byrow = T)
wfw10<-graph_from_adjacency_matrix(w.mat10, mode = "directed", weighted = T)
plot(wfw10, edge.width = E(wfw10)$weight, edge.arrow.size=1, layout=layout.circle(fw))

UncInd(w.mat10)$AMI
UncInd(w.mat10)$AMI/UncInd(w.mat10)$HR
AscInd(w.mat10)[1,4]# A/C
fluxind(w.mat10)$qC.w

