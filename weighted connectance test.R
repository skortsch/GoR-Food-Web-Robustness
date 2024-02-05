####understanding weighted connectance####
library(igraph)

# 5x5 test-matrix
mat<-matrix(0:1, ncol = 5, nrow = 5)
mat[,1]<-0
mat[,5]<-c(0,1,1,1,0)
mat[,2]<-c(1,0,0,0,0)
mat[,4]<-c(0,1,1,0,0)
mat[,3]<-c(1,1,0,0,0)
fw<-graph_from_adjacency_matrix(mat, mode = "directed")
plot(fw,edge.width = E(fw)$weight*1.2, edge.arrow.size=0.5,layout=layout.circle(fw))

w.mat<-mat

w.mat[1,]<-c(0,10,7,0,0)
w.mat[2,]<-c(0,0,6,3,3)
w.mat[3,]<-c(0,0,0,1,1)
w.mat
wfw<-graph_from_adjacency_matrix(w.mat, mode = "directed", weighted = T)


plot(wfw, edge.width = E(wfw)$weight*1.2, edge.arrow.size=0.5,layout=layout.circle(fw))

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

fluxind(mat)$qC.w #weighted connectance output, should be the "weighted weighted" connectance right?
fluxind(mat)$qC.uw
#fluxind(mat2)$qC.w # samma som bin w.con pga shannon samma då alla är likavärda
fluxind(w.mat)$qC.w# #weighted connectance output, weighted weighted connectance right? lägre w. con för att mindre jämnt fördelat.
fluxind(w.mat)$qC.uw


w.mat2<-w.mat
w.mat2[w.mat2>=10]<-24
wfw2<-graph_from_adjacency_matrix(w.mat2, mode = "directed", weighted = T)
plot(wfw2, edge.width = E(wfw2)$weight, edge.arrow.size=0.5, layout=layout.circle(fw))
fluxind(w.mat2)$qC.w #blir mindre då det är  mera skewed

w.mat3<-w.mat2
w.mat3[w.mat3>=10]<-50
wfw3<-graph_from_adjacency_matrix(w.mat3, mode = "directed", weighted = T)
plot(wfw3, edge.width = E(wfw3)$weight, edge.arrow.size=0.5, layout=layout.circle(fw))
fluxind(w.mat3)$qC.w #blir mindre då det är ännu mera skewed


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

