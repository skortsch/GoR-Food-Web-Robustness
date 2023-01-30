
#working dir
setwd("C:/LocalData/susakort/supervision/Patrick Ståhl/data")

#load data
#biom.data<-read.table("clipboard", dec=",", header=T)
#biom.dat<-as.data.frame(t(biom.data))
#write.csv(biom.dat, "biom.dat.csv", row.names=F)

w_fw<-readRDS("weighted_foodwebs.rds")

#biomass data yearly
biom.data<-read.csv("biom.dat.csv")
rownames(biom.data)<-c(1981, 1986, 1991, 1996, 2001, 2006, 2011)

#binary randomisations time by species matrix
mods_curve<-readRDS("model_curveball.rds")

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