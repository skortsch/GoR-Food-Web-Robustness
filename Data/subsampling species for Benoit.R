###Script to subsample species list for Benoit###
mean_biom<-matrix(0, 33, 34)
rownames(mean_biom)<-rownames(biomC[,,1000])
colnames(mean_biom)<-colnames(biomC[,,1000])
median_biom<-matrix(0, 33, 34)
rownames(median_biom)<-rownames(biomC[,,1000])
colnames(median_biom)<-colnames(biomC[,,1000])
sd_biom<-matrix(0, 33, 34)
rownames(sd_biom)<-rownames(biomC[,,1000])
colnames(sd_biom)<-colnames(biomC[,,1000])

for (i in 1:33){
  for (j in 1:34){
mean_biom[i,j]<-mean(biomC[i,j,])
median_biom[i,j]<-median(biomC[i,j,])
sd_biom[i,j]<-sd(biomC[i,j,])
  }
}

