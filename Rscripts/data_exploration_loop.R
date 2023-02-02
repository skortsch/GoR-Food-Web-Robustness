
#data exploration loop
setwd("D:/R_Git/GoR-Food-webs/Rscripts")
w_fw<-readRDS("../Data/weighted_foodwebs.rds")

library(igraph)

par(mfrow=c(4,2))
for(i in 1:length(w_fw)){
  fw_i<-w_fw[[i]]
  
  fw<-graph_from_adjacency_matrix(
    fw_i,
    mode = c("directed"),
    weighted = T
  )
  deg.sp<-degree(fw,mode = c("all"),loops = TRUE,normalized = FALSE)
  sum.link.weights.sp<-apply(fw_i, 1, sum)+apply(fw_i, 2, sum)
  names(sum.link.weights.sp)
  names(deg.sp)
  
  year<-c("1981", "1986", "1991", "1996", "2001", "2006", "2011")
  
  plot(deg.sp[-c(2:4)], log(sum.link.weights.sp+1)[-c(2:4)], xlab="degree", ylab="sum.link.weights", main = year[i])

}
