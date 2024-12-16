
library("igraph")
#https://www.data-imaginist.com/2017/introducing-tidygraph/

#Do pollinator degree and sum of link weights correlate?

#import data
w_fw<-readRDS("../Data/weighted_foodwebs.rds")

fw_81<-w_fw[[1]]
  
fw81<-graph_from_adjacency_matrix(
  fw_81,
  mode = c("directed"),
  weighted = TRUE,
)

deg.sp<-degree(fw81,mode = c("all"),loops = TRUE,normalized = FALSE)

sum.link.weights.sp<-apply(fw_81, 1, sum)+apply(fw_81, 2, sum)
names(sum.link.weights.sp)
names(deg.sp)

hist(sum.link.weights.sp)

plot(deg.sp, log(sum.link.weights.sp+1))

#
fw_86<-w_fw[[2]]

fw86<-graph_from_adjacency_matrix(
  fw_86,
  mode = c("directed"),
  weighted = TRUE,
)

deg.sp<-degree(fw86,mode = c("all"),loops = TRUE,normalized = FALSE)

sum.link.weights.sp<-apply(fw_86, 1, sum)+apply(fw_86, 2, sum)
names(sum.link.weights.sp)
names(deg.sp)

plot(deg.sp[-c(2:4)], log(sum.link.weights.sp+1)[-c(2:4)])

#
fw_2011<-w_fw[[7]]

fw2011<-graph_from_adjacency_matrix(
  fw_2011,
  mode = c("directed"),
  weighted = TRUE,
)

deg.sp<-degree(fw2011,mode = c("all"),loops = TRUE,normalized = FALSE)

sum.link.weights.sp<-apply(fw_2011, 1, sum)+apply(fw_2011, 2, sum)
names(sum.link.weights.sp)
names(deg.sp)

plot(deg.sp[-c(2:4)], log(sum.link.weights.sp+1)[-c(2:4)], xlab="degree", ylab="sum.link.weights")

