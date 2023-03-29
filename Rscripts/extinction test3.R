######## gammla kod från extinction tests ###############


#lägst till högsta
order_degree.low<-as.numeric(names(sort(deg.sp)))
#remove basal nodes from extinction order: autotrophs=2, mixotrophs=3, detritus=15
order_degree.low<-setdiff(order_degree.low, remove)

degree.extinc.low<-ExtinctionOrder(fw.1981, Order = order_degree.low, NetworkType = "Trophic")
ExtinctionPlot(History = degree.extinc.low$sims)


### topological roles extinctions ###

uwroles81<-uw.top.roles[[1]]

uwroles81[,4]<-cbind(1:25)

test<-sort(uwroles81[2:4], decreasing = T)
order_among.mod<-as.numeric(unlist(test[3]))

order_among.mod<-setdiff(order_among.mod, remove)

#among modle degree högsta till lägsta
among.mod.ext<-ExtinctionOrder(fw.1981, Order= order_among.mod, NetworkType = "Trophic")
ExtinctionPlot(History = among.mod.ext$sims)

#within module degree högst till lägsta
test<-sort(uwroles81[1:4], decreasing = T)
order_within.mod<-as.numeric(unlist(test[4]))
order_within.mod<-setdiff(order_within.mod, remove)

within.mod.ext<-ExtinctionOrder(fw.1981, Order = order_within.mod, NetworkType = "Trophic")
ExtinctionPlot(History = within.mod.ext$sims)

#null model extinction sim, random extinctions
rand<-RandomExtinctions(fw.1981, nsim = 1000, NetworkType = "Trophic")

CompareExtinctions(rand, degree.extinc)#order= högsta degree till lägsta
#CompareExtinctions(rand, degree.extinc.low)#Order= lägsta degree till högsta
CompareExtinctions(rand, among.mod.ext)#Order= högsta among module degree till lägsta
CompareExtinctions(rand, within.mod.ext)#Order= högsta within module degree till lägsta

###### se vad som händer efter varje steg?





############################################
#weighted sum of weights, högst till lägsta

fw_81<-w_fw[[1]]
fw_81
colnames(fw_81)<-c(1:25)
rownames(fw_81)<-c(1:25)
fw_81

fw81<-graph_from_adjacency_matrix(
  fw_81,
  mode = c("directed"),
  weighted = TRUE,
)
deg.sp<-degree(fw81,mode = c("all"),loops = TRUE,normalized = FALSE)
sum.link.weights.sp<-apply(fw_81, 1, sum)+apply(fw_81, 2, sum)


order_sum_link<-as.numeric(names(sort(sum.link.weights.sp, decreasing = T)))
order_sum_link<-setdiff(order_sum_link, remove)

sum.link.ext<-ExtinctionOrder(fw_81, Order =order_sum_link, NetworkType = "Trophic")
ExtinctionPlot(History = sum.link.ext$sims)

#sum.link.ext<-ExtinctionOrder(fw_81, Order =order_sum_link, NetworkType = "Trophic", IS=0.5)
#ExtinctionPlot(History = sum.link.ext$sims)

rand.w<-RandomExtinctions(fw_81,nsim = 100, NetworkType = "Trophic")
CompareExtinctions(rand.w, sum.link.ext)

#weighted among mod
wroles81<-w.top.roles[[1]]
wroles81[,4]<-cbind(1:25)

test<-sort(wroles81[2:4], decreasing = T)
w.order_among.mod<-as.numeric(unlist(test[3]))
w.order_among.mod<-setdiff(order_among.mod, remove)

w.among.mod.ext<-ExtinctionOrder(fw_81, Order = w.order_among.mod, NetworkType = "Trophic")
ExtinctionPlot(History = w.among.mod.ext$sims)

CompareExtinctions(rand.w, w.among.mod.ext)