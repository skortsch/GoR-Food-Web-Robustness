
setwd("C:/LocalData/susakort/GitHub/GoR-Food-webs/Data")

#load data
suffix<-"R1000"
insuffix<-"R1000"
load("GoRoff_v2.Rdata")
load(paste0("ConstSW5_GoRoff_", suffix, ".Rdata"))
load(paste0("GoRoff_FWTSConstS_",insuffix, ".Rdata")) # food web metrics for the 1000 food webs


##load Food web meta web data
load("metawebs_GoR.RData")
load("uw_fw_GoR.RData")
load("w_fw_GoR.RData")


index<-matrix(NA, 34, 1000)

index.arr <- list() #could be list of output for the 1000 food webs per year

for (i in 1:1000){
  #i<-1
  sp.pa.year<-1*((t(biomC[,,i])[, match(row.names(un_mw), row.names(biomC[,,1]))])>0) #reorder to use for making networks (1000 per year)
  
  for (j in 1:34){
  fw.sel<-which(sp.pa.year[j,]==1)  
  fw<-un_mw[fw.sel, fw.sel]
  
  #WHAT METRIC DO YOU WANT TO CALCULATE ADD HERE
  #Example with number of species per network (unweighted)
  
  index[j,i]<-dim(fw)[1] #number of species
  

  #YOU NEED TO DO THE SAME FOR WEIGHTED NETWORKS
  #YOU CAN USE THE BIOMASS INFO FOR THIS AND THE UNEWEIGHTED fw
  #INSTEAD OF MATRICES YOU CAN ASLO STORE YOUR OUTPUT IN LISTS
  #THESE LISTS WE CAN ADD TO ARRAYS AND THEN CALCULATE THE MEAN VALUES
  
  
  }
}

index.median<- apply(index, 1, median)

#To make the array
years<-34
dim.index<-5
no.fw<-1000

arr <- array(unlist(index), dim = c(years, dim.index, no.fw))
dimnames(arr) <- list(dimnames(biomC)[[2]],
                      colnames(dim.index), 
                      1:dim(biomC)[3])

#calc median value of array
apply(arr, 1:2, median, na.rm = TRUE)




