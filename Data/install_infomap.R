# Install  (if not installed) and load necessary packages
package.list=c("attempt", "cowplot", "igraph", "ggalluvial","magrittr","metafolio","tidyverse","vegan", "devtools")
loaded <-  package.list %in% .packages()
package.list <-  package.list[!loaded]
installed <-  package.list %in% .packages(TRUE)
if (!all(installed)) install.packages(package.list[!installed], repos="http://cran.rstudio.com/")

# Install infomapecology 
devtools::install_github('Ecological-Complexity-Lab/infomap_ecology_package', force=T)

library(infomapecology)

# Check the version.
packageDescription('infomapecology')

setwd('where your Infomap file and R script will live')
install_infomap()

# Check Infomap is running
setwd('where your Infomap file and R script now live')
check_infomap() # Make sure file can be run correctly. Should return TRUE



### test ####

a<-w.webs.1000[[1]] #i
fw<-a[[1]] #j
mon.layer<-create_monolayer_object(fw, directed = T, bipartite = F)

t<-run_infomap_monolayer(mon.layer, two_level = T, seed=123
                         ,flow_model = "directed", trials = 100, silent = T)

t$L #should be the modularity
print(t$modules, n=dim(fw)[1])# module1 col is which module a sp is in.

print(t$edge_list, n=dim(t$edge_list)[1])


####loop for weighted modularity####

mod<-matrix(NA,34,1000)
for (i in 1:1000){
  a<-w.webs.1000[[i]]
  for(j in 1:34){
    fw<-a[[j]]
    mon.layer<-create_monolayer_object(fw, directed = T, bipartite = F)
    output<-run_infomap_monolayer(mon.layer, two_level = T, seed=123
                                  ,flow_model = "directed", trials = 100,
                                  silent = T, verbose = F)
    mod[j,i]<-output$L
  }
}
mod
plot(apply(mod,1, median))

plot(apply(mod,1,median), apply(w.R50th30,1,median))
abline(lm(apply(w.R50th30,1,median)~apply(mod,1,median)))
summary(lm(apply(w.R50th30,1,median)~apply(mod,1,median)))
