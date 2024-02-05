Paths<-function(web){

basal <- rownames(subset(web, apply(web, 2, sum)==0) & apply(web, 1, sum)!= 0)
basaedge.list_web <- graph.adjacency(web, mode = "directed");
paths_prey <- shortest.paths(graph = edge.list_web, v= V(edge.list_web),
to = V(edge.list_web)[basal], mode = "in", weights = NULL, algorithm = "unweighted")
paths_prey[is.infinite(paths_prey)] <- NA
shortest_paths <- as.matrix(apply(paths_prey, 1, min, na.rm=TRUE))
longest_paths <- as.matrix(apply(paths_prey, 1, max, na.rm=TRUE))
meanShort<-shortest_paths
meanLong<-longest_paths

PATHS<-cbind(meanShort, meanLong)

PATHS
}

