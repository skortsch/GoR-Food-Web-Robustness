plotfw<-function(net, col=NULL, lab=NULL, size=NULL,
                 nylevel=5, maxsize=10, labcex=0.01,
                 ynum=6, ylab= "Trophic Level", ...){
  n <- vcount(net)
  if (!is.null(col)){
    V(net)$color <- col
  } else{
    V(net)$color <- rainbow(vcount(net))
  }
  if (!is.null(lab)){
    V(net)$name <- lab
  }
  if (!is.null(size)){
    V(net)$size <- size
  } else {
    V(net)$size <- maxsize/2
  }
  
  tl <- trophiclevels(net)
  dgpred <- tl
  
  bks <- c(0.9, seq(1.9, max(tl), length.out = nylevel))
  ynod <- cut(tl, breaks = bks, include.lowest = TRUE, 
              labels = 1:(length(bks)-1))
  xnod <- rep(0,n)
  for (i in 1:(length(bks)-1)){
    l <- sum(ynod==i)
    xnod[ynod==i] <- seq(-l,l,length.out = l)[order(tl[ynod==i])]
  }
  coo <- cbind(xnod,tl)
  
  #y axis with 1 and continuous axis from 2 to max TL.
  yax <- c(1, seq(2, max(tl), length.out = ynum-1))
  labax <- round(yax, 1)
  #rescale xax between -1 and 1
  laby <- (yax-min(yax))/(max(yax)-min(yax))*2 - 1
  
  plot(net, layout= coo, vertex.label.color="black", 
       vertex.label.cex=labcex, ...)
  axis(2, at = laby, labels = labax)
  mtext(ylab, side = 2, line=2)
  res <- data.frame(coo, "size"= V(net)$size, "color"= V(net)$color)
  names(res) <- c("x", "y", "size", "color")
  row.names(res) <- V(net)$name
  invisible(res)
}

trophiclevels<-function(net){
  mat <- get.adjacency(net, sparse=F)
  edge.list_web <- graph.adjacency(mat, mode = "directed");
  
  basal <- rownames(subset(mat, apply(mat, 2, sum)==0) & apply(mat, 1, sum)!= 0)
  paths_prey <- suppressWarnings(shortest.paths(graph = net, v= V(net), to = V(net)[basal], 
                                                mode = "in", weights = NULL, algorithm = "unweighted"))
  
  paths_prey[is.infinite(paths_prey)] <- NA
  shortest_paths <- suppressWarnings(as.matrix(apply(paths_prey, 1, min, na.rm=TRUE)))
  #for species with no prey apart of them
  shortest_paths[is.infinite(shortest_paths)] <- NA
  
  in_deg <- apply(mat, 2, sum) #==degree(net, mode = "in")
  # Shortest TL
  sTL <- 1 + shortest_paths  # Commonly, detritus have a TL value of 1. (Shortest path to basal = 0)
  
  S<-dim(mat)[1] #== vcount(net)
  # Creating the matrix  
  short_TL_matrix <- mat*matrix(rep(sTL, length(sTL)), ncol=length(sTL))
  
  prey_ave <- ifelse(in_deg==0, 0, 1 / in_deg)
  
  sumShortTL <- apply(short_TL_matrix, 2, sum, na.rm=T) # sum all shortest path
  
  # Short-weighted TL weigth by the number of prey
  SWTL <- 1 +(prey_ave * sumShortTL)
  
  #check that only basal species have TL of 1
  SWTL[!rownames(mat)%in%basal & SWTL == 1] <- NA
  
  #res<-data.frame("swTL"= SWTL, "lwTL"=LWTL)
  #rownames(res)<-rownames(mat)
  return(SWTL)
}
