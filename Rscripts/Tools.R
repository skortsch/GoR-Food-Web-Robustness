## Define a function for plotting a matrix ----------------
#x : the matrix
#horiz : add a horizontal line 
#pal : name of palette used from RColorBrewer
myHeatmap <- function(x, horiz=NULL, verti =NULL, pal="BuGn", colscale=TRUE, 
                      mary=0, col.l="black", breaks=NULL, labcol= "", 
                      rm.empty=FALSE, cex.x=0.7, cex.y=0.7, 
                      tck.x=-0.015, tck.y=-0.025, padj.x = -0.5, hadj.y=1, ...){
  require(RColorBrewer)
  min <- min(x, na.rm=TRUE)
  max <- max(x, na.rm=TRUE)
  yLabels <- rownames(x)
  xLabels <- colnames(x)
  title <-c()
  # check for additional function arguments
  if( length(list(...)) ){
    Lst <- list(...)
    if( !is.null(Lst$zlim) ){
      min <- Lst$zlim[1]
      max <- Lst$zlim[2]
    }
    if( !is.null(Lst$yLabels) ){
      yLabels <- c(Lst$yLabels)
    }
    if( !is.null(Lst$xLabels) ){
      xLabels <- c(Lst$xLabels)
    }
    if( !is.null(Lst$title) ){
      title <- Lst$title
    }
  }
  # check for null values
  if( is.null(xLabels) ){
    xLabels <- c(1:ncol(x))
  }
  if( is.null(yLabels) ){
    yLabels <- c(1:nrow(x))
  }
  if(colscale) {
    layout(matrix(data=c(1,2), nrow=1, ncol=2), widths=c(4,1), heights=c(1,1))
  }
  # Set color palette
  if (length(pal)==1){
    if (is.null(breaks)){
      ColorRamp<- brewer.pal(9,pal)
      ColorLevels <- seq(min, max, length=10)
    } else {
      ColorRamp<- brewer.pal(length(breaks)-1,pal)
      ColorLevels <- breaks
    }
  } else {
    if (is.null(breaks)){
      ColorRamp<- pal
      ColorLevels <- seq(min, max, length=length(pal)+1)
    } else {
      ColorRamp<- pal
      ColorLevels <- breaks
    }
  }
  
  
  # Reverse Y axis
  reverse <- nrow(x) : 1
  yLabels <- yLabels[reverse]
  x <- x[reverse,]
  
  # Data Map
  par(mar = c(2,5+mary,1.5,1))
  #par(mar = c(0.1,0.1,0.1,0.1))
  image(1:length(xLabels), 1:length(yLabels), t(x), col=ColorRamp, xlab="",
        ylab="", axes=FALSE, zlim=c(min,max), breaks=ColorLevels)
  
  if( !is.null(title) ){
    title(main=title)
  }
  if (!is.null(horiz) ){
    abline(h=horiz, lwd=3, col=col.l)
  }
  if (!is.null(verti) ){
    abline(v=verti, lwd=2)
  }
  if (rm.empty){
    axis(BELOW<-1, at=(1:length(xLabels))[xLabels!=""], 
         labels=xLabels[xLabels!=""], cex.axis=cex.x, tck=tck.x, padj = padj.x)
    axis(LEFT <-2, at=(1:length(yLabels))[yLabels!=""], 
         labels=yLabels[yLabels!=""], las= 1, cex.axis= cex.y, 
         tck=tck.y, hadj = hadj.y)
  } else {
    axis(BELOW<-1,  at=1:length(xLabels), labels=xLabels, 
         cex.axis=cex.x, tck=tck.x, padj = padj.x)
    axis(LEFT <-2,  at=1:length(yLabels), labels=yLabels, 
         las= 1, cex.axis=cex.y, tck=tck.y, hadj = hadj.y)
  }
  
  
  # Color Scale
  if (colscale) {
    par(mar = c(3,2.5,2.5,2))
    nColorLevels<-ColorLevels[-1]-(diff(ColorLevels)/2)
    image(1, nColorLevels,
          matrix(data=nColorLevels, ncol=length(nColorLevels),nrow=1),
          col=ColorRamp,
          xlab="",ylab="",
          xaxt="n")
    axis(2, at = median(nColorLevels), labels = labcol)
  }
  if(colscale) {
    layout(1)
  } else {
    return(ColorLevels)
  }
}


lunique <- function(dat){
  return(length(unique(dat)))
}

panel.cor.m <- function(x, y, digits=2)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r <- cor(x, y,method = "pearson", use = "complete.obs")
  r2 <- abs(cor(x, y,method = "pearson", use = "complete.obs"))
  txt <- format(c(r, 0.123456789), digits=digits)[1]
  cex <- 0.8/strwidth(txt)
  test <- cor.test(x,y,method = "pearson", use = "complete.obs")
  # borrowed from printCoefmat
  Signif <- symnum(test$p.value, corr = FALSE, na = FALSE,
                   cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1),
                   symbols = c("***", "**", "*", ".", " "))
  col1<-ifelse(r<0, "blue", "black")
  text(0.5, 0.5, txt, cex = cex * r2, col=col1)
  text(.8, .8, Signif, cex=cex, col=2)
}

require(RColorBrewer)
plot4C <- function (x, alpha = c(0.1, 0.05, 0.01, 0.001), col = brewer.pal(9, "RdBu")[9:1], txt=TRUE, ...) {
  res <- data.frame(matrix(1, length(x$varnames.Q), length(x$varnames.R)))
  names(res) <- x$varnames.R
  row.names(res) <- x$varnames.Q
  xrand <- x$tabG
  balpha <- sort(c(1,alpha,0))
  for (i in 1:nrow(res)) {
    for (j in 1:ncol(res)) {
      idx.var <- ncol(res) * (i - 1) + j
      catx <- as.numeric(cut(xrand$adj.pvalue[idx.var], breaks= balpha))
      if (xrand$obs[idx.var] < 0) {
        res[i, j] <- catx
      } else {
        res[i, j] <- 10-catx
      }
    }
  }
  x1 <- 1:ncol(res)
  y <- nrow(res):1
  opar <- par(mai = par("mai"), srt = par("srt"))
  on.exit(par(opar))
  table.prepare(x = x1, y = y, row.labels = row.names(res), 
                col.labels = names(res), clabel.row = 1, clabel.col = 1, 
                grid = FALSE, pos = "paint")
  xtot <- x1[col(as.matrix(res))]
  ytot <- y[row(as.matrix(res))]
  xdelta <- (max(x1) - min(x1))/(length(x1) - 1)/2
  ydelta <- (max(y) - min(y))/(length(y) - 1)/2
  z <- unlist(res)
  Signif <- c("***", "**", "*", ".", " ", ".", "*", "**", "***")
  rect(xtot - xdelta, ytot - ydelta, xtot + xdelta, ytot + 
         ydelta, col = col[z], border = "grey90")
  if(txt){
    text(xtot, ytot,  labels= Signif[z])
  }
  rect(min(xtot) - xdelta, min(ytot) - ydelta, max(xtot) + 
         xdelta, max(ytot) + ydelta, col = NULL)
}

table.prepare <- function (x, y, row.labels, col.labels, clabel.row, clabel.col,
                           grid, pos) {
  cexrow <- par("cex") * clabel.row
  cexcol <- par("cex") * clabel.col
  wx <- range(x)
  wy <- range(y)
  maxx <- max(x)
  maxy <- max(y)
  minx <- min(x)
  miny <- min(y)
  dx <- diff(wx)/(length(x))
  dy <- diff(wy)/(length(y))
  if (cexrow > 0) {
    ## ncar <- max(nchar(paste(" ", row.labels, " ", sep = "")))
    ## strx <- par("cin")[1] * ncar * cexrow/2 + 0.1
    strx <- max(strwidth(paste(" ", row.labels, " ", sep = ""), units = "inches", cex=cexrow))
  }
  else strx <- 0.1
  if (cexcol > 0) {
    ##ncar <- max(nchar(paste(" ", col.labels, " ", sep = "")))
    ##stry <- par("cin")[1] * ncar * cexcol/2 + 0.1
    stry <- max(strwidth(paste(" ", col.labels, " ", sep = ""), units = "inches", cex=cexcol))
  }
  else stry <- 0.1
  if (pos == "righttop") {
    par(mai = c(0.1, 0.1, stry, strx))
    xlim <- wx + c(-dx, 2 * dx)
    ylim <- wy + c(-2 * dy, 2 * dy)
    plot.default(0, 0, type = "n", xlab = "", ylab = "", 
                 xaxt = "n", yaxt = "n", xlim = xlim, ylim = ylim, 
                 xaxs = "i", yaxs = "i", frame.plot = FALSE)
    if (cexrow > 0) {
      for (i in 1:length(y)) {
        ynew <- seq(miny, maxy, le = length(y))
        ynew <- ynew[rank(y)]
        text(maxx + 2 * dx, ynew[i], row.labels[i], adj = 0, 
             cex = cexrow, xpd = NA)
        segments(maxx + 2 * dx, ynew[i], maxx + dx, y[i])
      }
    }
    if (cexcol > 0) {
      par(srt = 90)
      for (i in 1:length(x)) {
        xnew <- seq(minx, maxx, le = length(x))
        xnew <- xnew[rank(x)]
        text(xnew[i], maxy + 2 * dy, col.labels[i], adj = 0, 
             cex = cexcol, xpd = NA)
        segments(xnew[i], maxy + 2 * dy, x[i], maxy + 
                   dy)
      }
      par(srt = 0)
    }
    if (grid) {
      col <- "lightgray"
      for (i in 1:length(y)) segments(maxx + dx, y[i], 
                                      minx - dx, y[i], col = col)
      for (i in 1:length(x)) segments(x[i], miny - dy, 
                                      x[i], maxy + dy, col = col)
    }
    rect(minx - dx, miny - dy, maxx + dx, maxy + dy)
    return(invisible())
  }
  if (pos == "phylog") {
    par(mai = c(0.1, 0.1, stry, strx))
    xlim <- wx + c(-dx, 2 * dx)
    ylim <- wy + c(-dy, 2 * dy)
    plot.default(0, 0, type = "n", xlab = "", ylab = "", 
                 xaxt = "n", yaxt = "n", xlim = xlim, ylim = ylim, 
                 xaxs = "i", yaxs = "i", frame.plot = FALSE)
    if (cexrow > 0) {
      for (i in 1:length(y)) {
        ynew <- seq(miny, maxy, le = length(y))
        ynew <- ynew[rank(y)]
        text(maxx + 2 * dx, ynew[i], row.labels[i], adj = 0, 
             cex = cexrow, xpd = NA)
        segments(maxx + 2 * dx, ynew[i], maxx + dx, y[i])
      }
    }
    if (cexcol > 0) {
      par(srt = 90)
      xnew <- x[2:length(x)]
      x <- xnew
      for (i in 1:length(x)) {
        text(xnew[i], maxy + 2 * dy, col.labels[i], adj = 0, 
             cex = cexcol, xpd = NA)
        segments(xnew[i], maxy + 2 * dy, x[i], maxy + 
                   dy)
      }
      par(srt = 0)
    }
    minx <- min(x)
    if (grid) {
      col <- "lightgray"
      for (i in 1:length(y)) segments(maxx + dx, y[i], 
                                      minx - dx, y[i], col = col)
      for (i in 1:length(x)) segments(x[i], miny - dy, 
                                      x[i], maxy + dy, col = col)
    }
    rect(minx - dx, miny - dy, maxx + dx, maxy + dy)
    rect(-dx, miny - dy, minx - dx, maxy + dy)
    return(c(0, minx - dx))
  }
  if (pos == "leftbottom") {
    par(mai = c(stry, strx, 0.05, 0.05))
    xlim <- wx + c(-2 * dx, dx)
    ylim <- wy + c(-2 * dy, dy)
    plot.default(0, 0, type = "n", xlab = "", ylab = "", 
                 xaxt = "n", yaxt = "n", xlim = xlim, ylim = ylim, 
                 xaxs = "i", yaxs = "i", frame.plot = FALSE)
    if (cexrow > 0) {
      for (i in 1:length(y)) {
        ynew <- seq(miny, maxy, le = length(y))
        ynew <- ynew[rank(y)]
        w9 <- strwidth(row.labels[i], cex = cexrow)
        text(minx - w9 - 2 * dx, ynew[i], row.labels[i], 
             adj = 0, cex = cexrow, xpd = NA)
        segments(minx - 2 * dx, ynew[i], minx - dx, y[i])
      }
    }
    if (cexcol > 0) {
      par(srt = -90)
      for (i in 1:length(x)) {
        xnew <- seq(minx, maxx, le = length(x))
        xnew <- xnew[rank(x)]
        text(xnew[i], miny - 2 * dy, col.labels[i], adj = 0, 
             cex = cexcol, xpd = NA)
        segments(xnew[i], miny - 2 * dy, x[i], miny - 
                   dy)
      }
      par(srt = 0)
    }
    if (grid) {
      col <- "lightgray"
      for (i in 1:length(y)) segments(maxx + 2 * dx, y[i], 
                                      minx - dx, y[i], col = col)
      for (i in 1:length(x)) segments(x[i], miny - 2 * 
                                        dy, x[i], maxy + dy, col = col)
    }
    rect(minx - dx, miny - dy, maxx + dx, maxy + dy)
    return(invisible())
  }
  if (pos == "paint") {
    
    dx <- diff(wx)/(length(x) - 1)/2
    dy <- diff(wy)/(length(y) - 1)/2
    
    par(mai = c(0.2, strx, stry, 0.1))
    xlim <- wx + c(-dx, dx)
    ylim <- wy + c(-dy, dy)
    plot.default(0, 0, type = "n", xlab = "", ylab = "", 
                 xaxt = "n", yaxt = "n", xlim = xlim, ylim = ylim, 
                 xaxs = "i", yaxs = "i", frame.plot = TRUE)
    if (cexrow > 0) {
      ynew <- seq(miny, maxy, le = length(y))
      ynew <- ynew[rank(y)]
      ##w9 <- strwidth(row.labels, cex = cexrow)
      ##text(minx - w9 - 3 * dx/4, ynew, row.labels, adj = 0, cex = cexrow, xpd = NA)
      mtext(at =  ynew, side = 2, text = paste(row.labels," ", sep = ""), adj = 1, cex = cexrow, las = 1)
    }
    if (cexcol > 0) {
      xnew <- seq(minx, maxx, le = length(x))
      xnew <- xnew[rank(x)]
      ## par(srt = 90)
      ## text(xnew, maxy + 3 * dy/4, col.labels, adj = 0, cex = cexcol, xpd = NA)
      mtext(at = xnew, side = 3, text = paste(" ", col.labels, sep = ""), adj = 0, cex = cexcol, las = 2)
      par(srt = 0)
    }
    return(invisible())
  }
}

matchmat <- function(dat, col=NULL, row=NULL){
  temp <- c()
  if(!is.null(col)){
    for (i in as.character(col)){
      if (i %in% colnames(dat)){
        temp <- cbind(temp, dat[,i])
      } else {
        temp <- cbind(temp, rep(0, nrow(dat)))
      }
    }
    colnames(temp) <- col
  } else {
    temp <- dat
  }
  res <- c()
  if(!is.null(row)){
    for (j in as.character(row)){
      if (j %in% row.names(temp)){
        res <- rbind(res, temp[j,])
      } else {
        res <- rbind(res, rep(0, ncol(temp)))
      }
    }
    row.names(res) <- row
  } else {
    res <-  temp
  }
  return(res)
}

#colscale : compute the scale of colors uniformly from the range of values
# dat : vector of data
# pal : color palette
colscale <- function(dat, pal, met="ran", na.rm=TRUE){
  if (met=="ran"){
    br <- seq(min(dat, na.rm = na.rm), max(dat, na.rm = na.rm), length.out = length(pal)+1)
  }
  if (met%in%c("ran90","ran95")){
    xq <- quantile(dat, probs = c(0.1, 0.9), na.rm=na.rm)
    br <- c(min(dat, na.rm = na.rm), seq(xq[1], xq[2], length.out = length(pal)-1), max(dat, na.rm = na.rm))
    if (length(unique(br))<length(pal)+1){
      br <- seq(min(dat, na.rm = na.rm), max(dat, na.rm = na.rm), length.out = length(pal)+1)
    }
  } 
  if (met=="abs"){
    maxabs <- max(abs(dat), na.rm=na.rm)
    br <- seq(-maxabs, maxabs, length.out = length(pal)+1)
  }
  if (met=="abs90"){
    xq <- quantile(abs(dat), probs = 0.95, na.rm=na.rm)
    br <- seq(-xq, xq, length.out = length(pal)+1)
    br[1] <- min(dat, na.rm = na.rm)
    br[length(pal)+1] <- max(dat, na.rm = na.rm)
  }
  col <- pal[cut(dat, include.lowest = TRUE, breaks=br)]
  res <- list("br"=br, "col"=col)
  return(res)
}


#add.colscale : add the scale of colors in a plot
# dat : vector of data
# pal : color palette
add.colscale <- function(br, pal,posi = "topleft", rd=1, las=1, lab="", cex=0.8, ratio = 0.2, inset = 0.01){
  add.scatter(plot.colscale(br, pal, rd=rd, las=las, lab=lab, cex=cex, posi=posi), posi=posi, ratio = ratio, inset = inset)
}

#plot.colscale : function to be coupled with add.scatter
plot.colscale <- function(br, pal, rd=1, las=1, posi="left", lab="", 
                          cex=0.8, bg.col="white", n=4){
  opar=par("mar","xaxt","yaxt","plt", "las")
  on.exit(par(opar))
  axp <- ifelse(posi == "bottomright" || posi == "topright", 2, 4)
  adjp <- ifelse(posi == "bottomright" || posi == "topright", 1, 0)
  par(mar=c(0,2,2,2),plt=par("plt"), las=las)
  plot(0, xlim=c(0,1), ylim=c(0.1,length(pal)-0.1), type="n", 
       bty="n", xlab="", ylab="", xaxt="n",yaxt="n")
  for (i in 1:length(pal)){
    rect(0,i-1,1,i,col = pal[i])
  }
  showAx <- round(seq(0,length(pal), length.out=n))
  axis(axp, at=showAx, labels = round(br,rd)[showAx+1], cex.axis=cex)
  mtext(lab, 3, line = 0.5, adj = adjp, cex=cex)
}
#add.scatter : function modified from ade4 package
add.scatter <- function (func, posi = c("bottomleft", "bottomright", "topleft", 
                                        "topright"), ratio = 0.2, inset = 0.01, bg.col="white") {
  if (tolower(posi[1]) == "none") 
    return()
  if (ratio > 0.99) 
    ratio <- 0.99
  if (ratio < 0) 
    ratio <- 0.2
  if (length(inset) == 2) {
    inset.x <- inset[1]
    inset.y <- inset[2]
  }
  else {
    inset.x <- inset[1]
    inset.y <- inset[1]
  }
  inset[inset < 0] <- 0
  plotreg0 <- par("plt")
  plotreg <- plotreg0 + c(inset.x, -inset.x, inset.y, -inset.y)
  on.exit(par(plt = plotreg0))
  posi <- tolower(posi[1])
  if (posi == "bottomleft" || posi == "bottom") {
    x1 <- plotreg[1]
    y1 <- plotreg[3]
  }
  else if (posi == "topleft" || posi == "top") {
    x1 <- plotreg[1]
    y1 <- plotreg[4] - ratio
  }
  else if (posi == "bottomright") {
    x1 <- plotreg[2] - ratio
    y1 <- plotreg[3]
  }
  else if (posi == "topright") {
    x1 <- plotreg[2] - ratio
    y1 <- plotreg[4] - ratio
  }
  else if (posi == "right") {
    x1 <- plotreg[2] - ratio
    y1 <- mean(plotreg[c(3,4)]) - ratio
  }
  else if (posi == "left") {
    x1 <- plotreg[1]
    y1 <- mean(plotreg[c(3,4)]) - ratio
  }
  else stop("Unknown position required")
  x2 <- x1 + ratio
  y2 <- y1 + ratio
  # if (posi == "bottomright" || posi == "topright") {
  #   par(plt = c(x1+(ratio*0.5), x2, y1, y2), new = TRUE)
  # } else {
  #   par(plt = c(x1, x1 + ratio/2, y1, y2), new = TRUE)
  # }
  # plot.new()
  # polygon(c(-0.1, 1.1, 1.1, -0.1), c(-0.1, -0.1, 1.1, 1.1), 
  #         border = NA, col = bg.col)
  if (posi == "bottomright" || posi == "topright" || posi == "right" ) {
    par(plt = c(x1+(ratio*0.8), x2, y1+ (ratio/10), y2- (ratio/4)), new = TRUE)
  } else {
    par(plt = c(x1, x1 + ratio*0.2, y1+ (ratio/10), y2- (ratio/4)), new = TRUE)
    # par(plt = c(x1+(ratio*0.4), x1 + ratio*0.6, y1+ (ratio/10), y2- (ratio/4)), new = TRUE)
  }
  eval(func)
  return(invisible(match.call()))
}

#Sentense case
firstup <- function(x) {
  x <- ifelse(nchar(x)>1, paste0(toupper(substr(x, 1, 1)), tolower(substr(x, 2, nchar(x)))), x)
  return(x)
}


#trophiclevels: compute short weighted trophic levels SWTL 
# from Williams and Martinez 2004, de Santana et al. 2013
# Trophic levels more appropriate for topologies
# Also include values based on the longest path - long weighted trophic levels LWTL
trophiclevels<-function(net){
  mat <- get.adjacency(net, sparse=F)
  
  basal <- rownames(subset(mat, apply(mat, 2, sum)==0) & apply(mat, 1, sum)!= 0)
  paths_prey <- suppressWarnings(shortest.paths(graph = net, v= V(net), to = V(net)[basal], 
                                                mode = "in", weights = NULL, algorithm = "unweighted"))
  
  paths_prey[is.infinite(paths_prey)] <- NA
  shortest_paths <- suppressWarnings(as.matrix(apply(paths_prey, 1, min, na.rm=TRUE)))
  longest_paths <- suppressWarnings(as.matrix(apply(paths_prey, 1, max, na.rm=TRUE)))
  #for species with no prey apart of them
  shortest_paths[is.infinite(shortest_paths)] <- NA
  longest_paths[is.infinite(longest_paths)] <- NA
  
  # Shortest TL
  sTL <- 1 + shortest_paths 
  # Longest TL
  lTL <- 1 + longest_paths
  
  #Compute average prey trophic level
  W <- t(mat)
  rs <- rowSums(W)
  W <- W/matrix(rs, ncol = ncol(W), nrow = nrow(W))
  W[0 == rs, ] <- 0
  I <- diag(ncol(W))
  tl0<-rowSums(I - W)
  result <- tryCatch(solve(I - W), error = function(e) e)
  if ("error" %in% class(result)) {
    avtl <- rep(NA, ncol(pm))
    names(avtl) <- colnames(pm)
  }
  else {
    avtl <- rowSums(result)
  }
  
  # Short-weighted TL is the average of shortest TL and prey-averaged TL
  SWTL <- (sTL + avtl)/2
  
  # Long-weighted TL is the average of longest TL and prey-averaged TL
  LWTL <- (lTL + avtl)/2
  
  #check that only basal species have TL of 1
  SWTL[!rownames(mat)%in%basal & SWTL == 1] <- NA
  LWTL[!rownames(mat)%in%basal & LWTL == 1] <- NA
  
  res<-data.frame("swTL"= SWTL, "lwTL"=LWTL)
  rownames(res)<-rownames(mat)
  return(res)
}


#PlotFW
plotfw <- function(net, met=c("tl", "deg"), col=NULL, lab=NULL, yax=TRUE, 
                   labcex=0.01, size=NULL, coo=NULL, maxsize=10, 
                   intx=10, retcoo=FALSE, rescale=TRUE, nylevel=6, 
                   ewidth=0.3, esize=0.3,labcol="black",...){
  n <- vcount(net)
  if (is.null(col)){
    V(net)$color <- rainbow(vcount(net))
  } else{
    V(net)$color <- col
  }
  if (!is.null(lab)){
    V(net)$name <- lab
  }
  if (is.null(size)){
    V(net)$size <- (degree(net)/max(degree(net))*maxsize)+1
  } else {
    V(net)$size <- size
  }
  if (is.null(coo)){
    if ("tl"%in%met){
      tl <- trophiclevels(net)
      dgpred <- tl$swTL
      uny <- "Trophic Level"
    } else {
      dgin <- degree(net, mode="in")
      dgout <- degree(net, mode="out")
      dgpred <- (dgin-dgout)/degree(net)
      dgpred <- ifelse(is.nan(dgpred), 0, dgpred) #nodes disconnected, degree=0
      uny <- "(Din-Dout)/(Din+Dout)"
    }
    bks <- quantile(dgpred, probs = seq(0,1, length.out = nylevel))
    bks <- sort(unique(bks))
    ynod <- cut(dgpred, breaks = bks, include.lowest = TRUE, 
                labels = 1:(length(bks)-1))
    xnod <- rep(0,n)
    for (i in 1:(length(bks)-1)){
      l <- sum(ynod==i)
      xnod[ynod==i] <- seq(-l,l,length.out = l)[order(dgpred[ynod==i])]
    }
    coo <- cbind(xnod,dgpred)
  } else {
    xnod <- coo[,1]
    dgpred <- coo[,2]
    uny <- ""
  }
  if (max(dgpred)-min(dgpred)>10){
    xax <- seq(min(dgpred), max(dgpred))
    labax <- xax[xax%%intx==0]
    laby <- seq(-1, 1, length.out = length(xax))[xax%%intx==0]
  } else {
    if (rescale){
      xax <- seq(min(dgpred), max(dgpred), length.out = intx)
      labax <- round(xax, 1)
      laby <- seq(-1, 1, length.out = intx)
    } else {
      xax <- seq(min(dgpred), max(dgpred), length.out = intx)
      labax <- round(xax, 1)
      laby <- seq(min(dgpred), max(dgpred), length.out = intx)
    }
  }
  plot(net, layout= coo, vertex.label.color=labcol, 
       edge.width=ewidth, edge.arrow.size=esize,
       vertex.label.cex=labcex, rescale=rescale,...)
  if (yax){
    axis(2, at = laby, labels = labax)
    mtext(uny, side = 2, line=2)
  }
  if (retcoo){
    res <- data.frame(coo, V(net)$size, V(net)$color)
    names(res) <- c("x", "y", "size", "color")
    row.names(res) <- V(net)$name
    return(res)
  }
}

plotfwtl <- function(net, col=NULL, lab=NULL, yax=TRUE, 
                   labcex=0.01, size=NULL, rescale=TRUE, maxsize=10, 
                   intx=10, retcoo=FALSE, nylevel=6, 
                   ewidth=0.3, esize=0.3,labcol="black",...){
  n <- vcount(net)
  if (is.null(col)){
    V(net)$color <- rainbow(vcount(net))
  } else{
    V(net)$color <- col
  }
  if (!is.null(lab)){
    V(net)$label <- lab
  }
  if (is.null(size)){
    V(net)$size <- (degree(net)/max(degree(net))*maxsize)+1
  } else {
    V(net)$size <- size
  }
  
  tl <- trophiclevels(net)
  dgpred <- tl$swTL
  uny <- "Trophic Level"
  
  bks <- c(0.9, seq(1.9, max(dgpred), length.out = nylevel))
  ynod <- cut(dgpred, breaks = bks, include.lowest = TRUE, 
                labels = 1:(length(bks)-1))
  maxx <- max(table(ynod))
  xnod <- rep(0,n)
  for (i in 1:nylevel){
    l <- sum(ynod==i)
    #tranformation of l
    ltr <- (l/maxx)**(1/2)*maxx
    if (l>1) {
      xnod[ynod==i] <- seq(-ltr,ltr,length.out = l)#[order(dgpred[ynod==i])]
    } else {
      xnod[ynod==i] <- 0
    }
  }
  coo <- cbind(xnod,dgpred)

  xax <- c(1, seq(2, max(dgpred), length.out = intx))
  labax <- round(xax, 1)
  #rescale xax between -1 and 1
  laby <- (xax-min(xax))/(max(xax)-min(xax))*2 - 1
  
  plot(net, layout= coo, vertex.label.color=labcol, 
       edge.width=ewidth, edge.arrow.size=esize,
       vertex.label.cex=labcex, rescale=rescale,...)
  if (yax){
    axis(2, at = laby, labels = labax)
    mtext(uny, side = 2, line=2)
  }
  if (retcoo){
    res <- data.frame(coo, V(net)$size, V(net)$color)
    names(res) <- c("x", "y", "size", "color")
    row.names(res) <- V(net)$name
    return(res)
  }
}


#Omnivory Petchey code
Fraction.omnivores <- function(net, TLs=trophiclevels(net)[,1]) {
  web <- get.adjacency(net, sparse=F)
  non.int.TL <- web[,TLs %% 1 != 0]
  if(is.matrix(non.int.TL))
    frac.omniv <- sum(apply(non.int.TL, 2, sum) > 1)  / length(web[,1])
  if(is.vector(non.int.TL)) 
    frac.omniv <- (sum(non.int.TL) > 1)  / length(web[,1])
  frac.omniv
}

Level.omnivory <- function(net, TLs=trophiclevels(net)[,1]) {
  web <- get.adjacency(net, sparse=F)
  if( sum(is.na(TLs)) == length(TLs) )
    rr <- NA
  
  if( sum(is.na(TLs)) != length(TLs) ) {
    web.TLs <- matrix(rep(TLs, length(web[,1])), length(web[,1]), length(web[,1]))
    lo.pc <- numeric(length=length(web[,1]))
    for(i in 1:length(web[,1])) {
      tt <- web.TLs[web[,i]==1,i]
      if(length(tt)==0 | sum(!is.na(tt))==0 )
        lo.pc[i] = NA
      if(length(tt)>0 & sum(!is.na(tt))!=0)
        lo.pc[i] <- sd(tt)
    }
    rr <- mean(lo.pc, na.rm=T)
  }
  rr
}

#Persistence plot
# tab: matrice of abundance with species in row, time in column
persistplot <- function(tab, model="both", col=c("limegreen", "dodgerblue"),  main="Persistence plot",...){
  if (!tolower(model)%in%c("both", "loess", "poly")){
    stop("model should be either 'both', 'loess', 'poly'")
  }
  #calculate total abundance and number of years
  totA <- apply(tab, 1, sum, na.rm=TRUE)
  nbY <- apply(tab>0, 1, sum, na.rm=TRUE)
  #remove species not seen
  nbY <- nbY[totA>0]
  totA <- totA[totA>0]
  #plot
  plot(nbY, log(totA), main=main, xlab="number of years",
       ylab="log(total biomass)", ...)
  xseq <- seq(min(nbY), max(nbY), length.out = 100)
  #text(nbY, log(totA), labels = substr(names(nbY),1,6), xpd=NA)
  highY <- c()
  highT <- c()
  #Polynomial
  if (tolower(model)%in%c("both", "poly")){
    nbY2 <- nbY**2
    nbY3 <- nbY**3
    poly <-lm(log(totA) ~ nbY + nbY2 + nbY3)
    ypred <- predict(poly,list(nbY=xseq, nbY2=xseq**2, nbY3=xseq**3))
    lines(xseq, ypred, col = col[1], lwd = 2)
    infl <- which.min(diff(ypred))+1
    points(xseq[infl], ypred[infl], col=col[1], pch=3)
    segments(xseq[infl], par()$usr[3], xseq[infl], ypred[infl], lty=2, col=col[1])
    segments(par()$usr[1], ypred[infl], xseq[infl], ypred[infl], lty=2, col=col[1])
    highY <- unique(c(highY, names(nbY)[nbY>xseq[infl]]))
    highT <- unique(c(highT, names(totA)[log(totA)>ypred[infl]]))
  }
  #LOESS
  if (tolower(model)%in%c("both", "loess")){
    ncol <- ifelse(model=="both", 2, 1)
    lo <- loess(log(totA) ~ nbY)
    ypred = predict(lo,xseq)
    lines(xseq, ypred, col = col[ncol], lwd = 2)
    infl <- which.min(diff(ypred))+1
    points(xseq[infl], ypred[infl], col=col[ncol], pch=3)
    segments(xseq[infl], par()$usr[3], xseq[infl], ypred[infl], lty=2, col=col[ncol])
    segments(par()$usr[1], ypred[infl], xseq[infl], ypred[infl], lty=2, col=col[ncol])
    highY <- unique(c(highY, names(nbY)[nbY>xseq[infl]]))
    highT <- unique(c(highT, names(totA)[log(totA)>ypred[infl]]))
  }
  if (model=="both"){
    legend("bottomright", lwd = 2, legend=c("Polynomial", "LOESS"), 
           col=col, bty="n")
  }
  #show <- log(totA)>ypred[infl]&nbY<xseq[infl]
  #text(nbY[show], log(totA)[show], labels = short(names(totA))[show], pos = 3)
  res <- list("HighYear"=highY, "HighBiomass"=highT[!highT%in%highY])
  return(res)
}


PCheatmap <- function(dat, col=c("springgreen4","chartreuse3","yellow","darkgoldenrod1", "red"), 
                         npc=2, n=1, shortname=short(colnames(dat))){
  layout(matrix(data=c(2,1,4,3), nrow=2), widths=c(4,1), heights=c(1,4))
  require(ade4)
  dat <- t(dat)
  #remove years with NA values
  scdat <- scale(dat)
  scdat[dat==0] <- NA
  keepYr <- complete.cases(dat)
  dat2 <- dat[keepYr,]
  mvar <- dudi.pca(dat2, nf=npc, scannf = FALSE, center = TRUE, scale = TRUE)
  odPC <- order(mvar$co[,n])
  bk <- quantile(scdat, probs = seq(0,1,length.out = length(col)+1), na.rm=TRUE)
  bk <- bk + c(-0.1, rep(0,length(col)-1), 0.1) #to be sure to include all values
  par(mar=c(4,6,0,0))
  image(scdat[,odPC], axes=FALSE, col = col, breaks=bk)
  box()
  axis(1, at = (seq(0,nrow(dat), length.out = nrow(dat)))/nrow(dat), 
       labels = row.names(dat))
  axis(2, at = (seq(0,ncol(dat), length.out = ncol(dat)))/ncol(dat), 
       labels = shortname[odPC], las=1, xpd=NA)
  tsPC1 <- rep(NA, nrow(dat))
  tsPC1[keepYr] <- mvar$li[,n]
  par(mar=c(0,6,1,0), xaxs="i", yaxs="r")
  yr <- as.numeric(row.names(dat))
  plot(yr, tsPC1, type="b", pch=16,
       ylab="loadings", xaxt="n", xlim=range(yr)+c(-0.5, 0.5))
  abline(h=0, lty=3)
  par(mar=c(4,0,0,1), yaxs="i", xaxs="r")
  plot(mvar$co[odPC,n],1:ncol(dat), type="p",
       xlab="loadings", yaxt="n", pch=16, 
       ylim=c(0.5, ncol(dat)+0.5))
  abline(v=0, lty=3)
  labpc <- paste0("PC",n, ": ", round(mvar$eig[n]/sum(mvar$eig)*100), "%")
  mtext(labpc, side = 3, line = 0, xpd=NA)
  par(mar=c(2,5,2,2), yaxs="r", xaxs="r")
  plot(0, xlim=c(0,1), ylim=c(0.1,length(col)-0.1), type="n", 
       xaxt="n", yaxt="n", bty="n", xlab="", ylab="")
  for (i in 1:length(col)){
    rect(0,i-1,1,i,col = col[i])
  }
  axis(2, at=c(0,length(col)/2, length(col)), labels = c("Low", "Mean", "High"), las=1)
} 

BINheatmap <- function(dat, col=c("red","springgreen4"), 
                      npc=2, n=1, shortname=short(colnames(dat))){
  layout(matrix(data=c(2,1,4,3), nrow=2), widths=c(4,1), heights=c(1,4))
  require(ade4)
  bin <- t(dat>0)
  bin[bin] <- 1
  bin[bin==0] <- -1
  keepYr <- complete.cases(bin)
  keepSp <- apply(bin,2,sum)<nrow(bin)
  dat2 <- bin[keepYr,keepSp]
  mvar <- dudi.pca(dat2, nf=npc, scannf = FALSE, center = TRUE, scale = TRUE)
  odPC <- order(mvar$co[,n])
  odPC2 <- match(colnames(dat2)[odPC], colnames(bin))
  par(mar=c(4,6,0,0))
  image(bin[,c(which(!keepSp), odPC2)], axes=FALSE, col = col, breaks=c(-1.5,0,1.5))
  box()
  axis(1, at = (seq(0,nrow(bin), length.out = nrow(bin)))/nrow(bin), 
       labels = row.names(bin))
  axis(2, at = (seq(0,ncol(bin), length.out = ncol(bin)))/ncol(bin), 
       labels = shortname[c(which(!keepSp), odPC2)], las=1, xpd=NA)
  tsPC1 <- rep(NA, nrow(bin))
  tsPC1[keepYr] <- mvar$li[,n]
  par(mar=c(0,6,1,0), xaxs="i", yaxs="r")
  yr <- as.numeric(row.names(bin))
  plot(yr, tsPC1, type="b", pch=16,
       ylab="loadings", xaxt="n", xlim=range(yr)+c(-0.5, 0.5))
  abline(h=0, lty=3)
  par(mar=c(4,0,0,1), yaxs="i", xaxs="r")
  plot(c(rep(0, sum(!keepSp)), mvar$co[odPC,n]), 1:ncol(bin), type="p",
       xlab="loadings", yaxt="n", pch=16, 
       ylim=c(0.5, ncol(bin)+0.5))
  abline(v=0, lty=3)
  labpc <- paste0("PC",n, ": ", round(mvar$eig[n]/sum(mvar$eig)*100), "%")
  mtext(labpc, side = 3, line = 0, xpd=NA)
  par(mar=c(2,5,2,2), yaxs="r", xaxs="r")
  plot(0, xlim=c(0,1), ylim=c(0.1,length(col)-0.1), type="n", 
       xaxt="n", yaxt="n", bty="n", xlab="", ylab="")
  for (i in 1:length(col)){
    rect(0,i-1,1,i,col = col[i])
  }
  axis(2, at=c(1.5,0.5), labels = c("Presence", "Absence"), las=1)
} 


short <- function(dat, nchar = 3){
  #find the separator
  sepL <- c(length(grep("\\.", dat)), length(grep("_", dat)), length(grep(" ", dat)))
  sep <- c("\\.", "_", " ")[which.max(sepL)]
  #split
  return(unlist(lapply(strsplit(dat, sep), shorten)))
}

shorten <- function(dat, nchar = 3){
  if (length(dat)>1){
    new <- paste(substr(dat[1], 1, nchar), substr(dat[2], 1, nchar), sep = "_")
  } else {
    new <- substr(dat[1], 1, nchar*2)
  }
  return(new)
}

#Plotlogs is a visualization of time series
plotts<-function(dat, col=NULL, leg=NULL, xax=TRUE, shortname=short(colnames(dat)), ...){
  if (is.null(col)){
    col<-rainbow(ncol(dat))
  }
  limy<-range(dat, na.rm = TRUE)
  plot(dat[,1], ylim=limy, col=col[1], type="l", xaxt="n", ...)
  points(dat[,1],col=col[1], pch=16)
  if(xax){
    axis(1,1:nrow(dat), row.names(dat))
  }
  for (i in 2:ncol(dat)){
    lines(dat[,i], col=col[i])
    points(dat[,i],col=col[i], pch=16)
  }
  if (!is.null(leg)){
    nco <- ifelse(ncol(dat)>10, 2, 1)
    ordL <- order(apply(dat, 2, sum, na.rm=TRUE), decreasing = TRUE)
    legend(leg, legend = shortname[ordL], col=col[ordL], lty=1, ncol = nco, cex=0.8)
  }
}

cleanspl <- function(name){
  temp <- paste0(toupper(substr(name, 1, 1)), tolower(substr(name, 2, nchar(name))))
  temp <- gsub("\\.$", "", temp) #remove last .
  temp <- gsub(" $", "", temp) #remove empty space
  temp <- gsub("^ ", "", temp)
  temp <- gsub(" .$", "", temp) #remove last single character
  temp <- gsub(" spp$", "", temp)
  temp <- gsub(" sp$", "", temp)
  return(temp)
}

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
  
  return(res)
}


#Compute food web topology indicators
fwind <- function(net, ab=NULL, los=NULL, ef=NULL, fg=NULL, retnet = FALSE, loop=TRUE){
  if(sum(degree(net)==0)>0)
    {warning(paste("Presence of", sum(degree(net)==0), "singleton !"))}
  if(components(net)$no>1)
    {warning(paste("Presence of disconnected sub-graphs!"))}
  res <- c()
  
  #Species richness : number of nodes
  res$S <- vcount(net)
  
  #Number of links
  res$L <- ecount(net)
  
  #Number of links : number of nodes
  res$LD <- ecount(net)/vcount(net)
  
  #Connectance - Density
  #proportion of present edges from all possible edges in the network
  res$C <- edge_density(net, loops=loop)
  
  # Modularity : 
  # Community structureal measures
  res$Clust<- transitivity(net,vids=NULL)
  
  # Compartmentalisation estimation based on simulated annealing 
  # (similar procedure to Guimera et al. 2010)
  spingc <- spinglass.community(net) 		
  res$Mod <- spingc$modularity 	
  
  #Generality: the mean number of prey per consumer
  pred <- degree(net, mode="in")>0
  res$G <- sum(degree(net, mode="in")[pred])/sum(pred)
  res$Gsd <- sd(degree(net, mode="in"))/res$LD
  
  #Vulnerability: the mean number of consumers per prey
  prey <- degree(net, mode="out")>0
  res$V <- sum(degree(net, mode="out")[prey])/sum(prey)
  res$Vsd <- sd(degree(net, mode="out"))/res$LD 
  
  # Proportion of cannibals
  web <- get.adjacency(net, sparse=F)
  res$Can <- sum(diag(web))/dim(web)[1]				# Fraction of cannibals
  
  #Longest path - Diameter
  res$LongPath <- diameter(net)
  # Mean shortest path length
  # removing the diagonal 
  sp <- shortest.paths(net)   # isSymmetric(sp) # TRUE
  res$ShortPath <- mean(sp[upper.tri(sp)]) 
  # Mean path
  res$MeanPath<- average.path.length(net)
  
  #Degree Centrality
  res$DegC <- centr_degree(net, mode="all", normalized=T)$centralization
  
  #Closeness centrality 
  # based on distance to others in the graph
  res$CloC <- centr_clo(net, mode="all", normalized=T)$centralization
  
  # Betweenness centrality 
  # based on a broker position connecting others
  res$BetC <- centr_betw(net, directed=T, normalized=T)$centralization
  
  tlnodes <- trophiclevels(net)
  res$TL <-mean(tlnodes$swTL)
  
  # Level of omnivory
  res$PercOmni <- Fraction.omnivores(net, TLs=tlnodes$swTL) 
  res$Omni <- Level.omnivory(net, TLs=tlnodes$swTL)
  
  #Motifs
  triad.count <- triad.census(net)
  names(triad.count) <- c("empty", "single", "mutual", "s5", "s4", "s1", "d4",
                          "d3", "s2", "s3","d8", "d2", "d1", "d5", "d7", "d6")
  
  res$MotS1 <- as.numeric(triad.count["s1"])/ecount(net)
  res$MotS2 <- as.numeric(triad.count["s2"])/ecount(net)
  res$MotS3 <- as.numeric(triad.count["s3"])/ecount(net)
  res$MotS4 <- as.numeric(triad.count["s4"])/ecount(net)
  res$MotS5 <- as.numeric(triad.count["s5"])/ecount(net)
  
  #Node weighted indicators
  if (!is.null(ab)){
    if(vcount(net) != length(ab)){stop("wrong dimensions")}
    #Weighted average number of link
    res$wLD <- sum(degree(net)*ab)/(2*(sum(ab)))
    
    #Weighted connectance
    if (loop){
      res$wC <- sum(degree(net)*ab)/(2*sum(ab)*vcount(net))
    } else {
      res$wC <- sum(degree(net)*ab)/(2*sum(ab)*(vcount(net)-1))
    }
    
    #Weighted generality
    res$wG <- sum((degree(net, mode="in")*ab)[pred])/(sum(ab[pred]))
    
    #Weighted vulnerability
    res$wV <- sum((degree(net, mode="out")*ab)[prey])/(sum(ab[prey]))
    
    res$wTL <- sum(tlnodes$swTL*ab)/sum(ab)
  }
  
  if (!(is.null(los) & is.null(ef))){
    adj<-as.matrix(get.adjacency(net))
    #Try fluxing - sometimes results in error
    fluxI <- try(expr = fluxing(adj, ab, los, ef, ef.level="prey"), silent = TRUE)
    #Error in fluxing(adjI, abI, losI, efI) : 
    # model chosen is unable to determine fluxes accoringly to data
    if(length(fluxI)>1){
      fluxI <- fluxI*86.4 #conversion from J/sec to kJ/day 
      
      if (!is.null(fg)){
        # if DOM is a node, then include it in detrivory 
        idDet<-which(colnames(fluxI)%in%c("Detritus", "Dissolved organic matter")) 
        res$detritivory<-SumrowSums(fluxI[idDet,])
        # if DOM is a node, then include it in detrivory 
        idBact<-which(colnames(fluxI)%in%c("Bacteria")) 
        res$bacterivory<-SumrowSums(fluxI[idBact,])
        idPhy<-which(fg=="Phytoplankton")
        res$nphyto<-length(idPhy)
        res$phytoplanktivory<-SumrowSums(fluxI[idPhy,])
        idZoo<-which(fg=="Zooplankton")
        res$nZoo<-length(idZoo)
        res$zooplanktivory<-SumrowSums(fluxI[idZoo,])
        idBen<-which(fg=="Benthos")
        res$nBenthos<-length(idBen)
        res$benthivory<-SumrowSums(fluxI[idBen,])
        idFish<-which(fg=="Fish")
        res$nFish<-length(idFish)
        res$piscivory<-SumrowSums(fluxI[idFish,])
        res$totFluxes<-sum(fluxI)
      }
        
      #tnet indicators
      w_net<-suppressWarnings(as.tnet(fluxI, type="weighted one-mode tnet")) 
      #betweenness
      #weighted betweenness centrality
      bet_w<-betweenness_w(w_net, directed=NULL, alpha=1)
      #rownames(bet_w)<-names(V(netI))
      #mean weighted betweenness centrality
      res$qBetC<-mean(bet_w[,2])
      
      #weighted closeness centrality
      #Returns a data.frame with three columns: 
      #the first column contains the nodes' ids, 
      #the second column contains the closeness scores, 
      #and the third column contains the normalised closeness scores (i.e., divided by N-1).
      clo_w<-closeness_w(w_net, directed=NULL, gconly=FALSE, precomp.dist=NULL, alpha=1)
      #rownames(clo_w)<-names(V(netI))
      #mean weighted closeness centrality
      res$qCloC<-mean(clo_w[,2])
      
      res<-c(res, fluxind(fluxI))
      
    } else {
      if (!is.null(fg)){
        res$ndet<-NA
        res$detritivory<-NA
        res$nphyto<-NA
        res$phytoplanktivory<-NA
        res$nZoo <- NA
        res$zooplanktivory <- NA
        res$nBenthos<- NA
        res$benthivory<- NA
        res$nFish<-NA
        res$piscivory<-NA
        res$totFluxes<-NA
      }
      res$qBetC<-NA
      res$qCloC<-NA
      res$qLD.uw<-NA
      res$qC.uw<-NA
      res$qLD.w<-NA
      res$qC.w<-NA
      res$qL<-NA
      res$qG.uw<-NA
      res$qG.w<-NA
      res$qV.uw<-NA
      res$qV.w<-NA
      res$qGsd.uw<-NA
      res$qGsd.w<-NA
      res$qVsd.uw<-NA
      res$qVsd.w<-NA
    }
  }
  
  if (retnet) {
    res2 <- list("res"=res, "net"=web, "netF"=fluxI)
    return(res2)
  } else {
    return(res)
  }

}

fwindts <- function(net, biom, los=NULL, ef=NULL, fg=NULL, retnet=FALSE, loop=TRUE){
  if(any(V(net)$name!=row.names(biom))){
    stop("Names of net and biom do not match")
  }
  res <- c()
  if (retnet) {
    listnet <- list()
    listnetF <- list()
  }
  for (i in 1:ncol(biom)){
    abI <- biom[,i]
    spI <- row.names(biom)[abI>0]
    if(!is.null(los)){
      losI <- los[abI>0]
    } else {
      losI <- NULL
    }
    if(!is.null(ef)){
      efI <- ef[abI>0]
    } else {
      efI <- NULL
    }
    if(!is.null(fg)){
      fgI <- fg[abI>0]
    } else {
      fgI <- NULL
    }
    abI <- abI[abI>0]
    netI <- delete_vertices(net, V(net)$name[!V(net)$name %in% spI])
    #Compute indicators
    tmp <- fwind(netI, ab=abI, los=losI, ef=efI, fg=fgI, retnet=retnet, loop = loop)
    if (retnet){
      res <- rbind(res, tmp$res)
      listnet[[i]] <- tmp$net
      listnetF[[i]] <- tmp$netF
    } else {
      res <- rbind(res, tmp)
    }
  }
  res <- apply(res, 2, as.numeric)
  res <- as.data.frame(res)
  row.names(res) <-  colnames(biom)
  if (retnet) {
    res2 <- list("res"=res, "net"=listnet, "netF"= listnetF)
    return(res2)
  } else {
    return(res)
  }
}

#Moving window average (if fun=mean)
movwin <- function(x, w, fun = mean, ...){
  res <- c()
  for (i in 1:(length(x)-w+1)){
    res <- c(res, fun(x[i:(i+w-1)], ...))
  }
  return(res)
}

movmat <- function(mat, win, fun, ...){
  mat <- as.matrix(mat)
  res <- c()
  for (i in 1:nrow(mat)){
    res <- rbind(res, movwin(mat[i,], win, fun, ...))
  }
  #newnames <- paste0(substr(movwin(colnames(mat), win, fun=min), 3,4), 
  #                   "-", substr(movwin(colnames(mat), win, fun=max), 3, 4))
  newnames <- movwin(as.numeric(colnames(mat)), win, fun=mean)
  colnames(res) <- newnames
  row.names(res) <- row.names(mat)
  return(res)
}

firstword <- function(x){
  sapply(strsplit(x, split = " "), function(y)y[1])
}

SumrowSums<- function(dat){
  if (length(dim(dat))>1){
    out <-sum(rowSums(dat))
  } else {
    out<-sum(dat)
  }
  return(out)
}

add.PCaxis <- function(eig, x=1, y=2, cex=0.7){
  #Name
  gy <- strheight("PCN", cex = cex)/1.5
  gx <- strwidth("00 %", cex = cex)/1.5
  text(par("usr")[2] - gx, gy, col = "grey40", cex = cex, 
       labels = paste0("PC", x), xpd=NA)
  text(-gy, par("usr")[3] + gx, col = "grey40", cex = cex, 
       labels = paste0("PC", y), srt = 90, xpd=NA)
  #Variance
  var <- paste(round(eig[c(x,y)]/sum(eig),2)*100, "%")
  gx <- strwidth("00 %", cex = cex)/1.5
  gy <- strheight("00 %", cex = cex)/1.5
  text(par("usr")[2] - gx, -gy, col = "grey40", cex = cex, 
       labels = var[1], xpd=NA)
  text(+gy, par("usr")[3] + gx, col = "grey40", cex = cex, 
       labels = var[2], srt = 90, xpd=NA)
}

#Species accumulation curves
speacc <- function(species, id, nrep=100, plot=TRUE, ...){
  maxsp <- lunique(species)
  maxh <- lunique(id)
  idh <- sort(unique(id))
  spaccF <- matrix(NA, nrow = round(maxh/2), ncol=nrep)
  for (i in 1:round(maxh/2)){
    for (k in 1:nrep){
      ha <- idh[sample(seq_along(idh), i, replace = FALSE)]
      nsp <- lunique(species[id%in%ha])
      spaccF[i,k] <- nsp
    }
  }
  if (plot){
    limY <- c(min(apply(spaccF,1,mean)-apply(spaccF,1,sd)), max(apply(spaccF,1,mean)+apply(spaccF,1,sd)))
    plot(1:nrow(spaccF), apply(spaccF,1,mean), type="l", lwd=2, ylim=limY,
         ylab="number of species",xlab="number of stations", ...)
    lines(1:nrow(spaccF), apply(spaccF,1,mean)-apply(spaccF,1,sd), col="grey")
    lines(1:nrow(spaccF), apply(spaccF,1,mean)+apply(spaccF,1,sd), col="grey")
    legend("bottomright", legend = c("mean", "+/- sd"), lwd=c(2,1),
           col=c("black", "grey"), bty="n")
  } else {
    return(spaccF)
  }
}


#From ade4
scatterutil.eti.circ2 <- function (x, y, label, clabel, origin = c(0, 0), boxes = TRUE, col="black", ...) 
{
  if (is.null(label)) 
    return(invisible())
  if (any(is.na(label))) 
    return(invisible())
  if (any(label == "")) 
    return(invisible())
  xref <- x - origin[1]
  yref <- y - origin[2]
  for (i in 1:(length(x))) {
    cha <- as.character(label[i])
    cha <- paste(" ", cha, " ", sep = "")
    cex0 <- par("cex") * clabel
    xh <- strwidth(cha, cex = cex0)
    yh <- strheight(cha, cex = cex0) * 5/6
    if ((xref[i] > yref[i]) & (xref[i] > -yref[i])) {
      x1 <- x[i] + xh/2
      y1 <- y[i]
    }
    else if ((xref[i] > yref[i]) & (xref[i] <= (-yref[i]))) {
      x1 <- x[i]
      y1 <- y[i] - yh
    }
    else if ((xref[i] <= yref[i]) & (xref[i] <= (-yref[i]))) {
      x1 <- x[i] - xh/2
      y1 <- y[i]
    }
    else if ((xref[i] <= yref[i]) & (xref[i] > (-yref[i]))) {
      x1 <- x[i]
      y1 <- y[i] + yh
    }
    if (boxes) {
      rect(x1 - xh/2, y1 - yh, x1 + xh/2, y1 + yh, col = "white", 
           border = col[i], ...)
    }
    text(x1, y1, cha, cex = cex0, col=col[i], ...)
  }
}

#DFA - from comita package
#' Run a Dynamic Factor Analysis with diagonal and unequal matrix model
#'
#' @param dat matrix with 
#' @param npc number of processes
#' @return Result of multivariate analysis with 
#' \itemize{
#' \item \code{ts} scores of the time on the principal component
#' \item \code{co} scores of the variables on the principal component
#' \item \code{eig} eigen values
#' \item \code{dat} original dataset
#' \item \code{npc} number of selected principal component
#' \item \code{mita} method used for integrated trend analysis 
#' }
#' @note R code inspired from https://nwfsc-timeseries.github.io/atsa-labs/sec-dfa.html
#' Holmes, E. E., M. D. Scheuerell, and E. J. Ward. Applied time series analysis for fisheries and environmental data. NOAA Fisheries, Northwest Fisheries Science Center, 2725 Montlake Blvd E., Seattle, WA 98112.
#' @references Zuur AF, Fryer RJ, Jolliffe IT, Dekker R, Beukema JJ (2003) Estimating common trends in multivariate time series using dynamic factor analysis. Environmetrics 14:665–685
#' 
dfa <- function(dat, npc){
  #parameters 
  mm <- npc ## number of processes
  BB <- "identity"  # diag(mm) ## 'BB' is identity: 1's along the diagonal & 0's elsewhere
  uu <- "zero"  # matrix(0,mm,1)  ## 'uu' is a column vector of 0's
  ## 'CC' and 'cc' are for covariates
  CC <- "zero"  # matrix(0,mm,1) 
  cc <- "zero"  # matrix(0,1,wk_last)
  QQ <- "identity"  # diag(mm) ## 'QQ' is identity
  aa <- "zero" ## 'aa' is the offset/scaling
  ## 'DD' and 'd' are for covariates
  DD <- "zero"  # matrix(0,mm,1)
  dd <- "zero"  # matrix(0,1,wk_last)
  RR <- "diagonal and unequal"## 'RR' is var-cov matrix for obs errors
  ## 'ZZ' is loadings matrix
  #tab <- tab[,6:10]
  Z_vals <- c()
  for (i in 1:ncol(dat)){
    for (j in 1:npc){
      Z_vals <- c(Z_vals, ifelse(j<=i, paste0("z",i, j), 0))
    }
  }
  ZZ <- matrix(Z_vals, nrow=ncol(dat), ncol=npc, byrow=TRUE)
  ## list with specifications for model vectors/matrices
  mod_list <- list(B=BB, U=uu, C=CC, c=cc, Q=QQ, Z=ZZ, A=aa, D=DD, d=dd, R=RR)
  ## list with model inits
  init_list <- list(x0=matrix(rep(0,mm),mm,1))
  ## list with model control parameters
  con_list <- list(maxit=3000, allow.degen=TRUE)
  
  ## fit MARSS
  ret <- as.matrix(t(dat))
  mvar <- MARSS::MARSS(y=ret, model=mod_list, inits=init_list, control=con_list)
  
  ## get the estimated ZZ
  Z_est <- coef(mvar, type="matrix")$Z
  ## get the inverse of the rotation matrix
  H_inv <- stats::varimax(Z_est)$rotmat
  ## rotate factor loadings
  Z_rot = Z_est %*% H_inv   
  ## rotate processes
  proc_rot = solve(H_inv) %*% mvar$states
  
  ts=t(proc_rot)
  co=Z_rot
  #compute variance explained
  res <- list("ts"=t(proc_rot), "co"=Z_rot, "eig"= NULL, 
              "dat"=dat, "npc"=npc, "mita"="DFA")
  return(res)
}

#https://gist.github.com/Jfortin1/72ef064469d1703c6b30
lighten <- function(color, factor = 0.5) {
  if ((factor > 1) | (factor < 0)) stop("factor needs to be within [0,1]")
  col <- col2rgb(color)
  col <- col + (255 - col)*factor
  col <- rgb(t(col), maxColorValue=255)
  col
}
