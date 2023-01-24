### ----------------------------------------------------------------------------
#   Defines functions to plot colored PCoA.
#
#   Author : Camille Albouy, minor modif by Pierre Veron
#   Updated 2022/05/25 : add a function addColorSpace to implement the definition
#     of the color gradient and correction of a bug (discontinuity of the color 
#     gradient).
### ----------------------------------------------------------------------------

# Update : 
# The problems was in the function getDataColor : 
# 1. The split between two spaces was for x = 0, but this is not always the middle
#    of the range, so I changed it to (max + min)/2
# 2. Once the points are splitted, the bounding to [0,1] is done separately between both
#    halves of the set, so that the transformations are not consistent between
#    the two parts.
#    Example : before corrections 
#     |x     y   | would be splitted into a|-1 ,0.6| and c|1 , 0.3| 
#    a|-1    0.6 |                        b|-.2,0  |     d|.2, 0  |
#    b|-0.2  0   |
#    c|1     0.3 |   
#    d|0.2   0   |                        
#    
#                    and transformed in   a|0,  1| and   c|1, 1|
#                                         b|1,  0| and   d|0, 0|
#    we can see that points a and c have the same y coordinate after transformation
#    whereas they do not have the same y coordinate before (the transformation is 
#      not continuous). 
#
#   After correction, they are and transformed in   a|0  ,1| and     c|1  , 0.5|
#                                                   b|0.1,0| and     d|0.1, 0  |
#  The main function to use is addColorSpace   


#' addColotSpace
#'
#' @param data a dataframe withfirst two columns contain the first two axes
#' of the resulted PCoA 
#'
#' @return the same dataframe with a new column "col" containing the color in 
#' the colorspace
addColorSpace <- function(data) {
  roug <- c(1,0,0)
  jaun <- c(1,1,0)
  turq <- c(0,1,1)
  bleu <- c(0,0,1)
  mag  <- c(1,0,1)
  noir <- c(0,0,0)
  vert <- c(0,1,0)
  
  
  violetcl   <- biColorInter (col1 = bleu, col2 = roug, frac = 0.5)
  violet     <- biColorInter (col1 = bleu, col2 = roug, frac = 0.6)
  
  #vertcl     <- biColorInter (col1 = jaun, col2 =turq, frac = 0.2)
  vert     <- c(0.65,0.88,0.53)
  
  
  marron     <- biColorInter (col1 = roug, col2 = vert, frac = 0.5)
  rougefabli <- biColorInter ( col1 = roug, col2 = mag, frac = 0.25)
  
  vertcl <- c(0.37,0.68,0.37)
  
  
  jaunc <- c(1,0.996,0.58);bleuc <- c(0.58,0.788,1)
  rougec <- c(1,0.635,0.596);orange <- c(1,0.69,0.125)
  vert <- c(0.24,0.6,0.37)
  marron <- biColorInter(col1=rougec,col2=vert,frac=0.5)
  
  # add the color gradient
  Rich1_ratio1 <- getdatacolor(data=data,condition=">=0",col1=jaun,col2=orange,col3=bleuc, col4=vert)
  Rich1_ratio2 <- getdatacolor(data=data,condition= "<0",col1=roug,col2=jaun,col3=bleu,col4=bleuc)
  # Rich1_ratio1 <- getdatacolor(data=PCoA$li,condition=">=0",col1=orange,col3=bleuc,col2=jaun, col4=vert)
  # Rich1_ratio2 <- getdatacolor(data=PCoA$li,condition= "<0",col1=roug,col2=orange,col3=bleu,col4=bleuc)
  
  rbind(Rich1_ratio1,Rich1_ratio2)
}

get_pcoa_plot <- function(data=PCoA$li,xlim=c(-5,5),ylim=c(-7,5),ax=c(1,2),offset=0.4,
                          xlab="PC 1",ylab="PC 2",col_bgf="grey76",colsclass=c("#40AAD7","#08519F"),
                          Traits = as.factor(Vect_pos),col_vect,pos,label_el="",ellipse=F,
                          cex.txt=0.5,cex.lab=0.6,pos.labx=1.2,pos.laby=1.7,cex.axis=0.45,
                          lbl = NA, pch = 21,
                          col.text = "black",cex.pts = 0.8){
  
  plot(rnorm(100),col="white",ylim=ylim,xlim =xlim,type="n",axes=F,xlab="",ylab="",mgp=c(1.3,1,0),cex.lab=cex.lab)

  if(ellipse==T){
    s.class(data[ax[1]:ax[2]], Traits, wt = rep(1, length(Traits)), xax = 1, 
            yax = 2, cstar = 0, cellipse = 1.5, axesell = TRUE, 
            label = label_el, clabel = 0.7,cpoint = 0.7, pch = 20, 
            col = colsclass, xlim = NULL, ylim = NULL,lwd=0.5, 
            grid = FALSE, addaxes =F, origin = c(0,0), 
            include.origin = TRUE, sub = "", csub = 1, possub = "bottomleft", 
            cgrid = 1.3, pixmap = NULL, contour = NULL, area = NULL, add.plot = TRUE)
  }
  if (any(is.na(lbl))) {
  label1<- rownames(data)
  } else {
    label1 <- lbl
  }
  points(x=data[,ax[1]],y=data[,ax[2]],xlab="",pch=pch,col="black",bg=col_vect,cex=cex.pts,lwd=0.4)
  if (any(is.na(pos))) {
    text(x=data[,ax[1]],y=data[,ax[2]],label=label1,cex=cex.txt,adj = 0.5,col = col.text) 
  } else {
    text(x=data[,ax[1]],y=data[,ax[2]],label=label1,cex=cex.txt,
         pos = pos, offset = offset, col = col.text) 
  }
  
  # axis(side=1,cex.axis=cex.axis,at = seq(xlim[1],xlim[2],0.1),lwd=0.35,tcl=-0.25,bg="white",
  #      labels=as.character(round(seq(xlim[1],xlim[2],0.1),2)),mgp=c(0,0.05,0))
  # axis(side=2,cex.axis=cex.axis,las=2,at = seq(ylim[1],ylim[2],0.1),lwd=0.35,tcl=-0.25,bg="white",
  #      labels=as.character(round(seq(ylim[1],ylim[2],0.1),2)),mgp=c(0,0.5,0))
  
  axis(side=1,cex.axis=cex.axis,at = seq(xlim[1],xlim[2],0.1),lwd=0.35,tcl=-0.25,bg="white",
       labels=F, line = 0)
  axis(side=1,cex.axis=cex.axis,at = seq(xlim[1],xlim[2],0.1),lwd=0.35,tick = F,bg="white",
       labels=as.character(round(seq(xlim[1],xlim[2],0.1),2)), line = -1)
  
  axis(side=2,cex.axis=cex.axis,las=2,at = seq(ylim[1],ylim[2],0.1),lwd=0.35,tcl=-0.25,bg="white",
       labels=F, line = 0)
  axis(side=2,cex.axis=cex.axis,las=2,at = seq(ylim[1],ylim[2],0.1),lwd=0.35,tick = F,bg="white",
       labels=as.character(round(seq(ylim[1],ylim[2],0.1),2)), line = -0.6)
  
  
  mtext(text = xlab,side = 1,line=pos.labx,at = 0,cex=cex.lab)
  mtext(text = ylab,side = 2,line=pos.laby,at = 0,cex=cex.lab)
  abline(h=0,lty=2,lwd=0.5); abline(v=0,lty=2,lwd=0.5); 
  box(lwd=0.55)
} # end of get_pcoa_plot  

get_translate <- function (data=Data_linesST[[1]],x=0.0003,y=0){
  dd <- coordinates(data)[[1]][[1]]
  cbind(dd[,1]+x,dd[,2]+y)
  
}

s.class <- function (dfxy, fac, wt = rep(1, length(fac)), xax = 1, yax = 2, 
                     cstar = 1, cellipse = 1.5, axesell = TRUE, label = levels(fac), 
                     clabel = 1, cpoint = 1, pch = 20, col = rep(1, length(levels(fac))), 
                     xlim = NULL, ylim = NULL, grid = TRUE, addaxes = TRUE, origin = c(0, 
                                                                                       0), include.origin = TRUE, sub = "", csub = 1, possub = "bottomleft", 
                     cgrid = 1, pixmap = NULL, contour = NULL, area = NULL, add.plot = FALSE,lwd=0.2) 
{
  opar <- par(mar = par("mar"))
  par(mar = c(0.1, 0.1, 0.1, 0.1))
  on.exit(par(opar))
  dfxy <- data.frame(dfxy)
  if (!is.data.frame(dfxy)) 
    stop("Non convenient selection for dfxy")
  if (any(is.na(dfxy))) 
    stop("NA non implemented")
  if (!is.factor(fac)) 
    stop("factor expected for fac")
  dfdistri <- ade4:::fac2disj(fac) * wt
  coul <- col
  w1 <- unlist(lapply(dfdistri, sum))
  dfdistri <- t(t(dfdistri)/w1)
  coox <- as.matrix(t(dfdistri)) %*% dfxy[, xax]
  cooy <- as.matrix(t(dfdistri)) %*% dfxy[, yax]
  if (nrow(dfxy) != nrow(dfdistri)) 
    stop(paste("Non equal row numbers", nrow(dfxy), nrow(dfdistri)))
  coo <- scatterutil.base(dfxy = dfxy, xax = xax, yax = yax, 
                          xlim = xlim, ylim = ylim, grid = grid, addaxes = addaxes, 
                          cgrid = cgrid, include.origin = include.origin, origin = origin, 
                          sub = sub, csub = csub, possub = possub, pixmap = pixmap, 
                          contour = contour, area = area, add.plot = add.plot)
  if (cpoint > 0) 
    for (i in 1:ncol(dfdistri)) {
      pch <- rep(pch, length = nrow(dfxy))
      points(coo$x[dfdistri[, i] > 0], coo$y[dfdistri[, i] > 0], pch = pch[dfdistri[, i] > 0], cex = par("cex") * 
               cpoint, col = coul[i])
    }
  if (cstar > 0) 
    for (i in 1:ncol(dfdistri)) {
      scatterutil.star(coo$x, coo$y, dfdistri[, i], cstar = cstar, 
                       coul[i])
    }
  if (cellipse > 0) 
    for (i in 1:ncol(dfdistri)) {
      scatterutil.ellipseb(coo$x, coo$y, dfdistri[, i], 
                           cellipse = cellipse, axesell = axesell, coul[i],lwd=lwd)
    }
  if (clabel > 0) 
    scatterutil.eti(coox, cooy, label, clabel, coul = col)
  
  invisible(match.call())
}


scatterutil.ellipseb<- function (x, y, z, cellipse, axesell,lwd ,coul = rep(1, length(x))) 
{
  if (any(is.na(z))) 
    return(invisible())
  if (sum(z * z) == 0) 
    return(invisible())
  util.ellipse <- function(mx, my, vx, cxy, vy, coeff) {
    lig <- 100
    epsi <- 1e-10
    x <- 0
    y <- 0
    if (vx < 0) 
      vx <- 0
    if (vy < 0) 
      vy <- 0
    if (vx == 0 && vy == 0) 
      return(NULL)
    delta <- (vx - vy) * (vx - vy) + 4 * cxy * cxy
    delta <- sqrt(delta)
    l1 <- (vx + vy + delta)/2
    l2 <- vx + vy - l1
    if (l1 < 0) 
      l1 <- 0
    if (l2 < 0) 
      l2 <- 0
    l1 <- sqrt(l1)
    l2 <- sqrt(l2)
    test <- 0
    if (vx == 0) {
      a0 <- 0
      b0 <- 1
      test <- 1
    }
    if ((vy == 0) && (test == 0)) {
      a0 <- 1
      b0 <- 0
      test <- 1
    }
    if (((abs(cxy)) < epsi) && (test == 0)) {
      if (vx > vy) {
        a0 <- 1
        b0 <- 0
      }
      else {
        a0 <- 0
        b0 <- 1
      }
      test <- 1
    }
    if (test == 0) {
      a0 <- 1
      b0 <- (l1 * l1 - vx)/cxy
      norm <- sqrt(a0 * a0 + b0 * b0)
      a0 <- a0/norm
      b0 <- b0/norm
    }
    a1 <- 2 * pi/lig
    c11 <- coeff * a0 * l1
    c12 <- (-coeff) * b0 * l2
    c21 <- coeff * b0 * l1
    c22 <- coeff * a0 * l2
    angle <- 0
    for (i in 1:lig) {
      cosinus <- cos(angle)
      sinus <- sin(angle)
      x[i] <- mx + c11 * cosinus + c12 * sinus
      y[i] <- my + c21 * cosinus + c22 * sinus
      angle <- angle + a1
    }
    return(list(x = x, y = y, seg1 = c(mx + c11, my + c21, 
                                       mx - c11, my - c21), seg2 = c(mx + c12, my + c22, 
                                                                     mx - c12, my - c22)))
  }
  z <- z/sum(z)
  m1 <- sum(x * z)
  m2 <- sum(y * z)
  v1 <- sum((x - m1) * (x - m1) * z)
  v2 <- sum((y - m2) * (y - m2) * z)
  cxy <- sum((x - m1) * (y - m2) * z)
  ell <- util.ellipse(m1, m2, v1, cxy, v2, cellipse)
  if (is.null(ell)) 
    return(invisible())
  polygon(ell$x, ell$y, border = coul,lwd=lwd)
  if (axesell) 
    segments(ell$seg1[1], ell$seg1[2], ell$seg1[3], ell$seg1[4], 
             lty = 2, col = coul,lwd=lwd)
  if (axesell) 
    segments(ell$seg2[1], ell$seg2[2], ell$seg2[3], ell$seg2[4], 
             lty = 2, col = coul,lwd=lwd)
}




# mixted colors function 
biColorInter <- function(col1,col2,frac){
  
  #colors needs to be in RGB [0;1] space (will add for [0;255]
  #decomposing colors
  r1 <- col1[1]
  g1 <- col1[2]
  b1 <- col1[3]
  
  r2 <- col2[1]
  g2 <- col2[2]
  b2 <- col2[3]
  
  #deltas calculation
  Dr <- r2 - r1
  Dg <- g2 - g1
  Db <- b2 - b1
  
  #new r, g and b
  red = r1 + (Dr * frac)
  gre = g1 + (Dg * frac)
  blu = b1 + (Db * frac)
  
  unlist(c(red,gre,blu))
  
}#eo biColorInter

################### 4 colors interpolation

# col3         col4
#
#
# col1         col2


# four colors (RGB)


quadriColorInter <- function(colx1=roug,colx2=jaun,coly1=bleu,coly2=turq,fracx,fracy){
  
  finColX1 <- biColorInter(colx1,colx2,fracx)
  finColX2 <- biColorInter(coly1,coly2,fracx)
  
  fin <- biColorInter(finColX1,finColX2,fracy)
  
  rgb(fin[1],fin[2],fin[3])
  
}#eo quadriColorInter


getdatacolor <- function (data=ww ,condition = ">=0",col1, col2, col3, col4){
  med <- median(data[, 1])#(max(data[,1]) + min(data[,1])) / 2
  if(condition == ">=0"){
    w <- data[data[,1]>=med,]
    c <- w
    c[, 1] <- (c[, 1] - med) / (max(c) - med)
  } else{
    w <- data[data[,1]<med,]
    c <- w
    c[, 1] <- (c[, 1] - min(c[, 1])) / (med  - min(c[,1]))
  }
  c[, 2] <- (c[,2] - min(data[,2])) / (max(data[,2]) - min(data[,2]))

  col <- sapply(1:nrow(w),function(i){
    quadriColorInter (col1,col2,col3,col4,fracx=c[i,1],fracy=c[i,2])
  })
  
  cbind(w,col)
  
} # end of getdatacolor()


