### ----------------------------------------------------------------------------
#    Defines useful functions
### ----------------------------------------------------------------------------




#' short_species_name
#'
#' @param spec a name of species (character with two words)
#'
#' @return a char with shorter name
#' @export
#'
#' @examples
#' short_species_name("Chelon ramada")
#' [1] "C. ramada"
#' 
#' short_species_name("Chelon")
#' [1] "Chelon"
#' 
#' short_species_name("Stromis boa boa")
#' [1] "S. boa"
#' 
short_species_name <- function(spec) {
  if (length(spec) == 1) {
    a <- strsplit(spec, " ")[[1]]
    if (length(a) < 2) {
      spec
    } else {
      paste0(substr(a[1],1,1), ". ", a[2])
    }
  } else if (length(spec) > 1) {
    c(short_species_name(spec[1]), short_species_name(spec[-1]) )
  } else {
    spec
  }
}



#' addImg
#' Adds an image to an existing graph
#' Taken from :  https://ogeek.cn/qa/?qa=958668/
#' @param obj an image file imported as an array (e.g. png::readPNG, jpeg::readJPEG)
#' @param x mid x coordinate for image
#' @param y mid y coordinate for image
#' @param width width of image (in x coordinate units)
#' @param interpolate (passed to graphics::rasterImage) A logical vector 
#' (or scalar) indicating whether to apply linear interpolation to the image 
#' when drawing. 
#' @param angle rotation of the image
addImg <- function(
  obj, 
  x = NULL,
  y = NULL,
  width = NULL,
  interpolate = TRUE,
  angle = 0
){
  if(is.null(x) | is.null(y) | is.null(width)){stop("Must provide args 'x', 'y', and 'width'")}
  USR <- par()$usr # A vector of the form c(x1, x2, y1, y2) giving the extremes of the user coordinates of the plotting region
  PIN <- par()$pin # The current plot dimensions, (width, height), in inches
  DIM <- dim(obj) # number of x-y pixels for the image
  ARp <- DIM[1]/DIM[2] # pixel aspect ratio (y/x)
  WIDi <- width/(USR[2]-USR[1])*PIN[1] # convert width units to inches
  HEIi <- WIDi * ARp # height in inches
  HEIu <- HEIi/PIN[2]*(USR[4]-USR[3]) # height in units
  rasterImage(image = obj, 
              xleft = x-(width/2), xright = x+(width/2),
              ybottom = y-(HEIu/2), ytop = y+(HEIu/2), 
              interpolate = interpolate,
              angle = angle)
}


#' randomize
#' 
#' Randomly permutes the lines of a dataframe, but keeps the names of the rows. 
#' This function is designed to work for traits table.
#' Types of the columns must be numeric, character or factor. No need that all the 
#' columns are same type
#'
#' @param df a dataframe containing traits values for instance, with species name
#' zs rownames and traits names as colnames.
#'
#' @return a new dataframe with randomy permuted rows but rownames in the same order. 
#' Colnames and types do not change, provided they are all numeric, factor or char
#' 
randomize <- function(df) {
  n <- nrow(df)
  m <- ncol(df)
  ord <- sample(1:n, size = n)
  new <- t(sapply(ord, function(i) {
    sapply(1:m, function(j) {
      as.character(df[i,j])
    })
  }))
  rownames(new) <- rownames(df)
  colnames(new) <- colnames(df)
  new <- as.data.frame(new)
  classes <- sapply(df, class)
  for (col in colnames(df)) { # default is character
    if (classes[col] == "numeric") {
      new[, col] <- as.numeric(new[, col])
    } else if (classes[col] == "factor") {
      new[, col] <- as.factor(new[, col])
    }
  }
  new
}



#' circle_intersect
#' Calculates the area of the intersection between two circles
#' @param r1 radius of circle 1
#' @param r2 radius of circle 2
#' @param d distance between two centers
#'
#' @return the area of the intersection
#' @export
#'
#' @examples
circle_intersect <- function(r1, r2, d) {
  if (r2 > r1) { # r1 should be greater than r2
    r <- r1
    r1 <- r2
    r2 <- r
  }
  if (d <= r1-r2) {
    pi*r2^2
  } else if (d >= r2+r1) {
    0
  } else { # equation from https://mathworld.wolfram.com/Circle-CircleIntersection.html
    r <- r2
    R <- r1
    A <- r^2 * acos((d^2+r^2-R^2) / (2*d*r))
    B <- R^2 * acos((d^2+R^2-r^2) / (2*d*R))
    C <- -0.5 * sqrt((-d+r+R)*(d+r-R)*(d-r+R)*(d+r+R))
    A+B+C
  }
}

#' reciproq_circle_intersect 
#' Calculates the distance between two circles such that the area of intersection
#' equals the parameter area.
#' If area = 0 : the distance is arbitrarly set to 1.2 * (r1+r2) 
#' (disks are disjoint)
#' If area is larger than the area of the smallest circle, this is impossible, 
#' but case treated as (3). 
#' If area is equal to the area of the smallest circle, the smallest disk is inside the 
#' largest, returns |r1-r2| (arbitrarly)
#' If area is >0 but smaller than the area of the little circle, returns a value
#' found by dichotomy with a precision 0.01*r1
#' 
#' This function is used to draw Venn diagrams. 
#' @param r1 radius of circle 1
#' @param r2 radius of circle 2
#' @param area desired area
#'
#' @return the distance between the two centers of the circles
#' @export
#'
#' @examples
reciproq_circle_intersect <- function(r1, r2, area) {
  if (r2 > r1) { # r1 should be greater than r2
    r <- r1
    r1 <- r2
    r2 <- r
  }
  # find d so that area = circle_intersect(r1,r2,d)
  if (area <= 0) {
    1.2 * (r1+r2)
  } else if (area >= pi * r2^2) {
    r1-r2
  } else {
    d_low <- r1-r2
    d_hig <- r1+r2
    eps <- 0.01*r1
    while(d_hig - d_low  > eps) {
      d_med <- (d_hig + d_low)/2
      if (circle_intersect(r1,r2,d_med) > area) {
        d_low <- d_med
      } else {
        d_hig <- d_med
      }
    }
    (d_hig+d_low) /2
  }
}

#' my_venn_diag
#' Draws a Venn diagram of two sets.
#' The area of the circles are proportional to the size of the sets
#' A plot must be defined before executing this function, for instance with the
#' following command 
#' plot(c(), c(), xlim = c(-3,3), ylim = c(-3,3), xaxt = "n", yaxt = "n", axes = F)
#' @param setlist a list of two vectors 
#' @param colors a vector containing 2 colors in a format accepted by 
#' grDevices::col2rgb
#' Accepted format : Hexadecimal string, positive integer, named color
#' @param alpha a float between 0 and 1 : transparency of the circles
#' @param x position x of the Venn diagram 
#' @param y position y of the Venn diagram (the center of the diagram is set at 
#' the middle between the two centers of the circles)
#' @param size positive float : size factor, the radius of a circle that would 
#' contain 1 element
#' @param angle float in radians : orientation of the Venn diagram
#' @param print.size a boolean : if T the function displays the size of the sets
#' and the size of the intersection (unless the intersection is emtpy or equal to 
#' one of the sets)
#' @param text.lbl list of characters or NULL : 
#' * if NULL labels will be displayed if print.size with the number of elements 
#' in each part
#' * if list of characters : labels to be displayed 
#' text.lbl[[1]] and text.lbl[[2]] correspond to the exclusive sets of setlist
#' text.lbl[[3]] is the label for the intersection
#' @param text.offset a float for text position >1 for text outside the circles, 
#' <1 for text inside
#' @param text.cex size of font 
#'
#' @return a list containing the graphical parameters of the circles 
#' @export
#'
#' @examples
my_venn_diag <- function(setlist, 
                         colors = c("red", "blue"), 
                         alpha = 0.5, 
                         x = 0,
                         y = 0, 
                         size = 0.1, 
                         angle = 0, 
                         print.size = T,
                         text.lbl = NULL, 
                         text.offset = 1.3,
                         text.cex = 1.3) {
  inter <- length(intersect(setlist[[1]], setlist[[2]]))
  r1 <- size * sqrt(length(setlist[[1]]))
  r2 <- size * sqrt(length(setlist[[2]]))
  d <- reciproq_circle_intersect(r1, r2, pi * size^2 * inter)
  colors_alpha <- sapply(colors, function(col) {
    adjustcolor(col, alpha = alpha)
  })
  dx <- d/2 * cos(angle)
  dy <- d/2 * sin(angle)
  draw.circle(x = x - dx, y = y-dy, radius = r1, col = colors_alpha[1])
  draw.circle(x = x+ dx, y = y+dy, radius = r2, col = colors_alpha[2])
  if (is.null(text.lbl) & print.size) {
    text(x = x- dx - text.offset*r1*cos(angle), 
         y = y - dy - text.offset*r1*sin(angle), 
         label = length(setlist[[1]]),
         col = colors[1], cex =text.cex)
    text(x = x+ dx + text.offset*r2*cos(angle), 
         y = y + dy + text.offset*r2*sin(angle), 
         label = length(setlist[[2]]),
         col = colors[2], cex = text.cex)
    if (inter > 0 & inter < min(length(setlist[[1]]),length(setlist[[2]]))) {
      text(x = x + cos(angle) * (r1-r2)/2, 
           y = y + sin(angle) * (r1-r2)/2, 
           label = inter,
           col = "black", cex = text.cex)
    }
  } else if (!is.null(text.lbl)) {
    text.lbl <- lapply(text.lbl, function(t) {
      if (length(t) == 0) {
        ""
      } else {
        paste0(t, collapse = "\n")
      }
    })
    text(x = x- dx - text.offset*r1*cos(angle), 
         y = y - dy - text.offset*r1*sin(angle), 
         label = text.lbl[[1]],
         col = colors[1], cex =text.cex)
    text(x = x+ dx + text.offset*r2*cos(angle), 
         y = y + dy + text.offset*r2*sin(angle), 
         label = text.lbl[[2]],
         col = colors[2], cex = text.cex)
    text(x = x + cos(angle) * (r1-r2)/2, 
         y = y + sin(angle) * (r1-r2)/2, 
         label = text.lbl[[3]],
         col = "black", cex = text.cex)
  }
  list("circle1" = c("radius" = r1, "center.x" = x - dx, "center.y" = y-dy),
       "circle2" = c("radius" = r2, "center.x" = x + dx, "center.y" = y+dy))
}

#' my_clip
#' 
#' Bounds a vector between lower and upper bounds 
#' @param x a numeric vector
#' @param lower_bound float or NA
#' @param upper_bound float or NA
#'
#' @return a vector of same size as x with values bounded between lower_bound 
#' and upper bound. If one of bounds is set to NA, the vector is not bounded in
#' this direction. 
#' @export
#'
#' @examples
my_clip <- function(x, lower_bound = NA, upper_bound = NA) {
  out <- pmax(pmin(x, upper_bound, na.rm = T), lower_bound, na.rm = T)
  out[which(is.na(x))] <- NA
  out
}

#' graduation
#'
#' @param x float 
#' @param y float position of the bottom left corner of the graduation
#' @param width float width of the graduation
#' @param height float height of the graduation
#' @param minval value corresponding to the left of the rule
#' @param maxval value corresponding to the right of the rule
#' @param dgrad interval between two graduations
#' @param lbl vector of the labels for the graduations or NA (if NA, the label 
#' are automatically calculated)
#' @param startgrad value of the first graduation
#' @param text.offset float, vertical position of the text, relative to height
#' @param text.cex cex for the text
#' @param outer.param a list of named parameters that will be passed as arguments 
#' of the rect function (or an empty list for default parameters)
#' @param grad.param a  list of named parameters for the tick of the graduation
#' that will be passed as arguments to the line function (or an empty list for 
#' default parameters)
#'
#' @return
#' @export
#'
#' @examples
graduation <- function (x = 0,
                        y = 0, 
                        width = 10,
                        height = 1,
                        minval = 0,
                        maxval = 10,
                        dgrad = 1, 
                        lbl = NA,
                        startgrad = NA,
                        text.offset = -0.5,
                        text.cex = 0.5,
                        outer.param = list("col" = "white",
                                           "border" = "black",
                                           "lwd" = 2),
                        grad.param = list())  {
  default.outer.param <- list("xleft" = x, "xright" = x + width, "ybottom" = y, "ytop"= y+height)
  paramlist <- unique(c(names(default.outer.param), names(outer.param)))
  rect.param <- lapply(paramlist, function(param) {
    if (is.null(outer.param[[param]])) {
      default.outer.param[[param]]
    } else (
      outer.param[[param]]
    )
  })
  names(rect.param) <- paramlist 
  do.call(rect, rect.param)
  if (is.na(startgrad)) {
    startgrad <- minval
  }
  grads <- seq(startgrad, maxval, dgrad)
  ngrad <-length(grads)
  dx <- width / ngrad
  if (any(is.na(lbl))) {
    lbl <- grads
  }
  
  default.grad.param <- list("x" = NA, "y" = c(y, y+height))
  paramlist <- unique(c(names(default.grad.param), names(grad.param)))
  line.param <- lapply(paramlist, function(param) {
    if (is.null(grad.param[[param]])) {
      default.grad.param[[param]]
    } else (
      grad.param[[param]]
    )
  })
  names(line.param) <- paramlist
  for (i in 1:ngrad) {
    x_t <- x + grads[i] * width / (maxval-minval)
    line.param[["x"]] <- c(x_t, x_t)
    do.call(lines, line.param)
    text(x = x_t, y = y + text.offset * height, lbl[i], cex = text.cex, adj = c(0.5,0))
  }
}

triangle_mark <- function(x, y, width, ...) {
  x_p <- c(x, x-width/2, x+width/2, x)
  dy <- width * sqrt(3)/2
  y_p <- c(y, y + dy, y+dy, y)
  polygon(x = x_p, y = y_p, ...)
}


#' p_to_star
#' Turns a p-value into a string with following rule 
#' *** < 0.001 < ** < 0.01 < * < 0.05 < ""
#' @param p 
#'
#' @return a character: "", "*", "**" or "***"
#' @export
#'
#' @examples
p_to_star <- function(p) {
  if (p < 0.001) {
    "***"
  } else if (p < 0.01) {
    "**"
  } else if (p < 0.05) {
    "*"
  } else {
    ""
  }
}

#' ses_to_signif
#' Turns a SES value into a string with rule
#' --- < -3.29 < -- <-2.58 < - < -1.96 < empty < 1.96 < + < 2.58 < ++ < 3.29 < +++
#' @param x
#'
#' @return a character: "+", "++", "+++", "-", "--", "---" or ""
#' @export
#'
#' @examples
ses_to_signif <- function(x) {
  if (x < -3.29) {
    "---"
  } else if (x < -2.58) {
    "--"
  } else if (x < -1.96) {
    "-"
  } else if (x <= 1.96) {
    ""
  } else if (x <= 2.58) {
    "+"
  } else if (x <= 3.29) {
    "++"
  } else {
    "+++"
  }
}

#' draws a pie slice on an existing graph
#'
#' @param x 
#' @param y position of the center of the pie slice
#' @param r radius of the pie slice
#' @param nsteps (int) number of segments (the pie slice is approximated with a polygon)
#' @param angle vector with values of the begining and the end of the pie slice 
#' in radians 
#' @param lwd linewidth
#' @param ... other arguments passed to the polygon function (color...)
#'
#' @return
#' @export
#'
#' @examples
pie_slice <- function(x,y,r,nsteps=1000,angle = c(0,pi),lwd = 0.7, ...){  
  rs <- seq(angle[1],angle[2],len=nsteps) 
  xc <- c(x+r*cos(rs) ,x)
  yc <- c(y+r*sin(rs) ,y)
  polygon(xc,yc,lwd = lwd, ...) 
} 




pre_map <- function(GISdata, xlim = c(-5,0), ylim = c(43,47)) { # code to execute before plotting points on the map
  
  breaksbathy <- c(-5000,-4500,-4000,-3500,-3000,-2500,-2000,-1500,-1400,-1300,-1200,-1100,
                   -1000,-900,-800,-700,-600,-500,-450,-400,-350,-300,-250,-200,-150,-100,-90,
                   -80,-70,-60,-50,-40,-20,-10,-5,0)
  plot(c(), c(),col="white",type="n",axes=F,xlab="",ylab="",xlim=xlim,ylim=ylim,cex.main=1.3, main="",font=2, asp = 1)
  colour <- colorRampPalette(c("#4BA1EC","white"))(length(breaksbathy)-1) 
  image(GISdata$Bathy,col=colour,breaks=breaksbathy,axes=F,add=TRUE)
  plot(GISdata$Countries,add=TRUE,col="grey40",lwd=0.3,border="black")
  plot(GISdata$Countries[which(GISdata$Countries@data[,"NAME_SORT"]=="France"),],col="grey60",add=T,lwd=0.3)
  plot(GISdata$Rivers,add=T,col="grey10",lwd=0.3)
  
}

post_map <- function(xtickslab = T, ytickslab = T, add.northarrow = F, 
                     add.depth.legend = F, cex.axis = 0.85, dpt.lgd.x = c(0.13,0.16), 
                     dpt.lgd.y = c(0.25,0.48), dpt.title.x = 0.08, dpt.title.y = 0.52, 
                     dpt.title.cex = 0.55, line.xlab = -1, line.ylab=-0.6) { # code to execute after plotting a map
  breaksbathy <- c(-5000,-4500,-4000,-3500,-3000,-2500,-2000,-1500,-1400,-1300,-1200,-1100,
                   -1000,-900,-800,-700,-600,-500,-450,-400,-350,-300,-250,-200,-150,-100,-90,
                   -80,-70,-60,-50,-40,-20,-10,-5,0)
  long_labels <- (-6):(0)
  lat_labels <- 43:48
  colour <- colorRampPalette(c("#4BA1EC","white"))(length(breaksbathy)-1) 
  
  axis(side=1,cex.axis=0.6,lwd=0.35,tcl=-0.25,bg="white",labels=F,
       at=long_labels,line = 0)
  axis(side=2,las=2,cex.axis=0.6,lwd=0.35,tcl=-0.25,bg="white",labels=F,
       at=lat_labels, line = 0)
  
  if (xtickslab) {
    axis(side=1,cex.axis=cex.axis,lwd=0.35,tick = F,bg="white",labels=paste0(long_labels, "°"),
         at=long_labels,line = line.xlab)
  } 
  if (ytickslab) {
    axis(side=2,las=2,cex.axis=cex.axis,lwd=0.35,tick = F,bg="white",labels=paste0(lat_labels, "°"),
         at=lat_labels, line = line.ylab)
  } 
  
  if (add.northarrow) {
    prettymapr::addnortharrow(pos="topleft",padin=c(0.15,0.15),scale=0.35,cols=c("white","black"),text.col="black")
  }
  
  if (add.depth.legend) {
    colorlegend(posx=dpt.lgd.x,posy=dpt.lgd.y,col=rev(colour),lwd=0.3,zlim=c(0,-5200),cex= 0.8 * cex.axis,
                zlevels=6,dz=-1000,tcl=-0.2,offset=0.2)
    text(x=dpt.title.x,y=dpt.title.y,"Depth (m)",cex=dpt.title.cex, font = 2)
  }
  box(lwd = 0.5)
}

round_pcoa_interia <- function(s) { # "PC1 (39.7 %)" --> "PC1 (40%)"
  split <- strsplit(s, '(', fixed = T)[[1]]
  value <- as.numeric(strsplit(split[2], ' ', fixed = T)[[1]][1])
  paste0(split[1], "(", round(value), "%)")
}

