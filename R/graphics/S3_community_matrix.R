### ----------------------------------------------------------------------------
#   This script plots a community matrix 
#
#   Author : Pierre Veron, pierre.veron.2017@polytechnique.org
### ----------------------------------------------------------------------------

plot_community_matrix <- function() {
  # Sort species by taxonomic proximity ------------------------------------------
  pth_tax <- paste(Reftax$class, Reftax$order, Reftax$family, 
                   Reftax$genus, Reftax$species, sep = "/")
  Reftax <- Reftax[order(pth_tax), ]
  
  # Color matrix
  color_none <- "gray95"
  
  # Building local incidence matrix 
  Inc_mat_eDNA <- Incidence_matrix_eDNA
  Inc_mat_trawl <- Incidence_matrix_trawl
  
  colnames(Inc_mat_eDNA) <- gsub("\\.", " ", colnames(Inc_mat_eDNA))
  colnames(Inc_mat_trawl) <- gsub("\\.", " ", colnames(Inc_mat_trawl))
  
  # Sort incidence matrices sites
  Inc_mat_eDNA <- Inc_mat_eDNA[paste0(names(Name_stations), "e"), ]
  Inc_mat_eDNA <- Inc_mat_eDNA[,Reftax$taxname]
  Inc_mat_trawl <- Inc_mat_trawl[paste0(names(Name_stations), "t"), ]
  Inc_mat_trawl <- Inc_mat_trawl[,Reftax$taxname]
  Color_matrix <- Inc_mat_eDNA
  rownames(Color_matrix) <- gsub("e", "", rownames(Color_matrix))
  
  for (i_site in 1:nrow(Color_matrix)) {
    Color_matrix[i_site, which(Inc_mat_eDNA[i_site, ] == 1 & Inc_mat_trawl[i_site,] == 1)] <- "black"
    Color_matrix[i_site, which(Inc_mat_eDNA[i_site, ] == 1 & Inc_mat_trawl[i_site,] == 0)] <- color_edna
    Color_matrix[i_site, which(Inc_mat_eDNA[i_site, ] == 0 & Inc_mat_trawl[i_site,] == 1)] <- color_trawl
    Color_matrix[i_site, which(Inc_mat_eDNA[i_site, ] == 0 & Inc_mat_trawl[i_site,] == 0)] <- color_none
  }
  
  Color_matrix <- Color_matrix[, Reftax$taxname]
  
  # Draw plot
  out_file <- "figures/Fig_S3_community_matrix"
  ext <- ".pdf"
  
  image_width <- 30
  image_height <- 34
  res = 800
  
  pdf(file = paste0(out_file, ext),
      width = image_width / 2.54,
      height = image_height / 2.54)
  
  par(mar = c(0,0,0,0))
  
  dx <- 0.95
  dy <- 0.85
  
  plot(c(), c(), xlim = c(-8, 3.5+nrow(Color_matrix)), ylim = c(0.4,  ncol(Color_matrix)-4), xaxt = 'n', yaxt = 'n',
       bty="n", xlab = "", ylab = "")
  
  a <- sapply(1:nrow(Color_matrix), function (x) {
    b <- sapply(1:ncol(Color_matrix), function(y) {
      rect(xleft = x, xright = x + dx, ybottom = ncol(Color_matrix) - y, ytop = ncol(Color_matrix) - y + dy, border = "black",
           col = Color_matrix[x,y])
    })
  })
  
  a <- sapply(1:nrow(Reftax), function(i) {
    det_eDNA <- sum(Inc_mat_eDNA[, i]) > 0
    det_trawl <- sum(Inc_mat_trawl[, i]) > 0
    col <- if (det_eDNA & det_trawl) {
      "black"
    } else if (det_eDNA) {
      color_edna
    } else if (det_trawl) {
      color_trawl
    } else {
      color_none
    }
    rect(xleft = -0.6, xright = -0.6+dx, ybottom = ncol(Color_matrix) - i, ytop = ncol(Color_matrix) - i + dy, border = "black",
         col = col)
  })
  
  # Name of the taxa
  a <- sapply(1:ncol(Color_matrix), function(y) {
    text(x = 0.75-1.6, y = ncol(Color_matrix) - y + dy/2, labels = colnames(Color_matrix)[y], adj = c(1,0.35), cex = 0.75)
  })
  
  # Name of the orders
  a <- sapply(unique(Reftax$order), function(ord) {
    ytop <- ncol(Color_matrix) - min(which(Reftax$order == ord))
    ybot <-  ncol(Color_matrix) -max(which(Reftax$order == ord))
    if (ybot > 0) {
      lines(x = c(-8,-0.75), y = c(ybot, ybot), lty = 3, col = "gray", lwd = 0.8)
    }
    if (ytop < ncol(Color_matrix) - 1) {
      lines(x = c(-8,-0.75), y = 1 + c(ytop, ytop), lty = 3, col = "gray", lwd = 0.8)
    }
    if (ytop - ybot >= 1) {
      text(x = -5.3, y = (1+ytop+ybot)/2, labels = ord, adj = c(1, 0.5), srt = 0, cex = 0.8)
    }
    else {
      #
    }
  })
  
  # Name of the classes 
  a <- sapply(unique(Reftax$class), function(cl) {
    ytop <- ncol(Color_matrix) - min(which(Reftax$class == cl))
    ybot <-  ncol(Color_matrix) -max(which(Reftax$class == cl))
    if (ybot > 0) {
      lines(x = c(-10,-0.75), y = c(ybot, ybot), lty = 1, col = "black", lwd = 0.8)
    }
    if (ytop < ncol(Color_matrix) - 1) {
      lines(x = c(-10,-0.75), y = 1 + c(ytop, ytop), lty = 1, col = "black", lwd = 0.8)
    }
    if (ytop - ybot >= 1) {
      text(x = -8.6, y = (1+ytop+ybot)/2, labels = cl, adj = 0.5, srt = 90, cex = 1)
    }
    else {
      #
    }
  })
  
  # Name of the sites
  text(x = (1+nrow(Color_matrix))/2, y = -2.5, labels = "Site", adj = 0.5)
  
  a <- sapply(1:nrow(Color_matrix), function(x) {
    text(x = 0.5 + x, y = -0.85, labels = Name_stations[rownames(Color_matrix)[x]], adj = 0.5, cex = 1)
  })
  
  text(x = -0.1, y = -0.85, labels = "all", adj = 0.5)
  
  # Legend 
  legend("right",  legend = c("eDNA", "trawling", "both", "none"),
         fill = c(color_edna, color_trawl, "black", color_none), inset = c(0.01, 0), ncol = 1,
         title = "Detection by")
  dev.off()
}

plot_community_matrix()