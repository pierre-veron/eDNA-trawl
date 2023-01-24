### ----------------------------------------------------------------------------
#   This scripts calculates principal coordinates analysis on the Jaccard 
#    distance between all sites (eDNA and trawling). 
#   It plots this PCoA along axis 1, 2 and 3, with color gradient and with
#    a map. 
#   It adds ellipses corresponding to depth (> or < 100m)
#
#   Author : Pierre Veron, pierre.veron.2017@polytechnique.org
#   Created 2022/03/21
### ----------------------------------------------------------------------------


plot_pcoa_eDNA_and_trawl <- function() {
  # Calculate Jaccard distance and PCoA ------------------------------------------
  Bet <- betapart::beta.pair(Incidence_matrix, index.family = "jaccard")
  Jtu <- as.matrix(Bet$beta.jtu)
  Jne <- as.matrix(Bet$beta.jne)
  Jac <- as.matrix(Bet$beta.jac)
  
  # PCoA
  PCoA <- dudi.pco(ade4::quasieuclid(as.dist(Jac)),scannf=FALSE, nf =8) # nf = 8
  Eigen <- PCoA$eig/sum(PCoA$eig)*100
  
  PCoA$li$A2 <- -PCoA$li$A2
  stat_plot <- addColorSpace(PCoA$li)
  saveRDS(stat_plot,"output/PCoA_colors.RDS")
  
  
  # Build the dataframe to be plotted --------------------------------------------
  stat_plot$latitude <- as.vector(sapply(rownames(stat_plot), function(station) {
    if (grepl("e", station)) {
      station <- gsub("e", "", station)
    } else {
      station <- gsub("t", "", station)
    }
    as.numeric(EDNA_metadata[which(EDNA_metadata$station == station), "Latitude"][1])
  }))
  
  stat_plot$longitude <- as.vector(sapply(rownames(stat_plot), function(station) {
    if (grepl("e", station)) {
      station <- gsub("e", "", station)
    } else {
      station <- gsub("t", "", station)
    }
    as.numeric(EDNA_metadata[which(EDNA_metadata$station == station), "Longitude"][1])
  }))
  
  stat_plot$latitude_edna <- as.vector(sapply(rownames(stat_plot), function(station) {
    if (grepl("e", station)) {
      station <- gsub("e", "", station)
    } else {
      station <- gsub("t", "", station)
    }
    as.numeric(EDNA_metadata[which(EDNA_metadata$station == station), "Latitude"][1])
  }))
  stat_plot$longitude_edna <- as.vector(sapply(rownames(stat_plot), function(station) {
    if (grepl("e", station)) {
      station <- gsub("e", "", station)
    } else {
      station <- gsub("t", "", station)
    }
    as.numeric(EDNA_metadata[which(EDNA_metadata$station == station), "Longitude"][1])
    
  }))
  
  stat_plot$type <- as.factor(as.vector(sapply(rownames(stat_plot), function(station) {
    if (grepl("e", station)) {
      "eDNA"
    } else {
      "trawling"
    }
  })))
  
  
  stat_plot$region <- as.factor(as.vector(sapply(rownames(stat_plot), function(station) {
    station <- gsub("t", "", station)
    station <- gsub("e", "", station)
    
    if (EDNA_metadata[which(EDNA_metadata$station == station),]$avg_depth[1] > 100) {
      "> 100 m"
    } else {
      "< 100 m"
    }
  })))
  
  # rename the stations 
  for (i in 1:length(Name_stations)) {
    rownames(stat_plot) <- gsub(names(Name_stations[i]), Name_stations[i], rownames(stat_plot))
  }
  
  stat_plot$pch <- as.vector(sapply(stat_plot$type, function(typ) {
    if (typ == "eDNA") {
      21
    } else {
      22
    }
  }))
  

  # text color (set manually to have good contrast)
  stat_plot$textcol <- as.vector(sapply(stat_plot$col, function(color) {
    rgb <- col2rgb(color) / 256
    bright <- (0.2126*rgb[1] + 0.7152*rgb[2] + 0.0722*rgb[3])
    if (bright > 0.65) {
      "black"
    } else {
      "white"
    }
  }))
  
  stat_plot$type_region <- as.factor(paste(as.character(stat_plot$type), as.character(stat_plot$region)))
  
  lbl <- gsub("e","", gsub("t", "", rownames(stat_plot)))
  
  pdf("figures/Fig_6_PCoA_map.pdf", width = 17/2.54, height = 8/2.54) 
  
  par(mar=c(2,2.1,1,0),mfrow=c(1,2))
  get_pcoa_plot(data=stat_plot[,c(1,2)], lbl = lbl,
                xlim=c(-0.4,0.4),ylim=c(-0.3,0.5),ax=c(1,2),Traits=stat_plot$type_region,
                xlab="",ylab="",col_bgf="grey76",colsclass=rep("black", 4),
                cex.txt = 0.4, offset = 0,
                col_vect=as.character(stat_plot$col),
                pos=NA, 
                col.text = stat_plot$textcol,
                cex.lab=0.3,cex.axis=0.6,ellipse=T, pch = stat_plot$pch,
                label_el = NULL,  cex.pts = 1.15) # to display labels on ellipses, set label_el = levels(stat_plot$type_region),
  xlab <- paste0("PCoA Axis 1 (", round(Eigen[1],1)," %)")
  ylab <- paste0("PCoA Axis 2 (", round(Eigen[2],1)," %)")
  title(xlab = xlab, line = 0.8, cex.lab = 0.7)            
  title(ylab = ylab, line = 1.4, cex.lab = 0.7)
  
  # adding ellipses labels manually
  text(0.02, 0.4, "trawl, < 100 m", adj = 0, cex = 0.7)
  text(0.22, -0.2, "trawl, > 100 m", adj = 0, cex = 0.7)
  text(-0.2, 0.03, "eDNA, < 100 m", adj = 1, cex = 0.7)
  text(-0.23, -0.285, "eDNA, > 100 m", adj = 0, cex = 0.7)
  mtext(text="a",side=3,adj = 0, cex = 1)
  legend("bottomleft", legend = c("eDNA", "trawling"), pch = c(21,22), cex =0.7,box.lwd = 0.5 , bg = "white" )
  # map
  par(mar=c(2,1.6,1,0.3)) 
  plot(0,0, type='n', xaxt='n', yaxt='n', xlab='', ylab='', bty='n')
  mtext(text="b",side=3,adj = 0, cex = 1)
  par(new = "T")
  pre_map(GISdata, xlim = c(-5,0), ylim = c(43,48))
  # plot points 
  dx <- 0.03
  dy <- 0.03
  
  points(x=stat_plot[which(stat_plot$type == "trawling"),]$longitude_edna-dx,
         y=stat_plot[which(stat_plot$type == "trawling"),]$latitude_edna-dy,
         bg=stat_plot[which(stat_plot$type == "trawling"),]$col,
         pch=22,
         col="black",
         cex=0.9,
         lwd=0.35)
  points(x=stat_plot[which(stat_plot$type == "eDNA"),]$longitude_edna+dx,
         y=stat_plot[which(stat_plot$type == "eDNA"),]$latitude_edna+dy,
         bg=stat_plot[which(stat_plot$type == "eDNA"),]$col,
         pch=21,
         col="black",
         cex=0.9,
         lwd=0.35)
  
  dx <- rep(0,nrow(stat_plot[which(stat_plot$type == "eDNA"),]))
  dy <- dx
  dx[7] <- -0.05
  dx[14] <- -0.1
  dy[10] <- 0.45
  dx[12] <- 0.1
  dy[12] <- 0.47
  dy[7] <- 0.3
  dy[14] <- 0.35
  dx[10] <- -0.05
  text(x=stat_plot[which(stat_plot$type == "eDNA"),]$longitude-0.1 + dx,
       y=stat_plot[which(stat_plot$type == "eDNA"),]$latitude-0.38 + dy,
       label=gsub("e", "", rownames(stat_plot[which(stat_plot$type == "eDNA"),])),
       pos=3,
       cex=0.8,
       offset=0.1)
  post_map(cex.axis = 0.6, add.northarrow = T, add.depth.legend = T, dpt.title.x = 0.08, dpt.title.y = 0.52,)
  
  dev.off()
  
  # Fig S6 ---------------------------------------------------------------------
  
  pdf(file= "figures/Fig_S6_PCoA_axes2_3.pdf",width=14/2.54,height=6/2.54)
  par(mar=c(2,2,1,0.8),mfrow=c(1,2), mgp = c(3,1,0))
  get_pcoa_plot(data=stat_plot[,c(1,3)], lbl = lbl,
                xlim=c(-0.4,0.4),ylim=c(-0.4,0.4),ax=c(1,2),Traits=stat_plot$region,
                xlab="",ylab="",col_bgf="grey76",colsclass=c("black", "black"),
                cex.txt = 0.27, offset = 0,
                col_vect=as.character(stat_plot$col),
                pos=NA,
                col.text = stat_plot$textcol,
                label_el = levels(stat_plot$region),
                cex.lab=0.3,cex.axis=0.5,ellipse=T, pch = stat_plot$pch)
  xlab <- paste0("PCoA Axis 1 (", round(Eigen[1],1)," %)")
  ylab <- paste0("PCoA Axis 3 (", round(Eigen[3],1)," %)")
  title(xlab = xlab, line = 0.8, cex.lab = 0.7)            
  title(ylab = ylab, line = 1.4, cex.lab = 0.7)
  mtext(text="a",side=3,adj = 0, cex = 1)
  legend("bottomleft", legend = c("eDNA", "trawling"), pch = c(21,22), cex =0.5,box.lwd = 0.5  )
  
  # map
  par(mar=c(2,2,1,0.8))
  get_pcoa_plot(data=stat_plot[,c(2,3)], lbl = lbl,
                xlim=c(-0.3,0.5),ylim=c(-0.4,0.4),ax=c(1,2),Traits=stat_plot$region,
                xlab="",ylab="",col_bgf="grey76",colsclass=c("black", "black"),
                cex.txt = 0.27, offset = 0,
                col_vect=as.character(stat_plot$col),
                pos=NA,
                col.text = stat_plot$textcol,
                label_el = levels(stat_plot$region),
                cex.lab=0.3,cex.axis=0.5,ellipse=T, pch = stat_plot$pch)
  xlab <- paste0("PCoA Axis 2 (", round(Eigen[2],1)," %)")
  ylab <- paste0("PCoA Axis 3 (", round(Eigen[3],1)," %)")
  title(xlab = xlab, line = 0.8, cex.lab = 0.7)            
  title(ylab = ylab, line = 1.4, cex.lab = 0.7)

  mtext(text="b",side=3,adj = 0, cex = 1)
  dev.off()
}

plot_pcoa_eDNA_and_trawl()

