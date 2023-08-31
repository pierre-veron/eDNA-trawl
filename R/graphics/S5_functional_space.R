### ----------------------------------------------------------------------------
#   Plots the functional space by site (with a difference in eDNA and trawling)
#
#   Author: Pierre Veron, pierre.veron.2017@polytechnique.org
#   Created: 2022/04/12
### ----------------------------------------------------------------------------

plot_fct_space_by_site <- function() {
  # Load and select one fct space ------------------------------------------------
  fct_space_list <- readRDS("output/fct_space.rds")

  i_rand <- 1 # sample(1:length(fct_space_list), size = 1)
  fct_space <- fct_space_list[[i_rand]]
  colnames(fct_space ) <- gsub("PC","PCoA axis ",colnames(fct_space))
  
  # Parameters -------------------------------------------------------------------
  dims <- c(1,2) # which axis of the fct space to be plotted
  
  # Plot initialisation ----------------------------------------------------------
  out_file <- "figures/Fig_S5_functional_space"
  ext <- ".pdf"
  
  image_width <- 17
  image_height <- 20
  res = 800
  
  if (ext == ".tif") {
    tiff(filename = paste0(out_file, ext),
         width = image_width,
         height = image_height, 
         units = 'cm',
         res = res)
  } else if (ext == ".pdf") {
    pdf(file = paste0(out_file, ext), 
        width = image_width / 2.54,
        height = image_height / 2.54)
  } else  {
    jpeg(filename = paste0(out_file, ext),
         width = image_width,
         height = image_height, 
         units = 'cm',
         res = res)
  }
  
  par(mfrow = c(5,3))
  par(oma = c(4, 4, 0, 0))
  xlim <- c(floor(10*min(fct_space[, dims[1]]))/10, 0.1 * ceiling(10*max(fct_space[, dims[1]])))
  ylim <- c(floor(10*min(fct_space[, dims[2]]))/10, 0.1 * ceiling(10*max(fct_space[, dims[2]])))
  # Loop by site -----------------------------------------------------------------
  for (i_site in 1:(0.5*nrow(Incidence_matrix))) {
    # init subplot
    par(mar = 0.5+c(0,0,0,0))
    plot(c(), c(), xlim = xlim, ylim = ylim, xlab = "", ylab = "", xaxt = "n", yaxt = "n")
    axis(1, labels = F)
    axis(2, labels = F)
    
    # select sub-functional space
    Sub_sp <- as.data.frame(fct_space[, dims])
    site_names <- paste0(names(Name_stations[i_site]), c('e', 't'))
    Det_species <- Incidence_matrix[site_names, ]
    rownames(Det_species) <- c("eDNA", "trawl")
    colnames(Det_species) <- gsub("\\.", " ", colnames(Det_species))
    Det_species <- Det_species[, which(colSums(Det_species) > 0)]
    Sub_sp <- Sub_sp[colnames(Det_species), ]
    eDNA_species <- colnames(Det_species[, which(Det_species["eDNA",] > 0)])
    trawl_species <- colnames(Det_species[, which(Det_species["trawl",] > 0)])
    
    # define color 
    Sub_sp$color <- NA
    Sub_sp[which(rownames(Sub_sp) %in% eDNA_species), "color"] <- color_edna
    Sub_sp[which(rownames(Sub_sp) %in% trawl_species), "color"] <- color_trawl
    Sub_sp[which(rownames(Sub_sp) %in% intersect(eDNA_species, trawl_species)), "color"] <- "black"
    
    # define convex hulls
    Sub_sp_edna <- Sub_sp[eDNA_species, dims]
    ch_edna <- chull(Sub_sp_edna)
    
    ch_edna <- Sub_sp_edna[c(ch_edna, ch_edna[1]),]
    
    
    
    Sub_sp_trawl <- Sub_sp[which(rownames(Sub_sp) %in% trawl_species), dims]
    ch_trawl <- chull(Sub_sp_trawl)
    ch_trawl <- Sub_sp_trawl[c(ch_trawl, ch_trawl[1]),]
    
    # draw convex hull
    polygon(x = ch_edna[,1], y = ch_edna[,2], border = color_edna, 
            col = adjustcolor(color_edna, alpha = 0.3),
            density = NULL)
    polygon(x = ch_trawl[,1], y = ch_trawl[,2], border = color_trawl, 
            col = adjustcolor(color_trawl, alpha = 0.3),
            density = NULL)
    
    # draw points
    points(x = Sub_sp[, 1], y = Sub_sp[, 2], col = Sub_sp$color, pch = 16, cex = 0.8)
    
    # add some text
    text(x = 0.53, y = 0.32, labels = i_site, cex = 2)
    
    x <- -0.29
    y <- -0.48
    text(x = x, y = y, labels = "FD q0  \nFDiv  \nFEve", cex = 0.85, adj = c(0,0))
    
    indic <- round(RESULTS_fct_both$Alpha[site_names, c("FD_q0", "Mean.fdiv", "Mean.feve")],2)
    text(x = x + 0.17, y = y, labels = paste0(indic[1, ], collapse = "\n"),cex = 0.85, col = color_edna, adj = c(0,0))
    text(x = x + 0.34, y = y, labels = paste0(indic[2, ], collapse = "\n"),cex = 0.85, col = color_trawl, adj = c(0,0))
    if (i_site == 15) {
      par(xpd = NA)
      legend("bottomright", inset = c(0.04,-0.3), pch = 16, col = c(color_edna, color_trawl, "black"), 
             legend = c("eDNA", "trawling","both"), title = "detected by", bg = "white")
    }
  }
  
  mtext(colnames(fct_space)[dims[1]], side = 1, outer = TRUE, line = 2)
  mtext(colnames(fct_space)[dims[2]], side = 2, outer = TRUE, line = 2)
  
  dev.off()
  
}

plot_fct_space_by_site()