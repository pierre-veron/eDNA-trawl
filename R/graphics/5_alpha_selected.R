plot_alpha_div_selected <- function() {
  selected_sites <- c(2,5,8,15)
  out_file <- 'figures/Fig_5_alpha_selected'
  ext <- ".pdf"
  
  image_width <- 17
  image_height <- 13
  res = 800
  
  color_none <- "gray60"
  color_both <- "black"
  
  # Prepare incidence matrix -----------------------------------------------------
  Inc_mat <- Incidence_matrix
  colnames(Inc_mat) <- gsub("\\.", " ", colnames(Inc_mat))
  
  Inc_mat_eDNA <- Incidence_matrix_eDNA
  colnames(Inc_mat_eDNA) <- gsub("\\.", " ", colnames(Inc_mat_eDNA))
  
  Inc_mat_trawl <- Incidence_matrix_trawl
  colnames(Inc_mat_trawl) <- gsub("\\.", " ", colnames(Inc_mat_trawl))
  
  Inc_mat_merged <- Incidence_matrix_merged
  colnames(Inc_mat_merged) <- gsub("\\.", " ", colnames(Inc_mat_merged))
  # Sort species by taxonomic proximity ------------------------------------------
  pth_tax <- paste(Reftax$class, Reftax$order, Reftax$family, 
                   Reftax$genus, Reftax$species, sep = "/")
  Reftax <- Reftax[order(pth_tax), ]
  
  
  # Sort incidence matrices sites
  Inc_mat_eDNA <- Inc_mat_eDNA[paste0(names(Name_stations), "e"), ]
  Inc_mat_eDNA <- Inc_mat_eDNA[, Reftax$taxname]
  Inc_mat_trawl <- Inc_mat_trawl[paste0(names(Name_stations), "t"), ]
  Inc_mat_trawl <- Inc_mat_trawl[,Reftax$taxname]
  Color_matrix <- Inc_mat_eDNA
  rownames(Color_matrix) <- gsub("e", "", rownames(Color_matrix))
  
  # Create list of taxa
  tax_sites <- colSums(Inc_mat_merged[names(Name_stations)[selected_sites], ])
  tax_sites <- tax_sites[which(tax_sites > 0)]
  list_taxa <- sort(names(tax_sites))
  names(list_taxa) <- list_taxa
  abb_taxa <- sapply(list_taxa, function(tax) {
    if (grepl(" ", tax)) {
      words <- strsplit(tax, " ")[[1]]
      paste0(substr(words[1], 1,1), toupper(substr(words[2], 1,1)), substr(words[2], 2,2))
    } else {
      substr(tax, 1,3)
    }
  })
  
  # avoid duplicates
  abb_taxa["Chelon"] <- "Chn"
  abb_taxa["Chelidonichthys"] <- "Chs"
  abb_taxa["Scorpaena scrofa"] <- "SScr"
  abb_taxa["Scomber scombrus"] <- "SSco"
  abb_taxa["Sarda sarda"] <- "SaSa"
  abb_taxa["Scomberesox saurus"] <- "ScSa"
  for (i_site in 1:nrow(Color_matrix)) {
    Color_matrix[i_site, which(Inc_mat_eDNA[i_site, ] == 1 & Inc_mat_trawl[i_site,] == 1)] <- "black"
    Color_matrix[i_site, which(Inc_mat_eDNA[i_site, ] == 1 & Inc_mat_trawl[i_site,] == 0)] <- color_edna
    Color_matrix[i_site, which(Inc_mat_eDNA[i_site, ] == 0 & Inc_mat_trawl[i_site,] == 1)] <- color_trawl
    Color_matrix[i_site, which(Inc_mat_eDNA[i_site, ] == 0 & Inc_mat_trawl[i_site,] == 0)] <- color_none
  }
  
  Color_matrix <- Color_matrix[, Reftax$taxname]
  
  # Prepare phylo tree -----------------------------------------------------------
  iconic_species <- c("Leucoraja naevus", "Mola mola", "Pleuronectes platessa", 
                      "Engraulis", "Lamna nasus", "Hippocampus hippocampus", "Syngnathus acus",
                      "Capros aper", "Myctophum punctatum", "Trigla lyra", "Chelon ramada", 
                      "Merluccius", "Dicentrarchus labrax", "Lophius", "Trachurus trachurus",
                      "Chimaera monstrosa", "Anguilla")
  iconic_species <- as.vector(sapply(iconic_species, function(spec) {
    strsplit(spec, " ")[[1]][1]
  }))
  itree <- 30 # sample(1:length(Trees), 1)
  phyl_tree <-Trees[[itree]] # we take one tree randomly
  df_ph <- Reftax
  df_ph$no.rd.taxname <- df_ph$taxname
  # random draw
  Genuses <- df_ph[which(df_ph$rank == "genus"), "genus"]
  for (gen in Genuses) {
    df_ph[which(df_ph$genus == gen), "taxname"] <- 
      sample(Reftax_species[which(Reftax_species$genus == gen), "species"],1)
  }
  
  detected <- gsub(" ", "_", df_ph$taxname)
  phyl_tree <- keep.tip(phyl_tree, tip =  detected)
  
  

  phyl_symbols_size = sapply(phyl_tree$tip.label, function(lab){
    0.8
  })
  
  # Prepare functional space -----------------------------------------------------
  fct_space_list <- readRDS("output/fct_space.rds")
  
  i_rand <- 1 # 84 
  fct_space <- as.data.frame(fct_space_list[[i_rand]])
  
  # Plot initialisation ----------------------------------------------------------
  
  pdf(file = paste0(out_file, ext),
      width = image_width / 2.54,
      height = image_height / 2.54,
      title = "Alpha div selected")
  
  
  
  par(mfrow = c(3,length(selected_sites)), oma = c(0.2,0,0,0))
  
  # PART 1 : Taxonomic -----------------------------------------------------------
  for (site in selected_sites) {
    par(mar = c(0,0,0.4,0))
    station <- names(Name_stations)[which(Name_stations == site)]
    Color_taxa <- Color_matrix[station, ]
    Color_taxa <- Color_taxa[which(Color_taxa != color_none)]
    
    Orders <- unique(Reftax[which(Reftax$taxname %in% colnames(Color_taxa)),
                                  "order"])
    ylim <- c(0,27)
    xlim <- c(-6,6)
    plot(c(), c(), xlim = xlim, ylim = ylim, xaxt = 'n', yaxt = 'n', axes = F,
         xlab = "", ylab = "")
    x <- xlim[1]
    y <- par("usr")[4]
    dy <- 2.5 # size of the title 
    
    rect(xleft = -6, xright = 6, ybottom = y - 0.9*dy, ytop = y, col = "gray80", border = NA )
    text(x = 0, y = y - 0.45 * dy, labels = paste0("Site ", site), adj = c(0.5,0.5), cex = 1.3)
    
    y <- y - 1.4*dy
    ymax <- y
    ymin <- par("usr")[3]
    
    dx <- 4
    
    # split by orders and select species detected in that site from each order
    for (i_ord in 1:length(Orders)) {
      ord <- Orders[i_ord]
      Species <- Reftax[which(Reftax$order == ord), "taxname"]
      Color_taxa_order <- Color_taxa[which(colnames(Color_taxa) %in% Species)]
      if (length(Color_taxa_order) > 0) {
        rect(xleft = x+0.04*dx, xright = x+ 0.96 * dx, ytop = y, ybottom = y-length(Color_taxa_order),
             col = "white", border = "black", lwd = 1.5)
        for (i_tax in 1:length(Color_taxa_order)) {
          
          rect(xleft = x+0.04*dx, xright = x + 0.96 * dx, ytop = y, ybottom = y-1, 
               col = as.character(Color_taxa_order[i_tax]), 
               border = NA)
          text(x = x+0.05*dx, y = y-0.5, labels = short_species_name(colnames(Color_taxa_order[i_tax])),
               adj = c(0,0.5), cex = 0.52, col = "white", font = 3)
          y <- y-1
        }
        y <- y - 0.7
      }
      
      size_next_ord <- if (i_ord == length(Orders)) {
        0
      } else {
        ord <- Orders[i_ord + 1]
        Species <- Reftax[which(Reftax$order == ord), "taxname"]
        length(Color_taxa[which(colnames(Color_taxa) %in% Species)])
      }
      if (y-size_next_ord-0.7 < ymin) {
        x <- x+dx
        y <- ymax
      }
    }
    
    # legend 
    if (site == selected_sites[1]) {
      x <- 2.2
      y <- ymin
      rect(x,y, x+4.2,y+9.9, col = "white", border = "black", lwd = 1)
      text(x + 2.1, y + 8.8, "Detected by", font = 2, adj = 0.5, cex = 0.8)
      text(x +2.1, y +6.5, "eDNA", col = color_edna, adj = 0.5, cex = 0.8)
      text(x +2.1, y +4.7, "trawling", col = color_trawl, adj = 0.5, cex = 0.8)
      text(x +2.1, y +2.9, "both", col = "black", adj = 0.5, cex = 0.8)
      text(x +2.1, y +1.1, "none", col = color_none, adj = 0.5, cex = 0.8)
    }
    
  }
  
  # PART 2 : Phylogenetic ------------------------------------------------------
  for (site in selected_sites) {
    
    par(mar = c(0,0,0,0))
    #df_ph <- Rd_reftax
    station <- names(Name_stations)[which(Name_stations == site)]
    df_ph$color <- "gray"
    det_edna <- colnames(Inc_mat[,which(Inc_mat[paste0(station, 'e'),] == 1)])
    det_trawl <- colnames(Inc_mat[,which(Inc_mat[paste0(station, 't'),] == 1)])
    det_both <- intersect(det_edna, det_trawl)
    det_edna_only <- setdiff(det_edna, det_both)
    det_trawl_only <- setdiff(det_trawl, det_both)
    
    df_ph[which(df_ph$no.rd.taxname %in% det_edna_only),"color"] <- color_edna
    df_ph[which(df_ph$no.rd.taxname %in% det_trawl_only),"color"] <- color_trawl
    df_ph[which(df_ph$no.rd.taxname %in% det_both),"color"] <- color_both
    phyl_colors = sapply(phyl_tree$tip.label, function(lab){
      lab_sp <- gsub("_", " ", lab)
      as.character(df_ph[which(df_ph$taxname == lab_sp),"color"])
    })
    
    phyl_symbols = sapply(phyl_tree$tip.label, function(lab){
      lab_sp <- gsub("_", " ", lab)
      
      real_taxa <- df_ph[which(df_ph$taxname == lab_sp), "no.rd.taxname"]
      
      if (real_taxa  %in% union(det_trawl, det_edna)) {
        16
      } else {
        NA
      }
    })
    
    phyl_labels <- sapply(phyl_tree$tip.label, function(lab) {
      tax <- gsub("_", " ", lab)
      real_taxa <- df_ph[which(df_ph$taxname == tax), "no.rd.taxname"]
      
      
      if (real_taxa  %in% union(det_trawl, det_edna)) {
        abb_taxa[tax]
      } else {
        ""
      }
    })
    
    edges_eDNA <- get_edges_parent(phyl_tree, which(gsub("_", " ", phyl_tree$tip.label) %in% 
                                                      df_ph[which(df_ph$no.rd.taxname %in% det_edna), "taxname"]))
    edges_trawl <- get_edges_parent(phyl_tree, which(gsub("_", " ", phyl_tree$tip.label) %in%
                                                       df_ph[which(df_ph$no.rd.taxname %in% det_trawl), "taxname"]))
    phyl_edge_colors <- rep(color_none, nrow(phyl_tree$edge))
    phyl_edge_colors[intersect(edges_eDNA, edges_trawl)] <- "black"
    phyl_edge_colors[setdiff(edges_eDNA, edges_trawl)] <- color_edna
    phyl_edge_colors[setdiff(edges_trawl, edges_eDNA)] <- color_trawl
    
    
    # plot the tree
    marg <- 450
    plot.phylo(phyl_tree, type = 'fan', tip.color = phyl_colors, show.tip.label = F, no.margin = T,
               rotate.tree = 0, x.lim =c(-marg,marg), y.lim = c(-marg,marg),
               edge.color = phyl_edge_colors)
    tiplabels(pch = phyl_symbols, col = phyl_colors, cex = phyl_symbols_size)
    
    
    images_names <- c("shark", "mugiliformes", "perciformes", "hippocampus", 
                      "scorpaeniformes", "moroniformes", "caproiformes", "tetradontiformes",
                      "lophiiformes", "gadiformes", "clupeiformes", "chimaera", "ray",
                      "anguilliformes")
    angles <- c(-4, 11, 55,76,117,148, 178, 185, 195,235, 287, 312, 320,  303 )
    width <- rep(100, length(images_names))
    width[4] <- 35
    width[8] <- 80
    
    delta_radius <- rep(0, length(images_names))
    delta_radius[4] <- -18
    delta_radius[5] <- -10
    delta_radius[11] <- -15
    delta_radius[12] <- -10
    
    if (F) { # add images
      sapply(1:length(images_names), function(i) {
        img <- readPNG(here("data", "icons", paste0(images_names[i], ".png")))
        angle <- pi/180* angles[i]
        r <- (1.05*marg + delta_radius[i])
        addImg(img, x = r*cos(angle), y = r*sin(angle), width = width[i], angle = 0)
      })
    }
    

    # graduation line for clustering signal
    x <- -200
    y <- 130
    w <- 400
    
    lines(x = c(x-0.1*w,x+1.1*w), y = c(y,y), lwd = 0.7)
    lines(x = c(x,x), y = c(y,y+15), lwd = 0.7)
    lines(x = c(x+w,x+w), y = c(y,y+15), lwd = 0.7)
    text(x = x, y = y-50, labels = "clustered", cex = 0.65, adj = c(0.5,0))
    text(x = x+w, y = y-50, labels = "dispersed", cex = 0.65, adj = c(0.5,0))
    
    val_edna <- results_phylo_signal_alpha[site,"eDNA.D-value.Obs"]
    pval_edna <- results_phylo_signal_alpha[site,"eDNA.p_random"]
    triangle_mark(x = x + val_edna *w, y = y+15, width = 60, col = color_edna, border = color_edna)
    signif <- if (pval_edna <0.001) {"***"} else if (pval_edna < 0.01) {"**"} else if (pval_edna < 0.05) {"*"} else {""}
    text(x = x + val_edna *w, y = y + 100, adj = 0.5, labels=  signif, col = color_edna, cex = 1.4)
    
    val_trawl <- results_phylo_signal_alpha[site,"trawl.D-value.Obs"]
    pval_trawl <- results_phylo_signal_alpha[site,"trawl.p_random"]
    triangle_mark(x = x + val_trawl *w, y = y+ 15, width = 60, col = color_trawl, border = color_trawl)
    signif <- if (pval_trawl <0.001) {"***"} else if (pval_trawl < 0.01) {"**"} else if (pval_trawl < 0.05) {"*"} else {""}
    text(x = x + val_trawl *w, y = y + 100, adj = 0.5, labels=  signif, col = color_trawl, cex = 1.4)
  }
  
  # PART 3 : Functional ----------------------------------------------------------
  dims <- c(1,2)
  xlim <- range(fct_space[, dims[1]])
  ylim <- range(fct_space[, dims[2]])

  
  par(mar = c(0,0.5,0,0.5))
  for (site in selected_sites) {
    print(site)
    fct_space$color <- color_none
    fct_space$size <- 0.5
    
    
    station <- names(Name_stations)[which(Name_stations == site)]
    det_edna <- colnames(Inc_mat[,which(Inc_mat[paste0(station, 'e'),] == 1)])
    det_trawl <- colnames(Inc_mat[,which(Inc_mat[paste0(station, 't'),] == 1)])
    det_both <- intersect(det_edna, det_trawl)
    det_edna_only <- setdiff(det_edna, det_both)
    det_trawl_only <- setdiff(det_trawl, det_both)
    
    fct_space[det_edna_only, "color"] <- color_edna
    fct_space[det_trawl_only, "color"] <- color_trawl
    fct_space[det_both, "color"] <- color_both
    fct_space[union(det_edna, det_trawl), "size"] <- 1
    
    # define convex hulls
    Sub_sp_edna <- fct_space[det_edna, dims]
    ch_edna <- chull(Sub_sp_edna)
    ch_edna <- Sub_sp_edna[c(ch_edna, ch_edna[1]),]
    
    Sub_sp_trawl <- fct_space[det_trawl, dims]
    ch_trawl <- chull(Sub_sp_trawl)
    ch_trawl <- Sub_sp_trawl[c(ch_trawl, ch_trawl[1]),]
    
    
    
    # plot
    plot(c(), c(), xlim = xlim, ylim = ylim, xlab = "", ylab = "", xaxt = "n", yaxt = "n",
         axes =T)
    
    
    # draw convex hull
    polygon(x = ch_edna[,1], y = ch_edna[,2], border = color_edna, 
            col = adjustcolor(color_edna, alpha = 0.3),
            density = NULL, lwd = 0.5)
    polygon(x = ch_trawl[,1], y = ch_trawl[,2], border = color_trawl, 
            col = adjustcolor(color_trawl, alpha = 0.3),
            density = NULL, lwd = 0.5)
    
    # first plot undetected points
    df <- fct_space[which(fct_space$color == color_none),]

    
    df <- fct_space[which(fct_space$color != color_none),]

    text(x = df[, dims[1]], y = df[, dims[2]], adj = 0.5, labels = abb_taxa[rownames(df)], col = df$color, cex = 0.4 )
    dx_ch <- 0
    dy_ch <- 0
    col_ch <- df[rownames(ch_edna)[-nrow(ch_edna)], "color"]
    text(x = ch_edna[-nrow(ch_edna),1] + dx_ch, y = ch_edna[-nrow(ch_edna), 2] + dy_ch, labels = abb_taxa[rownames(ch_edna)[-nrow(ch_edna)]], cex = 0.4, adj = 0.5, col = col_ch)
    
    col_ch <- df[rownames(ch_trawl)[-nrow(ch_trawl)], "color"]
    text(x = ch_trawl[-nrow(ch_trawl),1] + dx_ch, y = ch_trawl[-nrow(ch_trawl), 2] + dy_ch, labels = abb_taxa[rownames(ch_trawl)[-nrow(ch_trawl)]], cex = 0.4, adj = 0.5, col = col_ch)
    axis(1, labels = F, tck = 0.03)
    axis(2, labels = F, tck = 0.03)
    if (site == selected_sites[1]) {
      mtext(colnames(fct_space)[dims[1]], side = 1, outer = FALSE, line = -1.35, cex = 0.6, adj = 0.4)
      mtext(colnames(fct_space)[dims[2]], side = 2, outer = FALSE, line = -1.35, cex = 0.6, adj = 0.92)
    }
  }
  
  dev.off()
}

plot_alpha_div_selected()