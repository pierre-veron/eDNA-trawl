plot_figure1 <- function() {
  
  # General graphical parameters -------------------------------------------------
  out_file <- "figures/Fig_1_gamma_div"
  ext <- ".pdf"
  image_width <- 17
  image_height <- 15
  res = 800
  
  color_rays <- "#44CF6C" 
  color_sharks <- "#9EE6B2"
  color_actino <- c("#B19206", "#F9DC5C", "#FBEA9D")
  color_holo <- "#877EBB"
  color_both <- "grey26"
  
  pch_edna <- 21
  pch_trawl <- 22
  
  # Prepare map ------------------------------------------------------------------
  Points_edna <- EDNA_metadata[which(EDNA_metadata$replicate == 1), 
                               c("station", "Latitude", "Longitude")]
  Points_edna$name <- Name_stations[Points_edna$station] 
  
  Points_trawl <- TRAWL_metadata
  Points_trawl$name <- Name_stations[Points_trawl$station]
  
  # Prepare taxo tree ------------------------------------------------------------
  df <- Reftax
  df$pathString <- paste("Eugn.", # for Eugnasthomata, the common clade to all detected species  
                         df$class,
                         df$subclass, 
                         df$infraclass,
                         df$order,                       
                         df$family, 
                         df$genus,
                         df$species,
                         sep = "/")
  df$end_label <- ''
  for (irow in 1:nrow(df)) {
    pth <- strsplit(df$pathString[irow], split = "/")[[1]]
    if (pth[8] == 'NA') { # no species but a genus
      df$pathString[irow] <- paste(pth[1], pth[2], pth[3],pth[4], pth[5], pth[5], pth[6], pth[7], sep = "/")
      df$end_label[irow] <- pth[7]
    } else { # there is a species 
      df$end_label[irow] <- pth[8]
    }
  }
  
  df$pathString <- gsub("NA", "", df$pathString)
  
  df$color <- ""
  
  df$det_by_eDNA <- as.vector(sapply(df$taxname, function(tax) {
    tax <- gsub(" ", ".", tax)
    sum(Incidence_matrix[grep("e", rownames(Incidence_matrix)),tax]) > 0
  }))
  
  df$det_by_trawl <- as.vector(sapply(df$taxname, function(tax) {
    tax <- gsub(" ", ".", tax)
    sum(Incidence_matrix[grep("t", rownames(Incidence_matrix)),tax]) > 0
  }))
  
  df$det_by <- as.vector(sapply(rownames(df), function(i) {
    if ((df[i, "det_by_eDNA"]) & (df[i, "det_by_trawl"])) {
      "both"
    } else if (df[i, "det_by_eDNA"]) {
      "eDNA"
    } else if (df[i, "det_by_trawl"]) {
      "trawl"
    }
  }))
  df[which(df$det_by == "eDNA"),'color'] <- color_edna
  df[which(df$det_by == "trawl"),'color'] <- color_trawl
  df[which(df$det_by == "both"),'color'] <- color_both
  
  # create tree
  tax_tree <- as.Node(df)
  tax_tree <- as.phylo(tax_tree)
  
  # define lists of parameters that are readable by plot.phylo
  tax_colors = sapply(tax_tree$tip.label, function(lab){
    lab_sp <- gsub("_", " ", lab)
    as.character(df[which(df$end_label == lab_sp),"color"])
  })
  
  tax_symbols = sapply(tax_tree$tip.label, function(lab){
    16
  })
  
  tax_symbols_size = sapply(tax_tree$tip.label, function(lab){
    0.8
  })
  
  tax_ring_cat <- sapply(tax_tree$tip.label, function(lab) {
    lab_sp <- gsub("_", " ", lab)
    infra <- df[which(df$taxname == lab_sp), "infraclass"]
    if (!is.na(infra)) {
      infra
    } else {
      df[which(df$taxname == lab_sp), "subclass"]
    }
  })
  
  Orders_actino <- unique(df[which(df$class == "Actinopterygii"), "order"])
  colors_orders_actino <- sapply(1:length(Orders_actino), function(i) {
    color_actino[1+i%%length(color_actino)]
  })
  names(colors_orders_actino) <- Orders_actino
  
  tax_ring_col <- sapply(tax_tree$tip.label, function(lab) {
    lab_sp <- gsub("_", " ", lab)
    taxa <- df[which(df$taxname == lab_sp), ]
    if (taxa$class == "Actinopterygii") {
      colors_orders_actino[taxa$order]
    } else if (taxa$class == "Holocephali"){
      color_holo
    } else if (taxa$infraclass == "Batoidea") {
      color_rays
    } else if (taxa$infraclass == "Selachii") {
      color_sharks
    } else {
      "white"
    }
  })
  
  
  edges_eDNA <- get_edges_parent(tax_tree, which(gsub("_", " ", tax_tree$tip.label) %in%
                                                   df[which(df$det_by_eDNA), "end_label"]))
  edges_trawl <- get_edges_parent(tax_tree, which(gsub("_", " ", tax_tree$tip.label) %in%
                                                    df[which(df$det_by_trawl), "end_label"]))
  tax_edge_colors <- rep("black", nrow(tax_tree$edge))
  tax_edge_colors[intersect(edges_eDNA, edges_trawl)] <- "black"
  tax_edge_colors[setdiff(edges_eDNA, edges_trawl)] <- color_edna
  tax_edge_colors[setdiff(edges_trawl, edges_eDNA)] <- color_trawl
  
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
  df_ph <- df
  
  # random draw
  Genuses <- df_ph[which(df_ph$rank == "genus"), "genus"]
  for (gen in Genuses) {
    df_ph[which(df_ph$genus == gen), "taxname"] <- 
      Reftax_species[which(Reftax_species$genus == gen), "species"][1]
  }
  
  detected <- gsub(" ", "_", df_ph$taxname)
  phyl_tree <- keep.tip(phyl_tree, tip = detected)
  
  df_ph$color <- ""
  df_ph[which(df_ph$det_by == "eDNA"),'color'] <- color_edna
  df_ph[which(df_ph$det_by == "trawl"),'color'] <- color_trawl
  df_ph[which(df_ph$det_by == "both"),'color'] <- color_both
  
  phyl_colors = sapply(phyl_tree$tip.label, function(lab){
    lab_sp <- gsub("_", " ", lab)
    as.character(df_ph[which(df_ph$taxname == lab_sp),"color"])
  })
  
  phyl_labels = sapply(phyl_tree$tip.label, function(lab){
    lab_sp <- gsub("_", " ", lab)
    label <- ""
    for (spec in iconic_species) {
      if (length(grep(spec, lab_sp)) > 0) {
        label <- spec # short_species_name(lab_sp)
      }
    }
    label
  })
  
  phyl_symbols = sapply(phyl_tree$tip.label, function(lab){ 16})
  phyl_symbols_size = sapply(phyl_tree$tip.label, function(lab){0.8})
  
  edges_eDNA <- get_edges_parent(phyl_tree, which(gsub("_", " ", phyl_tree$tip.label) %in%
                                                    df_ph[which(df_ph$det_by_eDNA), "taxname"]))
  edges_trawl <- get_edges_parent(phyl_tree, which(gsub("_", " ", phyl_tree$tip.label) %in%
                                                     df_ph[which(df_ph$det_by_trawl), "taxname"]))
  phyl_edge_colors <- rep("black", nrow(phyl_tree$edge))
  phyl_edge_colors[intersect(edges_eDNA, edges_trawl)] <- "black"
  phyl_edge_colors[setdiff(edges_eDNA, edges_trawl)] <- color_edna
  phyl_edge_colors[setdiff(edges_trawl, edges_eDNA)] <- color_trawl
  
  # Prepare functional space -----------------------------------------------------
  path_fct_space <- "output/fct_space.rds"
  
  Fct_space <- as.data.frame(readRDS(path_fct_space)[[4]])
  colnames(Fct_space) <- gsub("PC","PCoA axis ",colnames(Fct_space))
  
  Fct_space[df$taxname,"det_by"] <- df$det_by
  
  Fct_space[,1] <- -Fct_space[,1] # just reversing the first axis for nice graph
  
  Fct_space$color <- ""
  Fct_space[which(Fct_space$det_by == "eDNA"), 'color'] <- color_edna
  Fct_space[which(Fct_space$det_by == "trawl"), 'color'] <- color_trawl
  Fct_space[which(Fct_space$det_by == "both"), 'color'] <- color_both
  
  Fct_space$pch <- 16
  Fct_space$cex <- 0.8
  
  # Plot initialisation ----------------------------------------------------------
  pdf(file = paste0(out_file, ext),
      width = image_width / 2.54,
      height = image_height / 2.54)
  
  
  # Map --------------------------------------------------------------------------
  par(mar=c(2,2,1,0.8),mfrow=c(2,2), mgp = c(3,1,0))

  plot(0,0, type='n', xaxt='n', yaxt='n', xlab='', ylab='', bty='n')
  mtext(text="a",side=3,adj = 0, cex = 1)
  par(new = "T")
  pre_map(GISdata = GISdata, xlim = c(-6,0), ylim = c(43,48))
  # plot points 
  points(x=Points_trawl$Longitude,
         y=Points_trawl$Latitude,
         bg=color_trawl,
         pch=pch_trawl,
         col="black",
         cex=0.75,
         lwd=0.35)
  
  points(x=Points_edna$Longitude,
         y=Points_edna$Latitude,
         bg=color_edna,
         pch=pch_edna,
         col="black",
         cex=0.7,
         lwd=0.35)
  
  text(x=Points_edna$Longitude-0.1,
       y=Points_edna$Latitude+0.1,
       label=Points_edna$name,
       pos=3,
       cex=0.95,
       offset=0.1)
  post_map(add.northarrow = T, add.depth.legend = T,dpt.lgd.x = c(0.14,0.16), 
           dpt.lgd.y = c(0.35,0.58), dpt.title.x = 0.06, dpt.title.y = 0.62, dpt.title.cex = 0.6)
  
  legend(x = -0.01, y = 0.27, legend = c("eDNA", "trawling"), 
         pch = c(pch_edna, pch_trawl), 
         pt.bg = c(color_edna, color_trawl),
         box.lwd = 0.5,
         cex = 0.95, title = "Sampling sites", bg = "white")
  
  # Taxonomic tree ---------------------------------------------------------------
  par(mar=c(2,2,1,0.8), mgp = c(3,1,0))
  
  plot(0,0, type='n', xaxt='n', yaxt='n', xlab='', ylab='', bty='n')
  mtext(text="b",side=3,adj = 0, cex = 1)
  par(new = "T")
  plot.phylo(tax_tree, type = 'fan', tip.color = "white", show.tip.label = F, no.margin = T,
             rotate.tree = 90, x.lim = c(-130,130), y.lim = c(-130,130),
             edge.color = "white")
  #tiplabels(pch = tax_symbols, col = tax_colors, cex = tax_symbols_size)
  ring(13, tax_tree, col = tax_colors, offset = 15)
  
  # add images
  images_names <- c("ray" ,"shark", "clupeiformes", "anguilliformes", "gobiiformes",
                    "argentiniformes", "pleuronectiformes", "tetradontiformes",
                    "beryciformes", "caproiformes", "perciformes", "myctophiformes",
                    "scorpaeniformes", "mugiliformes", "gadiformes", "moroniformes",
                    "syngnathiformes", "lophiiformes", "carangiformes", "chimaera")
  
  width <- rep(35, length(images_names))
  delta_radius <- rep(0, length(images_names))
  width[8] <- 25
  width[11] <- 35
  width[18] <- 25
  width[17] <- 38
  delta_radius[3] <- -5
  delta_radius[11] <- -15
  delta_radius[13] <- -10
  delta_radius[14] <- -10
  delta_radius[17] <- 10
  delta_radius[18] <- -12
  
  angles <- c(50, 75, 98, 120, 150 ,158, 180, 205, 218, 263, 272, 285, 302,
              314, 338, 3,13.5, 20, 32, 37) #positions of images in degree

  sapply(1:length(images_names), function(i) {
    #img <- readPNG(here("data", "icons", paste0(images_names[i], ".png")))
    imgw <- readPNG(here("data", "icons", paste0(images_names[i],"_W", ".png")))
    angle <- pi/180* angles[i]
    r <- (126.5 + delta_radius[i])
    addImg(imgw, x = r*cos(angle), y = r*sin(angle), width = width[i], angle = 0)
    #addImg(img, x = r*cos(angle), y = r*sin(angle), width = width[i], angle = 0)
    
 })
     #})
  angle = pi/180*55
  text(x = 145*cos(angle), y = 145*sin(angle), labels = "Elasmobranchii", col = "black",
       cex = 0.75, srt = -30) 
  angle = pi/180*225
  text(x = 145*cos(angle), y = 145*sin(angle), labels = "Actinopt.", col = "black",
       cex = 0.75, srt = -45) 
  angle = pi/180*35
  text(x = 155*cos(angle), y = 155*sin(angle), labels = "C. monstrosa", col = "black",
       cex = 0.6, srt = 0, font = 3) 
  # include venn diagram
  #venn <- readPNG(here("figures", "venn_diagram", "manual_venn_diag.png"))
  set_edna <- df[which(df$det_by_eDNA), "taxname"]
  set_trawl <- df[which(df$det_by_trawl), "taxname"]
  par(xpd = NA)
  vd <- my_venn_diag(list(set_edna, set_trawl), colors = c(color_edna, color_trawl), 
                     size = 8, alpha = 0.5, angle = pi/4, 
                     text.offset = 1.1, x = 0, y = 0, text.cex = 1)
  angle <- -pi/4
  par(xpd = F)
  
  
  # Phylogenetic tree ------------------------------------------------------------
  par(mar=c(2,2,1,0.8), mgp = c(3,1,0))
  
  marg <- 500
  plot(0,0, type='n', xaxt='n', yaxt='n', xlab='', ylab='', bty='n')
  mtext(text="c",side=3,adj = 0, cex = 1)
  par(new = "T")
  age <- max(node.depth.edgelength(phyl_tree))
  
  # draw a black phyl tree for the scale
  plot.phylo(phyl_tree, type = 'fan', tip.color = NA, show.tip.label = F, no.margin = T,
             rotate.tree = 0, x.lim =c(-marg,marg), y.lim = c(-marg,marg),
             edge.color = NA, plot = F)
  
  scales <- c(50,100,150,200,250,300,350,400)
  for (s in scales) {draw.circle(x=0,y=0,radius=age- s,border="gray",lty=3)}
  scales <- c(200,300,400)
  for (s in scales) {
    rect(xleft = -35, xright = 35, ybottom = age-s-20, ytop = age-s+20, border = "white",
         col = "white")
    text(x = 0, y = age - s, labels = paste0(s, " Ma"), col = "black", cex = 0.6)
  }
  par(new = "T")
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
  
  sapply(1:length(images_names), function(i) {
    img <- readPNG(here("data", "icons", paste0(images_names[i], ".png")))
    angle <- pi/180* angles[i]
    r <- (1.05*marg + delta_radius[i])
    addImg(img, x = r*cos(angle), y = r*sin(angle),width=width[i],angle=0)
    #addImg(1-img, x = r*cos(angle), y = r*sin(angle),width=width[i]+20,angle=0)
  })
  
  # Functional space -------------------------------------------------------------
  par(mar=c(2,2,1,0.8), mgp = c(3,1,0))
  
  plot(0,0, type='n', xaxt='n', yaxt='n', xlab='', ylab='', bty='n')
  mtext(text="d",side=3,adj = 0, cex = 1)
  par(new = "T", mar=c(1.5,1.5,1,0.8) )
  
  # the two following functions aim at encapsulating the functions to plot
  plot1 <- function (expression) { Hmisc::subplot(expression,x=c(-1.1,1.1),
                   y = c(-1.05,0)) }
  
  plot2 <- function(expression) { Hmisc::subplot(expression, x = c(-1.1,1.1),
                   y = c(0,1.05)) }
  
  plot(0,0, xlim=c(-1,1),ylim=c(-1,1),type='n',xaxt='n',yaxt='n',xlab='',ylab='',bty='n')
  
  to_plot <- function(dm, show.xlab = T, show.legend = F){ # dm : dimensions to plot
    xlim <- c(floor(10*min(Fct_space[, dm[1]]))/10, 0.1 * ceiling(10*max(Fct_space[, dm[1]])))
    ylim <- c(floor(10*min(Fct_space[, dm[2]]))/10, 0.1 * ceiling(10*max(Fct_space[, dm[2]])))
    xlab <- round_pcoa_interia(colnames(Fct_space)[dm[1]]) #paste0("Dim ", dm[1], " (", round(Eigen_fct[dm[1]],1), " %)")
    ylab <- round_pcoa_interia(colnames(Fct_space)[dm[2]]) #paste0("Dim ", dm[2], " (", round(Eigen_fct[dm[2]],1), " %)")
    
    # build convex hull
    Sub_fct_space <- Fct_space[which(Fct_space$det_by %in% c("eDNA", "both")), dm]
    ch_edna <- chull(Sub_fct_space)
    ch_edna <- Sub_fct_space[c(ch_edna, ch_edna[1]),]
    
    Sub_fct_space <- Fct_space[which(Fct_space$det_by %in% c("trawl", "both")), dm]
    ch_trawl <- chull(Sub_fct_space)
    ch_trawl <- Sub_fct_space[c(ch_trawl, ch_trawl[1]),]
    Sub_fct_space <- Fct_space[, dm]
    ch_all <- chull(Sub_fct_space)
    ch_all <- Sub_fct_space[c(ch_all, ch_all[1]),]
    if (show.xlab) {
      xaxt = NULL
    } else {
      xaxt = "n"
    }
    
    plot(c(), c(), xlim = xlim, ylim = ylim, xlab = "", ylab = "", cex.axis = 0.7, 
         xaxt = "n", yaxt = "n")
    if (show.xlab) {
      axis(1, tck = 0.04, labels = F)
    }
    axis(2, tck = 0.04, labels = F)
    mtext(ylab, side = 2, line = 0.5, cex = 0.5)
    if (show.xlab) {
      mtext(xlab, side = 1, line = 0.5, cex = 0.5)
    }
    polygon(x = ch_edna[,1], y = ch_edna[,2], border = color_edna, 
            col = adjustcolor(color_edna, alpha = 0.3),
            density = NULL)
    polygon(x = ch_trawl[,1], y = ch_trawl[,2], border = color_trawl, 
            col = adjustcolor(color_trawl, alpha = 0.3),
            density = NULL)
    print(ch_all)
    if (2 %in% dm) {
      dx_ch <- 0.05*c(-1,1,-1,0,0.5)
      dy_ch <- 0.07*c(-1,0,1,1,1)
    } else {
      dx_ch <- 0.05*c(-1,1,-1,0,0,0,0.5)
      dy_ch <- 0.07*c(-1,0,-1,1,0.5,1,0.5)
    }
    col_ch <- as.vector(sapply(rownames(ch_all)[-nrow(ch_all)], function (taxnm) {
      as.character(switch(df[which(df$taxname == taxnm), "det_by"],
                          "eDNA" = color_edna, 
                          "trawl" = color_trawl,
                          "both" = color_both))
    }))
    text(x = ch_all[-nrow(ch_all),1] + dx_ch, y = ch_all[-nrow(ch_all), 2] + dy_ch, labels = short_species_name(rownames(ch_all)[-nrow(ch_all)]), cex = 0.6, adj = 0, col = col_ch, font = 3)
    points(x = Fct_space[, dm[1]], y = Fct_space[, dm[2]], col = Fct_space$color, 
           pch = Fct_space$pch, cex = Fct_space$cex)
    
    if (show.legend) {
    }
    
  }
  plot1(to_plot(c(1,2), show.legend = T))
  plot2(to_plot(c(1,3), show.xlab = F))
 
   par(xpd = NA)
  legend("topleft", legend = c("eDNA", "trawling", "both"), 
         pch = 16, col = c(color_edna, color_trawl, color_both), cex = 0.8, 
         bg = "white", inset = c(-0.2,-0.1))
  par(xpd = F)
  
  dev.off()
}

plot_figure1()
