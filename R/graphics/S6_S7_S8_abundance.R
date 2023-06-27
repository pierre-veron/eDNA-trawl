plot_abundance_fig <- function() {
  
  # Finding common taxa ----------------------------------------------------------
  EDNA_2019_RelatReads <- EDNA_2019_reads
  for (col in colnames(EDNA_2019_RelatReads)) {
    EDNA_2019_RelatReads[, col] <- EDNA_2019_RelatReads[, col] / sum(EDNA_2019_RelatReads[, col])
  }
  filters <- colnames(EDNA_2019_RelatReads)
  df <- EDNA_2019_RelatReads
  df$fish <- rownames(df)
  df <- as.data.frame(df %>% pivot_longer(cols = !fish , names_to = "filter", values_to = "relat_reads"))
  df <- df[which(df$relat_reads > 0),]
  
  
  df$station <- as.vector(sapply(df$filter, function(filt) {
    as.character(EDNA_metadata[which(EDNA_metadata$filter == filt), "station"])
  }))
  
  df$repl <- as.vector(sapply(df$filter, function(filt) {
    as.character(EDNA_metadata[which(EDNA_metadata$filter == filt), "replicate"])
  }))
  
  df$nb_trawl <- as.vector(sapply(1:nrow(df), function(i) {
    fish <- df[i, "fish"]
    stat <- df[i, "station"]
    n <- Trawl_2019_abundance[fish, stat]
    if (is.na(n)) {
      n <- 0
    }
    n
  }))
  
  df$reads <- as.vector(sapply(1:nrow(df), function(i) {
    fish <- df[i, "fish"]
    filt <- df[i, "filter"]
    n <- EDNA_2019_reads[fish, filt]
    n
  }))
  
  
  Common_taxa <- df[which(df$nb_trawl > 0),]
  Common_taxa$log_nb_trawl <- log(Common_taxa$nb_trawl)
  
  # Fig S7 - CDF of abundance ----------------------------------------------------
  cdf_nb <- ecdf(Common_taxa$nb_trawl)
  cdf_read <- ecdf(Common_taxa$reads)
  
  out_file <- "figures/Fig_S7_abundance_distribution"
  ext <- ".pdf"
  
  image_width <- 17
  image_height <- 9
  
  pdf(file = paste0(out_file, ext),
      width = image_width / 2.54,
      height = image_height / 2.54)
  
  par(mfrow =  c(1,2), mar = c(0.5,1,1,1), oma = c(2,2,0,0))
  
  plot(cdf_nb,   main = "Trawl", xlab = "",
       ylab = "", xaxt = "n" , yaxt = "n")
  axis(side = 1, at = 2e4*(0:6), line = 0, tck = -0.015, labels = F)
  axis(side = 1, at = 2e4*(0:6), line = -0.6, tick = F, lty = 0, 
       labels = 2*(0:6))
  
  mtext("Number of individuals (1e4)", side = 1, line = 1.5)
  mtext("CDF", side = 2, line = 2)
  
  axis(side = 2, at = 0.2*(0:5), line = 0, tck = -0.015, labels = F)
  axis(side = 2, at = 0.2*(0:5), line = -0.3, tick = F, lty = 0, labels = T, las = T)
  
  plot(cdf_read,  main = "eDNA", xlab = "",
       ylab = "", xaxt = "n" , yaxt = "n")
  
  axis(side = 1, at = 1e5*(0:5), line = 0, tck = -0.015, labels = F)
  axis(side = 1, at = 1e5*(0:5), line = -0.6, tick = F, lty = 0, 
       labels = 0:5)
  
  mtext("Number of reads (1e5)", side = 1, line = 1.5)
  axis(side = 2, at = 0.2*(0:5), line = 0, tck = -0.015, labels = F)
  
  dev.off()
  
  
  # Fig S8 - Correlation all filters ---------------------------------------------
  # building a dataframe that contains the average and difference of the relative
  # number of reads in each sampling site, between the two filters 
  # and the trawling catches 
  
  Common_taxa_merged_filt <- Common_taxa
  Common_taxa_merged_filt <- as.data.frame(Common_taxa_merged_filt %>% 
                                             group_by(fish, station) %>% 
                                             summarise_at(vars(relat_reads, nb_trawl, log_nb_trawl), list(group = mean)))
  
  colnames(Common_taxa_merged_filt) <- c("fish", "station", "av_relat_reads", "nb_trawl", "log_nb_trawl")
  
  Common_taxa_merged_filt$delta_relat_reads <- 
    as.vector(sapply(1:nrow(Common_taxa_merged_filt), function(i) {
      fish <- Common_taxa_merged_filt[i, "fish"]
      station <- Common_taxa_merged_filt[i, "station"]
      relat_reads <- Common_taxa[which(Common_taxa$fish == fish & Common_taxa$station == station), "relat_reads"]
      if (length(relat_reads) > 1) {
        abs(relat_reads[1] - relat_reads[2])
      } else {
        0
      }
    }))
  
  
  out_file <- "figures/Fig_S8_abundance"
  ext <- ".pdf"
  
  image_width <- 17
  image_height <- 15
  res = 400
  
  
  pdf(file = paste0(out_file, ext),
      width = image_width / 2.54,
      height = image_height / 2.54)
  
  par(mfrow = c(1,1), mar = c(2.5,2.8,1,1))
  plot(Common_taxa_merged_filt$log_nb_trawl, Common_taxa_merged_filt$av_relat_reads, pch = 20, 
       xlab = "",
       ylab = "", 
       ylim = c(0, max( Common_taxa_merged_filt$av_relat_reads + 0.5 * Common_taxa_merged_filt$delta_relat_reads)), 
       xaxt = "n" , yaxt = "n")
  Hmisc::errbar(x = Common_taxa_merged_filt$log_nb_trawl, y =  Common_taxa_merged_filt$av_relat_reads,
                yplus = Common_taxa_merged_filt$av_relat_reads + 0.5 * Common_taxa_merged_filt$delta_relat_reads,
                yminus =  Common_taxa_merged_filt$av_relat_reads - 0.5 * Common_taxa_merged_filt$delta_relat_reads, 
                add = T,cap = 0)
  lreg <- lm(Common_taxa_merged_filt$av_relat_reads ~ Common_taxa_merged_filt$log_nb_trawl)
  print(summary(lreg))
  abline(lreg, col = "blue", lt = 2)
  
  
  axis(side = 1, at = 2*(0:6), line = 0, tck = -0.015, labels = F)
  axis(side = 1, at = 2*(0:6), line = -0.6, tick = F, lty = 0, labels = T)
  
  axis(side = 2, at = 0.1*(0:4), line = 0, tck = -0.015, labels = F)
  axis(side = 2, at = 0.1*(0:4), line = -0.3, tick = F, lty = 0, labels = T, las = T)
  mtext("Log nb individuals trawl", side = 1, line = 1.4)
  mtext("Relative abundance eDNA", side = 2, line = 1.9)
  dev.off()
  
  
  
  # Fig S9 - Correlation per filter ----------------------------------------------
  out_file <- "figures/Fig_S9_abundance_per_filter"
  ext <- ".pdf"
  
  image_width <- 17
  image_height <- 20
  res = 400
  
  color_repl <- c("#AD7A99", "#3E8914")
  
  pdf(file = paste0(out_file, ext),
      width = image_width / 2.54,
      height = image_height / 2.54)
  
  par(mfrow = c(5,3))
  par(oma = c(4, 4, 0, 0))
  for (i in 1:length(Name_stations)) {
    Station_fish <- Common_taxa[which(Common_taxa$station == names(Name_stations)[i]),]
    
    if (i %in% c(1, 4, 7, 10, 13)) {
      mar_dy <- 0
    } else {
      mar_dy <- 0
    }
    if (i %in% c(13, 14, 15)) {
      mar_dx <- 0
    } else {
      mar_dx <- 0
    }
    if (i == 7) {
      ylab <- "Filter-relative number of reads"
    } else {
      ylab <- ""
    }
    if (i == 14) {
      xlab <- "Log nb of individuals trawl"
    } else {
      xlab <- ""
    }
    par(mar = c(1+mar_dx,1+mar_dy,1,1))
    plot(c(), c(), xlim = c(0,12),
         ylim = c(0,0.5),
         xlab = xlab, ylab = ylab,
         cex.lab = 1.2,
         main= "",
         xaxt = "n",
         yaxt = "n")
    for (irep in c(1,2)) {
      Filter_fish <- Station_fish[which(Station_fish$repl == irep), ]
      lreg <- lm(Filter_fish$relat_reads ~ Filter_fish$log_nb_trawl)
      
      lreg_details <- c("intercept" = as.numeric(lreg$coefficients[1]), 
                        "slope" = as.numeric(lreg$coefficients[2]),
                        "r.squared" = summary(lreg)[["r.squared"]], 
                        "p.intercept" = summary(lreg)$coefficients[1,4],
                        "p.slope" = summary(lreg)$coefficients[2,4], 
                        "n" = nrow(Filter_fish))
      
      points(Filter_fish$log_nb_trawl, Filter_fish$relat_reads, pch = 19,col = color_repl[irep])
      abline(lreg, col = color_repl[irep], lty = 2)
      
      if (lreg_details["p.slope"] < 0.001) {
        lev_signif <- "***"
      } else if (lreg_details["p.slope"] < 0.01) {
        lev_signif <- "**"
      } else if (lreg_details["p.slope"] < 0.05){
        lev_signif <- "*"
      } else {
        lev_signif <- "ns"
      }
      lreg_label <- paste0("r² = ", round(lreg_details["r.squared"], 2),
                           ", ", lev_signif)
      if (i %in% c(7, 9, 15)) {
        text(x = 8.00, y = 0.42 - irep * 0.06, labels = lreg_label, adj = 0, cex = 0.9, col = color_repl[irep])
      } else {
        text(x = 0.00, y = 0.42 - irep * 0.06, labels = lreg_label, adj = 0, cex = 0.9, col = color_repl[irep])
      }
      
    }
    
    text(x = 0, y = 0.45, labels = paste0("Site ",i ), cex = 1.3, adj = 0)
    if (i %in% c(1, 4, 7, 10, 13)) {
      axis(2, at = 0.1*(0:5),
           las = 2, tck  =F, line = -0.3,lty = 0)
    } 
    axis(2,at =  0.1*(0:5), 
         labels = F, tck  =-0.03)
    if (i %in% c(13, 14, 15)) {
      axis(1, at =  2*(0:6),
           tck  =F, line = -0.3, lty = 0)
    }
    axis(1, at =  2*(0:6),
         labels = F, tck  =-0.03)
  }
  
  mtext("Log nb of individuals trawl", side = 1, outer = TRUE, line = 2)
  mtext("Filter-relative number of reads", side = 2, outer = TRUE, line = 2)
  
  dev.off()
  
  
  # Without Trachurus
  Common_taxa = Common_taxa[which(Common_taxa$fish != "Trachurus"),]
  # Fig S7 - CDF of abundance ----------------------------------------------------
  cdf_nb <- ecdf(Common_taxa$nb_trawl)
  cdf_read <- ecdf(Common_taxa$reads)
  
  out_file <- "figures/without_Trachurus_Fig_S7_abundance_distribution"
  ext <- ".pdf"
  
  image_width <- 17
  image_height <- 9
  
  pdf(file = paste0(out_file, ext),
      width = image_width / 2.54,
      height = image_height / 2.54)
  
  par(mfrow =  c(1,2), mar = c(0.5,1,1,1), oma = c(2,2,0,0))
  
  plot(cdf_nb,   main = "Trawl", xlab = "",
       ylab = "", xaxt = "n" , yaxt = "n")
  axis(side = 1, at = 2e4*(0:6), line = 0, tck = -0.015, labels = F)
  axis(side = 1, at = 2e4*(0:6), line = -0.6, tick = F, lty = 0, 
       labels = 2*(0:6))
  
  mtext("Number of individuals (1e4)", side = 1, line = 1.5)
  mtext("CDF", side = 2, line = 2)
  
  axis(side = 2, at = 0.2*(0:5), line = 0, tck = -0.015, labels = F)
  axis(side = 2, at = 0.2*(0:5), line = -0.3, tick = F, lty = 0, labels = T, las = T)
  
  plot(cdf_read,  main = "eDNA", xlab = "",
       ylab = "", xaxt = "n" , yaxt = "n")
  
  axis(side = 1, at = 1e5*(0:5), line = 0, tck = -0.015, labels = F)
  axis(side = 1, at = 1e5*(0:5), line = -0.6, tick = F, lty = 0, 
       labels = 0:5)
  
  mtext("Number of reads (1e5)", side = 1, line = 1.5)
  axis(side = 2, at = 0.2*(0:5), line = 0, tck = -0.015, labels = F)
  
  dev.off()
  
  
  # Fig S8 - Correlation all filters ---------------------------------------------
  # building a dataframe that contains the average and difference of the relative
  # number of reads in each sampling site, between the two filters 
  # and the trawling catches 
  
  Common_taxa_merged_filt <- Common_taxa
  Common_taxa_merged_filt <- as.data.frame(Common_taxa_merged_filt %>% 
                                             group_by(fish, station) %>% 
                                             summarise_at(vars(relat_reads, nb_trawl, log_nb_trawl), list(group = mean)))
  
  colnames(Common_taxa_merged_filt) <- c("fish", "station", "av_relat_reads", "nb_trawl", "log_nb_trawl")
  
  Common_taxa_merged_filt$delta_relat_reads <- 
    as.vector(sapply(1:nrow(Common_taxa_merged_filt), function(i) {
      fish <- Common_taxa_merged_filt[i, "fish"]
      station <- Common_taxa_merged_filt[i, "station"]
      relat_reads <- Common_taxa[which(Common_taxa$fish == fish & Common_taxa$station == station), "relat_reads"]
      if (length(relat_reads) > 1) {
        abs(relat_reads[1] - relat_reads[2])
      } else {
        0
      }
    }))
  
  
  out_file <- "figures/without_Trachurus_Fig_S8_abundance"
  ext <- ".pdf"
  
  image_width <- 17
  image_height <- 15
  res = 400
  
  
  pdf(file = paste0(out_file, ext),
      width = image_width / 2.54,
      height = image_height / 2.54)
  
  par(mfrow = c(1,1), mar = c(2.5,2.8,1,1))
  plot(Common_taxa_merged_filt$log_nb_trawl, Common_taxa_merged_filt$av_relat_reads, pch = 20, 
       xlab = "",
       ylab = "", 
       ylim = c(0, max( Common_taxa_merged_filt$av_relat_reads + 0.5 * Common_taxa_merged_filt$delta_relat_reads)), 
       xaxt = "n" , yaxt = "n")
  Hmisc::errbar(x = Common_taxa_merged_filt$log_nb_trawl, y =  Common_taxa_merged_filt$av_relat_reads,
                yplus = Common_taxa_merged_filt$av_relat_reads + 0.5 * Common_taxa_merged_filt$delta_relat_reads,
                yminus =  Common_taxa_merged_filt$av_relat_reads - 0.5 * Common_taxa_merged_filt$delta_relat_reads, 
                add = T,cap = 0)
  lreg <- lm(Common_taxa_merged_filt$av_relat_reads ~ Common_taxa_merged_filt$log_nb_trawl)
  print("without Trachurus")
  print(summary(lreg))
  abline(lreg, col = "blue", lt = 2)
  
  
  axis(side = 1, at = 2*(0:6), line = 0, tck = -0.015, labels = F)
  axis(side = 1, at = 2*(0:6), line = -0.6, tick = F, lty = 0, labels = T)
  
  axis(side = 2, at = 0.1*(0:4), line = 0, tck = -0.015, labels = F)
  axis(side = 2, at = 0.1*(0:4), line = -0.3, tick = F, lty = 0, labels = T, las = T)
  mtext("Log nb individuals trawl", side = 1, line = 1.4)
  mtext("Relative abundance eDNA", side = 2, line = 1.9)
  dev.off()
  
  
  
  # Fig S9 - Correlation per filter ----------------------------------------------
  out_file <- "figures/without_Trachurus_Fig_S9_abundance_per_filter"
  ext <- ".pdf"
  
  image_width <- 17
  image_height <- 20
  res = 400
  
  color_repl <- c("#AD7A99", "#3E8914")
  
  pdf(file = paste0(out_file, ext),
      width = image_width / 2.54,
      height = image_height / 2.54)
  
  par(mfrow = c(5,3))
  par(oma = c(4, 4, 0, 0))
  for (i in 1:length(Name_stations)) {
    Station_fish <- Common_taxa[which(Common_taxa$station == names(Name_stations)[i]),]
    
    if (i %in% c(1, 4, 7, 10, 13)) {
      mar_dy <- 0
    } else {
      mar_dy <- 0
    }
    if (i %in% c(13, 14, 15)) {
      mar_dx <- 0
    } else {
      mar_dx <- 0
    }
    if (i == 7) {
      ylab <- "Filter-relative number of reads"
    } else {
      ylab <- ""
    }
    if (i == 14) {
      xlab <- "Log nb of individuals trawl"
    } else {
      xlab <- ""
    }
    par(mar = c(1+mar_dx,1+mar_dy,1,1))
    plot(c(), c(), xlim = c(0,12),
         ylim = c(0,0.5),
         xlab = xlab, ylab = ylab,
         cex.lab = 1.2,
         main= "",
         xaxt = "n",
         yaxt = "n")
    for (irep in c(1,2)) {
      Filter_fish <- Station_fish[which(Station_fish$repl == irep), ]
      lreg <- lm(Filter_fish$relat_reads ~ Filter_fish$log_nb_trawl)
      
      lreg_details <- c("intercept" = as.numeric(lreg$coefficients[1]), 
                        "slope" = as.numeric(lreg$coefficients[2]),
                        "r.squared" = summary(lreg)[["r.squared"]], 
                        "p.intercept" = summary(lreg)$coefficients[1,4],
                        "p.slope" = summary(lreg)$coefficients[2,4], 
                        "n" = nrow(Filter_fish))
      
      points(Filter_fish$log_nb_trawl, Filter_fish$relat_reads, pch = 19,col = color_repl[irep])
      abline(lreg, col = color_repl[irep], lty = 2)
      
      if (lreg_details["p.slope"] < 0.001) {
        lev_signif <- "***"
      } else if (lreg_details["p.slope"] < 0.01) {
        lev_signif <- "**"
      } else if (lreg_details["p.slope"] < 0.05){
        lev_signif <- "*"
      } else {
        lev_signif <- "ns"
      }
      lreg_label <- paste0("r² = ", round(lreg_details["r.squared"], 2),
                           ", ", lev_signif)
      if (i %in% c(7, 9, 15)) {
        text(x = 8.00, y = 0.42 - irep * 0.06, labels = lreg_label, adj = 0, cex = 0.9, col = color_repl[irep])
      } else {
        text(x = 0.00, y = 0.42 - irep * 0.06, labels = lreg_label, adj = 0, cex = 0.9, col = color_repl[irep])
      }
      
    }
    
    text(x = 0, y = 0.45, labels = paste0("Site ",i ), cex = 1.3, adj = 0)
    if (i %in% c(1, 4, 7, 10, 13)) {
      axis(2, at = 0.1*(0:5),
           las = 2, tck  =F, line = -0.3,lty = 0)
    } 
    axis(2,at =  0.1*(0:5), 
         labels = F, tck  =-0.03)
    if (i %in% c(13, 14, 15)) {
      axis(1, at =  2*(0:6),
           tck  =F, line = -0.3, lty = 0)
    }
    axis(1, at =  2*(0:6),
         labels = F, tck  =-0.03)
  }
  
  mtext("Log nb of individuals trawl", side = 1, outer = TRUE, line = 2)
  mtext("Filter-relative number of reads", side = 2, outer = TRUE, line = 2)
  
  dev.off()
}

plot_abundance_fig()