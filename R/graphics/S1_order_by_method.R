### ----------------------------------------------------------------------------
#   This script draws the distributions of taxonomic groups by detection method
#
#   Author : Pierre Veron
#   pierre.veron.2017@polytechnique.org
#   Created 2022/04/25
### ----------------------------------------------------------------------------

plot_distribution_clades <- function() {
  rank <- "order"
  df <- Reftax
  
  # General graphical parameters -------------------------------------------------
  out_file <- "figures/Fig_S1_order_by_method"
  ext <- ".pdf"
  
  image_width <- 20
  image_height <- 14
  res = 800
  
  pdf(file = paste0(out_file, ext),
      width = image_width / 2.54,
      height = image_height / 2.54)
  
  df <- df[order(df[,'class']), ]
  
  # counting the number of sites where each taxa is detected
  df$site_edna <- as.vector(sapply(df$taxname, function(tax) {
    sum(Incidence_matrix_eDNA[, gsub(" ", ".", tax)])
  }))
  
  df$site_trawl <- as.vector(sapply(df$taxname, function(tax) {
    sum(Incidence_matrix_trawl[, gsub(" ", ".", tax)])
  }))
  
  # group by order (rank = "order")
  clades <- unique(df[, rank])
  
  freq_edna <- as.vector(sapply(clades, function(cl) {
    sum(df[which(df[, rank] == cl), "site_edna"])
  }))
  freq_edna <- freq_edna / sum(freq_edna)
  names(freq_edna) <- clades
  
  freq_trawl <- as.vector(sapply(clades, function(cl) {
    sum(df[which(df[, rank] == cl), "site_trawl"])
  }))
  freq_trawl <- freq_trawl / sum(freq_trawl)
  names(freq_trawl) <- clades
  
  
  # Sort by freq_trawl
  sorting <- order(-freq_edna)
  freq_trawl <- freq_trawl[sorting]
  freq_edna <- freq_edna[sorting]
  clades <- clades[sorting]
  
  if (rank == "order") {
    short_clades <- gsub("formes", ".", clades)
  } else {
    short_clades <- clades
  }
  
  # Plot -------------------------------------------------------------------------
  par(mar = c(5,3,0.5,0.5))
  plot(c(), c(), xlim = c(0.5,length(clades)+0.5),
       ylim = max(freq_trawl, freq_edna) * c(0,1.1), 
       xaxt = 'n',
       yaxt = 'n',
       xlab = "",
       ylab = "")
  
  axis(side = 1, at = 1:length(clades), line = 0, tck = -0.01, labels = F)
  axis(side = 1, at = 1:length(clades), line = -0.6, tick = F, lty = 0, 
       labels = short_clades, cex.axis = 0.8, las = 2)
  
  axis(side = 2, at = seq(-0.16,0.16, 0.04), line = 0, tck = -0.01, labels = F)
  axis(side = 2, at = seq(-0.16,0.16, 0.04), line = -0.6, tick = F, lty = 0, 
       labels = T, las = 2, cex.axis = 0.8)
  
  mtext("Frequency", side = 2, line = 2)
  
  
  # grid
  for (a in  seq(-0.16,0.16, 0.02)) {
    abline(a = a, b = 0, col = "lightgray", lty = 1, lwd  = 1)
  }

  
  rect(xleft=1:length(clades) - 0.4, ybottom=0, xright=1:length(clades) , 
       ytop = freq_edna, col = color_edna, border = NA)
  rect(xleft=1:length(clades) , ybottom=0, xright=1:length(clades) + 0.4, 
       ytop = freq_trawl, col = color_trawl, border = NA)
  
  
  
  # Chi square test --------------------------------------------------------------
  
  freq_edna <- as.vector(sapply(clades, function(cl) {
    sum(df[which(df[, rank] == cl), "site_edna"])
  }))
  names(freq_edna) <- clades
  
  freq_trawl <- as.vector(sapply(clades, function(cl) {
    sum(df[which(df[, rank] == cl), "site_trawl"])
  }))
  names(freq_trawl) <- clades
  
  distrib <- t(data.frame( "trawling" = freq_trawl, "eDNA" = freq_edna))
  print("Chi-square test on distribution among clades between eDNA and trawling:")
  chitest <- chisq.test(distrib)
  
  text(x = 0.5, y = 0.17, 
       labels = bquote("Pearson's"~chi^2~"test with"~.(chitest$parameter)~"degrees of freedom:"~chi^2~"="~.(signif(chitest$statistic, 2))~", p = "~.(signif(chitest$p.value,2))),
       adj = 0)
  legend("topright", inset = 0.02, fill = c(color_edna, color_trawl), legend = c("eDNA", "trawling"), bg = "white")
  dev.off()
}

plot_distribution_clades()