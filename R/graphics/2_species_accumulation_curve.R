### ----------------------------------------------------------------------------
#   Draws a species accumulation curve for eDNA and for trawling and compares them
### ----------------------------------------------------------------------------

plot_SAC <- function() {
  out_file <- "figures/Fig_2_species_accumulation_curve"
  ext <- ".pdf"
  image_width <- 17
  image_height <- 15
  res = 800
  
  
  Data_Atl <- Incidence_matrix
  
  # Calculates the spec accum
  data_e <-  Data_Atl[grep("e",rownames(Data_Atl)),]
  data_t <- Data_Atl[grep("t",rownames(Data_Atl)),]
  data_m <- Incidence_matrix_merged
  
  speac_e <- specaccum(comm= data_e )
  speac_t <- specaccum(comm= data_t )
  speac_m <- specaccum(comm = data_m)
  
  # interpolation
  sr_edna <- speac_e$richness[1] #SR for 1 site 
  sr_trawl <- speac_t$richness
  site_encadr <- speac_t$sites[c(max(which(sr_trawl < sr_edna)), min(which(sr_trawl > sr_edna)))]
  rich_encadr <- sr_trawl[c(max(which(sr_trawl < sr_edna)), min(which(sr_trawl > sr_edna)))]
  
  site_for_equal_sr <- site_encadr[1] +  (sr_edna - rich_encadr[1])*(site_encadr[2] - site_encadr[1]) / (rich_encadr[2] - rich_encadr[1])
  
  # Fit 
  models <- c("asymp" , "gompertz", "michaelis-menten", "logis")
  fit.speac <- lapply(list("eDNA" = speac_e, 
                           "trawling" = speac_t, 
                           "merged" = speac_m), 
                      function(speac) {
                        fitmodels <- lapply(models, function(modl) {
                          fitspecaccum(speac, modl)})
                        names(fitmodels) <- models
                        fitmodels
                      }
  )
  
  # Plot
  
  pdf(file = paste0(out_file, ext), 
      width = image_width / 2.54,
      height = image_height / 2.54)
  
  
  par(mar = c(3,3,0.2,0.2))
  
  plot(speac_m, col = "black",xlab = "",  ylab = "",lwd = 1.7, xaxt = "n", yaxt = "n", xlim = c(1,15), lty = 0)
  points(x = speac_m$sites, y = speac_m$richness, pch = 16, col = "black")
  plot(fit.speac$merged$asymp, lty = 2, col = "black", add = T)
  val_asym_m <- coef(fit.speac$merged$asymp)["Asym"]
  
  plot(speac_e, add = T, col = color_edna, lwd = 1.7, lty = 0)
  points(x = speac_e$sites, y = speac_e$richness, pch = 16, col = color_edna)
  plot(fit.speac$eDNA$asymp, lty = 2, col =color_edna, add = T)
  val_asym_e <- coef(fit.speac$eDNA$asymp)["Asym"]
  
  plot(speac_t,add=T,col=color_trawl, lwd = 1.7, lty = 0)
  points(x = speac_t$sites, y = speac_t$richness, pch = 16, col = color_trawl)
  plot(fit.speac$trawling$asymp, lty = 2, col =color_trawl, add = T)
  val_asym_t <- coef(fit.speac$trawling$asymp)["Asym"]
  
  
  text(x = 6, y = 30, label = paste0("trawling\nmax = ", round(val_asym_t, 1)), col = color_trawl)
  text(x = 12, y = 80, label =paste0("eDNA\nmax = ", round(val_asym_e, 1)), col = color_edna)
  text(x = 4, y = 95, label = paste0("combined eDNA and trawling\nmax = ", round(val_asym_m, 1)), col = "black")
  
  axis(at = 1:15,side=1,tcl=-0.4,bg="white", labels=F, line = 0)
  axis(side=1,tick = F,bg="white",labels=T, line = -0.5)
  
  axis(side=2,tcl=-0.4,bg="white", labels=F, line = 0)
  axis(side=2,tick = F,bg="white",labels=T, line = -0.35, las = 1)
  
  mtext("Number of sites", side = 1, line = 1.8, cex = 1.3)
  mtext("Taxa richness", side = 2, line = 1.9, cex = 1.3)
  dev.off()
  fit.speac
}

fit.speac <- plot_SAC()