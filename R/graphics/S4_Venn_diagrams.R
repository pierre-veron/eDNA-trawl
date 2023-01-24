### ----------------------------------------------------------------------------
#   This script plots a Venn Diagram by sampling site between eDNA and trawling.
#
#   Author : Pierre Veron, pierre.veron.2017@polytechnique.org
#   Created 2022/04/20
### ----------------------------------------------------------------------------

plot_venn_diag <- function() {
  out_file <- "figures/Fig_S4_Venn_diagrams"
  
  ext <- ".pdf"
  
  image_width <- 17
  image_height <- 20
  res = 800
  
  
  pdf(file = paste0(out_file, ext), 
      width = image_width / 2.54,
      height = image_height / 2.54)
  
  par(mfrow = c(5,3))
  par(oma = c(0, 0, 0, 0))
  for (i_site in 1:length(Name_stations)) {
    station <- names(Name_stations)[i_site]
    Set_edna <- colnames(Incidence_matrix)[which(Incidence_matrix[paste0(station, "e"),] > 0)]
    Set_trawl <- colnames(Incidence_matrix)[which(Incidence_matrix[paste0(station, "t"),] > 0)]
    par(mar = c(0,0,0,0))
    plot(c(), c(), xlim = c(-3,3), ylim = c(-3,3), xaxt = "n", yaxt = "n", axes = F)
    vd <-my_venn_diag(list(Set_edna, Set_trawl), colors = c(color_edna, color_trawl),
                      alpha = 0.5, size = 0.25, angle = pi/8, text.offset =1.4)
    text(x = -2.5, y = 2.5, label = Name_stations[station], cex = 1.7)
    if (i_site == 1) {
      text(-2,-2, "eDNA", col = color_edna, cex = 1.3)
      text(2,-1,"trawling", col = color_trawl, cex = 1.3)
    }
  }
  dev.off()
}

plot_venn_diag()
