### ----------------------------------------------------------------------------
#   Creates a map with richness indices of diversity with + and - for SES
#    
#
#   Author : Pierre Veron, pierre.veron.2017@polytechnique.org
#   Created 2022/04
### ----------------------------------------------------------------------------

plot_figure <- function() {
  # Graphical settings -----------------------------------------------------------
  circle_factor <- 0.32 # size of the pies for median indices 
  
  # PART 1: Alpha indices --------------------------------------------------------
  imgSignif <- list("+" = readPNG(here("data", "icons", "char", "plus.png")),
                    "-" = readPNG(here("data", "icons", "char", "minus.png")))
  
  
  
  
  ## Define indices to plot ------------------------------------------------------
  indices_rich <- c('taxo' = 'SR.mean', # which indices to plot
                    'func' = 'FD.mean',
                    'phyl' = 'PD.mean')
  
  indic_names <-c('taxo' = 'SR', # which indices to plot
                  'func' = 'FD',
                  'phyl' = 'PD')
  
  ## Initiate plot ---------------------------------------------------------------
  image_width <- 10
  image_height <- 10
  res = 800
  
  pdf(file = "figures/Fig_4_alpha_div_merged.pdf",
      width = image_width / 2.54,
      height = image_height / 2.54)
  par(mfrow = c(1,1))
  
  
  ## Draw the 3 maps alpha -------------------------------------------------------
  
  par(mar = c(1, 1.6, 0.15, 0.15))
  
  plot(0,0, type='n', xaxt='n', yaxt='n', xlab='', ylab='', bty='n')
  
  par(new = "T")
  
  
  Points <- Alpha_merged
  
  
  Points$div_t <- Points[, indices_rich['taxo']]
  scale_t <- circle_factor / median(Points$div_t)
  
  
  #functio indice
  Points$div_f <- Points[, indices_rich["func"]]
  scale_f <- circle_factor / median(Points$div_f)
  
  Points$signif_f <- as.vector(sapply(Points$station, function(site) {
    ses <- RESULTS_SES_fct_merged[site, "FD_q0"]
    if (ses >  2.58) {
      "+"
    } else if (ses > 1.96) {
      "+"
    } else if (ses < -2.58) {
      "-"
    } else if (ses < -1.96) {
      "-"
    } else {
      ""
    }
  }))
  
  #phylo indice
  Points$div_p <- Points[, indices_rich["phyl"]]
  scale_p <- circle_factor / median(Points$div_p)
  
  Points$signif_p <- as.vector(sapply(Points$station, function(site) {
    ses <- RESULTS_SES_phyl_merged[site, "PD.mean"]
    if (ses >  2.58) {
      "+"
    } else if (ses > 1.96) {
      "+"
    } else if (ses < -2.58) {
      "-"
    } else if (ses < -1.96) {
      "-"
    } else {
      ""
    }
  }))
  
  
  pre_map(GISdata = GISdata, xlim = c(-5,0), ylim = c(43,47))
  # plot points 
  
  
  
  eps <- 0 # space between pie slices
  dx <- 0.05
  
  shift <- list("X0491" = c(-1,46.4), 
                "X0481" = c(-0.75,43.5),
                "X0477" = c(-2.5,43.5),
                "X0450" = c(-4.2,45.5),
                "X0496" = c(-1.5,47.2)) # where to plot a point (if not set, the point
  # will be plotted at its real position)
  
  for (i in 1:nrow(Points)) {
    Pt <- Points[i, ]
    if (!is.null(shift[[Pt$station]])) {
      x <- shift[[Pt$station]][1]
      y <- shift[[Pt$station]][2]
      lines(x = c(x, Pt$Longitude), y = c(y, Pt$Latitude),
            col = "black",
            lwd = 1.5)
      
    }
  }
  
  for (i in 1:nrow(Points)) {
    
    Pt <- Points[i, ]
    angles <- pi * c(2/3+eps, 1-eps, 1/3+eps, 2/3-eps, eps, 1/3 - eps)
    color <- "#635666"
    dy <- 1
    
    
    # taxo
    if (!is.null(shift[[Pt$station]])) {
      x <- shift[[Pt$station]][1]
      y <- shift[[Pt$station]][2]
    } else {
      x <- Pt$Longitude
      y <- Pt$Latitude
    }
    points(x = Pt$Longitude, y = Pt$Latitude, pch = 16, col = "black")
    
    pie_slice(x = x - dx, 
              y = y + dy* 0.025, 
              r = scale_t * Pt$div_t,
              col = color,
              angle = angles[c(1,2)])
    
    pie_slice(x = x, 
              y = y + 0.5*dx*dy + dy*0.025, 
              r = scale_p * Pt$div_p,
              col = color,
              angle = angles[c(3,4)])
    
    if (Pt$signif_p %in% c("+", "-")) {
      addImg(imgSignif[[Pt$signif_p]], 
             x = x + 1.2* 0.245 * cos(mean(angles[c(3,4)])),
             y = y + dy* 0.025  + 0.85 *0.245 * sin(mean(angles[c(3,4)])),
             width = 0.15)
    }
    pie_slice(x = x + dx, 
              y = y + dy *0.025, 
              r = scale_f * Pt$div_f,
              col = color,
              angle = angles[c(5,6)])
    
    if (Pt$signif_f %in% c("+", "-")) {
      addImg(imgSignif[[Pt$signif_f]], 
             x = x + 1.1* 0.245 * cos(mean(angles[c(5,6)])),
             y = y + dy* 0.025  + 0.7 *0.245 * sin(mean(angles[c(5,6)])),
             width = 0.15)
    }
    
    
    
  }
  
  
  ## Legend --------------------------------------------------------------------
  x <- -3.6
  y <- 43.9
  rect(x - 1.49,y-1.5,x+0.6, y + 0.85, col = "white")
  
  
  dy <- 1
  pie_slice(x = x - 0.445 - dx, 
            y = y + dy* 0.025 - 0.15, 
            r = 1.7*scale_t * median(Points$div_t),
            col = color,
            angle = angles[c(1,2)])
  
  
  
  pie_slice(x = x - 0.445, 
            y = y + 0.5*dx+ dy*0.025 - 0.15, 
            r = 1.7*scale_p * median(Points$div_p),
            col = color,
            angle = angles[c(3,4)])
  
  pie_slice(x = x -0.445 + dx, 
            y = y  + dy *0.025 - 0.15, 
            r = 1.7*scale_f * median(Points$div_f),
            col = color,
            angle = angles[c(5,6)])
  
  
  text(x - 0.445,
       y = y + 0.63,
       label = "Richness",
       cex = 0.95,
       adj = 0.5,
       col = "black",
       font = 2)
  
  
  r <- 0.4 
  angles <- c(5*pi/6)
  text(x = x + r * cos(angles) - 0.445,
       y = y + r * sin(angles) - 0.15,
       label = indic_names["taxo"],
       adj = 0.5,
       col = "white",
       cex = 0.65)
  
  angles <- c(pi/2)
  text(x = x + r * cos(angles) - 0.445,
       y = y + r * sin(angles) - 0.15,
       label = indic_names["phyl"],
       adj = 0.5,
       col = "white",
       cex = 0.65)
  angles <- c(pi/6)
  text(x = x + r * cos(angles) - 0.445,
       y = y + r * sin(angles) - 0.15,
       label = indic_names["func"],
       adj = 0.5,
       col = "white",
       cex = 0.65)
  
  text(x - 1.4,
       y - 0.75,
       label = "For phylo and
functio:
+ : overdispersed
- : clustered",
       cex = 0.8,
       adj = 0)
  
  post_map(xtickslab = T, ytickslab = T, cex.axis = 0.85)
  dev.off()
}

plot_figure()