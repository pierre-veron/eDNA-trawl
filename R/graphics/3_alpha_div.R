### ----------------------------------------------------------------------------
#   Creates three different maps (richness, regularity and divergence) with 
#    diversity indices plotted with pie slices, and ++/-- for SES significance. 
#
#   Author : Pierre Veron, pierre.veron.2017@polytechnique.org
#   Created 2022/04
### ----------------------------------------------------------------------------

# Figure overview :
#   
# Richness         |  Regularity       |  Divergence
# i_comp = 1       |  i_comp = 2       |  i_comp = 3
# i_panel= 1       |  i_panel = 2      |  i_panel = 3
#                  |                   |




plot_figure <- function() {
  
  # load the SES. of the diversity metrics.
  Alpha = read.delim("diversity_indices/alpha.csv", sep = ";", h=T)
  
  # load the SES. of the diversity metrics.
  Alpha_SES = read.delim("diversity_indices/alpha_ses.csv", sep = ";", h=T)
  
  Alpha$eSR.mean = as.numeric(gsub(",", ".", Alpha$eSR.mean, fixed = T))
  Alpha$tSR.mean = as.numeric(gsub(",", ".", Alpha$tSR.mean, fixed = T))
  
  mean.SR.eDNA = mean(Alpha$eSR.mean)
  mean.SR.trawl = mean(Alpha$tSR.mean)
  
  # remove the comma and convert to numeric
  Alpha_SES.mat = as.data.frame(do.call(cbind, lapply(1:dim(Alpha_SES[1:13])[2], function(x){
    as.numeric(gsub(",", ".", Alpha_SES[,x], fixed = T))
  })))
  names(Alpha_SES.mat) = names(Alpha_SES)[1:13]
  Alpha_SES.mat = data.frame(X = Alpha_SES.mat[,1], eSR = Alpha$eSR.mean, Alpha_SES.mat[,2:7], tSR = Alpha$tSR.mean, Alpha_SES.mat[,8:13])
  
  mean.ses = apply(Alpha_SES.mat, 2, mean)
  sd.ses = apply(Alpha_SES.mat, 2, sd)
  
  col.e = c(2:8)
  col.t = c(9:15)
  metric = gsub("^e", "", names(Alpha_SES.mat)[col.e], perl = T)
  
  ## Non paramteric Wilcoxon test with paired test.
  Wilcox.DF.paired = as.data.frame(do.call(rbind, lapply(1:length(col.e), function(x){
    a = wilcox.test(Alpha_SES.mat[,col.e[x]], Alpha_SES.mat[,col.t[x]], paired = T)
    c(a$statistic, round(a$p.value, digits = 4))
  })))
  names(Wilcox.DF.paired) = c("Statistic", "p.value")
  
  col.e1 = c(2:8)
  col.t1 = c(9:15)
  Wilcox.DF.paired = round(data.frame(Wilcox.DF.paired, Mean.ses.eDNA = mean.ses[col.e1], Sd.ses.eDNA = sd.ses[col.e1], Mean.ses.Trawl =  mean.ses[col.t1], Sd.ses.Trawl =  sd.ses[col.t1]), digits = 4)
  Wilcox.DF.paired = data.frame(Metric = metric, Wilcox.DF.paired)
  
  Wilcox.DF.paired
  
  A = Wilcox.DF.paired[,c(1,4,5)]
  B = Wilcox.DF.paired[,c(1,6,7)]
  names(A) = c("Metric", "Mean", "Sd")
  names(B) = c("Metric", "Mean", "Sd")
  
  
  DF1 = rbind(A, B)
  DF1 = data.frame(Method = c(rep("eDNA", 7), rep("Trawling", 7)),  DF1)
  DF1$Se = DF1$Sd/sqrt(15)
  # Remove the SR which is not a SES.
  DF1.SES = DF1[-c(1,8),]
  
  DF1.SES$Metric = rep(c("1_PD.ses", "5_VPD.ses", "3_MPD.ses", "2_FD.ses", "6_FEve.ses", "4_FDiv.ses"),2)
  DF1.SES$Metricnb = rep(c("1", "5", "3", "2", "6", "4"),2)
  DF1.SES <- DF1.SES[order(DF1.SES$Metricnb),]

  # Graphical settings -----------------------------------------------------------
  circle_factor <- 0.25 # size of the pies for median indices 
  
  # PART 1: Alpha indices --------------------------------------------------------
  imgSignif <- list("+" = readPNG(here("data", "icons", "char", "plus.png")),
                    "-" = readPNG(here("data", "icons", "char", "minus.png")))
  
  ## Define points ---------------------------------------------------------------
  Points_edna <- EDNA_metadata[which(EDNA_metadata$replicate == 1), 
                               c("station", "Latitude", "Longitude")]
  Points_trawl <- Points_edna
  
  Points_edna$name <- Name_stations[Points_edna$station] 
  Points_trawl$name <- Name_stations[Points_trawl$station]
  Points_edna <- as.data.frame(Points_edna)
  colnames(Points_trawl) <- c("station", "Latitude", "Longitude", "name")
  
  rownames(Points_trawl) <- Points_trawl$station
  
  Points_trawl[Points_edna$station, c("Latitude", "Longitude")] <- Points_edna[, c("Latitude", "Longitude")]
  rownames(Points_edna) <- paste0(Points_edna$station, "e")
  Points_edna$type <- "eDNA"
  rownames(Points_trawl) <- paste0(Points_trawl$station, "t")
  Points_trawl$type <- "trawl"
  Points <- rbind(Points_edna, Points_trawl)
  
  ## Define indices to plot ------------------------------------------------------
  indices_rich <- c('taxo' = 'sp_richness', # which indices to plot
                    'func' = 'FD_q0','phyl' = 'PD.mean')
  
  indices_eve <- c('taxo' = NULL, # which indices to plot
                   'func' = 'Mean.feve', 'phyl' = 'VPD.mean')
  
  indices_div <- c('taxo' = NULL,
                   'func' = 'Mean.fdiv',
                   'phyl' = 'MPD.mean')
  
  indices_list <- list("Richness" = indices_rich,"Divergence" = indices_div,"Regularity" = indices_eve)
                       
  indic_names_list <- list("Richness" = c('taxo' = 'SR','func' = 'FD','phyl' = 'PD'), 
                           "Divergence" = c('taxo' = '', 'func' = 'FDiv','phyl' = 'MPD'),
                           "Regularity" = c('taxo' = '','func' = 'FEve', 'phyl' = 'VPD'))
                                                                           
  
  ## Initiate plot -------------------------------------------------------------
  image_width <- 20
  image_height <- 20
  res = 800
  
  pdf(file="figures/Fig_3_alpha_div.pdf",width=image_width/2.54,height=image_height/2.54)
  par(mfrow = c(2,2))
  
  par(mar=c(3.2, 3.1, 1.3, 0.15))
  a <- barplot(DF1.SES$Mean,col=c("#ab3902", "#3E9BFE"),axes=F,lwd=0.3,
               ylim=c(min(DF1.SES$Mean - DF1.SES$Se)-0.25,max(DF1.SES$Mean + DF1.SES$Se+0.25)))
  
  
  sapply(1:length(a[,1]),function(i){segments(y0=(DF1.SES$Mean - DF1.SES$Se)[i],
                                              x0=a[i,1],y1=(DF1.SES$Mean + DF1.SES$Se)[i],
                                              x1=a[i,1],lwd=0.5)})
  sapply(1:length(a[,1]),function(i){segments(y0=(DF1.SES$Mean - DF1.SES$Se)[i],
                                              x0=a[i,1]-0.2,y1=(DF1.SES$Mean - DF1.SES$Se)[i],
                                              x1=a[i,1]+0.2,lwd=0.5)})
  sapply(1:length(a[,1]),function(i){segments(y0=(DF1.SES$Mean + DF1.SES$Se)[i],
                                              x0=a[i,1]-0.2,y1=(DF1.SES$Mean + DF1.SES$Se)[i],
                                              x1=a[i,1]+0.2,lwd=0.5)})
  axis(side = 2,las=2,cex.axis=0.75,lwd=0.2,mgp=c(0.35,0.45,0),tcl=-0.15)
  axis(side = 1,las=2, at=c(1.3,3.7,6.1, 8.5, 10.9, 13.3), 
       label=c('PD.ses', 'FD.ses', 'MPD.ses', 'FDiv.ses', 'VPD.ses', 'FEve.ses'),las=1, 
       cex.axis=0.7,mgp=c(0.25,0.35,0),lwd=0.2 )
  abline(h=-1.96,lty=2)
  legend("bottomright",legend=c("eDNA","Trawling"),pt.cex=1.5,pch=15,
         col=c("#ab3902", "#3E9BFE"),title="Methods",cex=0.8,box.lwd=0.2)
  mtext("Biodiversity indices",side = 1,line = 2, at=7.5,cex=0.8)
  mtext("SES values",side = 2,line = 1.25, at=-0.5,cex=0.8)
  mtext("a",side=2,line = 1.25, at=1.53,las=2,cex=1)
  box(lwd=0.2)
  
  
  label_map <- c("b", "c", "d")
  
  ## Draw the 3 panels -------------------------------------------------------
  for (i_panel in 1:3) {
    if (i_panel == 1) {
      par(mar = c(1.2, 1.55, 1.3, 0.15))
    } else {
      par(mar = c(1.2, 1.55, 1.3, 0.15))
    }
    plot(0,0, type='n', xaxt='n', yaxt='n', xlab='', ylab='', bty='n')
    mtext(text=label_map[i_panel],side=3,adj = 0, cex = 1.1)
    par(new = "T")
    
    indices <- indices_list[[i_panel]]
    
    # taxo indice
    df_alpha_t <- RESULTS_tax_trawl$Alpha
    rownames(df_alpha_t) <- paste0(rownames(df_alpha_t), "t")
    
    df_alpha_e <- RESULTS_tax_eDNA$Alpha
    rownames(df_alpha_e) <- paste0(rownames(df_alpha_e), "e")
    
    if (!is.null(indices["taxo"])) {
      df_alpha <- rbind(df_alpha_e, df_alpha_t)
      Points$div_t <- df_alpha[rownames(Points), indices['taxo']]
      
      scale_t <- circle_factor / median(Points$div_t)
      
    } else {
      Points$div_t <- 0
      scale_t <- 0
    }
    
    #functio indice
    df_alpha <- RESULTS_fct_both$Alpha
    df_alpha <- df_alpha[which(rownames(df_alpha) != "total_site"),]
    Points$div_f <- (df_alpha[rownames(Points), indices['func']])
    
    scale_f <- circle_factor / median(Points$div_f)
    
    Points$signif_f <- as.vector(sapply(rownames(Points), function(site) {
      ses <- RESULTS_SES_fct_both[site, indices['func']]
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
    
    # indice phylo
    df_alpha <- RESULTS_phyl_both$Alpha
    df_alpha <- df_alpha[which(rownames(df_alpha) != "total_site"),]
    Points$div_p <- (df_alpha[rownames(Points), indices["phyl"]])
    
    scale_p <- circle_factor / median(Points$div_p)
    
    Points$signif_p <- as.vector(sapply(rownames(Points), function(site) {
      ses <- RESULTS_SES_phyl_both[site,  indices["phyl"]]
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
      if (Pt$type == "eDNA"){
        angles <- pi * c(2/3+eps, 1-eps, 1/3+eps, 2/3-eps, eps, 1/3 - eps)
        color <- color_edna
        dy <- 1
      }else {
        angles <- pi * c(1+eps, 4/3-eps, 4/3+eps, 5/3-eps, 5/3+eps, 2-eps)
        color <- color_trawl
        dy <- -1
      }
      
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
    
    # Points labels 
    if (F) {
      Points_labels <- Points %>% 
        filter(type == "eDNA") %>%
        select(station, Latitude, Longitude)
      rownames(Points_labels) <- Name_stations[Points_labels$station]
      Points_labels$dx <- 0
      Points_labels$dy <- 0
      Points_labels[c(1,2,4,5,6,7,11,14), "dx"] <- 0.35 * cos(3.5)
      Points_labels[c(1,2,4,5,6,7,11,14), "dy"] <- 0.35 * sin(3.5)
      
      Points_labels[c(13), "dx"] <- 0.38 * cos(3.5)
      Points_labels[c(13), "dy"] <- 0.38 * sin(3.5)
      
      Points_labels[c(9,12,15), "dx"] <- 0.15 * cos(5*pi/4)
      Points_labels[c(9,12,15), "dy"] <- 0.15 * sin(5*pi/4)
      
      Points_labels[c(3,8), "dx"] <- 0.13 * cos(0)
      Points_labels[c(3,8), "dy"] <- 0.13 * sin(0)
      
      Points_labels[10, "dx"] <- 0.40
      
      text(x = Points_labels$Longitude + Points_labels$dx,
           y = Points_labels$Latitude + Points_labels$dy,
           labels = rownames(Points_labels),
           col = "black")
    }
    
    ## Legend --------------------------------------------------------------------
    x <- -3.6
    y <- 43.9
    rect(x - 1.49,y-1.4,x+0.6, y + 0.85, col = "white")
    
    
    dy <- 1
    color <- color_edna
    angles <- pi * c(2/3+eps, 1-eps, 1/3+eps, 2/3-eps, eps, 1/3 - eps)
    pie_slice(x = x - dx, 
              y = y + dy* 0.025, 
              r = 1.7*scale_t * median(Points$div_t),
              col = color,
              angle = angles[c(1,2)])
    
    
    pie_slice(x = x, 
              y = y + 0.5*dx+ dy*0.025, 
              r = 1.7*scale_p * median(Points$div_p),
              col = color,
              angle = angles[c(3,4)])
    
    pie_slice(x = x + dx, 
              y = y  + dy *0.025, 
              r = 1.7*scale_f * median(Points$div_f),
              col = color,
              angle = angles[c(5,6)])
    
    
    
    dy <- -1
    color <- color_trawl
    angles <- pi * c(1+eps, 4/3-eps, 4/3+eps, 5/3-eps, 5/3+eps, 2-eps)
    pie_slice(x = x - dx, 
              y = y + dy* 0.025, 
              r = 1.7*scale_t * median(Points$div_t),
              col = color,
              angle = angles[c(1,2)])
    
    
    pie_slice(x = x, 
              y = y -0.5*dx + dy*0.025, 
              r = 1.7*scale_p * median(Points$div_p),
              col = color,
              angle = angles[c(3,4)])
    
    pie_slice(x = x + dx, 
              y = y + dy *0.025, 
              r = 1.7*scale_f * median(Points$div_f),
              col = color,
              angle = angles[c(5,6)])
    text(x - 0.445,
         y = y + 0.63,
         label = names(indices_list)[i_panel],
         cex = 1.1,
         adj = 0.5,
         col = "black",
         font = 2)
    
    text(x - 1,
         y + 0.15,
         label = "eDNA",
         cex = 0.9,
         adj = 0.5,
         col = color_edna,
         font = 2)
    
    text(x - 1,
         y - 0.15,
         label = "trawling",
         cex = 0.9,
         adj = 0.5,
         col = color_trawl,
         font = 2)
    r <- 0.35 
    angles <- c(5*pi/6, 7*pi/6)
    text(x = x + r * cos(angles),
         y = y + r * sin(angles),
         label = indic_names_list[[i_panel]]["taxo"],
         adj = 0.5,
         col = "white",
         cex = 0.62)
    
    angles <- c(pi/2, 3*pi/2)
    text(x = x + r * cos(angles),
         y = y + r * sin(angles),
         label = indic_names_list[[i_panel]]["phyl"],
         adj = 0.5,
         col = "white",
         cex = 0.62)
    angles <- c(pi/6, 11*pi/6)
    text(x = x + r * cos(angles),
         y = y + r * sin(angles),
         label = indic_names_list[[i_panel]]["func"],
         adj = 0.5,
         col = "white",
         cex = 0.62)
    
    text(x - 1.4,
         y - 1,
         label = if (i_panel == 1) {"+ : overdispersed
- : clustered"} else {
  "\n+ : overdispersed
- : clustered"
},
cex = 0.8,
adj = 0)
    
    post_map(xtickslab = (i_panel == 3|i_panel == 2), ytickslab = (i_panel == 1 |i_panel == 2), cex.axis = 0.85, line.xlab = -0.8)
    
  }
  
dev.off()
}

plot_figure()
