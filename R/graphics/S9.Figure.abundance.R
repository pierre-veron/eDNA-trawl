
# Fig S9 - Correlation per filter ----------------------------------------------



S9.Figure.abundance = function(Data.input = NULL, Res.bestModelInput = NULL,   out_file = "figures/Fig_S9_abundance_per_filter", ext = ".pdf", y.var = NULL, y.text = NULL, x.text = NULL, family  = NULL, Filter.Station = NULL){

  
  image_width <- 17
  image_height <- 20
  res = 400
  
  color_repl <- "dodgerblue"
  
  pdf(file = paste0(out_file, ext),
      width = image_width / 2.54,
      height = image_height / 2.54)
  
  par(mfrow = c(5,3))
  par(oma = c(4, 4, 0, 0))
  UniqFilter = unique(Res.bestModelInput$filter)
  
  for (i in 1:length(UniqFilter)) {
    Filter_fish <- Data.input[which(Data.input$filter == UniqFilter[i]),]
    
    
    if (i %in% c(1, 4, 7, 10, 13, 16, 19, 22, 25, 28)) {
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
      ylab <- "Number of reads"
    } else {
      ylab <- ""
    }
    if (i == 14) {
      xlab <- "Allometrically scaled abundance"
    } else {
      xlab <- ""
    }
    par(mar = c(1+mar_dx,1+mar_dy,1,1))
      
      Filter_bestModel = Res.bestModelInput[which(Res.bestModelInput$filter == UniqFilter[i]),]
      
      x.var = Filter_bestModel$nb.col.APTs
      #lreg <- lm(Filter_fish$relat_reads ~ Filter_fish$log_nb_trawl)
      
      lreg_details <- c("intercept" = as.numeric(Filter_bestModel$Intercept), 
                        "slope" = as.numeric(Filter_bestModel$Slope),
                        "r.squared" = Filter_bestModel$R2, 
                        "p.intercept" = NA,
                        "p.slope" = Filter_bestModel$P.value.slope, 
                        "n" = Filter_bestModel$nb.taxa)
      plot(Filter_fish[,x.var], Filter_fish[,y.var], pch = 19, col = "black", xlab = xlab, ylab = ylab, cex.lab = 1.2, main= "")
      
      #points(Filter_fish[,x.var], Filter_fish[,y.var], pch = 19,col = color_repl)
      #abline(lreg, col = color_repl[irep], lty = 2)
      expl.var  = Filter_fish[,x.var]
      #range.minmax = c(min(Filter_fish[,x.var]), max(Filter_fish[,x.var]))
      newx = sort(runif(100,min(Filter_fish[,x.var]), max(Filter_fish[,x.var])))
      
      #newx = seq(range.minmax[1], range.minmax[2],  round((range.minmax[2] - range.minmax[1])/100)) 
      
      newdata = data.frame(expl.var = newx)
      
      
      if(family == "poisson"){

      #lreg <- glm(Filter_fish[,y.var] ~ expl.var, family = "poisson")
      lreg <- glm.nb(Filter_fish[,y.var] ~ expl.var)
      conf_interval <- predict(lreg, newdata, se.fit= TRUE, response = T)
      conf_interval.DF = data.frame(cbind(conf_interval$fit, conf_interval$se.fit))
      conf_interval.DF$lower = conf_interval$fit - conf_interval$se.fit
      conf_interval.DF$upper = conf_interval$fit + conf_interval$se.fit
      conf_interval.DF$X1 = exp(conf_interval.DF$X1)
      conf_interval.DF$lower = exp(conf_interval.DF$lower)
      conf_interval.DF$upper = exp(conf_interval.DF$upper)
      } else {
      lreg <- lm(Filter_fish[,y.var] ~ expl.var)
      conf_interval <- predict(lreg, newdata, se.fit= TRUE, response = T)
      conf_interval.DF = data.frame(cbind(conf_interval$fit, conf_interval$se.fit))
      conf_interval.DF$lower = conf_interval$fit - conf_interval$se.fit
      conf_interval.DF$upper = conf_interval$fit + conf_interval$se.fit
      }

      # To match what can be done with GLS, I used the standard error of the model fit, instead of the confidence intreval at 0.95.
# conf_interval <- predict(lreg, newdata= data.frame(log_nb_trawl = newx), interval="confidence", level = 0.95)       



### Plot the confidence interval at 95%

if(Filter_bestModel$P.value.slope <= 0.05){

lines(newx, conf_interval.DF[,1], col=color_repl, lty=1, lwd = 2) 
#lines(newx, conf_interval.DF[,3], col=color_repl, lty=2) 
#lines(newx, conf_interval.DF[,4], col=color_repl, lty=2) 
} else {
lines(newx, conf_interval.DF[,1], col=color_repl, lty=2, lwd = 1) 
}


      if (lreg_details["p.slope"] < 0.001) {
        lev_signif <- "***"
      } else if (lreg_details["p.slope"] < 0.01) {
        lev_signif <- "**"
      } else if (lreg_details["p.slope"] < 0.05){
        lev_signif <- "*"
      } else {
        lev_signif <- "ns"
      }
      lreg_label <- paste0("R² = ", round(lreg_details["r.squared"], 2),
                           ", ", lev_signif, "\n", "b = ", Res.bestModelInput$best.allom.b[i])
      legend("topright", legend = lreg_label, adj = 0, cex = 0.9, col = color_repl, bty = "n")
      #if (i %in% c(7, 9, 15)) {
      #  text(x = 8.00, y = 0.42 - irep * 0.06, labels = lreg_label, adj = 0, cex = 0.9, col = color_repl)
      #} else {
      #  text(x = 0.00, y = 0.42 - irep * 0.06, labels = lreg_label, adj = 0, cex = 0.9, col = color_repl)
      #}
      
      legend("topleft", legend = paste0("Site ", Filter.Station[i, 3], ", filtre ",  Filter.Station[i,4]), cex = 1 , adj = 0, bty = "n")
    
    #text(x = 0, y = 0.45, labels = paste0("Site ", Filter.Station[i, 3]), cex = 1.3, adj = 0)
    #if (i %in% c(1, 4, 7, 10, 13, 16, 19, 22, 25, 28)) {
    #  axis(2, at = 0.1*(0:5),
    #       las = 2, tck  =F, line = -0.3,lty = 0)
    #} 
    #axis(2,at =  0.1*(0:5), 
    #     labels = F, tck  =-0.03)
    #if (i %in% c(28, 29, 30)) {
    #  axis(1, at =  2*(0:6),
    #       tck  =F, line = -0.3, lty = 0)
    #}
    #axis(1, at =  2*(0:6),
    #     labels = F, tck  =-0.03)
  }
  
  mtext(x.text, side = 1, outer = TRUE, line = 2)
  mtext(y.text, side = 2, outer = TRUE, line = 2)
  
  dev.off()
  }


# For debug
# rm(newx, conf_interval.DF,conf_interval, Filter_fish, lreg_details, lreg, Filter_bestModel)






