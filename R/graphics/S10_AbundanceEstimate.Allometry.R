### Script to estimate the relationships between eDNA copies and Abundance/biomass of fish, and using the best allometric scaling coefficient.


# Correspond to the section called "Qualitative comparison of abundance estimates" in the material and method of the manuscript.



# Additional packages
require(nlme)
require(MASS)
require(AICcmodavg)


### Cox Snell R2 function
R2.Cox.Snell = function(Fitted.model, Null.model, n){
R2.res = 1-exp((-2/n)*(logLik(Fitted.model)[1] - logLik(Null.model)[1]))
return(R2.res)
}


# Load an additional dataset 
Common_taxa.DF = readRDS("data/Abundance.Biomass.Common.eDNA.Trawl.rds")



#### Compute the "allometrically scaled abundance per tow" based on Yates et al. 2022 Envir DNA approach



Nb.scaling.fac = seq(0, 1, 0.01)
av.weight = Common_taxa.DF$biomass.kg/Common_taxa.DF$nb_trawl

APT.DF = data.frame(do.call(cbind, lapply(1:length(Nb.scaling.fac), function(x){
av.weight^Nb.scaling.fac[x] * Common_taxa.DF$nb_trawl
})))
names(APT.DF) = paste("APTb.", Nb.scaling.fac, sep = "")



### A single table with all the APT estimated for the 101 allometric b components.
Common_taxa.DF.APTs = data.frame(Common_taxa.DF, APT.DF)

#### 100 random selection of one of the two replicates filter for each site.
Unique.Sp.ID = unique(Common_taxa.DF.APTs$ID)

Overall.100rand.DF = lapply(1:100, function(x){
do.call(rbind, lapply(1:length(Unique.Sp.ID), function(j){
Common_taxa.DF.APTs[sample(which(Common_taxa.DF.APTs$ID == Unique.Sp.ID[j]),1),]
}))
})






####################################################################################
############## Global approach using Average relative eDNA number of reads between 
############## filters per site (Pierre Veron 's original approach) ################
####################################################################################


# Load the dataset where the number of reads were summed for the two filters per site, and average among filters
Aggreg.Read.mean.DF2 = readRDS("data/Abundance.Biomass.Common.eDNA.Trawl.Agg.Filters.rds")

### Defined the column with the APT scale values estimate for a b parameter from 0 to 1, by 0.01 increment.
Allo.scal.pos = seq(9, 109)

# Fit a linear model
lm.overal.av.rel.read = lm(Aggreg.Read.mean.DF2$av_relat_reads ~ Aggreg.Read.mean.DF2$log_nb_trawl)
summary(lm.overal.av.rel.read)
# Multiple R-squared:  0.1377,    Adjusted R-squared:  0.1337


# Gls model accounting for the heterosdasticity, fitted with ML optimisation
gls.overal.av.rel.read.ML = gls(av_relat_reads ~ log_nb_trawl, data = Aggreg.Read.mean.DF2, weights = varExp(form =~ log_nb_trawl), method = "ML")
AIC(lm.overal.av.rel.read, gls.overal.av.rel.read.ML)
### the GLS model is better

# Gls model accounting for the heterosdasticity, fitted with REML optimisation
gls.overal.av.rel.read = gls(av_relat_reads ~ log_nb_trawl, data = Aggreg.Read.mean.DF2, weights = varExp(form =~ log_nb_trawl))
summary(gls.overal.av.rel.read)


# Gls model null model 1
gls.overal.av.rel.read.null.null = gls(av_relat_reads ~ 1, data = Aggreg.Read.mean.DF2)

# Gls model null model 2
gls.overal.av.rel.read.null = gls(av_relat_reads ~ 1, data = Aggreg.Read.mean.DF2, weights = varExp(form =~ log_nb_trawl))
AIC(gls.overal.av.rel.read, gls.overal.av.rel.read.null, gls.overal.av.rel.read.null.null)


### Cox Snell R2.
R2.Cox.Snell(gls.overal.av.rel.read, gls.overal.av.rel.read.null.null, dim(Aggreg.Read.mean.DF2)[1]) ### 0.2261845



###################################################################################################
#### Test for a potential residual spatial autocorrelation in the models.
head(TRAWL_metadata)

# Get the geographic coordinates of the sites
Aggreg.Read.mean.DF3 = merge(Aggreg.Read.mean.DF2, TRAWL_metadata, by.x = 3, by.y = 1)

### Because multiple species are caught at the same location, and each sampling points should have different coordinates, we jittered a bit the coordinates.
Aggreg.Read.mean.DF4 = Aggreg.Read.mean.DF3
Aggreg.Read.mean.DF4$Longitude =jitter(Aggreg.Read.mean.DF3$Longitude, amount = 0.00001)
Aggreg.Read.mean.DF4$Latitude =jitter(Aggreg.Read.mean.DF3$Latitude, amount = 0.00001)

# Gls model accounting for both the heteroscedasticity and the spatal autocorrelation of the rsiduals using a corExp structure with a nugget effect.
gls.overal.av.rel.read.Space.corExp = gls(av_relat_reads ~ log_nb_trawl, data = Aggreg.Read.mean.DF4, weights = varExp(form =~ log_nb_trawl), correlation = corExp(form =~ Longitude + Latitude, nugget = TRUE))
summary(gls.overal.av.rel.read.Space.corExp)

# Gls model accounting for both the heteroscedasticity and the spatal autocorrelation of the rsiduals using a corGaus structure with a nugget effect.
gls.overal.av.rel.read.Space.corGaus = gls(av_relat_reads ~ log_nb_trawl, data = Aggreg.Read.mean.DF4, weights = varExp(form =~ log_nb_trawl), correlation = corGaus(form =~ Longitude + Latitude, nugget = TRUE))

summary(gls.overal.av.rel.read.Space.corGaus)


# Gls model accounting for both the heteroscedasticity and the spatal autocorrelation of the rsiduals using a corSpher structure with a nugget effect.
gls.overal.av.rel.read.Space.corSpher = gls(av_relat_reads ~ log_nb_trawl, data = Aggreg.Read.mean.DF4, weights = varExp(form =~ log_nb_trawl), correlation = corSpher(form =~ Longitude + Latitude, nugget = TRUE))
summary(gls.overal.av.rel.read.Space.corSpher)


### Model selection of the best stuctrure of the spatial autocorrelation.
AIC(gls.overal.av.rel.read, gls.overal.av.rel.read.Space.corExp, gls.overal.av.rel.read.Space.corGaus, gls.overal.av.rel.read.Space.corSpher)
### The three structures accounting for the spatial autocorrelation provided the same results and did not improve the AIC significantly compare to the model without spatial autocorrelation structure.

# We choose the corExp
# refit with the ML optimisation to make it comparable with LM
gls.overal.av.rel.read.Space.corExp.ML = gls(av_relat_reads ~ log_nb_trawl, data = Aggreg.Read.mean.DF4, weights = varExp(form =~ log_nb_trawl), correlation = corExp(form =~ Longitude + Latitude, nugget = TRUE), method  = "ML")

summary(gls.overal.av.rel.read.Space.corExp.ML)

AIC(lm.overal.av.rel.read, gls.overal.av.rel.read.ML, gls.overal.av.rel.read.Space.corExp.ML)
### So the best model includes only a variance structure in the residuals for the heteroscedaticity and not for the spatial autocorrelation.

# Null model with spatial correlation and variance structure
gls.overal.av.rel.read.Space.corExp.ML.null = gls(av_relat_reads ~ 1, data = Aggreg.Read.mean.DF4, weights = varExp(form =~ log_nb_trawl), correlation = corExp(form =~ Longitude + Latitude, nugget = TRUE), method  = "ML")

# Null model with variance structure only
gls.overal.av.rel.read.Space.corExp.ML.null.null = gls(av_relat_reads ~ 1, data = Aggreg.Read.mean.DF4, weights = varExp(form =~ log_nb_trawl), method  = "ML")

# Null model with intercept only
gls.overal.av.rel.read.Space.corExp.ML.null.null.null = gls(av_relat_reads ~ 1, data = Aggreg.Read.mean.DF4, method  = "ML")


R2.Cox.Snell(gls.overal.av.rel.read.Space.corExp.ML,gls.overal.av.rel.read.Space.corExp.ML.null,dim(Aggreg.Read.mean.DF4)[1]) ### 0.03306472 explained by the fixed effect only 

R2.Cox.Snell(gls.overal.av.rel.read.Space.corExp.ML,gls.overal.av.rel.read.Space.corExp.ML.null.null,dim(Aggreg.Read.mean.DF4)[1]) ### 0.04654265 explained by the fixed effect and spatial autocorrelation

R2.Cox.Snell(gls.overal.av.rel.read.Space.corExp.ML,gls.overal.av.rel.read.Space.corExp.ML.null.null.null,dim(Aggreg.Read.mean.DF4)[1]) ### 0.2700532explained by the fixed effect and spatial autocorrelation




### Modif Fg. S8 from Pierre Veron's OriginalFig.

# Fig S8 - Correlation all filters ---------------------------------------------
  # building a dataframe that contains the average and difference of the relative
  # number of reads in each sampling site, between the two filters 
  # and the trawling catches 


  Common_taxa_merged_filt <- Common_taxa.DF
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
  
  
  out_file <- "figures/Fig_S8_abundance.NEW"
  ext <- ".pdf"
  
  image_width <- 17
  image_height <- 15
  res = 400
  
  
  # pdf(file = paste0(out_file, ext),
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
  #abline(lreg, col = "dodgerblue", lt = 1, lwd = 2)
  
newx = sort(Common_taxa_merged_filt$log_nb_trawl)
# To match what can be done with GLS, I used the standard error of the model fit, instead of the confidence interval at 0.95.
# conf_interval <- predict(lreg, newdata= data.frame(log_nb_trawl = newx), interval="confidence", level = 0.95)       

conf_interval <- predict(lreg, newdata= data.frame(log_nb_trawl = newx), se.fit= TRUE)

conf_interval.DF = data.frame(cbind(conf_interval$fit, conf_interval$se.fit))
conf_interval.DF$lower = conf_interval$fit - conf_interval$se.fit
conf_interval.DF$upper = conf_interval$fit + conf_interval$se.fit
           
conf_interval.DF = conf_interval.DF[order(conf_interval.DF[,1]),]                       

### Plot the confidence interval at 95%
lines(newx, conf_interval.DF[,1], col="dodgerblue", lty=1, lwd = 2)           
lines(newx, conf_interval.DF[,3], col="dodgerblue", lty=2)
lines(newx, conf_interval.DF[,4], col="dodgerblue", lty=2)
  
  
#### with a gls model
gls.reg = gls(av_relat_reads ~ log_nb_trawl, data = Common_taxa_merged_filt, weights = varExp(form=~log_nb_trawl))
  
newx = data.frame(log_nb_trawl = sort(Common_taxa_merged_filt$log_nb_trawl))
# conf_interval = intervals(gls.reg, level = 0.95)

conf_interval <- AICcmodavg::predictSE(gls.reg, newdata= newx, se.fit = TRUE)  
conf_interval.DF = data.frame(cbind(conf_interval$fit, conf_interval$se.fit))
conf_interval.DF$lower = conf_interval$fit - conf_interval$se.fit
conf_interval.DF$upper = conf_interval$fit + conf_interval$se.fit
  
conf_interval.DF = conf_interval.DF[order(conf_interval.DF[,1]),]                       
 

lines(newx[,1], conf_interval.DF[,1], col="orange", lty=1, lwd = 2)           
lines(newx[,1], conf_interval.DF[,3], col="orange", lty=2)
lines(newx[,1], conf_interval.DF[,4], col="orange", lty=2)
  
  
  
  legend("bottomright", legend = c("LM", "GLS"), lty = 1, lwd = 2, col = c("dodgerblue", "orange"))
  
  axis(side = 1, at = 2*(0:6), line = 0, tck = -0.015, labels = F)
  axis(side = 1, at = 2*(0:6), line = -0.6, tick = F, lty = 0, labels = T)
  
  axis(side = 2, at = 0.1*(0:4), line = 0, tck = -0.015, labels = F)
  axis(side = 2, at = 0.1*(0:4), line = -0.3, tick = F, lty = 0, labels = T, las = T)
  mtext("Log number of individuals per trawl", side = 1, line = 1.4)
  mtext("Relative abundance of eDNA reads per taxa", side = 2, line = 1.9)
  
dev.off()
  

gls.reg = gls(av_relat_reads ~ log_nb_trawl, data = Common_taxa_merged_filt, weights = varExp(form=~log_nb_trawl))

gls.reg.null.null = gls(av_relat_reads ~ 1, data = Common_taxa_merged_filt)
AIC(gls.reg, gls.reg.null.null)
summary(gls.reg)
summary(gls.reg.null.null)


## Make more sens
R2.Cox.Snell(gls.reg, gls.reg.null.null, dim(Common_taxa_merged_filt)[1])
# 0.2261845


gls.reg.null = gls(av_relat_reads ~ 1, data = Common_taxa_merged_filt, weights = varExp(form=~log_nb_trawl))

R2.Cox.Snell(gls.reg, gls.reg.null, dim(Common_taxa_merged_filt)[1])

gls.reg.ML = gls(av_relat_reads ~ log_nb_trawl, data = Common_taxa_merged_filt, weights = varExp(form=~log_nb_trawl), method = "ML")

gls.reg.null.ML = gls(av_relat_reads ~ 1, data = Common_taxa_merged_filt, weights = varExp(form=~log_nb_trawl), method = "ML")

# Cox Snell R2 related to the fixed effect only.
R2.Cox.Snell(gls.reg.ML, gls.reg.null.ML , dim(Common_taxa_merged_filt)[1])






###### Plot the Same Figure S8 but without the Trachurus sp. as requested by the Reviewer for ICES
###### (This figures was used only to address a reviewer's comment and did not appear in the main doucment or the supp. mat.)

# Remove the Trachurus
unique(Common_taxa_merged_filt$fish)
 
Common_taxa_merged_filt.noTrach = Common_taxa_merged_filt[-which(Common_taxa_merged_filt$fish == "Trachurus"),]


 out_file <- "figures/Fig_S8_abundance.NEW.revision"
  ext <- ".pdf"
  
  image_width <- 34
  image_height <- 18
  res = 400
  
  
  pdf(file = paste0(out_file, ext),
      width = image_width / 2.54,
      height = image_height / 2.54)      


  par(mfrow = c(1,2), mar = c(2.5,2.8,3,1))

### Figure With Trachurus 


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
  #abline(lreg, col = "dodgerblue", lt = 1, lwd = 2)
  
newx = sort(Common_taxa_merged_filt$log_nb_trawl)
# To match what can be done with GLS, I used the standard error of the model fit, instead of the confidence interval at 0.95.
# conf_interval <- predict(lreg, newdata= data.frame(log_nb_trawl = newx), interval="confidence", level = 0.95)       

conf_interval <- predict(lreg, newdata= data.frame(log_nb_trawl = newx), se.fit= TRUE)

conf_interval.DF = data.frame(cbind(conf_interval$fit, conf_interval$se.fit))
conf_interval.DF$lower = conf_interval$fit - conf_interval$se.fit
conf_interval.DF$upper = conf_interval$fit + conf_interval$se.fit
           
conf_interval.DF = conf_interval.DF[order(conf_interval.DF[,1]),]                       

### Plot the confidence interval at 95%
lines(newx, conf_interval.DF[,1], col="dodgerblue", lty=1, lwd = 2)           
lines(newx, conf_interval.DF[,3], col="dodgerblue", lty=2)
lines(newx, conf_interval.DF[,4], col="dodgerblue", lty=2)
  
  
#### with a gls model
  gls.reg = gls(av_relat_reads ~ log_nb_trawl, data = Common_taxa_merged_filt, weights = varExp(form=~log_nb_trawl))
  
  newx = data.frame(log_nb_trawl = sort(Common_taxa_merged_filt$log_nb_trawl))
  # conf_interval = intervals(gls.reg, level = 0.95)

conf_interval <- AICcmodavg::predictSE(gls.reg, newdata= newx, se.fit = TRUE)  
conf_interval.DF = data.frame(cbind(conf_interval$fit, conf_interval$se.fit))
conf_interval.DF$lower = conf_interval$fit - conf_interval$se.fit
conf_interval.DF$upper = conf_interval$fit + conf_interval$se.fit
  
conf_interval.DF = conf_interval.DF[order(conf_interval.DF[,1]),]                       
 

lines(newx[,1], conf_interval.DF[,1], col="orange", lty=1, lwd = 2)           
lines(newx[,1], conf_interval.DF[,3], col="orange", lty=2)
lines(newx[,1], conf_interval.DF[,4], col="orange", lty=2)
  
  
  
  legend("bottomright", legend = c("LM", "GLS"), lty = 1, lwd = 2, col = c("dodgerblue", "orange"))
  
  axis(side = 1, at = 2*(0:6), line = 0, tck = -0.015, labels = F)
  axis(side = 1, at = 2*(0:6), line = -0.6, tick = F, lty = 0, labels = T)
  
  axis(side = 2, at = 0.1*(0:4), line = 0, tck = -0.015, labels = F)
  axis(side = 2, at = 0.1*(0:4), line = -0.3, tick = F, lty = 0, labels = T, las = T)
  mtext("Log number of individuals per trawl", side = 1, line = 1.4)
  mtext("Relative abundance of eDNA reads per taxa", side = 2, line = 1.9)
  
  mtext(substitute(paste("With", italic("Trachurus"))),  , side = 3, line = 0.9)


### Figure Without Trachurus 
  plot(Common_taxa_merged_filt.noTrach$log_nb_trawl, Common_taxa_merged_filt.noTrach$av_relat_reads, pch = 20, 
       xlab = "",
       ylab = "", 
       ylim = c(0, max( Common_taxa_merged_filt.noTrach$av_relat_reads + 0.5 * Common_taxa_merged_filt.noTrach$delta_relat_reads)), 
       xaxt = "n" , yaxt = "n")
  Hmisc::errbar(x = Common_taxa_merged_filt.noTrach$log_nb_trawl, y =  Common_taxa_merged_filt.noTrach$av_relat_reads,
                yplus = Common_taxa_merged_filt.noTrach$av_relat_reads + 0.5 * Common_taxa_merged_filt.noTrach$delta_relat_reads,
                yminus =  Common_taxa_merged_filt.noTrach$av_relat_reads - 0.5 * Common_taxa_merged_filt.noTrach$delta_relat_reads, 
                add = T,cap = 0)
  
  
  lreg <- lm(Common_taxa_merged_filt.noTrach$av_relat_reads ~ Common_taxa_merged_filt.noTrach$log_nb_trawl)
  print(summary(lreg))
  #abline(lreg, col = "dodgerblue", lt = 1, lwd = 2)
  
newx = sort(Common_taxa_merged_filt.noTrach$log_nb_trawl)
# To match what can be done with GLS, I used the standard error of the model fit, instead of the confidence interval at 0.95.
# conf_interval <- predict(lreg, newdata= data.frame(log_nb_trawl = newx), interval="confidence", level = 0.95)       

conf_interval <- predict(lreg, newdata= data.frame(log_nb_trawl = newx), se.fit= TRUE)

conf_interval.DF = data.frame(cbind(conf_interval$fit, conf_interval$se.fit))
conf_interval.DF$lower = conf_interval$fit - conf_interval$se.fit
conf_interval.DF$upper = conf_interval$fit + conf_interval$se.fit
           
conf_interval.DF = conf_interval.DF[order(conf_interval.DF[,1]),]                       

### Plot the confidence interval at 95%
lines(newx, conf_interval.DF[,1], col="dodgerblue", lty=1, lwd = 2)           
lines(newx, conf_interval.DF[,3], col="dodgerblue", lty=2)
lines(newx, conf_interval.DF[,4], col="dodgerblue", lty=2)
  
  
#### with a gls model
  gls.reg = gls(av_relat_reads ~ log_nb_trawl, data = Common_taxa_merged_filt.noTrach, weights = varExp(form=~log_nb_trawl))
  
  newx = data.frame(log_nb_trawl = sort(Common_taxa_merged_filt.noTrach$log_nb_trawl))
  # conf_interval = intervals(gls.reg, level = 0.95)

conf_interval <- AICcmodavg::predictSE(gls.reg, newdata= newx, se.fit = TRUE)  
conf_interval.DF = data.frame(cbind(conf_interval$fit, conf_interval$se.fit))
conf_interval.DF$lower = conf_interval$fit - conf_interval$se.fit
conf_interval.DF$upper = conf_interval$fit + conf_interval$se.fit
  
conf_interval.DF = conf_interval.DF[order(conf_interval.DF[,1]),]                       
 

lines(newx[,1], conf_interval.DF[,1], col="orange", lty=2, lwd = 2)           
lines(newx[,1], conf_interval.DF[,3], col="orange", lty=2)
lines(newx[,1], conf_interval.DF[,4], col="orange", lty=2)
  
  
  
  legend("bottomright", legend = c("LM", "GLS"), lty = 1, lwd = 2, col = c("dodgerblue", "orange"))
  
  axis(side = 1, at = 2*(0:6), line = 0, tck = -0.015, labels = F)
  axis(side = 1, at = 2*(0:6), line = -0.6, tick = F, lty = 0, labels = T)
  
  axis(side = 2, at = 0.1*(0:4), line = 0, tck = -0.015, labels = F)
  axis(side = 2, at = 0.1*(0:4), line = -0.3, tick = F, lty = 0, labels = T, las = T)
  mtext("Log number of individuals per trawl", side = 1, line = 1.4)
  mtext("Relative abundance of eDNA reads per taxa", side = 2, line = 1.9)
  
  gls.reg.null.null = gls(av_relat_reads ~ 1, data = Common_taxa_merged_filt.noTrach)
  r2 = R2.Cox.Snell(gls.reg, gls.reg.null.null, dim(Common_taxa_merged_filt)[1])
  
  
  mtext(expression("Without" ~italic("Trachurus")),  , side = 3, line = 0.9)
  
  
dev.off()
 

lreg <- lm(Common_taxa_merged_filt.noTrach$av_relat_reads ~ Common_taxa_merged_filt.noTrach$log_nb_trawl)
summary(lreg)

gls.reg = gls(av_relat_reads ~ log_nb_trawl, data = Common_taxa_merged_filt.noTrach, weights = varExp(form=~log_nb_trawl))
summary(gls.reg)

gls.reg.null.null = gls(av_relat_reads ~ 1, data = Common_taxa_merged_filt.noTrach)
AIC(gls.reg, gls.reg.null.null)
summary(gls.reg)

## Make more sens
R2.Cox.Snell(gls.reg, gls.reg.null.null, dim(Common_taxa_merged_filt.noTrach)[1])
# 0.09981877







#####################################################################################################
#####################################################################################################
### Local models by filter 
#####################################################################################################


#######################################################################################################
#### NEW Figure S9 with new color and avoid extrapolation, and make the distinction between significant and non significant relationships



  # Fig S9 - Correlation per filter ----------------------------------------------
  out_file <- "figures/Fig_S9_abundance_per_filter.NEW"
  
  ext <- ".pdf"
  
  image_width <- 17
  image_height <- 20
  res = 400
  
  color_repl <- c("dodgerblue", "orange")
  
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
      expl.var = Filter_fish$log_nb_trawl
      lreg <- lm(Filter_fish$relat_reads ~ expl.var)
      
      lreg_details <- c("intercept" = as.numeric(lreg$coefficients[1]), 
                        "slope" = as.numeric(lreg$coefficients[2]),
                        "r.squared" = summary(lreg)[["r.squared"]], 
                        "p.intercept" = summary(lreg)$coefficients[1,4],
                        "p.slope" = summary(lreg)$coefficients[2,4], 
                        "n" = nrow(Filter_fish))
      
      points(Filter_fish$log_nb_trawl, Filter_fish$relat_reads, pch = 19,col = color_repl[irep])
      
      newx = sort(runif(100,min(Filter_fish$log_nb_trawl), max(Filter_fish$log_nb_trawl)))
      newdata = data.frame(expl.var = newx)
      conf_interval <- predict(lreg, newdata, se.fit= TRUE, response = T)
      conf_interval.DF = data.frame(cbind(conf_interval$fit, conf_interval$se.fit))
      conf_interval.DF$lower = conf_interval$fit - conf_interval$se.fit
      conf_interval.DF$upper = conf_interval$fit + conf_interval$se.fit
      
      if(summary(lreg)$coefficients[2,4] <= 0.05){

      lines(newx, conf_interval.DF[,1], col=color_repl[irep], lty=1, lwd = 2) 
      } else {
      lines(newx, conf_interval.DF[,1], col=color_repl[irep], lty=2, lwd = 1) 
      }

            
      #abline(lreg, col = color_repl[irep], lty = 2)
      
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
  
  mtext("Log number of individuals per trawl", side = 1, outer = TRUE, line = 2)
  mtext("Relative abundance of eDNA reads per taxa and per filter", side = 2, outer = TRUE, line = 2)
  
  dev.off()











#####################################################################################################
#####################################################################################################
### Local model by filter and estimate the best Allometric coefficient by filter
#####################################################################################################

Common_taxa.DF.APTs[1:5, 1:10]

UniqFilter = unique(Common_taxa.DF.APTs$filter)
Allo.scal.pos = c(11:111)


######################################################################################################
### Nb read


###### Model fit with a GLM Poisson
PerFilter.Nb_read = lapply(1:length(UniqFilter), function(x){
DFtemp = Common_taxa.DF.APTs[Common_taxa.DF.APTs$filter == UniqFilter[x],]
res1 = do.call(rbind, lapply(1:length(Allo.scal.pos), function(j){
lm1= glm(DFtemp$reads ~ DFtemp[, Allo.scal.pos[j]], family=poisson()) # glm with a poisson model
R2 = (summary(lm1)$null.deviance - summary(lm1)$deviance)/summary(lm1)$null.deviance
c(AIC(lm1), R2, summary(lm1)$coefficients[1,1], summary(lm1)$coefficients[2,1], summary(lm1)$coefficients[2,4], dim(DFtemp)[1])
}))
colnames(res1) = c("AIC", "R2", "Intercept", "Slope", "P.value.slope", "nb.taxa")
res1
})

head(PerFilter.Nb_read[[1]])

BestAllo.Perfilter.nb_reads = unlist(lapply(1:length(PerFilter.Nb_read), function(x){
which.min(PerFilter.Nb_read[[x]][,1])
}))

NB.taxa.perFilter = as.vector(unlist(lapply(1:length(PerFilter.Nb_read), function(x){
PerFilter.Nb_read[[x]][1,6]
})))


BestModel.Perfilter.nb_reads = do.call(rbind,lapply(1:length(PerFilter.Nb_read), function(x){
PerFilter.Nb_read[[x]][which.min(PerFilter.Nb_read[[x]][,1]),]
}))

BestModel.Perfilter.nb_reads.DF = data.frame(unique(Common_taxa.DF.APTs[,c(3,5, 6)]) ,  BestModel.Perfilter.nb_reads, best.allom.b = Nb.scaling.fac[BestAllo.Perfilter.nb_reads], nb.col.APTs = Allo.scal.pos[BestAllo.Perfilter.nb_reads])


mean(BestModel.Perfilter.nb_reads[,2]) # average R2 of the best model for each filter = 0.2927086





######################################################################################################
#### Model fit with a GLM. nb (negative bionomial instead of a poissons model (retained in the paper)


### Nb read
PerFilter.Nb_read.glm.nb = lapply(1:length(UniqFilter), function(x){
DFtemp = Common_taxa.DF.APTs[Common_taxa.DF.APTs$filter == UniqFilter[x],]
res1 = do.call(rbind, lapply(1:length(Allo.scal.pos), function(j){
lm1= tryCatch(glm.nb(DFtemp$reads ~ DFtemp[, Allo.scal.pos[j]]), error=function(e) "error") # glm with a poisson model
model.type = 1
if(class(lm1)[1] == "character"){
lm1= glm(DFtemp$reads ~ DFtemp[, Allo.scal.pos[j]], family=poisson()) # glm with a poisson model
model.type = 0
}

R2 = (summary(lm1)$null.deviance - summary(lm1)$deviance)/summary(lm1)$null.deviance
c(AIC(lm1), R2, summary(lm1)$coefficients[1,1], summary(lm1)$coefficients[2,1], summary(lm1)$coefficients[2,4], dim(DFtemp)[1], model.type)
}))
colnames(res1) = c("AIC", "R2", "Intercept", "Slope", "P.value.slope", "nb.taxa", "model.type")
res1
})


BestAllo.PerFilter.Nb_read.glm.nbs = unlist(lapply(1:length(PerFilter.Nb_read.glm.nb), function(x){
which.min(PerFilter.Nb_read.glm.nb[[x]][,1])
}))

NB.taxa.perFilter = as.vector(unlist(lapply(1:length(PerFilter.Nb_read.glm.nb), function(x){
PerFilter.Nb_read.glm.nb[[x]][1,6]
})))


BestModel.PerFilter.Nb_read.glm.nbs = do.call(rbind,lapply(1:length(PerFilter.Nb_read.glm.nb), function(x){
PerFilter.Nb_read.glm.nb[[x]][which.min(PerFilter.Nb_read.glm.nb[[x]][,1]),]
}))

BestModel.PerFilter.Nb_read.glm.nbs.DF = data.frame(unique(Common_taxa.DF.APTs[,c(3,5, 6)]) ,  BestModel.PerFilter.Nb_read.glm.nbs, best.allom.b = Nb.scaling.fac[BestAllo.PerFilter.Nb_read.glm.nbs], nb.col.APTs = Allo.scal.pos[BestAllo.PerFilter.Nb_read.glm.nbs])


mean(BestModel.PerFilter.Nb_read.glm.nbs[,2]) # average R2 of the best model for each filter. 0.1787705

### These model makes much more sens, despite a lower average R2, the higher R2s for the Poisson models were driven by several extreme values only.




######################################################################################################
######################################################################################################
#### Build a new figure for each relationships between the best allometric coefficient  and the raw number of read for each filter and sites separately. 
######################################################################################################

### This figure will be the Figure 10



UniqFilter = unique(BestModel.Perfilter.nb_reads.DF$filter)
  Station = rep(1:15, each = 2)
  Filtre = rep(1:2, 15)
  Filter.station = cbind(unique(cbind(BestModel.Perfilter.nb_reads.DF$filter, BestModel.Perfilter.nb_reads.DF$station)), Station.code = Station, Filtre.repl = Filtre)



######################################################################################################
#### Based on glm.nb  (retained in the mansucript).
#part 1
BestModel.PerFilter.Nb_read.glm.nbs.DF.part1 = BestModel.PerFilter.Nb_read.glm.nbs.DF[1:15, ]

Common_taxa.DF.APTs.part1 = Common_taxa.DF.APTs[Common_taxa.DF.APTs$filter %in% BestModel.PerFilter.Nb_read.glm.nbs.DF.part1$filter,]


S9.Figure.abundance(Data.input = Common_taxa.DF.APTs.part1, Res.bestModelInput = BestModel.PerFilter.Nb_read.glm.nbs.DF.part1,   out_file = "figures/Fig_S10.Part1", ext = ".pdf", y.var = 8, y.text = "Number of eDNA reads per taxa", x.text = "allometrically scaled abundance", family = "poisson", Filter.Station = Filter.station[1:15,])



#part 2
BestModel.PerFilter.Nb_read.glm.nbs.DF.part2 = BestModel.PerFilter.Nb_read.glm.nbs.DF[16:30, ]

Common_taxa.DF.APTs.part2 = Common_taxa.DF.APTs[Common_taxa.DF.APTs$filter %in% BestModel.PerFilter.Nb_read.glm.nbs.DF.part2$filter,]


S9.Figure.abundance(Data.input = Common_taxa.DF.APTs.part2, Res.bestModelInput = BestModel.PerFilter.Nb_read.glm.nbs.DF.part2,   out_file = "figures/Fig_S10.Part2", ext = ".pdf", y.var = 8, y.text = "Number of eDNA reads per taxa", x.text = "allometrically scaled abundance", family = "poisson", Filter.Station = Filter.station[16:30,])


BestModel.PerFilter.Nb_read.glm.nbs.DF = as.data.frame(BestModel.PerFilter.Nb_read.glm.nbs.DF)






######################################################################################################
###### Identify the outlier species driving the allometric relationships fitted by the glm.nb model.

Identification.outlier.PerFilter = do.call(rbind, lapply(1:dim(BestModel.PerFilter.Nb_read.glm.nbs.DF)[1], function(x){
DFtemp = Common_taxa.DF.APTs[which(Common_taxa.DF.APTs$station == BestModel.PerFilter.Nb_read.glm.nbs.DF[x,"station"] & Common_taxa.DF.APTs$filter == 
BestModel.PerFilter.Nb_read.glm.nbs.DF[x, "filter"]),]
DFtemp[which.max(DFtemp[, BestModel.PerFilter.Nb_read.glm.nbs.DF[x, "nb.col.APTs"]]), c("fish", "reads", "relat_reads")]
}))
Identification.outlier.PerFilter.DF=data.frame(BestModel.PerFilter.Nb_read.glm.nbs.DF[,c(1:2)],Filter.station[,c(3:4)], Identification.outlier.PerFilter)

BestModel.PerFilter.Nb_read.glm.nbs.DF2 = data.frame(BestModel.PerFilter.Nb_read.glm.nbs.DF, Filter.station[,c(3:4)], Identification.outlier.PerFilter)

### Oultier taxa of all the significant and positive relationships
BestModel.PerFilter.Nb_read.glm.nbs.DF2[which(BestModel.PerFilter.Nb_read.glm.nbs.DF2$Slope > 0 & BestModel.PerFilter.Nb_read.glm.nbs.DF2$P.value.slope <= 0.05), c("Station.code", "Filtre.repl", "R2", "nb.taxa", "best.allom.b", "fish", "reads", "relat_reads")]
# Trachurus comes 5 times and Capros Aper 1 

### Oultier taxa of all the positive relationships
Oultier.All.Positive = BestModel.PerFilter.Nb_read.glm.nbs.DF2[which(BestModel.PerFilter.Nb_read.glm.nbs.DF2$Slope > 0), c("Station.code", "Filtre.repl", "R2", "nb.taxa", "best.allom.b", "fish", "reads", "relat_reads")]

table(Oultier.All.Positive$fish)


### Oultier taxa of all the significant and negative relationships
BestModel.PerFilter.Nb_read.glm.nbs.DF2[which(BestModel.PerFilter.Nb_read.glm.nbs.DF2$Slope < 0 & BestModel.PerFilter.Nb_read.glm.nbs.DF2$P.value.slope <= 0.05), c("Station.code", "Filtre.repl", "R2", "nb.taxa", "best.allom.b", "fish", "reads", "relat_reads")]


### Oultier taxa of all the negative relationships
Oultier.All.Negative = BestModel.PerFilter.Nb_read.glm.nbs.DF2[which(BestModel.PerFilter.Nb_read.glm.nbs.DF2$Slope < 0), c("Station.code", "Filtre.repl", "R2", "nb.taxa", "best.allom.b", "fish", "reads", "relat_reads")]

table(Oultier.All.Negative$fish)









######################################################################################################
### Relative number of reads.
PerFilter.rel_read = lapply(1:length(UniqFilter), function(x){
DFtemp = Common_taxa.DF.APTs[Common_taxa.DF.APTs$filter == UniqFilter[x],]
res1 = do.call(rbind, lapply(1:length(Allo.scal.pos), function(j){
lm1= lm(DFtemp$relat_reads ~ DFtemp[, Allo.scal.pos[j]]) # glm with a poisson model
#R2 = (summary(lm1)$null.deviance - summary(lm1)$deviance)/summary(lm1)$null.deviance
c(AIC(lm1), summary(lm1)$r.squared, summary(lm1)$coefficients[1,1], summary(lm1)$coefficients[2,1], summary(lm1)$coefficients[2,4], dim(DFtemp)[1])
}))
colnames(res1) = c("AIC", "R2", "Intercept", "Slope", "P.value.slope", "nb.taxa")
res1
})

BestAllo.Perfilter.rel_reads = unlist(lapply(1:length(PerFilter.rel_read), function(x){
which.min(PerFilter.rel_read[[x]][,1])
}))

table(BestAllo.Perfilter.rel_reads)


BestModel.Perfilter.rel_reads = do.call(rbind,lapply(1:length(PerFilter.rel_read), function(x){
PerFilter.rel_read[[x]][which.min(PerFilter.rel_read[[x]][,1]),]
}))

mean(BestModel.Perfilter.rel_reads[,2])


BestModel.Perfilter.rel_reads.DF = data.frame(unique(Common_taxa.DF.APTs[,c(3,5, 6)]) ,  BestModel.Perfilter.rel_reads, best.allom.b = Nb.scaling.fac[BestAllo.Perfilter.rel_reads], nb.col.APTs = Allo.scal.pos[BestAllo.Perfilter.rel_reads])

table(BestModel.Perfilter.rel_reads.DF$best.allom.b)



######## A figure S11 with individuals relationships with the relative number of reads.

# part 1
BestModel.Perfilter.rel_reads.part1 = BestModel.Perfilter.rel_reads.DF[1:15, ]

Common_taxa.DF.APTs.part1 = Common_taxa.DF.APTs[Common_taxa.DF.APTs$filter %in% BestModel.Perfilter.rel_reads.part1$filter,]


S9.Figure.abundance(Data.input = Common_taxa.DF.APTs.part1, Res.bestModelInput = BestModel.Perfilter.rel_reads.part1,   out_file = "figures/Fig_S11.Part1", ext = ".pdf", y.var = 4, y.text = "Relative abundance of eDNA reads per taxaa", x.text = "allometrically scaled abundance", family = "gaussian", Filter.Station = Filter.station[1:15,])


#part 2
BestModel.Perfilter.rel_reads.part2 = BestModel.Perfilter.rel_reads.DF[16:30, ]

Common_taxa.DF.APTs.part2 = Common_taxa.DF.APTs[Common_taxa.DF.APTs$filter %in% BestModel.Perfilter.rel_reads.part2$filter,]


S9.Figure.abundance(Data.input = Common_taxa.DF.APTs.part2, Res.bestModelInput = BestModel.Perfilter.rel_reads.part2,   out_file = "figures/Fig_S11.Part2", ext = ".pdf", y.var = 4, y.text = "Relative abundance of eDNA reads per taxa", x.text = "allometrically scaled abundance", family = "gaussian", Filter.Station = Filter.station[16:30,])





##### Export a single table with bets allometric coefficients with GLM.NB (number of reads) and the best allometric coefficient with the relative number of read fitted with a lm model.

d1 = data.frame(Response.Var = rep("Nb. eDNA reads", dim(BestModel.PerFilter.Nb_read.glm.nbs.DF)[1]), BestModel.PerFilter.Nb_read.glm.nbs.DF[,c(1:2)], Filter.station[,c(3:4)], round(BestModel.PerFilter.Nb_read.glm.nbs.DF[,c(4:9, 11)], digits = 3))


d2 = data.frame(Response.Var = rep("Relative nb. eDNA reads", dim(BestModel.Perfilter.rel_reads.DF)[1]), BestModel.Perfilter.rel_reads.DF[,c(1:2)], Filter.station[,c(3:4)], round(BestModel.Perfilter.rel_reads.DF[,c(4:9, 11)], digits = 3))
names(d2) = names(d1)

write.table(data.frame(rbind(d1, d2)), file = "output/RESULTS_TABLE.S3.Best_AllometricScaledAbundance.PerFilter.GLM.NB.And.RelNb.eDNAreads.txt", sep ="\t", row.names = F, quote = F)






#################################################################################################
#################################################################################################
#### Basal comparison with the approach from Pierre Veron relat-read ~ log10(nb_individuals)

PerFilter.rel_read.vs.log.nb.indiv = do.call(rbind, lapply(1:length(UniqFilter), function(x){
DFtemp = Common_taxa.DF.APTs[Common_taxa.DF.APTs$filter == UniqFilter[x],]
lm1= lm(DFtemp$relat_reads ~ DFtemp$log_nb_trawl) # glm with a poisson model
#R2 = (summary(lm1)$null.deviance - summary(lm1)$deviance)/summary(lm1)$null.deviance
res1 = c(AIC(lm1), summary(lm1)$r.squared, summary(lm1)$coefficients[1,1], summary(lm1)$coefficients[2,1], summary(lm1)$coefficients[2,4], dim(DFtemp)[1])
names(res1) = c("AIC", "R2", "Intercept", "Slope", "P.value.slope", "nb.taxa")
res1
}))

head(PerFilter.rel_read.vs.log.nb.indiv, n = 10)



length(which((BestModel.Perfilter.rel_reads[,2] - PerFilter.rel_read.vs.log.nb.indiv[,2]) > 0))/30
# so 66.66667% of the case the best allometric scaling model explained better that the relative proportion of read vs log10 (nb individuals)

(mean(BestModel.Perfilter.rel_reads[,2]) - mean(PerFilter.rel_read.vs.log.nb.indiv[,2]))*100




##### Test for haw many replicate the R2 is better with the best allometric scaling coefficient fiited with a GLM.NB model compared to the initial Pierre Veron's Approach
length(which((BestModel.PerFilter.Nb_read.glm.nbs.DF[,5] - PerFilter.rel_read.vs.log.nb.indiv[,2]) > 0))/30 # 43%

(mean(BestModel.PerFilter.Nb_read.glm.nbs.DF[,5]) - mean(PerFilter.rel_read.vs.log.nb.indiv[,2]))*100 # -2.85% less with 





