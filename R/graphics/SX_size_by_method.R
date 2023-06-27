### ----------------------------------------------------------------------------
#   This script draws the distributions of fish size by method
#
#   Author : Pierre Veron
#   pierre.veron.2017@polytechnique.org
#   Created 2022/04/25
### ----------------------------------------------------------------------------

# Defining proper categories for the traits to be plotted ----------------------
Temp_traits = Traits
Temp_traits$Max_length <- factor(sapply(Traits$Max_length, function(x) {
  if (x <25) {
    "0 - 25"
  } else if (x < 50) {
    "25 - 50"
  } else if (x < 75) {
    "50 - 75"
  } else if (x < 100) {
    "75 - 100"
  } else if (x < 500) {
    "100 - 500"
  } else {
    "err"
  }
}), levels = c("0 - 25","25 - 50","50 - 75","75 - 100","100 - 500"))


# counting the number of sites where each taxa is detected
temp_reftax <- Reftax
temp_reftax$site_edna <- as.vector(sapply(temp_reftax$taxname, function(tax) {
  sum(Incidence_matrix_eDNA[, gsub(" ", ".", tax)])
}))

temp_reftax$site_trawl <- as.vector(sapply(temp_reftax$taxname, function(tax) {
  sum(Incidence_matrix_trawl[, gsub(" ", ".", tax)])
}))

# General graphical parameters -------------------------------------------------
out_file <- "figures/Fig_SX_fish_size_by_method"
ext <- ".pdf"

image_width <- 17
image_height <- 15
res = 800

pdf(file = paste0(out_file, ext),
    width = image_width / 2.54,
    height = image_height / 2.54)

par(mfrow = c(1,1), oma = c(0,2,1,0), mar = c(3.2,1,0.5,1))

trait <- "Max_length"
if (!is.factor(Temp_traits[, trait])) {
  classes <- quantile(Temp_traits[, trait], seq(0,1,length.out = n_class+1))
  classes[1] <- 0.99*classes[1]
  classes[length(classes)] <- 1.01*classes[length(classes)] # to include all values with < and > 
  names(classes) <- NULL
  classes
}




# Random and plotting ----------------------------------------------------------
n_rd_draw <- 100
List_rd_distri <- lapply(1:n_rd_draw, function(i) {
  
  Rd_reftax <- temp_reftax
  
  Genus <- Rd_reftax[which(Rd_reftax$rank == "genus"), "taxname"]
  Species <- as.vector(sapply(Genus, function(gen) {
    sample(Reftax_species[which(Reftax_species$genus == gen), "species"],1)
  }))
  Rd_reftax[which(Rd_reftax$rank == "genus"), "taxname"] <- Species
  
  # Calculate distributions 
  if (is.factor(Temp_traits[, trait])) { #  for factors
    lev <- levels(Temp_traits[, trait])
    val_traits <- as.numeric(Temp_traits[Rd_reftax$taxname, trait])
    val_traits_edna <- rep(val_traits, Rd_reftax$site_edna)
    distri_edna <- sapply(1:length(lev), function(j) {
      length(which(val_traits_edna == j))
    }) 
    
    val_traits_trawl <- rep(val_traits, Rd_reftax$site_trawl)
    distri_trawl <- sapply(1:length(lev), function(j) {
      length(which(val_traits_trawl == j))
    })
    
    dsty_edna <- NULL
    dsty_trawl <- NULL
    
  } else { # for numeric values 
    val_traits <- Temp_traits[Rd_reftax$taxname, trait]
    val_traits_edna <- rep(val_traits, Rd_reftax$site_edna)
    dsty_edna <- density(val_traits_edna)
    val_traits_trawl <- rep(val_traits, Rd_reftax$site_trawl)
    dsty_trawl <- density(val_traits_trawl)
    
    distri_edna <- sapply(1:n_class, function(i_cl) {
      length(which(val_traits_edna >= classes[i_cl] & val_traits_edna < classes[i_cl+1]))
    })
    
    distri_trawl <- sapply(1:n_class, function(i_cl) {
      length(which(val_traits_trawl >= classes[i_cl] & val_traits_trawl < classes[i_cl+1]))
    })
    
  }
  list("eDNA" = distri_edna, "trawl" = distri_trawl, "dsty.eDNA" = dsty_edna, "dsty.trawl" = dsty_trawl)
})

# Calculate mean and std of densities
Mat_distri_eDNA <- sapply(List_rd_distri, function(rd_distri) {
  rd_distri$eDNA
})

av_distri_eDNA <- rowMeans(Mat_distri_eDNA)
sd_distri_eDNA <- rowSds(Mat_distri_eDNA)

Mat_distri_trawl <- sapply(List_rd_distri, function(rd_distri) {
  rd_distri$trawl
})

av_distri_trawl <- rowMeans(Mat_distri_trawl)
sd_distri_trawl <- rowSds(Mat_distri_trawl)

# Chi square test
chitest <- chisq.test(t(data.frame("eDNA" = av_distri_eDNA,
                                   "trawl"= av_distri_trawl)))


# Plotting -------------------------------------------------------------------

lev <- levels(Temp_traits[, trait])

dsty_edna <- av_distri_eDNA / sum(av_distri_eDNA)
sd_dsty_eDNA <- sd_distri_eDNA / sum(av_distri_eDNA)

dsty_trawl <- av_distri_trawl / sum(av_distri_trawl)
sd_dsty_trawl <- sd_distri_trawl / sum(av_distri_trawl)

xaxislength = length(lev)+0.2
plot(c(), c(), xlim = c(1 - 0.4 + 0.02*xaxislength,length(lev)+0.4 - 0.02*xaxislength),
     ylim = 1 * c(0,0.8), 
     xaxt = 'n',
     yaxt = 'n',
     xlab = "",
     ylab = "")

axis(side = 1, at = 1:length(lev), line = 0, tck = -0.01, labels = F)
axis(side = 1, at = 1:length(lev), line = -0.6, tick = F, lty = 0, labels = lev, cex.axis = 1)

axis(side = 2, at = seq(-1,1, 0.2), line = 0, tck = -0.01, labels = F)

mtext("Frequency", side = 2, line = 1.75)
axis(side = 2, at = seq(-1,1, 0.2), line = -0.6, tick = F, lty = 0, labels = T, las = 2, cex.axis = 1)
# grid
for (a in  seq(-1,1, 0.2)) {
  abline(a = a, b = 0, col = "lightgray", lty = 1, lwd  = 1)
}
abline(a = 0, b = 0, col = "black", lty = 1, lwd  = 1)



rect(xleft=1:length(lev) - 0.4, ybottom=0, xright=1:length(lev) , ytop = dsty_edna, 
     col = color_edna, border = NA)
rect(xleft=1:length(lev) , ybottom=0, xright=1:length(lev) + 0.4, ytop = dsty_trawl, 
     col = color_trawl, border = NA)
errbar(x = 1:length(lev) - 0.2, y = dsty_edna, yplus =  dsty_edna+sd_dsty_eDNA, yminus =  dsty_edna-sd_dsty_eDNA, add = T, pch = NA)
errbar(x = 1:length(lev) + 0.2, y = dsty_trawl, yplus = dsty_trawl+sd_dsty_trawl, yminus = dsty_trawl-sd_dsty_trawl, add = T, pch = NA)
mtext("Fish maximum size (cm)", side = 1, line = 2)


legend("topright", inset = c(0.02,0.02), fill = c(color_edna, color_trawl), legend = c("eDNA", "trawling"), bg = "white", cex = 1.6)
dev.off()
