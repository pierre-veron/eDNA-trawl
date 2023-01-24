### ----------------------------------------------------------------------------
#    Load used packages
### ----------------------------------------------------------------------------


lib_vect = c("tidyverse", "here", "gapminder", "rmarkdown", "taxize", "stringr",
             "data.tree", "treemap", "readxl", "geosphere", "ape", "cluster",
             "wesanderson", "Hmisc", "hillR","raster","rgeos","rgdal","sp","maptools",
             "shape","ade4","picante","phytools","parallel","RCurl","XML",
             "taxize","PhyloMeasures","fishtree","plotrix",'sf',"gtools",
             "mFD", "betapart", "usedist", "caper", "readr", "phangorn", "png", "unikn",
             "dismo","parallel","worrms","patchwork",
             "reshape2", "matrixStats", "fitdistrplus", "latex2exp")

sapply(lib_vect,library,character.only=TRUE)
select <- dplyr::select
rm(lib_vect)
