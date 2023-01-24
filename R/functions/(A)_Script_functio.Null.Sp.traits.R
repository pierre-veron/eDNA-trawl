####NULL MODELS SES FONCTIO 

###Set Working directory
setwd("/home/romane/Documents/STAGE_PFE_IFREMER")

###library
lib_vect <- c("readr","mFD","stringr","parallel","worrms","patchwork","reshape2", "hillR")
sapply(lib_vect,library,character.only=TRUE)


#fonction source
source ("fonctions_source/Functional.Beta.MPFD.MNND.R")
source ("fonctions_source/get_q_functio.R")


#data importation
tab_eDNA<- read.csv(file="Excel_CSV/Resultats DE200249_IFREMER_Poissons_Quantif_TRIEE.csv", header=F) #eDNA SPYGEN table that has already been sorted


sp_tr_all_sp <- read_csv("Excel_CSV/Tab_sp_tr_ALL_SPECIES_csv.csv", #traits table already sorted and standardized
                         col_types = cols(Average.Depth = col_number(), 
                                          Body_shape = col_factor(levels = c()), 
                                          Depth.Range = col_number(), 
                                          Env_2 = col_factor(levels = c("pelagic","bathypelagic","demersal","bathydemersal"), ordered=TRUE), 
                                          Max_length = col_factor(levels = c("6-26cm", "26-79cm", "79-100cm", "100-459cm"), ordered=TRUE), 
                                          Repro = col_factor(levels = c()), 
                                          Trophic_Level = col_number())) 
#### Data preparation for the pipeline 
cat("########   Modification and division of eDNA table to get a table of copies ########","\n")

#(entry:tab_eDNA)

#modif eDNA table
spyID <- tab_eDNA[2,7:ncol(tab_eDNA)]
site <- tab_eDNA[3,7:ncol(tab_eDNA)]
colnames(tab_eDNA)[7:ncol(tab_eDNA)] <- str_c(spyID, paste="_", site)
colnames(tab_eDNA)[1:4]<-tab_eDNA[4,1:4]

tab_eDNA[5:6]<-NULL #delete the DB and the pres/abs 

#subdivision of eDNA table in a tables only containing the number of copies : 
tab_nb_copies <- tab_eDNA[,tab_eDNA[4,] %in% "Nb Copies"]
tab_nb_copies <- cbind(tab_eDNA[,1:4],tab_nb_copies)
tab_nb_copies <- tab_nb_copies[-c(1,2,3,4),]  #delete all the useless rows
tab_nb_copies[5:ncol(tab_nb_copies)]<-apply(tab_nb_copies[5:ncol(tab_nb_copies)],2,as.numeric)


# Classification computation
cat("######## Get the taxonomic classification of each species  ########", "\n")
#matching each species with its classification
species_list <- sp_tr_all_sp$Species

Classif <- lapply(1: length(species_list), function(x){
  error_message<- "Erreur : (204) No Content"
  tryCatch({
    cat(x,":",species_list[x],"\n")
    data.frame(wm_classification(wm_name2id(species_list[x])))
  },error=function(error_message) { return(NA) })
})

Classif_end <- do.call(rbind,lapply(Classif,function(w){
  if(is.na(w)){
    return(rep(NA,4))
  } else {
    yop <- sapply(c("Class","Order","Family", "Genus"),function(x) grep(x,w$rank))
    return(w$scientificname[yop])
  }}))

### Complete the NA
Genus <- do.call(rbind,strsplit(species_list[which(is.na(Classif_end[,1]))]," "))[,1]
Cl_Na <- do.call(rbind,lapply(Genus,function(w){Classif_end[grep(w,Classif_end[,4])[1],]}))
Classif_end[which(is.na(Classif_end[,1])),] <- Cl_Na

colnames(Classif_end) <- c("Class","Order","Family", "Genus")
Classif_end<-as.data.frame(Classif_end)
#sp_tr_all_classif <- cbind(Classif_end,sp_tr_all_sp)


### part of PIPELINE APPLICATION
nb_rd=100
tab_eDNA=tab_nb_copies
Classif=Classif_end

##################################################################
###FUNCTION 2 : Classification of species and 100 Randomisation of species belonging to the same genus/family that have been detetected at these 2 levels by spygen
## (entry : sp_tr_all_sp, tab_eDNA) 

sp_tr_all_classif <- na.omit(cbind(Classif,sp_tr_all_sp))
##searching for fish detected at the FAMILY level by spygen in the traits table and randomisation 100 times
sp_tr_fam <- sp_tr_all_classif[which(sp_tr_all_classif$Family %in% tab_eDNA$Taxon),]
fam_list <- unique(sp_tr_fam$Family)

list_fam_random <- list()
for (i in (1:nb_rd)){ #randomisation 100 times
  list_fam_random[[i]] <- do.call(rbind,lapply(seq(1,length(fam_list)), function(x){
    DFtempt <- sp_tr_fam[which(sp_tr_fam$Family==fam_list[x]),]
    DFtempt <- DFtempt[sample(1:nrow(DFtempt), 1, replace=FALSE), ]
    row.names (DFtempt) <- fam_list[x]; DFtempt 
  }))
} #end of i

##searching for fish detected at GENUS level detectetd by Spygen and randomisation 
sp_tr_genus <- sp_tr_all_classif[which(sp_tr_all_classif$Genus %in%tab_eDNA$Taxon & !sp_tr_all_classif$Species %in% tab_eDNA$Taxon ),] #condition "&" to avoid randomisation of fish already detetected at the species level
genus_list <- unique(sp_tr_genus$Genus)

list_genus_random<-list()
for (i in (1:nb_rd)){ #randomisation
  list_genus_random[[i]]<- do.call(rbind,lapply(seq(1,length(genus_list)), function(x){
    DFtempt <- sp_tr_genus[which(sp_tr_genus$Genus==genus_list[x]),] 
    DFtempt <-  DFtempt[sample(1:nrow(DFtempt), 1, replace=FALSE), ]
    row.names (DFtempt) <- genus_list[x]; DFtempt 
  }))
}

cat("######## Randomisation  ########", "\n")
#Grouping the lists of randomisation of family and genus and adding the rest of fish already detected at the species level by spygen
sp_tr_species_only <- sp_tr_all_classif[which(sp_tr_all_classif$Species %in% tab_eDNA$Taxon),]
row.names(sp_tr_species_only) <- sp_tr_species_only$Species #put species as rownames (to match with tab_eDNA in function 3) 

list_sp_tr_final <- lapply(seq(1,length(list_fam_random)), function(x){
  Test <- rbind(list_fam_random[[x]],list_genus_random[[x]], sp_tr_species_only)
  Test <- Test[order(rownames(Test)),]
  subset(Test, select= c(-Class, -Order, -Family, -Genus, -Species))
})

#(exit : list_sp_tr_final)



###################################################
#RANDOMISATION 10 times list_sp_tr_final to get 1000 traits data frame as null model

Sp_tr.NullModel = do.call(c, lapply(seq(1, length(list_sp_tr_final)), function(x){
  lapply(seq(1, 10), function(i){   #randomisation 10 times
    tr.temp = list_sp_tr_final[[x]]
    rand.sp = sample(rownames(tr.temp))
    rownames(tr.temp) = rand.sp
    tr.temp
  })
}))
class(Sp_tr.NullModel)
length(Sp_tr.NullModel)

#WARNING : DON'T SAVE SEVERAL NULL MODEL BECAUSE EACH TIME THE RANDOMISATION IS DIFFERENT 
#saveRDS(Sp_tr.NullModel, file="fichiers_RDS/Sp_tr.NullModel")
saveRDS(Sp_tr.NullModel, file="fichiers_RDS/Sp_tr.NullModel2")


