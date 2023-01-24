### ----------------------------------------------------------------------------
#    This script exports all the calculated indices in a .csv file
### ----------------------------------------------------------------------------

# Data from separated methods --------------------------------------------------
df_t_f <- RESULTS_fct_both$Alpha[paste0(names(Name_stations), "t"), ]
colnames(df_t_f) <- paste0("t", colnames(df_t_f))
rownames(df_t_f) <- Name_stations[gsub("t", "", rownames(df_t_f))]

df_e_f <- RESULTS_fct_both$Alpha[paste0(names(Name_stations), "e"), ]
colnames(df_e_f) <- paste0("e", colnames(df_e_f))
rownames(df_e_f) <- Name_stations[gsub("e", "", rownames(df_e_f))]

df_t_p <- RESULTS_phyl_both$Alpha[paste0(names(Name_stations), "t"), ]
colnames(df_t_p) <- paste0("t", colnames(df_t_p))
rownames(df_t_p) <- Name_stations[gsub("t", "", rownames(df_t_p))]

df_e_p <- RESULTS_phyl_both$Alpha[paste0(names(Name_stations), "e"), ]
colnames(df_e_p) <- paste0("e", colnames(df_e_p))
rownames(df_e_p) <- Name_stations[gsub("e", "", rownames(df_e_p))]


# Spatial
Spat <- EDNA_metadata[which(EDNA_metadata$replicate == 1), c("station", "Latitude", "Longitude",
                                                "max_depth")]
Spat <- as.data.frame(Spat)
rownames(Spat) <- Name_stations[Spat$station]

df <- cbind(Spat, df_e_f, df_t_f, df_e_p, df_t_p)


Alpha <- df


Alpha <- Alpha[, c("station", "Latitude", "Longitude", "max_depth","eMean.sp_richn", 
                   "ePD.mean", "ePD.sd", "eVPD.mean", "eVPD.sd", "eMPD.mean", 
                   "eMPD.sd", "eFD_q0", "eFD_q0.1","eMean.feve", "eSd.feve",
                   "eMean.fdiv", "eSd.fdiv","eMean.fric", "eSd.fric","tMean.sp_richn", 
                   "tPD.mean", "tPD.sd", "tVPD.mean", "tVPD.sd", "tMPD.mean", 
                   "tMPD.sd", "tFD_q0", "tFD_q0.1","tMean.feve", "tSd.feve",
                   "tMean.fdiv", "tSd.fdiv","tMean.fric", "tSd.fric")]

colnames(Alpha) <- c("station", "Latitude", "Longitude", "Depth", "eSR.mean", 
                     "ePD.mean", "ePD.sd", "eVPD.mean", "eVPD.sd", "eMPD.mean", 
                     "eMPD.sd", "eFD.mean", "eFD.sd","eFEve.mean", "eFEve.sd",
                     "eFDiv.mean", "eFDiv.sd","eFRic.mean", "eFRic.sd","tSR.mean", 
                     "tPD.mean", "tPD.sd", "tVPD.mean", "tVPD.sd", "tMPD.mean", 
                     "tMPD.sd", "tFD.mean", "tFD.sd","tFEve.mean", "tFeve.sd",
                     "tFDiv.mean", "tFDiv.sd","tFRic.mean", "tFRic.sd")



write.csv2(Alpha, "diversity_indices/alpha.csv")


Alpha_eDNA_SES <- as.data.frame(cbind(RESULTS_SES_phyl_both[paste0(Alpha$station, 'e'), c("PD.mean", "VPD.mean", "MPD.mean")],
                                      RESULTS_SES_fct_both[paste0(Alpha$station, 'e'), c("FD_q0", "Mean.feve", "Mean.fdiv")]))
colnames(Alpha_eDNA_SES) <- paste0("e", c("PD", "VPD", "MPD","FD", "FEve", "FDiv"), ".ses")
rownames(Alpha_eDNA_SES) <- rownames(Alpha)


Alpha_trawl_SES <- as.data.frame(cbind(RESULTS_SES_phyl_both[paste0(Alpha$station, 't'), c("PD.mean", "VPD.mean", "MPD.mean")],
                                       RESULTS_SES_fct_both[paste0(Alpha$station, 't'), c("FD_q0", "Mean.feve", "Mean.fdiv")]))
colnames(Alpha_trawl_SES) <- paste0("t", c("PD", "VPD", "MPD","FD", "FEve", "FDiv"), ".ses")
rownames(Alpha_trawl_SES) <- rownames(Alpha)

Alpha_SES <- cbind(Alpha_eDNA_SES, Alpha_trawl_SES, Alpha)

write.csv2(Alpha_SES, "diversity_indices/alpha_ses.csv")
rm(df, df_e_f, df_e_p, df_t_f, df_t_p, Alpha_eDNA_SES, Alpha_trawl_SES)

# Data from merged method ------------------------------------------------------
df_m_f <- RESULTS_fct_merged$Alpha[names(Name_stations),]
df_m_p <- RESULTS_phyl_merged$Alpha[names(Name_stations),]

df <- cbind(Spat, df_m_f, df_m_p)



Alpha_merged <- df





Alpha_merged <- Alpha_merged[, c("station", "Latitude", "Longitude", "max_depth", "Mean.sp_richn", 
                                 "PD.mean", "PD.sd", "VPD.mean", "VPD.sd", "MPD.mean", 
                                 "MPD.sd", "FD_q0", "FD_q0.1","Mean.feve", "Sd.feve",
                                 "Mean.fdiv", "Sd.fdiv","Mean.fric", "Sd.fric")]

colnames(Alpha_merged) <- c("station", "Latitude", "Longitude", "Depth", "SR.mean", 
                            "PD.mean", "PD.sd", "VPD.mean", "VPD.sd", "MPD.mean", 
                            "MPD.sd", "FD.mean", "FD.sd","FEve.mean", "FEve.sd",
                            "FDiv.mean", "FDiv.sd","FRic.mean", "FRic.sd")
write.csv2(Alpha, "diversity_indices/alpha_merged.csv")

Alpha_merged_SES <- as.data.frame(cbind(RESULTS_SES_phyl_merged[Alpha_merged$station, c("PD.mean", "VPD.mean", "MPD.mean")],
                                        RESULTS_SES_fct_merged[Alpha_merged$station, c("FD_q0", "Mean.feve", "Mean.fdiv")]))
colnames(Alpha_merged_SES) <- paste0( c("PD", "VPD", "MPD","FD", "FEve", "FDiv"), ".ses")
rownames(Alpha_merged_SES) <- rownames(Alpha_merged)


Alpha_merged_SES <- cbind(Alpha_merged_SES, Alpha_merged)
write.csv2(Alpha_merged_SES, "diversity_indices/alpha_merged_ses.csv")


rm(Spat, df, df_m_f, df_m_p)
