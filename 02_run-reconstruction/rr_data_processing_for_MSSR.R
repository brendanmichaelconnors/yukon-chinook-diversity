#------------------------------------------------------------------------------#
# rr_data_processing_for_MSSR.R
#
# Wrangle run-reconstruction data into format required for multi-stock SS 
# spawner-recruit modelling
#------------------------------------------------------------------------------#

rr_modAll <- read.csv("./02_run-reconstruction/rr_outputs/ensemble.csv") # ensemble run-reconstruction output
pop_age_comp <- readRDS("./01_inputs/data/pop_age_comp_2008_2019.rds") # age composition by population from merging GSI and ASL data
agg_data <- read.csv("./01_inputs/data/jtc_age_and_harvest_data.csv") # aggregate JTC age comp data and CDN harvest rate

stockIds <- read.csv("02_run-reconstruction/data/stockIDs.csv")
stockIds$stock <- c("LowerMainstem",
                    "WhiteDonjek",
                    "Pelly",
                    "Stewart",
                    "Carmacks",
                    "Teslin",
                    "MiddleMainstem",
                    "UpperLakesAndMainstem")  

rr_modAll$sd <- rr_modAll[,8] # approximate SD of border passage

xx <- merge(rr_modAll,stockIds,by="stock")

RR_output <- arrange(xx,id_MSSR,year)

agg_data <- subset(agg_data, year>1984 & year <2020)

CDN_u <- agg_data$yk_rv_har/(agg_data$spwn+agg_data$yk_rv_har)# CDN exploitation rate

# --- set values for list --------------------------------------------------------------
  ns <- 8
  
  nt <- 35
  
  na <- 4
  
  a_max <- 7
  
  C_tot_t_obs <- agg_data$total_har
    
  tau_C_obs <- round(C_tot_t_obs*0.10) 
    
  v <- rep(1,8)
    
  S_obs <- round(RR_output$X50.*(1-rep(CDN_u,8)),0)
  
  S_obs_t <- rep(1:nt,8)
  
  S_obs_s <- sort(rep(1:8,nt))
  
  S_obs_n <- length(S_obs_s)
  
  tau_S_obs <- round((RR_output$sd),0)
  
  x_tas_obs_agg <- round(agg_data[,15:18]*100); colnames(x_tas_obs_agg) <- c("a4", "a5", "a6", "a7"); rownames(x_tas_obs_agg) <- seq(1985,2019)
  
  x_tas_obs <- array(NA, dim=c(35,4,8)); x_tas_obs [24:34,,] <- pop_age_comp[]; colnames(x_tas_obs) <- c("a4", "a5", "a6", "a7"); rownames(x_tas_obs) <- seq(1985,2019)
						
  ESS_ts_agg <- matrix(NA,35,1); ESS_ts_agg[,1] <- (rowSums(x_tas_obs_agg)); rownames(ESS_ts_agg) <- seq(1985,2019); colnames(ESS_ts_agg) <- c("aggregate");  ESS_ts_agg <- as.data.frame(ESS_ts_agg)

  ESS_ts <- matrix(NA,11,8); for(i in 1:8){ESS_ts[,i] <- (rowSums(x_tas_obs[24:34,,i]))}; rownames(ESS_ts) <- seq(2008,2018); ESS_ts <- as.data.frame(ESS_ts)

  R_wish <- matrix(0,ns,ns)
  
  diag(R_wish) <- rep(1,8)
  
  df_wish <- 9
 
# --- create list ----------------------------------------------------------------------------------
Yukon_chinook_data_for_BS_ScenAA.RR <- list("ns" = ns, 
                                "nt" = nt, 
                                "na" = na, 
                                "a_max" = a_max,
                                "C_tot_t_obs" = C_tot_t_obs,
                                "tau_C_obs" = tau_C_obs,
                                "v" = v,
                                "S_obs" = S_obs, 
                                "S_obs_t" = S_obs_t,
                                "S_obs_s" = S_obs_s,
                                "S_obs_n" = S_obs_n,
                                "tau_S_obs" = tau_S_obs,
                                "x_tas_obs" = x_tas_obs,
                                "ESS_ts" = ESS_ts,
                                "x_tas_obs_agg" = x_tas_obs_agg,
                                "ESS_ts_agg" = ESS_ts_agg,
                                "R_wish" = R_wish, 
                                "df_wish" = df_wish) 

saveRDS(Yukon_chinook_data_for_BS_ScenAA.RR,"./02_run-reconstruction/rr_outputs/Yukon_data_for_SS_SRA.RDS")

# --- create new outputs for MS-SRA --------------------------------------------------------------------
esc_data <- matrix(NA,280,8)
esc_data <- as.data.frame(esc_data)
colnames(esc_data) <- c("stock","year","mean","sd","lwr95CI","upr95CI","cv","obs")
esc_data[,1] <- c(rep("LwrMain",35),
                  rep("White-Donjek",35),
                  rep("Pelly",35),
                  rep("Stewart",35),
                  rep("Carmacks",35),
                  rep("Teslin",35),
                  rep("MidMain",35),
                  rep("UprLksMain",35))

esc_data[,2] <- rep(1985:2019,times=8)
esc_data[,3] <- S_obs
esc_data[,4] <- tau_S_obs
esc_data[,5] <- round(RR_output$X2.5.*(1-rep(CDN_u,8)),0)
esc_data[,6] <- round(RR_output$X97.5.*(1-rep(CDN_u,8)),0)
esc_data[,7] <- tau_S_obs/S_obs
esc_data[,8] <- rep(1,280)

write_csv(esc_data,"./05_integrated-SS-SRA/fit-integrated-model/inputs/esc-data.csv")  


age_data <- matrix(NA,315,7)
age_data <- as.data.frame(age_data)
colnames(age_data) <- c("stock","year","a4","a5","a6","a7","n")

age_data[,1] <- c(rep("aggregate",35),
                  rep("LwrMain",35),
                  rep("White-Donjek",35),
                  rep("Pelly",35),
                  rep("Stewart",35),
                  rep("Carmacks",35),
                  rep("Teslin",35),
                  rep("MidMain",35),
                  rep("UprLksMain",35))

age_data[,3:6] <- rbind(x_tas_obs_agg,
                              x_tas_obs[,,1],
                              x_tas_obs[,,2],
                              x_tas_obs[,,3],
                              x_tas_obs[,,4],
                              x_tas_obs[,,5],
                              x_tas_obs[,,6],
                              x_tas_obs[,,7],
                              x_tas_obs[,,8])

age_data[,7] <- rowSums(age_data[,3:6])

write_csv(age_data,"./05_integrated-SS-SRA/fit-integrated-model/inputs/age-data.csv")  


catch_data <- matrix(NA,35,7)
catch_data <- as.data.frame(catch_data)
colnames(catch_data) <- c("year","mean","sd","median","lwr95CI","upr95CI","cv")
catch_data[,1] <- rep(1985:2019)
catch_data[,2] <- C_tot_t_obs
catch_data[,3] <- round(C_tot_t_obs*0.10)
catch_data[,4] <- C_tot_t_obs
catch_data[,5] <- C_tot_t_obs-(2*catch_data[,3])
catch_data[,6] <- C_tot_t_obs+(2*catch_data[,3])
catch_data[,7] <- catch_data[,3]/C_tot_t_obs

write_csv(catch_data,"./05_integrated-SS-SRA/fit-integrated-model/inputs/catch-data.csv")  

