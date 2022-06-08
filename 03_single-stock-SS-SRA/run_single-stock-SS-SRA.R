# ------------------------------------------------------------------------------------- #
# run_single-stock-SS-SRA.R
#
# Fit single stock state-space spawner-recruitment models
# ------------------------------------------------------------------------------------- #
input_data <- readRDS("./02_run-reconstruction/rr_outputs/Yukon_data_for_SS_SRA.RDS")
age_and_harvest <- read.csv("./01_inputs/data/jtc_age_and_harvest_data.csv")

stock <- c("LowerMainstem",
           "WhiteDonjek",
           "Pelly",
           "Stewart",
           "Carmacks",
           "Teslin",
           "MiddleMainstem",
           "UpperLakesAndMainstem") 
           

# prepare age data (age)
age <- input_data$x_tas_obs_agg/100
age$year <- row.names(age)
age <- age[,c(5,1:4)]
rownames(age) <- c()
age$truecount <- 1

# individual ages
age_comps <- input_data$x_tas_obs
agg_age <- as.matrix(input_data$x_tas_obs_agg)
age_ess<- matrix(100,35,8)

for(i in 1:8){age_ess[24:34,i] <- rowSums(age_comps[24:34,,i])
              age_comps[,,i] <- age_comps[,,i]/(rowSums(age_comps[,,i]))
              age_comps[1:23,,i] <-  agg_age[1:23,]/100    
              age_comps[35,,i] <-  agg_age[35,]/100 
              
    }

dimnames(age_comps) <- list(seq(1985,2019),c("a4","a5","a6","a7"),stock)

dimnames(age_ess) <- list(seq(1985,2019),stock)

# prepare spawner escapement data (esc)
population <- rep(stock,each=35,times=1)
true <- rep(1,each=35,times=8)

esc  <- cbind(input_data$S_obs, population, (input_data$tau_S_obs/input_data$S_obs), true)
esc <- as.data.frame(esc)
colnames(esc)<-c("spawn","population","CV","trueCount")
esc$spawn <- as.numeric(as.character(esc$spawn))
esc$CV <- as.numeric(as.character(esc$CV))
esc$trueCount <- as.numeric(as.character(esc$trueCount))

# prepare harvest data (harv)
er2 <- age_and_harvest %>%
  filter(year %in% c(1985:2019)) %>%
  select(er) %>%
  as.vector()
er2 <- as.vector(er2)
er <- rbind(er2, er2, er2, er2, er2, er2, er2, er2)
harv <- cbind(esc, er)
harv$harv <- (harv$spawn/(1-harv$er))*harv$er
harv <- harv[,c(2,4,6)]

# array to store posteriors
populations <- unique(esc$population)

posteriors <- array(NA,
                    dim = c(10000, 562, 8),
                    dimnames = list(NULL,NULL,populations))

# loop through populations and fit SS-SRA to each
for(i in populations){
  out <- jags_fit(i, 
                  age = age, 
                  esc = esc, 
                  har = harv)
  
  out.mcmc <- as.mcmc(out)
  posteriors[,,i] <- as.matrix(out.mcmc, chain=FALSE)
  
}

# 3-D array for main text figures
# dimensions are: posterior samples X parameters X population
saveRDS(posteriors,"./01_inputs/posteriors/single-pop-SS-SRA-posteriors.RDS")  


# 2-D matrix for supplement B
supp_B_post <- cbind(c(seq(1,2000),seq(1,2000),seq(1,2000),seq(1,2000),seq(1,2000)),
                     c(rep(1,2000),rep(2,2000),rep(3,2000),rep(4,2000),rep(5,2000)), 
                     posteriors[,,1],
                     posteriors[,,2],
                     posteriors[,,3],
                     posteriors[,,4],
                     posteriors[,,5],
                     posteriors[,,6],
                     posteriors[,,7],
                     posteriors[,,8])

colnames(supp_B_post) <-c("ITER", "CHAIN",pop_1,pop_2,pop_3,pop_4,pop_5,pop_6,pop_7,pop_8)

saveRDS(supp_B_post,"./05_integrated-SS-SRA/posterior-samples/seperated-posterior.RDS")  
