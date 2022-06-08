# ------------------------------------------------------------------------------------- #
# trib_vs_RR.R
# 
# Comparison of multi-population run-reconstruction estimates of CDN-origin Yukon Chinook 
# escapement at watershed scale and tributary level indices of escapement for individual assessment 
# projects
#---------------------------------------------------------------------------------------#

# load "data", manipulate, merge and standardize
indices <- read.csv("./02_run-reconstruction/data/DATA_InputFiles_HM_Canadian_Model_Data_2017.csv") # tributary indices from H. Hamazaki for JTC IMEG EG RR
rr <- readRDS("./02_run-reconstruction/rr_outputs/Yukon_data_for_SS_SRA.RDS") # run-reconstruction output

# Order of populations in rr 
# Lower Mainstem
# White-Donjek
# Pelly
# Stewart
# Carmacks
# Teslin
# Middle Mainstem
# Upper Lakes and Mainstem

indices <- indices[,c(1,3:13)] # subset for just indices 
long_indices <- gather(indices, key="system", value = "index",-1) # indices to long format

# assign a "stock" or population/drainage to each index
esc_indices <- long_indices %>%
  filter(Year > 1984) %>%
  mutate(stock = case_when(system == "a.Tincup" ~ "White-Donjek",
                          system == "f.Tatchun" ~ "Carmacks",
                          system == "a.Little.Sal" ~ "Carmacks",
                          system == "a.Nisutlin" ~ "Teslin",
                          system == "a.Ross" ~ "Pelly",
                          system == "a.Wolf" ~ "Teslin",
                          system == "w.Blind" ~ "Pelly",
                          system == "s.whitehorse" ~ "Upper Lakes and Mainstem",
                          system == "s.Big.Sal" ~ "Carmacks",
                          system == "a.Big.Sal" ~ "Carmacks",
                          system == "s.Teslin" ~ "Teslin"),
         tributary = case_when(system == "a.Tincup" ~ "Tincup (aerial)",
                          system == "f.Tatchun" ~ "Tatchun (foot)",
                          system == "a.Little.Sal" ~ "Little Salmon (aerial)",
                          system == "a.Nisutlin" ~ "Nisutlin (aerial)",
                          system == "a.Ross" ~ "Ross (aerial)",
                          system == "a.Wolf" ~ "Wolf (aerial)",
                          system == "w.Blind" ~ "Blind (weir)",
                          system == "s.whitehorse" ~ "Whitehorse (fishway)",
                          system == "s.Big.Sal" ~ "Big Salmon (sonar)",
                          system == "a.Big.Sal" ~ "Big Salmon (aerial)",
                          system == "s.Teslin" ~ "Teslin (sonar)")
                          )

# combine with run-reconstruction output
ms_rr <- cbind(rr$S_obs,rep(1985:2019,8),c(rep("Lower Mainstem",35),
                                           rep("White-Donjek",35),
                                           rep("Pelly",35),
                                           rep("Stewart",35),
                                           rep("Carmacks",35),
                                           rep("Teslin",35),
                                           rep("Middle Mainstem",35),
                                           rep("Upper Lakes and Mainstem",35)))

# a little more housekeeping
colnames(ms_rr) <- c("escape","Year","stock")
ms_rr <- as_tibble(ms_rr)
ms_rr$escape <- as.numeric(as.character(ms_rr$escape))
ms_rr$Year <- as.numeric(as.character(ms_rr$Year))

# output run-reconstructions as a csv
#write.csv(ms_rr, file="CDN-Yukon-RR_02Mar2020.csv")

# merge the indices and run-reconstructions by year and stock
rr_escInd <- full_join(esc_indices,ms_rr, by = c("Year","stock"))

mean_index <- rr_escInd %>%
	group_by(tributary) %>%
	summarize(Mean=mean(index, na.rm = T))
	
scale_this <- function(x){
  (x - mean(x, na.rm=TRUE)) / sd(x, na.rm=TRUE)
}

mean_index <- mean_index[-12,]	
mean_index$Mean <- round(mean_index$Mean)
# standardize the index and escapement estimates
scaled_data <- rr_escInd %>% 
  group_by(system)%>%
  mutate(scaled_index = scale_this(log(index)),
         scaled_esc = scale_this(log(escape))) %>%
  droplevels() 

scaled_data <- subset(scaled_data, system != "NA")

# plot indices vs run-reconstructions
ggplot(scaled_data, aes(x = scaled_esc, y = scaled_index)) +
  geom_smooth(method="lm", color="grey") +
  geom_point(size=1)+ 
  xlab("Standardized run reconstruction (log)") +
  ylab("Standardized tributary index (log)") +
  coord_cartesian(xlim = c(-3,3), ylim=c(-3,3)) +
  facet_wrap(~tributary, nrow=4) + 
  theme_bw()  +
  geom_text(
  	data = mean_index,
 	mapping = aes(x = -Inf, y = -Inf, label = Mean),
  	hjust = -0.4,
  	vjust = -7,
  	size = 3) +
  	geom_point(aes(x = scaled_esc, y = scaled_index, color = Year), size=1) +
  theme(strip.background = element_blank())

  
ggsave("./04_figures/figures/figureS11.jpeg", height = 5, width = 5.5)
