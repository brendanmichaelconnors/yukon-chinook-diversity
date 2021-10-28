# ------------------------------------------------------------------------------------- #
# make_figures.R
#
# Manuscript figures
# ------------------------------------------------------------------------------------- #

# load data/outputs ----
border_passage <- read.csv("01_inputs/data/borderCounts.csv") # border passage
gsi <- read.csv("01_inputs/data/20Oct2021-update.gsiSamplesAllProbs.csv") # GSI samples
ASL <- read.csv("01_inputs/data/ASL_GSI_modified.csv") # ASL data
jtc <- read.csv("01_inputs/data/jtc_age_and_harvest_data.csv") # aggregate JTC age comp data and CDN harvest rate
ensemble <- read.csv("./02_run-reconstruction/rr_outputs/ensemble.csv") # ensemble run-reconstruction output
RR <-load("./02_run-reconstruction/rr_outputs/rpt.fullCor.Rdata") # fully estimated var-covar run-reconstruction model output
SSSR <- readRDS("./01_inputs/posteriors/single-pop-SS-SRA-posteriors.RDS") # posterior samples from SS-SRA models

# pre-process data ----
pop_colors <- viridis(8)
pop_colors_ord <- pop_colors[c(1,4,2,3,5,6,7,8)]

ensemble <- ensemble %>%
  rename(lwr_2.5 = X2.5.,
         lwr_25 = X25.,
         med = X50.,
         upr_75 = X75.,
         upr_97.5 = X97.5.) %>%
  mutate(population = case_when(stock == "LowerMainstem" ~ "Lower Mainstem",
                                stock == "WhiteDonjek" ~ "White-Donjek",
                                stock == "Pelly" ~ "Pelly",
                                stock == "Stewart" ~ "Stewart",
                                stock == "Carmacks" ~ "Carmacks",
                                stock == "Teslin" ~ "Teslin",
                                stock == "MiddleMainstem" ~ "Middle Mainstem",
                                stock == "UpperLakesAndMainstem" ~ "Upper Lakes and Mainstem")) %>%
  select(population, year, lwr_2.5: upr_97.5)

ensemble <- ensemble %>%
  mutate(Color = case_when(population == "Lower Mainstem" ~ "#440154FF",
                           population == "Stewart" ~ "#46337EFF",
                           population == "Pelly" ~ "#365C8DFF",
                           population == "White-Donjek" ~ "#277F8EFF",
                           population == "Middle Mainstem" ~ "#1FA187FF",
                           population == "Carmacks" ~ "#4AC16DFF",
                           population == "Upper Lakes and Mainstem" ~ "#9FDA3AFF",
                           population == "Teslin" ~ "#FDE725FF"))

ensemble[,3:7] <- ensemble[,3:7]/1000

ensemble$pops_f <- factor(ensemble$population, levels = c("Lower Mainstem", "White-Donjek", "Stewart", "Pelly", "Middle Mainstem", "Carmacks", "Upper Lakes and Mainstem", "Teslin"))

# posterior wrangling
colnames(SSSR) <- pop

population <- rep(c("Lower Mainstem",
                    "White-Donjek",
                    "Pelly",
                    "Stewart",
                    "Carmacks",
                    "Teslin",
                    "Middle Mainstem",
                    "Upper Lakes and Mainstem"),each=10000)

Post.YlC <- as.data.frame(SSSR[1:10000,1:562,1])
Post.YS <- as.data.frame(SSSR[1:10000,1:562,2])
Post.YP <- as.data.frame(SSSR[1:10000,1:562,3])
Post.YWD <- as.data.frame(SSSR[1:10000,1:562,4])
Post.Ym <- as.data.frame(SSSR[1:10000,1:562,5])
Post.YC <- as.data.frame(SSSR[1:10000,1:562,6])
Post.Yu <- as.data.frame(SSSR[1:10000,1:562,7])
Post.YT <- as.data.frame(SSSR[1:10000,1:562,8])

Posteriors <- do.call("rbind",list(Post.YlC,Post.YS,Post.YP,Post.YWD,Post.Ym,Post.YC,Post.Yu,Post.YT)) 
Posteriors.df <- cbind(population,Posteriors)

# combine alpha and beta population dataframes
alpha_beta_per <- Posteriors.df %>%
  dplyr::select(population, alpha,beta,phi,sigma.R,U.msy) %>%
  mutate(equilibrium = log(alpha)/beta,
         lnalpha = log(alpha)) %>%
  group_by(population) %>%
  dplyr::summarise(alpha_low = exp(quantile(lnalpha, probs=0.05)),
                   alpha_med = exp(quantile(lnalpha, probs=0.5)),
                   alpha_upr = exp(quantile(lnalpha, probs=0.95)),
                   equil_low = quantile(equilibrium, probs=0.05),
                   equil_med = quantile(equilibrium, probs=0.5),
                   equil_upr = quantile(equilibrium, probs=0.95),
                   phi_med = quantile(phi, probs=0.5),
                   phi_low= quantile(phi, probs=0.05),
                   phi_upr = quantile(phi, probs=0.95),
                   sigma.R_med = quantile(sigma.R, probs=0.5),
                   sigma.R_low= quantile(sigma.R, probs=0.05),
                   sigma.R_upr = quantile(sigma.R, probs=0.95),
                   U.msy = quantile(U.msy, probs=0.5)) %>%
  as.data.frame()

# create resid dataframe
logresid_df <- Posteriors.df[,c(1,47:84)] %>%
  group_by(population) %>%
  sample_n(1000,replace=TRUE) %>%
  select(1,2,13,24,34:39,3:12,14:23,25:33) %>%
  as.data.frame()

logresid_perc <- logresid_df %>%
  group_by(population) %>%
  gather(key = year, value = value, "log.resid[1]":"log.resid[31]") %>%
  ungroup() %>% 
  group_by(population,year) %>%
  dplyr::summarise(lwr = quantile(value, probs = c(0.05)),
                   med = quantile(value, probs = c(0.5)),
                   upr = quantile(value, probs = c(0.95))) %>%
  ungroup() %>%
  mutate(year2 = rep(c(1,10:19,2,20:31,3:9),times=8),
         Adjyear = year2 + 6,
         BroodYear = year2 + 1981,
         Year = 1982 + Adjyear) %>%
  dplyr::select(1,6,7,8,9,3,4,5) %>%
  arrange(population,Adjyear) %>%
  as.data.frame()

logresid_perc_scale <- logresid_df %>%
  group_by(population) %>%
  gather(key = year, value = value, "log.resid[1]":"log.resid[31]") %>%
  mutate(value=(value-mean(value))/sd(value))%>%
  ungroup() %>% 
  group_by(population,year) %>%
  dplyr::summarise(lwr = quantile(value, probs = c(0.05)),
                   med = quantile(value, probs = c(0.5)),
                   upr = quantile(value, probs = c(0.95))) %>%
  ungroup() %>%
  mutate(year2 = rep(c(1,10:19,2,20:31,3:9),times=8),
         Adjyear = year2 + 6,
         BroodYear = year2 + 1981,
         Year = 1982 + Adjyear) %>%
  dplyr::select(1,6,7,8,9,3,4,5) %>%
  arrange(population,Adjyear) %>%
  as.data.frame()

# create spawner dataframe
spawn_df <- Posteriors.df[,c(1,492:526)] %>%
  group_by(population) %>%
  sample_n(1000,replace=TRUE) %>%
  as.data.frame()

spawn_df2 <- spawn_df[,c(1,2,13,24,31:36,3:12,14:23,25:30)]

spawn <- spawn_df2 %>%
  group_by(population) %>%
  gather(key=cat, value = value, c("S[1]":"S[35]")) %>%
  mutate(param = "S") %>%
  group_by(population, cat) %>%
  dplyr::mutate(lwr = quantile(value, probs = 0.05),
                med = quantile(value, probs = 0.5),
                upr = quantile(value, probs = 0.95)) %>%
  distinct(population, cat, param, lwr, med, upr) %>%
  as.data.frame()

# create recruit dataframe
rec_df <- Posteriors.df[,c(1,453:490)] %>%
  group_by(population) %>%
  sample_n(1000,replace=TRUE) %>%
  as.data.frame()

rec_df2 <- rec_df[,c(1,2,13,24,34:39,3:12,14:23,25:33)]

recruit <- rec_df2 %>%
  group_by(population) %>%
  gather(key=cat, value = value, c("R[1]":"R[38]")) %>%
  mutate(param = "R") %>%
  group_by(population, cat) %>%
  dplyr::mutate(lwr = quantile(value, probs = 0.05),
                med = quantile(value, probs = 0.5),
                upr = quantile(value, probs = 0.95)) %>%
  distinct(population, cat, param, lwr, med, upr) %>%
  as.data.frame()

# create equilibrium values and subsample 1000 values per population
alphabeta_df <-  Posteriors.df %>%
  dplyr::select(population,alpha,beta) %>%
  mutate(equilibrium = log(alpha)/beta,
         alpha = alpha) %>%
  group_by(population) %>%
  sample_n(1000,replace=TRUE) %>%
  as.data.frame()


# range of equilibrium values
alphabeta_df2_summ <- alphabeta_df %>% 
  group_by(population) %>%
  summarise(min = min(equilibrium),
            max = max(equilibrium)) %>%
  as.data.frame()

# create vector of abundances for each dataframe
sp_YC <- as.vector(seq(0,20,length.out=1000))
sp_YlC <- as.vector(seq(0,20,length.out=1000))
sp_Ym <- as.vector(seq(0,40,length.out=1000))
sp_YP <- as.vector(seq(0,25,length.out=1000))
sp_YS <- as.vector(seq(0,25,length.out=1000))
sp_YT <- as.vector(seq(0,25,length.out=1000))
sp_Yu <- as.vector(seq(0,20,length.out=1000))
sp_YWD <- as.vector(seq(0,13,length.out=1000))

# repeat substock abundance vector 1000 times for each value of alpha and beta
spw_YC <- rep(c(sp_YC),(1000))
spw_YlC <- rep(c(sp_YlC),(1000))
spw_Ym <- rep(c(sp_Ym),(1000))
spw_YP <- rep(c(sp_YP),(1000))
spw_YS <- rep(c(sp_YS),(1000))
spw_YT <- rep(c(sp_YT),(1000)) 
spw_Yu <- rep(c(sp_Yu),(1000))
spw_YWD <- rep(c(sp_YWD),(1000))

# combine abundances into single vector
spw_df <- c(spw_YC,spw_YlC,spw_Ym,spw_YP,spw_YS,spw_YT,spw_Yu,spw_YWD)

# repeat each alpha beta row 1000x
alphabeta_df3 <- alphabeta_df[rep(seq_len(nrow(alphabeta_df)),each=1000),]

# add column of abundances to alphabeta_df2 dataframe
alphabeta_df3$abund <- NA
alphabeta_df3$abund <- spw_df

# add on 'y' value that has the predicted ricker estimates per abundance value
alphabeta_df3$y <- NA
alphabeta_df3$y <- alphabeta_df3$alpha*alphabeta_df3$abund*exp(-alphabeta_df3$beta*alphabeta_df3$abund)


# add scenario simulation number
scen <- rep(1:1000, each = 1000)
alphabeta_df3$scenario <- NA
alphabeta_df3$scenario <- rep(c(scen),c(1))

# create percentile dataframe for alphabeta_df3
alphabeta_df3_perc <- alphabeta_df3 %>% 
  group_by(population, abund) %>%
  dplyr::summarise(q_05 = quantile(y, probs=0.05),
                   q_95 = quantile(y, probs=0.95),
                   q_50 = quantile(y, probs=0.5)) %>%
  as.data.frame()

# bind spawners and recruits together
RS_df <- rbind(recruit, spawn) %>%
  as.data.frame() %>%
  mutate(year = as.numeric(c(rep(c(1:38),each=8),rep(c(1:35),each=8)))) %>%
  mutate(Adjyear = case_when(param == "R" ~ year-7,
                             param == "S" ~ year)) %>%
  filter(Adjyear >= 1) %>%
  mutate(BroodYear = Adjyear + 1981) %>%
  dplyr::select(population, BroodYear, param, lwr, med, upr) %>%
  as.data.frame()

Rec_df <- RS_df %>% filter(param == "R") %>%
  dplyr::rename(R_lwr = lwr,
                R_med = med,
                R_upr = upr) %>%
  dplyr::select(-param)

Spawn_df <- RS_df %>% filter(param == "S") %>%
  dplyr::rename(S_lwr = lwr,
                S_med = med,
                S_upr = upr) %>%
  dplyr::select(-param)

RS_df2 <- full_join(Rec_df, Spawn_df)

alpha_beta_per$population <- factor(alpha_beta_per$population, levels = c("Lower Mainstem","Stewart","Pelly",
                                                                          "White-Donjek","Middle Mainstem",
                                                                          "Carmacks","Upper Lakes and Mainstem","Teslin"))

# create dummy data for figure 8 plotting
dummy <- data.frame(population = c("Lower Mainstem","White-Donjek","Stewart","Pelly","Middle Mainstem","Carmacks","Upper Lakes and Mainstem","Teslin"), 
                    BroodYear = 1981, 
                    R_med = rep(5,8),
                    R_lwr = rep(5,8),
                    R_upr = c(20,18,15,10,15,14,9,15),
                    S_lwr = rep(5,8),
                    S_med = rep(5,8),
                    S_upr = rep(5,8))

RS_df2 <- full_join(RS_df2, dummy)

alphabeta_df3_perc$population <- factor(alphabeta_df3_perc$population, levels = c("Lower Mainstem","White-Donjek","Stewart","Pelly",
                                                                                  "Middle Mainstem",
                                                                                  "Carmacks","Upper Lakes and Mainstem","Teslin"))

RS_df2$population <- factor(RS_df2$population, levels = c("Lower Mainstem","White-Donjek","Stewart","Pelly",
                                                          "Middle Mainstem",
                                                          "Carmacks","Upper Lakes and Mainstem","Teslin"))



# Figure 3: multi-panel daily border passage + corr run-timing + annual border passage ----

# daily border passage
load("./02_run-reconstruction/rr_outputs/rpt.fullCor.Rdata")
rr_mod <- rpt
geo_stk_num <- c(1,2,4,3,7,5,8,6)
avg_run_time <- matrix(NA,126,8)

for(s in geo_stk_num){
  for( d in 1: rr_mod $nD )
  {
    avg_run_time[d,s] <- mean(rr_mod $rho_dst[d,s, ])
  }
}

run_time <- cbind(seq(160,285,1),avg_run_time)
colnames(run_time) <- c("day","Lower Mainstem", "White-Donjek", "Pelly", "Stewart", "Carmacks", "Teslin", "Middle Mainstem", "Upper Mainstem")
run_time <- as.data.frame(run_time)

run_time_long <- run_time%>%
  pivot_longer(cols=2:9, names_to = "population", values_to="passage")

run_time_long$population <- factor(run_time_long$population, levels = c("Teslin","Upper Mainstem","Carmacks","Middle Mainstem",
                                                                            "Pelly", "Stewart","White-Donjek","Lower Mainstem"))
a <- ggplot(run_time_long, aes(x = day, y = population, height=passage, fill=population,color=population)) +
      geom_density_ridges(stat="identity",scale = 1.4) +
      scale_fill_manual(values = rev(c("#440154FF","#277F8EFF" ,"#46337EFF" ,"#365C8DFF"  ,"#1FA187FF","#4AC16DFF","#9FDA3AFF","#FDE725FF" )))+
      scale_color_manual(values = rev(c("#440154FF","#277F8EFF" ,"#46337EFF" ,"#365C8DFF"  ,"#1FA187FF","#4AC16DFF","#9FDA3AFF","#FDE725FF" )))+
      coord_cartesian(xlim = c(175, 265)) +
      xlab("Day of year") +
      theme_bw() +
      theme(axis.title = element_text(size=9),
          axis.text = element_text(size=6),
          axis.title.y = element_blank(),
          panel.grid.minor = element_blank(),
          legend.position = "none",
          plot.margin = unit(c(1,1,1,0), units = "lines")) 


# correlation in run-timing 
run_corr <- cov2cor(rpt$cov_ss)
run_corr_2 <- run_corr[c(6,8,5,7,3,4,2,1),c(6,8,5,7,3,4,2,1)]
colnames(run_corr_2) <- c("Teslin","Upper Mainstem","Carmacks","Middle Mainstem","Pelly","Stewart","White-Donjek","Lower Mainstem")
rownames(run_corr_2) <- c("Teslin","Upper Mainstem","Carmacks","Middle Mainstem","Pelly","Stewart","White-Donjek","Lower Mainstem")

b <- ggcorrplot(run_corr_2, type = "upper",
                outline.col = "white",
                lab = "TRUE",
                legend.title = "Correlation",
                show.legend = FALSE,
                lab_size = 2,
                tl.cex = 6)  

# border passage by population 
c <- ggplot(ensemble, aes(x = as.factor(year), y = med, fill = pops_f)) + 
          geom_bar(stat = "identity",colour = "black", width=1) +
          scale_fill_manual(values = c("#440154FF", "#277F8EFF","#46337EFF", "#365C8DFF", "#1FA187FF", "#4AC16DFF","#9FDA3AFF", "#FDE725FF")) +
          geom_errorbar(data=ensemble, aes(x= as.factor(year), ymin= lwr_2.5, ymax=upr_97.5, width=.1)) +
          xlab("Year") +
          ylab("Border passage (000s)") +
          coord_cartesian(ylim = c(0,40))+
          scale_x_discrete(breaks = c("1984","1988","1992","1996","2000","2004","2008","2012","2016")) +
          scale_y_continuous(position = "right") +
          theme_bw() +
          facet_wrap(~pops_f, ncol=2) +
          theme(strip.text = element_text(size=6),
                axis.text.x = element_text(angle=45, hjust=1),
                axis.text = element_text(size=6),
                axis.title = element_text(size=9),
                legend.position = "none",
                panel.grid.minor = element_blank())

# combine three plots
g <- ggarrange(ggarrange(a,b,nrow =2, labels = c("a", "b")),c,
               ncol = 2, labels = c("","c"), hjust=0.7)

jpeg("./04_figures/figures/figure3.jpeg", width = 6, height = 5, units = "in", res = 600)
print(g)
dev.off()

# Figure 8: Multi-panel SR and trade-offs

source("./04_figures/tradeoffs.R")

# Figure 4: spawner-recruit ----

a <- ggplot() +
  geom_ribbon(data = alphabeta_df3_perc, aes(x = abund, ymin = q_05, ymax = q_95), 
              fill = "grey80", alpha=0.5, linetype=2, colour="gray46") +
  geom_line(data = alphabeta_df3_perc, aes(x = abund, y = q_50), color="black", size = 0.75) +
  geom_errorbar(data = RS_df2, aes(x= S_med, y = R_med, ymin = R_lwr, ymax = R_upr), width=0,colour="grey42", width=0.2, size=0.3) +
  geom_errorbarh(data = RS_df2, aes(x= S_med, y = R_med, xmin = S_lwr, xmax = S_upr), 
                 height=0,colour = "grey42", width=0.2, size = 0.3) + 
  geom_point(data = RS_df2, aes(x= S_med, y = R_med, color=BroodYear, width=0.9), size=1) +
  xlab("Spawners (000s)") +
  ylab("Recruits (000s)") +
  facet_wrap(~population, ncol=2, scales="free") +
  theme_bw() +
  theme(strip.text = element_text(size=6),
        axis.title = element_text(size=9),
        axis.text = element_text(size=6),
        legend.justification = c(0,0),
        legend.position = c(0.82,0.365),
        legend.key.size = unit(4, "pt"),
        legend.background = element_blank(),
        legend.text = element_text(size = 4),
        legend.title = element_text(size = 5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

pop_colors <- viridis(8,alpha=0.5)[c(6,1,5,3,2,8,7,4)]
pop_colors_text <- viridis(8)[c(6,1,5,3,2,8,7,4)]

b <- ggplot(data = alpha_beta_per, aes(x = equil_med, y = alpha_med, fill=population,label=population)) +
  geom_errorbar(data=alpha_beta_per, mapping=aes(x=equil_med, ymin=alpha_low, ymax=alpha_upr), width=0,color = pop_colors, size =0.5) +
  geom_errorbarh(data=alpha_beta_per, mapping=aes(x=alpha_med, xmin=equil_low, xmax=equil_upr), height=0,color = pop_colors, size = 0.5) +
  geom_point(size=2, pch=21, colour="black") +
  xlab("Equilibrium size (000s)") +
  ylab("Intrinsic productivity")  +
  scale_fill_viridis(option = "viridis", discrete = TRUE) +
  coord_cartesian(xlim = c(0, 30), ylim = c(2,8)) +
  scale_y_continuous(breaks = c(2,4,6,8,10)) +
  scale_x_continuous(breaks = c(0,5,10,15,20,25,30)) +
  theme_bw() +
  theme(axis.title = element_text(size=9),
        axis.text = element_text(size=6),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none",
        plot.margin = unit(c(1,2,3,1), units = "lines")) +
  geom_text_repel(color = pop_colors_text,
                  point.padding=1, 
                  nudge_x =c(0,0,0,0,0,0,0,0), 
                  nudge_y=c(0,0,0,0,0,0,0,0),
                  size=2)

source("./04_figures/tradeoffs.R")

c <- ggplot(t3.3, aes(x=U*100, y=med_V2)) +
  geom_line(linetype = "solid", size=0.5) +
  geom_ribbon(aes(ymin=low_V2, ymax=upp_V2), alpha=0.2, lty=1) +
  geom_line(aes(x=U*100,y=over.med*100), lty=2, size=0.5) +
  geom_ribbon(aes(ymin=over.low*100, ymax=over.up*100), alpha=0.2, lty=2) +
  geom_line(aes(x=U*100,y=ext.med*100), lty=3, size=0.5) +
  geom_ribbon(aes(ymin=ext.low*100, ymax=ext.up*100), alpha=0.2, lty=3) +
  scale_y_continuous(limits = c(0,30), breaks=c(0,10,20,30,40)) +
  scale_y_continuous(sec.axis = sec_axis(~./1, name = "Populations at risk (%)"),limits = c(0,100)) +
  scale_x_continuous(limits = c(0,100), breaks=c(0,20,40,60,80,100)) +
  ylab("Harvest (000s)") +
  xlab("Harvest rate (%)") +
  theme_bw() +
  theme(axis.title = element_text(size= 9),
        axis.text = element_text(size = 7),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none",
        plot.margin = unit(c(1,0,3,0.5), units = "lines")) +
  annotate("text", x = c(9,9,9),
           y = c(95,90,85),
           label = c("Harvest", "Overfished", "Extinct"),
           color="black", 
           size=2,
           hjust=0) +
  annotate("segment", x = c(0,0,0),
           xend=c(7,7,7),
           y = c(95,90,85),
           yend = c(95,90,85),
           lty = c(1,2,3),
           color="black", 
           size=0.3)

g <- ggarrange(a,ggarrange(b,c,nrow =2, labels = c("b", "c")),
               ncol = 2, labels = c("a"))

jpeg("./04_figures/figures/figure4.jpeg", width = 6, height = 6, units = "in", res = 600)
print(g)
dev.off()

# Figure 5: multi-panel productivity correlation ----
logresid_perc_scale$population <- factor(logresid_perc_scale$population, levels = c("Lower Mainstem","White-Donjek","Stewart","Pelly",
                                                                                    "Middle Mainstem",
                                                                                    "Carmacks","Upper Lakes and Mainstem","Teslin"))

a <- ggplot(logresid_perc_scale, aes(x=BroodYear, y = med , color=population), show.legend = F) +
  geom_point(size=1,show.legend = F) + 
  scale_color_manual(values = pop_colors_ord, name=population) +
  geom_line(show.legend = F) + 
  geom_ribbon(aes(ymin = lwr, ymax = upr), show.legend = F, alpha=0.4) +
  scale_color_manual(values = pop_colors_ord) +
  geom_hline(yintercept = 0, lty = "dotted") +
  xlab("Brood year") +
  ylab("Productivity index") +
  scale_x_continuous(breaks=c(1985,1990,1995, 2000, 2005, 2010)) +
  facet_wrap(~population, ncol=2) + 
  theme(legend.position = "none") +
  theme_bw() +
  theme(axis.text.x = element_text(angle=45, hjust = 1, size=6),
        strip.text = element_text(size=6),
        axis.title = element_text(size=9),
        axis.text = element_text(size=6),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

resids <- logresid_perc[,c(1,2,7)]
cor_resids <- spread(resids, population, med)[,-1]
corr_matrix <- cor(cor_resids[,c(6,7,1,3,4,5,8,2)])

b <- ggcorrplot(cor(cor_resids[,c(6,7,1,3,4,5,8,2)]), type = "upper",
                outline.col = "white",
                lab = "TRUE",
                legend.title = "Correlation",
                show.legend = FALSE,
                lab_size = 2,
                tl.cex = 6) +
  theme(axis.title = element_text(size=3),
        axis.text = element_text(size=3))

phis <- matrix(NA,1000,8); colnames(phis)<-unique(spawn_df$population)
for (i in unique(spawn_df$population)){
  for(j in 1:1000){
    pop <- subset(Posteriors.df, population == i)[sample(1:1000,1),]
    phis[j,i] <-pop$phi
  }
}
colnames(phis)<-c(as.character(unique(spawn_df$population)))
phis.df <- as.data.frame(phis)                 
phis.long <-pivot_longer(data=phis.df,cols=1:8, names_to = "pop", values_to = "phi")
phis.long$pop <- factor(phis.long$pop, levels = c("Teslin","Upper Lakes and Mainstem","Carmacks",
                                                  "Middle Mainstem","Pelly","Stewart",
                                                  "White-Donjek","Lower Mainstem"))

c <- ggplot(data = phis.long, aes(x = pop, y = phi, fill=pop)) +
  geom_violin() +
  ylab("Lag-1 correlation in recruitment") +
  theme_bw() +
  coord_flip()+
  ylim(-0.5,1) +
  scale_fill_manual(values=rev(pop_colors_ord))+
  theme(axis.title.y = element_blank(),
        axis.text = element_text(size=6),
        axis.title.x = element_text(size=9),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none",
        plot.margin = unit(c(0.25,0.25,0.7,0.5), units = "lines"))

g <- ggarrange(a,ggarrange(b,c,nrow =2, labels = c("b", "c")),
               ncol = 2, labels = c("a"))

jpeg("./04_figures/figures/figure5.jpeg", width = 6, height = 5, units = "in", res = 600)
print(g)
dev.off()

# Figure 6: CV in run-size ----

indCVs <- matrix(NA,1000,8); colnames(indCVs)<-unique(spawn_df$population)
for (i in unique(spawn_df$population)){
  for(j in 1:1000){
    pop <- subset(spawn_df, population == i)[sample(1:1000,1),]
    indCVs[j,i] <-sd(as.numeric(pop[-1]))/mean(as.numeric(pop[-1]))
  }
}

aggCV <- matrix(NA,1000,1)
for(j in 1:1000){
  rsamp <- sample(1:1000,1);rsamps<-c(rsamp,rsamp+1000,rsamp+2000,rsamp+3000,
                                      rsamp+4000,rsamp+5000,rsamp+6000,rsamp+7000)
  agg<-colSums(spawn_df[rsamps,2:36])
  aggCV[j] <-sd(agg)/mean(agg)
}

CVs<-cbind(indCVs,aggCV)
colnames(CVs)<-c(as.character(unique(spawn_df$population)),"Aggregate")
CVs.df <- as.data.frame(CVs)                 
CV.long <-pivot_longer(data=CVs.df,cols=1:9, names_to = "pop", values_to = "cv")
CV.long$pop <- factor(CV.long$pop, levels = c("Aggregate","Lower Mainstem","White-Donjek","Stewart","Pelly",
                                              "Middle Mainstem",
                                              "Carmacks","Upper Lakes and Mainstem","Teslin"))
weights <- ensemble%>%
  group_by(population)%>%
  summarize(sum_size = sum(med))%>%
  mutate(weights = sum_size / sum(sum_size))

ensemble$population <-lapply(ensemble$population, as.character)
run_size_Carmacks <- rbind(ensemble[ensemble$population == "Upper Lakes and Mainstem",],
                           ensemble[ensemble$population == "Teslin",],
                           ensemble[ensemble$population == "Carmacks",]
)
run_size_Teslin <- rbind(ensemble[ensemble$population == "Teslin",])

PE_running_avg_full <- matrix(NA, 25,1)
PE_running_avg_Carmacks <- PE_running_avg_full
CV_running_avg_full <- PE_running_avg_full
CV_running_avg_Carmacks <- PE_running_avg_full
CV_running_avg_Teslin <- PE_running_avg_full

for (i in 1:25){
  PE_running_avg_full[i,] <- PE_running(ensemble, 1984 + i, 1994 + i)[1]
  PE_running_avg_Carmacks[i,] <- PE_running(run_size_Carmacks, 1984 + i, 1994 + i)[1]
  CV_running_avg_full[i,] <- PE_running(ensemble, 1984 + i, 1994 + i)[2]
  CV_running_avg_Carmacks[i,] <- PE_running(run_size_Carmacks, 1984 + i, 1994 + i)[2]
  CV_running_avg_Teslin[i,] <- PE_running(run_size_Teslin, 1984 + i, 1994 + i)[2]
}

full <- cbind(seq(1995,2019),CV_running_avg_full, rep("Lower",length(CV_running_avg_full)))
middle <- cbind(seq(1995,2019),CV_running_avg_Carmacks, rep("Middle",length(CV_running_avg_full)))
upper <- cbind(seq(1995,2019),CV_running_avg_Teslin, rep("Headwaters",length(CV_running_avg_full)))

rt_cvs <- rbind(full,middle,upper)
colnames(rt_cvs) <- c("Year", "CV", "Location")
rt_cvs<-as.data.frame(rt_cvs)
rt_cvs$Year <- as.numeric(as.character(rt_cvs$Year))
rt_cvs$CV <- as.numeric(as.character(rt_cvs$CV))
rt_cvs$Location <- factor(rt_cvs$Location , levels = c("Headwaters","Middle", "Lower"))

a <- ggplot(data = CV.long, aes(x = pop, y = cv, fill=pop)) +
  geom_violin() +
  ylab("CV in run size") +
  theme_bw() +
  coord_flip()+
  ylim(0,1) +
  scale_fill_manual(values=c("dark grey", pop_colors_ord))+
  theme(axis.title.y = element_blank(),
        axis.text = element_text(size=7),
        axis.title.x = element_text(size=9),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none",
        plot.margin = unit(c(0.25,0.25,0.7,0.5), units = "lines"))

b <- ggplot(data = rt_cvs, aes(x = Year, y = CV, color=Location)) +
  geom_line(size=1.2) +
  ylab("CV in run size") +
  scale_color_grey(start=0.8, end=0.2) +
  scale_y_continuous(breaks = c(0,0.15,0.3,0.45,0.6,0.75), limits=c(0,0.85),position = "right") +
  theme_bw() +
  theme(axis.title = element_text(size=9),
        axis.text.x = element_text(angle=45, hjust = 1, size=7),
        axis.text.y = element_text(size=7),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.justification = c(0,0),
        legend.position = c(0.52,0.05),
        legend.key.size = unit(7, "pt"),
        legend.background = element_blank(),
        legend.text = element_text(size = 7),
        legend.title = element_text(size = 8),
        plot.margin = unit(c(0.25,0.25,0,1), units = "lines"))

g <- ggarrange(a,b, labels = c("a", "b"),widths = c(1.3, 1), ncol = 2)

jpeg("./04_figures/figures/figure6.jpeg", width = 6, height = 2.5, units = "in", res = 600)
print(g)
dev.off()

# Figure S1: border passage and GSI samples -----
bc <- border_passage %>%
  drop_na(count) %>%
  group_by(year, count_type) %>%
  mutate(year_count = sum(count)) %>%
  mutate(max_count = max(count)) %>%
  ungroup() %>%
  mutate(prop = (count/year_count)) %>%
  mutate(prop2 = (count/max_count)) %>%
  mutate(year_gear = paste0(year,"_",count_type)) %>%
  filter(year_gear %not_in% c("2005_fishWheel","2006_fishWheel","2007_fishWheel")) %>%
  filter(year %not_in% c("1988","1989","1990","1998")) %>%
  as.data.frame()

gsi.prop <- gsi %>%
  mutate(year_gear = paste0(year,"_",gear)) %>%
  filter(year_gear %not_in% c("2008_Fish Wheel","2010_Fish Wheel","2011_Fish Wheel","2012_Fish Wheel")) %>%
  group_by(year, sample_num) %>%
  filter(prob == max(prob)) %>%
  ungroup() %>%
  group_by(year_gear) %>%
  mutate(total_count = n()) %>%
  ungroup() %>%
  group_by(year_gear, region_name) %>%
  mutate(pop_count = n(),
         pop_prop = (pop_count/total_count)*100) %>%
  ungroup() %>%
  distinct(year, region_name, pop_prop) %>%
  spread(key=region_name, pop_prop) %>%
  as.data.frame()

gsi.count <- gsi %>%
  mutate(year_gear = paste0(year,"_",gear)) %>%
  filter(year_gear %not_in% c("2008_Fish Wheel","2010_Fish Wheel","2011_Fish Wheel","2012_Fish Wheel")) %>%
  group_by(year, sample_num) %>%
  filter(prob == max(prob)) %>%
  ungroup() %>%
  group_by(year_gear) %>%
  mutate(total_count = n()) %>%
  ungroup() %>%
  group_by(year_gear, region_name) %>%
  mutate(pop_count = n(),
         pop_prop = (pop_count/total_count)*100) %>%
  ungroup() %>%
  distinct(year, region_name, pop_count) %>%
  spread(key=region_name, pop_count) %>%
  as.data.frame()

gsi.join <- full_join(gsi.count, gsi.prop) %>% 
  dplyr::select(1, 3, 6, 5, 9, 4, 2, 8, 7) %>%
  as.data.frame()
#write.csv(gsi.join, file = "gsi.count.csv")

gsi %>% 
  group_by(year, sample_num) %>%
  filter(prob == max(prob)) %>% 
  ungroup() %>%
  group_by(year) %>%
  summarise(med = median(prob)) %>% 
  ungroup() %>%
  mutate(max = max(med),
         min = min(med)) %>%
  summarise(quants = quantile(med, probs = c(0.05, 0.5, 0.95))) %>%
  as.data.frame()

gsi2 <- gsi %>%
  mutate(year_gear = paste0(year,"_",gear)) %>%
  filter(year_gear %not_in% c("2008_Fish Wheel","2010_Fish Wheel","2011_Fish Wheel","2012_Fish Wheel")) %>%
  group_by(year, sample_num) %>%
  filter(prob == max(prob)) %>%
  ungroup() %>%
  group_by(year) %>%
  mutate(year_count = n()) %>%
  ungroup() %>%
  group_by(year, julian) %>%
  mutate(julian_count = n()) %>%
  ungroup() %>%
  group_by(year) %>%
  mutate(max_count = max(julian_count)) %>%
  ungroup() %>%
  mutate(julian_prop = (julian_count/year_count)) %>%
  mutate(julian_prop2 = (julian_count/max_count)) %>%
  distinct(year, year_count, julian, julian_prop,julian_prop2) %>%
  as.data.frame()

gsi2 %>% distinct(year, gear) %>% arrange(gear, year)
bc %>% distinct(year, count_type) %>% arrange(count_type, year)
bc %>% group_by(year) %>% summarise(sum = sum(prop)) %>% as.data.frame()
gsi2 %>% group_by(year) %>% summarise(sum = sum(julian_prop)) %>% as.data.frame()

a <- gsi2 %>% distinct(year, year_count) %>% arrange(year) %>% 
  mutate(xpos = 200,
         ypos = 0,
         vjustvar = 0) %>%
  as.data.frame()

g <- ggplot() +
  geom_vline(xintercept=c(180,210,240), color="light grey",lwd=0.25,lty=2) +
  geom_bar(data = bc, aes(x=as.numeric(julian), y=prop2),stat = "identity") +
  geom_bar(data = gsi2, aes(x=as.numeric(julian), y=(julian_prop2)*-1), 
           fill="red", alpha=0.85, stat = "identity") +
  scale_x_continuous(limits=c(175,250), breaks = c(180,210,240)) +
  xlab("Day of year") +
  ylab("Run and GSI sample sizes") +
  facet_wrap(~year, scales = "free_y") +
  theme_bw() +
  theme(axis.text.x = element_text(size=8, angle = 45, hjust = 1),
        strip.text = element_text(size=9),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        aspect.ratio=1,
        panel.spacing.x=unit(0.1, "lines"),
        panel.spacing.y=unit(0.3, "lines"),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank())+
  geom_text(data = a, 
            mapping = aes(x = 237, y = -0.5, label = year_count, hjust = 1, vjust = 2),
            size=3, color = "red")
jpeg("04_figures/figures/figureS1.jpeg", width = 6, height = 8, units = "in", res = 600)
print(g)
dev.off()


# Figure S11: Age composition ----
ASL$new_AGE <- as.numeric(ASL$new_AGE)

data <- ASL %>% 
  mutate(year_gear = paste0(year,"_",gear)) %>%
  filter(year_gear %not_in% c("2008_Fish Wheel","2010_Fish Wheel","2011_Fish Wheel","2012_Fish Wheel")) %>%
  drop_na(new_AGE) %>%
  group_by(gear, year, sample_num) %>%
  filter(prob == max(prob)) %>%
  ungroup() %>%
  mutate(new_AGE = case_when(new_AGE == 8 ~ 7,
                             new_AGE == 3 ~ 4,
                             new_AGE %in% c(4,5,6,7) ~ new_AGE)) %>%
  group_by(gear, year, region_name) %>%
  mutate(count = n()) %>%
  ungroup() %>%
  group_by(gear, year, new_AGE, region_name) %>%
  mutate(age_count = n()) %>%
  ungroup() %>%
  mutate(age_prop = age_count/count) %>%
  dplyr::select(year, gear, region_name, new_AGE, count, age_count, age_prop) %>%
  arrange(gear, region_name, new_AGE, year) %>%
  mutate(Population = case_when(region_name == "Yukon Carmacks" ~ "Carmacks",
                                region_name == "Yukon Lower Canadian" ~ "Lower Mainstem",
                                region_name == "Yukon mainstem" ~ "Middle Mainstem",
                                region_name == "Yukon Pelly" ~ "Pelly",
                                region_name == "Yukon Stewart" ~ "Stewart",
                                region_name == "Yukon upper" ~ "Upper Lakes and Mainstem",
                                region_name == "Yukon White-Donjek" ~ "White-Donjek",
                                region_name == "Yukon Teslin" ~ "Teslin")) %>%
  rename(Fish.Age = new_AGE) %>%
  as.data.frame()

data$Fish.Age <- as.factor(data$Fish.Age)

jtc2 <- jtc %>%
  dplyr::select(year, age_4, age_5, age_6, age_7) %>%
  gather(age, age_prop, age_4:age_7) %>%
  mutate(Population = "Aggregate") %>%
  mutate(Fish.Age = case_when(age == "age_4" ~ "4",
                              age == "age_5" ~ "5",
                              age == "age_6" ~ "6",
                              age == "age_7" ~ "7")) %>%
  dplyr::select(year, Population, Fish.Age, age_prop) %>%
  as.data.frame()

data2 <- full_join(data, jtc2) %>% arrange(year) %>% as.data.frame()

data2$Population <- factor(data2$Population, levels = c("Lower Mainstem","White-Donjek","Stewart","Pelly",
                                                        "Middle Mainstem","Carmacks","Upper Lakes and Mainstem","Teslin","Aggregate"))

g <- ggplot(data2, aes(x=year, y=age_prop, fill=Fish.Age)) +
  geom_bar(stat= "identity", position = "fill") +
  scale_x_continuous(limits = c(1985,2018), breaks = seq(1985,2020,5)) +
  scale_fill_viridis(option = "viridis", discrete = TRUE) +
  xlab("Year") +
  ylab("Proportion") +
  facet_wrap(~Population) +
  labs(fill="Age class")+
  theme_bw() +
  theme(axis.text.x = element_text(angle=45, hjust = 1, vjust=1, size=6),
        strip.text = element_text(size=6),
        axis.title = element_text(size=9),
        axis.text = element_text(size=6),
        legend.key.size = unit(10, "pt"),
        legend.background = element_blank(),
        legend.text = element_text(size = 7),
        legend.title = element_text(size = 9),
        panel.grid.minor = element_blank())


jpeg("04_figures/figures/figures11.jpeg", width = 6, height = 5, units = "in", res = 200)
print(g)
dev.off()

# summary stats for main text ----

# border passage
ensemble %>%
  group_by(pops_f) %>%
  summarize(mean_passage = mean(med))

# spawner-recruitment parameters
alpha_beta_per 

# median lag-1 autocorrelation in rec resids
median(alpha_beta_per$phi_med)

# correlation in recruitment resids
mean(corr_matrix[row(corr_matrix)!=col(corr_matrix)])

# median % change in productivity between periods
early_resids <- logresid_perc_scale %>%
                  group_by(population) %>%
                  filter(BroodYear %in% (1990:1995)) %>%
                  summarise(mean(med))

late_resids <- logresid_perc_scale %>%
  group_by(population) %>%
  filter(BroodYear %in% (2005:2010)) %>%
  summarise(mean(med))

median(-1-(late_resids$`mean(med)`/early_resids$`mean(med)`))

  
# harvest-risk tradeoffs

t3.3$med_V2[which(t3.3$med_V2==max(t3.3$med_V2))] # median max-mixed stock harvest
t3.3$low_V2[which(t3.3$low_V2==max(t3.3$low_V2))] # lower CI max-mixed stock harvest
t3.3$upp_V2[which(t3.3$upp_V2==max(t3.3$upp_V2))] # upper CI max-mixed stock harvest

t3.3$U[which(t3.3$med_V2==max(t3.3$med_V2))] # median max-mixed stock harvest rate
t3.3$U[which(t3.3$low_V2==max(t3.3$low_V2))] # lower CI max-mixed stock harvest rate
t3.3$U[which(t3.3$upp_V2==max(t3.3$upp_V2))] # upper CI max-mixed stock harvest rate

t3.3$over.med[which(t3.3$med_V2==max(t3.3$med_V2))] # median overfished at max-mixed stock harvest 
t3.3$over.low[which(t3.3$low_V2==max(t3.3$low_V2))] # lower CI overfished at max-mixed stock harvest
t3.3$over.up[which(t3.3$upp_V2==max(t3.3$upp_V2))] # upper CI overfished at max-mixed stock harvest

t3.3$ext.med[which(t3.3$med_V2==max(t3.3$med_V2))] # median overfished at max-mixed stock harvest 
t3.3$ext.low[which(t3.3$low_V2==max(t3.3$low_V2))] # lower CI overfished at max-mixed stock harvest
t3.3$ext[which(t3.3$upp_V2==max(t3.3$upp_V2))] # upper CI overfished at max-mixed stock harvest

# run-size CVs
colMeans(indCVs)
colMeans(aggCV)

# run-size variance dampening
var_damp <- matrix(NA,1000,1)
for(i in 1:1000){
  sub_samp <- CVs.df[sample(1:1000,1),]
  weight_mean_CV <- sum(sub_samp[1:8]*weights$weights)
  var_damp[i,1] <- as.numeric(weight_mean_CV/sub_samp[9])
}

quantile(as.vector(var_damp), probs=c(0.025,0.5,0.975))

#var_damp_unweight <- matrix(NA,1000,1) 
#for(i in 1:1000){
#  sub_samp <- CVs.df[sample(1:1000,1),]
#  var_damp_unweight[i,1] <- as.numeric((sum(sub_samp[1:8])/8)/sub_samp[9])
#}

#quantile(as.vector(var_damp_unweight), probs=c(0.025,0.5,0.975))

# difference in CVs between sections of river
rt_cvs %>%
  group_by(Location) %>%
  summarise(mean(CV))