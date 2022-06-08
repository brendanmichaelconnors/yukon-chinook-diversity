#------------------------------------------------------------------------------#
# Link GSI and ASL data
# just for test fishery data (Eagle 2005 +) 
#------------------------------------------------------------------------------#

# 1. load data ----
require(tidyverse)

ASL_EAGLE <- read.csv("01_inputs/data/ASL_Output_Chinook_Eagle_2005-2019.csv", stringsAsFactors = FALSE) %>% filter(sampleYear != 2019)
GSI <- read.csv("01_inputs/data/20Oct2021-update.gsiSamplesAllProbs.csv", stringsAsFactors = FALSE) %>% filter(year != 2019)

head(ASL_EAGLE)
head(GSI)

# 2. merge ASL and GSI datasets ----
ASL_EAGLE.1 <- ASL_EAGLE %>%
  mutate(AGEx = ageFresh + ageSalt + 1) %>%
  mutate(AGE = case_when(AGEx == 3 ~ 4,
                         AGEx == 8 ~ 7,
                         AGEx %in% c(4,5,6,7) ~ AGEx)) %>%
  dplyr::select(-AGEx) %>%
  dplyr::select(gear, sampleYear, sampleDate, Genetic.Sample.Number, sexID, length, ageFresh, ageSalt, totalAge, AGE) %>%
  rename(year = sampleYear, Length = length, Fresh.Water.Age = ageFresh, 
         Salt.Water.Age = ageSalt, Total.Age = totalAge, sample_num = Genetic.Sample.Number) %>%
  mutate(gear = case_when(gear == "Drift Gillnet" ~ "Test Fishery",
                          gear == "Set Gillnet" ~ "Test Fishery"),
         Sex = case_when(sexID == 1 ~ "male",
                         sexID == 2 ~ "female")) %>%
  dplyr::select(-sexID) %>%
  as.data.frame()

eagle_GSI_ASL <- GSI %>%
  filter(year %in% c(2008:2018),
         gear == "Test Fishery") %>%
  left_join(., ASL_EAGLE.1, by = c("year", "sample_num")) %>%
  group_by(year, julian, sample_num) %>%
  filter(prob == max(prob)) %>%
  arrange(year, sample_num) %>%
  as.data.frame()

write.csv(eagle_GSI_ASL, file = "01_inputs/data/ASL_GSI.csv")

# 3. calculte age comps ----
eagle_GSI_ASL_age_comps <- eagle_GSI_ASL %>%
  drop_na(AGE) %>%
  select(year, region, AGE) %>%
  group_by(year, region, AGE) %>%
  summarize(AGE_COUNT = n()) %>%
  mutate(Population = case_when(region == 36 ~ "Lower Mainstem",
                                region == 34 ~ "Pelly",
                                region == 35 ~ "Stewart",
                                region == 33 ~ "Middle Mainstem",
                                region == 38 ~ "White-Donjek",
                                region == 32 ~ "Carmacks",
                                region == 30 ~ "Upper Lakes and Mainstem",
                                region == 31 ~ "Teslin")) %>%
  mutate(Population_order = case_when(region == 36 ~ 1,
                                region == 34 ~ 3,
                                region == 35 ~ 4,
                                region == 33 ~ 7,
                                region == 38 ~ 2,
                                region == 32 ~ 5,
                                region == 30 ~ 8,
                                region == 31 ~ 6)) %>%
  spread(key = AGE, value = AGE_COUNT) %>%
  replace_na(list("4" = 0, "5" = 0, "6" = 0, "7" = 0)) %>%
  arrange(Population_order, year) %>%
  gather(key = AGE, value = AGE_COUNT, "4":"7") %>%
  arrange(Population_order,  AGE, year) %>%
  as.data.frame()

age_comp <- array(eagle_GSI_ASL_age_comps$AGE_COUNT,dim=c(11,4,8), dimnames = list(seq(2008,2018), 
                                                    seq(4,7),
                                                    c("Lower Mainstem" ,
                                                      "White-Donjek",
                                                      "Pelly",
                                                      "Stewart" ,
                                                      "Carmacks", 
                                                      "Teslin",
                                                      "Middle Mainstem",
                                                      "Upper Lakes and Mainstem")))


saveRDS(age_comp, file = "01_inputs/data/pop_age_comp_2008_2019.rds")

