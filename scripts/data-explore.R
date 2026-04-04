#### DATA EXPLORATION ####

## packages
library(tidyverse)

## data
data <- read_csv("data/gen_data_final_fonseca.csv")

## how many are shifting in the wrong direction?
data %>% 
  select(shift_vel_sign) %>% 
  mutate(across(everything(), as.factor)) %>%
  summary()
  
## i want to see if Nucleotide_diversity is the same in multiple entriies of the same spp
data %>% 
  select(spp, Nucleotide_diversity) %>% 
  group_by(spp) %>% 
  summarise(n = n(),
            mean = mean(Nucleotide_diversity),
            sd = sd(Nucleotide_diversity))
