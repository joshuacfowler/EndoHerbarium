library(ggplot2)
library(tidyverse)

# Poa autumnalis records
occur <- read.csv(file = "~joshuacfowler/Downloads/SymbOutput_2019-01-22_134351_DwC-A/occurrences.csv")
View(occur)
colnames(occur)
occur <- occur %>% 
  filter(year>=100)

ggplot(data = occur)+
  stat_count(aes(x = year)) + 
  labs(title ="POAU Collection counts by year from NANSH", x = "county", y = "count")+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))


ggplot(data = occur) +
  stat_count(aes(x = stateProvince))+
  labs(title ="POAU Collection counts by State from NANSH", x = "county", y = "count")+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))


occurtexas <- occur %>% 
  filter(stateProvince == "Texas")

ggplot(data = occurtexas)+
  stat_count(aes(x = county)) +
  labs(title ="POAU Collection counts in TX from NANSH", x = "county", y = "count")+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

years_locations <- occur %>% 
  select(year, eventDate, county, stateProvince)
