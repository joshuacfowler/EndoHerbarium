library(ggplot2)
library(tidyverse)

# Poa autumnalis specimen records
occurPOAU <- read.csv(file = "~joshuacfowler/Downloads/SymbOutput_2019-01-22_134351_DwC-A/occurrences.csv")
View(occurPOAU)
occurPOAU <- occurPOAU %>% 
  filter(year>=100)

ggplot(data = occurPOAU)+
  stat_count(aes(x = year)) + 
  labs(title ="POAU Collection counts by year from NANSH", x = "county", y = "count")+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))


ggplot(data = occurPOAU) +
  stat_count(aes(x = stateProvince))+
  labs(title ="POAU Collection counts by State from NANSH", x = "county", y = "count")+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))


occurPOAUtexas <- occurPOAU %>% 
  filter(stateProvince == "Texas")

ggplot(data = occurPOAUtexas)+
  stat_count(aes(x = county)) +
  labs(title ="POAU Collection counts in TX from NANSH", x = "county", y = "count")+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

years_locations_POAU <- occurPOAU %>% 
  select(year, eventDate, county, stateProvince)


ggplot(data = occurPOAUtexas)+
  geom_tile(aes(x =county, y = year))

# Agrostis hyemalis specimen records
occurAGHY <- read.csv(file = "~joshuacfowler/Downloads/SymbOutput_2019-01-29_131843_DwC-A/occurrences.csv")
View(occurAGHY)
occurAGHY <- occurAGHY %>% 
  filter(year>=100)

ggplot(data = occurAGHY)+
  stat_count(aes(x = year)) + 
  labs(title ="AGHY Collection counts by year from NANSH", x = "county", y = "count")+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))


ggplot(data = occurAGHY) +
  stat_count(aes(x = stateProvince))+
  labs(title ="AGHY Collection counts by State from NANSH", x = "county", y = "count")+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))


occurAGHYtexas <- occurAGHY %>% 
  filter(stateProvince == "Texas")

ggplot(data = occurAGHYtexas)+
  stat_count(aes(x = county)) +
  labs(title ="AGHY Collection counts in TX from NANSH", x = "county", y = "count")+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

years_locations_AGHY <- occurAGHY %>% 
  select(year, eventDate, county, stateProvince)
View(years_locations_AGHY)


ggplot(data = occurAGHYtexas)+
  geom_tile(aes(x =county, y = year))

# Agrostis perenans specimen records
occurAGPE <- read.csv(file = "~joshuacfowler/Downloads/SymbOutput_2019-01-29_162445_DwC-A/occurrences.csv")
View(occurAGPE)
occurAGPE <- occurAGPE %>% 
  filter(year>=1500)

ggplot(data = occurAGPE)+
  stat_count(aes(x = year)) + 
  labs(title ="AGPE Collection counts by year from NANSH", x = "county", y = "count")+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))


ggplot(data = occurAGPE) +
  stat_count(aes(x = stateProvince))+
  labs(title ="AGPE Collection counts by State from NANSH", x = "county", y = "count")+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))


occurAGPEtexas <- occurAGPE %>% 
  filter(stateProvince == "Texas")

ggplot(data = occurAGPEtexas)+
  stat_count(aes(x = county)) +
  labs(title ="AGPE Collection counts in TX from NANSH", x = "county", y = "count")+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

years_locations_AGPE <- occurAGPE %>% 
  select(year, eventDate, county, stateProvince)
View(years_locations_AGPE)

ggplot(data = occurAGPEtexas)+
  geom_tile(aes(x =county, y = year))

# This file download is funky and doesn't work.
# Elymus virginicus specimen records
occurELVI <- read.csv2(file = "~joshuacfowler/Downloads/SymbOutput_2019-01-29_163142_DwC-A/occurrences.csv", sep = ",")
View(occurELVI)
occurELVI <- occurELVI %>% 
  filter(year>=1500)

ggplot(data = occurELVI)+
  stat_count(aes(x = year)) + 
  labs(title ="ELVI Collection counts by year from NANSH", x = "county", y = "count")+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))


ggplot(data = occurELVI) +
  stat_count(aes(x = stateProvince))+
  labs(title ="ELVI Collection counts by State from NANSH", x = "county", y = "count")+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))


occurELVItexas <- occurELVI %>% 
  filter(stateProvince == "Texas")

ggplot(data = occurELVItexas)+
  stat_count(aes(x = county)) +
  labs(title ="ELVI Collection counts in TX from NANSH", x = "county", y = "count")+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

years_locations_ELVI <- occurELVI %>% 
  select(year, eventDate, county, stateProvince)
View(years_locations_ELVI)


ggplot(data = occurELVItexas)+
  geom_tile(aes(x =county, y = year))

# Poa sylvestris specimen records
occurPOSY <- read.csv2(file = "~joshuacfowler/Downloads/SymbOutput_2019-01-29_163447_DwC-A/occurrences.csv", sep = ",")
View(occurPOSY)
occurPOSY <- occurPOSY %>% 
  filter(year>=100)

ggplot(data = occurPOSY)+
  stat_count(aes(x = year)) + 
  labs(title ="POSY Collection counts by year from NANSH", x = "county", y = "count")+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))


ggplot(data = occurPOSY) +
  stat_count(aes(x = stateProvince))+
  labs(title ="POSY Collection counts by State from NANSH", x = "county", y = "count")+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))


occurPOSYtexas <- occurPOSY %>% 
  filter(stateProvince == "Texas")

ggplot(data = occurPOSYtexas)+
  stat_count(aes(x = county)) +
  labs(title ="POSY Collection counts in TX from NANSH", x = "county", y = "count")+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

years_locations_POSY <- occurPOSY %>% 
  select(year, eventDate, county, stateProvince)
View(years_locations_POSY)

ggplot(data = occurPOSYtexas)+
  geom_tile(aes(x =county, y = year))




# UT Austin collections ---------------------------------------------------
########## UT Austin collections
# Poa autumnalis
UTAoccurPOAU <- read.csv(file = "~joshuacfowler/Downloads/webreq_DwC-A/occurrences.csv")
View(UTAoccurPOAU)
occurPOAU <- UTAoccurPOAU %>% 
  filter(year>=100)

ggplot(data = UTAoccurPOAU)+
  stat_count(aes(x = year)) + 
  labs(title ="POAU Collection counts by year from UTA", x = "county", y = "count")+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

occurPOAUtexas <- occurPOAU %>% 
  filter(stateProvince == "Texas")

ggplot(data = occurPOAUtexas)+
  stat_count(aes(x = county)) +
  labs(title ="POAU Collection counts in TX from NANSH", x = "county", y = "count")+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

UTA_years_locations_POAU <- UTAoccurPOAU %>% 
  select(catalogNumber, year, eventDate, county, stateProvince)

View(UTA_years_locations_POAU)
write.csv(UTA_years_locations_POAU)

ggplot(data = UTAoccurPOAU)+
  geom_tile(aes(x =county, y = year))+
  labs(title ="POAU Collection dates by county in TX from UTA")+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))


# Poa sylvestris - only one record
UTAoccurPOSY <- read.csv(file = "~joshuacfowler/Downloads/webreq_DwC-A (1)/occurrences.csv")
View(UTAoccurPOSY)
occurPOSY <- UTAoccurPOSY %>% 
  filter(year>=100)

ggplot(data = UTAoccurPOSY)+
  stat_count(aes(x = year)) + 
  labs(title ="POSY Collection counts by year from UTA", x = "county", y = "count")+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))


ggplot(data = UTAoccurPOSY) +
  stat_count(aes(x = stateProvince))+
  labs(title ="POSY Collection counts by State from NANSH", x = "county", y = "count")+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))




UTA_years_locations_POSY <- UTAoccurPOSY %>% 
  select(catalogNumber, year, eventDate, county, stateProvince)

View(UTA_years_locations_POSY)
write.csv(UTA_years_locations_POSY)

ggplot(data = UTAoccurPOSY)+
  geom_tile(aes(x =county, y = year))+
  labs(title ="POSY Collection dates by county in TX from UTA")+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

# Elymus virginicus 
UTAoccurELVI <- read.csv(file = "~joshuacfowler/Downloads/webreq_DwC-A (2)/occurrences.csv")
View(UTAoccurELVI)
occurELVI <- UTAoccurELVI %>% 
  filter(year>=100)

ggplot(data = UTAoccurELVI)+
  stat_count(aes(x = year)) + 
  labs(title ="ELVI Collection counts by year from UTA", x = "county", y = "count")+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))


ggplot(data = occurELVI) +
  stat_count(aes(x = stateProvince))+
  labs(title ="ELVI Collection counts by State from NANSH", x = "county", y = "count")+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))


occurELVItexas <- occurELVI %>% 
  filter(stateProvince == "Texas")

ggplot(data = occurELVItexas)+
  stat_count(aes(x = county)) +
  labs(title ="ELVI Collection counts in TX from UTA", x = "county", y = "count")+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

UTA_years_locations_ELVI <- UTAoccurELVI %>% 
  select(catalogNumber, year, eventDate, county, stateProvince)

View(UTA_years_locations_ELVI)

write.csv(UTA_years_locations_ELVI)
ggplot(data = UTAoccurELVI)+
  geom_tile(aes(x =county, y = year))+
  labs(title ="ELVI Collection dates by county in TX from UTA")+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))



# Agrostis hyemalis 
UTAoccurAGHY <- read.csv(file = "~joshuacfowler/Downloads/webreq_DwC-A (3)/occurrences.csv")
View(UTAoccurAGHY)
occurAGHY <- UTAoccurAGHY %>% 
  filter(year>=100)

ggplot(data = UTAoccurAGHY)+
  stat_count(aes(x = year)) + 
  labs(title ="AGHY Collection counts by year from UTA", x = "county", y = "count")+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))


ggplot(data = occurAGHY) +
  stat_count(aes(x = stateProvince))+
  labs(title ="AGHY Collection counts by State from NANSH", x = "county", y = "count")+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))


occurAGHYtexas <- occurAGHY %>% 
  filter(stateProvince == "Texas")

ggplot(data = occurAGHYtexas)+
  stat_count(aes(x = county)) +
  labs(title ="AGHY Collection counts in TX from UTA", x = "county", y = "count")+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))


ggplot(data = UTAoccurAGHY)+
  geom_tile(aes(x =county, y = year))+
  labs(title ="AGHY Collection dates by county in TX from UTA")+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

UTA_years_locations_AGHY <- UTAoccurAGHY %>% 
  select(catalogNumber, year, eventDate, county, stateProvince)

View(UTA_years_locations_AGHY)
write.csv(UTA_years_locations_AGHY)