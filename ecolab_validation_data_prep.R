# Title: Contemporary endhopyte prevalence surveys
# Purpose: Imports survey data collected between 2018 and 2020 as part of the Texas Ecolab program, then saves data for model validation
# Authors: Joshua Fowler
# Date Updated: Mar 2, 2023
library(tidyverse)
library(readxl)
library(lubridate)

# ecolab scores, primarily for A. hyemalis, scores from 2018 and a few sites from 2015 in first csv
scores_2018 <- read_excel(path ="~/Dropbox/Endophyte removal experiment 2018-2019/ecolab2.0.xlsx", sheet = "Plants Scored"  ) %>% 
  rename(SiteID = siteID, PlantID = individualID, SpeciesID = speciesID, seeds.scored = s.scored, seeds.eplus = s.eplus, addit.seed.scored = `add. s.scored`,addit.seed.eplus = `add. s.eplus` ) %>%
  mutate(PlantID = as.character(PlantID)) %>% 
  dplyr::select(-Notes, -date, -`total s.scored`, -`total eplus`, -vwc)

locations_2018 <- read_excel(path = "~/Dropbox/Endophyte removal experiment 2018-2019/ecolab2.0.xlsx") %>% 
  select(-contact, -email) %>% 
  rename(SiteID = siteID) %>% 
  rename(lat = gps_lat, lon = gps_lon) %>% 
  select(SiteID, county, Year, lat, lon)

scores_2019 <- read_excel(path = "~/Dropbox/Josh&Tom - shared/ecolab2019/Endophyte_field_collections_2019.xlsx", sheet = "Endophyte_scoring_2019") %>% 
  rename(PlantID = Plantid, addit.seed.eplus = addit.seed_eplus) %>% 
  dplyr::select(SiteID, SpeciesID, PlantID, seeds.scored, seeds.eplus, addit.seed.scored, addit.seed.eplus) %>% 
  mutate(seeds.eplus = as.numeric(seeds.eplus),
         seeds.scored = as.numeric(seeds.scored)) %>% 
  mutate(Year = 2019) 
  

locations_2019 <- read_excel("~/Dropbox/Josh&Tom - shared/ecolab2019/Endophyte_field_collections_2019.xlsx", sheet = "Collection Locations") %>% 
  rename(lon = `longitude (dd)`, lat = `latitude (dd)`) %>% 
  mutate(Year = 2019) %>% 
  dplyr::select(-`elevation (m)`, -time)

scores_2020 <- read_excel(path = "~/Dropbox/Josh&Tom - shared/ecolab2020/Ecolab_2020_AvaJohnson/Endophyte_field_collections_2020.xlsx", sheet = "Endophyte_scoring_2020") %>% 
  rename(PlantID = Plantid, addit.seed.eplus = addit.seed_eplus) %>% 
  dplyr::select(SiteID, SpeciesID, PlantID,seeds.scored, seeds.eplus, addit.seed.scored, addit.seed.eplus) %>% 
  mutate(PlantID = as.character(PlantID)) %>% 
  mutate(SiteID = case_when(SiteID == "Mont20-1" ~ "MONT20-1",
                            SiteID == "RL-1" ~ "RL1",
                            TRUE ~ SiteID)) %>% 
  mutate(Year = 2020) %>% 
  mutate(seeds.eplus = as.numeric(seeds.eplus),
         seeds.scored = as.numeric(seeds.scored))
  
  
locations_2020 <- read_excel(path = "~/Dropbox/Josh&Tom - shared/ecolab2020/Ecolab_2020_AvaJohnson/Endophyte_field_collections_2020.xlsx", sheet = "Collection Locations") %>%
  rename(lon = `longitude (dd)`, lat = `latitude (dd)`) %>% 
  mutate(Year = 2020) %>% 
  dplyr::select(-notes, -time)

ecolab_2018 <- scores_2018 %>% 
  left_join(locations_2018) 
ecolab_2019 <- scores_2019 %>% 
  full_join(locations_2019) 
ecolab_2020 <- scores_2020 %>% 
  full_join(locations_2020) 

ecolab <- ecolab_2018 %>% 
  full_join(ecolab_2019) %>% 
  full_join(ecolab_2020) %>% 
  mutate(total_seeds = case_when(is.na(addit.seed.scored) ~ seeds.scored,
                                 !is.na(addit.seed.scored) ~ seeds.scored + addit.seed.scored),
         total_eplus = case_when(is.na(addit.seed.eplus) ~ seeds.eplus,
                        !is.na(addit.seed.eplus) ~ seeds.eplus + addit.seed.eplus),
        endo_status = case_when(total_eplus>0 ~ 1, total_eplus<=0 ~ 0))

ecolab_prevalence <- ecolab %>% 
  filter(!is.na(seeds.scored) & seeds.scored>0) %>% 
  filter(!is.na(lat)) %>% 
  group_by(SiteID,  SpeciesID, lat, lon, Year) %>% 
  summarize(sample_size = n(),
            endo_prev = sum(endo_status)/sample_size) %>% 
  filter(SpeciesID == "AGHY")


# scores from Sneck et al. 2017, for E. virginicus
sneck_surveys <- read_csv(file = "~/Dropbox/Josh&Tom - shared/Endo_Herbarium/sneckelvielca.csv") %>% 
  filter(!is.na(Site)) %>% 
  mutate(SpeciesID = case_when(.$Species ==  "Both" ~ "ELCA;ELVI", 
                             .$Species == "EV" ~ "ELVI",
                             .$Species == "EC" ~ "ELCA")) %>% 
  separate(`E+ (%)`, c("Eprev1","Eprev2"), sep = ";") %>% 
  separate(`Plant N`, c("PlantNo1", "PlantNo2", sep = ";")) %>% 
  mutate(endo_prev = case_when(SpeciesID == "ELVI" ~ .01*(as.numeric(Eprev1)),
                               SpeciesID == "ELCA;ELVI" ~ .01*(as.numeric(Eprev2))),
         sample_size = case_when(SpeciesID == "ELVI" ~ as.numeric(PlantNo1),
                                 SpeciesID == "ELCA;ELVI" ~ as.numeric(PlantNo2))) %>% 
  dplyr::select( `Site name`,  SpeciesID, endo_prev, sample_size, `Lat.`, `Long.`, `Sampling Date`) %>% 
  filter(grepl("ELVI", SpeciesID)) %>% 
  rename(SiteID = `Site name`, date = `Sampling Date`, lat = `Lat.`, lon = `Long.`) %>% 
  mutate(SpeciesID = "ELVI",
         Year = year(mdy(date))) %>% 
  select(-date)


# now combining the AGHY and ELVI scores
contemp_surveys <- ecolab_prevalence %>% 
  full_join(sneck_surveys)

write_csv(contemp_surveys, file = "contemp_surveys.csv")


ecolab_random_sample <- ecolab %>% 
  filter(SpeciesID == "AGHY") %>% 
  filter(!is.na(seeds.scored) & seeds.scored>0) %>% 
  filter(!is.na(lat)) %>% 
  group_by(SiteID) %>% 
  slice_sample(n = 1) 
write_csv(ecolab_random_sample, file = "contemp_random_sample.csv")

ggplot(ecolab_random_sample)+
  geom_point(aes(x = lon, y = endo_status))


# getting the county for each survey
library(ggmap)
register_google(key = "AIzaSyDuAdpozRqmb8Sms-XfivxLi3tzlifJdMw")
contemp_surveys_county <- ecolab_prevalence %>% 
  mutate(coords = as.character(paste(lat, lon))) %>% 
  geocode(coords, output = "more")
  
str(geocode("33.683441 -93.987111", output = "more")$address)
  



