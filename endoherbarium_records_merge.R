# Title: Digital Herbarium Records and Sample Integration
# Purpose: Imports and merges downloaded digitized herbarium records with the Endo_Herbarium database
# Authors: Joshua Fowler and Lani Dufresne
# Updated: Mar 23, 2020

library(tidyverse)
library(fuzzyjoin)
library(readxl)
library(ggmap)
library(lubridate)

# Read in digitized herbarium records
# UT Austin downloaded from TORCH
AGHY_UTAustin <- read_csv(file = "~/Dropbox/Josh&Tom - shared/Endo_Herbarium/DigitizedHerbariumRecords/UTAustinAGHYrecords/occurrences.csv") %>% 
  mutate(new_id = gsub("[a-zA-Z ]", "", catalogNumber)) %>%   #creates a new id that we will use to merge with the 
  mutate(municipality = as.character(municipality)) %>% 
  dplyr::select(id, catalogNumber, country, stateProvince, county, municipality, locality, decimalLatitude, decimalLongitude, coordinateUncertaintyInMeters, new_id, eventDate, day, month, year) %>% 
  mutate(Spp_code = "AGHY")
ELVI_UTAustin <- read_csv(file = "~/Dropbox/Josh&Tom - shared/Endo_Herbarium/DigitizedHerbariumRecords/UTAustinELVIrecords/occurrences.csv") %>% 
  mutate(new_id = gsub("[a-zA-Z ]", "", catalogNumber)) %>%    #creates a new id that we will use to merge with the 
  mutate(municipality = as.character(municipality)) %>% 
  dplyr::select(id, catalogNumber, country, stateProvince, county, municipality, locality, decimalLatitude, decimalLongitude, coordinateUncertaintyInMeters, new_id, eventDate, day, month, year) %>% 
  mutate(Spp_code = "ELVI")
  
UTAustin_torch <- rbind(AGHY_UTAustin, ELVI_UTAustin)

# Texas A&M digitized records (Includes both AGHY and ELVI)
AM_records <- read_csv(file = "~/Dropbox/Josh&Tom - shared/Endo_Herbarium/DigitizedHerbariumRecords/TexasA&M_digitized_records/Fowler Data.csv") %>% 
  unite(Institution_specimen_id, CollectionCode, id, sep = "") %>%
  separate(DateCollected, into = c("year", "month", "day"), remove = FALSE) %>% 
  mutate(year = as.numeric(year), month = as.numeric(month), day = as.numeric(day)) %>% 
  select(-contains("X"))


# BRIT digitized records downloaded from TORCH (includes Vanderbilt and U of Louisiana Monroe) 
# This was downloaded on Apr 5, and we can get more specimens transcribed and download again.
AGHY_BRIT <- read_csv(file = "~/Dropbox/Josh&Tom - shared/Endo_Herbarium/DigitizedHerbariumRecords/BRIT_records/BRIT_AGHY_TORCH_records/SymbOutput_2020-04-05_162933_DwC-A/occurrences.csv") %>% 
  filter(!is.na(county), !is.na(eventDate)) %>%
  dplyr::select(id, catalogNumber, country, stateProvince, county, municipality, locality, decimalLatitude, decimalLongitude, coordinateUncertaintyInMeters, eventDate, day, month, year) %>% 
  mutate(Spp_code = "AGHY")
ELVI_BRIT <- read_csv(file = "~/Dropbox/Josh&Tom - shared/Endo_Herbarium/DigitizedHerbariumRecords/BRIT_records/BRIT_ELVI_TORCH_records/SymbOutput_2020-04-05_163112_DwC-A/occurrences.csv") %>% 
  filter(!is.na(county), !is.na(eventDate)) %>% 
  dplyr::select(id, catalogNumber, country, stateProvince, county, municipality, locality, decimalLatitude, decimalLongitude, coordinateUncertaintyInMeters, eventDate, day, month, year) %>% 
  mutate(Spp_code = "ELVI")


BRIT_torch <- rbind(AGHY_BRIT, ELVI_BRIT)



# Read in our transcribed datasheet from google sheets
# I am reading these in as csv files that I saved from the excel file because excel does weird things with the date entries.
specimen_info <- read_csv(file = "~joshuacfowler/Dropbox/Josh&Tom - shared/Endo_Herbarium/Endo_Herbarium_specimen.csv") %>% 
  mutate(eventDate = Date_collected) %>%
  mutate(TEXT_Date_collected = gsub("'", "", TEXT_Date_collected)) %>%  #TEXT_Date_collected is saved with a starting ' to preserve the date as text while saving from excel, so I remove that here
  separate(TEXT_Date_collected, into = c("month", "day", "year"), remove = FALSE) %>% 
  mutate(year = as.numeric(year), month = as.numeric(month), day = as.numeric(day))%>% 
  mutate(new_id = case_when(grepl("LL00", Institution_specimen_id) ~ gsub("[a-zA-Z ]", "", Institution_specimen_id),
                            !grepl("LL00", Institution_specimen_id) ~ Institution_specimen_id)) %>% 
  filter(!is.na(Specimen_id))

# This is the sample info and we will filter for only those that we have scored so far.
sample_info <- read_csv(file = "~joshuacfowler/Dropbox/Josh&Tom - shared/Endo_Herbarium/Endo_Herbarium_sample.csv") %>%  
  filter(!is.na(Endo_status_liberal), !is.na(Specimen_id)) %>% 
  mutate(Specimen_id_temp = Specimen_id) %>% 
  separate(Specimen_id_temp, c("Herbarium_id", "Spp_code", "Specimen_no"), "_")



endo_herb <- specimen_info %>% 
  merge(sample_info, by = c("Specimen_id" = "Specimen_id")) %>% 
  select(-contains("X1"),-contains("X2"))




# These are the matches from BRIT for transcription using the file that Jason Best shared of their database, Mar 25th, 2020
BRIT_matches <- read_csv(file = "~/Dropbox/Josh&Tom - shared/Endo_Herbarium/DigitizedHerbariumRecords/BRIT_records/id matches for Torch transcription/fowler_matches_Mar25.csv") %>% 
  select("catalogNumber") %>% 
  mutate(transcibed = NA)


data_matches <- endo_herb %>%
  filter(Herbarium_id == "BRIT") %>% 
  rename( catalogNumber = Institution_specimen_id)

matches <- BRIT_matches %>% 
  merge(data_matches, by = "catalogNumber")
# write_csv(matches, path = "~/Dropbox/Josh&Tom - shared/Endo_Herbarium/DigitizedHerbariumRecords/BRIT_records/id matches for Torch transcription/scored_BRIT_matches_Mar25.csv")


# Now we merge all of our records together with the endo_herb database

# Merge in the AM_records that we have so far
endo_herb1 <- endo_herb %>% 
  left_join(AM_records, by = c( "new_id" = "Institution_specimen_id")) %>% 
  mutate(County = case_when(is.na(County) ~ county,
                               is.na(county) ~ County)) %>% 
  mutate(State = case_when(is.na(State) ~ state,
                               is.na(state) ~ State)) %>% 
  mutate(Country = case_when(is.na(Country.x) ~ Country.y,
                               is.na(Country.y) ~ Country.x)) %>% 
  mutate(Municipality = case_when(is.na(Municipality) ~ Municpality,
                                      is.na(Municpality) ~ Municipality)) %>% 
  mutate(Locality = case_when(is.na(Locality_info) ~ Locality,
                                  is.na(Locality) ~ Locality_info)) %>% 
  mutate(year = case_when(is.na(year.x) ~ year.y,
                                       is.na(year.y) ~ year.x)) %>% 
  mutate(month = case_when(is.na(month.x) ~ month.y,
                          is.na(month.y) ~ month.x)) %>% 
  mutate(day = case_when(is.na(day.x) ~ day.y,
                          is.na(day.y) ~ day.x)) %>% 
  select(Sample_id, Institution_specimen_id, new_id, Country, State, County, Municipality, Locality, year, month, day, tissue_type, seed_scored, seed_eplus, Endo_status_liberal, Endo_status_conservative)

# Merge in the BRIT records that we have so far
endo_herb2 <- endo_herb1 %>% 
  left_join(BRIT_torch, by = c("new_id" = "catalogNumber")) %>% 
  mutate(County = case_when(is.na(County) ~ county,
                            is.na(county) ~ County)) %>% 
  mutate(State = case_when(is.na(State) ~ stateProvince,
                           is.na(stateProvince) ~ State)) %>% 
  mutate(Country = case_when(is.na(Country) ~ country,
                             is.na(country) ~ Country)) %>% 
  mutate(Municipality = case_when(is.na(Municipality) ~ municipality,
                                  is.na(municipality) ~ Municipality)) %>% 
  mutate(Locality = case_when(is.na(Locality) ~ locality,
                              is.na(locality) ~ Locality)) %>% 
  mutate(year = case_when(is.na(year.x) ~ year.y,
                          is.na(year.y) ~ year.x)) %>% 
  mutate(month = case_when(is.na(month.x) ~ month.y,
                           is.na(month.y) ~ month.x)) %>% 
  mutate(day = case_when(is.na(day.x) ~ day.y,
                         is.na(day.y) ~ day.x)) %>% 
  select(Sample_id, Institution_specimen_id, new_id, Country, State, County, Municipality, Locality, year, month, day, tissue_type, seed_scored, seed_eplus, Endo_status_liberal, Endo_status_conservative)

# Merge in the UT Austin records that we have so far
endo_herb3 <- endo_herb2 %>% 
  left_join(UTAustin_torch, by = c("new_id" = "new_id")) %>% 
  mutate(County = case_when(is.na(County) ~ county,
                            is.na(county) ~ County)) %>% 
  mutate(State = case_when(is.na(State) ~ stateProvince,
                           is.na(stateProvince) ~ State)) %>% 
  mutate(Country = case_when(is.na(Country) ~ country,
                             is.na(country) ~ Country)) %>% 
  mutate(Municipality = case_when(is.na(Municipality) ~ municipality,
                                  is.na(municipality) ~ Municipality))  %>% 
  mutate(Locality = case_when(is.na(Locality) ~ locality,
                              is.na(locality) ~ Locality)) %>% 
  mutate(year = case_when(is.na(year.x) ~ year.y,
                          is.na(year.y) ~ year.x)) %>% 
  mutate(month = case_when(is.na(month.x) ~ month.y,
                           is.na(month.y) ~ month.x)) %>% 
  mutate(day = case_when(is.na(day.x) ~ day.y,
                         is.na(day.y) ~ day.x)) %>% 
  select(Sample_id, Institution_specimen_id, new_id, Country, State, County, Municipality, Locality, year, month, day, tissue_type, seed_scored, seed_eplus, Endo_status_liberal, Endo_status_conservative)


# Now I am going to link these county/locality records to a gps point with ggmap
# This requires and API key which you can set up through google, look at ?register_google. 
# There are restrictions to the total number of queries that you can do per day and month, and if you go over, it costs money, so we will save the output. I believe we have a free trial for year.
# One other note, is that we are only using the level of county/city/state which I believe should be pretty accurate through google. I'm not sure that it could accurately do a more detailed locality string
# I have found software that could do this, but would require a human to double check the output.

endo_herb_georef <-endo_herb3 %>% 
  unite("location_string" , sep = ", " , Municipality,County,State,Country, remove = FALSE, na.rm = TRUE) %>% 
  filter(Endo_status_liberal <= 1) %>% 
  mutate_geocode(location_string)
# write_csv(endo_herb_georef, path = "~/Dropbox/Josh&Tom - shared/Endo_Herbarium/DigitizedHerbariumRecords/endo_herb_georef.csv")
endo_herb_georef <- read_csv(file = "~/Dropbox/Josh&Tom - shared/Endo_Herbarium/DigitizedHerbariumRecords/endo_herb_georef.csv") %>% 
  filter(Country != "Canada")


# Now we can explore the data
plot(endo_herb_georef$lon, endo_herb_georef$lat)
hist(endo_herb_georef$year)
hist(endo_herb_georef$lon)

long_date_mod <- glm(Endo_status_liberal == 1 ~ lon*year, data = endo_herb_georef, family = binomial)

plot(endo_herb_georef$lon, endo_herb_georef$Endo_status_liberal)
plot(endo_herb_georef$lat, endo_herb_georef$Endo_status_liberal)
plot(endo_herb_georef$year, endo_herb_georef$Endo_status_liberal)


newdat1900 <- data.frame(lon = seq(-120,-60,1), year = 1900)
newdat1950 <- data.frame(lon = seq(-120,-60,1), year = 1950)
newdat2000 <- data.frame(lon = seq(-120,-60,1), year = 2000)

y_pred1900 <- predict(long_date_mod, newdata = newdat1900, type = "response")
y_pred1950 <- predict(long_date_mod, newdata = newdat1950, type = "response")
y_pred2000 <- predict(long_date_mod, newdata = newdat2000, type = "response")

plot(endo_herb_georef$lon, endo_herb_georef$Endo_status_liberal)
  lines(newdat1900$lon, y_pred1900, col = "red3")
  lines(newdat1950$lon, y_pred1950, col = "gray39")
  lines(newdat1950$lon, y_pred2000, col = "royalblue3")
  
  






count_of_complete_records <- endo_herb_georef %>% 
  filter(!is.na(lon), !is.na(year), !is.na(Endo_status_liberal))
dim(count_of_complete_records)



#  Messing with Lani's data
  
lani_endo <- read_csv(file = "~/Downloads/endonew.csv") %>% 
  filter(`01liberal` <=1)
plot(lani_endo$decimalLongitude, lani_endo$`01liberal`)
long_date_mod <- glm(`01liberal` == 1 ~ decimalLongitude*Date, data = lani_endo, family = binomial)



ggplot(data = lani_endo)+
  geom_point(aes(x = decimalLongitude, y = `01liberal`)) +
  stat_smooth(aes(x = decimalLongitude, y = `01liberal`), method = glm, method.args = list(family = "binomial"))
  
ggplot(data = lani_endo)+
  geom_point(aes(x = decimalLatitude, y = `01liberal`)) +
  stat_smooth(aes(x = decimalLatitude, y = `01liberal`), method = glm, method.args = list(family = "binomial"))


ggplot(data = lani_endo)+
  geom_point(aes(x = Date, y = `01liberal`)) +
  stat_smooth(aes(x = Date, y = `01liberal`), method = glm, method.args = list(family = "binomial"))

# 





