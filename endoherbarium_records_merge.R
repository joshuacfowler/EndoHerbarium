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



endo_herbarium <- specimen_info %>% 
  merge(sample_info, by = c("Specimen_id" = "Specimen_id")) %>% 
  select(-contains("X1"),-contains("X2"))




# These are the matches from BRIT for transcription using the file that Jason Best shared of their database, Mar 25th, 2020
BRIT_matches <- read_csv(file = "~/Dropbox/Josh&Tom - shared/Endo_Herbarium/DigitizedHerbariumRecords/BRIT_records/id matches for Torch transcription/fowler_matches_Mar25.csv") %>% 
  select("catalogNumber") %>% 
  mutate(transcibed = NA)


data_matches <- endo_herbarium %>%
  filter(Herbarium_id == "BRIT") %>% 
  rename( catalogNumber = Institution_specimen_id)

matches <- BRIT_matches %>% 
  merge(data_matches, by = "catalogNumber")
# write_csv(matches, path = "~/Dropbox/Josh&Tom - shared/Endo_Herbarium/DigitizedHerbariumRecords/BRIT_records/id matches for Torch transcription/scored_BRIT_matches_Mar25.csv")


# Now we merge all of our records together with the endo_herbarium database

# Merge in the AM_records that we have so far
endo_herbarium1 <- endo_herbarium %>% 
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
endo_herbarium2 <- endo_herbarium1 %>% 
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
endo_herbarium3 <- endo_herbarium2 %>% 
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
endo_herbarium_georeferenced <- mutate_geocode(endo_herbarium3, county)




# I'm gonna filter for just AGHY to be able to merge this with the digitized records for AGHY
# This is still a little bit messy, but it has both datasets merged together. 
# The next step would be to actually ccombine some of the columns in each dataset that are the same. Sometimes this cause issues if there are conflicting values in one or the other.
# You can either add them to the merge and if they are distinct, they will fit together nicely, or if you can create new columns by combining the two existing columns.
AGHY_specimen_info <- specimen_info %>% 
  left_join(AGHY_UTAustin, by = c("new_id" = "new_id")) %>% 
  mutate(SPP_CODE = case_when(is.na(Spp_code.x) ~ Spp_code.y,
                              !is.na(Spp_code.x) ~ Spp_code.x),
         COUNTRY = case_when(is.na(Country) ~ country,
                            !is.na(Country) ~ Country),
         STATE = case_when(is.na(State) ~ stateProvince,
                           !is.na(State) ~ State),
         COUNTY = case_when(is.na(county) ~ County,
                            !is.na(county) ~ county),
         MUNICIPALITY = Municipality,
         LOCALITY = case_when(is.na(Locality_info) ~ locality,
                                  !is.na(Locality_info) ~ Locality_info),
         YEAR = case_when(is.na(year.x) ~ year.y,
                          !is.na(year.x) ~ year.x),
         MONTH = case_when(is.na(month.x) ~ month.y,
                          !is.na(month.x) ~ month.x),
         DAY = case_when(is.na(day.x) ~ day.y,
                          !is.na(day.x) ~ day.x)) %>% 
  dplyr::select(Specimen_id, Institution_specimen_id, catalogNumber, new_id, SPP_CODE, COUNTRY, STATE, COUNTY, MUNICIPALITY, LOCALITY, YEAR, MONTH, DAY ) %>% 
  filter(SPP_CODE == "AGHY")
  

data <- specimen_info %>% 
  full_join(sample_info, by = c("Specimen_id" = "Specimen_id")) %>% 
  filter(Specimen_id == contains("BRIT"))
  group_by(newSPP)
  summarise(n())





model <-  glm(Endo_status_liberal ~ YEAR, data = AGHY_data, family = "binomial")



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





