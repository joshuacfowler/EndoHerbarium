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
# UTAustin
AGHY_UTAustin <- read_csv(file = "~/Dropbox/Josh&Tom - shared/Endo_Herbarium/DigitizedHerbariumRecords/UTAustinAGHYrecords/occurrences.csv") %>% 
  mutate(new_id = gsub("[a-zA-Z ]", "", catalogNumber)) %>%    #creates a new id that we will use to merge with the 
  dplyr::select(id, catalogNumber, country, stateProvince, county, municipality, locality, decimalLatitude, decimalLongitude, coordinateUncertaintyInMeters, new_id, eventDate, day, month, year) %>% 
  mutate(Spp_code = "AGHY")
ELVI_UTAustin <- read_csv(file = "~/Dropbox/Josh&Tom - shared/Endo_Herbarium/DigitizedHerbariumRecords/UTAustinELVIrecords/occurrences.csv") %>% 
  mutate(new_id = gsub("[a-zA-Z ]", "", catalogNumber)) %>%    #creates a new id that we will use to merge with the 
  dplyr::select(id, catalogNumber, country, stateProvince, county, municipality, locality, decimalLatitude, decimalLongitude, coordinateUncertaintyInMeters, new_id, eventDate, day, month, year) %>% 
  mutate(Spp_code = "ELVI")

# Texas A&M digitized records
AM_records <- read_xlsx(path = "~/Dropbox/Josh&Tom - shared/Endo_Herbarium/DigitizedHerbariumRecords/TexasA&M_digitized_records/Fowler Data.xlsx")



# Read in our transcribed datasheet
# I am reading these in as csv files that I saved from the excel file because excel does weird things with the date entries.
specimen_info <- read_csv(file = "~joshuacfowler/Dropbox/Josh&Tom - shared/Endo_Herbarium/Endo_Herbarium_specimen.csv") %>% 
  mutate(eventDate = Date_collected) %>% 
  separate(Date_collected, into = c("month", "day", "year"),sep = "/") %>% 
  mutate(year = as.numeric(year), month = as.numeric(month), day = as.numeric(day))%>% 
  mutate(new_id = case_when(grepl("LL00", Institution_specimen_id) ~ gsub("[a-zA-Z ]", "", Institution_specimen_id),
                            !grepl("LL00", Institution_specimen_id) ~ Institution_specimen_id)) %>% 
  filter(!is.na(Specimen_id))
  
sample_info <- read_csv(file = "~joshuacfowler/Dropbox/Josh&Tom - shared/Endo_Herbarium/Endo_Herbarium_sample.csv") %>%  
  filter(!is.na(Endo_status_liberal), !is.na(Specimen_id)) %>% 
  mutate(Specimen_id_temp = Specimen_id) %>% 
  separate(Specimen_id_temp, c("Herbarium_id", "Spp_code", "Specimen_no"), "_")



data <- specimen_info %>% 
  merge(sample_info, by = c("Specimen_id" = "Specimen_id"))




# These are the matches from BRIT for transcription, Mar 25th, 2020
BRIT_matches <- read_csv(file = "~/Dropbox/Josh&Tom - shared/Endo_Herbarium/DigitizedHerbariumRecords/BRIT_catalogNumberMatches/fowler_matches_Mar25.csv") %>% 
  select("catalogNumber") %>% 
  mutate(transcibed = NA)


data_matches <- data %>%
  filter(Herbarium_id == "BRIT") %>% 
  rename( catalogNumber = Institution_specimen_id)

matches <- BRIT_matches %>% 
  merge(data_matches, by = "catalogNumber")
write_csv(matches, path = "~/Dropbox/Josh&Tom - shared/Endo_Herbarium/DigitizedHerbariumRecords/BRIT_catalogNumberMatches/scored_BRIT_matches_Mar25.csv")



# we can see what id's aren't found in both datasets
intersect(AGHY_UTAustin$new_id, specimen_info$new_id)
setdiff(AGHY_UTAustin$new_id, specimen_info$new_id)

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

# now we can merge in the specimen and the sample info based on the Specimen id
AGHY_data <- AGHY_specimen_info %>% 
  left_join(sample_info, by = c("Specimen_id" = "Specimen_id")) %>% 
  dplyr::select(Specimen_id, YEAR, MONTH, DAY, COUNTRY, STATE, COUNTY, MUNICIPALITY, LOCALITY, Endo_status_liberal, Endo_status_conservative) %>% 
  filter(!is.na(Endo_status_conservative)) %>% 
  mutate(DATE = as.Date(ymd(paste(YEAR, MONTH, DAY, sep = "-"))))
write_csv(AGHY_data, path = "AGHY_data_forgeoreference.csv")





# Simple model to look at date effects on infections
model <-  glm(Endo_status_liberal ~ YEAR, data = AGHY_data, family = "binomial")

yrep <- predict(model, list(YEAR = AGHY_data$YEAR), type = "response" )

plot(AGHY_data$YEAR, AGHY_data$Endo_status_liberal, pch = 16, xlab = "WEIGHT (g)", ylab = "VS")
lines(AGHY_data$YEAR, yrep)
# Not significant, but maybe a slight decrease

model <-  glm(Endo_status_conservative ~ YEAR, data = AGHY_data, family = "binomial")

yrep <- predict(model, list(YEAR = AGHY_data$YEAR), type = "response" )

plot(AGHY_data$YEAR, AGHY_data$Endo_status_conservative, pch = 16, xlab = "WEIGHT (g)", ylab = "VS")
lines(AGHY_data$YEAR, yrep)





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





