# Title: Digital Herbarium Records and Sample Integration
# Purpose: Imports and merges downloaded digitized herbarium records with the Endo_Herbarium database
# Authors: Joshua Fowler and Lani Dufresne
# Updated: Nov. 10, 2020

library(tidyverse)
library(readxl)
library(ggmap)
library(lubridate)
library(rstan)
library(car)
library(sf) # for making maps
library(here) # for making maps (lets you set filepaths)
library(ggmap) # for making maps
library(rnaturalearth)# for making maps
library(rnaturalearthdata)# for making maps
library(rgeos)# for making maps
library(prism)
library(raster) ##working with raster data
library(reshape2)
library(viridis)

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
AGPE_UTAustin <-   read_csv(file = "~/Dropbox/Josh&Tom - shared/Endo_Herbarium/DigitizedHerbariumRecords/UTAustinAGPErecords/occurrences.csv") %>% 
  mutate(new_id = gsub("[a-zA-Z ]", "", catalogNumber)) %>%    #creates a new id that we will use to merge with the 
  mutate(municipality = as.character(municipality)) %>% 
  dplyr::select(id, catalogNumber, country, stateProvince, county, municipality, locality, decimalLatitude, decimalLongitude, coordinateUncertaintyInMeters, new_id, eventDate, day, month, year) %>% 
  mutate(Spp_code = "ELVI")

UTAustin_torch <- rbind(AGHY_UTAustin, ELVI_UTAustin, AGPE_UTAustin)

# Texas A&M digitized records (Includes both AGHY and ELVI)
AM_records <- read_csv(file = "~/Dropbox/Josh&Tom - shared/Endo_Herbarium/DigitizedHerbariumRecords/TexasA&M_digitized_records/Fowler Data.csv") %>% 
  unite(Institution_specimen_id, CollectionCode, id, sep = "") %>%
  separate(DateCollected, into = c("year", "month", "day"), remove = FALSE) %>% 
  mutate(year = as.numeric(year), month = as.numeric(month), day = as.numeric(day)) %>% 
  select(-contains("X"))


# BRIT digitized records downloaded from TORCH (includes Vanderbilt and U of Louisiana Monroe) 
# This was downloaded on Jul17 and we can get more specimens transcribed and download again.
AGHY_BRIT <- read_csv(file = "~/Dropbox/Josh&Tom - shared/Endo_Herbarium/DigitizedHerbariumRecords/BRIT_records/BRIT_AGHY_TORCH_records/SymbOutput_2020-11-10_144536_DwC-A/occurrences.csv",
                      col_types = cols(id = col_double(),
                                       institutionCode = col_character(),
                                       collectionCode = col_logical(),
                                       ownerInstitutionCode = col_logical(),
                                       basisOfRecord = col_character(),
                                       occurrenceID = col_character(),
                                       catalogNumber = col_character(),
                                       otherCatalogNumbers = col_character(),
                                       higherClassification = col_character(),
                                       kingdom = col_character(),
                                       phylum = col_character(),
                                       class = col_logical(),
                                       order = col_character(),
                                       family = col_character(),
                                       scientificName = col_character(),
                                       taxonID = col_double(),
                                       scientificNameAuthorship = col_character(),
                                       genus = col_character(),
                                       subgenus = col_logical(),
                                       specificEpithet = col_character(),
                                       verbatimTaxonRank = col_character(),
                                       infraspecificEpithet = col_character(),
                                       taxonRank = col_character(),
                                       identifiedBy = col_character(),
                                       dateIdentified = col_character(),
                                       identificationReferences = col_character(),
                                       identificationRemarks = col_character(),
                                       taxonRemarks = col_character(),
                                       identificationQualifier = col_character(),
                                       typeStatus = col_logical(),
                                       recordedBy = col_character(),
                                       associatedCollectors = col_character(),
                                       recordNumber = col_character(),
                                       eventDate = col_character(),
                                       year = col_double(),
                                       month = col_double(),
                                       day = col_double(),
                                       startDayOfYear = col_double(),
                                       endDayOfYear = col_logical(),
                                       verbatimEventDate = col_character(),
                                       occurrenceRemarks = col_character(),
                                       habitat = col_character(),
                                       substrate = col_character(),
                                       verbatimAttributes = col_character(),
                                       fieldNumber = col_logical(),
                                       informationWithheld = col_logical(),
                                       dataGeneralizations = col_logical(),
                                       dynamicProperties = col_character(),
                                       associatedTaxa = col_character(),
                                       reproductiveCondition = col_character(),
                                       establishmentMeans = col_character(),
                                       cultivationStatus = col_logical(),
                                       lifeStage = col_logical(),
                                       sex = col_logical(),
                                       individualCount = col_logical(),
                                       preparations = col_logical(),
                                       country = col_character(),
                                       stateProvince = col_character(),
                                       county = col_character(),
                                       municipality = col_character(),
                                       locality = col_character(),
                                       locationRemarks = col_logical(),
                                       localitySecurity = col_double(),
                                       localitySecurityReason = col_logical(),
                                       decimalLatitude = col_double(),
                                       decimalLongitude = col_double(),
                                       geodeticDatum = col_character(),
                                       coordinateUncertaintyInMeters = col_double(),
                                       verbatimCoordinates = col_character(),
                                       georeferencedBy = col_character(),
                                       georeferenceProtocol = col_character(),
                                       georeferenceSources = col_character(),
                                       georeferenceVerificationStatus = col_character(),
                                       georeferenceRemarks = col_character(),
                                       minimumElevationInMeters = col_double(),
                                       maximumElevationInMeters = col_double(),
                                       minimumDepthInMeters = col_logical(),
                                       maximumDepthInMeters = col_logical(),
                                       verbatimDepth = col_logical(),
                                       verbatimElevation = col_character(),
                                       disposition = col_logical(),
                                       language = col_logical(),
                                       recordEnteredBy = col_character(),
                                       modified = col_datetime(format = ""),
                                       `sourcePrimaryKey-dbpk` = col_logical(),
                                       collId = col_double(),
                                       recordId = col_character(),
                                       references = col_character())) %>% 
  filter(!is.na(county), !is.na(eventDate)) %>% 
  separate(eventDate, into = c("year", "month", "day"), remove = FALSE) %>% 
  mutate_at(c("day", "month", "year"), as.numeric) %>% 
  dplyr::select(id, catalogNumber, country, stateProvince, county, municipality, locality, decimalLatitude, decimalLongitude, coordinateUncertaintyInMeters, eventDate, day, month, year) %>% 
  mutate(Spp_code = "AGHY")

ELVI_BRIT <- read_csv(file = "~/Dropbox/Josh&Tom - shared/Endo_Herbarium/DigitizedHerbariumRecords/BRIT_records/BRIT_ELVI_TORCH_records/SymbOutput_2020-11-10_144847_DwC-A/occurrences.csv",
                      col_types = cols(id = col_double(),
                                       institutionCode = col_character(),
                                       collectionCode = col_logical(),
                                       ownerInstitutionCode = col_logical(),
                                       basisOfRecord = col_character(),
                                       occurrenceID = col_character(),
                                       catalogNumber = col_character(),
                                       otherCatalogNumbers = col_character(),
                                       higherClassification = col_character(),
                                       kingdom = col_character(),
                                       phylum = col_character(),
                                       class = col_logical(),
                                       order = col_character(),
                                       family = col_character(),
                                       scientificName = col_character(),
                                       taxonID = col_double(),
                                       scientificNameAuthorship = col_character(),
                                       genus = col_character(),
                                       subgenus = col_logical(),
                                       specificEpithet = col_character(),
                                       verbatimTaxonRank = col_character(),
                                       infraspecificEpithet = col_character(),
                                       taxonRank = col_character(),
                                       identifiedBy = col_character(),
                                       dateIdentified = col_character(),
                                       identificationReferences = col_character(),
                                       identificationRemarks = col_character(),
                                       taxonRemarks = col_character(),
                                       identificationQualifier = col_character(),
                                       typeStatus = col_logical(),
                                       recordedBy = col_character(),
                                       associatedCollectors = col_character(),
                                       recordNumber = col_character(),
                                       eventDate = col_character(),
                                       year = col_double(),
                                       month = col_double(),
                                       day = col_double(),
                                       startDayOfYear = col_double(),
                                       endDayOfYear = col_logical(),
                                       verbatimEventDate = col_character(),
                                       occurrenceRemarks = col_character(),
                                       habitat = col_character(),
                                       substrate = col_character(),
                                       verbatimAttributes = col_character(),
                                       fieldNumber = col_logical(),
                                       informationWithheld = col_logical(),
                                       dataGeneralizations = col_logical(),
                                       dynamicProperties = col_character(),
                                       associatedTaxa = col_character(),
                                       reproductiveCondition = col_character(),
                                       establishmentMeans = col_character(),
                                       cultivationStatus = col_logical(),
                                       lifeStage = col_logical(),
                                       sex = col_logical(),
                                       individualCount = col_logical(),
                                       preparations = col_logical(),
                                       country = col_character(),
                                       stateProvince = col_character(),
                                       county = col_character(),
                                       municipality = col_character(),
                                       locality = col_character(),
                                       locationRemarks = col_logical(),
                                       localitySecurity = col_double(),
                                       localitySecurityReason = col_logical(),
                                       decimalLatitude = col_double(),
                                       decimalLongitude = col_double(),
                                       geodeticDatum = col_character(),
                                       coordinateUncertaintyInMeters = col_double(),
                                       verbatimCoordinates = col_character(),
                                       georeferencedBy = col_character(),
                                       georeferenceProtocol = col_character(),
                                       georeferenceSources = col_character(),
                                       georeferenceVerificationStatus = col_character(),
                                       georeferenceRemarks = col_character(),
                                       minimumElevationInMeters = col_double(),
                                       maximumElevationInMeters = col_double(),
                                       minimumDepthInMeters = col_logical(),
                                       maximumDepthInMeters = col_logical(),
                                       verbatimDepth = col_logical(),
                                       verbatimElevation = col_character(),
                                       disposition = col_logical(),
                                       language = col_logical(),
                                       recordEnteredBy = col_character(),
                                       modified = col_datetime(format = ""),
                                       `sourcePrimaryKey-dbpk` = col_logical(),
                                       collId = col_double(),
                                       recordId = col_character(),
                                       references = col_character())) %>% 
  filter(!is.na(county), !is.na(eventDate)) %>% 
  separate(eventDate, into = c("year", "month", "day"), remove = FALSE) %>% 
  mutate_at(c("day", "month", "year"), as.numeric) %>% 
  dplyr::select(id, catalogNumber, country, stateProvince, county, municipality, locality, decimalLatitude, decimalLongitude, coordinateUncertaintyInMeters, eventDate, day, month, year) %>% 
  mutate(Spp_code = "ELVI")


AGPE_BRIT <- read_csv(file = "~/Dropbox/Josh&Tom - shared/Endo_Herbarium/DigitizedHerbariumRecords/BRIT_records/BRIT_AGPE_TORCH_records/SymbOutput_2020-11-10_144221_DwC-A/occurrences.csv",
                      col_types = cols(id = col_double(),
                                       institutionCode = col_character(),
                                       collectionCode = col_logical(),
                                       ownerInstitutionCode = col_logical(),
                                       basisOfRecord = col_character(),
                                       occurrenceID = col_character(),
                                       catalogNumber = col_character(),
                                       otherCatalogNumbers = col_character(),
                                       higherClassification = col_character(),
                                       kingdom = col_character(),
                                       phylum = col_character(),
                                       class = col_logical(),
                                       order = col_character(),
                                       family = col_character(),
                                       scientificName = col_character(),
                                       taxonID = col_double(),
                                       scientificNameAuthorship = col_character(),
                                       genus = col_character(),
                                       subgenus = col_logical(),
                                       specificEpithet = col_character(),
                                       verbatimTaxonRank = col_character(),
                                       infraspecificEpithet = col_character(),
                                       taxonRank = col_character(),
                                       identifiedBy = col_character(),
                                       dateIdentified = col_character(),
                                       identificationReferences = col_character(),
                                       identificationRemarks = col_character(),
                                       taxonRemarks = col_character(),
                                       identificationQualifier = col_character(),
                                       typeStatus = col_logical(),
                                       recordedBy = col_character(),
                                       associatedCollectors = col_character(),
                                       recordNumber = col_character(),
                                       eventDate = col_character(),
                                       year = col_double(),
                                       month = col_double(),
                                       day = col_double(),
                                       startDayOfYear = col_double(),
                                       endDayOfYear = col_logical(),
                                       verbatimEventDate = col_character(),
                                       occurrenceRemarks = col_character(),
                                       habitat = col_character(),
                                       substrate = col_character(),
                                       verbatimAttributes = col_character(),
                                       fieldNumber = col_logical(),
                                       informationWithheld = col_logical(),
                                       dataGeneralizations = col_logical(),
                                       dynamicProperties = col_character(),
                                       associatedTaxa = col_character(),
                                       reproductiveCondition = col_character(),
                                       establishmentMeans = col_character(),
                                       cultivationStatus = col_logical(),
                                       lifeStage = col_logical(),
                                       sex = col_logical(),
                                       individualCount = col_logical(),
                                       preparations = col_logical(),
                                       country = col_character(),
                                       stateProvince = col_character(),
                                       county = col_character(),
                                       municipality = col_character(),
                                       locality = col_character(),
                                       locationRemarks = col_logical(),
                                       localitySecurity = col_double(),
                                       localitySecurityReason = col_logical(),
                                       decimalLatitude = col_double(),
                                       decimalLongitude = col_double(),
                                       geodeticDatum = col_character(),
                                       coordinateUncertaintyInMeters = col_double(),
                                       verbatimCoordinates = col_character(),
                                       georeferencedBy = col_character(),
                                       georeferenceProtocol = col_character(),
                                       georeferenceSources = col_character(),
                                       georeferenceVerificationStatus = col_character(),
                                       georeferenceRemarks = col_character(),
                                       minimumElevationInMeters = col_double(),
                                       maximumElevationInMeters = col_double(),
                                       minimumDepthInMeters = col_logical(),
                                       maximumDepthInMeters = col_logical(),
                                       verbatimDepth = col_logical(),
                                       verbatimElevation = col_character(),
                                       disposition = col_logical(),
                                       language = col_logical(),
                                       recordEnteredBy = col_character(),
                                       modified = col_datetime(format = ""),
                                       `sourcePrimaryKey-dbpk` = col_logical(),
                                       collId = col_double(),
                                       recordId = col_character(),
                                       references = col_character())) %>% 
  filter(!is.na(county), !is.na(eventDate)) %>% 
  separate(eventDate, into = c("year", "month", "day"), remove = FALSE) %>% 
  mutate_at(c("day", "month", "year"), as.numeric) %>% 
  dplyr::select(id, catalogNumber, country, stateProvince, county, municipality, locality, decimalLatitude, decimalLongitude, coordinateUncertaintyInMeters, eventDate, day, month, year) %>% 
  mutate(Spp_code = "AGPE")


BRIT_torch <- rbind(AGHY_BRIT, ELVI_BRIT, AGPE_BRIT)



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
  mutate(Spp_code = case_when(is.na(Spp_code.x) ~ Spp_code.y,
                            is.na(Spp_code.y) ~ Spp_code.x)) %>% 
  dplyr::select(-contains("X1"),-contains("X2"))




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
  mutate(Spp_code = case_when(is.na(Spp_code.x) ~ Spp_code.y,
                              is.na(Spp_code.y) ~ Spp_code.x)) %>% 
  dplyr::select(Sample_id, Institution_specimen_id, Spp_code, new_id, Country, State, County, Municipality, Locality, year, month, day, tissue_type, seed_scored, seed_eplus, Endo_status_liberal, Endo_status_conservative)

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
  mutate(Spp_code = case_when(is.na(Spp_code.x) ~ Spp_code.y,
                            is.na(Spp_code.y) ~ Spp_code.x)) %>% 
  dplyr::select(Sample_id, Institution_specimen_id, Spp_code, new_id, Country, State, County, Municipality, Locality, year, month, day, tissue_type, seed_scored, seed_eplus, Endo_status_liberal, Endo_status_conservative)

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
  mutate(Spp_code = case_when(is.na(Spp_code.x) ~ Spp_code.y,
                              is.na(Spp_code.y) ~ Spp_code.x)) %>% 
  dplyr::select(Sample_id, Institution_specimen_id, Spp_code, new_id, Country, State, County, Municipality, Locality, year, month, day, tissue_type, seed_scored, seed_eplus, Endo_status_liberal, Endo_status_conservative)


# Now I am going to link these county/locality records to a gps point with ggmap
# This requires and API key which you can set up through google, look at ?register_google. 
# There are restrictions to the total number of queries that you can do per day and month, and if you go over, it costs money, so we will save the output. I believe we have a free trial for year.
# One other note, is that we are only using the level of county/city/state which I believe should be pretty accurate through google. I'm not sure that it could accurately do a more detailed locality string
# I have found software that could do this, but would require a human to double check the output.

# endo_herb_georef <-endo_herb3 %>%
#   unite("location_string" , sep = ", " , Municipality,County,State,Country, remove = FALSE, na.rm = TRUE) %>%
#   filter(Endo_status_liberal <= 1) %>%
# mutate_geocode(location_string) # Uncomment this to run the geocoding.
# write_csv(endo_herb_georef, path = "~/Dropbox/Josh&Tom - shared/Endo_Herbarium/DigitizedHerbariumRecords/endo_herb_georef.csv")
endo_herb_georef <- read_csv(file = "~/Dropbox/Josh&Tom - shared/Endo_Herbarium/DigitizedHerbariumRecords/endo_herb_georef.csv") %>% 
  filter(Country != "Canada")


# Now we can explore the data
plot(endo_herb_georef$lon, endo_herb_georef$lat)
plot(endo_herb_georef$lon, endo_herb_georef$Endo_status_liberal)
plot(endo_herb_georef$lat, endo_herb_georef$Endo_status_liberal)
plot(endo_herb_georef$year, endo_herb_georef$Endo_status_liberal)

# counts of scores
endo_herb_AGHY <- endo_herb_georef %>% 
  filter(grepl("AGHY", Sample_id)) %>% 
  filter(!is.na(lon) & !is.na(year))

 binned_AGHY <- endo_herb_AGHY %>% 
   mutate(binned_lon = cut(lon, breaks = 12), binned_year = cut(year, breaks = 3)) %>%  
   group_by(binned_lon, binned_year) %>%   
   summarise(mean_lon = mean(lon),
             mean_year = mean(year),
             mean_endo = mean(Endo_status_liberal),
             sample = n())

endo_herb_ELVI <- endo_herb_georef %>% 
  filter(grepl("ELVI", Sample_id)) %>% 
  filter(!is.na(lon) & !is.na(year))

binned_ELVI <- endo_herb_ELVI %>% 
  mutate(binned_lon = cut(lon, breaks = 12), binned_year = cut(year, breaks = 3)) %>%  
  group_by(binned_lon, binned_year) %>%   
  summarise(mean_lon = mean(lon),
            mean_year = mean(year),
            mean_endo = mean(Endo_status_liberal),
            sample = n())

# AGHY
plot(endo_herb_AGHY$lon, endo_herb_AGHY$lat)
plot(endo_herb_AGHY$year,endo_herb_AGHY$Endo_status_liberal)
plot(endo_herb_AGHY$lon,endo_herb_AGHY$Endo_status_liberal)
hist(endo_herb_AGHY$year)
hist(endo_herb_AGHY$lon)

long_date_mod <- glm(Endo_status_liberal ~ lon * year, data = subset(endo_herb_AGHY), family = binomial)
long_date_mod <- glm(Endo_status_liberal ~ year * lon, data = subset(endo_herb_AGHY), family = binomial)

anova(long_date_mod, test = "Chisq")
summary(long_date_mod)

Anova(long_date_mod, test = "LR")

newdat1920 <- data.frame(lon = seq(-120,-60,1), year = 1920)
newdat1950 <- data.frame(lon = seq(-120,-60,1), year = 1950)
newdat2000 <- data.frame(lon = seq(-120,-60,1), year = 2000)
newdat <- rbind(newdat1920, newdat1950, newdat2000)
y_pred <- predict(long_date_mod, newdata = newdat, type = "response")
y_CI <- predict(long_date_mod, newdata = newdat, interval = "confidence", type = "link", se.fit=TRUE)
linkinv <- family(long_date_mod)$linkinv ## inverse-link function

newpred <- newdat
newpred$pred0 <- y_CI$fit
newpred$pred <- linkinv(y_CI$fit)
alpha <- 0.95
sc <- abs(qnorm((1-alpha)/2))  ## Normal approx. to likelihood
alpha2 <- 0.5
sc2 <- abs(qnorm((1-alpha2)/2))  ## Normal approx. to likelihood
newpred <- transform(newpred,
                    lwr=linkinv(pred0-sc*y_CI$se.fit),
                    upr=linkinv(pred0+sc*y_CI$se.fit),
                    lwr2=linkinv(pred0-sc2*y_CI$se.fit),
                    upr2=linkinv(pred0+sc2*y_CI$se.fit))



AGHY_herb <- ggplot() +
  geom_point(data = binned_AGHY,aes(x = mean_lon, y = mean_endo, size = sample, color = binned_year)) + 
  geom_line(data = newpred, aes(x = lon, y = pred, group = year, color = as.character(year))) +
  geom_ribbon(data = newpred, aes(x = lon, ymin = lwr, ymax = upr, group = year, fill = as.factor(year)), alpha = .2) +
  theme_classic() + labs(y = "Mean Endophyte Prevalence", x = "Longitude", color = "year")+
  scale_colour_manual(breaks = c("1920", "1950", "2000"),
                      values = c("#fc8d59", "#636363", "#91bfdb", "#fc8d69", "#635363", "#81bfdb")) +
  scale_fill_manual(breaks = c("1920", "1950", "2000"),
                    values = c("#fc8d59", "#636363", "#91bfdb")) +
  xlim(-105,-60) + guides(fill = FALSE)

 AGHY_herb
ggsave(AGHY_herb, filename = "~/Documents/AGHYherb.tiff", width = 4, height = 3)





### Plot for binned over longitude, all times merged

longbin_AGHY <- endo_herb_AGHY %>% 
  mutate(binned_lon = cut(lon, breaks = 25)) %>%  
  group_by(binned_lon) %>%   
  summarise(mean_lon = mean(lon),
            mean_endo = mean(Endo_status_liberal),
            sample = n())

long_mod <- glm(Endo_status_liberal ~ lon, data = subset(endo_herb_AGHY), family = binomial)
summary(long_mod)
anova(long_mod, test = "Chisq")

newdatlong <- data.frame(lon = seq(-110,-60,1))
y_pred <- predict(long_mod, newdata = newdatlong, type = "response")
y_CI <- predict(long_mod, newdata = newdatlong, interval = "confidence", type = "link", se.fit=TRUE)
linkinv <- family(long_mod)$linkinv ## inverse-link function


newpred <- newdatlong
newpred$pred0 <- y_CI$fit
newpred$pred <- linkinv(y_CI$fit)
alpha <- 0.95
sc <- abs(qnorm((1-alpha)/2))  ## Normal approx. to likelihood
alpha2 <- 0.5
sc2 <- abs(qnorm((1-alpha2)/2))  ## Normal approx. to likelihood
newpred <- transform(newpred,
                     lwr=linkinv(pred0-sc*y_CI$se.fit),
                     upr=linkinv(pred0+sc*y_CI$se.fit),
                     lwr2=linkinv(pred0-sc2*y_CI$se.fit),
                     upr2=linkinv(pred0+sc2*y_CI$se.fit))



AGHY_herb_long <- ggplot() +
  geom_point(data = longbin_AGHY,aes(x = mean_lon, y = mean_endo, size = sample)) + 
  geom_line(data = newpred, aes(x = lon, y = pred)) +
  geom_ribbon(data = newpred, aes(x = lon, ymin = lwr, ymax = upr, alpha = .5)) +
  theme_classic() + labs(y = "Mean Endophyte Prevalence", x = "Longitude")+
  xlim(-110,-60) + guides(alpha = FALSE)

AGHY_herb_long
ggsave(AGHY_herb_long, filename = "~/Documents/AGHYherb_long.tiff", width = 4, height = 3)




#### Plot for change over time, binned by longitude

binned_AGHY_time <- endo_herb_AGHY %>% 
  mutate(binned_lon = cut(lon, breaks = 4), binned_year = cut(year, breaks = 12)) %>% 
  group_by(binned_lon, binned_year) %>%   
  summarise(mean_lon = mean(lon),
            mean_year = mean(year),
            mean_endo = mean(Endo_status_liberal),
            sample = n()) %>% 
  mutate(lon = case_when(mean_lon <= -97.3 ~ -102,
                         mean_lon > -97.3 & mean_lon < -88.1 ~ -92,
                         mean_lon > -88.1 & mean_lon < -78.9 ~ -83,
                         mean_lon >= -78.9 ~ -73))

newdat73 <- data.frame(lon = -73, year = seq(1880,2015,1))
newdat83 <- data.frame(lon = -83, year = seq(1880,2015,1))
newdat92 <- data.frame(lon = -92, year = seq(1880,2015,1))
newdat102 <- data.frame(lon = -102, year = seq(1880,2015,1))
newdat <- rbind(newdat73, newdat83,newdat92, newdat102)
y_pred <- predict(long_date_mod, newdata = newdat, type = "response")
y_CI <- predict(long_date_mod, newdata = newdat, interval = "confidence", type = "link", se.fit=TRUE)
linkinv <- family(long_date_mod)$linkinv ## inverse-link function

newpred <- newdat
newpred$pred0 <- y_CI$fit
newpred$pred <- linkinv(y_CI$fit)
alpha <- 0.95
sc <- abs(qnorm((1-alpha)/2))  ## Normal approx. to likelihood
alpha2 <- 0.5
sc2 <- abs(qnorm((1-alpha2)/2))  ## Normal approx. to likelihood
newpred <- transform(newpred,
                     lwr=linkinv(pred0-sc*y_CI$se.fit),
                     upr=linkinv(pred0+sc*y_CI$se.fit),
                     lwr2=linkinv(pred0-sc2*y_CI$se.fit),
                     upr2=linkinv(pred0+sc2*y_CI$se.fit))




AGHY_herb_time <- ggplot() +
  geom_point(data = binned_AGHY_time,aes(x = mean_year, y = mean_endo, size = sample))+
  geom_line(data = newpred, aes(x = year, y = pred,  group = lon)) +
  geom_ribbon(data = newpred, aes(x = year, ymin = lwr, ymax = upr, group = lon), alpha = .2) +
  theme_classic() + labs(y = "Mean Endophyte Prevalence", x = "Year", color = "Longitude")+
  facet_wrap(~lon, nrow = 1)+
  xlim(1880,2015) + guides(fill = FALSE)
AGHY_herb_time
ggsave(AGHY_herb_time, filename = "~/Documents/AGHY_herb_time.tiff", width = 6, height = 3)




# Making a map of scored AGHY
usa <- st_as_sf(maps::map("state", fill=TRUE, plot =FALSE))
AGHY_herb_map <- ggplot()+
  geom_sf(data = usa, fill = "white") +
  geom_point(data = endo_herb_AGHY, aes(x = lon, y = lat), lwd = 2,alpha = .5) +
  theme_minimal() + lims(x = c(-105,-72)) + labs(x = c("Longitude"), y = c("Latitude"))

AGHY_herb_map
ggsave(AGHY_herb_map, filename = "~/Documents/AGHY_herb_map.tiff")



# Get annual prism climate data to make map of climate change magnitude
# Should have these PRISM files saved after running the first time
get_prism_annual(type = c("ppt"), years=1895:2016, keepZip = TRUE)
get_prism_annual(type = c("tmean"), years=1895:2016, keepZip = TRUE)

ls_prism_data(name=TRUE)

new_file<-c(1:(length(1895:2016))) ##change to corresponding file numbers
RS <- prism_stack(ls_prism_data()[new_file,1]) ##raster file
to_slice <- grep("_2016",RS[,1],value=T)##search through names


df <- data.frame(rasterToPoints(RS)) ##creates a dataframe of points
year.df <- melt(df, c("x", "y")) %>% 
  separate(variable, into = c("PRISM", "Variable", "Label", "Resolution", "Year", "file")) %>% 
  rename("lon" = "x", "lat" = "y")


normals <- year.df %>% 
  mutate(time.period = case_when(Year>=1895 & Year <=1956 ~ "start",
                                 Year>=1957 & Year <=2016 ~ "end")) %>% 
  group_by(lon, lat, Variable, time.period) %>% 
  summarize(normal = mean(value)) 


normals_diff <- normals %>% 
  pivot_wider(names_from = "time.period", values_from = "normal") %>% 
  group_by(lon,lat,Variable) %>% 
  mutate(normal_diff = end - start)


ggplot()+
  geom_raster(data = subset(normals_diff, Variable == "tmean"), aes(x = lon, y = lat, fill = normal_diff)) +
  scale_fill_distiller("Change in Mean Temperature (\u00B0C)", palette = "Spectral", limits = c(-1.15,1.15)) +
  theme_minimal() +coord_fixed(ratio=1.3)

ggplot()+
  geom_raster(data = subset(normals_diff, Variable == "ppt"), aes(x = lon, y = lat, fill = normal_diff)) +
  scale_fill_distiller("Change in Annual Precipitation (mm.)", palette = "Spectral", trans = "reverse", limits = c(250,-250)) +
  theme_minimal() +coord_fixed(ratio=1.3)




# Map of collections with change in temp and ppt

tempchangemap <- ggplot()+
  geom_raster(data = subset(normals_diff, Variable == "tmean"), aes(x = lon, y = lat, fill = normal_diff)) +
  geom_sf(data = usa, fill = "transparent") +
  geom_point(data = endo_herb_AGHY, aes(x = lon, y = lat), shape = 4, lwd = 1.5,alpha = .6) +
  theme_minimal() + labs(x = c("Longitude"), y = c("Latitude"), fill = c("Mean Temp. Change (\u00B0C)"))+
  scale_fill_distiller( palette = "Spectral", limits = c(-1.15,1.15)) + lims(x = c(-108,-65))

tempchangemap
ggsave(tempchangemap, filename = "~/Documents/tempchangemap2.tiff")



pptchangemap <- ggplot()+
  geom_raster(data = subset(normals_diff, Variable == "ppt"), aes(x = lon, y = lat, fill = normal_diff)) +
  geom_sf(data = usa, fill = "transparent") +
  geom_point(data = endo_herb_AGHY, aes(x = lon, y = lat), shape = 4, lwd = 1.5,alpha = .6) + 
  theme_minimal() + labs(x = c("Longitude"), y = c("Latitude"), fill = c("Annual Precip. Change (mm.)"))+
  scale_fill_distiller("Annual Precip. Change (mm.)", palette = "Spectral", trans = "reverse", limits = c(150,-150)) + lims(x = c(-108,-65))
  
pptchangemap
ggsave(pptchangemap, filename = "~/Documents/pptchangemap2.tiff")





#










# ELVI
plot(endo_herb_ELVI$lon, endo_herb_ELVI$lat)
hist(endo_herb_ELVI$year)
hist(endo_herb_ELVI$lon)
  
long_date_mod <- glm(Endo_status_liberal == 1 ~ lon*year, data = subset(endo_herb_ELVI, lon < -90), family = binomial)
summary(long_date_mod)
newdat1920 <- data.frame(lon = seq(-120,-60,1), year = 1920)
newdat1950 <- data.frame(lon = seq(-120,-60,1), year = 1950)
newdat2000 <- data.frame(lon = seq(-120,-60,1), year = 2000)
newdat <- rbind(newdat1920, newdat1950, newdat2000)
y_pred1920 <- predict(long_date_mod, newdata = newdat1920, type = "response")
y_pred1950 <- predict(long_date_mod, newdata = newdat1950, type = "response")
y_pred2000 <- predict(long_date_mod, newdata = newdat2000, type = "response")


y_pred <- predict(long_date_mod, newdata = newdat, type = "response")
newpred <- cbind(newdat, y_pred)

plot(endo_herb_ELVI$lon, endo_herb_ELVI$Endo_status_liberal)
lines(newdat1920$lon, y_pred1920, col = "red3")
lines(newdat1950$lon, y_pred1950, col = "gray39")
lines(newdat1950$lon, y_pred2000, col = "royalblue3")
 

anova(long_date_mod, test = "Chisq")

ggplot() +
  geom_point(data = binned_ELVI,aes(x = mean_lon, y = mean_endo, size = sample, color = binned_year)) + 
  geom_line(data = newpred, aes(x = lon, y = y_pred, group = year, color = as.character(year))) + 
  theme_classic() + 
  scale_colour_manual(values = c("#fc8d59", "#636363", "#91bfdb", "#fc8d69", "#635363", "#81bfdb")) +
  xlim(-100,-85)


# messing around with latitude

lat_date_mod <- glm(Endo_status_liberal == 1 ~ lat*year, data = subset(endo_herb_ELVI, lat < -90), family = binomial)
lat_date_mod <- glm(Endo_status_liberal == 1 ~ lat*year, data = endo_herb_AGHY, family = binomial)

newdat1900 <- data.frame(lat = seq(25,100,1), year = 1900)
newdat1950 <- data.frame(lat = seq(25,100,1), year = 1950)
newdat2000 <- data.frame(lat = seq(25,100,1), year = 2000)

y_pred1900 <- predict(lat_date_mod, newdata = newdat1900, type = "response")
y_pred1950 <- predict(lat_date_mod, newdata = newdat1950, type = "response")
y_pred2000 <- predict(lat_date_mod, newdata = newdat2000, type = "response")

plot(endo_herb_AGHY$lat, endo_herb_AGHY$Endo_status_liberal)
lines(newdat1900$lat, y_pred1900, col = "red3")
lines(newdat1950$lat, y_pred1950, col = "gray39")
lines(newdat1950$lat, y_pred2000, col = "royalblue3")



anova(lat_date_mod, test = "Chisq")



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

# Messing around with Stan models

invlogit<-function(x){exp(x)/(1+exp(x))}
logit = function(x) { log(x/(1-x)) }

#########################################################################################################
# Bernoulli GLM for endo ~  lon + year + lon*year   -------------------------
########################################################################################
## here is the Stan model ##
## run this to optimize computer system settings
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
set.seed(123)

## MCMC settings
ni <- 1000
nb <- 500
nc <- 3
endo_herb_forstan <- endo_herb_georef %>% 
  mutate(SPP_CODE = case_when(grepl("AGHY", Sample_id) ~ 0,
                              grepl("ELVI", Sample_id) ~ 1)) %>% 
  filter(!is.na(year), !is.na(lon), !is.na(SPP_CODE))

endo_herb_list <- list(endo = as.integer(endo_herb_forstan$Endo_status_liberal),
                  lon = endo_herb_forstan$lon,
                  year = endo_herb_forstan$year,
                  spp = endo_herb_forstan$SPP_CODE,
                  N = length(endo_herb_forstan$Endo_status_liberal),
                  K = 3L)

sink("endoherb.stan")
cat("
    data { 
     int<lower=0> N;                       // number of observations
     int<lower=0> K;                       // number of predictors
     int<lower = 0, upper = 1> endo[N];    // infection status (response)
     real<lower = 0> year[N];                // year of collection
     real<upper = 0> lon[N];             // longitude
     int<lower = 0, upper = 1> spp[N];    // Spp specific intercept
    }
    
    parameters {
    vector[K] beta;
    real<lower = 0> sigma;
    }
    
    transformed parameters{
    real mu[N];                             //Linear Predictor
    for(n in 1:N){
    mu[n] = beta[1] + beta[2]*lon[n] + beta[3]*year[n];
    }
    }
    model {
    // Priors
    beta ~ normal(0,sigma); // prior for predictor intercepts
    // Likelihood
        endo ~ bernoulli_logit(mu);
    }
    
         generated quantities{
    }
  
      ", fill = T)
sink()

stanmodel <- stanc("endoherb.stan")

sm<- stan(file = "endoherb.stan", data = endo_herb_list,
              iter = ni, warmup = nb, chains = nc, save_warmup = FALSE)
# saveRDS(sm, file = "endoherb_fit.rds")


