# Title: Digital Herbarium Records and Sample Integration
# Purpose: Imports and merges downloaded digitized herbarium records with the Endo_Herbarium database
# Authors: Joshua Fowler Ella Segal, Lani Dufresne, and Tom Miller
# Updated: Nov. 26, 2020

library(tidyverse) # for data manipulation and ggplot
# library(slider) # add on to tidyverse to calculate sliding windows for climate data
library(fuzzyjoin) # add on to tidyverse to merge tables on nearest values
library(readxl)
library(lubridate)
library(rstan)
library(brms) #Stan package for bayesian models similar to lme4
library(modelr)
library(bbmle) #use for AIC
# library(tidybayes)
library(car)
library(sf) # for making maps
# library(here) # for making maps (lets you set filepaths)
library(ggmap) # for making maps
library(prism)
library(raster) ##working with raster data
# library(reshape2)
# library(viridis)
# library(gganimate)
# library(gifski)
library(INLA)
################################################################################
############ Read in digitized herbarium records ############################### 
################################################################################

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
  mutate(Spp_code = "AGPE")

UTAustin_torch <- rbind(AGHY_UTAustin, ELVI_UTAustin, AGPE_UTAustin)

# Texas A&M digitized records (Includes both AGHY and ELVI)
AM_records <- read_csv(file = "~/Dropbox/Josh&Tom - shared/Endo_Herbarium/DigitizedHerbariumRecords/TexasA&M_digitized_records/Fowler Data.csv") %>% 
  unite(Institution_specimen_id, CollectionCode, id, sep = "") %>%
  separate(DateCollected, into = c("year", "month", "day"), remove = FALSE) %>% 
  mutate(year = as.numeric(year), month = as.numeric(month), day = as.numeric(day)) %>% 
  dplyr::select(-contains("X"))


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

# Lani's year and county info for her AGPE samples
AGPE_meta <- read_csv(file = "~/Dropbox/Josh&Tom - shared/Endo_Herbarium/DigitizedHerbariumRecords/agpe_working.csv") %>% 
  dplyr::select(ID1,ID2) %>% 
  rename(Specimen_id = ID1, year = ID2)

# Reading in the Louisiana state university digitized records (all three species, using TORCH search term: "Elymus virginicus; Elymus virginicus f. australis; Elymus virginicus f. hirsutiglumis; Elymus virginicus f. jejunus; Elymus virginicus f. lasiolepis; Elymus virginicus f. submuticus; Elymus virginicus subsp. interruptus; Elymus virginicus subsp. villosus; Elymus virginicus var. arcuatus; Elymus virginicus var. australis; Elymus virginicus var. glabriflorus; Elymus virginicus var. glaucus; Elymus virginicus var. halophilus; Elymus virginicus var. hirsutiglumis; Elymus virginicus var. intermedius; Elymus virginicus var. jejunus; Elymus virginicus var. jenkensii; Elymus virginicus var. jenkinsii; Elymus virginicus var. micromeris; Elymus virginicus var. submuticus; Elymus virginicus var. minor; Elymus virginicus var. virginicus; Agrostis hyemalis; Agrostis hyemalis f. exaristata; Agrostis hyemalis f. tuckermanii; Agrostis hyemalis var. elata; Agrostis hyemalis var. geminata; Agrostis hyemalis var. hyemalis; Agrostis hyemalis var. keweenawensis; Agrostis hyemalis var. laxiflora; Agrostis hyemalis var. nutkaensis; Agrostis hyemalis var. oreophila; Agrostis hyemalis var. scabra; Agrostis hyemalis var. subrepens; Agrostis hyemalis var. tenuis; Agrostis hiemalis; Agrostis hiemalis var. geminata; Agrostis hiemalis var. laxiflora; Agrostis hiemalis var. subrepens; Agrostis perennans; Agrostis perennans f. aestivalis; Agrostis perennans f. atherophora; Agrostis perennans f. chaetophora; Agrostis perennans f. chaotophora; Agrostis perennans f. perennans; Agrostis perennans var. aestivalis; Agrostis perennans var. elata; Agrostis perennans var. humilis; Agrostis perennans var. perennans")

LSU_records <- read_csv(file = "~/Dropbox/Josh&Tom - shared/Endo_Herbarium/DigitizedHerbariumRecords/LouisianaStateUniversity_digitized_records/SymbOutput_2022-04-15_085400_DwC-A/occurrences.csv",
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
                                         collID = col_double(),
                                         recordID = col_character(),
                                         references = col_character())) %>% 
  mutate(Spp_code = case_when(genus == "Elymus" ~ "ELVI",
                              genus == "Agrostis" & specificEpithet == "hyemalis" | specificEpithet == "hiemalis" ~ "AGHY",
                              genus == "Agrostis" & specificEpithet == "perennans" ~ "AGPE")) %>%  # there are some Agrostis scabra, that may need to be sorted out cause they could be part of hyemalis
  dplyr::select(id, Spp_code, catalogNumber, country, stateProvince, county, municipality, locality, decimalLatitude, decimalLongitude, coordinateUncertaintyInMeters, eventDate, day, month, year) 

# Reading in the Oklahoma state university digitized records
OKLA_records <- read_csv(file = "~/Dropbox/Josh&Tom - shared/Endo_Herbarium/DigitizedHerbariumRecords/OklahomaStateUniversity_digitized_records/SymbOutput_2022-04-15_090328_DwC-A/occurrences.csv",
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
                                          collID = col_double(),
                                          recordID = col_character(),
                                          references = col_character())) %>% 
  mutate(new_id = paste0("OKLA", otherCatalogNumbers)) %>% #this is the id that I recorded in endo_herbarium, but with OKLA on the front. These sheets had a lot that were digitized, but on this number not the barcode
  mutate(Spp_code = case_when(genus == "Elymus" ~ "ELVI",
                              genus == "Agrostis" & specificEpithet == "hyemalis" | specificEpithet == "hiemalis" ~ "AGHY",
                              genus == "Agrostis" & specificEpithet == "perennans" ~ "AGPE")) %>%    # there are some Agrostis scabra, that may need to be sorted out cause they could be part of hyemalis
  dplyr::select(id, catalogNumber, new_id, genus, specificEpithet, Spp_code, catalogNumber, country, stateProvince, county, municipality, locality, decimalLatitude, decimalLongitude, coordinateUncertaintyInMeters, eventDate, day, month, year) 

# Reading in the University of Oklahoma digitized records
OKL_records <- read_csv(file = "~/Dropbox/Josh&Tom - shared/Endo_Herbarium/DigitizedHerbariumRecords/UniversityofOklahoma_digitized_records/SymbOutput_2022-04-15_091522_DwC-A/occurrences.csv",
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
                                         collID = col_double(),
                                         recordID = col_character(),
                                         references = col_character())) %>% 
  mutate(Spp_code = case_when(genus == "Elymus" ~ "ELVI",
                              genus == "Agrostis" & specificEpithet == "hyemalis" | specificEpithet == "hiemalis" ~ "AGHY",
                              genus == "Agrostis" & specificEpithet == "perennans" ~ "AGPE")) %>%  # there are some Agrostis scabra, that may need to be sorted out cause they could be part of hyemalis
  dplyr::select(id, catalogNumber, otherCatalogNumbers, Spp_code, catalogNumber, country, stateProvince, county, municipality, locality, decimalLatitude, decimalLongitude, coordinateUncertaintyInMeters, eventDate, day, month, year) #For the specimens that are digitized, they are found in the otherCatalogNumbers column. For OKL, I recorded a string of id's. We need to take the first, which is usually of the form "OKL######"


# Reading in the University of Kansas digitized records (downloaded from TORCH, but may need to cross check with the KU botany search)
KANU_records <- read_csv(file = "~/Dropbox/Josh&Tom - shared/Endo_Herbarium/DigitizedHerbariumRecords/UniversityofKansas_digitized_records/SymbOutput_2022-04-15_101720_DwC-A/occurrences.csv",
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
                                          collID = col_double(),
                                          recordID = col_character(),
                                          references = col_character())) %>% 
  mutate(new_id = paste0("KANU00", catalogNumber)) %>% #this is the id that I recorded in endo_herbarium, but with KANU00 on the front.
  mutate(Spp_code = case_when(genus == "Elymus" ~ "ELVI",
                              genus == "Agrostis" & specificEpithet == "hyemalis" | specificEpithet == "hiemalis" ~ "AGHY",
                              genus == "Agrostis" & specificEpithet == "perennans" ~ "AGPE")) %>%   # there are some Agrostis scabra, that may need to be sorted out cause they could be part of hyemalis
  dplyr::select(id, new_id, catalogNumber, otherCatalogNumbers, Spp_code, catalogNumber, country, stateProvince, county, municipality, locality, decimalLatitude, decimalLongitude, coordinateUncertaintyInMeters, eventDate, day, month, year) #For the specimens that are digitized, they are found in the otherCatalogNumbers column
  
# Reading in the Missouri Botanic Garden digitized records # Mobot we need to find the specimens based on collector name and number
MOBOT_records <- read_csv(file = "~/Dropbox/Josh&Tom - shared/Endo_Herbarium/DigitizedHerbariumRecords/MOBOT_digitized_records/SymbOutput_2022-04-15_100805_DwC-A/occurrences.csv",
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
                                           collID = col_double(),
                                           recordID = col_character(),
                                           references = col_character())) %>%  
  mutate(coll_lastname = word(recordedBy, -1),
         new_id = paste0(coll_lastname, recordNumber)) %>% 
  mutate(Spp_code = case_when(genus == "Elymus" ~ "ELVI",
                                genus == "Agrostis" & specificEpithet == "hyemalis" | specificEpithet == "hiemalis" ~ "AGHY",
                                genus == "Agrostis" & specificEpithet == "perennans" ~ "AGPE")) %>%  # there are some Agrostis scabra, that may need to be sorted out cause they could be part of hyemalis
  dplyr::select(id, new_id, recordedBy, coll_lastname, catalogNumber, otherCatalogNumbers, Spp_code, catalogNumber, country, stateProvince, county, municipality, locality, decimalLatitude, decimalLongitude, coordinateUncertaintyInMeters, eventDate, day, month, year) #For the specimens that are digitized, they are found in the otherCatalogNumbers column
  
################################################################################
############ Read in endophyte scores and our transcribed specimen data ############################### 
################################################################################

# Read in our transcribed datasheet from google sheets which I have backed up in dropbox
# I am reading these in as csv files that I saved from the excel file because excel does weird things with the date entries.
specimen_info <- read_csv(file = "~joshuacfowler/Dropbox/Josh&Tom - shared/Endo_Herbarium/Endo_Herbarium_specimen.csv") %>% 
  dplyr::select(-contains("...")) %>% 
  mutate(eventDate = Date_collected) %>% 
  mutate(TEXT_Date_collected = gsub("'", "", TEXT_Date_collected)) %>%  #TEXT_Date_collected is saved with a starting ' to preserve the date as text while saving from excel, so I remove that here
  separate(TEXT_Date_collected, into = c("month", "day", "year"), remove = FALSE) %>% 
  mutate(year = as.numeric(year), month = as.numeric(month), day = as.numeric(day)) %>% 
  mutate(new_id = case_when(grepl("LL00", Institution_specimen_id) ~ gsub("[a-zA-Z ]", "", Institution_specimen_id),
                            grepl("MO",Original_herbarium_id) ~ substr(Institution_specimen_id, 2, nchar(Institution_specimen_id)),
                            grepl("OKL_", Specimen_id) ~ str_split(Institution_specimen_id, fixed(";"), simplify = TRUE)[,1], 
                            TRUE ~ Institution_specimen_id)) %>% 
  mutate(Specimen_id_temp = Specimen_id) %>% 
  separate(Specimen_id_temp, c("Herbarium_id", "Spp_code", "Specimen_no"), "_") %>% 
  filter(!is.na(Specimen_id)) %>% 
  left_join(AGPE_meta, by = c("Specimen_id")) %>% 
  mutate(year = case_when(Spp_code != "AGPE" ~ year.x,
                          Spp_code == "AGPE" & is.na(year.y) ~ year.x,
                          Spp_code == "AGPE" & is.na(year.x) ~ year.y)) %>% 
  dplyr::select(-year.x, - year.y)

# This is the sample info and we will filter for only those that we have scored so far.
sample_info <- read_csv(file = "~joshuacfowler/Dropbox/Josh&Tom - shared/Endo_Herbarium/Endo_Herbarium_sample.csv") %>%  
  dplyr::select(-contains("...")) %>% 
  mutate(Specimen_id_temp = Specimen_id) %>% 
  separate(Specimen_id_temp, c("Herbarium_id", "Spp_code", "Specimen_no"), "_")

testsample_info<- read_csv(file = "~joshuacfowler/Dropbox/Josh&Tom - shared/Endo_Herbarium/Endo_Herbarium_sample.csv") %>% 
  group_by(Specimen_id) %>% 
  summarize(count = n()) %>% 
  filter(count == 2)

endo_scores <- specimen_info %>% 
  merge(sample_info, by = c("Specimen_id" = "Specimen_id", "Herbarium_id" = "Herbarium_id", "Spp_code" = "Spp_code", "Specimen_no" = "Specimen_no")) %>% 
  mutate(Specimen_no = as.numeric(Specimen_no)) %>% 
  dplyr::select(-contains("X1"),-contains("X2"))


# Read in the fitness and size data that we have collected

aghy_fitness <- read_csv(file = "~joshuacfowler/Dropbox/Josh&Tom - shared/Endo_Herbarium/AGHY_fitness_data.csv") %>% 
  filter(!is.na(surface_area_cm2)) %>% 
  dplyr::select(-seed, -stem)
  
elvi_fitness <- read_csv(file = "~joshuacfowler/Dropbox/Josh&Tom - shared/Endo_Herbarium/ELVI_fitness_data.csv") %>% 
  filter(!is.na(Area)) %>%
  pivot_longer(cols = c(Length1, Length2, Length3, Length4, Length5,
                        Length_w_fuzzy_1, Length_w_fuzzy_2, Length_w_fuzzy_3, Length_w_fuzzy_4, Length_w_fuzzy_5),
               names_to = c("Infl_measurement", "Infl_no"), names_pattern = "(.+)(.+)") %>% 
  filter(!is.na(value)) %>% 
  pivot_wider(names_from = Infl_measurement, values_from = value) %>% 
  group_by(Specimen_id, Institution_specimen_id, new_id, Herbarium_id, Spp_code, Specimen_no, Area) %>% 
  summarize(mean_infl_length = mean(Length),
            mean_inflplusawn_length = mean(Length_w_fuzzy_),
            infl_count = n())

fitness_data <- full_join(aghy_fitness, elvi_fitness, by = c("Specimen_id" = "Specimen_id",
                                                         "Institution_specimen_id" = "Institution_specimen_id", 
                                                         "new_id"= "new_id", 
                                                         "Herbarium_id" = "Herbarium_id", 
                                                         "Spp_code" = "Spp_code", 
                                                         "Specimen_no" = "Specimen_no", 
                                                         "surface_area_cm2" = "Area"))  %>% 
  dplyr::select( -Institution_specimen_id, -new_id)


endo_herb <- endo_scores %>% 
  left_join(fitness_data, by = c("Specimen_id" = "Specimen_id",
                                 "Herbarium_id" = "Herbarium_id", 
                                 "Spp_code" = "Spp_code", 
                                 "Specimen_no" = "Specimen_no"))
  

# These are the matches from BRIT for transcription using the file that Jason Best shared of their database, Mar 25th, 2020
BRIT_matches <- read_csv(file = "~/Dropbox/Josh&Tom - shared/Endo_Herbarium/DigitizedHerbariumRecords/BRIT_records/id matches for Torch transcription/fowler_matches_Mar25.csv") %>% 
  dplyr::select("catalogNumber") %>% 
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
  dplyr::select(Sample_id, Institution_specimen_id, Spp_code, new_id, Country, State, County, Municipality, Locality, year, month, day, tissue_type, seed_scored_1, seed_eplus_1, Endo_status_liberal_1, Endo_status_conservative_1, Date_scored_1, scorer_id_1, seed_scored_2, seed_eplus_2, Endo_status_liberal_2, Endo_status_conservative_2, Date_scored_2, scorer_id_2, surface_area_cm2, mean_infl_length, mean_inflplusawn_length, infl_count)

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
  dplyr::select(Sample_id, Institution_specimen_id, Spp_code, new_id, Country, State, County, Municipality, Locality, year, month, day, tissue_type, seed_scored_1, seed_eplus_1, Endo_status_liberal_1, Endo_status_conservative_1, Date_scored_1, scorer_id_1, seed_scored_2, seed_eplus_2, Endo_status_liberal_2, Endo_status_conservative_2, Date_scored_2, scorer_id_2, surface_area_cm2, mean_infl_length, mean_inflplusawn_length, infl_count)

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
  dplyr::select(Sample_id, Institution_specimen_id, Spp_code, new_id, Country, State, County, Municipality, Locality, year, month, day, tissue_type, seed_scored_1, seed_eplus_1, Endo_status_liberal_1, Endo_status_conservative_1, Date_scored_1, scorer_id_1, seed_scored_2, seed_eplus_2, Endo_status_liberal_2, Endo_status_conservative_2, Date_scored_2, scorer_id_2, surface_area_cm2, mean_infl_length, mean_inflplusawn_length, infl_count) %>% 
filter(!duplicated(Sample_id))


# Merge in the LSU records that we have so far
endo_herb4 <- endo_herb3 %>% 
  left_join(LSU_records, by = c("new_id" = "catalogNumber")) %>% 
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
  dplyr::select(Sample_id, Institution_specimen_id, Spp_code, new_id, Country, State, County, Municipality, Locality, year, month, day, tissue_type, seed_scored_1, seed_eplus_1, Endo_status_liberal_1, Endo_status_conservative_1, Date_scored_1, scorer_id_1, seed_scored_2, seed_eplus_2, Endo_status_liberal_2, Endo_status_conservative_2, Date_scored_2, scorer_id_2, surface_area_cm2, mean_infl_length, mean_inflplusawn_length, infl_count) %>% 
  filter(!duplicated(Sample_id))

# Merge in the OKL records that we have so far. This one is weird, and need to figure out which column to match on.
endo_herb5 <- endo_herb4 %>% 
  left_join(OKL_records, by = c("new_id" = "otherCatalogNumbers")) %>% 
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
  dplyr::select(Sample_id, Institution_specimen_id, Spp_code, new_id, Country, State, County, Municipality, Locality, year, month, day, tissue_type, seed_scored_1, seed_eplus_1, Endo_status_liberal_1, Endo_status_conservative_1, Date_scored_1, scorer_id_1, seed_scored_2, seed_eplus_2, Endo_status_liberal_2, Endo_status_conservative_2, Date_scored_2, scorer_id_2, surface_area_cm2, mean_infl_length, mean_inflplusawn_length, infl_count) %>%
  filter(!duplicated(Sample_id))

# Merge in the OKLA records that we have so far
endo_herb6 <- endo_herb5 %>% 
  left_join(OKLA_records, by = c("new_id" = "new_id")) %>% 
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
  dplyr::select(Sample_id, Institution_specimen_id, Spp_code, new_id, Country, State, County, Municipality, Locality, year, month, day, tissue_type, seed_scored_1, seed_eplus_1, Endo_status_liberal_1, Endo_status_conservative_1, Date_scored_1, scorer_id_1, seed_scored_2, seed_eplus_2, Endo_status_liberal_2, Endo_status_conservative_2, surface_area_cm2, mean_infl_length, mean_inflplusawn_length, infl_count) %>% 
  filter(!duplicated(Sample_id))


# Merge in the KANU records that we have so far
endo_herb7 <- endo_herb6 %>% 
  left_join(KANU_records, by = c("new_id" = "new_id")) %>% 
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
  dplyr::select(Sample_id, Institution_specimen_id, Spp_code, new_id, Country, State, County, Municipality, Locality, year, month, day, tissue_type, seed_scored_1, seed_eplus_1, Endo_status_liberal_1, Endo_status_conservative_1, Date_scored_1, scorer_id_1, seed_scored_2, seed_eplus_2, Endo_status_liberal_2, Endo_status_conservative_2, surface_area_cm2, mean_infl_length, mean_inflplusawn_length, infl_count) %>% 
  filter(!duplicated(Sample_id))



# Merge in the MOBOT records that we have so far
ci_str_detect <- function(x, y){str_detect(x, regex(y, ignore_case = TRUE))}
endo_herb8 <- endo_herb7 %>% 
  fuzzy_left_join(MOBOT_records, match_fun = ci_str_detect, by = c("new_id" = "new_id")) %>%  
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
  mutate(new_id = new_id.x) %>% 
  dplyr::select(Sample_id, Institution_specimen_id, Spp_code, new_id, Country, State, County, Municipality, Locality, year, month, day, tissue_type, seed_scored, seed_eplus, Endo_status_liberal, Endo_status_conservative, surface_area_cm2, mean_infl_length, mean_inflplusawn_length, infl_count) %>% 
  filter(!duplicated(Sample_id))




specimen_counts <- endo_herb8 %>% 
  mutate(Spp_code = case_when(grepl("AGHY", Sample_id) ~ "AGHY",
                              grepl("ELVI", Sample_id) ~ "ELVI",
                              grepl("AGPE", Sample_id) ~ "AGPE")) %>% 
  group_by(Spp_code, tissue_type) %>% 
  summarize(n())
# Now I am going to link these county/locality records to a gps point with ggmap
# This requires and API key which you can set up through google, look at ?register_google.
# There are restrictions to the total number of queries that you can do per day and month, and if you go over, it costs money, so we will save the output. I believe we have a free trial for year.
# One other note, is that we are only using the level of county/city/state which I believe should be pretty accurate through google. I'm not sure that it could accurately do a more detailed locality string
# I have found software that could do this, but would require a human to double check the output.
#
# register_google()
# endo_herb_georef <-endo_herb8 %>%
#   unite("location_string" , sep = ", " , Municipality,County,State,Country, remove = FALSE, na.rm = TRUE) %>%
#   # filter(Endo_status_liberal <= 1) %>%
#   mutate_geocode(location_string) # Uncomment this to run the geocoding.
# write_csv(endo_herb_georef, file = "~/Dropbox/Josh&Tom - shared/Endo_Herbarium/DigitizedHerbariumRecords/endo_herb_georef.csv")
endo_herb_georef <- read_csv(file = "~/Dropbox/Josh&Tom - shared/Endo_Herbarium/DigitizedHerbariumRecords/endo_herb_georef.csv") %>%
  filter(Country != "Canada") %>%
  mutate(Spp_code = case_when(grepl("AGHY", Sample_id) ~ "AGHY",
                              grepl("ELVI", Sample_id) ~ "ELVI",
                              grepl("AGPE", Sample_id) ~ "AGPE"))

################################################################################
############ Visualizing samples and analysis ############################### 
################################################################################


specimen_counts <- endo_herb_georef %>% 
  mutate(Spp_code = case_when(grepl("AGHY", Sample_id) ~ "AGHY",
                              grepl("ELVI", Sample_id) ~ "ELVI",
                              grepl("AGPE", Sample_id) ~ "AGPE")) %>% 
  group_by(Spp_code, tissue_type) %>% 
  summarize(n())
# Now we can explore the data
hist(endo_herb_georef$year)
plot(endo_herb_georef$lon, endo_herb_georef$lat)
plot(endo_herb_georef$lon, endo_herb_georef$Endo_status_liberal)
plot(endo_herb_georef$lat, endo_herb_georef$Endo_status_liberal)
plot(endo_herb_georef$year, endo_herb_georef$Endo_status_liberal)

# counts of scores
table(endo_herb_georef$Spp_code, endo_herb_georef$Endo_status_liberal) # Spp code has NA's, have to look in sample id
table(endo_herb_georef$Spp_code) # Spp code has NA's, have to look in sample id

endo_herb_AGHY <- endo_herb_georef %>% 
  filter(Spp_code == "AGHY") %>% 
  filter(!is.na(lon) & !is.na(year)) 


# messing aruond with INLA model fitting
# generate mesh
coords <- cbind(endo_herb_AGHY$lon, endo_herb_AGHY$lat)
max.edge = diff(range(coords[,2]))/(3*5)
mesh <- inla.mesh.2d(loc = coords, max.edge = c(1,2)*max.edge,
            cutoff = max.edge/5)
mesh$n # the number of mesh vertices
plot(mesh)
points(coords, col = "red")


formula <- formula(Endo_status_liberal ~ year + f(lon, model = "iid"))
I <- inla(Endo_status_liberal ~ year + f(lon, model = "iid") , family = "binomial",
          data = endo_herb_AGHY, verbose = TRUE)

 binned_AGHY <- endo_herb_AGHY %>% 
   mutate(binned_lon = cut(lon, breaks = 20), binned_year = cut(year, breaks = 2)) %>%  
   group_by(binned_lon, binned_year) %>%   
   summarise(mean_lon = mean(lon),
             mean_year = mean(year),
             mean_endo = mean(Endo_status_liberal),
             sample = n())
 binned3_AGHY <- endo_herb_AGHY %>% 
   mutate(binned_lon = cut(lon, breaks = 10), binned_year = cut(year, breaks = 3)) %>%  
   group_by(binned_lon, binned_year) %>%   
   summarise(mean_lon = mean(lon),
             mean_year = mean(year),
             mean_endo = mean(Endo_status_liberal),
             sample = n())
 

endo_herb_ELVI <- endo_herb_georef %>% 
  filter(Spp_code == "ELVI") %>% 
  filter(!is.na(lon) & !is.na(year))

binned_ELVI <- endo_herb_ELVI %>% 
  mutate(binned_lon = cut(lon, breaks =10), binned_year = cut(year, breaks = 2)) %>%  
  group_by(binned_lon, binned_year) %>%   
  summarise(mean_lon = mean(lon),
            mean_year = mean(year),
            mean_endo = mean(Endo_status_liberal),
            sample = n())

binned_ELVI_lat <- endo_herb_ELVI %>% 
  mutate(binned_lat = cut(lat, breaks =10), binned_year = cut(year, breaks = 2)) %>%  
  group_by(binned_lat, binned_year) %>%   
  summarise(mean_lat = mean(lat),
            mean_year = mean(year),
            mean_endo = mean(Endo_status_liberal),
            sample = n())


# AGPE need to get year records from Lani/ many of these are only partially digitized
endo_herb_AGPE <- endo_herb_georef %>%
  filter(Spp_code == "AGPE") %>% 
  filter(!is.na(lon) & !is.na(year))

binned_AGPE <- endo_herb_AGPE %>%
  mutate(binned_lon = cut(lon, breaks = 6), binned_year = cut(year, breaks = 2)) %>%
  group_by(binned_lon, binned_year) %>%
  summarise(mean_lon = mean(lon),
            mean_year = mean(year),
            mean_endo = mean(Endo_status_liberal),
            sample = n())

############################################################################
##### Spatial and temporal trends ##########################################
############################################################################
# AGHY
plot(endo_herb_AGHY$lon, endo_herb_AGHY$lat)
plot(endo_herb_AGHY$year,endo_herb_AGHY$Endo_status_liberal)
plot(endo_herb_AGHY$lon,endo_herb_AGHY$Endo_status_liberal)
hist(endo_herb_AGHY$year)
hist(endo_herb_AGHY$lon)

heatmap_aghy <- endo_herb_AGHY %>%
  mutate(binned_lon = cut(lon, breaks = 25), binned_lat = cut(lat, breaks = 25), binned_year = cut(year, breaks = 2)) %>%  
  group_by(binned_lon, binned_lat, binned_year) %>%   
  summarise(mean_lon = mean(lon),
            mean_lat = mean(lat),
            mean_year = mean(year),
            mean_endo = mean(Endo_status_liberal),
            sample = n())
head(heatmap_aghy)

ggplot(heatmap_aghy, aes(x = binned_lon, y = binned_year, fill = mean_endo)) +
  geom_raster()
ggplot(heatmap_aghy, aes(x = binned_lon, y = binned_lat, fill = mean_endo, alpha = sample)) +
  facet_wrap(~binned_year)+
  geom_raster()

ggplot(heatmap_aghy, aes(x = binned_lon, y = binned_lat, z = mean_endo)) +
  facet_wrap(~binned_year)+
  geom_contour()
# model selection
aghy_models <- list()
aghy_models[[1]] <- glm(Endo_status_liberal ~ lon, data = endo_herb_AGHY, family = binomial)
aghy_models[[2]] <- glm(Endo_status_liberal ~ lat, data = endo_herb_AGHY, family = binomial)
aghy_models[[3]] <- glm(Endo_status_liberal ~ year, data = endo_herb_AGHY, family = binomial)
aghy_models[[4]] <- glm(Endo_status_liberal ~ lon + lat, data = endo_herb_AGHY, family = binomial)
aghy_models[[5]] <- glm(Endo_status_liberal ~ lon * lat, data = endo_herb_AGHY, family = binomial)
aghy_models[[6]] <- glm(Endo_status_liberal ~ lon + year, data = endo_herb_AGHY, family = binomial)
aghy_models[[7]] <- glm(Endo_status_liberal ~ lon * year, data = endo_herb_AGHY, family = binomial)
aghy_models[[8]] <- glm(Endo_status_liberal ~ lat + year, data = endo_herb_AGHY, family = binomial)
aghy_models[[9]] <- glm(Endo_status_liberal ~ lon * year, data = endo_herb_AGHY, family = binomial)
aghy_models[[10]] <- glm(Endo_status_liberal ~ lat + lon + year, data = endo_herb_AGHY, family = binomial)
aghy_models[[11]] <- glm(Endo_status_liberal ~ lat * lon + year, data = endo_herb_AGHY, family = binomial)
aghy_models[[12]] <- glm(Endo_status_liberal ~ lat + lon * year, data = endo_herb_AGHY, family = binomial)
aghy_models[[13]] <- glm(Endo_status_liberal ~ lon + lat * year, data = endo_herb_AGHY, family = binomial)
aghy_models[[14]] <- glm(Endo_status_liberal ~ lat * lon * year, data = endo_herb_AGHY, family = binomial)
aghy_models[[15]] <- glm(Endo_status_liberal ~ 1, data = endo_herb_AGHY, family = binomial)

AICtab(aghy_models) # The best model according to AIC is model 5, longitude by latitude interaction



long_date_mod <- glm(Endo_status_liberal ~ lon+year  , data = subset(endo_herb_AGHY), family = binomial)
anova(long_date_mod, test = "Chisq")
summary(long_date_mod)

Anova(long_date_mod, test = "LR")



newdat1930 <- data.frame(lon = seq(-120,-60,1), year = 1930)
newdat2000 <- data.frame(lon = seq(-120,-60,1), year = 1990)
newdat <- rbind(newdat1930, newdat2000)
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
  geom_point(data = endo_herb_AGHY, aes(x = lon, y = Endo_status_liberal), col = "darkgray", pch = "|")+
  geom_point(data = binned_AGHY,aes(x = mean_lon, y = mean_endo, size = sample, color = binned_year)) +
  geom_line(data = newpred, aes(x = lon, y = pred, group = year, color = as.character(year))) +
  geom_ribbon(data = newpred, aes(x = lon, ymin = lwr, ymax = upr,  group = year, fill = as.factor(year)), alpha = .2)+ 
  theme_classic() + labs(y = "Mean Endophyte Prevalence", x = "Longitude", color = "year") +
  scale_colour_manual(breaks = c( "1930", "1990"),
                      values = c( "#636363","#41b056",    "#636363", "#41b056")) +
  scale_fill_manual(breaks = c(  "1930",  "1990"),
                    values = c(   "#636363", "#41b056")) +
  xlim(-107,-68) + guides(fill = FALSE)
AGHY_herb
ggsave(AGHY_herb, filename = "~/Documents/AGHYherbgreen.tiff", width = 4, height = 3)





## here is the brms Stan model version ##
## run this to optimize computer system settings
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
set.seed(123)

## MCMC settings
ni <- 2000
nb <- 1000
nc <- 3

long_date_mod_bayesian <- brm(formula = Endo_status_liberal ~  lon + year + lon*year, 
                              data = endo_herb_AGHY,
                              family = bernoulli(link = "logit"),
                              prior = c(set_prior("normal(0,5)", class = "b")),
                              # prior = c(set_prior("normal(0,5)", class = "b", coef = "lon"),
                              #           set_prior("normal(0,5)", class = "b", coef = "year"),
                              #           set_prior("normal(0,1)", class = "b", coef = "lon:year")),
                              warmup = nb, 
                              iter = ni, 
                              chains = nc)
# plot the traceplots
mcmc_plot(long_date_mod_bayesian, 
         type = "trace")
# summary of the parameters
summary(long_date_mod_bayesian)

# figure out which priors the brm function is using
prior_summary(long_date_mod_bayesian)

plot(conditional_effects(long_date_mod_bayesian, effects = "lon:year")) +

endo_herb_AGHY %>% 
  mutate(binned_year = cut(year, breaks = 3)) %>%  
  group_by(year = binned_year) %>% 
  data_grid(lon = seq_range(lon, n = 100))+
  add_fitted_draws(long_date_mod_bayesian) %>% 
  ggplot(aes(x = lon, y = Endo_status_liberal)) +
  stat_lineribbon(aes(y = .value)) 
  geom_point(data = binned_AGHY, aes(x = binned_lon, y = binned_year)) +
  scale_fill_brewer(palette = "Greys") +
  scale_color_brewer(palette = "Set2")


AGHY_herb_bayes <- ggplot() +
  geom_point(data = binned_AGHY,aes(x = mean_lon, y = mean_endo, size = sample, color = binned_year)) + 
  geom_line(data = newpred, aes(x = lon, y = pred, group = year, color = as.character(year))) +
  geom_ribbon(data = newpred, aes(x = lon, ymin = lwr, ymax = upr, group = year, fill = as.factor(year)), alpha = .2) +
  theme_classic() + labs(y = "Mean Endophyte Prevalence", x = "Longitude", color = "year")+
  scale_colour_manual(breaks = c("1920", "1950", "2000"),
                      values = c("#fc8d59", "#636363", "#91bfdb", "#fc8d69", "#635363", "#81bfdb")) +
  scale_fill_manual(breaks = c("1920", "1950", "2000"),
                    values = c("#fc8d59", "#636363", "#91bfdb")) +
  xlim(-105,-60) + guides(fill = FALSE)

AGHY_herb_bayes

ggsave(AGHY_herb, filename = "~/Documents/AGHYherb.tiff", width = 4, height = 3)


# ELVI
plot(endo_herb_ELVI$lon, endo_herb_ELVI$lat)
plot(endo_herb_ELVI$year,endo_herb_ELVI$Endo_status_liberal)
plot(endo_herb_ELVI$lon,endo_herb_ELVI$Endo_status_liberal)
hist(endo_herb_ELVI$year)
hist(endo_herb_ELVI$lon)

# model selection
elvi_models <- list()
elvi_models[[1]] <- glm(Endo_status_liberal ~ lon, data = endo_herb_ELVI, family = binomial)
elvi_models[[2]] <- glm(Endo_status_liberal ~ lat, data = endo_herb_ELVI, family = binomial)
elvi_models[[3]] <- glm(Endo_status_liberal ~ year, data = endo_herb_ELVI, family = binomial)
elvi_models[[4]] <- glm(Endo_status_liberal ~ lon + lat, data = endo_herb_ELVI, family = binomial)
elvi_models[[5]] <- glm(Endo_status_liberal ~ lon * lat, data = endo_herb_ELVI, family = binomial)
elvi_models[[6]] <- glm(Endo_status_liberal ~ lon + year, data = endo_herb_ELVI, family = binomial)
elvi_models[[7]] <- glm(Endo_status_liberal ~ lon * year, data = endo_herb_ELVI, family = binomial)
elvi_models[[8]] <- glm(Endo_status_liberal ~ lat + year, data = endo_herb_ELVI, family = binomial)
elvi_models[[9]] <- glm(Endo_status_liberal ~ lon * year, data = endo_herb_ELVI, family = binomial)
elvi_models[[10]] <- glm(Endo_status_liberal ~ lat + lon + year, data = endo_herb_ELVI, family = binomial)
elvi_models[[11]] <- glm(Endo_status_liberal ~ lat * lon + year, data = endo_herb_ELVI, family = binomial)
elvi_models[[12]] <- glm(Endo_status_liberal ~ lat + lon * year, data = endo_herb_ELVI, family = binomial)
elvi_models[[13]] <- glm(Endo_status_liberal ~ lon + lat * year, data = endo_herb_ELVI, family = binomial)
elvi_models[[14]] <- glm(Endo_status_liberal ~ lat * lon * year, data = endo_herb_ELVI, family = binomial)
elvi_models[[15]] <- glm(Endo_status_liberal ~ 1, data = endo_herb_ELVI, family = binomial)

AICtab(elvi_models) # The best model according to AIC is model 2, latitude only




long_date_mod <- glm(Endo_status_liberal ~ lon*year, data = subset(endo_herb_ELVI), family = "binomial")
anova(long_date_mod, test = "Chisq")
summary(long_date_mod)

Anova(long_date_mod, test = "LR")
# plot(long_date_mod)


newdat1930 <- data.frame(lon = seq(-110,-60,1), year = 1930)
newdat2000 <- data.frame(lon = seq(-110,-60,1), year = 2000)
newdat <- rbind(newdat1930, newdat2000)
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



ELVI_herb <- ggplot() +
  geom_point(data = endo_herb_ELVI, aes(x = lon, y = Endo_status_liberal), col = "darkgray", pch = "|")+
  geom_point(data = binned_ELVI,aes(x = mean_lon, y = mean_endo, size = sample, color = binned_year)) +
  geom_line(data = newpred, aes(x = lon, y = pred, group = year, color = as.character(year))) +
  geom_ribbon(data = newpred, aes(x = lon, ymin = lwr, ymax = upr, group = year, fill = as.factor(year)), alpha = .2) +
  theme_classic() + labs(y = "Mean Endophyte Prevalence", x = "Latitude", color = "year")+
  scale_colour_manual(breaks = c("1930", "2000"),
                      values = c(  "#636363", "#4c90dd","#636363", "#4c90dd")) +
  scale_fill_manual(breaks = c("1930", "2000"),
                    values = c( "#636363", "#4c90dd")) +
  xlim(-110,-65) + guides(fill = FALSE)

ELVI_herb
 ggsave(ELVI_herb, filename = "~/Documents/ELVIherb_blue.tiff", width = 4, height = 3)

 
 
#ELVI latitude model
 
 lat_date_mod <- glm(Endo_status_liberal ~ lat*year, data = subset(endo_herb_ELVI), family = "binomial")
 anova(lat_date_mod, test = "Chisq")
 summary(lat_date_mod)
 
 Anova(lat_date_mod, test = "LR")
 # plot(lat_date_mod)
 
 
 newdat1930 <- data.frame(lat = seq(26,48,1), year = 1930)
 newdat2000 <- data.frame(lat = seq(26,48,1), year = 2000)
 newdat <- rbind(newdat1930, newdat2000)
 y_pred <- predict(lat_date_mod, newdata = newdat, type = "response")
 y_CI <- predict(lat_date_mod, newdata = newdat, interval = "confidence", type = "link", se.fit=TRUE)
 linkinv <- family(lat_date_mod)$linkinv ## inverse-link function
 
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
 
 
 
 ELVI_herblat <- ggplot() +
   geom_point(data = endo_herb_ELVI, aes(x = lat, y = Endo_status_liberal), col = "darkgray", pch = "|")+
   geom_point(data = binned_ELVI_lat,aes(x = mean_lat, y = mean_endo, size = sample, color = binned_year)) +
   geom_line(data = newpred, aes(x = lat, y = pred, group = year, color = as.character(year))) +
   geom_ribbon(data = newpred, aes(x = lat, ymin = lwr, ymax = upr, group = year, fill = as.factor(year)), alpha = .2) +
   theme_classic() + labs(y = "Mean Endophyte Prevalence", x = "Latitude", color = "year")+
   scale_colour_manual(breaks = c("1930", "2000"),
                       values = c(  "#636363", "#4c90dd","#636363", "#4c90dd")) +
   scale_fill_manual(breaks = c("1930", "2000"),
                     values = c( "#636363", "#4c90dd")) +
   xlim(26,48) + guides(fill = FALSE)
 
 ELVI_herblat
 ggsave(ELVI_herblat, filename = "~/Documents/ELVIherblat_blue.tiff", width = 4, height = 3)
 
 

# AGPE
plot(endo_herb_AGPE$lon, endo_herb_AGPE$lat)
plot(endo_herb_AGPE$year,endo_herb_AGPE$Endo_status_liberal)
plot(endo_herb_AGPE$lon,endo_herb_AGPE$Endo_status_liberal)
hist(endo_herb_AGPE$year)
hist(endo_herb_AGPE$lon)

# model selection
agpe_models <- list()
agpe_models[[1]] <- glm(Endo_status_liberal ~ lon, data = endo_herb_AGPE, family = binomial)
agpe_models[[2]] <- glm(Endo_status_liberal ~ lat, data = endo_herb_AGPE, family = binomial)
agpe_models[[3]] <- glm(Endo_status_liberal ~ year, data = endo_herb_AGPE, family = binomial)
agpe_models[[4]] <- glm(Endo_status_liberal ~ lon + lat, data = endo_herb_AGPE, family = binomial)
agpe_models[[5]] <- glm(Endo_status_liberal ~ lon * lat, data = endo_herb_AGPE, family = binomial)
agpe_models[[6]] <- glm(Endo_status_liberal ~ lon + year, data = endo_herb_AGPE, family = binomial)
agpe_models[[7]] <- glm(Endo_status_liberal ~ lon * year, data = endo_herb_AGPE, family = binomial)
agpe_models[[8]] <- glm(Endo_status_liberal ~ lat + year, data = endo_herb_AGPE, family = binomial)
agpe_models[[9]] <- glm(Endo_status_liberal ~ lon * year, data = endo_herb_AGPE, family = binomial)
agpe_models[[10]] <- glm(Endo_status_liberal ~ lat + lon + year, data = endo_herb_AGPE, family = binomial)
agpe_models[[11]] <- glm(Endo_status_liberal ~ lat * lon + year, data = endo_herb_AGPE, family = binomial)
agpe_models[[12]] <- glm(Endo_status_liberal ~ lat + lon * year, data = endo_herb_AGPE, family = binomial)
agpe_models[[13]] <- glm(Endo_status_liberal ~ lon + lat * year, data = endo_herb_AGPE, family = binomial)
agpe_models[[14]] <- glm(Endo_status_liberal ~ lat * lon * year, data = endo_herb_AGPE, family = binomial)
agpe_models[[15]] <- glm(Endo_status_liberal ~ 1, data = endo_herb_AGPE, family = binomial)

AICtab(agpe_models) # The best model according to AIC is model 1, longitude only



long_date_mod <- glm(Endo_status_liberal ~ lon* year, data = endo_herb_AGPE, family = binomial)
anova(long_date_mod, test = "Chisq")
summary(long_date_mod)

Anova(long_date_mod, test = "LR")



newdat1930 <- data.frame(lon = seq(-100,-70,1), year = 1930)
newdat2000 <- data.frame(lon = seq(-100,-70,1), year = 2000)
newdat <- rbind(newdat1930, newdat2000)
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



AGPE_herb <- ggplot() +
  geom_point(data = endo_herb_AGPE, aes(x = lon, y = Endo_status_liberal), col = "darkgray", pch = "|")+
  geom_point(data = binned_AGPE,aes(x = mean_lon, y = mean_endo, size = sample, color = binned_year)) + 
  geom_line(data = newpred, aes(x = lon, y = pred, group = year, color = as.character(year))) +
  geom_ribbon(data = newpred, aes(x = lon, ymin = lwr, ymax = upr, group = year, fill = as.factor(year)), alpha = .2) +
  theme_classic() + labs(y = "Mean Endophyte Prevalence", x = "Longitude", color = "year")+
  scale_colour_manual(breaks = c("1930", "2000"),
                      values = c(  "#636363","#eda44e","#636363", "#eda44e")) +
  scale_fill_manual(breaks = c("1930", "2000"),
                    values = c( "#636363", "#eda44e")) +
  xlim(-100,-70) + guides(fill = FALSE)

AGPE_herb
ggsave(AGPE_herb, filename = "~/Documents/AGPEherb_orange.tiff", width = 4, height = 3)




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
  geom_point(data = subset(endo_herb_georef, Spp_code == "AGHY"), aes(x = lon, y = Endo_status_liberal),color = "darkgrey", pch = "|") +
  geom_point(data = longbin_AGHY,aes(x = mean_lon, y = mean_endo, size = sample)) +
  geom_line(data = newpred, aes(x = lon, y = pred)) +
  geom_ribbon(data = newpred, aes(x = lon, ymin = lwr, ymax = upr, alpha = .5)) +
  theme_classic() + labs(y = "Mean Endophyte Prevalence", x = "Longitude")+
  xlim(-110,-70) + guides(alpha = FALSE)

AGHY_herb_long
ggsave(AGHY_herb_long, filename = "~/Documents/AGHYherb_long.tiff", width = 4, height = 3)



### Plot for binned over longitude, all times merged

longbin_ELVI <- endo_herb_ELVI %>%
  mutate(binned_lon = cut(lon, breaks = 25)) %>%  
  group_by(binned_lon) %>%   
  summarise(mean_lon = mean(lon),
            mean_endo = mean(Endo_status_liberal),
            sample = n())

long_mod <- glm(Endo_status_liberal ~ lon, data = (subset(endo_herb_ELVI, lat < 45)), family = binomial)
summary(long_mod)
anova(long_mod, test = "Chisq")

newdatlong <- data.frame(lon = seq(-100,-65,1))
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



ELVI_herb_lon <- ggplot() +
  geom_point(data = subset(endo_herb_georef, Spp_code == "ELVI" ), aes(x = lon, y = Endo_status_liberal), color = "darkgrey", shape = "|")+
  geom_point(data = longbin_ELVI,aes(x = mean_lon, y = mean_endo, size = sample)) + 
  geom_line(data = newpred, aes(x = lon, y = pred)) +
  geom_ribbon(data = newpred, aes(x = lon, ymin = lwr, ymax = upr, alpha = .5)) +
  theme_classic() + labs(y = "Mean Endophyte Prevalence", x ="Longitude")+
  xlim(-100,-70) + guides(alpha = FALSE)

ELVI_herb_lon
ggsave(ELVI_herb_lon, filename = "~/Documents/ELVIherb_lon.tiff", width = 4, height = 3)




### Plot for binned over longitude, all times merged

longbin_AGPE <- endo_herb_AGPE %>% 
  mutate(binned_lon = cut(lon, breaks = 25)) %>%  
  group_by(binned_lon) %>%   
  summarise(mean_lon = mean(lon),
            mean_endo = mean(Endo_status_liberal),
            sample = n())

long_mod <- glm(Endo_status_liberal ~ lon, data = subset(endo_herb_AGPE), family = binomial)
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



AGPE_herb_long <- ggplot() +
  geom_point(data = longbin_AGPE,aes(x = mean_lon, y = mean_endo, size = sample)) + 
  geom_line(data = newpred, aes(x = lon, y = pred)) +
  geom_ribbon(data = newpred, aes(x = lon, ymin = lwr, ymax = upr, alpha = .5)) +
  theme_classic() + labs(y = "Mean Endophyte Prevalence", x = "Longitude")+
  xlim(-100,-70) + guides(alpha = FALSE)

AGPE_herb_long
ggsave(AGPE_herb_long, filename = "~/Documents/AGPEherb_long.tiff", width = 4, height = 3)



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


############################################################################
##### Pulling in Climate Data ###############################################
############################################################################

# Making a map of scored AGHY
usa <- st_as_sf(maps::map("state", fill=TRUE, plot =FALSE))
texas <- usa %>% 
  filter(ID == "texas")
AGHY_herb_map <- ggplot()+
  geom_sf(data = usa, fill = "white") +
  geom_point(data = endo_herb_AGHY, aes(x = lon, y = lat), lwd = 2,alpha = .5) +
  theme_minimal() + lims(x = c(-105,-72)) + labs(x = c("Longitude"), y = c("Latitude"))

AGHY_herb_map
ggsave(AGHY_herb_map, filename = "~/Documents/AGHY_herb_map.tiff")

AGHY_endocolor_map <- ggplot()+
  geom_sf(data = usa, fill = "white") +
  geom_point(data = endo_herb_AGHY, aes(x = lon, y = lat, color = as.factor(Endo_status_liberal)), lwd = 2,alpha = .5) +
  theme_minimal() + scale_color_manual(values = c("#636363", "#2c7fb8")) + lims(x = c(-105,-72)) + labs(x = c("Longitude"), y = c("Latitude"), color = "Endophyte presence")

AGHY_endocolor_map
ggsave(AGHY_endocolor_map, filename = "~/Documents/AGHY_endocolor_map.tiff")

endo_herb_georef_seed <- endo_herb_georef %>% 
  filter(tissue_type == "seed") %>% 
  filter(!is.na(Spp_code)) %>% 
  mutate(year_bin = case_when(year < 1900 ~ "pre-1900",
                              1900<=year & year<1930 ~ "1900-1930",
                              1930<=year & year<1960 ~ "1930-1960",
                              1960<=year & year<1990 ~ "1960-1990",
                              1990<=year ~ "post-1990")) %>% 
  mutate(year_bin = fct_relevel(year_bin, "pre-1900", "1900-1930", "1930-1960","1960-1990","post-1990")) %>% 
  mutate(Sample_id_temp = Sample_id) %>% 
  separate(Sample_id_temp, into = c("Herbarium_code", "Species", "Sample_no"))

endocollections_map <- ggplot()+
  geom_sf(data = usa, fill = "white") +
  geom_point(data = subset(endo_herb_georef_seed, !is.na(Endo_status_liberal)), aes(x = lon, y = lat, color = Spp_code, pch = Spp_code), lwd = 1.5) +
  theme_minimal() + 
  scale_color_manual(values = c(   "#41b056", "#eda44e", "#4c90dd")) +
  scale_shape_manual(values = c(1, 2, 3)) +
  lims(x = c(-105,-72)) + 
  labs(x = c("Longitude"), y = c("Latitude"), color = "Species", pch = "Species") 

endocollections_map
ggsave(endocollections_map, filename = "~/Documents/endocollections_map.tiff")


endocollections_map_byherb <- ggplot()+
  geom_sf(data = usa, fill = "white") +
  geom_point(data = subset(endo_herb_georef_seed, Herbarium_code != "AKL" & Herbarium_code != "OLKA"), aes(x = lon, y = lat, color = Herbarium_code, pch = Spp_code), alpha = .5,lwd = 1.5) +
  theme_minimal() + 
  scale_color_brewer(palette = "Paired") +
  scale_shape_manual(values = c(1, 2, 3)) +
  lims(x = c(-105,-72)) + 
  labs(x = c("Longitude"), y = c("Latitude"), color = "Herbarium", pch = "Species") 
  
  endocollections_map_byherb
ggsave(endocollections_map_byherb, filename = "~/Documents/endocollections_map_byherb.tiff")


endocollections_map_byyear <- ggplot()+
  geom_sf(data = usa, fill = "white") +
  geom_point(data = endo_herb_georef_seed, aes(x = lon, y = lat, color = Spp_code, pch = Spp_code), lwd = 1.5) +
  facet_wrap(~year_bin)+
  theme_minimal() + 
  scale_color_manual(values = c(   "#41b056", "#eda44e", "#4c90dd")) +
  scale_shape_manual(values = c(1, 2, 3)) +
  lims(x = c(-105,-72)) + 
  labs(x = c("Longitude"), y = c("Latitude"), color = "Species", pch = "Species")

endocollections_map_byyear
ggsave(endocollections_map_byyear, filename = "~/Documents/endocollections_map_byyear.tiff", width = 8, height = 6)



# making an animated map
anim_herb_map <- endocollections_map +
  ggtitle('{closest_state}')+
  transition_states(year,
                    transition_length = 2,
                    state_length = 1)+
  shadow_mark()

gganimate::animate(anim_herb_map,
        nframes = length(unique(endo_herb_georef$year)),
        detail = 3,
        fps = 3,
        renderer = gifski_renderer(),
        height = 800, width = 1000, 
        res = 400, 
        end_pause = 3)
gganimate::anim_save("anim_herb_map.gif", nframes = length(unique(endo_herb_georef$year)),
                     detail = 3,
                     fps = 3,
                     anim_herb_map,
                     height = 800, 
                     width = 1000)




# Get annual prism climate data to make map of climate change magnitude
# There are several steps that take a bit of time, 1: downloading PRISM raster files, 2: extracting climate data from raster files
# Should have these PRISM files saved after running the first time in a folder called prismtmp (or whatever you want to name it)

# Uncomment to download raster files
get_prism_annual(type = c("ppt"), years=1895:2017, keepZip = F)
get_prism_annual(type = c("tmean"), years=1895:2017, keepZip = F)
get_prism_annual(type = c("tmin"), years=1895:2017, keepZip = F)
get_prism_annual(type = c("tmax"), years=1895:2017, keepZip = F)



ls_prism_data(name=TRUE)

# Grab the prism data and compile the files
climate_data <- ls_prism_data() %>%  
  prism_stack(.)  

# Extract projection coordinates from raster stack
climate_crs <- climate_data@crs@projargs

# Now we will extract the climate data from the raster files and save to a dataframe.
# This takes a while, so uncomment following lines and then you should save the output as an R data object

# point_df <- data.frame(rasterToPoints(climate_data)) ##creates a dataframe of points (This step takes a bit of time)
# year.df <- melt(point_df, c("x", "y")) %>%
#   separate(variable, into = c("PRISM", "Variable", "Label", "Resolution", "Year", "file")) %>%
#   rename("lon" = "x", "lat" = "y")
# saveRDS(year.df, file = "PRISM_climate_year_df.Rda")

# year.df <- readRDS(file = "PRISM_climate_year_df.Rda")


# Creating maps to visualize sample locations and overall change in climate
normals <- year.df %>% 
  mutate(time.period = case_when(Year>=1895 & Year <=1956 ~ "start",
                                 Year>=1957 & Year <=2017 ~ "end")) %>% 
  group_by(lon, lat, Variable, time.period) %>% 
  summarize(normal = mean(value))


normals_diff <- normals %>% 
  pivot_wider(names_from = "time.period", values_from = "normal") %>% 
  group_by(lon,lat,Variable) %>% 
  mutate(normal_diff = end - start)

ggplot()+
  geom_raster(data = subset(normals, Variable == "tmean" & time.period == "start"), aes(x = lon, y = lat, fill = normal)) +
  scale_fill_distiller("Change in Mean Temperature (\u00B0C)", palette = "Spectral") +#, limits = c(-1.15,1.15)) +
  theme_minimal() +coord_fixed(ratio=1.3) + 

ggplot()+
  geom_raster(data = subset(normals_diff, Variable == "tmean"), aes(x = lon, y = lat, fill = normal_diff)) +
  scale_fill_distiller("Change in Mean Temperature (\u00B0C)", palette = "Spectral", limits = c(-1.15,1.15)) +
  theme_minimal() +coord_fixed(ratio=1.3)

ggplot()+
  geom_raster(data = subset(normals_diff, Variable == "ppt"), aes(x = lon, y = lat, fill = normal_diff)) +
  scale_fill_distiller("Change in Annual Precipitation (mm.)", palette = "Spectral", trans = "reverse", limits = c(250,-250)) +
  theme_minimal() +coord_fixed(ratio=1.3)


# Map of collections with ppt or tmean

pptmap <- ggplot()+
  geom_raster(data = subset(normals, Variable == "ppt" & time.period == "end"), aes(x = lon, y = lat, fill = normal)) +
  geom_sf(data = usa, fill = "transparent", lwd = .2) +
  geom_point(data = subset(endo_herb_georef, !is.na(Spp_code) &Spp_code != "AGPE"), aes(x = lon, y = lat, shape = Spp_code), lwd = 1,alpha = .5) +
  theme_minimal() + labs(x = c("Longitude"), y = c("Latitude"), fill = c("Annual Precip. (mm.)"))+
  scale_shape_manual(values = c(1, 2, 3)) +
  scale_fill_viridis_b("Annual Precip. (mm.)",  limits = c( 0,2000)) + lims(x = c(-105,-65), y = c(25,45))

pptmap
ggsave(pptmap, filename = "~/Documents/pptmap.tiff")

tmeanmap <- ggplot()+
  geom_raster(data = subset(normals, Variable == "tmean" & time.period == "end"), aes(x = lon, y = lat, fill = normal)) +
  geom_sf(data = usa, fill = "transparent") +
  geom_point(data = endo_herb_georef, aes(x = lon, y = lat),shape = 1, lwd = .5,alpha = .5) +
  theme_minimal() + labs(x = c("Longitude"), y = c("Latitude"), fill = c("Mean Temp. (\u00B0C)"))+
  scale_fill_viridis("Mean Temp. (\u00B0C)", option = "magma") + lims(x = c(-108,-65))

tmeanmap
ggsave(tmeanmap, filename = "~/Documents/tmeanmap.tiff")





# texas precip maps

texaspptmap <- ggplot()+
  geom_raster(data = subset(normals, Variable == "ppt" & time.period == "end"), aes(x = lon, y = lat, fill = normal)) +
  geom_sf(data = texas, fill = "transparent", lwd = 1) +
  theme_minimal() + labs(x = c("Longitude"), y = c("Latitude")) + #, fill = c("Mean Annual Precip. (mm.)"))
  scale_fill_viridis_b(limits = c(0,2000))  + lims(x = c(-110,-90), y = c(26,40))

texaspptmap
ggsave(texaspptmap, filename = "~/Documents/texaspptmap.tiff")


pptchangemap <- ggplot()+
  geom_raster(data = subset(normals_diff, Variable == "ppt"), aes(x = lon, y = lat, fill = normal_diff)) +
  geom_sf(data = texas, fill = "transparent", lwd = 1) +
  theme_minimal() + labs(x = c("Longitude"), y = c("Latitude")) +#, fill = c("Change in Mean Annual Precip. (mm.)"))+
  scale_fill_viridis( option = "plasma", direction = -1, limits = c( -200,200)) + lims(x = c(-110,-90), y = c(26,40))

pptchangemap
ggsave(pptchangemap, filename = "~/Documents/pptchangenmap.tiff")





# Map of collections with change in temp and ppt

tempchangemap <- ggplot()+
  geom_raster(data = subset(normals_diff, Variable == "tmean"), aes(x = lon, y = lat, fill = normal_diff)) +
  geom_sf(data = usa, fill = "transparent") +
  geom_point(data = endo_herb_AGHY, aes(x = lon, y = lat), shape = 4, lwd = 2,alpha = .75) +
  theme_minimal() + labs(x = c("Longitude"), y = c("Latitude"), fill = c("Mean Temp. Change (\u00B0C)"))+
  scale_fill_distiller( palette = "Spectral", limits = c(-1.15,1.15)) + lims(x = c(-108,-65))

tempchangemap
ggsave(tempchangemap, filename = "~/Documents/tempchangemap2.tiff")



pptchangemap <- ggplot()+
  geom_raster(data = subset(normals_diff, Variable == "ppt"), aes(x = lon, y = lat, fill = normal_diff)) +
  geom_sf(data = usa, fill = "transparent") +
  geom_point(data = endo_herb_AGHY, aes(x = lon, y = lat), shape = 4, lwd = 2,alpha = .75) +
  theme_minimal() + labs(x = c("Longitude"), y = c("Latitude"), fill = c("Annual Precip. Change (mm.)"))+
  scale_fill_distiller("Annual Precip. Change (mm.)", palette = "Spectral", trans = "reverse", limits = c(150,-150)) + lims(x = c(-108,-65))
  
pptchangemap
ggsave(pptchangemap, filename = "~/Documents/pptchangemap2.tiff")


# # Now I will download monthly prism data 
# Since the monthly data come in a pretty gnarly shape, I will extract it at just the locations of our collections from the raster files.
# # For monthly, uncomment the following
# get_prism_monthlys(type = c("ppt"), years=1895:2017, mon = 1:12, keepZip = TRUE)
# get_prism_monthlys(type = c("tmean"), years=1895:2017, mon = 1:12, keepZip = TRUE)
# 
# Grab the prism data and compile the files
climate_data <- ls_prism_data() %>%  
  prism_stack(.)  

# Extract project coordinates from raster stack
climate_crs <- climate_data@crs@projargs

# We can get the lat/longs from our endophyte scores data frame "endo_herb_georef"
collection_sites <- endo_herb_georef %>% 
  dplyr::select(Sample_id, lat, lon) %>% 
  drop_na()

# Convert these locations to format that can be matched to Prism climate data
coordinates(collection_sites) <- c('lon', 'lat')
proj4string(collection_sites) <- CRS(climate_crs)

# Extract the  climate data from the raster stack for those sites 
climate_collections <- data.frame(coordinates(collection_sites), 
                           collection_sites$Sample_id, 
                           extract(climate_data, collection_sites))

# Reshape data. Col 1:3 are lat, long, and site ID. Col 4:ncol are climate data
# Column headers include date and climate type info
climate_collections <- climate_collections %>% 
  pivot_longer(cols = starts_with("PRISM"), names_to = "header", values_to = "value")


# The column header includes the date and data type, but also some other metadata that we don't need
# Here, I remove the extra info from the column header
climate_collections$header <- gsub('PRISM_', '', climate_collections$header) %>% 
  gsub('stable_4kmM3_', '', .) %>% 
  gsub('stable_4kmM2_', '', .) %>%
  gsub('_bil', '', .)

# Split header into type (precipitation or temperature), year, and month
climate_collections <- separate(climate_collections, 'header', 
                         into = c('variable', 'year'), 
                         sep = '_')
saveRDS(climate_collections, file = "climate_collections.Rda")


# Reshape data-- make a separate column for temperature and precipitation
climate_collections1 <- unique(climate_collections)
climate_collections1 <- climate_collections1 %>% 
  pivot_wider(names_from = variable, values_from = value) %>% 
  mutate(year = as.numeric(year)) %>%  # Make year numeric variables
  rename(Sample_id = collection_sites.Sample_id) %>% 
  arrange(Sample_id) # Order data by sample id


decadalmean <- climate_collections1 %>% 
  group_by(Sample_id, lon, lat) %>% 
  arrange(Sample_id,lon, lat, year) %>% 
  mutate(decade_ppt = slider::slide_dbl(ppt, mean, .before = 10, .after = 0),
         fiveyear_ppt = slider::slide_dbl(ppt, mean, .before = 5, .after = 0),
         decade_tmean = slider::slide_dbl(tmean, mean, .before = 10, .after = 0),
         fiveyear_tmean = slider::slide_dbl(tmean, mean, .before = 5, .after = 0),
         decade_tmax = slider::slide_dbl(tmax, mean, .before = 10, .after = 0),
         fiveyear_tmax = slider::slide_dbl(tmax, mean, .before = 5, .after = 0),
         decade_tmin = slider::slide_dbl(tmin, mean, .before = 10, .after = 0),
         fiveyear_tmin = slider::slide_dbl(tmin, mean, .before = 5, .after = 0)) %>% 
  ungroup() 


############################################################################
##### Analyzing climate and prevalence data ################################
############################################################################


# merge the climate variables with the endophyte scores

climate_endo_herb_georef <- endo_herb_georef %>% 
  left_join(decadalmean,
             by = c("Sample_id", "lon", "lat", "year")) # joins the dataframes based on the lat,long, and year 
# I need to update this file name
saveRDS(climate_endo_herb_georef, file = "climate__endo_herb_georef.Rda")
climate_endo_herb_georef <- readRDS(file = "climate__endo_herb_georef.Rda")


climate_AGHY <- climate_endo_herb_georef %>% 
  filter(grepl("AGHY", Sample_id))
binned_climate_AGHY <- climate_AGHY %>% 
  mutate(binned_ppt = cut(decade_ppt, breaks = 12), binned_year = cut(year, breaks = 3)) %>%  
  group_by(binned_ppt, binned_year) %>%   
  summarise(mean_ppt = mean(decade_ppt),
            mean_year = mean(year),
            mean_endo = mean(Endo_status_liberal),
            sample = n())

climate_ELVI <- climate_endo_herb_georef %>% 
  filter(grepl("ELVI", Sample_id))
binned_climate_ELVI <- climate_ELVI %>% 
  mutate(binned_ppt = cut(decade_ppt, breaks = 12), binned_year = cut(year, breaks = 3)) %>%  
  group_by(binned_ppt, binned_year) %>%   
  summarise(mean_ppt = mean(decade_ppt),
            mean_year = mean(year),
            mean_endo = mean(Endo_status_liberal),
            sample = n())
binned_tmean_ELVI <- climate_ELVI %>% 
  mutate(binned_tmean = cut(decade_tmean, breaks = 12), binned_year = cut(year, breaks = 3)) %>%  
  group_by(binned_tmean, binned_year) %>%   
  summarise(mean_tmean = mean(decade_tmean),
            mean_year = mean(year),
            mean_endo = mean(Endo_status_liberal),
            sample = n())



# AGHY
# model selection

aghy_climate_mod <- list()
aghy_climate_mod[[1]] <- 
ppt_mod <- glm(decade_tmax ~ lon +lat * year, data = climate_AGHY)
anova(ppt_mod, test = "Chisq")
summary(ppt_mod)
ggplot(data = climate_AGHY) +
  geom_point(aes(x = year, y = decade_ppt, color = lon))

AGHY_climate_mod <- glm(Endo_status_liberal ~  year +decade_ppt + decade_tmean, data = subset(climate_AGHY, climate_AGHY$lat < 35), family = binomial)
anova(AGHY_climate_mod, test = "Chisq")
summary(AGHY_climate_mod)

Anova(AGHY_climate_mod, test = "LR")


newdat1920 <- data.frame(decade_tmean = mean(climate_AGHY$decade_tmean, na.rm = T), decade_ppt = seq_range(min(climate_AGHY$decade_ppt, na.rm = T):max(climate_AGHY$decade_ppt, na.rm = T), n = length(climate_AGHY$decade_ppt)), year = 1920)
newdat1950 <- data.frame(decade_tmean = mean(climate_AGHY$decade_tmean, na.rm = T), decade_ppt = seq_range(min(climate_AGHY$decade_ppt, na.rm = T):max(climate_AGHY$decade_ppt, na.rm = T), n = length(climate_AGHY$decade_ppt)), year = 1950)
newdat2000 <- data.frame(decade_tmean = mean(climate_AGHY$decade_tmean, na.rm = T), decade_ppt = seq_range(min(climate_AGHY$decade_ppt, na.rm = T):max(climate_AGHY$decade_ppt, na.rm = T), n = length(climate_AGHY$decade_ppt)), year = 2000)
newdat <- rbind(newdat1920, newdat1950, newdat2000)
y_pred <- predict(AGHY_climate_mod, newdata = newdat, type = "response")
y_CI <- predict(AGHY_climate_mod, newdata = newdat, interval = "confidence", type = "link", se.fit=TRUE)
linkinv <- family(AGHY_climate_mod)$linkinv ## inverse-link function

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



AGHY_ppt <- ggplot() +
  geom_point(data = binned_climate_AGHY,aes(x = mean_ppt, y = mean_endo, size = sample , color = as.factor(binned_year))) +
  geom_line(data = newpred, aes(x = decade_ppt, y = pred, group = year, color = as.character(year))) +
  geom_ribbon(data = newpred, aes(x = decade_ppt, ymin = lwr, ymax = upr, group = year, fill = as.factor(year)), alpha = .2) +
  theme_classic() + labs(y = "Mean Endophyte Prevalence", x = "Ppt", color = "year")+
  scale_colour_manual(breaks = c("1920", "1950", "2000"),
                      values = c("#fc8d59", "#636363", "#91bfdb", "#fc8d69", "#635363", "#81bfdb")) +
  scale_fill_manual(breaks = c("1920", "1950", "2000"),
                    values = c("#fc8d59", "#636363", "#91bfdb")) + guides(fill = FALSE)

AGHY_ppt
ggsave(AGHY_ppt, filename = "~/Documents/AGHYppt.tiff", width = 4, height = 3)

#  ELVI


ELVI_climate_mod <- glm(Endo_status_liberal ~  year + decade_ppt * decade_tmean, data = climate_ELVI, family = binomial)
anova(ELVI_climate_mod, test = "Chisq")
summary(ELVI_climate_mod)

Anova(ELVI_climate_mod, test = "LR")


newdat1920 <- data.frame(decade_tmean = mean(climate_ELVI$decade_tmean, na.rm = T), decade_ppt = seq_range(min(climate_ELVI$decade_ppt, na.rm = T):max(climate_ELVI$decade_ppt, na.rm = T), n = length(climate_ELVI$decade_ppt)), year = 1920)
newdat1950 <- data.frame(decade_tmean = mean(climate_ELVI$decade_tmean, na.rm = T), decade_ppt = seq_range(min(climate_ELVI$decade_ppt, na.rm = T):max(climate_ELVI$decade_ppt, na.rm = T), n = length(climate_ELVI$decade_ppt)), year = 1950)
newdat2000 <- data.frame(decade_tmean = mean(climate_ELVI$decade_tmean, na.rm = T), decade_ppt = seq_range(min(climate_ELVI$decade_ppt, na.rm = T):max(climate_ELVI$decade_ppt, na.rm = T), n = length(climate_ELVI$decade_ppt)), year = 2000)
newdat <- rbind(newdat1920, newdat1950, newdat2000)
y_pred <- predict(ELVI_climate_mod, newdata = newdat, type = "response")
y_CI <- predict(ELVI_climate_mod, newdata = newdat, interval = "confidence", type = "link", se.fit=TRUE)
linkinv <- family(ELVI_climate_mod)$linkinv ## inverse-link function

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



ELVI_ppt <- ggplot() +
  geom_point(data = binned_climate_ELVI,aes(x = mean_ppt, y = mean_endo, size = sample , color = as.factor(binned_year))) +
  geom_line(data = newpred, aes(x = decade_ppt, y = pred, group = year, color = as.character(year))) +
  geom_ribbon(data = newpred, aes(x = decade_ppt, ymin = lwr, ymax = upr, group = year, fill = as.factor(year)), alpha = .2) +
  theme_classic() + labs(y = "Mean Endophyte Prevalence", x = "Ppt", color = "year")+
  scale_colour_manual(breaks = c("1920", "1950", "2000"),
                      values = c("#fc8d59", "#636363", "#91bfdb", "#fc8d69", "#635363", "#81bfdb")) +
  scale_fill_manual(breaks = c("1920", "1950", "2000"),
                    values = c("#fc8d59", "#636363", "#91bfdb")) + guides(fill = FALSE)

ELVI_ppt
ggsave(ELVI_ppt, filename = "~/Documents/ELVIppt.tiff", width = 4, height = 3)



newdat1920 <- data.frame(decade_tmean = seq_range(min(climate_ELVI$decade_tmean, na.rm = T):max(climate_ELVI$decade_tmean, na.rm = T),n = length(climate_ELVI$decade_tmean)), decade_ppt = mean(climate_ELVI$decade_ppt, na.rm = T), year = 1920)
newdat1950 <- data.frame(decade_tmean = seq_range(min(climate_ELVI$decade_tmean, na.rm = T):max(climate_ELVI$decade_tmean, na.rm = T),n = length(climate_ELVI$decade_tmean)), decade_ppt = mean(climate_ELVI$decade_ppt, na.rm = T), year = 1950)
newdat2000 <- data.frame(decade_tmean = seq_range(min(climate_ELVI$decade_tmean, na.rm = T):max(climate_ELVI$decade_tmean, na.rm = T),n = length(climate_ELVI$decade_tmean)), decade_ppt = mean(climate_ELVI$decade_ppt, na.rm = T), year = 2000)
newdat <- rbind(newdat1920, newdat1950, newdat2000)
y_pred <- predict(ELVI_climate_mod, newdata = newdat, type = "response")
y_CI <- predict(ELVI_climate_mod, newdata = newdat, interval = "confidence", type = "link", se.fit=TRUE)
linkinv <- family(ELVI_climate_mod)$linkinv ## inverse-link function

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



ELVI_tmean <- ggplot() +
  geom_point(data = binned_tmean_ELVI,aes(x = mean_tmean, y = mean_endo, size = sample , color = as.factor(binned_year))) +
  geom_line(data = newpred, aes(x = decade_tmean, y = pred, group = year, color = as.character(year))) +
  geom_ribbon(data = newpred, aes(x = decade_tmean, ymin = lwr, ymax = upr, group = year, fill = as.factor(year)), alpha = .2) +
  theme_classic() + labs(y = "Mean Endophyte Prevalence", x = "Tmean", color = "year")+
  scale_colour_manual(breaks = c("1920", "1950", "2000"),
                      values = c("#fc8d59", "#636363", "#91bfdb", "#fc8d69", "#635363", "#81bfdb")) +
  scale_fill_manual(breaks = c("1920", "1950", "2000"),
                    values = c("#fc8d59", "#636363", "#91bfdb")) + guides(fill = FALSE)

ELVI_tmean
ggsave(ELVI_tmean, filename = "~/Documents/ELVItmean.tiff", width = 4, height = 3)



############################################################################
##### Analyzing fitness data ###############################################
############################################################################

aghy_fitness_data <- endo_herb_AGHY %>% 
  filter(!is.na(surface_area_cm2))
plot(aghy_fitness_data$year, aghy_fitness_data$surface_area_cm2)
plot(aghy_fitness_data$lon, aghy_fitness_data$surface_area_cm2)
plot( aghy_fitness_data$Endo_status_liberal, aghy_fitness_data$surface_area_cm2)


area_mod <- glm(surface_area_cm2~ lat+lon+year+Endo_status_liberal, data = aghy_fitness_data, family = gaussian)
area_mod <- glm(surface_area_cm2~ year, data = aghy_fitness_data, family = gaussian)

anova(area_mod, test = "Chisq")
summary(area_mod)


# ELVI

elvi_fitness_data <- endo_herb_ELVI %>% 
  filter(!is.na(surface_area_cm2))

plot(elvi_fitness_data$year, elvi_fitness_data$mean_infl_length)
plot(elvi_fitness_data$lon, elvi_fitness_data$mean_infl_length)
plot(elvi_fitness_data$lat, elvi_fitness_data$mean_infl_length)

plot(elvi_fitness_data$Endo_status_liberal, elvi_fitness_data$mean_infl_length)

plot(elvi_fitness_data$year, elvi_fitness_data$infl_count)
plot(elvi_fitness_data$lon, elvi_fitness_data$infl_count)
plot(elvi_fitness_data$lat, elvi_fitness_data$infl_count)

plot(elvi_fitness_data$Endo_status_liberal, elvi_fitness_data$infl_count)

plot(elvi_fitness_data$year,  elvi_fitness_data$surface_area_cm2)
plot(elvi_fitness_data$lon, elvi_fitness_data$surface_area_cm2)
plot(elvi_fitness_data$Endo_status_liberal, elvi_fitness_data$surface_area_cm2)


plot(elvi_fitness_data$surface_area_cm2, elvi_fitness_data$mean_infl_length)
plot(elvi_fitness_data$surface_area_cm2, elvi_fitness_data$infl_count)


# model selection
elvi_models <- list()
elvi_models[[1]] <- glm(mean_infl_length ~ lon, data = elvi_fitness_data, family = gaussian)
elvi_models[[2]] <- glm(mean_infl_length ~ lat, data = elvi_fitness_data, family = gaussian)
elvi_models[[3]] <- glm(mean_infl_length ~ year, data = elvi_fitness_data, family = gaussian)
elvi_models[[4]] <- glm(mean_infl_length ~ Endo_status_liberal, data = elvi_fitness_data, family = gaussian)
elvi_models[[5]] <- glm(mean_infl_length ~ lon + lat, data = elvi_fitness_data, family = gaussian)
elvi_models[[6]] <- glm(mean_infl_length ~ lon * lat, data = elvi_fitness_data, family = gaussian)
elvi_models[[7]] <- glm(mean_infl_length ~ lon + year, data = elvi_fitness_data, family = gaussian)
elvi_models[[8]] <- glm(mean_infl_length ~ lon * year, data = elvi_fitness_data, family = gaussian)
elvi_models[[9]] <- glm(mean_infl_length ~ lon + Endo_status_liberal, data = elvi_fitness_data, family = gaussian)
elvi_models[[10]] <- glm(mean_infl_length ~ lon * Endo_status_liberal, data = elvi_fitness_data, family = gaussian)
elvi_models[[11]] <- glm(mean_infl_length ~ lat + year, data = elvi_fitness_data, family = gaussian)
elvi_models[[12]] <- glm(mean_infl_length ~ lat * year, data = elvi_fitness_data, family = gaussian)
elvi_models[[13]] <- glm(mean_infl_length ~ lat + Endo_status_liberal, data = elvi_fitness_data, family = gaussian)
elvi_models[[14]] <- glm(mean_infl_length ~ lat * Endo_status_liberal, data = elvi_fitness_data, family = gaussian)
elvi_models[[15]] <- glm(mean_infl_length ~ lat + lon + year, data = elvi_fitness_data, family = gaussian)
elvi_models[[16]] <- glm(mean_infl_length ~ lat * lon + year, data = elvi_fitness_data, family = gaussian)
elvi_models[[17]] <- glm(mean_infl_length ~ lat + lon * year, data = elvi_fitness_data, family = gaussian)
elvi_models[[18]] <- glm(mean_infl_length ~ lon + lat * year, data = elvi_fitness_data, family = gaussian)
elvi_models[[19]] <- glm(mean_infl_length ~ lat * lon * year, data = elvi_fitness_data, family = gaussian)

elvi_models[[20]] <- glm(mean_infl_length ~ lat + lon + Endo_status_liberal, data = elvi_fitness_data, family = gaussian)
elvi_models[[21]] <- glm(mean_infl_length ~ lat * lon + Endo_status_liberal, data = elvi_fitness_data, family = gaussian)
elvi_models[[22]] <- glm(mean_infl_length ~ lat + lon * Endo_status_liberal, data = elvi_fitness_data, family = gaussian)
elvi_models[[23]] <- glm(mean_infl_length ~ lon + lat * Endo_status_liberal, data = elvi_fitness_data, family = gaussian)
elvi_models[[24]] <- glm(mean_infl_length ~ lat * lon * Endo_status_liberal, data = elvi_fitness_data, family = gaussian)

elvi_models[[25]] <- glm(mean_infl_length ~ Endo_status_liberal + lon + year, data = elvi_fitness_data, family = gaussian)
elvi_models[[26]] <- glm(mean_infl_length ~ Endo_status_liberal * lon + year, data = elvi_fitness_data, family = gaussian)
elvi_models[[27]] <- glm(mean_infl_length ~ Endo_status_liberal + lon * year, data = elvi_fitness_data, family = gaussian)
elvi_models[[28]] <- glm(mean_infl_length ~ lon + Endo_status_liberal * year, data = elvi_fitness_data, family = gaussian)
elvi_models[[29]] <- glm(mean_infl_length ~ Endo_status_liberal * lon * year, data = elvi_fitness_data, family = gaussian)


elvi_models[[30]] <- glm(mean_infl_length ~ lat + Endo_status_liberal + year, data = elvi_fitness_data, family = gaussian)
elvi_models[[31]] <- glm(mean_infl_length ~ lat * Endo_status_liberal + year, data = elvi_fitness_data, family = gaussian)
elvi_models[[32]] <- glm(mean_infl_length ~ lat + Endo_status_liberal * year, data = elvi_fitness_data, family = gaussian)
elvi_models[[33]] <- glm(mean_infl_length ~ Endo_status_liberal + lat * year, data = elvi_fitness_data, family = gaussian)
elvi_models[[34]] <- glm(mean_infl_length ~ lat * Endo_status_liberal * year, data = elvi_fitness_data, family = gaussian)

elvi_models[[35]] <- glm(mean_infl_length ~ lat + lon + year+ Endo_status_liberal, data = elvi_fitness_data, family = gaussian)
elvi_models[[36]] <- glm(mean_infl_length ~ lat * lon + year+ Endo_status_liberal, data = elvi_fitness_data, family = gaussian)
elvi_models[[37]] <- glm(mean_infl_length ~ lat + lon * year+ Endo_status_liberal, data = elvi_fitness_data, family = gaussian)
elvi_models[[38]] <- glm(mean_infl_length ~ lon + lat * year+ Endo_status_liberal, data = elvi_fitness_data, family = gaussian)
elvi_models[[39]] <- glm(mean_infl_length ~ lat * lon * year+ Endo_status_liberal, data = elvi_fitness_data, family = gaussian)
elvi_models[[40]] <- glm(mean_infl_length ~ lat * lon * year* Endo_status_liberal, data = elvi_fitness_data, family = gaussian)

elvi_models[[41]] <- glm(mean_infl_length ~ 1, data = elvi_fitness_data, family = gaussian)

AICtab(elvi_models) # best model according to AIC is 14, lat * endo_status





infl_mod <- glm(mean_infl_length ~ lon*Endo_status_liberal, data = elvi_fitness_data, family = gaussian)
anova(infl_mod, test = "Chisq")
summary(infl_mod)

# plot(infl_mod)


newdat1<- data.frame(lon = seq(-100,-70,1), Endo_status_liberal = 1)
newdat0 <- data.frame(lon = seq(-100,-70,1), Endo_status_liberal = 0)
newdat <- rbind(newdat0, newdat1)
y_pred <- predict(infl_mod, newdata = newdat, type = "response")
y_CI <- predict(infl_mod, newdata = newdat, interval = "confidence", type = "link", se.fit=TRUE)
linkinv <- family(infl_mod)$linkinv ## inverse-link function

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



ELVI_outcome <- ggplot() +
  geom_point(data = elvi_fitness_data, aes(x = lon, y =mean_infl_length, pch = as.factor(Endo_status_liberal)), color = "#4c90dd", lwd = 2)+
  geom_line(data = newpred, aes(x = lon, y = pred, group = Endo_status_liberal, linetype = as.factor(Endo_status_liberal)), color = "#4c90dd") +
  geom_ribbon(data = newpred, aes(x = lon, ymin = lwr, ymax = upr, group = Endo_status_liberal, fill = as.factor(Endo_status_liberal)), alpha = .2) +
  theme_classic() + labs(y = "Mean Inflorescence Length", x = "Longitude", linetype = "Endophyte_status", shape = "Endophyte_status")+
  scale_linetype_manual(breaks = c("1", "0"),
              values = c("solid", "dashed"))+
  scale_shape_manual(breaks = c("1", "0"),
                     values = c(20,1))+
  scale_fill_manual(breaks = c("1", "0"),
                    values = c(  "#4c90dd", "#636363")) +
  xlim(-100,-70) + guides(fill = FALSE)

ELVI_outcome
ggsave(ELVI_outcome, filename = "~/Documents/ELVIinfl.tiff", width = 4, height = 3)




# Endophyte effect plot
infl_mod <- glm(mean_infl_length ~ lon*Endo_status_liberal, data = (subset(elvi_fitness_data)), family = gaussian)
anova(infl_mod, test = "Chisq")
summary(infl_mod)

# plot(infl_mod)




newdat1<- data.frame(lon = seq(-100,-70,1), Endo_status_liberal = 1)
newdat0 <- data.frame(lon = seq(-100,-70,1),  Endo_status_liberal = 0)
newdat <- rbind(newdat1,newdat0)
y_pred <- predict(infl_mod, newdata = newdat, type = "response")
y_CI <- predict(infl_mod, newdata = newdat, interval = "confidence", type = "link", se.fit=TRUE)
linkinv <- family(infl_mod)$linkinv ## inverse-link function

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
newpred_endo <- newpred %>% 
  dplyr::select(lon, Endo_status_liberal, pred, lwr, upr) %>% 
  pivot_wider(names_from = Endo_status_liberal, values_from = c(pred, lwr, upr)) %>% 
  mutate(endo_effect = pred_1-pred_0,
         endo_lwr = lwr_1 - lwr_0,
         endo_upr = upr_1 - upr_0)

ELVI_outcomelon <- ggplot() +
  geom_line(data = newpred_endo, aes(x= lon, y = endo_effect)) +
  geom_ribbon(data = newpred_endo, aes(x = lon, ymin = endo_lwr, ymax = endo_upr), alpha = .2) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  theme_classic() + labs(y = "Endophyte Effect on Mean Inflorescence Length", x = "Longitude", linetype = "Endophyte_status", shape = "Endophyte_status")+
  scale_color_manual(breaks = c("1930", "2000"),
                     values = c(  "#636363", "#4c90dd")) +
  scale_fill_manual(breaks = c("1930", "2000"),
                    values = c(  "#636363", "#4c90dd")) 
ELVI_outcomelon

ggsave(ELVI_outcomelon, filename = "~/Documents/ELVIinfleffect.tiff", width = 4, height = 3)

# model for inflorescence count
inflcount_mod <- glm(infl_count ~  lon*Endo_status_liberal, data = elvi_fitness_data, family = poisson)

anova(inflcount_mod, test = "Chisq")
summary(inflcount_mod)




newdat1<- data.frame(lon = seq(-100,-70,1), Endo_status_liberal = 1)
newdat0 <- data.frame(lon = seq(-100,-70,1), Endo_status_liberal = 0)
newdat <- rbind(newdat0, newdat1)
y_pred <- predict(inflcount_mod, newdata = newdat, type = "response")
y_CI <- predict(inflcount_mod, newdata = newdat, interval = "confidence", type = "link", se.fit=TRUE)
linkinv <- family(inflcount_mod)$linkinv ## inverse-link function

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




ELVI_inflcount_outcome <- ggplot() +
  geom_point(data = elvi_fitness_data, aes(x = lon, y =infl_count, pch = as.factor(Endo_status_liberal)), color = "#4c90dd", lwd = 2)+
  geom_line(data = newpred, aes(x = lon, y = pred, group = Endo_status_liberal, linetype = as.factor(Endo_status_liberal)), color = "#4c90dd") +
  geom_ribbon(data = newpred, aes(x = lon, ymin = lwr, ymax = upr, group = Endo_status_liberal, fill = as.factor(Endo_status_liberal)), alpha = .2) +
  theme_classic() + labs(y = "Inflorescence Count", x = "Longitude", linetype = "Endophyte_status", shape = "Endophyte_status")+
  scale_linetype_manual(breaks = c("1", "0"),
                        values = c("solid", "dashed"))+
  scale_shape_manual(breaks = c("1", "0"),
                     values = c(20,1))+
  scale_fill_manual(breaks = c("1", "0"),
                    values = c(  "#4c90dd", "#636363")) +
  xlim(-100,-70) + guides(fill = FALSE)

ELVI_inflcount_outcome
ggsave(ELVI_inflcount_outcome, filename = "~/Documents/ELVIinflcount.tiff", width = 4, height = 3)


# Endophyte effect plot
inflcount_mod <- glm(infl_count ~ lon*Endo_status_liberal, data = (subset(elvi_fitness_data)), family = poisson)
anova(inflcount_mod, test = "Chisq")
summary(inflcount_mod)

# plot(infl_mod)




newdat1<- data.frame(lon = seq(-100,-70,1), Endo_status_liberal = 1)
newdat0 <- data.frame(lon = seq(-100,-70,1),  Endo_status_liberal = 0)
newdat <- rbind(newdat1,newdat0)
y_pred <- predict(inflcount_mod, newdata = newdat, type = "response")
y_CI <- predict(inflcount_mod, newdata = newdat, interval = "confidence", type = "link", se.fit=TRUE)
linkinv <- family(inflcount_mod)$linkinv ## inverse-link function

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
newpred_endo <- newpred %>% 
  dplyr::select(lon, Endo_status_liberal, pred, lwr, upr) %>% 
  pivot_wider(names_from = Endo_status_liberal, values_from = c(pred, lwr, upr)) %>% 
  mutate(endo_effect = pred_1-pred_0,
         endo_lwr = lwr_1 - lwr_0,
         endo_upr = upr_1 - upr_0)

ELVI_outcomecountlon <- ggplot() +
  geom_line(data = newpred_endo, aes(x= lon, y = endo_effect)) +
  geom_ribbon(data = newpred_endo, aes(x = lon, ymin = endo_lwr, ymax = endo_upr), alpha = .2) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  theme_classic() + labs(y = "Endophyte Effect on Inflorescence Count", x = "Longitude", linetype = "Endophyte_status", shape = "Endophyte_status")+
  scale_color_manual(breaks = c("1930", "2000"),
                     values = c(  "#636363", "#4c90dd")) +
  scale_fill_manual(breaks = c("1930", "2000"),
                    values = c(  "#636363", "#4c90dd")) 
ELVI_outcomecountlon

ggsave(ELVI_outcomecountlon, filename = "~/Documents/ELVIinflcounteffect.tiff", width = 4, height = 3)






# Looking at change over time in interaction outcome

infl_mod <- glm(mean_infl_length ~ lon*year*Endo_status_liberal, data = (subset(elvi_fitness_data)), family = gaussian)
anova(infl_mod, test = "Chisq")
summary(infl_mod)

# plot(infl_mod)



newdat1_1930<- data.frame(lon = seq(-100,-70,1), year = 1930, Endo_status_liberal = 1)
newdat0_1930 <- data.frame(lon = seq(-100,-70,1), year = 1930, Endo_status_liberal = 0)
newdat1_2000<- data.frame(lon = seq(-100,-70,1), year = 2000, Endo_status_liberal = 1)
newdat0_2000 <- data.frame(lon = seq(-100,-70,1), year = 2000, Endo_status_liberal = 0)
newdat <- rbind(newdat0_1930, newdat1_1930, newdat0_2000, newdat1_2000)
y_pred <- predict(infl_mod, newdata = newdat, type = "response")
y_CI <- predict(infl_mod, newdata = newdat, interval = "confidence", type = "link", se.fit=TRUE)
linkinv <- family(infl_mod)$linkinv ## inverse-link function

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
newpred_endo <- newpred %>% 
  dplyr::select(year,lon, Endo_status_liberal, pred, lwr, upr) %>% 
  pivot_wider(names_from = Endo_status_liberal, values_from = c(pred, lwr, upr)) %>% 
  mutate(endo_effect = pred_1-pred_0,
         endo_lwr = lwr_1 - lwr_0,
         endo_upr = upr_1 - upr_0)

ELVI_outcomeyear <- ggplot() +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50")+
  geom_line(data = newpred_endo, aes(x= lon, y = endo_effect, group = year, color = as.factor(year))) +
  geom_ribbon(data = newpred_endo, aes(x = lon, ymin = endo_lwr, ymax = endo_upr, group = year, fill = as.factor(year)), alpha = .2) +
  theme_classic() + labs(y = "Endophyte Effect on Mean Inflorescence Length", x = "Longitude", linetype = "Endophyte_status", shape = "Endophyte_status")+
  scale_color_manual(breaks = c("1930", "2000"),
                    values = c(  "black", "#4c90dd")) +
  scale_fill_manual(breaks = c("1930", "2000"),
                  values = c(  "#636363", "#4c90dd")) 
ELVI_outcomeyear
  
ggsave(ELVI_outcomeyear, filename = "~/Documents/ELVIoutcomeyear.tiff", width = 4.5, height = 3)
# geom_line(data = newpred, aes(x = lon, y = pred, group = Endo_status_liberal, linetype = as.factor(Endo_status_liberal)), color = "#4c90dd") 
#   geom_ribbon(data = newpred, aes(x = year, ymin = lwr, ymax = upr, group = Endo_status_liberal, fill = as.factor(Endo_status_liberal)), alpha = .2) +
#   theme_classic() + labs(y = "Mean Inflorescence Length", x = "Year", linetype = "Endophyte_status", shape = "Endophyte_status")+
#   scale_linetype_manual(breaks = c("1", "0"),
#                         values = c("solid", "dashed"))+
#   scale_shape_manual(breaks = c("1", "0"),
#                      values = c(20,1))+
#   scale_fill_manual(breaks = c("1", "0"),
#                     values = c(  "#4c90dd", "#636363")) +
#   xlim(-100,-70) + guides(fill = FALSE)

ELVI_outcomeyear
ggsave(ELVI_outcomeyear, filename = "~/Documents/ELVIinflyear.tiff", width = 4, height = 3)



# model for inflorescence count
inflcount_mod <- glm(infl_count ~  year*Endo_status_liberal, data = elvi_fitness_data, family = poisson)

anova(inflcount_mod, test = "Chisq")
summary(inflcount_mod)



newdat1<- data.frame(year = seq(1900,2016,1), Endo_status_liberal = 1)
newdat0 <- data.frame(year= seq(1900,2016,1), Endo_status_liberal = 0)
newdat <- rbind(newdat0, newdat1)
y_pred <- predict(inflcount_mod, newdata = newdat, type = "response")
y_CI <- predict(inflcount_mod, newdata = newdat, interval = "confidence", type = "link", se.fit=TRUE)
linkinv <- family(inflcount_mod)$linkinv ## inverse-link function

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




ELVI_inflcount_outcomeyear <- ggplot() +
  geom_point(data = elvi_fitness_data, aes(x = year, y =infl_count, pch = as.factor(Endo_status_liberal)), color = "#4c90dd", lwd = 2)+
  geom_line(data = newpred, aes(x = year, y = pred, group = Endo_status_liberal, linetype = as.factor(Endo_status_liberal)), color = "#4c90dd") +
  geom_ribbon(data = newpred, aes(x = year, ymin = lwr, ymax = upr, group = Endo_status_liberal, fill = as.factor(Endo_status_liberal)), alpha = .2) +
  theme_classic() + labs(y = "Inflorescence Count", x = "Year", linetype = "Endophyte_status", shape = "Endophyte_status")+
  scale_linetype_manual(breaks = c("1", "0"),
                        values = c("solid", "dashed"))+
  scale_shape_manual(breaks = c("1", "0"),
                     values = c(20,1))+
  scale_fill_manual(breaks = c("1", "0"),
                    values = c(  "#4c90dd", "#636363"))+
  xlim(1910,2016) + ylim(0, 14) + guides(fill = FALSE)

ELVI_inflcount_outcomeyear
ggsave(ELVI_inflcount_outcomeyear, filename = "~/Documents/ELVIinflcountyear.tiff", width = 4, height = 3)



# Looking at change over time in interaction outcome

inflcount_modyear <- glm(infl_count ~ lon*year*Endo_status_liberal, data = (subset(elvi_fitness_data)), family = poisson)
anova(inflcount_modyear, test = "Chisq")
summary(inflcount_modyear)

# plot(infl_mod)



newdat1_1930<- data.frame(lon = seq(-100,-70,1), year = 1930, Endo_status_liberal = 1)
newdat0_1930 <- data.frame(lon = seq(-100,-70,1), year = 1930, Endo_status_liberal = 0)
newdat1_2000<- data.frame(lon = seq(-100,-70,1), year = 2000, Endo_status_liberal = 1)
newdat0_2000 <- data.frame(lon = seq(-100,-70,1), year = 2000, Endo_status_liberal = 0)
newdat <- rbind(newdat0_1930, newdat1_1930, newdat0_2000, newdat1_2000)
y_pred <- predict(inflcount_modyear, newdata = newdat, type = "response")
y_CI <- predict(inflcount_modyear, newdata = newdat, interval = "confidence", type = "link", se.fit=TRUE)
linkinv <- family(inflcount_modyear)$linkinv ## inverse-link function

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
newpred_endo <- newpred %>% 
  dplyr::select(year,lon, Endo_status_liberal, pred, lwr, upr) %>% 
  pivot_wider(names_from = Endo_status_liberal, values_from = c(pred, lwr, upr)) %>% 
  mutate(endo_effect = pred_1-pred_0,
         endo_lwr = lwr_1 - lwr_0,
         endo_upr = upr_1 - upr_0)

ELVI_countoutcomeyear <- ggplot() +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50")+
  geom_line(data = newpred_endo, aes(x= lon, y = endo_effect, group = year, color = as.factor(year))) +
  geom_ribbon(data = newpred_endo, aes(x = lon, ymin = endo_lwr, ymax = endo_upr, group = year, fill = as.factor(year)), alpha = .2) +
  theme_classic() + labs(y = "Endophyte Effect on Mean Inflorescence Length", x = "Longitude", linetype = "Endophyte_status", shape = "Endophyte_status")+
  scale_color_manual(breaks = c("1930", "2000"),
                     values = c(  "black", "#4c90dd")) +
  scale_fill_manual(breaks = c("1930", "2000"),
                    values = c(  "#636363", "#4c90dd")) +
  coord_cartesian(ylim = c(-25, 5)) 
ELVI_countoutcomeyear

ggsave(ELVI_countoutcomeyear, filename = "~/Documents/ELVIcountoutcomeyear.tiff", width = 4.5, height = 3)



#########################################################################################################
# Bernoulli GLM for endo ~  lon + year + lon*year   -------------------------
########################################################################################
# Messing around with Stan models

invlogit<-function(x){exp(x)/(1+exp(x))}
logit = function(x) { log(x/(1-x)) }


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


