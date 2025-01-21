# Title: Digital Herbarium Records and Endophyte Sample Merging
# Purpose: Imports and merges downloaded digitized herbarium records with the Endo_Herbarium database
# Authors: Joshua Fowler, Mallory Tucker, Ella Segal, Lani Dufresne, and Tom Miller
# Updated: Nov 4, 2022

library(tidyverse) # for data manipulation and ggplot
# library(slider) # add on to tidyverse to calculate sliding windows for climate data
library(fuzzyjoin) # add on to tidyverse to merge tables on nearest values
library(readxl)
library(lubridate)
library(ggmap)
library(prism) # to import prism raster files
library(terra)
library(sf)
library(INLA) # using this here to work with the shape boundary
library(inlabru) # using this here for some plotting stuff

################################################################################
############ Read in digitized herbarium records ############################### 
################################################################################

# UT Austin downloaded from TORCH
AGHY_UTAustin <- read_csv(file = "~/Dropbox/Josh&Tom - shared/Endo_Herbarium/DigitizedHerbariumRecords/UTAustinAGHYrecords/occurrences.csv") %>%
  mutate(new_id = gsub("[a-zA-Z ]", "", catalogNumber)) %>%   #creates a new id that we will use to merge with the
  mutate(municipality = as.character(municipality)) %>% 
  dplyr::select(id, catalogNumber, country, stateProvince, county, municipality, locality, decimalLatitude, decimalLongitude, coordinateUncertaintyInMeters, recordedBy, recordNumber, new_id, eventDate, day, month, year) %>%
  mutate(Spp_code = "AGHY")
ELVI_UTAustin <- read_csv(file = "~/Dropbox/Josh&Tom - shared/Endo_Herbarium/DigitizedHerbariumRecords/UTAustinELVIrecords/occurrences.csv") %>% 
  mutate(new_id = gsub("[a-zA-Z ]", "", catalogNumber)) %>%    #creates a new id that we will use to merge with the 
  mutate(municipality = as.character(municipality)) %>% 
  dplyr::select(id, catalogNumber, country, stateProvince, county, municipality, locality, decimalLatitude, decimalLongitude, coordinateUncertaintyInMeters, recordedBy, recordNumber, new_id, eventDate, day, month, year) %>% 
  mutate(Spp_code = "ELVI")
AGPE_UTAustin <-   read_csv(file = "~/Dropbox/Josh&Tom - shared/Endo_Herbarium/DigitizedHerbariumRecords/UTAustinAGPErecords/occurrences.csv") %>% 
  mutate(new_id = gsub("[a-zA-Z ]", "", catalogNumber)) %>%    #creates a new id that we will use to merge with the 
  mutate(municipality = as.character(municipality)) %>% 
  dplyr::select(id, catalogNumber, country, stateProvince, county, municipality, locality, decimalLatitude, decimalLongitude, coordinateUncertaintyInMeters, recordedBy, recordNumber, new_id, eventDate, day, month, year) %>% 
  mutate(Spp_code = "AGPE")

UTAustin_torch <- rbind(AGHY_UTAustin, ELVI_UTAustin, AGPE_UTAustin) %>% 
  mutate(primary_collector = case_when(str_detect(recordedBy, "Jr.") ~ word(recordedBy, start = 1, end = 2, sep = fixed(",")),
                                       !str_detect(recordedBy, "Jr.") ~ word(recordedBy, start = 1, sep = fixed(","))),
         primary_collector = word(primary_collector, sep = fixed("with")),
         primary_collector = word(primary_collector, sep = fixed("and")),
         primary_collector = word(primary_collector, sep = fixed("&")),
         primary_collector = word(primary_collector, sep = fixed("asst")),
         primary_collector = word(primary_collector, sep = fixed("assit")),
         primary_collector = word(primary_collector, sep = fixed("et")),
         primary_collector = word(primary_collector, sep = fixed("|")),
         primary_collector = str_replace_all(primary_collector, "�", " "),
         primary_collector = case_when(primary_collector == "Jes s Vald s" ~ "Jesus Valdes", 
                                       primary_collector == "Pedro Tenorio L." ~ "Pedro Tenorio", 
                                       primary_collector == "Roger W. S" ~ "Roger W. Sanders",
                                       TRUE ~ primary_collector)) %>% 
  mutate(collector_firstname = word(primary_collector),
         collector_lastname = case_when(str_detect(primary_collector,"Jr.") ~ word(str_trim(primary_collector), start = -2, end = -1),
                                 !str_detect(primary_collector,"Jr.") ~ word(str_trim(primary_collector), -1)))

# Texas A&M digitized records (Includes both AGHY and ELVI)
AM_records <- read_csv(file = "~/Dropbox/Josh&Tom - shared/Endo_Herbarium/DigitizedHerbariumRecords/TexasA&M_digitized_records/Fowler Data.csv") %>% 
  unite(Institution_specimen_id, CollectionCode, id, sep = "") %>% 
  separate(DateCollected, into = c("year", "month", "day"), remove = FALSE) %>% 
  mutate(DateCollected = as_date(DateCollected)) %>% 
  mutate(primary_collector = word(collector, sep = fixed("|")),
         primary_collector = case_when(primary_collector == "Michael M. " ~ "Michael M. MacRoberts", TRUE ~ primary_collector),
         collector_firstname = word(primary_collector),
         collector_lastname = case_when(str_detect(primary_collector,"Jr.") ~ word(str_trim(primary_collector), start = -2, end = -1),
                                        str_detect(primary_collector, "III") ~ word(str_trim(primary_collector),  start = -2, end = -1),
                                        !str_detect(primary_collector,"Jr.") | !str_detect(primary_collector, "III") ~ word(str_trim(primary_collector), -1))) %>%
  mutate(collector_lastname = gsub( "[()]", "", collector_lastname)) %>% 
  mutate(year = as.numeric(year), month = as.numeric(month), day = as.numeric(day)) %>% 
  mutate(year = case_when(grepl("Spring, 1942", VerbatimDateCollected) ~ 1942,
                          TRUE ~ year)) %>% 
  mutate(Spp_code = case_when(species == "hyemalis" ~ "AGHY", species == "hiemalis" ~ "AGHY",
                              species == "virginicus" ~ "ELVI", species == "virginiana" ~ "ELVI")) %>% 
  dplyr::select(-contains("..."), -`Found?`) 

# Same records but downloaded from symbiota now!
AM_symbiota_records <- read_csv(file = "~/Dropbox/Josh&Tom - shared/Endo_Herbarium/DigitizedHerbariumRecords/TexasA&M_digitized_records/SymbOutput_2024-02-12_094323_DwC-A/occurrences.csv") %>% 
  mutate(primary_collector = recordedBy) %>% 
  mutate(primary_collector = case_when(primary_collector == "Michael M." ~ "Michael M. MacRoberts",
                                       primary_collector == "John C. Fr\x8emont" ~ "John C. Fremont",
                                      TRUE ~ primary_collector),
         primary_collector = gsub("\\--.*", "", primary_collector), # removing some notes in the collector info for two records
         collector_firstname = word(primary_collector),
         collector_lastname = case_when(str_detect(primary_collector,"Jr.") ~ word(str_trim(primary_collector), start = -2, end = -1),
                                        str_detect(primary_collector, "III") ~ word(str_trim(primary_collector),  start = -2, end = -1),
                                        !str_detect(primary_collector,"Jr.") | !str_detect(primary_collector, "III") ~ word(str_trim(primary_collector), -1))) %>% 
  mutate(Spp_code = case_when(specificEpithet == "hyemalis" ~ "AGHY", specificEpithet == "hiemalis" ~ "AGHY",
                              specificEpithet == "perennans" ~ "AGPE", specificEpithet == "schiedeana" ~ "AGPE",
                              specificEpithet == "virginicus" ~ "ELVI", specificEpithet == "striatus" ~ "ELVI")) %>% 
  mutate(Institution_specimen_id = word(otherCatalogNumbers,2, sep  = "; AccNo: "))

AM_merged_records <- AM_symbiota_records %>% 
  full_join(AM_records, by = c("Institution_specimen_id" = "Institution_specimen_id", 
                               "Spp_code" = "Spp_code", 
                               "year" = "year",
                               "month" = "month", 
                               "day" = "day",
                               "eventDate" = "DateCollected" ,
                               "primary_collector" = "primary_collector",
                               "collector_firstname" = "collector_firstname",
                               "collector_lastname" = "collector_lastname",
                               "recordNumber" = "collectorNumber")) 
  # group_by(Institution_specimen_id)  
  # coalesce()
  # dplyr::select(Institution_specimen_id, year, month, day, scientificName, primary_collector)


# BRIT digitized records downloaded from TORCH (includes Vanderbilt and U of Louisiana Monroe) 
# This was downloaded on Jul17 and we can get more specimens transcribed and download again.
AGHY_BRIT <- read_csv(file = "~/Dropbox/Josh&Tom - shared/Endo_Herbarium/DigitizedHerbariumRecords/BRIT_records/BRIT_AGHY_TORCH_records/SymbOutput_2024-04-03_142217_DwC-A/occurrences.csv",
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
  filter(!is.na(county), !is.na(eventDate)) %>% 
  filter(scientificName != "Agrostis perennans") %>% 
  separate(eventDate, into = c("year", "month", "day"), remove = FALSE) %>% 
  mutate_at(c("day", "month", "year"), as.numeric) %>% 
  mutate(decimalLongitude = case_when(decimalLongitude>0 ~ decimalLongitude*-1,
                                      TRUE ~ decimalLongitude)) %>% 
  mutate(recordedBy = case_when(recordedBy == "R. K Godfrey &" ~ "R. K Godfrey", TRUE ~ recordedBy)) %>% 
  dplyr::select(id, catalogNumber, country, stateProvince, county, municipality, locality, decimalLatitude, decimalLongitude, coordinateUncertaintyInMeters, recordedBy, recordNumber, eventDate, day, month, year) %>% 
  mutate(Spp_code = "AGHY")

ELVI_BRIT <- read_csv(file = "~/Dropbox/Josh&Tom - shared/Endo_Herbarium/DigitizedHerbariumRecords/BRIT_records/BRIT_ELVI_TORCH_records/SymbOutput_2024-04-03_141739_DwC-A/occurrences.csv",
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
  filter(!is.na(county), !is.na(eventDate)) %>% 
  separate(eventDate, into = c("year", "month", "day"), remove = FALSE) %>% 
  mutate_at(c("day", "month", "year"), as.numeric) %>% 
  dplyr::select(id, catalogNumber, country, stateProvince, county, municipality, locality, decimalLatitude, decimalLongitude, coordinateUncertaintyInMeters,recordedBy, recordNumber, eventDate, day, month, year) %>% 
  mutate(Spp_code = "ELVI")


AGPE_BRIT <- read_csv(file = "~/Dropbox/Josh&Tom - shared/Endo_Herbarium/DigitizedHerbariumRecords/BRIT_records/BRIT_AGPE_TORCH_records/SymbOutput_2024-04-03_142041_DwC-A/occurrences.csv",
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
                                       references = col_character()))  %>% 
  filter(!is.na(county), !is.na(eventDate)) %>% 
  separate(eventDate, into = c("year", "month", "day"), remove = FALSE) %>% 
  mutate_at(c("day", "month", "year"), as.numeric) %>% 
  dplyr::select(id, catalogNumber, country, stateProvince, county, municipality, locality, decimalLatitude, decimalLongitude, coordinateUncertaintyInMeters, recordedBy, recordNumber, eventDate, day, month, year) %>% 
  mutate(Spp_code = "AGPE")


BRIT_torch <- rbind(AGHY_BRIT, ELVI_BRIT, AGPE_BRIT) %>% 
  mutate(primary_collector = recordedBy,
         primary_collector = case_when(recordedBy == "Vernon Legrett and Debbie Legrett" ~ "Vernon Leggett, Debbie Leggett", TRUE ~ primary_collector),
         primary_collector = case_when(recordedBy == "John and Connie Taylor" | recordedBy == "John & Connie Taylor" ~ "John Taylor", TRUE ~ primary_collector),
         primary_collector = case_when(recordedBy == "S. & G. Jones" ~ "S. Jones", TRUE ~ primary_collector),
         primary_collector = case_when(recordedBy == "B.E. Dutton and David E. Taylor" ~ "B.E. Dutton", TRUE ~ primary_collector),
         primary_collector = case_when(recordedBy == "Dr. Sanders" | recordedBy == "Dana R. Sanders, Sr." ~ "Dana Sanders", TRUE ~ primary_collector),
         primary_collector = case_when(primary_collector == "Alfred Traveers" ~ "Alfred Traverse", TRUE ~ primary_collector),
         primary_collector = case_when(primary_collector == "A. D. McKeller" ~ "A. D. McKellar", TRUE ~ primary_collector),
         primary_collector = str_replace_all(primary_collector, "�", " "),
         primary_collector = str_replace_all(primary_collector, ", Jr.", ""), 
         primary_collector = str_replace_all(primary_collector, "Dr. ", "")) %>%
  mutate(collector_firstname = word(primary_collector),
         collector_lastname = case_when(str_detect(primary_collector,"Jr.") ~ word(str_trim(primary_collector), start = -2, end = -1),
                                        !str_detect(primary_collector,"Jr.") ~ word(str_trim(primary_collector), -1))) %>% 
  mutate(collector_firstname = case_when(primary_collector == "B.R " ~ "B.R & M.H",
                                         TRUE ~ collector_firstname),
         collector_lastname = case_when(primary_collector == "B.R " ~ "MacRoberts",
                                        TRUE ~ collector_lastname))

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
  mutate(primary_collector = case_when(str_detect(recordedBy, "Jr.") ~ word(recordedBy, start = 1, end = 2, sep = fixed(",")),
                                       !str_detect(recordedBy, "Jr.") ~ word(recordedBy, start = 1, sep = fixed(","))),
         primary_collector = case_when(recordedBy == "J. D. Schneidau (Jr.)" ~ "John D. Schneidau", 
                                       recordedBy == "John D. Schneidau (Jr.)" ~ "John D. Schneidau", TRUE ~ primary_collector),
         collector_firstname = word(primary_collector),
         collector_lastname = case_when(str_detect(primary_collector,"Jr.") ~ word(str_trim(primary_collector), start = -2, end = -1),
                                        !str_detect(primary_collector,"Jr.") ~ word(str_trim(primary_collector), -1))) %>% 
  mutate(Spp_code = case_when(genus == "Elymus" ~ "ELVI",
                              genus == "Agrostis" & specificEpithet == "hyemalis" | specificEpithet == "hiemalis" ~ "AGHY",
                              genus == "Agrostis" & specificEpithet == "perennans" ~ "AGPE")) %>%# there are some Agrostis scabra, that may need to be sorted out cause they could be part of hyemalis
  mutate(county = case_when(catalogNumber == "LSU00074817" ~ "Rensselaer", TRUE ~ county)) %>% 
  dplyr::select(id, Spp_code, catalogNumber, country, stateProvince, county, municipality, locality, decimalLatitude, decimalLongitude, coordinateUncertaintyInMeters, recordedBy, primary_collector, collector_lastname, collector_firstname, recordNumber, eventDate, day, month, year) 

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
  mutate(primary_collector = case_when(str_detect(recordedBy, "Jr") ~ word(recordedBy, start = 1, end = 2, sep = fixed(",")),
                                       !str_detect(recordedBy, "Jr") ~ word(recordedBy, start = 1, sep = fixed(","))),
         primary_collector = case_when(recordedBy == "Featherly & Cornelius" ~ "Henry I. Featherly", TRUE ~ primary_collector),
         primary_collector = word(primary_collector, sep = fixed("&")),
         collector_firstname = word(primary_collector),
         collector_lastname = case_when(str_detect(primary_collector,"Jr") ~ word(str_trim(primary_collector), start = -2, end = -1),
                                        !str_detect(primary_collector,"Jr") ~ word(str_trim(primary_collector), -1))) %>% 
  mutate(Spp_code = case_when(genus == "Elymus" ~ "ELVI",
                              genus == "Agrostis" & specificEpithet == "hyemalis" | specificEpithet == "hiemalis" ~ "AGHY",
                              genus == "Agrostis" & specificEpithet == "perennans" ~ "AGPE")) %>%    # there are some Agrostis scabra, that may need to be sorted out cause they could be part of hyemalis
  dplyr::select(id, catalogNumber, new_id, genus, specificEpithet, Spp_code, catalogNumber, country, stateProvince, county, municipality, locality, decimalLatitude, decimalLongitude, coordinateUncertaintyInMeters,   primary_collector, collector_lastname, collector_firstname, recordedBy, recordNumber, eventDate, day, month, year) 

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
  mutate(primary_collector = case_when(str_detect(recordedBy, "Jr") ~ word(recordedBy, start = 1, end = 2, sep = fixed(",")),
                                       !str_detect(recordedBy, "Jr") ~ word(recordedBy, start = 1, sep = fixed(","))),
         collector_firstname = word(primary_collector),
         collector_lastname = case_when(str_detect(primary_collector,"Jr") ~ word(str_trim(primary_collector), start = -2, end = -1),
                                        !str_detect(primary_collector,"Jr") ~ word(str_trim(primary_collector), -1))) %>% 
  mutate(Spp_code = case_when(genus == "Elymus" ~ "ELVI",
                              genus == "Agrostis" & specificEpithet == "hyemalis" | specificEpithet == "hiemalis" ~ "AGHY",
                              genus == "Agrostis" & specificEpithet == "perennans" ~ "AGPE")) %>%  # there are some Agrostis scabra, that may need to be sorted out cause they could be part of hyemalis
  dplyr::select(id, catalogNumber, otherCatalogNumbers, Spp_code, catalogNumber, country, stateProvince, county, municipality, locality, decimalLatitude, decimalLongitude, coordinateUncertaintyInMeters, primary_collector, collector_lastname, collector_firstname,  recordedBy, recordNumber, eventDate, day, month, year) #For the specimens that are digitized, they are found in the otherCatalogNumbers column. For OKL, I recorded a string of id's. We need to take the first, which is usually of the form "OKL######"


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
  mutate(primary_collector = case_when(str_detect(recordedBy, "Jr") ~ word(recordedBy, start = 1, end = 2, sep = fixed(",")),
                                       str_detect(recordedBy, "et al.") ~ word(recordedBy, start = 1, end = 2, sep = fixed(" ")),
                                       str_detect(recordedBy, "II") ~ word(recordedBy, start = 1, end = 2, sep = fixed(" ")),
                                       !str_detect(recordedBy, "Jr") & !str_detect(recordedBy, "II") ~ word(recordedBy, start = 1, sep = fixed(","))),
         primary_collector = word(primary_collector, sep = fixed(";")),
         primary_collector = word(primary_collector, sep = fixed("&")),
         collector_firstname = word(primary_collector),
         collector_lastname = case_when(str_detect(primary_collector,"Jr") ~ word(str_trim(primary_collector), start = -2, end = -1),
                                        str_detect(primary_collector, "II") ~ word(str_trim(primary_collector), start = -2, end = -1),
                                        !str_detect(primary_collector,"Jr") & !str_detect(primary_collector, "II") ~ word(str_trim(primary_collector), -1))) %>% 
  mutate(Spp_code = case_when(genus == "Elymus" ~ "ELVI",
                              genus == "Agrostis" & specificEpithet == "hyemalis" | specificEpithet == "hiemalis" ~ "AGHY",
                              genus == "Agrostis" & specificEpithet == "perennans" ~ "AGPE")) %>%   # there are some Agrostis scabra, that may need to be sorted out cause they could be part of hyemalis
  dplyr::select(id, new_id, catalogNumber, otherCatalogNumbers, Spp_code, catalogNumber, country, stateProvince, county, municipality, locality, decimalLatitude, decimalLongitude, coordinateUncertaintyInMeters, primary_collector, collector_lastname, collector_firstname,  recordedBy, recordNumber, eventDate, day, month, year) #For the specimens that are digitized, they are found in the otherCatalogNumbers column

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
                                           references = col_character()))  %>% 
  mutate(primary_collector = recordedBy,
         primary_collector = case_when(primary_collector == "Kuster, M." ~ "M. Kuster", 
                                       primary_collector == "Leis, S." ~ "S. Leis", 
                                       primary_collector == "Sutherland, C." ~ "C. Sutherland",
                                       primary_collector == "Serviss, Brett E." ~ "Brett E. Serviss", TRUE ~ primary_collector),
         collector_firstname = word(primary_collector),
         collector_lastname = case_when(str_detect(primary_collector,"Jr.") ~ word(str_trim(primary_collector), start = -2, end = -1),
                                        !str_detect(primary_collector,"Jr.") ~ word(str_trim(primary_collector), -1))) %>% 
  mutate(new_id = case_when(collector_lastname == "Shimek" ~ paste0(collector_lastname, year),
                            primary_collector == "Edith C. Bicknell" ~ paste0("EdithCBicknell", catalogNumber), 
                            TRUE ~ paste0(collector_lastname, recordNumber))) %>% 
  mutate(Spp_code = case_when(genus == "Elymus" ~ "ELVI",
                              genus == "Agrostis" & specificEpithet == "hyemalis" | specificEpithet == "hiemalis" ~ "AGHY",
                              genus == "Agrostis" & specificEpithet == "perennans" ~ "AGPE")) %>%  # there are some Agrostis scabra, that may need to be sorted out cause they could be part of hyemalis
  dplyr::select(id, new_id, recordedBy, primary_collector, collector_firstname, collector_lastname, recordID, catalogNumber, otherCatalogNumbers, Spp_code, catalogNumber, country, stateProvince, county, municipality, locality, decimalLatitude, decimalLongitude, coordinateUncertaintyInMeters, eventDate, day, month, year) #For the specimens that are digitized, they are found in the otherCatalogNumbers column

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
                          Spp_code == "AGPE" & is.na(year.x) ~ year.y, TRUE ~ year.x)) %>% 
  mutate(hand_georef_lon = case_when(as.numeric(decimalLongitude)>0 ~ (as.numeric(decimalLongitude))*-1,TRUE ~ as.numeric(decimalLongitude)),
         hand_georef_lat = as.numeric(decimalLatitude)) %>%
  mutate(primary_collector = case_when(str_detect(Collected_by, ", Jr") ~ word(Collected_by, start = 1, sep = fixed(",")),
                                       str_detect(Collected_by, " Jr") ~ word(Collected_by, start = 1, end = 2, sep = fixed(" ")),
                                       TRUE ~ Collected_by),
         primary_collector = word(primary_collector, sep = fixed(";")),
         primary_collector = word(primary_collector, sep = fixed("with")),
         primary_collector = word(primary_collector, sep = fixed("and")),
         primary_collector = word(primary_collector, sep = fixed("&")),
         primary_collector = str_replace(primary_collector,"\n", ""),
         primary_collector = str_replace(primary_collector, "Hoagl", "Hoagland"),
         primary_collector = str_replace(primary_collector, "R. Lence Weldon, Susan Harris, Elizabeth Harris, Tony L. Burgess", "R. Lence Weldon"),
         primary_collector = str_replace(primary_collector, "G. Davidse, D. Bell, ", "G. Davidse"),
         primary_collector = str_replace(primary_collector, "F.L. Johnson, R. Rudman,", "F.L. Johnson"),
         collector_firstname = word(primary_collector),
         collector_lastname = case_when(str_detect(primary_collector,"Jr") ~ word(str_trim(primary_collector), start = -2, end = -1),
                                        !str_detect(primary_collector,"Jr") ~ word(str_trim(primary_collector), -1))) %>% 
  dplyr::select(-year.x, - year.y)

# This is the sample info and we will filter for only those that we have scored so far.
sample_info <- read_csv(file = "~joshuacfowler/Dropbox/Josh&Tom - shared/Endo_Herbarium/Endo_Herbarium_sample.csv") %>%  
  dplyr::select(-contains("...")) %>% 
  rename(notes_uncertainty_1 = `notes (uncertainty, maturity of seeds, weird fungal morphology)`) %>% 
  mutate(across(everything(), as.character)) %>% 
  pivot_longer(cols = c(-Specimen_id, -Sample_id, -tissue_type)) %>% 
  mutate(score_number = str_sub(name, -1), name = sub("_[^_]+$", "", name)) %>% 
  filter(score_number == 1 | score_number == 2 & !is.na(value)) %>% 
  filter(!is.na(Specimen_id)) %>% 
  pivot_wider(id_cols = c(Specimen_id, Sample_id, tissue_type, score_number), names_from = name, values_from = value) %>%  
  mutate(across(c(seed_scored, seed_eplus, Endo_status_liberal, Endo_status_conservative), as.numeric)) %>% 
  mutate(Specimen_id_temp = Specimen_id) %>% 
  separate(Specimen_id_temp, c("Herbarium_id", "Spp_code", "Specimen_no"), "_")


testsample_info<- read_csv(file = "~joshuacfowler/Dropbox/Josh&Tom - shared/Endo_Herbarium/Endo_Herbarium_sample.csv") %>% 
  group_by(Specimen_id) %>% 
  summarize(count = n()) %>% 
  filter(count == 2)

endo_scores <- specimen_info  %>% 
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
                            is.na(county) ~ County, TRUE ~ county)) %>% 
  mutate(State = case_when(is.na(State) ~ state,
                           is.na(state) ~ State, TRUE ~ state)) %>% 
  mutate(Country = case_when(is.na(Country.x) ~ Country.y,
                             is.na(Country.y) ~ Country.x, TRUE ~ Country.y)) %>% 
  mutate(Municipality = case_when(is.na(Municipality) ~ Municpality,
                                  is.na(Municpality) ~ Municipality, TRUE ~ Municpality)) %>% 
  mutate(Locality = case_when(is.na(Locality_info) ~ Locality,
                              is.na(Locality) ~ Locality_info, TRUE ~ Locality)) %>% 
  mutate(year = case_when(is.na(year.x) ~ year.y,
                          is.na(year.y) ~ year.x, TRUE ~ year.y)) %>% 
  mutate(month = case_when(is.na(month.x) ~ month.y,
                           is.na(month.y) ~ month.x, TRUE ~ month.y)) %>% 
  mutate(day = case_when(is.na(day.x) ~ day.y,
                         is.na(day.y) ~ day.x, TRUE ~ day.y)) %>% 
  mutate(Spp_code = case_when(is.na(Spp_code.x) ~ Spp_code.y,
                              is.na(Spp_code.y) ~ Spp_code.x, TRUE ~ Spp_code.y)) %>% 
  mutate(primary_collector = case_when(is.na(primary_collector.x) ~ primary_collector.y,
                                       is.na(primary_collector.y) ~ primary_collector.x,  TRUE ~ primary_collector.y),
         collector_lastname = case_when(is.na(collector_lastname.x) ~ collector_lastname.y,
                                        is.na(collector_lastname.y) ~ collector_lastname.x,TRUE ~ collector_lastname.y),
         collector_firstname = case_when(is.na(collector_firstname.x) ~ collector_firstname.y,
                                        is.na(collector_firstname.y) ~ collector_firstname.x, TRUE ~ collector_firstname.y))  %>% 
  dplyr::select(Sample_id, Institution_specimen_id, Spp_code, new_id, primary_collector, collector_lastname, collector_firstname, Country, State, County, Municipality, Locality, hand_georef_lat, hand_georef_lon, year, month, day, tissue_type, seed_scored, seed_eplus, Endo_status_liberal, Endo_status_conservative, Date_scored, scorer_id, score_number, surface_area_cm2, mean_infl_length, mean_inflplusawn_length, infl_count)

# Merge in the BRIT records that we have so far
endo_herb2 <- endo_herb1 %>% 
  left_join(BRIT_torch, by = c("new_id" = "catalogNumber")) %>% 
  mutate(County = case_when(is.na(County) ~ county,
                            is.na(county) ~ County, TRUE ~ county)) %>% 
  mutate(State = case_when(is.na(State) ~ stateProvince,
                           is.na(stateProvince) ~ State, TRUE ~ stateProvince)) %>% 
  mutate(Country = case_when(is.na(Country) ~ country,
                             is.na(country) ~ Country, TRUE ~ country)) %>% 
  mutate(Municipality = case_when(is.na(Municipality) ~ municipality,
                                  is.na(municipality) ~ Municipality, TRUE ~ municipality)) %>% 
  mutate(Locality = case_when(is.na(Locality) ~ locality,
                              is.na(locality) ~ Locality, TRUE ~ locality)) %>% 
  mutate(hand_georef_lat = case_when(is.na(hand_georef_lat) ~ decimalLatitude,
                                     is.na(decimalLatitude) ~ hand_georef_lat,TRUE ~ decimalLatitude),
         hand_georef_lon = case_when(is.na(hand_georef_lon) ~ decimalLongitude,
                                     is.na(decimalLongitude) ~ hand_georef_lon,TRUE ~ decimalLatitude)) %>% 
  mutate(year = case_when(is.na(year.x) ~ year.y,
                          is.na(year.y) ~ year.x, TRUE ~ year.y)) %>% 
  mutate(month = case_when(is.na(month.x) ~ month.y,
                           is.na(month.y) ~ month.x, TRUE ~ day.y)) %>% 
  mutate(day = case_when(is.na(day.x) ~ day.y,
                         is.na(day.y) ~ day.x, TRUE ~ day.y)) %>% 
  mutate(Spp_code = case_when(is.na(Spp_code.x) ~ Spp_code.y,
                              is.na(Spp_code.y) ~ Spp_code.x, TRUE ~ Spp_code.y)) %>% 
  mutate(primary_collector = case_when(is.na(primary_collector.x) ~ primary_collector.y,
                                       is.na(primary_collector.y) ~ primary_collector.x, TRUE ~ primary_collector.y),
         collector_lastname = case_when(is.na(collector_lastname.x) ~ collector_lastname.y,
                                        is.na(collector_lastname.y) ~ collector_lastname.x, TRUE ~ collector_lastname.y),
         collector_firstname = case_when(is.na(collector_firstname.x) ~ collector_firstname.y,
                                         is.na(collector_firstname.y) ~ collector_firstname.x, TRUE ~ collector_firstname.y)) %>% 
  dplyr::select(Sample_id, Institution_specimen_id, Spp_code, new_id, primary_collector, collector_lastname, collector_firstname, Country, State, County, Municipality, Locality, hand_georef_lat, hand_georef_lon, year, month, day, tissue_type, seed_scored, seed_eplus, Endo_status_liberal, Endo_status_conservative, Date_scored, scorer_id, score_number, surface_area_cm2, mean_infl_length, mean_inflplusawn_length, infl_count)

# Merge in the UT Austin records that we have so far
endo_herb3 <- endo_herb2 %>% 
  left_join(UTAustin_torch, by = c("new_id" = "new_id")) %>% 
  mutate(County = case_when(is.na(County) ~ county,
                            is.na(county) ~ County, TRUE ~ county)) %>% 
  mutate(State = case_when(is.na(State) ~ stateProvince,
                           is.na(stateProvince) ~ State, TRUE ~ stateProvince)) %>% 
  mutate(Country = case_when(is.na(Country) ~ country,
                             is.na(country) ~ Country, TRUE ~ country)) %>% 
  mutate(Municipality = case_when(is.na(Municipality) ~ municipality,
                                  is.na(municipality) ~ Municipality, TRUE ~ municipality))  %>% 
  mutate(Locality = case_when(is.na(Locality) ~ locality,
                              is.na(locality) ~ Locality, TRUE ~ locality)) %>% 
  mutate(hand_georef_lat = case_when(is.na(hand_georef_lat) ~ decimalLatitude,
                                     is.na(decimalLatitude) ~ hand_georef_lat, TRUE ~ decimalLatitude),
         hand_georef_lon = case_when(is.na(hand_georef_lon) ~ decimalLongitude,
                                     is.na(decimalLongitude) ~ hand_georef_lon, TRUE ~ decimalLongitude)) %>% 
  mutate(year = case_when(is.na(year.x) ~ year.y,
                          is.na(year.y) ~ year.x, TRUE ~ year.y)) %>% 
  mutate(month = case_when(is.na(month.x) ~ month.y,
                           is.na(month.y) ~ month.x, TRUE ~ month.y)) %>% 
  mutate(day = case_when(is.na(day.x) ~ day.y,
                         is.na(day.y) ~ day.x, TRUE ~ day.y)) %>% 
  mutate(Spp_code = case_when(is.na(Spp_code.x) ~ Spp_code.y,
                              is.na(Spp_code.y) ~ Spp_code.x, TRUE ~ Spp_code.y)) %>% 
  mutate(primary_collector = case_when(is.na(primary_collector.x) ~ primary_collector.y,
                                       is.na(primary_collector.y) ~ primary_collector.x, TRUE ~ primary_collector.y),
         collector_lastname = case_when(is.na(collector_lastname.x) ~ collector_lastname.y,
                                        is.na(collector_lastname.y) ~ collector_lastname.x, TRUE ~ collector_lastname.y),
         collector_firstname = case_when(is.na(collector_firstname.x) ~ collector_firstname.y,
                                         is.na(collector_firstname.y) ~ collector_firstname.x, TRUE ~ collector_firstname.y)) %>% 
  dplyr::select(Sample_id, Institution_specimen_id, Spp_code, new_id, primary_collector, collector_lastname, collector_firstname, Country, State, County, Municipality, Locality, hand_georef_lat, hand_georef_lon, year, month, day, tissue_type, seed_scored, seed_eplus, Endo_status_liberal, Endo_status_conservative, Date_scored, scorer_id, score_number, surface_area_cm2, mean_infl_length, mean_inflplusawn_length, infl_count)


# Merge in the LSU records that we have so far
endo_herb4 <- endo_herb3 %>% 
  left_join(LSU_records, by = c("new_id" = "catalogNumber")) %>% 
  mutate(County = case_when(is.na(County) ~ county,
                            is.na(county) ~ County, TRUE ~ county)) %>% 
  mutate(State = case_when(is.na(State) ~ stateProvince,
                           is.na(stateProvince) ~ State, TRUE ~ stateProvince)) %>% 
  mutate(Country = case_when(is.na(Country) ~ country,
                             is.na(country) ~ Country, TRUE ~ country)) %>%  
  mutate(Municipality = case_when(is.na(Municipality) ~ municipality,
                                  is.na(municipality) ~ Municipality, TRUE ~ municipality))  %>% 
  mutate(Locality = case_when(is.na(Locality) ~ locality,
                              is.na(locality) ~ Locality, TRUE ~ locality)) %>% 
  mutate(hand_georef_lat = case_when(is.na(hand_georef_lat) ~ decimalLatitude,
                                     is.na(decimalLatitude) ~ hand_georef_lat, TRUE ~ decimalLatitude),
         hand_georef_lon = case_when(is.na(hand_georef_lon) ~ decimalLongitude,
                                     is.na(decimalLongitude) ~ hand_georef_lon, TRUE ~ decimalLongitude)) %>% 
  mutate(year = case_when(is.na(year.x) ~ year.y,
                          is.na(year.y) ~ year.x, TRUE ~ year.y)) %>% 
  mutate(month = case_when(is.na(month.x) ~ month.y,
                           is.na(month.y) ~ month.x, TRUE ~ month.y)) %>% 
  mutate(day = case_when(is.na(day.x) ~ day.y,
                         is.na(day.y) ~ day.x, TRUE ~ day.y)) %>% 
  mutate(Spp_code = case_when(is.na(Spp_code.x) ~ Spp_code.y,
                              is.na(Spp_code.y) ~ Spp_code.x, TRUE ~ Spp_code.y)) %>% 
  mutate(primary_collector = case_when(is.na(primary_collector.x) ~ primary_collector.y,
                                       is.na(primary_collector.y) ~ primary_collector.x, TRUE ~ primary_collector.y),
         collector_lastname = case_when(is.na(collector_lastname.x) ~ collector_lastname.y,
                                        is.na(collector_lastname.y) ~ collector_lastname.x, TRUE ~ collector_lastname.y),
         collector_firstname = case_when(is.na(collector_firstname.x) ~ collector_firstname.y,
                                         is.na(collector_firstname.y) ~ collector_firstname.x, TRUE ~ collector_firstname.y)) %>% 
  dplyr::select(Sample_id, Institution_specimen_id, Spp_code, new_id, primary_collector, collector_lastname, collector_firstname, Country, State, County, Municipality, Locality, hand_georef_lat, hand_georef_lon, year, month, day, tissue_type, seed_scored, seed_eplus, Endo_status_liberal, Endo_status_conservative, Date_scored, scorer_id, score_number, surface_area_cm2, mean_infl_length, mean_inflplusawn_length, infl_count)

# Merge in the OKL records that we have so far. This one is weird, and need to figure out which column to match on.
endo_herb5 <- endo_herb4 %>% 
  left_join(OKL_records, by = c("new_id" = "otherCatalogNumbers")) %>% 
  mutate(County = case_when(is.na(County) ~ county,
                            is.na(county) ~ County,
                            TRUE ~ county)) %>%
  mutate(State = case_when(is.na(State) ~ stateProvince,
                           is.na(stateProvince) ~ State,
                           TRUE ~ stateProvince)) %>%
  mutate(Country = case_when(is.na(Country) ~ country,
                             is.na(country) ~ Country,
                             TRUE ~ country)) %>%
  mutate(Municipality = case_when(is.na(Municipality) ~ municipality,
                                  is.na(municipality) ~ Municipality,
                                  TRUE ~ municipality))  %>%
  mutate(Locality = case_when(is.na(Locality) ~ locality,
                              is.na(locality) ~ Locality,
                              TRUE ~ locality)) %>%
  mutate(hand_georef_lat = case_when(is.na(hand_georef_lat) ~ decimalLatitude,
                                     is.na(decimalLatitude) ~ hand_georef_lat),
         hand_georef_lon = case_when(is.na(hand_georef_lon) ~ decimalLongitude,
                                     is.na(decimalLongitude) ~ hand_georef_lon)) %>% 
  mutate(year = case_when(is.na(year.x) ~ year.y,
                          is.na(year.y) ~ year.x,
                          TRUE ~ year.y)) %>%
  mutate(month = case_when(is.na(month.x) ~ month.y,
                           is.na(month.y) ~ month.x,
                           TRUE ~ month.y)) %>%
  mutate(day = case_when(is.na(day.x) ~ day.y,
                         is.na(day.y) ~ day.x,
                         TRUE ~ day.y)) %>%
  mutate(Spp_code = case_when(is.na(Spp_code.x) ~ Spp_code.y,
                              is.na(Spp_code.y) ~ Spp_code.x)) %>%
  mutate(primary_collector = case_when(is.na(primary_collector.x) ~ primary_collector.y,
                                       is.na(primary_collector.y) ~ primary_collector.x,
                                       TRUE ~ primary_collector.y),
         collector_lastname = case_when(is.na(collector_lastname.x) ~ collector_lastname.y,
                                        is.na(collector_lastname.y) ~ collector_lastname.x,
                                        TRUE ~ collector_lastname.y),
         collector_firstname = case_when(is.na(collector_firstname.x) ~ collector_firstname.y,
                                         is.na(collector_firstname.y) ~ collector_firstname.x,
                                         TRUE ~ collector_firstname.y)) %>% 
  dplyr::select(Sample_id, Institution_specimen_id, Spp_code, new_id, primary_collector, collector_lastname, collector_firstname, Country, State, County, Municipality, Locality, hand_georef_lat, hand_georef_lon, year, month, day, tissue_type, seed_scored, seed_eplus, Endo_status_liberal, Endo_status_conservative, Date_scored, scorer_id, score_number, surface_area_cm2, mean_infl_length, mean_inflplusawn_length, infl_count)

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
  mutate(hand_georef_lat = case_when(is.na(hand_georef_lat) ~ decimalLatitude,
                                     is.na(decimalLatitude) ~ hand_georef_lat),
         hand_georef_lon = case_when(is.na(hand_georef_lon) ~ decimalLongitude,
                                     is.na(decimalLongitude) ~ hand_georef_lon)) %>% 
  mutate(year = case_when(is.na(year.x) ~ year.y,
                          is.na(year.y) ~ year.x)) %>% 
  mutate(month = case_when(is.na(month.x) ~ month.y,
                           is.na(month.y) ~ month.x)) %>% 
  mutate(day = case_when(is.na(day.x) ~ day.y,
                         is.na(day.y) ~ day.x)) %>% 
  mutate(Spp_code = case_when(is.na(Spp_code.x) ~ Spp_code.y,
                              is.na(Spp_code.y) ~ Spp_code.x)) %>% 
  mutate(primary_collector = case_when(is.na(primary_collector.x) ~ primary_collector.y,
                                       is.na(primary_collector.y) ~ primary_collector.x),
         collector_lastname = case_when(is.na(collector_lastname.x) ~ collector_lastname.y,
                                        is.na(collector_lastname.y) ~ collector_lastname.x),
         collector_firstname = case_when(is.na(collector_firstname.x) ~ collector_firstname.y,
                                         is.na(collector_firstname.y) ~ collector_firstname.x)) %>% 
  dplyr::select(Sample_id, Institution_specimen_id, Spp_code, new_id, primary_collector, collector_lastname, collector_firstname, Country, State, County, Municipality, Locality, hand_georef_lat, hand_georef_lon, year, month, day, tissue_type, seed_scored, seed_eplus, Endo_status_liberal, Endo_status_conservative, Date_scored, scorer_id, score_number, surface_area_cm2, mean_infl_length, mean_inflplusawn_length, infl_count)


# Merge in the KANU records that we have so far
endo_herb7 <- endo_herb6 %>% 
  left_join(KANU_records, by = c("new_id" = "new_id")) %>% 
  mutate(County = case_when(is.na(County) ~ county,
                            is.na(county) ~ County,
                            TRUE ~ County)) %>% 
  mutate(State = case_when(is.na(State) ~ stateProvince,
                           is.na(stateProvince) ~ State,
                           TRUE ~ State)) %>% 
  mutate(Country = case_when(is.na(Country) ~ country,
                             is.na(country) ~ Country,
                             TRUE ~ Country)) %>%  
  mutate(Municipality = case_when(is.na(Municipality) ~ municipality,
                                  is.na(municipality) ~ Municipality,
                                  TRUE ~ Municipality))  %>% 
  mutate(Locality = case_when(is.na(Locality) ~ locality,
                              is.na(locality) ~ Locality,
                              TRUE ~ Locality)) %>% 
  mutate(hand_georef_lat = case_when(is.na(hand_georef_lat) ~ decimalLatitude,
                                     is.na(decimalLatitude) ~ hand_georef_lat),
         hand_georef_lon = case_when(is.na(hand_georef_lon) ~ decimalLongitude,
                                     is.na(decimalLongitude) ~ hand_georef_lon)) %>% 
  mutate(year = case_when(is.na(year.x) ~ year.y,
                          is.na(year.y) ~ year.x,
                          TRUE ~ year.x)) %>% 
  mutate(month = case_when(is.na(month.x) ~ month.y,
                           is.na(month.y) ~ month.x,
                           TRUE ~ month.x)) %>% 
  mutate(day = case_when(is.na(day.x) ~ day.y,
                         is.na(day.y) ~ day.x,
                         TRUE ~ day.x)) %>% 
  mutate(Spp_code = case_when(is.na(Spp_code.x) ~ Spp_code.y,
                              is.na(Spp_code.y) ~ Spp_code.x,
                              TRUE ~ Spp_code.x)) %>% 
  mutate(primary_collector = case_when(is.na(primary_collector.x) ~ primary_collector.y,
                                       is.na(primary_collector.y) ~ primary_collector.x,
                                       TRUE ~ primary_collector.x),
         collector_lastname = case_when(is.na(collector_lastname.x) ~ collector_lastname.y,
                                        is.na(collector_lastname.y) ~ collector_lastname.x,
                                        TRUE ~ primary_collector.x),
         collector_firstname = case_when(is.na(collector_firstname.x) ~ collector_firstname.y,
                                         is.na(collector_firstname.y) ~ collector_firstname.x,
                                         TRUE ~ primary_collector.x)) %>% 
  dplyr::select(Sample_id, Institution_specimen_id, Spp_code, new_id, primary_collector, collector_lastname, collector_firstname, Country, State, County, Municipality, Locality, hand_georef_lat, hand_georef_lon, year, month, day, tissue_type, seed_scored, seed_eplus, Endo_status_liberal, Endo_status_conservative, Date_scored, scorer_id, score_number, surface_area_cm2, mean_infl_length, mean_inflplusawn_length, infl_count)



# Merge in the MOBOT records that we have so far
ci_str_detect <- function(x, y){str_detect(x, regex(y, ignore_case = TRUE))}
endo_herb8 <- endo_herb7 %>% 
  fuzzy_left_join(MOBOT_records, match_fun = ci_str_detect, by = c("new_id" = "new_id")) %>% 
  mutate(County = case_when(is.na(County) ~ county,
                            is.na(county) ~ County, TRUE ~ County)) %>% 
  mutate(State = case_when(is.na(State) ~ stateProvince,
                           is.na(stateProvince) ~ State, TRUE ~ State)) %>% 
  mutate(Country = case_when(is.na(Country) ~ country,
                             is.na(country) ~ Country, TRUE ~ Country)) %>%  
  mutate(Municipality = case_when(is.na(Municipality) ~ municipality,
                                  is.na(municipality) ~ Municipality, TRUE ~ Municipality))  %>% 
  mutate(Locality = case_when(is.na(Locality) ~ locality,
                              is.na(locality) ~ Locality, TRUE ~ Locality)) %>% 
  mutate(hand_georef_lat = case_when(is.na(hand_georef_lat) ~ decimalLatitude,
                                     is.na(decimalLatitude) ~ hand_georef_lat),
         hand_georef_lon = case_when(is.na(hand_georef_lon) ~ decimalLongitude,
                                     is.na(decimalLongitude) ~ hand_georef_lon)) %>% 
  mutate(year = case_when(is.na(year.x) ~ year.y,
                          is.na(year.y) ~ year.x, TRUE ~ year.y)) %>% 
  mutate(month = case_when(is.na(month.x) ~ month.y,
                           is.na(month.y) ~ month.x, TRUE ~ month.y)) %>% 
  mutate(day = case_when(is.na(day.x) ~ day.y,
                         is.na(day.y) ~ day.x, TRUE ~ day.y)) %>% 
  mutate(Spp_code = case_when(is.na(Spp_code.x) ~ Spp_code.y,
                              is.na(Spp_code.y) ~ Spp_code.x)) %>% 
  mutate(primary_collector = case_when(is.na(primary_collector.x) ~ primary_collector.y,
                                       is.na(primary_collector.y) ~ primary_collector.x, TRUE ~ primary_collector.y),
         collector_lastname = case_when(is.na(collector_lastname.x) ~ collector_lastname.y,
                                        is.na(collector_lastname.y) ~ collector_lastname.x, TRUE ~ collector_lastname.y),
         collector_firstname = case_when(is.na(collector_firstname.x) ~ collector_firstname.y,
                                         is.na(collector_firstname.y) ~ collector_firstname.x, TRUE ~ collector_firstname.y)) %>% 
  mutate(new_id = new_id.x) %>%
  mutate(primary_collector = case_when(primary_collector == "Hein, S." ~ "S. Hein", TRUE ~ primary_collector),
         collector_lastname = case_when(primary_collector == "S. Hein" ~ "Hein", TRUE ~ collector_lastname),
         collector_firstname = case_when(primary_collector == "S. Hein" ~ "S.", TRUE ~ collector_firstname),
         collector_lastname = str_replace_all(collector_lastname, "�", " "),
         collector_lastname = str_replace_all(collector_lastname, "xa0", " "),
         collector_lastname = str_replace_all(collector_lastname, ", Jr.", " "),
         collector_lastname = str_replace_all(collector_lastname, ", Jr", " "),
         collector_lastname = str_replace_all(collector_lastname, " Jr.", " ")) %>% 
  mutate(collector_full_string = paste(collector_firstname, collector_lastname),
         collector_first_initial  = str_sub(collector_firstname, 1,1),
         collector_string = paste(collector_first_initial, collector_lastname)) %>% 
  mutate(collector_string = case_when(collector_full_string == "Robert Jones" ~ "Robert Jones",
                                      collector_full_string == "Ronald Jones" ~ "Ronald Jones",
                                      collector_full_string == "Latimore Smith" ~ "Latimore Smith",
                                      collector_full_string == "Logan Smith" ~ "Logan Smith",
                                      collector_full_string == "NA NA" ~ NA,
                                      TRUE ~ collector_string)) %>% 
  mutate(scorer_id = case_when(scorer_id == "NA" ~ NA,
                               TRUE ~ scorer_id)) %>% 
  dplyr::select(Sample_id, Institution_specimen_id, Spp_code, new_id, primary_collector, collector_lastname, collector_firstname, collector_full_string, collector_string, Country, State, County, Municipality, Locality, hand_georef_lat, hand_georef_lon, year, month, day, tissue_type, seed_scored, seed_eplus, Endo_status_liberal, Endo_status_conservative, Date_scored, scorer_id, score_number, surface_area_cm2, mean_infl_length, mean_inflplusawn_length, infl_count)


endo_herb_checking <- endo_herb8 %>% 
  filter(tissue_type == "seed") %>% 
  filter(!is.na(Endo_status_liberal) & is.na(year)) %>% 
  dplyr::select(Sample_id, primary_collector, collector_lastname, collector_firstname, County, year, month, day, County, seed_scored, score_number)

specimen_counts <- endo_herb8 %>% 
  mutate(Spp_code = case_when(grepl("AGHY", Sample_id) ~ "AGHY",
                              grepl("ELVI", Sample_id) ~ "ELVI",
                              grepl("AGPE", Sample_id) ~ "AGPE")) %>% 
  group_by(Spp_code, tissue_type) %>% 
  summarize(n())




 # checking in on the accuracy of splitting the collector names

collector_string_count <- endo_herb8 %>% 
  group_by(collector_string) %>% 
  summarize(n()) 
 
# ggplot(collector_string_count)+
#    geom_bar(aes(y = `n()`, x = collector_string), stat = "identity")

length(unique(collector_string_count$collector_string))


# Creating scorer and collector levels
scorer_levels <- levels(as.factor(endo_herb8$scorer_id))
scorer_no <- paste0("Scorer",1:nlevels(as.factor(endo_herb8$scorer_id)))

endo_herb8$scorer_factor <- scorer_no[match(as.factor(endo_herb8$scorer_id), scorer_levels)]

# updating the collector levels

collector_levels <- levels(as.factor(endo_herb8$collector_string))
collector_no <- paste0("Collector",1:nlevels(as.factor(endo_herb8$collector_string)))

endo_herb8$collector_factor <- collector_no[match(as.factor(endo_herb8$collector_string), collector_levels)]

collector_checking <- endo_herb8 %>% 
  filter(is.na(collector_string) & !is.na(Endo_status_liberal) & !is.na(County))



location_checking <- endo_herb8 %>% 
  filter(!is.na(Endo_status_liberal) & is.na(Country))


# ggplot(unique_lastnames)+
#   geom_histogram(aes(x = no_records))




# Now I am going to link these county/locality records to a gps point with ggmap
# This requires and API key which you can set up through google, look at ?register_google.
# There are restrictions to the total number of queries that you can do per day and month, and if you go over, it costs money, so we will save the output. I believe we have a free trial for year.
# One other note, is that we are only using the level of county/city/state which I believe should be pretty accurate through google. I'm not sure that it could accurately do a more detailed locality string
# I have found software that could do this, but would require a human to double check the output.
#
register_google()
endo_herb_georef <-endo_herb8 %>%
  mutate(county_label = case_when(State == "Louisiana" | State == "LA" ~ " Parish",
                                  is.na(County) ~ NA,
                                  County == "" ~ NA,
                                  Country == "Canada" ~ NA,
                                  Country == "Mexico" ~ NA, 
                                  TRUE ~ " County")) %>% 
  unite("County_fixed", sep = "", County, county_label, remove = FALSE, na.rm = TRUE) %>% 
  unite("location_string" , sep = ", " , Municipality,County_fixed,State,Country, remove = FALSE, na.rm = TRUE) %>%
  # mutate(location_string = case_when(location_string == "NA NA" ~ NA, TRUE ~ location_string)) %>% 
  filter(Endo_status_liberal <= 1) %>% 
  mutate_geocode(location_string, output = "more") # Uncomment this to run the geocoding.

endo_herb_georef_1 <- endo_herb_georef %>% 
  mutate(sample_temp = Sample_id) %>% 
  separate(sample_temp, into = c("Herb_code", "spp_code", "specimen_code", "tissue_code")) %>% 
  filter(scorer_factor != "Scorer26")

write_csv(endo_herb_georef_1, file = "~/Dropbox/Josh&Tom - shared/Endo_Herbarium/DigitizedHerbariumRecords/endo_herb_georef.csv") # saving a version of the file used in ensuing scripts

endo_herb_georef_dryad <- endo_herb_georef_1 %>% 
  mutate(Spp_code = spp_code) %>% 
  dplyr::select(c(Sample_id, Institution_specimen_id, Spp_code, new_id, primary_collector, collector_lastname, collector_firstname, collector_full_string, collector_string, location_string, 
                  Country, State, County, Municipality, Locality, hand_georef_lat, hand_georef_lon, year, month, day, tissue_type, seed_scored, seed_eplus, 
                  Endo_status_liberal, Endo_status_conservative, Date_scored, scorer_id, score_number, scorer_factor, collector_factor, lat, lon))
write_csv(endo_herb_georef_dryad, file = "~/Desktop/endo_herb_georef.csv") # saving a version of the file used in ensuing scripts







endo_herb_georef <- read_csv(file = "~/Dropbox/Josh&Tom - shared/Endo_Herbarium/DigitizedHerbariumRecords/endo_herb_georef.csv") %>%
  # filter(Country != "Canada") %>%
  mutate(Spp_code = case_when(grepl("AGHY", Sample_id) ~ "AGHY",
                              grepl("ELVI", Sample_id) ~ "ELVI",
                              grepl("AGPE", Sample_id) ~ "AGPE")) 


           
specimen_counts <- endo_herb_georef %>% 
  mutate(Spp_code = case_when(grepl("AGHY", Sample_id) ~ "AGHY",
                              grepl("ELVI", Sample_id) ~ "ELVI",
                              grepl("AGPE", Sample_id) ~ "AGPE")) %>% 
  group_by(Spp_code, tissue_type) %>% 
  summarize(n())


scored_counts <- endo_herb_georef %>% 
  filter(!is.na(Endo_status_liberal), !is.na(lon), score_number == 1) %>% 
  mutate(Spp_code = case_when(grepl("AGHY", Sample_id) ~ "AGHY",
                              grepl("ELVI", Sample_id) ~ "ELVI",
                              grepl("AGPE", Sample_id) ~ "AGPE")) %>% 
  group_by(Spp_code) %>% 
  summarize(n())

scored_counts <- endo_herb_georef %>% 
  filter(!is.na(Endo_status_liberal), score_number == 1) %>% 
  mutate(Spp_code = case_when(grepl("AGHY", Sample_id) ~ "AGHY",
                              grepl("ELVI", Sample_id) ~ "ELVI",
                              grepl("AGPE", Sample_id) ~ "AGPE")) %>% 
  group_by(Spp_code, tissue_type) %>% 
  summarize(n())

hist(endo_herb_georef$year)
plot(endo_herb_georef$lon, endo_herb_georef$lat)
plot(endo_herb_georef$lon, endo_herb_georef$Endo_status_liberal)
plot(endo_herb_georef$lat, endo_herb_georef$Endo_status_liberal)
plot(filter(endo_herb_georef, year>200)$year, filter(endo_herb_georef, year>200)$Endo_status_liberal)


endo_herb_AGHY <- endo_herb_georef %>% 
  filter(!is.na(Endo_status_liberal)) %>% 
  filter(Spp_code == "AGHY") %>% 
  filter(!is.na(lon) & !is.na(year)) 
plot(endo_herb_AGHY$lon, endo_herb_AGHY$lat)

endo_herb_ELVI <- endo_herb_georef %>% 
  filter(!is.na(Endo_status_liberal)) %>% 
  filter(Spp_code == "ELVI") %>% 
  filter(!is.na(lon) & !is.na(year))
plot(endo_herb_ELVI$lon, endo_herb_ELVI$lat)

endo_herb_AGPE <- endo_herb_georef %>% 
  filter(!is.na(Endo_status_liberal)) %>% 
  filter(Spp_code == "AGPE") %>% 
  filter(!is.na(lon) & !is.na(year))
plot(endo_herb_AGPE$lon, endo_herb_AGPE$lat)

### Looking at how many samples are scored without locations and in need of transcription
#There is a lot of overlap in the needs collector info and need latitude info

# need to fix Fontana Kansas lat lon
need_collector <- endo_herb_georef %>% 
  filter(Spp_code %in% c("AGHY", "AGPE", "ELVI" )) %>% 
  filter(!is.na(Endo_status_liberal)) %>% 
  filter(is.na(collector_string)) %>% 
  dplyr::select(Sample_id, Institution_specimen_id, Spp_code, new_id, primary_collector, collector_string, location_string, Country, year, lon, lat)

write_csv(need_collector, file = "~/Downloads/endo_herbarium_need_collector.csv")


need_latlon <- endo_herb_georef %>% 
  filter(Spp_code %in% c("AGHY", "AGPE", "ELVI" )) %>% 
  filter(!is.na(Endo_status_liberal), is.na(lon)) %>% 
  dplyr::select(Sample_id, Institution_specimen_id, Spp_code, new_id, primary_collector, collector_string, location_string, Country, year, lon, lat) %>% 
  left_join(specimen_info %>% select(Specimen_id, photo_transcribing_notes), by = c("Sample_id" = "Specimen_id"))


write_csv(need_latlon, file = "~/Downloads/endo_herbarium_need_latlon.csv")

# How many from each collection
need_latlon_summary <- need_latlon %>% 
  mutate(collection = word(Sample_id, 1, sep = "_")) %>% 
  group_by(collection) %>% 
  summarize(count = n())

# summary of how many scored specimens from each collection
collections_summary <- endo_herb_georef %>% 
  mutate(collection = word(Sample_id, 1, sep = "_")) %>% 
  group_by(collection, Spp_code) %>% 
  summarize(count = n())



endo_herb_georef %>% filter(!is.na(Endo_status_liberal) & is.na(year)) %>% view()
endo_herb_georef %>% filter(!is.na(Endo_status_liberal) & !is.na(County) & is.na(lat)) %>% view()

####################################################################################################
######### Connecting the herbarium records to climate data from PRISM ##############################
####################################################################################################
# making a folder to store prism data
prism_set_dl_dir(paste0(getwd(),"/prism_download"))

# getting monthly data for mean temp and precipitation
# takes a long time the first time, but can skip when you have raster files saved on your computer.
get_prism_monthlys(type = "tmean", years = 1895:2020, mon = 1:12, keepZip = FALSE)
get_prism_monthlys(type = "ppt", years = 1895:2020, mon = 1:12, keepZip = FALSE)


# pulling out values to get normals for old and new time periods
tmean_annual_recent_norm <- terra::mean(terra::rast(pd_stack(prism_archive_subset(type = "tmean", temp_period = "monthly", year = 1990:2020))))
tmean_spring_recent_norm <- terra::mean(terra::rast(pd_stack(prism_archive_subset(type = "tmean", temp_period = "monthly", year = 1990:2020, mon = 1:4))))
tmean_summer_recent_norm <- terra::mean(terra::rast(pd_stack(prism_archive_subset(type = "tmean", temp_period = "monthly", year = 1990:2020, mon = 5:8))))
tmean_autumn_recent_norm <- terra::mean(terra::rast(pd_stack(prism_archive_subset(type = "tmean", temp_period = "monthly", year = 1990:2020, mon = 9:12))))


tmean_annual_old_norm <- terra::mean(terra::rast(pd_stack(prism_archive_subset(type = "tmean", temp_period = "monthly", year = 1895:1925))))
tmean_spring_old_norm <- terra::mean(terra::rast(pd_stack(prism_archive_subset(type = "tmean", temp_period = "monthly", year = 1895:1925, mon = 1:4))))
tmean_summer_old_norm <- terra::mean(terra::rast(pd_stack(prism_archive_subset(type = "tmean", temp_period = "monthly", year = 1895:1925, mon = 5:8))))
tmean_autumn_old_norm <- terra::mean(terra::rast(pd_stack(prism_archive_subset(type = "tmean", temp_period = "monthly", year = 1895:1925, mon = 9:12))))

# calculating standard deviation in temp

tmean_annual_recent_sd <- terra::stdev(terra::rast(pd_stack(prism_archive_subset(type = "tmean", temp_period = "monthly", year = 1990:2020))))
tmean_spring_recent_sd <- terra::stdev(terra::rast(pd_stack(prism_archive_subset(type = "tmean", temp_period = "monthly", year = 1990:2020, mon = 1:4))))
tmean_summer_recent_sd <- terra::stdev(terra::rast(pd_stack(prism_archive_subset(type = "tmean", temp_period = "monthly", year = 1990:2020, mon = 5:8))))
tmean_autumn_recent_sd <- terra::stdev(terra::rast(pd_stack(prism_archive_subset(type = "tmean", temp_period = "monthly", year = 1990:2020, mon = 9:12))))

tmean_annual_old_sd <- terra::stdev(terra::rast(pd_stack(prism_archive_subset(type = "tmean", temp_period = "monthly", year = 1895:1925))))
tmean_spring_old_sd <- terra::stdev(terra::rast(pd_stack(prism_archive_subset(type = "tmean", temp_period = "monthly", year = 1895:1925, mon = 1:4))))
tmean_summer_old_sd <- terra::stdev(terra::rast(pd_stack(prism_archive_subset(type = "tmean", temp_period = "monthly", year = 1895:1925, mon = 5:8))))
tmean_autumn_old_sd <- terra::stdev(terra::rast(pd_stack(prism_archive_subset(type = "tmean", temp_period = "monthly", year = 1895:1925, mon = 9:12))))

# calculating coefficient of variation in temp
tmean_annual_recent_cv <- tmean_annual_recent_sd/tmean_annual_recent_norm
tmean_spring_recent_cv <- tmean_spring_recent_sd/tmean_spring_recent_norm
tmean_summer_recent_cv <- tmean_summer_recent_sd/tmean_summer_recent_norm
tmean_autumn_recent_cv <- tmean_autumn_recent_sd/tmean_autumn_recent_norm

tmean_annual_old_cv <- tmean_annual_old_sd/tmean_annual_old_norm
tmean_spring_old_cv <- tmean_spring_old_sd/tmean_spring_old_norm
tmean_summer_old_cv <- tmean_summer_old_sd/tmean_summer_old_norm
tmean_autumn_old_cv <- tmean_autumn_old_sd/tmean_autumn_old_norm



# calculating the cumulative precipitation for each year and for each season within the year
ppt_annual_recent <- ppt_spring_recent <- ppt_summer_recent <- ppt_autumn_recent <- ppt_winter_recent<- list()
for(y in 1990:2020){
ppt_annual_recent[[y]] <- sum(terra::rast(pd_stack(prism_archive_subset(type = "ppt", temp_period = "monthly", year = y))))
ppt_spring_recent[[y]] <- sum(terra::rast(pd_stack(prism_archive_subset(type = "ppt", temp_period = "monthly", year = y, mon = 1:4))))
ppt_summer_recent[[y]] <- sum(terra::rast(pd_stack(prism_archive_subset(type = "ppt", temp_period = "monthly", year = y, mon = 5:8))))
ppt_autumn_recent[[y]] <- sum(terra::rast(pd_stack(prism_archive_subset(type = "ppt", temp_period = "monthly", year = y, mon = 9:12)))) 
}

ppt_annual_old <- ppt_spring_old <- ppt_summer_old <- ppt_autumn_old <- ppt_winter_old<- list()
for(y in 1895:1925){
  ppt_annual_old[[y]] <- sum(terra::rast(pd_stack(prism_archive_subset(type = "ppt", temp_period = "monthly", year = y))))
  ppt_spring_old[[y]] <- sum(terra::rast(pd_stack(prism_archive_subset(type = "ppt", temp_period = "monthly", year = y, mon = 1:4))))
  ppt_summer_old[[y]] <- sum(terra::rast(pd_stack(prism_archive_subset(type = "ppt", temp_period = "monthly", year = y, mon = 5:8))))
  ppt_autumn_old[[y]] <- sum(terra::rast(pd_stack(prism_archive_subset(type = "ppt", temp_period = "monthly", year = y, mon = 9:12)))) # including December here because the year needs to wrap around...
}


# Taking the mean of the cumulative precipation values
ppt_annual_recent_norm <- terra::mean(terra::rast(unlist(ppt_annual_recent)))
ppt_spring_recent_norm <- terra::mean(terra::rast(unlist(ppt_spring_recent)))
ppt_summer_recent_norm <- terra::mean(terra::rast(unlist(ppt_summer_recent)))
ppt_autumn_recent_norm <- terra::mean(terra::rast(unlist(ppt_autumn_recent)))

ppt_annual_old_norm <- terra::mean(terra::rast(unlist(ppt_annual_old)))
ppt_spring_old_norm <- terra::mean(terra::rast(unlist(ppt_spring_old)))
ppt_summer_old_norm <- terra::mean(terra::rast(unlist(ppt_summer_old)))
ppt_autumn_old_norm <- terra::mean(terra::rast(unlist(ppt_autumn_old)))

#calculating the standard devation in precip
ppt_annual_recent_sd <- terra::stdev(terra::rast(unlist(ppt_annual_recent)))
ppt_spring_recent_sd <- terra::stdev(terra::rast(unlist(ppt_spring_recent)))
ppt_summer_recent_sd <- terra::stdev(terra::rast(unlist(ppt_summer_recent)))
ppt_autumn_recent_sd <- terra::stdev(terra::rast(unlist(ppt_autumn_recent)))

ppt_annual_old_sd <- terra::stdev(terra::rast(unlist(ppt_annual_old)))
ppt_spring_old_sd <- terra::stdev(terra::rast(unlist(ppt_spring_old)))
ppt_summer_old_sd <- terra::stdev(terra::rast(unlist(ppt_summer_old)))
ppt_autumn_old_sd <- terra::stdev(terra::rast(unlist(ppt_autumn_old)))

# calculating the coefficient of variation in precip
ppt_annual_recent_cv <- ppt_annual_recent_sd/ppt_annual_recent_norm
ppt_spring_recent_cv <- ppt_spring_recent_sd/ppt_spring_recent_norm
ppt_summer_recent_cv <- ppt_summer_recent_sd/ppt_summer_recent_norm
ppt_autumn_recent_cv <- ppt_autumn_recent_sd/ppt_autumn_recent_norm

ppt_annual_old_cv <- ppt_annual_old_sd/ppt_annual_old_norm
ppt_spring_old_cv <- ppt_spring_old_sd/ppt_spring_old_norm
ppt_summer_old_cv <- ppt_summer_old_sd/ppt_summer_old_norm
ppt_autumn_old_cv <- ppt_autumn_old_sd/ppt_autumn_old_norm





# calculating the difference over time
tmean_annual_difference <- terra::diff(terra::rast(list(tmean_annual_old_norm, tmean_annual_recent_norm)))
tmean_spring_difference <- terra::diff(terra::rast(list(tmean_spring_old_norm, tmean_spring_recent_norm)))
tmean_summer_difference <- terra::diff(terra::rast(list(tmean_summer_old_norm, tmean_summer_recent_norm)))
tmean_autumn_difference <- terra::diff(terra::rast(list(tmean_autumn_old_norm, tmean_autumn_recent_norm)))

tmean_annual_sd_difference <- terra::diff(terra::rast(list(tmean_annual_old_sd, tmean_annual_recent_sd)))
tmean_spring_sd_difference <- terra::diff(terra::rast(list(tmean_spring_old_sd, tmean_spring_recent_sd)))
tmean_summer_sd_difference <- terra::diff(terra::rast(list(tmean_summer_old_sd, tmean_summer_recent_sd)))
tmean_autumn_sd_difference <- terra::diff(terra::rast(list(tmean_autumn_old_sd, tmean_autumn_recent_sd)))

tmean_annual_cv_difference <- terra::diff(terra::rast(list(tmean_annual_old_cv, tmean_annual_recent_cv)))
tmean_spring_cv_difference <- terra::diff(terra::rast(list(tmean_spring_old_cv, tmean_spring_recent_cv)))
tmean_summer_cv_difference <- terra::diff(terra::rast(list(tmean_summer_old_cv, tmean_summer_recent_cv)))
tmean_autumn_cv_difference <- terra::diff(terra::rast(list(tmean_autumn_old_cv, tmean_autumn_recent_cv)))


ppt_annual_difference <- terra::diff(terra::rast(list(ppt_annual_old_norm, ppt_annual_recent_norm)))
ppt_spring_difference <- terra::diff(terra::rast(list(ppt_spring_old_norm, ppt_spring_recent_norm)))
ppt_summer_difference <- terra::diff(terra::rast(list(ppt_summer_old_norm, ppt_summer_recent_norm)))
ppt_autumn_difference <- terra::diff(terra::rast(list(ppt_autumn_old_norm, ppt_autumn_recent_norm)))

ppt_annual_sd_difference <- terra::diff(terra::rast(list(ppt_annual_old_sd, ppt_annual_recent_sd)))
ppt_spring_sd_difference <- terra::diff(terra::rast(list(ppt_spring_old_sd, ppt_spring_recent_sd)))
ppt_summer_sd_difference <- terra::diff(terra::rast(list(ppt_summer_old_sd, ppt_summer_recent_sd)))
ppt_autumn_sd_difference <- terra::diff(terra::rast(list(ppt_autumn_old_sd, ppt_autumn_recent_sd)))

ppt_annual_cv_difference <- terra::diff(terra::rast(list(ppt_annual_old_cv, ppt_annual_recent_cv)))
ppt_spring_cv_difference <- terra::diff(terra::rast(list(ppt_spring_old_cv, ppt_spring_recent_cv)))
ppt_summer_cv_difference <- terra::diff(terra::rast(list(ppt_summer_old_cv, ppt_summer_recent_cv)))
ppt_autumn_cv_difference <- terra::diff(terra::rast(list(ppt_autumn_old_cv, ppt_autumn_recent_cv)))


# changing the crs of these rasters to match the pixels of our mesh
# converting the lat long to epsg 6703km in km
# define a crs
epsg6703km <- paste(
  "+proj=aea +lat_0=23 +lon_0=-96 +lat_1=29.5",
  "+lat_2=45.5 +x_0=0 +y_0=0 +datum=NAD83",
  "+units=km +no_defs"
)

raster::crs(tmean_annual_difference) <- epsg6703km
raster::crs(tmean_spring_difference) <- epsg6703km
raster::crs(tmean_summer_difference) <- epsg6703km
raster::crs(tmean_autumn_difference) <- epsg6703km

raster::crs(tmean_annual_sd_difference) <- epsg6703km
raster::crs(tmean_spring_sd_difference) <- epsg6703km
raster::crs(tmean_summer_sd_difference) <- epsg6703km
raster::crs(tmean_autumn_sd_difference) <- epsg6703km

raster::crs(tmean_annual_cv_difference) <- epsg6703km 
raster::crs(tmean_spring_cv_difference) <- epsg6703km 
raster::crs(tmean_summer_cv_difference) <- epsg6703km 
raster::crs(tmean_autumn_cv_difference) <- epsg6703km 



raster::crs(ppt_annual_difference) <- epsg6703km
raster::crs(ppt_spring_difference) <- epsg6703km
raster::crs(ppt_summer_difference) <- epsg6703km
raster::crs(ppt_autumn_difference) <- epsg6703km

raster::crs(ppt_annual_sd_difference) <- epsg6703km
raster::crs(ppt_spring_sd_difference) <- epsg6703km
raster::crs(ppt_summer_sd_difference) <- epsg6703km
raster::crs(ppt_autumn_sd_difference) <- epsg6703km

raster::crs(ppt_annual_cv_difference) <- epsg6703km 
raster::crs(ppt_spring_cv_difference) <- epsg6703km 
raster::crs(ppt_summer_cv_difference) <- epsg6703km 
raster::crs(ppt_autumn_cv_difference) <- epsg6703km 



# ppt_autumn_cv_difference %>% st_transform(epsg6703km)
#   st_as_sf(coords = c("lon", "lat"), crs = 4326, remove = FALSE)
#   st_transform(epsg6703km) %>% 
#   mutate(
#     coords.x1 = st_coordinates(.)[, 1],
#     coords.x2 = = st_coordinates(.)[, 1],
# )


######### Next extracting these climate values at our points. ####
# we can extract them at the georeferenced coordinates and also across the grid of points we are using for prediction

# this is the coordinates of our samples

endo_herb <- endo_herb_georef %>% 
  filter(!is.na(Endo_status_liberal)) %>%
  filter(!is.na(Spp_code)) %>% 
  filter(!is.na(lon) & !is.na(year)) %>% 
  filter(lon>-110 ) %>% 
  filter(Country != "Canada" & !is.na(County)) 
  
coords_df <- endo_herb %>% 
  distinct(lat,lon) %>% 
  filter(!is.na(lat), !is.na(lon))
coords <- cbind(coords_df$lon, coords_df$lat)

prism_diff_df <- tibble(lon = coords[,1], lat = coords[,2],
                        tmean_annual_diff = unlist(terra::extract(tmean_annual_difference, coords)),
                        tmean_spring_diff = unlist(terra::extract(tmean_spring_difference, coords)),
                        tmean_summer_diff = unlist(terra::extract(tmean_summer_difference, coords)),
                        tmean_autumn_diff = unlist(terra::extract(tmean_autumn_difference, coords)),
                        ppt_annual_diff = unlist(terra::extract(ppt_annual_difference, coords)),
                        ppt_spring_diff = unlist(terra::extract(ppt_spring_difference, coords)),
                        ppt_summer_diff = unlist(terra::extract(ppt_summer_difference, coords)),
                        ppt_autumn_diff = unlist(terra::extract(ppt_autumn_difference, coords)),
                        
                        tmean_annual_sd_diff = unlist(terra::extract(tmean_annual_sd_difference, coords)),
                        tmean_spring_sd_diff = unlist(terra::extract(tmean_spring_sd_difference, coords)),
                        tmean_summer_sd_diff = unlist(terra::extract(tmean_summer_sd_difference, coords)),
                        tmean_autumn_sd_diff = unlist(terra::extract(tmean_autumn_sd_difference, coords)),
                        ppt_annual_sd_diff = unlist(terra::extract(ppt_annual_sd_difference, coords)),
                        ppt_spring_sd_diff = unlist(terra::extract(ppt_spring_sd_difference, coords)),
                        ppt_summer_sd_diff = unlist(terra::extract(ppt_summer_sd_difference, coords)),
                        ppt_autumn_sd_diff = unlist(terra::extract(ppt_autumn_sd_difference, coords)),
                        
                        tmean_annual_cv_diff = unlist(terra::extract(tmean_annual_cv_difference, coords)),
                        tmean_spring_cv_diff = unlist(terra::extract(tmean_spring_cv_difference, coords)),
                        tmean_summer_cv_diff = unlist(terra::extract(tmean_summer_cv_difference, coords)),
                        tmean_autumn_cv_diff = unlist(terra::extract(tmean_autumn_cv_difference, coords)),
                        ppt_annual_cv_diff = unlist(terra::extract(ppt_annual_cv_difference, coords)),
                        ppt_spring_cv_diff = unlist(terra::extract(ppt_spring_cv_difference, coords)),
                        ppt_summer_cv_diff = unlist(terra::extract(ppt_summer_cv_difference, coords)),
                        ppt_autumn_cv_diff = unlist(terra::extract(ppt_autumn_cv_difference, coords)))
write_csv(prism_diff_df, file = "prism_diff_df.csv")

prism_diff_df <- read_csv(file = "prism_diff_df.csv")
# this is the coordinates we are using for the spatial model prediction currently
# this is the set of data for which we want predictions

fit_lists <- readRDS(file = "~/Dropbox/Josh&Tom - shared/Endo_Herbarium/Model_Output/fit_lists_withscorer.Rds")

mesh_lists <- readRDS(file = "~/Dropbox/Josh&Tom - shared/Endo_Herbarium/Model_Output/mesh_list_withscorer.Rds")

bdry_polygon_list <- readRDS(file = "~/Dropbox/Josh&Tom - shared/Endo_Herbarium/Model_Output/boundary_polygon_list.Rds")



# THese are the pixels across space for which we will extract climate

vrt_aghy <- readRDS(file = "aghy_distribution_df.rds")
vrt_agpe <- readRDS(file = "agpe_distribution_df.rds")
vrt_elvi <- readRDS(file = "elvi_distribution_df.rds")


dd_crs <- "GEOGCRS[\"NAD83\",\n    DATUM[\"North American Datum 1983\",\n        ELLIPSOID[\"GRS 1980\",6378137,298.257222101004,\n            LENGTHUNIT[\"metre\",1]]],\n    PRIMEM[\"Greenwich\",0,\n        ANGLEUNIT[\"degree\",0.0174532925199433]],\n    CS[ellipsoidal,2],\n        AXIS[\"geodetic latitude (Lat)\",north,\n            ORDER[1],\n            ANGLEUNIT[\"degree\",0.0174532925199433]],\n        AXIS[\"geodetic longitude (Lon)\",east,\n            ORDER[2],\n            ANGLEUNIT[\"degree\",0.0174532925199433]],\n    ID[\"EPSG\",4269]]"

vrt_aghy_dd <- vrt_aghy  %>% 
  spTransform(dd_crs)

vrt_agpe_dd <- vrt_agpe  %>% 
  spTransform(dd_crs)

vrt_elvi_dd <- vrt_elvi  %>% 
  spTransform(dd_crs)



aghy_prism_diff_pred_df <- tibble(lon = vrt_aghy_dd@coords[,1], lat = vrt_aghy_dd@coords[,2],
                             tmean_annual_diff = unlist(terra::extract(tmean_annual_difference, vrt_aghy_dd@coords)),
                             tmean_spring_diff = unlist(terra::extract(tmean_spring_difference, vrt_aghy_dd@coords)),
                             tmean_summer_diff = unlist(terra::extract(tmean_summer_difference, vrt_aghy_dd@coords)),
                             tmean_autumn_diff = unlist(terra::extract(tmean_autumn_difference, vrt_aghy_dd@coords)),
                             ppt_annual_diff = unlist(terra::extract(ppt_annual_difference, vrt_aghy_dd@coords)),
                             ppt_spring_diff = unlist(terra::extract(ppt_spring_difference, vrt_aghy_dd@coords)),
                             ppt_summer_diff = unlist(terra::extract(ppt_summer_difference, vrt_aghy_dd@coords)),
                             ppt_autumn_diff = unlist(terra::extract(ppt_autumn_difference, vrt_aghy_dd@coords)),
                             
                             tmean_annual_sd_diff = unlist(terra::extract(tmean_annual_sd_difference, vrt_aghy_dd@coords)),
                             tmean_spring_sd_diff = unlist(terra::extract(tmean_spring_sd_difference, vrt_aghy_dd@coords)),
                             tmean_summer_sd_diff = unlist(terra::extract(tmean_summer_sd_difference, vrt_aghy_dd@coords)),
                             tmean_autumn_sd_diff = unlist(terra::extract(tmean_autumn_sd_difference, vrt_aghy_dd@coords)),
                             ppt_annual_sd_diff = unlist(terra::extract(ppt_annual_sd_difference, vrt_aghy_dd@coords)),
                             ppt_spring_sd_diff = unlist(terra::extract(ppt_spring_sd_difference, vrt_aghy_dd@coords)),
                             ppt_summer_sd_diff = unlist(terra::extract(ppt_summer_sd_difference, vrt_aghy_dd@coords)),
                             ppt_autumn_sd_diff = unlist(terra::extract(ppt_autumn_sd_difference, vrt_aghy_dd@coords)),
                             
                             tmean_annual_cv_diff = unlist(terra::extract(tmean_annual_cv_difference, vrt_aghy_dd@coords)),
                             tmean_spring_cv_diff = unlist(terra::extract(tmean_spring_cv_difference, vrt_aghy_dd@coords)),
                             tmean_summer_cv_diff = unlist(terra::extract(tmean_summer_cv_difference, vrt_aghy_dd@coords)),
                             tmean_autumn_cv_diff = unlist(terra::extract(tmean_autumn_cv_difference, vrt_aghy_dd@coords)),
                             ppt_annual_cv_diff = unlist(terra::extract(ppt_annual_cv_difference, vrt_aghy_dd@coords)),
                             ppt_spring_cv_diff = unlist(terra::extract(ppt_spring_cv_difference, vrt_aghy_dd@coords)),
                             ppt_summer_cv_diff = unlist(terra::extract(ppt_summer_cv_difference, vrt_aghy_dd@coords)),
                             ppt_autumn_cv_diff = unlist(terra::extract(ppt_autumn_cv_difference, vrt_aghy_dd@coords))) %>% 
  st_as_sf(coords = c("lon", "lat"), crs = epsg6703km, remove = FALSE) %>% 
  st_transform(epsg6703km) %>% 
  mutate(
    coords.x1 = st_coordinates(.)[, 1],
    coords.x2 = st_coordinates(.)[, 2]
  ) %>% 
  as.data.frame() %>% 
  mutate(origin = "prism") %>% 
  mutate(species = "A. hyemalis")


write_csv(aghy_prism_diff_pred_df, file = "aghy_prism_diff_pred_df.csv")

agpe_prism_diff_pred_df <- tibble(lon = vrt_agpe_dd@coords[,1], lat = vrt_agpe_dd@coords[,2],
                                  tmean_annual_diff = unlist(terra::extract(tmean_annual_difference, vrt_agpe_dd@coords)),
                                  tmean_spring_diff = unlist(terra::extract(tmean_spring_difference, vrt_agpe_dd@coords)),
                                  tmean_summer_diff = unlist(terra::extract(tmean_summer_difference, vrt_agpe_dd@coords)),
                                  tmean_autumn_diff = unlist(terra::extract(tmean_autumn_difference, vrt_agpe_dd@coords)),
                                  ppt_annual_diff = unlist(terra::extract(ppt_annual_difference, vrt_agpe_dd@coords)),
                                  ppt_spring_diff = unlist(terra::extract(ppt_spring_difference, vrt_agpe_dd@coords)),
                                  ppt_summer_diff = unlist(terra::extract(ppt_summer_difference, vrt_agpe_dd@coords)),
                                  ppt_autumn_diff = unlist(terra::extract(ppt_autumn_difference, vrt_agpe_dd@coords)),
                                  
                                  tmean_annual_sd_diff = unlist(terra::extract(tmean_annual_sd_difference, vrt_agpe_dd@coords)),
                                  tmean_spring_sd_diff = unlist(terra::extract(tmean_spring_sd_difference, vrt_agpe_dd@coords)),
                                  tmean_summer_sd_diff = unlist(terra::extract(tmean_summer_sd_difference, vrt_agpe_dd@coords)),
                                  tmean_autumn_sd_diff = unlist(terra::extract(tmean_autumn_sd_difference, vrt_agpe_dd@coords)),
                                  ppt_annual_sd_diff = unlist(terra::extract(ppt_annual_sd_difference, vrt_agpe_dd@coords)),
                                  ppt_spring_sd_diff = unlist(terra::extract(ppt_spring_sd_difference, vrt_agpe_dd@coords)),
                                  ppt_summer_sd_diff = unlist(terra::extract(ppt_summer_sd_difference, vrt_agpe_dd@coords)),
                                  ppt_autumn_sd_diff = unlist(terra::extract(ppt_autumn_sd_difference, vrt_agpe_dd@coords)),
                                  
                                  tmean_annual_cv_diff = unlist(terra::extract(tmean_annual_cv_difference, vrt_agpe_dd@coords)),
                                  tmean_spring_cv_diff = unlist(terra::extract(tmean_spring_cv_difference, vrt_agpe_dd@coords)),
                                  tmean_summer_cv_diff = unlist(terra::extract(tmean_summer_cv_difference, vrt_agpe_dd@coords)),
                                  tmean_autumn_cv_diff = unlist(terra::extract(tmean_autumn_cv_difference, vrt_agpe_dd@coords)),
                                  ppt_annual_cv_diff = unlist(terra::extract(ppt_annual_cv_difference, vrt_agpe_dd@coords)),
                                  ppt_spring_cv_diff = unlist(terra::extract(ppt_spring_cv_difference, vrt_agpe_dd@coords)),
                                  ppt_summer_cv_diff = unlist(terra::extract(ppt_summer_cv_difference, vrt_agpe_dd@coords)),
                                  ppt_autumn_cv_diff = unlist(terra::extract(ppt_autumn_cv_difference, vrt_agpe_dd@coords))) %>% 
  st_as_sf(coords = c("lon", "lat"), crs = 4326, remove = FALSE) %>% 
  st_transform(epsg6703km) %>% 
  mutate(
    coords.x1 = st_coordinates(.)[, 1],
    coords.x2 = st_coordinates(.)[, 2]
  ) %>% 
  as.data.frame() %>% 
  mutate(origin = "prism") %>% 
  mutate(species = "A. perennans")


write_csv(agpe_prism_diff_pred_df, file = "agpe_prism_diff_pred_df.csv")


elvi_prism_diff_pred_df <- tibble(lon = vrt_elvi_dd@coords[,1], lat = vrt_elvi_dd@coords[,2],
                                  tmean_annual_diff = unlist(terra::extract(tmean_annual_difference, vrt_elvi_dd@coords)),
                                  tmean_spring_diff = unlist(terra::extract(tmean_spring_difference, vrt_elvi_dd@coords)),
                                  tmean_summer_diff = unlist(terra::extract(tmean_summer_difference, vrt_elvi_dd@coords)),
                                  tmean_autumn_diff = unlist(terra::extract(tmean_autumn_difference, vrt_elvi_dd@coords)),
                                  ppt_annual_diff = unlist(terra::extract(ppt_annual_difference, vrt_elvi_dd@coords)),
                                  ppt_spring_diff = unlist(terra::extract(ppt_spring_difference, vrt_elvi_dd@coords)),
                                  ppt_summer_diff = unlist(terra::extract(ppt_summer_difference, vrt_elvi_dd@coords)),
                                  ppt_autumn_diff = unlist(terra::extract(ppt_autumn_difference, vrt_elvi_dd@coords)),
                                  
                                  tmean_annual_sd_diff = unlist(terra::extract(tmean_annual_sd_difference, vrt_elvi_dd@coords)),
                                  tmean_spring_sd_diff = unlist(terra::extract(tmean_spring_sd_difference, vrt_elvi_dd@coords)),
                                  tmean_summer_sd_diff = unlist(terra::extract(tmean_summer_sd_difference, vrt_elvi_dd@coords)),
                                  tmean_autumn_sd_diff = unlist(terra::extract(tmean_autumn_sd_difference, vrt_elvi_dd@coords)),
                                  ppt_annual_sd_diff = unlist(terra::extract(ppt_annual_sd_difference, vrt_elvi_dd@coords)),
                                  ppt_spring_sd_diff = unlist(terra::extract(ppt_spring_sd_difference, vrt_elvi_dd@coords)),
                                  ppt_summer_sd_diff = unlist(terra::extract(ppt_summer_sd_difference, vrt_elvi_dd@coords)),
                                  ppt_autumn_sd_diff = unlist(terra::extract(ppt_autumn_sd_difference, vrt_elvi_dd@coords)),
                                  
                                  tmean_annual_cv_diff = unlist(terra::extract(tmean_annual_cv_difference, vrt_elvi_dd@coords)),
                                  tmean_spring_cv_diff = unlist(terra::extract(tmean_spring_cv_difference, vrt_elvi_dd@coords)),
                                  tmean_summer_cv_diff = unlist(terra::extract(tmean_summer_cv_difference, vrt_elvi_dd@coords)),
                                  tmean_autumn_cv_diff = unlist(terra::extract(tmean_autumn_cv_difference, vrt_elvi_dd@coords)),
                                  ppt_annual_cv_diff = unlist(terra::extract(ppt_annual_cv_difference, vrt_elvi_dd@coords)),
                                  ppt_spring_cv_diff = unlist(terra::extract(ppt_spring_cv_difference, vrt_elvi_dd@coords)),
                                  ppt_summer_cv_diff = unlist(terra::extract(ppt_summer_cv_difference, vrt_elvi_dd@coords)),
                                  ppt_autumn_cv_diff = unlist(terra::extract(ppt_autumn_cv_difference, vrt_elvi_dd@coords))) %>% 
  st_as_sf(coords = c("lon", "lat"), crs = 4326, remove = FALSE) %>% 
  st_transform(epsg6703km) %>% 
  mutate(
    coords.x1 = st_coordinates(.)[, 1],
    coords.x2 = st_coordinates(.)[, 2]
  ) %>% 
  as.data.frame() %>% 
  mutate(origin = "prism") %>% 
  mutate(species = "E. virginicus")

write_csv(elvi_prism_diff_pred_df, file = "elvi_prism_diff_pred_df.csv")


# Plotting the change in climate at our pixel values

prism_diff_pred_df <- aghy_prism_diff_pred_df %>% 
  full_join(agpe_prism_diff_pred_df) %>% 
  full_join(elvi_prism_diff_pred_df) %>% 
  pivot_longer(contains("_diff"))
  

ggplot(filter(prism_diff_pred_df, grepl("tmean", name) & !grepl("_cv_", name) ))+
  geom_point(aes(x = lon, y = lat, color = value), shape = 15) +
  facet_wrap(~species + name)
  scale_fill_b




ggplot(prism_diff_pred_df)+
  geom_point(aes(x = lon, y = lat, color = ppt_annual_diff))




epsg6703degree <- "+proj=aea +lat_0=23 +lon_0=-96 +lat_1=29.5 +lat_2=45.5 +x_0=0 +y_0=0 +datum=NAD83 +units=km +no_defs"
vrt_list <- list()

species_codes <- c("AGHY", "AGPE", "ELVI")
species_names <- c("A. hyemalis", "A. perennans", "E. virginicus")


for (s in 1:3){
  vrt <- NA
  mesh <- mesh_lists[[s]]
  bdry_polygon <- bdry_polygon_list[[s]]
  
  vrt_list[[species_codes[s]]] <- inlabru::fm_pixels(mesh, mask = bdry_polygon, format = "sp", dims = c(40,40))
}

vrt_df <- vrt_df %>% 
  st_transform(crs = epsg6703degree) %>% 
  mutate(lat= st_coordinates(vrt_df)[,1],
         lon = st_coordinates(vrt_df)[,2])

ggplot()+
  gg(mesh_lists[[1]])+
  gg(vrt_list[[1]], color = "red")

ggplot()+
  gg(mesh_lists[[2]])+
  gg(vrt_list[[2]], color = "red")

ggplot()+
  gg(mesh_lists[[3]])+
  gg(vrt_list[[3]], color = "red")


vrt_list[[1]]@coords

prism_diff_pred_df <- tibble(lon = vrt_list[[1]]@coords[,1], lat = vrt_list[[1]]@coords[,2],
                             tmean_annual_diff = unlist(terra::extract(tmean_annual_difference, vrt_list[[1]]@coords)),
                             tmean_spring_diff = unlist(terra::extract(tmean_spring_difference, vrt_list[[1]]@coords)),
                             tmean_summer_diff = unlist(terra::extract(tmean_summer_difference, vrt_list[[1]]@coords)),
                             tmean_autumn_diff = unlist(terra::extract(tmean_autumn_difference, vrt_list[[1]]@coords)),
                             ppt_annual_diff = unlist(terra::extract(ppt_annual_difference, vrt_list[[1]]@coords)),
                             ppt_spring_diff = unlist(terra::extract(ppt_spring_difference, vrt_list[[1]]@coords)),
                             ppt_summer_diff = unlist(terra::extract(ppt_summer_difference, vrt_list[[1]]@coords)),
                             ppt_autumn_diff = unlist(terra::extract(ppt_autumn_difference, vrt_list[[1]]@coords)),
                             
                             tmean_annual_sd_diff = unlist(terra::extract(tmean_annual_sd_difference, vrt_list[[1]]@coords)),
                             tmean_spring_sd_diff = unlist(terra::extract(tmean_spring_sd_difference, vrt_list[[1]]@coords)),
                             tmean_summer_sd_diff = unlist(terra::extract(tmean_summer_sd_difference, vrt_list[[1]]@coords)),
                             tmean_autumn_sd_diff = unlist(terra::extract(tmean_autumn_sd_difference, vrt_list[[1]]@coords)),
                             ppt_annual_sd_diff = unlist(terra::extract(ppt_annual_sd_difference, vrt_list[[1]]@coords)),
                             ppt_spring_sd_diff = unlist(terra::extract(ppt_spring_sd_difference, vrt_list[[1]]@coords)),
                             ppt_summer_sd_diff = unlist(terra::extract(ppt_summer_sd_difference, vrt_list[[1]]@coords)),
                             ppt_autumn_sd_diff = unlist(terra::extract(ppt_autumn_sd_difference, vrt_list[[1]]@coords)),
                             
                             tmean_annual_cv_diff = unlist(terra::extract(tmean_annual_cv_difference, vrt_list[[1]]@coords)),
                             tmean_spring_cv_diff = unlist(terra::extract(tmean_spring_cv_difference, vrt_list[[1]]@coords)),
                             tmean_summer_cv_diff = unlist(terra::extract(tmean_summer_cv_difference, vrt_list[[1]]@coords)),
                             tmean_autumn_cv_diff = unlist(terra::extract(tmean_autumn_cv_difference, vrt_list[[1]]@coords)),
                             ppt_annual_cv_diff = unlist(terra::extract(ppt_annual_cv_difference, vrt_list[[1]]@coords)),
                             ppt_spring_cv_diff = unlist(terra::extract(ppt_spring_cv_difference, vrt_list[[1]]@coords)),
                             ppt_summer_cv_diff = unlist(terra::extract(ppt_summer_cv_difference, vrt_list[[1]]@coords)),
                             ppt_autumn_cv_diff = unlist(terra::extract(ppt_autumn_cv_difference, vrt_list[[1]]@coords))) %>% 
  na.omit()
write_csv(prism_diff_pred_df, file = "prism_diff_pred_df.csv")

ggplot(prism_diff_pred_df)+
  geom_point(aes(x = lon, y = lat, color = ppt_annual_diff))
ggplot(prism_diff_pred_df)+
  geom_point(aes(x = lon, y = lat, color = ppt_annual_cv_diff))




###



  
ppt_annual_old <- terra::rast(pd_stack(prism_archive_subset(type = "ppt", temp_period = "monthly", year = 1895:1925)))


# ppt isn't right as is, need to get cumulative values for year or for season
ppt_spring_old <- terra::rast(pd_stack(prism_archive_subset(type = "ppt", temp_period = "monthly", year = 1895:1925, mon = 1:3)))


# calculating the change in annual values
tmean_annual_recent_avg <- terra::mean(tmean_annual_recent)
tmean_annual_old_avg <- terra::mean(tmean_annual_old)
tmean_spring_recent_avg <- terra::mean(tmean_spring_recent)
tmean_spring_old_avg <- terra::mean(tmean_spring_old)


ppt_annual_recent_avg <- terra::mean(ppt_annual_recent)
ppt_annual_old_avg <- terra::mean(ppt_annual_old)
ppt_spring_recent_avg <- terra::mean(ppt_spring_recent)
ppt_spring_old_avg <- terra::mean(ppt_spring_old)

tmean_spring_change_stack <- terra::rast(list(tmean_spring_old_avg, tmean_spring_recent_avg))
tmean_annual_change_stack <- terra::rast(list(tmean_annual_old_avg, tmean_annual_recent_avg))

ppt_spring_change_stack <- terra::rast(list(ppt_spring_old_avg, ppt_spring_recent_avg))
ppt_annual_change_stack <- terra::rast(list(ppt_annual_old_avg, ppt_annual_recent_avg))


tmean_spring_change <- terra::diff(tmean_spring_change_stack)
tmean_annual_change <- terra::diff(tmean_annual_change_stack)

ppt_spring_change <- terra::diff(ppt_spring_change_stack)
ppt_annual_change <- terra::diff(ppt_annual_change_stack)

plot(tmean_spring_change)
plot(tmean_annual_change)

plot(ppt_spring_change)
plot(ppt_annual_change)


# extracting the climate values at our herbarium specimens' coordinates for each year

tmean_stack <- terra::rast(pd_stack(prism_archive_subset(type = "tmean", temp_period = "annual", year = 1970:2020)))
ppt_stack <- terra::rast(pd_stack(prism_archive_subset(type = "ppt", temp_period = "monthly", year = 1895:2019)))


coords_df <- endo_herb_georef %>% 
  # distinct(lat,lon) %>% 
  filter(!is.na(lat), !is.na(lon))
coords <- cbind(coords_df$lon, coords_df$lat)

tmean_extract <- terra::extract(tmean_stack, coords)
tmean_extract$lon <-  coords_df$lon
tmean_extract$lat <-  coords_df$lat
tmean_extract$species <- coords_df$Spp_code
tmean_extract$endo <- coords_df$Endo_status_liberal

ggplot(tmean_extract)+geom_point(aes(x = PRISM_tmean_stable_4kmM3_2020))



lag_multiple <- function(x,name, n_vec){
  #apply lag function across a set of columns
  purrr::map(n_vec, lag, x = x) %>% 
    # name the columns according to "name" and the number of lag periods
    set_names(paste0(name,"_", "lag", n_vec)) %>% 
    as_tibble()
}

tmean_longer <- tmean_extract %>% 
  pivot_longer(cols = c(-lat, -lon)) %>% 
  mutate(year = year(ym(stringr::str_sub(name,start = 26, end = 31))),
         month = month(ym(stringr::str_sub(name,start = 26, end = 31)))) %>% 
  group_by(lat,lon) %>% 
  arrange(year, month) %>% 
  mutate(lag_multiple(value,"temp", 0:140)) %>% 
  ungroup() %>% 
  mutate(mean_1month_temp = temp_lag0,
         mean_2month_temp = rowMeans(select(., temp_lag0:temp_lag1), na.rm = TRUE),
         mean_3month_temp = rowMeans(select(., temp_lag0:temp_lag2), na.rm = TRUE),
         mean_4month_temp = rowMeans(select(., temp_lag0:temp_lag3), na.rm = TRUE),
         mean_5month_temp = rowMeans(select(., temp_lag0:temp_lag4), na.rm = TRUE),
         mean_6month_temp = rowMeans(select(., temp_lag0:temp_lag5), na.rm = TRUE),
         mean_7month_temp = rowMeans(select(., temp_lag0:temp_lag6), na.rm = TRUE),
         mean_8month_temp = rowMeans(select(., temp_lag0:temp_lag7), na.rm = TRUE),
         mean_9month_temp = rowMeans(select(., temp_lag0:temp_lag8), na.rm = TRUE),
         mean_10month_temp = rowMeans(select(., temp_lag0:temp_lag9), na.rm = TRUE),
         mean_11month_temp = rowMeans(select(., temp_lag0:temp_lag10), na.rm = TRUE),
         mean_annual_temp = rowMeans(select(., temp_lag0:temp_lag11), na.rm = TRUE),
         mean_24month_temp = rowMeans(select(., temp_lag0:temp_lag23), na.rm = TRUE),
         mean_decade_temp = rowMeans(select(., temp_lag0:temp_lag139), na.rm = TRUE)) %>% 
  select(lat, lon, year, month, mean_1month_temp, mean_2month_temp, mean_3month_temp, mean_annual_temp, mean_24month_temp, mean_decade_temp)

write_csv(tmean_longer, file = "~/Dropbox/Josh&Tom - shared/Endo_Herbarium/DigitizedHerbariumRecords/PRISM_tmean_longer.csv")

ppt_extract <- terra::extract(ppt_stack, coords)
ppt_extract$lon = coords_df$lon
ppt_extract$lat = coords_df$lat

ppt_longer <- ppt_extract %>% 
  pivot_longer(cols = c(-lat, -lon)) %>% 
  mutate(year =   year(ym(stringr::str_sub(name,start = 24, end = 29))),
         month = month(ym(stringr::str_sub(name,start = 24, end = 29)))) %>% 
  group_by(lat,lon) %>% 
  arrange(year, month) %>% 
  mutate(lag_multiple(value,"ppt", 0:140)) %>% 
  ungroup() %>% 
  mutate(total_1month_ppt = ppt_lag0,
         total_2month_ppt = rowSums(select(., ppt_lag0:ppt_lag1), na.rm = TRUE),
         total_3month_ppt = rowSums(select(., ppt_lag0:ppt_lag2), na.rm = TRUE),
         total_4month_ppt = rowSums(select(., ppt_lag0:ppt_lag3), na.rm = TRUE),
         total_5month_ppt = rowSums(select(., ppt_lag0:ppt_lag4), na.rm = TRUE),
         total_6month_ppt = rowSums(select(., ppt_lag0:ppt_lag5), na.rm = TRUE),
         total_7month_ppt = rowSums(select(., ppt_lag0:ppt_lag6), na.rm = TRUE),
         total_8month_ppt = rowSums(select(., ppt_lag0:ppt_lag7), na.rm = TRUE),
         total_9month_ppt = rowSums(select(., ppt_lag0:ppt_lag8), na.rm = TRUE),
         total_10month_ppt = rowSums(select(., ppt_lag0:ppt_lag9), na.rm = TRUE),
         total_11month_ppt = rowSums(select(., ppt_lag0:ppt_lag10), na.rm = TRUE),
         total_annual_ppt = rowSums(select(., ppt_lag0:ppt_lag11), na.rm = TRUE),
         total_24month_ppt = rowSums(select(., ppt_lag0:ppt_lag23), na.rm = TRUE),
         total_decade_ppt = rowSums(select(., ppt_lag0:ppt_lag139), na.rm = TRUE)) %>% 
  select(lat, lon, year, month, total_1month_ppt, total_2month_ppt, total_annual_ppt, total_24month_ppt, total_decade_ppt)
write_csv(ppt_longer, file = "~/Dropbox/Josh&Tom - shared/Endo_Herbarium/DigitizedHerbariumRecords/PRISM_ppt_longer.csv")


test <- tmean_longer %>% 
  filter(lon == unique(tmean_longer$lon)[1]) %>%    
  mutate(date = ym(paste0(year,"-",month))) %>% 
  ggplot() +
  geom_path(aes(x = date, y = mean_1month_temp))+
  geom_path(aes(x = date, y = mean_annual_temp), color = "red")+
  geom_path(aes(x = date, y = mean_decade_temp), color = "purple")

test
test <- ppt_longer %>% 
  filter(lon == unique(ppt_longer$lon)[1]) %>%    
  mutate(date = ym(paste0(year,"-",month))) %>% 
ggplot() +
  geom_path(aes(x = date, y = total_1month_ppt))+
  geom_path(aes(x = date, y = total_annual_ppt), color = "blue")+
  geom_path(aes(x = date, y = total_decade_ppt), color = "purple")


# merging the climate data with the herbarium scores
endo_herb_georef_climate <- endo_herb_georef %>% 
  left_join(tmean_longer) %>% 
  left_join(ppt_longer)
write_csv(endo_herb_georef_climate, file = "~/Dropbox/Josh&Tom - shared/Endo_Herbarium/DigitizedHerbariumRecords/endo_herb_georef_climate.csv")


