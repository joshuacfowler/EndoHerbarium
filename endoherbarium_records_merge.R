# Title: Digital Herbarium Records and Sample Integration
# Purpose: Imports and merges downloaded digitized herbarium records with the Endo_Herbarium database
# Authors: Joshua Fowler
# Updated: July 17, 2019

library(tidyverse)
library(googlesheets)
library(readxl)

# Read in digitized herbarium records
dig_AGHY_UTAustin <- read_csv(file = "~joshuacfowler/Dropbox/Josh&Tom - shared/Endo_Herbarium/DigitizedHerbariumRecords/UTAustinAGHYrecords/occurrences.csv")
# dig_AGHY_BRIT <- read_csv(file = "")

# Read in our transcribed datasheet
specimen_info <- read_csv(file = "~joshuacfowler/Dropbox/Josh&Tom - shared/Endo_Herbarium/testfromgsheets_specimen_info.csv")

# I'm gonna clean up the id variable here to be able to merge this with the digitized records for AGHY
AGHY_specimen_info <- specimen_info %>% 
  filter(str_detect(Specimen_id, 'AGHY')) %>% 
  left_join(dig_AGHY_UTAustin, by = c("Institution_specimen_id" = "catalogNumber"))

  

# Read in our scores for the samples
# We will need to clean this up to remove the samples that we haven't done yet
sample_info <- read_xlsx(path = "~joshuacfowler/Dropbox/Josh&Tom - shared/Endo_Herbarium/Endo_Herbarium_testfilefromgooglesheets.xlsx", sheet = "Sample info")

# We can fix the id's to be able to merge the two datasets so that we have a full data set of specimen info (hopefully)
AGHY_spec <- 