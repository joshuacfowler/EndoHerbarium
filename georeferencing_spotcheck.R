# Title: Spotchecking Georeferencing
# Purpose: Imports georeferenced herbarium records with endophyte scores from endoherbarium_records_merge.R and generates an interactive map to compare named location and lat, long
# Authors: Joshua Fowler 
# Updated: May 29, 2022

library(tidyverse) # for data manipulation and ggplot
library(leaflet) # for making interactive maps

################################################################################
############ Read in georeferenced herbarium records ############################### 
################################################################################

endo_herb_georef <- read_csv(file = "~/Dropbox/Josh&Tom - shared/Endo_Herbarium/DigitizedHerbariumRecords/endo_herb_georef.csv") %>% 
  filter(!is.na(Endo_status_liberal))

# setting endophyte color palette
spp_pal <- colorFactor(palette = "Accent",
                       domain = endo_herb_georef$Spp_code)

map <- leaflet() %>% 
  addTiles() %>% 
  addCircleMarkers(~lon, ~ lat, data = endo_herb_georef,
                   color = ~spp_pal(Spp_code),
                   popup = ~paste("Species:", Spp_code, "<br>",
                                  "Specimen_id:", Sample_id, "<br>",
                                  "Lat,Long:", lat, ",", lon, "<br>",
                                  "Locality String:", location_string))


map
