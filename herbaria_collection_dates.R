## Title: Range limit project herbarium record collection date exploration
## Purpose: Takes recent herbarium data and identifies collection dates based on range position
## Authors: Josh

# getwd()
# setwd("C:/Users/MillerLab/Desktop/Josh/EndoHerbarium")

library(tidyverse)
library(shiny)
library(leaflet)
library(viridis)
library(raster)
library(maps)
# Load in our downloaded data
POAU <- read.csv("SymbOutput_2019-04-09_073934_DwC-A_POAU/occurrences.csv")
SPNI <- read.csv("SymbOutput_2019-04-09_074053_DwC-A_SPNI/occurrences.csv")
AGPE <- read.csv("SymbOutput_2019-04-09_074419_DwC-A_AGPE/occurrences.csv")
AGHY <- read.csv("SymbOutput_2019-04-09_074328_DwC-A_AGHY/occurrences.csv")
ELCA <- read.csv("SymbOutput_2019-04-09_074622_DwC-A_ELCA/occurrences.csv")
ELVI <- read.csv("SymbOutput_2019-04-09_074502_DwC-A_ELVI/occurrences.csv")

# now I'm going to group the data by location
# POAU
POAU_loc <- POAU %>% 
  filter(!is.na(decimalLatitude), !is.na(decimalLongitude)) %>% 
  mutate(bin = cut(decimalLongitude, breaks = 20)) %>% 
  group_by(bin) %>% 
  summarize(mean_month = mean(month, na.rm = TRUE),
            mean_long = mean(decimalLongitude, na.rm = TRUE),
            samplesize = n())
POAU_loc

POAU_dates <- ggplot(data = POAU_loc)+
  geom_point(aes(x = mean_long, y = mean_month, size = samplesize))+
  geom_smooth(aes(x = mean_long, y = mean_month), method = "glm")
POAU_dates

POAU_state <- POAU %>%
  group_by(stateProvince) %>% 
  summarize(mean_month = mean(month, na.rm = TRUE),
            samplesize = n())
level_order <- c("New Mexico", "Texas", "TEXAS", "Oklahoma", "Arkansas", "Louisiana", "Mississippi","Tennessee", "Kentucky", "Ohio", "Alabama", "Georgia", "Florida", "North Carolina", "South Carolina", "Virginia", "Maryland")

POAU_dates_states <- ggplot(data = POAU_state)+
  geom_point(aes(x = factor(stateProvince, levels = level_order), y = mean_month, size = samplesize))
POAU_dates_states


# I'm gonna try to make an interactive county map

POAU_county<- POAU %>%
  group_by(stateProvince, county) %>% 
  summarize(mean_month = mean(month, na.rm = TRUE),
            samplesize = n())

USA <- getData("GADM", country = "usa", level = 2)
temp <- merge(USA, POAU_county,
                    by.x = c("NAME_1", "NAME_2"), by.y = c("stateProvince", "county"),
                    all.x = TRUE)
# create a color palette
mypal <- colorNumeric(palette = "viridis", domain = temp$mean_month, na.color = "grey")

leaflet() %>% 
  addProviderTiles("OpenStreetMap.Mapnik") %>%
  setView(lat = 39.8283, lng = -98.5795, zoom = 4) %>% 
  addPolygons(data = temp, stroke = FALSE, smoothFactor = 0.2, fillOpacity = 0.8,
              fillColor = ~mypal(temp$mean_month),
              popup = paste("Region: ", temp$NAME_2, "<br>",
                            "Mean Month: ", temp$mean_month, "<br>",
                            "Samples: ", temp$samplesize, "<br>")) %>%
  addLegend(position = "bottomleft", pal = mypal, values = temp$mean_month,
            title = "Mean Month",
            opacity = 1)

# SPNI
SPNI_loc <- SPNI %>% 
  filter(!is.na(decimalLatitude), !is.na(decimalLongitude)) %>% 
  mutate(bin = cut(decimalLongitude, breaks = 40)) %>% 
  group_by(bin) %>% 
  summarize(mean_month = mean(month, na.rm = TRUE),
            mean_long = mean(decimalLongitude, na.rm = TRUE),
            samplesize = n())
SPNI_loc

SPNI_dates <- ggplot(data = SPNI_loc)+
  geom_point(aes(x = mean_long, y = mean_month, size = samplesize))+
  geom_smooth(aes(x = mean_long, y = mean_month), method = "glm")
SPNI_dates

SPNI_state <- SPNI %>%
  group_by(stateProvince) %>% 
  summarize(mean_month = mean(month, na.rm = TRUE),
            samplesize = n())
level_order <- c("New Mexico", "Texas", "TEXAS", "Oklahoma", "Arkansas", "Louisiana", "Mississippi","Tennessee", "Kentucky", "Ohio", "West virginia", "Indiana", "Michigan", "Alabama", "Georgia", "Pennsylvania", "Florida", "North Carolina", "South Carolina", "Virginia", "Maryland")

SPNI_dates_states <- ggplot(data = SPNI_state)+
  geom_point(aes(x = factor(stateProvince, levels = level_order), y = mean_month, size = samplesize))
SPNI_dates_states

# I'm gonna try to make an interactive county map

SPNI_county<- SPNI %>%
  group_by(stateProvince, county) %>% 
  summarize(mean_month = mean(month, na.rm = TRUE),
            samplesize = n())

USA <- getData("GADM", country = "usa", level = 2)
temp <- merge(USA, SPNI_county,
              by.x = c("NAME_1", "NAME_2"), by.y = c("stateProvince", "county"),
              all.x = TRUE)
# create a color palette
mypal <- colorNumeric(palette = "viridis", domain = temp$mean_month, na.color = "grey")

leaflet() %>% 
  addProviderTiles("OpenStreetMap.Mapnik") %>%
  setView(lat = 39.8283, lng = -98.5795, zoom = 4) %>% 
  addPolygons(data = temp, stroke = FALSE, smoothFactor = 0.2, fillOpacity = 0.8,
              fillColor = ~mypal(temp$mean_month),
              popup = paste("Region: ", temp$NAME_2, "<br>",
                            "Mean Month: ", temp$mean_month, "<br>",
                            "Samples: ", temp$samplesize, "<br>")) %>%
  addLegend(position = "bottomleft", pal = mypal, values = temp$mean_month,
            title = "Mean Month",
            opacity = 1)




# AGHY
AGHY_loc <- AGHY %>% 
  filter(!is.na(decimalLatitude), !is.na(decimalLongitude)) %>% 
  mutate(bin = cut(decimalLongitude, breaks = 20)) %>% 
  group_by(bin) %>% 
  summarize(mean_month = mean(month, na.rm = TRUE),
            mean_long = mean(decimalLongitude, na.rm = TRUE),
            samplesize = n())
AGHY_loc

AGHY_dates <- ggplot(data = AGHY_loc)+
  geom_point(aes(x = mean_long, y = mean_month, size = samplesize))+
  geom_smooth(aes(x = mean_long, y = mean_month), method = "glm")
AGHY_dates

AGHY_state <- AGHY %>%
  group_by(stateProvince) %>% 
  summarize(mean_month = mean(month, na.rm = TRUE),
            samplesize = n())
level_order <- c("New Mexico", "Texas", "TEXAS", "Oklahoma", "Arkansas", "Louisiana", "Mississippi","Tennessee", "Kentucky", "Ohio", "Alabama", "Georgia", "Florida", "North Carolina", "South Carolina", "Virginia", "Maryland")

AGHY_dates_states <- ggplot(data = AGHY_state)+
  geom_point(aes(x = factor(stateProvince, levels = level_order), y = mean_month, size = samplesize))
AGHY_dates_states

# I'm gonna try to make an interactive county map

AGHY_county<- AGHY %>%
  group_by(stateProvince, county) %>% 
  summarize(mean_month = mean(month, na.rm = TRUE),
            samplesize = n())

USA <- getData("GADM", country = "usa", level = 2)
temp <- merge(USA, AGHY_county,
              by.x = c("NAME_1", "NAME_2"), by.y = c("stateProvince", "county"),
              all.x = TRUE)
# create a color palette
mypal <- colorNumeric(palette = "viridis", domain = temp$mean_month, na.color = "grey")

leaflet() %>% 
  addProviderTiles("OpenStreetMap.Mapnik") %>%
  setView(lat = 39.8283, lng = -98.5795, zoom = 4) %>% 
  addPolygons(data = temp, stroke = FALSE, smoothFactor = 0.2, fillOpacity = 0.8,
              fillColor = ~mypal(temp$mean_month),
              popup = paste("Region: ", temp$NAME_2, "<br>",
                            "Mean Month: ", temp$mean_month, "<br>",
                            "Samples: ", temp$samplesize, "<br>")) %>%
  addLegend(position = "bottomleft", pal = mypal, values = temp$mean_month,
            title = "Mean Month",
            opacity = 1)


# AGPE
AGPE_loc <- AGPE %>% 
  filter(!is.na(decimalLatitude), !is.na(decimalLongitude)) %>% 
  mutate(bin = cut(decimalLongitude, breaks = 100)) %>% 
  group_by(bin) %>% 
  summarize(mean_month = mean(month, na.rm = TRUE),
            mean_long = mean(decimalLongitude, na.rm = TRUE),
            samplesize = n())
AGPE_loc

AGPE_dates <- ggplot(data = AGPE_loc)+
  geom_point(aes(x = mean_long, y = mean_month, size = samplesize))+
  geom_smooth(aes(x = mean_long, y = mean_month), method = "glm")
AGPE_dates

AGPE_state <- AGPE %>%
  group_by(stateProvince) %>% 
  summarize(mean_month = mean(month, na.rm = TRUE),
            samplesize = n())
level_order <- c("New Mexico", "Texas", "TEXAS", "Oklahoma", "Arkansas", "Louisiana", "Mississippi","Tennessee", "Kentucky", "Ohio", "Alabama", "Georgia", "Florida", "North Carolina", "South Carolina", "Virginia", "Maryland")

AGPE_dates_states <- ggplot(data = AGPE_state)+
  geom_point(aes(x = factor(stateProvince, levels = level_order), y = mean_month, size = samplesize))
AGPE_dates_states

# I'm gonna try to make an interactive county map

AGPE_county<- AGPE %>%
  group_by(stateProvince, county) %>% 
  summarize(mean_month = mean(month, na.rm = TRUE),
            samplesize = n())

USA <- getData("GADM", country = "usa", level = 2)
temp <- merge(USA, AGPE_county,
              by.x = c("NAME_1", "NAME_2"), by.y = c("stateProvince", "county"),
              all.x = TRUE)
# create a color palette
mypal <- colorNumeric(palette = "viridis", domain = temp$mean_month, na.color = "grey")

leaflet() %>% 
  addProviderTiles("OpenStreetMap.Mapnik") %>%
  setView(lat = 39.8283, lng = -98.5795, zoom = 4) %>% 
  addPolygons(data = temp, stroke = FALSE, smoothFactor = 0.2, fillOpacity = 0.8,
              fillColor = ~mypal(temp$mean_month),
              popup = paste("Region: ", temp$NAME_2, "<br>",
                            "Mean Month: ", temp$mean_month, "<br>",
                            "Samples: ", temp$samplesize, "<br>")) %>%
  addLegend(position = "bottomleft", pal = mypal, values = temp$mean_month,
            title = "Mean Month",
            opacity = 1)

# ELVI
ELVI_loc <- ELVI %>% 
  filter(!is.na(decimalLatitude), !is.na(decimalLongitude)) %>% 
  mutate(bin = cut(decimalLongitude, breaks = 100)) %>% 
  group_by(bin) %>% 
  summarize(mean_month = mean(month, na.rm = TRUE),
            mean_long = mean(decimalLongitude, na.rm = TRUE),
            samplesize = n())
ELVI_loc

ELVI_dates <- ggplot(data = ELVI_loc)+
  geom_point(aes(x = mean_long, y = mean_month, size = samplesize))+
  geom_smooth(aes(x = mean_long, y = mean_month), method = "glm")
ELVI_dates

ELVI_state <- ELVI %>%
  group_by(stateProvince) %>% 
  summarize(mean_month = mean(month, na.rm = TRUE),
            samplesize = n())
level_order <- c("New Mexico", "Texas", "TEXAS", "Oklahoma", "Arkansas", "Louisiana", "Mississippi","Tennessee", "Kentucky", "Ohio", "Alabama", "Georgia", "Florida", "North Carolina", "South Carolina", "Virginia", "Maryland")

ELVI_dates_states <- ggplot(data = ELVI_state)+
  geom_point(aes(x = factor(stateProvince, levels = level_order), y = mean_month, size = samplesize))
ELVI_dates_states

# I'm gonna try to make an interactive county map

ELVI_county<- ELVI %>%
  group_by(stateProvince, county) %>% 
  summarize(mean_month = mean(month, na.rm = TRUE),
            samplesize = n())

USA <- getData("GADM", country = "usa", level = 2)
temp <- merge(USA, ELVI_county,
              by.x = c("NAME_1", "NAME_2"), by.y = c("stateProvince", "county"),
              all.x = TRUE)
# create a color palette
mypal <- colorNumeric(palette = "viridis", domain = temp$mean_month, na.color = "grey")

leaflet() %>% 
  addProviderTiles("OpenStreetMap.Mapnik") %>%
  setView(lat = 39.8283, lng = -98.5795, zoom = 4) %>% 
  addPolygons(data = temp, stroke = FALSE, smoothFactor = 0.2, fillOpacity = 0.8,
              fillColor = ~mypal(temp$mean_month),
              popup = paste("Region: ", temp$NAME_2, "<br>",
                            "Mean Month: ", temp$mean_month, "<br>",
                            "Samples: ", temp$samplesize, "<br>")) %>%
  addLegend(position = "bottomleft", pal = mypal, values = temp$mean_month,
            title = "Mean Month",
            opacity = 1)


# ELCA
ELCA_loc <- ELCA %>% 
  filter(!is.na(decimalLatitude), !is.na(decimalLongitude)) %>% 
  mutate(bin = cut(decimalLongitude, breaks = 100)) %>% 
  group_by(bin) %>% 
  summarize(mean_month = mean(month, na.rm = TRUE),
            mean_long = mean(decimalLongitude, na.rm = TRUE),
            samplesize = n())
ELCA_loc

ELCA_dates <- ggplot(data = ELCA_loc)+
  geom_point(aes(x = mean_long, y = mean_month, size = samplesize))+
  geom_smooth(aes(x = mean_long, y = mean_month), method = "glm")
ELCA_dates

ELCA_state <- ELCA %>%
  group_by(stateProvince) %>% 
  summarize(mean_month = mean(month, na.rm = TRUE),
            samplesize = n())
level_order <- c("New Mexico", "Texas", "TEXAS", "Oklahoma", "Arkansas", "Louisiana", "Mississippi","Tennessee", "Kentucky", "Ohio", "Alabama", "Georgia", "Florida", "North Carolina", "South Carolina", "Virginia", "Maryland")

ELCA_dates_states <- ggplot(data = ELCA_state)+
  geom_point(aes(x = factor(stateProvince, levels = level_order), y = mean_month, size = samplesize))
ELCA_dates_states

# I'm gonna try to make an interactive county map

ELCA_county<- ELCA %>%
  group_by(stateProvince, county) %>% 
  summarize(mean_month = mean(month, na.rm = TRUE),
            samplesize = n())

USA <- getData("GADM", country = "usa", level = 2)
temp <- merge(USA, ELCA_county,
              by.x = c("NAME_1", "NAME_2"), by.y = c("stateProvince", "county"),
              all.x = TRUE)
# create a color palette
mypal <- colorNumeric(palette = "viridis", domain = temp$mean_month, na.color = "grey")

leaflet() %>% 
  addProviderTiles("OpenStreetMap.Mapnik") %>%
  setView(lat = 39.8283, lng = -98.5795, zoom = 4) %>% 
  addPolygons(data = temp, stroke = FALSE, smoothFactor = 0.2, fillOpacity = 0.8,
              fillColor = ~mypal(temp$mean_month),
              popup = paste("Region: ", temp$NAME_2, "<br>",
                            "Mean Month: ", temp$mean_month, "<br>",
                            "Samples: ", temp$samplesize, "<br>")) %>%
  addLegend(position = "bottomleft", pal = mypal, values = temp$mean_month,
            title = "Mean Month",
            opacity = 1)


# Now I'm gonna make a map with all the species together

allspp <- POAU %>% 
  rbind(SPNI) %>% 
  rbind(AGHY) %>% 
  rbind(AGPE) %>% 
  rbind(ELVI) %>% 
  rbind(ELCA)

all_county<- allspp %>%
  group_by(stateProvince, county) %>% 
  summarize(mean_month = mean(month, na.rm = TRUE),
            species = length(unique(specificEpithet)),
            samplesize = n())

USA <- getData("GADM", country = "usa", level = 2)
temp <- merge(USA, all_county,
              by.x = c("NAME_1", "NAME_2"), by.y = c("stateProvince", "county"),
              all.x = TRUE)
# create a color palette
mypal <- colorNumeric(palette = "viridis", domain = temp$species, na.color = "grey")

leaflet() %>% 
  addProviderTiles("OpenStreetMap.Mapnik") %>%
  setView(lat = 39.8283, lng = -98.5795, zoom = 4) %>% 
  addPolygons(data = temp, stroke = FALSE, smoothFactor = 0.2, fillOpacity = 0.8,
              fillColor = ~mypal(temp$species),
              popup = paste("Region: ", temp$NAME_2, "<br>",
                            "Mean Month: ", temp$mean_month, "<br>",
                            "Samples: ", temp$samplesize, "<br>",
                            "# of Species", temp$species, "<br>")) %>%
  addLegend(position = "bottomleft", pal = mypal, values = temp$mean_month,
            title = "Mean Month",
            opacity = 1)


# make a table of mean month by species
all_county1<- allspp %>%
  group_by(stateProvince, specificEpithet) %>% 
  summarize(mean_month = mean(month, na.rm = TRUE),
            samplesize = n())
table(all_county1$stateProvince, all_county1$specificEpithet,  all_county1$mean_month)

View(all_county1)
all_county2<- allspp %>%
  group_by(specificEpithet) %>% 
  summarize(mean_month = mean(month, na.rm = TRUE),
            samplesize = n())
table(all_county2$mean_month, all_county2$specificEpithet)
unique(all_county2$specificEpithet)
View(all_county2)
