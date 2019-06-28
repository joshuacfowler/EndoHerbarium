library(ggplot2)
library(dplyr)
library(stringr)
library(readxl)
library(measurements)
library(sp)
library(tidyr)
library(reshape2)
library(tidyverse)

world_endo <- read_xlsx("~joshuacfowler/Downloads/Rudgers_endophyte_world_database2017.xlsx")
endosurveys <- read.csv("~joshuacfowler/Downloads/appendixB.csv")
sneck <- read.csv("~joshuacfowler/Downloads/sneckelvielca.csv")
ecolab_sites <- read_excel(path = "~joshuacfowler/Dropbox/Student Projects/Ceyda Kural/ecolab.xlsx")
ecolab_surveys <- read_excel(path ="~joshuacfowler/Dropbox/Student Projects/Ceyda Kural/ecolab.xlsx", sheet = "Plants Scored"  )

ecolab_sites_gpsfix <- ecolab_sites %>% 
  mutate(gps_lat = case_when(.$siteID == "GBP15" ~ 29.726, .$siteID == "KP15" ~  29.7220, .$siteID == "KP18"~ 29.7220, .$siteID == "LP15" ~ 30.0910, TRUE ~.$gps_lat )) %>% 
  mutate(gps_lon = case_when(.$siteID == "GBP15" ~ -95.687,.$siteID == "KP15" ~ -95.4180, .$siteID == "KP18"~ -95.4180, .$siteID == "LP15" ~  -97.4991, TRUE ~.$gps_lon ))

ecolab <- ecolab_sites_gpsfix %>% 
  merge(ecolab_surveys, by= "siteID") %>% 
  filter(s.scored>0, !is.na(s.eplus)) %>% 
  mutate(endostat = case_when(.$s.eplus==0 ~ 0, .$s.eplus >0 ~ 1)) %>% 
  group_by(siteID,speciesID, gps_lat, gps_lon) %>% 
  summarize(popprev = mean(endostat))


ecolab <- read_csv("~joshuacfowler/Downloads/endosurveys_withPRISM.csv") %>% 
  rename(gps_lon = Long, gps_lat = Lat, popprev = endoprev) %>% 
  mutate(speciesID = case_when(spp == "Agrostis_hyemalis" ~ "AGHY",
                               spp == "Poa_autumnalis" ~ "POAU",
                               spp == "Poa_sylvestris" ~ "POSY",
                               spp == "Agrostis_perennans" ~ "AGPE",
                               spp == "Elymus_virginicus" ~ "ELVI",
                               spp == "Elymus_canadensis" ~ "ELCA",
                               spp == "Sphenopholis_nitida" ~ "SPNI",
                               spp == "Sphenopholis_obtusata" ~ "SPOB",
                               spp == "Cinna_arundinacea" ~ "CIAR",
                               spp == "Lolium_perenne" ~ "LOPE",
                               spp == "Festuca_subverticillata" ~ "FESU"
                               ))
# View(ecolab)
ggplot(ecolab) +
  geom_point(aes(x = gps_lon, y = popprev))+ facet_wrap(~spp)

ggplot(ecolab) +
  geom_point(aes(x = gps_lat, y = popprev))+ facet_wrap(~speciesID)

ggplot(ecolab) +
  geom_point(aes(x = PRISMppt30y, y = popprev))+ facet_wrap(~speciesID)

ggplot(ecolab) +
  geom_point(aes(x = PRISMtmean3, y = popprev))+ facet_wrap(~speciesID)


# Longitude

ggplot(subset(ecolab, speciesID %in% "AGHY")) + geom_vline(xintercept =-98) +
  geom_smooth(aes(x = gps_lon, y = popprev), method = "lm") +
  geom_point(aes(x= gps_lon, y = popprev), lwd = 3) + ggtitle("AGHY") 

ggplot(subset(ecolab, speciesID %in% "AGPE")) + geom_vline(xintercept =-98) +
  geom_smooth(aes(x = gps_lon, y = popprev), method = "lm") +
  geom_point(aes(x= gps_lon, y = popprev), lwd = 3) + ggtitle("AGPE") 

ggplot(subset(ecolab, speciesID %in% "POAU"))+ 
  geom_smooth(aes(x = gps_lon, y = popprev), method = "lm")+
  geom_point(aes(x= gps_lon, y = popprev), lwd = 3) + ggtitle("POAU") 

ggplot(subset(ecolab, speciesID %in% "POSY"))+ 
  geom_smooth(aes(x = gps_lon, y = popprev), method = "lm")+
  geom_point(aes(x= gps_lon, y = popprev), lwd = 3) + ggtitle("POSY") 

ggplot(subset(ecolab, speciesID %in% "ELVI"))+ 
  geom_smooth(aes(x = gps_lon, y = popprev), method = "lm")+
  geom_point(aes(x= gps_lon, y = popprev), lwd = 3) + ggtitle("ELVI") 

ggplot(subset(ecolab, speciesID %in% "ELCA"))+ 
  geom_smooth(aes(x = gps_lon, y = popprev), method = "lm")+
  geom_point(aes(x= gps_lon, y = popprev), lwd = 3) + ggtitle("ELCA") 

ggplot(subset(ecolab, speciesID %in% "CIAR"))+ 
  geom_smooth(aes(x = gps_lon, y = popprev), method = "lm")+
  geom_point(aes(x= gps_lon, y = popprev), lwd = 3) + ggtitle("CIAR") 

ggplot(subset(ecolab, speciesID %in% "SPOB"))+ 
  geom_smooth(aes(x = gps_lon, y = popprev), method = "lm")+
  geom_point(aes(x= gps_lon, y = popprev), lwd = 3) + ggtitle("SPOB") 

ggplot(subset(ecolab, speciesID %in% "FESU"))+ 
  geom_smooth(aes(x = gps_lon, y = popprev), method = "lm")+
  geom_point(aes(x= gps_lon, y = popprev), lwd = 3) + ggtitle("FESU") 


# Latitude
ggplot(subset(ecolab, speciesID %in% "AGHY")) +
  geom_smooth(aes(x = gps_lat, y = popprev), method = "lm") +
  geom_point(aes(x= gps_lat, y = popprev), lwd = 3) + ggtitle("AGHY") 

ggplot(subset(ecolab, speciesID %in% "AGPE")) +
  geom_smooth(aes(x = gps_lat, y = popprev), method = "lm") +
  geom_point(aes(x= gps_lat, y = popprev), lwd = 3) + ggtitle("AGPE") 

ggplot(subset(ecolab, speciesID %in% "POAU"))+ 
  geom_smooth(aes(x = gps_lat, y = popprev), method = "lm")+
  geom_point(aes(x= gps_lat, y = popprev), lwd = 3) + ggtitle("POAU") 

ggplot(subset(ecolab, speciesID %in% "POSY"))+ 
  geom_smooth(aes(x = gps_lat, y = popprev), method = "lm")+
  geom_point(aes(x= gps_lat, y = popprev), lwd = 3) + ggtitle("POSY") 

ggplot(subset(ecolab, speciesID %in% "ELVI"))+ 
  geom_smooth(aes(x = gps_lat, y = popprev), method = "lm")+
  geom_point(aes(x= gps_lat, y = popprev), lwd = 3) + ggtitle("ELVI") 

ggplot(subset(ecolab, speciesID %in% "ELCA"))+ 
  geom_smooth(aes(x = gps_lat, y = popprev), method = "lm")+
  geom_point(aes(x= gps_lat, y = popprev), lwd = 3) + ggtitle("ELCA") 

ggplot(subset(ecolab, speciesID %in% "CIAR"))+ 
  geom_smooth(aes(x = gps_lat, y = popprev), method = "lm")+
  geom_point(aes(x= gps_lat, y = popprev), lwd = 3) + ggtitle("CIAR") 

ggplot(subset(ecolab, speciesID %in% "SPOB"))+ 
  geom_smooth(aes(x = gps_lat, y = popprev), method = "lm")+
  geom_point(aes(x= gps_lat, y = popprev), lwd = 3) + ggtitle("SPOB") 

ggplot(subset(ecolab, speciesID %in% "FESU"))+ 
  geom_smooth(aes(x = gps_lat, y = popprev), method = "lm")+
  geom_point(aes(x= gps_lat, y = popprev), lwd = 3) + ggtitle("FESU") 












sneckspp <- sneck %>% 
  mutate(sppcode = case_when(.$Species ==  "Both" ~ "ELVI;ELCA", 
                             .$Species == "EV" ~ "ELVI",
                             .$Species == "EC" ~ "ELCA"))

snecksplit <- sneckspp %>% 
  separate(E....., c("Eprev1","Eprev2"), sep = ";") %>% 
  select(Site, Site.name, Species, Lat., Long., Temp..., Ppt...mm., SPEI, Elev.m., Sampling.Date,Eprev1, Eprev2, sppcode) %>% 
  melt(id.vars = c("Site", "Site.name", "Species", "Lat.", "Long.", "Temp...", "Ppt...mm.", "SPEI", "Elev.m.", "Sampling.Date", "sppcode")) %>% 
  mutate(endoprev = case_when(.$sppcode == "ELVI" & .$variable == "Eprev1" ~ .$value,.$sppcode == "ELCA" & .$variable == "Eprev1" ~ .$value,.$sppcode == "ELVI;ELCA" & .$variable == "Eprev1" ~ .$value,.$sppcode == "ELVI;ELCA" & .$variable == "Eprev2" ~ .$value)) %>% 
  mutate(spp = case_when(.$sppcode == "ELVI" ~ "ELVI", .$sppcode == "ELCA" ~ "ELCA", .$sppcode == "ELVI;ELCA" & .$variable == "Eprev1" ~ "ELCA", .$sppcode == "ELVI;ELCA" & .$variable == "Eprev2"~ "ELVI")) %>% 
  mutate(endoprev = as.numeric(endoprev)*.01) %>% 
  filter(!is.na(endoprev))

         
# View(snecksplit)

ggplot(snecksplit)+
  geom_point(aes(x = Long., y = Lat., color = endoprev), lwd = 3) +facet_wrap(~spp)+ scale_colour_gradient(low = "lightblue", high = "blue")

ggplot(snecksplit)+
  geom_point(aes(x = Long., y = endoprev), lwd = 3) +facet_wrap(~spp)+ scale_colour_gradient(low = "lightblue", high = "blue")


ggplot(subset(snecksplit, spp %in% "ELVI"))+
  geom_point(aes(x = Long., y = Lat., lwd = endoprev, color = Ppt...mm.))+ scale_colour_gradient(low = "lightblue", high = "blue")
ggplot(subset(snecksplit, spp %in% "ELVI"))+
  geom_point(aes(x = Long., y = Lat., lwd = endoprev, color = Temp...))+ scale_colour_gradient(low = "lightblue", high = "blue")
ggplot(subset(snecksplit, spp %in% "ELVI"))+
  geom_point(aes(x = Long., y = Lat., lwd = endoprev, color = SPEI))+ scale_colour_gradient(low = "lightblue", high = "blue")



ggplot(subset(snecksplit, spp %in% "ELCA"))+
  geom_point(aes(x = Long., y = Lat., lwd = endoprev, color = Ppt...mm.))+ scale_colour_gradient(low = "lightblue", high = "blue")




ggplot(subset(snecksplit, spp %in% "ELCA"))+ geom_vline(xintercept = -94) +
  geom_point(aes(x = Long., y = Lat., color = endoprev), lwd = 3)+ scale_colour_gradient(low = "lightblue", high = "blue") + ggtitle("ELCA")

ggplot(subset(snecksplit, spp %in% "ELCA"))+ geom_vline(xintercept = -94) +
  geom_smooth(aes(x = Long., y = endoprev), method = "lm")+
  geom_point(aes(x = Long., y = endoprev), lwd = 3)+ scale_colour_gradient(low = "lightblue", high = "blue") + ggtitle("ELCA")
ggplot(subset(snecksplit, spp %in% "ELCA"))+
  geom_smooth(aes(x = Lat., y = endoprev), method = "lm")+
  geom_point(aes(x = Lat., y = endoprev), lwd = 3)+ scale_colour_gradient(low = "lightblue", high = "blue") + ggtitle("ELCA")


ggplot(subset(snecksplit, spp %in% "ELVI"))+ geom_vline(xintercept = -101)+
  geom_point(aes(x = Long., y = Lat., color = endoprev), lwd = 3)+ scale_colour_gradient(low = "lightblue", high = "blue") + ggtitle("ELVI")

ggplot(subset(snecksplit, spp %in% "ELVI"))+
  geom_smooth(aes(x = Lat., y = endoprev), method = "lm")+
  geom_point(aes(x= Lat., y = endoprev), lwd = 3)+ scale_colour_gradient(low = "lightblue", high = "blue") + ggtitle("ELVI") 

ggplot(subset(snecksplit, spp %in% "ELVI"))+
  geom_smooth(aes(x = Long., y = endoprev), method = "lm")+
  geom_point(aes(x= Long., y = endoprev), lwd = 3)+ scale_colour_gradient(low = "lightblue", high = "blue") + ggtitle("ELVI") 


ggplot(subset(snecksplit, spp %in% "ELVI")) + geom_vline(xintercept= -101)+
  geom_smooth(aes(x = Long., y = SPEI), method = "lm")+
  geom_point(aes(x = Long., y = SPEI), lwd = 3)+ ggtitle("ELVI") 

# SPEI (positive = wetter)
ggplot(subset(snecksplit, spp %in% "ELVI"))+
  geom_smooth(aes(x = SPEI, y = endoprev), method = "lm")+
  geom_point(aes(x = SPEI, y = endoprev), lwd = 3)+ ggtitle("ELVI")


ggplot(subset(snecksplit, spp %in% "ELVI")) + geom_vline(xintercept= -101)+
  geom_smooth(aes(x = Long., y = Ppt...mm.), method = "lm")+
  geom_point(aes(x = Long., y = Ppt...mm.), lwd = 3)+ ggtitle("ELVI") 

ggplot(subset(snecksplit, spp %in% "ELVI")) + geom_vline(xintercept= -101)+
  geom_smooth(aes(x = Long., y = Temp...), method = "lm")+
  geom_point(aes(x = Long., y = Temp...), lwd = 3)+ ggtitle("ELVI") 


# Jenn Rudger's 2009 data

endos <- endosurveys %>% 
  mutate(Lat = as.factor(Latitude), Long = as.factor(Longitude?)) %>% 
  mutate(plantsppcode = str_c(Plant.genus, Plant.species, sep = "_"))


# coming frequency measurements
endos$endofrequency = ifelse(!is.na(endos$Tiller.Freq....) & is.na(endos$Seed.Freq....), endos$Tiller.Freq..., 
                             ifelse(is.na(endos$Tiller.Freq....) & !is.na(endos$Seed.Freq....), endos$Seed.Freq...., 
                                    ifelse(!is.na(endos$Tiller.Freq....) & !is.na(endos$Seed.Freq....), (endos$Tiller.Freq.... + endos$Seed.Freq....)/2, NA)))
# Converting lat long coordinates
endos$Lat = gsub(' ', "", endos$Lat)
endos$Lat = gsub('N', "", endos$Lat)
endos$Lat = gsub('?', " ", endos$Lat)
endos$Lat = gsub('"', "", endos$Lat)
endos$Lat = gsub("'", " ", endos$Lat)

endos$Long = gsub(' ', "", endos$Long)
endos$Long = gsub('W', "", endos$Long)
endos$Long = gsub('?', " ", endos$Long)
endos$Long = gsub('"', "", endos$Long)
endos$Long = gsub("'", " ", endos$Long)
endos$Long = gsub("[:?:]", "", endos$Long)

endos1 <- endos %>% 
  separate(Lat, paste("lat",c("d","m","s"), sep="_") ) %>%
  separate(Long, paste("long",c("d","m","s"), sep="_" ) ) %>%
  mutate(lat_s = case_when(.$lat_s == "" ~ "00", !is.na(.$lat_s)~.$lat_s)) %>% 
  mutate(long_s = case_when(.$long_s == "" ~"00", !is.na(.$long_s)~.$long_s))  %>% 
  mutate_each(funs(as.numeric), lat_d,lat_m,lat_s, long_d, long_m, long_s) %>% 
  mutate(lat_dec=lat_d + lat_m/60 + lat_s/60^2,
            long_dec=long_d + long_m/60 + long_s/60^2)




# View(endos1)
ggplot(data = endos1) +
  geom_point(aes(x = long_dec, y = lat_dec, color = endofrequency)) + facet_wrap(~plantsppcode) + scale_x_reverse() + scale_colour_gradient(low = "lightblue", high = "blue")

ggplot(subset(endos1, plantsppcode %in% c("Poa_autumnalis"))) +
  geom_point(aes(x = long_dec, y = lat_dec, color = endofrequency), lwd = 3) + geom_vline(xintercept = 96) + scale_x_reverse()+ scale_colour_gradient(low = "lightblue", high = "blue") + ggtitle("POAU")

ggplot(subset(endos1, plantsppcode %in% c("Poa_autumnalis"))) +
  geom_smooth(aes(x = long_dec, y = endofrequency), method = "lm")+
  geom_point(aes(x = long_dec, y = endofrequency), lwd = 3) + geom_vline(xintercept = 96) + scale_x_reverse()+ scale_colour_gradient(low = "lightblue", high = "blue") + ggtitle("POAU")


ggplot(subset(endos1, plantsppcode %in% c("Poa_sylvestris"))) +
  geom_point(aes(x = long_dec, y = lat_dec, color = endofrequency), lwd = 3) + scale_x_reverse()+ scale_colour_gradient(low = "lightblue", high = "blue")

ggplot(subset(endos1, plantsppcode %in% c("Poa_sylvestris"))) +
  geom_smooth(aes(x = long_dec, y = endofrequency), method = "lm")+
  geom_point(aes(x = long_dec, y = endofrequency), lwd = 3) + scale_x_reverse()+ scale_colour_gradient(low = "lightblue", high = "blue")


ggplot(subset(endos1, plantsppcode %in% c("Festuca_subverticillata"))) +
  geom_point(aes(x = long_dec, y = lat_dec, color = endofrequency), lwd = 3) + scale_x_reverse()+ scale_colour_gradient(low = "lightblue", high = "blue")

ggplot(subset(endos1, plantsppcode %in% c("Festuca_subverticillata"))) +
  geom_smooth(aes(x = long_dec, y = endofrequency), method = "lm")+
  geom_point(aes(x = long_dec, y = lat_dec, color = endofrequency), lwd = 3) + scale_x_reverse()+ scale_colour_gradient(low = "lightblue", high = "blue")


ggplot(subset(endos1, plantsppcode %in% c("Elymus_virginicus"))) + ggtitle("ELVI")+
  geom_point(aes(x = long_dec, y = lat_dec, color = endofrequency), lwd = 3) + geom_vline(xintercept = 100) + scale_x_reverse()+ scale_colour_gradient(low = "lightblue", high = "blue")

ggplot(subset(endos1, plantsppcode %in% c("Elymus_virginicus"))) + ggtitle("ELVI")+
  geom_smooth(aes(x = long_dec, y = endofrequency), method = "lm")+
  geom_point(aes(x = long_dec, y = endofrequency), lwd = 3) + geom_vline(xintercept = 100) + scale_x_reverse()+ scale_colour_gradient(low = "lightblue", high = "blue")


ggplot(subset(endos1, plantsppcode %in% c("Agrostis_hyemalis"))) + geom_vline(xintercept =98) +
  geom_point(aes(x = long_dec, y = lat_dec, color = endofrequency), lwd = 3) + scale_x_reverse()+ scale_colour_gradient(low = "lightblue", high = "blue") + ggtitle("AGHY")

ggplot(subset(endos1, plantsppcode %in% c("Agrostis_hyemalis"))) + geom_vline(xintercept =98) +
  geom_smooth(aes(x = long_dec, y = endofrequency), method = "lm")+
  geom_point(aes(x = long_dec, y = endofrequency), lwd = 3) + scale_x_reverse()+ scale_colour_gradient(low = "lightblue", high = "blue") + ggtitle("AGHY")


ggplot(subset(endos1, plantsppcode %in% c("Agrostis_perennans"))) + ggtitle("AGPE") + geom_vline(xintercept = 95)+
  geom_point(aes(x = long_dec, y = lat_dec, color = endofrequency), lwd = 3) + scale_x_reverse()+ scale_colour_gradient(low = "lightblue", high = "blue")

ggplot(subset(endos1, plantsppcode %in% c("Agrostis_perennans"))) + ggtitle("AGPE") + geom_vline(xintercept = 95)+
  geom_smooth(aes(x = long_dec, y = endofrequency), method = "lm")+
  geom_point(aes(x = long_dec, y = endofrequency), lwd = 3) + scale_x_reverse()+ scale_colour_gradient(low = "lightblue", high = "blue")

ggplot(subset(endos1, plantsppcode %in% c("Agrostis_perennans"))) + ggtitle("AGPE") + geom_vline(xintercept = 95)+
  geom_smooth(aes(x = lat_dec, y = endofrequency), method = "lm")+
  geom_point(aes(x = lat_dec, y = endofrequency), lwd = 3) + scale_x_reverse()+ scale_colour_gradient(low = "lightblue", high = "blue")


ggplot(subset(endos1, plantsppcode %in% c("Sphenopholis_nitida"))) + ggtitle("SPNI") + geom_vline(xintercept = 95)+
  geom_point(aes(x = long_dec, y = lat_dec, color = endofrequency), lwd = 3) + scale_x_reverse()+ scale_colour_gradient(low = "lightblue", high = "blue")

ggplot(subset(endos1, plantsppcode %in% c("Sphenopholis_nitida"))) + ggtitle("SPNI") + geom_vline(xintercept = 95)+
  geom_smooth(aes(x = long_dec, y = endofrequency), method = "lm")+
  geom_point(aes(x = long_dec, y = endofrequency), lwd = 3) + scale_x_reverse()+ scale_colour_gradient(low = "lightblue", high = "blue")


ggplot(subset(endos1, plantsppcode %in% c("Brachyelytrum_erectum"))) +
  geom_point(aes(x = long_dec, y = lat_dec, color = endofrequency), lwd = 3) + scale_x_reverse()+ scale_colour_gradient(low = "lightblue", high = "blue")

ggplot(subset(endos1, plantsppcode %in% c("Brachyelytrum_erectum"))) +
  geom_smooth(aes(x = long_dec, y = endofrequency), method = "lm")+
  geom_point(aes(x = long_dec, y = endofrequency), lwd = 3) + scale_x_reverse()+ scale_colour_gradient(low = "lightblue", high = "blue")


ggplot(subset(endos1, plantsppcode %in% c("Cinna_arundinacea"))) +
  geom_point(aes(x = long_dec, y = lat_dec, color = endofrequency), lwd = 3) + scale_x_reverse()+ scale_colour_gradient(low = "lightblue", high = "blue")

ggplot(subset(endos1, plantsppcode %in% c("Cinna_arundinacea"))) +
  geom_smooth(aes(x = long_dec, y = endofrequency), method = "lm")+
  geom_point(aes(x = long_dec, y = endofrequency), lwd = 3) + scale_x_reverse()+ scale_colour_gradient(low = "lightblue", high = "blue")


ggplot(data = endos1) +
  geom_point(aes(x = long_dec, y = lat_dec, color = plantsppcode)) + scale_x_reverse()



# Merging the differnt surveys together
ecolab1 <- ecolab %>% 
  select(gps_lat, gps_lon, speciesID, siteID, popprev) %>% 
  mutate(popprev = as.numeric(popprev)) %>% 
  rename(Lat = gps_lat, Long = gps_lon, spp = speciesID, LocID = siteID, endoprev = popprev) %>% 
  mutate(data_origin = "ecolab")

snecksplit1 <- snecksplit %>% 
  select(Lat., Long., spp, Site, endoprev) %>% 
  rename(Lat = Lat., Long = Long., spp = spp, LocID = Site, endoprev = endoprev) %>%
  mutate(data_origin = "sneck")

endos2 <- endos1 %>% 
  select(lat_dec, long_dec, plantsppcode, No., endofrequency) %>% 
  mutate(long_dec = -(long_dec)) %>% 
  rename(Lat = lat_dec, Long = long_dec, spp = plantsppcode, LocID = No., endoprev = endofrequency) %>%
  mutate(data_origin = "rudgers2009") %>% 
  filter(!is.na(Long)) %>% 
  filter(!is.na(endoprev)) %>% 
  mutate(endoprev = .$endoprev*.01)
  

surveys <- ecolab1 %>% 
  merge(snecksplit1, by = c("Lat", "Long","spp","LocID", "endoprev", "data_origin"), all = TRUE) %>% 
  merge(endos2, by = c("Lat", "Long","spp","LocID", "endoprev", "data_origin"), all = TRUE) %>% 
  mutate(spp = case_when(.$spp == "ELVI" ~ "Elymus_virginicus", .$spp == "ELCA" ~ "Elymus_canadensis", .$spp == "AGHY" ~ "Agrostis_hyemalis", .$spp == "POAU" ~ "Poa_autumnalis", .$spp == "LOLI" ~ "Lolium_perenne", TRUE ~ as.character(.$spp)))
write.csv(surveys, file = "endosurveys")
# View(surveys)
# write.csv(surveys, file = "")
# I used the output file from endosurveys and used QGIS to extract the PRISM data (PPT and Temp) for the GPS locations
endosurveys_withPRISM <- read.csv("~joshuacfowler/Desktop/endosurveys_withPRISM.csv")

endoPRISM <- endosurveys_withPRISM %>% 
  rename(PPT = PRISMppt30y, Temp = PRISMtmean3) %>% 
  rename(data_origin = data_origi) #I'm not sure but it cut off the label



ggplot(endoPRISM)+
  geom_point(aes(x = Long, y = endoprev, color = data_origin))+facet_wrap(~spp)

ggplot(endoPRISM)+
  geom_point(aes(x = Long, y = endoprev, color = data_origin))
ggplot(endoPRISM)+
  geom_point(aes(x = Lat, y = endoprev, color = data_origin))

ggplot(endoPRISM)+
  geom_point(aes(x = Long, y = Lat, color = data_origin), position = position_jitter(width = .05, height = .05)) + ggtitle("Locations and data sources of Endophyte surveys")


ggplot(subset(endoPRISM, spp %in% "Poa_autumnalis")) +
  geom_smooth(aes(x = Long, y = endoprev), method = "glm")+
  geom_point(aes(x = Long, y = endoprev, color = data_origin))+ ggtitle("Poa_autumnalis")
  

ggplot(subset(endoPRISM, spp %in% "Poa_autumnalis")) +
  geom_smooth(aes(x = Lat, y = endoprev), method = "glm")+
  geom_point(aes(x = Lat, y = endoprev, color = data_origin))+ ggtitle("Poa_autumnalis")


ggplot(subset(endoPRISM, spp %in% "Elymus_virginicus")) +
  geom_smooth(aes(x = Long, y = endoprev), method = "glm")+
  geom_point(aes(x = Long, y = endoprev, color = data_origin))+ ggtitle("Elymus_virginicus")

ggplot(subset(endoPRISM, spp %in% "Elymus_virginicus")) +
  geom_smooth(aes(x = Lat, y = endoprev), method = "glm")+
  geom_point(aes(x = Lat, y = endoprev, color = data_origin))+ ggtitle("Elymus_virginicus")


ggplot(subset(endoPRISM, spp %in% "Elymus_canadensis")) +
  geom_smooth(aes(x = Long, y = endoprev), method = "glm")+
  geom_point(aes(x = Long, y = endoprev, color = data_origin))+ ggtitle("Elymus_canadensis")

ggplot(subset(endoPRISM, spp %in% "Elymus_canadensis")) +
  geom_smooth(aes(x = Lat, y = endoprev), method = "glm")+
  geom_point(aes(x = Lat, y = endoprev, color = data_origin))+ ggtitle("Elymus_canadensis")

ggplot(subset(endoPRISM, spp %in% "Agrostis_hyemalis")) +
  geom_smooth(aes(x = Long, y = endoprev), method = "glm")+
  geom_point(aes(x = Long, y = endoprev, color = data_origin))+ ggtitle("Agrostis_hyemalis")

ggplot(subset(endoPRISM, spp %in% "Agrostis_hyemalis")) +
  geom_smooth(aes(x = Lat, y = endoprev), method = "glm")+
  geom_point(aes(x = Lat, y = endoprev, color = data_origin))+ ggtitle("Agrostis_hyemalis")

ggplot(subset(endoPRISM, spp %in% "Agrostis_perennans")) +
  geom_smooth(aes(x = Long, y = endoprev), method = "glm")+
  geom_point(aes(x = Long, y = endoprev, color = data_origin))+ ggtitle("Agrostis_perennans")

ggplot(subset(endoPRISM, spp %in% "Agrostis_perennans")) +
  geom_smooth(aes(x = Lat, y = endoprev), method = "glm")+
  geom_point(aes(x = Lat, y = endoprev, color = data_origin))+ ggtitle("Agrostis_perennans")

ggplot(subset(endoPRISM, spp %in% "Sphenopholis_nitida")) +
  geom_smooth(aes(x = Long, y = endoprev), method = "glm")+
  geom_point(aes(x = Long, y = endoprev, color = data_origin))+ ggtitle("Sphenopholis_nitida")

ggplot(subset(endoPRISM, spp %in% "Sphenopholis_nitida")) +
  geom_smooth(aes(x = Lat, y = endoprev), method = "glm")+
  geom_point(aes(x = Lat, y = endoprev, color = data_origin))+ ggtitle("Sphenopholis_nitida")





ggplot(endoPRISM)+
  geom_point(aes(x = Long, y = PPT, color = data_origin))
ggplot(endoPRISM)+
  geom_point(aes(x = Lat, y = PPT, color = data_origin))
ggplot(endoPRISM)+
  geom_point(aes(x = Long, y = Temp, color = data_origin))
ggplot(endoPRISM)+
  geom_point(aes(x = Lat, y = Temp, color = data_origin))

ggplot(endoPRISM)+
  geom_point(aes(x = Long, y = Lat, shape = data_origin, color = PPT))

ggplot(endoPRISM)+
  geom_point(aes(x = Long, y = Lat, shape = data_origin, color = Temp))



ggplot(subset(endoPRISM, spp %in% "Poa_autumnalis")) +
  geom_smooth(aes(x = PPT, y = endoprev), method = "glm")+
  geom_point(aes(x = PPT, y = endoprev, color = data_origin))+ ggtitle("Poa_autumnalis")


ggplot(subset(endoPRISM, spp %in% "Poa_autumnalis")) +
  geom_smooth(aes(x = Temp, y = endoprev), method = "glm")+
  geom_point(aes(x = Temp, y = endoprev, color = data_origin))+ ggtitle("Poa_autumnalis")


ggplot(subset(endoPRISM, spp %in% "Elymus_virginicus")) +
  geom_smooth(aes(x = PPT, y = endoprev), method = "glm")+
  geom_point(aes(x = PPT, y = endoprev, color = data_origin))+ ggtitle("Elymus_virginicus")

ggplot(subset(endoPRISM, spp %in% "Elymus_virginicus")) +
  geom_smooth(aes(x = Temp, y = endoprev), method = "glm")+
  geom_point(aes(x = Temp, y = endoprev, color = data_origin))+ ggtitle("Elymus_virginicus")


ggplot(subset(endoPRISM, spp %in% "Elymus_canadensis")) +
  geom_smooth(aes(x = PPT, y = endoprev), method = "glm")+
  geom_point(aes(x = PPT, y = endoprev, color = data_origin))+ ggtitle("Elymus_canadensis")

ggplot(subset(endoPRISM, spp %in% "Elymus_canadensis")) +
  geom_smooth(aes(x = Temp, y = endoprev), method = "glm")+
  geom_point(aes(x = Temp, y = endoprev, color = data_origin))+ ggtitle("Elymus_canadensis")

ggplot(subset(endoPRISM, spp %in% "Agrostis_hyemalis")) +
  geom_smooth(aes(x = PPT, y = endoprev), method = "glm")+
  geom_point(aes(x = PPT, y = endoprev, color = data_origin))+ ggtitle("Agrostis_hyemalis")

ggplot(subset(endoPRISM, spp %in% "Agrostis_hyemalis")) +
  geom_smooth(aes(x = Temp, y = endoprev), method = "glm")+
  geom_point(aes(x = Temp, y = endoprev, color = data_origin))+ ggtitle("Agrostis_hyemalis")

ggplot(subset(endoPRISM, spp %in% "Agrostis_perennans")) +
  geom_smooth(aes(x = PPT, y = endoprev), method = "glm")+
  geom_point(aes(x = PPT, y = endoprev, color = data_origin))+ ggtitle("Agrostis_perennans")

ggplot(subset(endoPRISM, spp %in% "Agrostis_perennans")) +
  geom_smooth(aes(x = Temp, y = endoprev), method = "glm")+
  geom_point(aes(x = Temp, y = endoprev, color = data_origin))+ ggtitle("Agrostis_perennans")

ggplot(subset(endoPRISM, spp %in% "Sphenopholis_nitida")) +
  geom_smooth(aes(x = PPT, y = endoprev), method = "glm")+
  geom_point(aes(x = PPT, y = endoprev, color = data_origin))+ ggtitle("Sphenopholis_nitida")

ggplot(subset(endoPRISM, spp %in% "Sphenopholis_nitida")) +
  geom_smooth(aes(x = Temp, y = endoprev), method = "glm")+
  geom_point(aes(x = Temp, y = endoprev, color = data_origin))+ ggtitle("Sphenopholis_nitida")



# Jenn Rudger's "world" database


type2_endo <- world_endo %>% 
  mutate(sppname = str_c(plant_genus, plant_species, sep = "_")) %>% 
  filter(sppname == "Agrostis_hyemalis"|sppname == "Agrostis_perennans"|sppname == "Brachelytrum_erectrum"|sppname == "Cinna_arundinacea"|
           sppname == "Elymus_canadensis"|sppname == "Elymus_virginicus"|sppname == "Festuca_paradoxa"|sppname == "Festuca_suverticillata"|
           sppname == "Poa_chapmania"|sppname == "Poa_sylvestris"|sppname == "Poa_autumnalis"|sppname == "Sphenopholis_longiflora"|
           sppname == "Sphenopholis_nitida")
type2_endo1 <- type2_endo %>% 
  group_by(sppname) %>% 
  summarise(meanpopinfect = mean(percentage_of_populations_infected),
            numpop = as.integer(n()))

ggplot(data = type2_endo) +
  geom_line(aes(sppname, mean_percentage_individuals_infected), lwd = 2)+ theme(axis.text.x = element_text(angle = 90, hjust = 1))



ggplot(data = type2_endo1) + ggtitle("Percentage of populations with infection")+
  geom_col(aes(sppname, meanpopinfect, fill = numpop))+ theme(axis.text.x = element_text(angle = 90, hjust = 1))
           

View(type2_endo)
