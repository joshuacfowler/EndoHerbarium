# Species distributions models  Epichloe endophyte host species 
rm(list = ls())
# Load packages
library("prism")
library("raster")
library("dismo")
library("usdm")
library("car")
library("FactoMineR")
library("factoextra")
library("corrplot")
library("ENMTools")
library("vip")
library("pdp")
library("fastshap")
library("CalibratR")
library("maptools")
library("rgeos")
library("leaflet")
library("tidyverse")
# install.extras()
library("geodata")
library("terra")
library("rgdal")
library("sp")
library("ENMeval")
library("mecofun")
library("mgcv")
library("randomForest")
library("spThin")
# if( !("rJava" %in% rownames(installed.packages()))  ){
#   install.packages("rJava",repos="http://cran.r-project.org")
# }

# if(Sys.info()["sysname"] != "Windows" ){
#   dyn.load('/Library/Java/JavaVirtualMachines/temurin-17.jdk/Contents/Home/lib/server/libjvm.dylib')
# }
library("rJava")

set.seed(13)

# Climatic data----
## Data from PRISM---- 
# making a folder to store prism data
prism_set_dl_dir("/Users/jm200/Library/CloudStorage/Dropbox/Miller Lab/Herbaruim Project/Data/Clim")
# prism_set_dl_dir("/Users/jmoutouama/Dropbox/Miller Lab/Herbaruim Project/Data/Clim")

# getting monthly data for mean temp and precipitation
# takes a long time the first time, but can skip when you have raster files saved on your computer.
# get_prism_monthlys(type = "tmean", years = 1895:2020, mon = 1:12, keepZip = FALSE)
# get_prism_monthlys(type = "ppt", years = 1895:2020, mon = 1:12, keepZip = FALSE)

# pulling out values to get normals for old and new time periods
tmean_annual_recent_norm <- terra::mean(terra::rast(pd_stack(prism_archive_subset(type = "tmean", temp_period = "monthly", year = 1990:2020))))
tmean_spring_recent_norm <- terra::mean(terra::rast(pd_stack(prism_archive_subset(type = "tmean", temp_period = "monthly", year = 1990:2020, mon = 1:4))))
tmean_summer_recent_norm <- terra::mean(terra::rast(pd_stack(prism_archive_subset(type = "tmean", temp_period = "monthly", year = 1990:2020, mon = 5:8))))
tmean_autumn_recent_norm <- terra::mean(terra::rast(pd_stack(prism_archive_subset(type = "tmean", temp_period = "monthly", year = 1990:2020, mon = 9:12))))

# calculating standard deviation in temp

tmean_annual_recent_sd <- terra::stdev(terra::rast(pd_stack(prism_archive_subset(type = "tmean", temp_period = "monthly", year = 1990:2020))))
tmean_spring_recent_sd <- terra::stdev(terra::rast(pd_stack(prism_archive_subset(type = "tmean", temp_period = "monthly", year = 1990:2020, mon = 1:4))))
tmean_summer_recent_sd <- terra::stdev(terra::rast(pd_stack(prism_archive_subset(type = "tmean", temp_period = "monthly", year = 1990:2020, mon = 5:8))))
tmean_autumn_recent_sd <- terra::stdev(terra::rast(pd_stack(prism_archive_subset(type = "tmean", temp_period = "monthly", year = 1990:2020, mon = 9:12))))

# calculating the cumulative precipitation for each year and for each season within the year
ppt_annual_recent <- ppt_spring_recent <- ppt_summer_recent <- ppt_autumn_recent <- ppt_winter_recent<- list()
for(y in 1990:2020){
  ppt_annual_recent[[y]] <- sum(terra::rast(pd_stack(prism_archive_subset(type = "ppt", temp_period = "monthly", year = y))))
  ppt_spring_recent[[y]] <- sum(terra::rast(pd_stack(prism_archive_subset(type = "ppt", temp_period = "monthly", year = y, mon = 1:4))))
  ppt_summer_recent[[y]] <- sum(terra::rast(pd_stack(prism_archive_subset(type = "ppt", temp_period = "monthly", year = y, mon = 5:8))))
  ppt_autumn_recent[[y]] <- sum(terra::rast(pd_stack(prism_archive_subset(type = "ppt", temp_period = "monthly", year = y, mon = 9:12)))) 
}


# Taking the mean of the cumulative precipation values
ppt_annual_recent_norm <- terra::mean(terra::rast(unlist(ppt_annual_recent)))
ppt_spring_recent_norm <- terra::mean(terra::rast(unlist(ppt_spring_recent)))
ppt_summer_recent_norm <- terra::mean(terra::rast(unlist(ppt_summer_recent)))
ppt_autumn_recent_norm <- terra::mean(terra::rast(unlist(ppt_autumn_recent)))

#calculating the standard devation in precip
ppt_annual_recent_sd <- terra::stdev(terra::rast(unlist(ppt_annual_recent)))
ppt_spring_recent_sd <- terra::stdev(terra::rast(unlist(ppt_spring_recent)))
ppt_summer_recent_sd <- terra::stdev(terra::rast(unlist(ppt_summer_recent)))
ppt_autumn_recent_sd <- terra::stdev(terra::rast(unlist(ppt_autumn_recent)))


## Variance inflation factor (VIF)  to have  a measure of multicollinearity among the  variables----  
US_worldclim_recent_norm<-terra::rast(list(tmean_spring_recent_norm,tmean_summer_recent_norm,tmean_autumn_recent_norm,
                                           tmean_spring_recent_sd,tmean_summer_recent_sd,tmean_autumn_recent_sd,
                                           ppt_spring_recent_norm,ppt_summer_recent_norm,ppt_autumn_recent_norm,
                                           ppt_spring_recent_sd,ppt_summer_recent_sd,ppt_autumn_recent_sd))

US_worldclim_recent_norm_stack<-stack(US_worldclim_recent_norm)

plot(US_worldclim_recent_norm_stack)
(vif <- vifcor(US_worldclim_recent_norm_stack, th=0.7))


## Saved the rasters 
# writeRaster(tmean_spring_recent_norm, '/Users/jm200/Library/CloudStorage/Dropbox/Miller Lab/Herbaruim Project/Data/Final variable/recent/tmean_spring_recent_norm.tif')
# writeRaster(ppt_spring_recent_norm, '/Users/jm200/Library/CloudStorage/Dropbox/Miller Lab/Herbaruim Project/Data/Final variable/recent/ppt_spring_recent_norm.tif')
# writeRaster(ppt_summer_recent_norm, '/Users/jm200/Library/CloudStorage/Dropbox/Miller Lab/Herbaruim Project/Data/Final variable/recent/ppt_summer_recent_norm.tif')

# Final climatic variables to use for the SDM
US_worldclim_recent_norm_stack_vif<-terra::rast(list(tmean_summer_recent_norm,tmean_spring_recent_sd,tmean_summer_recent_sd,ppt_spring_recent_sd,ppt_summer_recent_sd))
names(US_worldclim_recent_norm_stack_vif) <- c("tmean_summer","tmean_spring_sd", "tmean_summer_sd","ppt_spring_sd","ppt_summer_sd")
US_worldclim_recent_norm_stack_final<-stack(US_worldclim_recent_norm_stack_vif)
plot(US_worldclim_recent_norm_stack_final)

#res(US_worldclim_recent_norm_stack_final[[1]])

# Agrostis hyemalis-----
## Import occurrence data from database and merge them later to existing online data
endo_herb_georef<-read.csv(url("https://www.dropbox.com/scl/fi/lwiv7o156kd11x2phsu93/filtered_endo_herb_georef.csv?rlkey=bxw5gjwq1mdh3sg1vqa7u1zgs&dl=1"), header=T)
# unique(endo_herb_georef$Country)
names(endo_herb_georef)
## Subset only AGHY records that in the US
endo_herb_georef %>% 
  dplyr::select(lon, lat,year,Spp_code,Country)%>% 
  dplyr::filter(Spp_code=="AGHY" & Country=="United States" & year %in% (1990:2020)) %>% 
  dplyr::rename(country=Country)->aghy_occgeoref_1990_2020 

## Download occurrence data from gbif for *Agrostis hyemalis*
# dir.create("/Users/jm200/Library/CloudStorage/Dropbox/Miller Lab/Herbaruim Project/Data/Occurence")
# aghy_occ_raw <- gbif(genus="Agrostis",species="hyemalis",download=TRUE) 
load("/Users/jm200/Library/CloudStorage/Dropbox/Miller Lab/Herbaruim Project/Data/Occurence/aghy_occ_raw.rdata")
# load("/Users/jmoutouama/Dropbox/Miller Lab/Herbaruim Project/Data/Occurence/aghy_occ_raw.rdata")
head(aghy_occ_raw) # to view the first few records the occurrence dataset use
# names(aghy_occ_raw)

# Clean occurrence data from gbif
aghy_occ <- subset(aghy_occ_raw,(!is.na(lat))&(!is.na(lon))&!is.na(country)) # here we remove erroneous coordinates, where either the latitude or longitude is missing
cat(nrow(aghy_occ_raw)-nrow(aghy_occ), "records are removed") # Show the number of records that are removed from the dataset. 
aghy_occ %>% 
  dplyr::select(country,lon, lat, year,basisOfRecord)%>% 
  filter(country=="United States")->aghy_occ_clean # We want only occurrence data for the US
dups_aghy <- duplicated(aghy_occ_clean[c("lat","lon")]  )
aghy_occ_unique <- aghy_occ_clean[!dups_aghy,] # Remove duplicated data based on latitude and longitude
cat(nrow(aghy_occ_clean)-nrow(aghy_occ_unique), "records are removed")
table(aghy_occ_unique$basisOfRecord) # show the frequency table of “basisOfRecord”
hist(aghy_occ_unique$year,main="",xlim=c(1990,2020),ylim=c(0,500),col="orange2",xlab="years")

aghy_occ_unique %>% 
  dplyr::select(lon, lat, year,country)%>% 
  filter(year %in% (1990:2020) & as.numeric(lon >=-110.374160))->aghy_occ_1990_2020 #  filter from 1990 to 2020 

summary(aghy_occ_1990_2020$year) # show a quick summary of years in the data

## Merge data from the database to GBIF  data
aghy_occ_final_1990_2020<-rbind(aghy_occ_1990_2020[,-3],aghy_occgeoref_1990_2020[,-c(3,4)])
dupfinal1990_2020 <- duplicated(aghy_occ_final_1990_2020[,c(1,2)])
aghy_occ_unique_final_1990_2020 <- aghy_occ_final_1990_2020[!dupfinal1990_2020,]
# unique(aghy_occ_unique_final_1895_2020$country)

## Cross-check coordinates by visual 
# library(maptools)
data(wrld_simpl)
coordinates(aghy_occ_unique_final_1990_2020) <- ~lon+lat
crs(aghy_occ_unique_final_1990_2020) <- crs(wrld_simpl)
ovr_aghy_1990_2020 <- over(aghy_occ_unique_final_1990_2020, wrld_simpl)
head(ovr_aghy_1990_2020)
ctr_aghy_1990_2020 <- ovr_aghy_1990_2020$NAME
i_aghy_1990_2020 <- which(is.na(ctr_aghy_1990_2020)) # points (identified by their record numbers) do not match any country (that is, they are in an ocean)
aghy_occ_unique_final_1990_2020[!i_aghy_1990_2020,]
j_aghy_1990_2020 <- which(ctr_aghy_1990_2020 != aghy_occ_unique_final_1990_2020$country)
cbind(ctr_aghy_1990_2020, aghy_occ_unique_final_1990_2020$country)[j_aghy_1990_2020,]
plot(aghy_occ_unique_final_1990_2020,col="grey")
plot(wrld_simpl, add=T, border='black', lwd=2)

## Spatial thinning of species occurrence records
r_aghy_1990_2020 <- raster(aghy_occ_unique_final_1990_2020) # create a RasterLayer with the extent of aghy_occ_unique_final_1990_2020
res(r_aghy_1990_2020) <- 0.04166667 # set the resolution of the cells
r_aghy_1990_2020 <- extend(r_aghy_1990_2020, extent(r_aghy_1990_2020)+1) # expand (extend) the extent of the RasterLayer a little
aghy_occ_sampled_final_1990_2020 <- gridSample(aghy_occ_unique_final_1990_2020, r_aghy_1990_2020, n=1)
p_aghy_1990_2020 <- rasterToPolygons(r_aghy_1990_2020)
# points(aghy_occ_sampled_final_1990_2020, cex=1, col='red', pch='x')
aghy_occ_sampled_final_1990_2020<-as.data.frame(aghy_occ_sampled_final_1990_2020)
# write_csv(aghy_occ_sampled_final_1990_2020,"/Users/jm200/Library/CloudStorage/Dropbox/Miller Lab/Herbaruim Project/Data/Occurence/aghy_thinned.csv")

## Projection on a map to double check if there is any point out of the US
plot(US_worldclim_recent_norm_stack_final,1)
plot(wrld_simpl, add=TRUE)
points(aghy_occ_sampled_final_1990_2020, col='#0072B2',pch=19)

## Maxent model for *Agrostis hyemalis*----
# aghy_thinned<-read_csv(url("https://www.dropbox.com/scl/fi/xviem0o588sns5c6gtg4s/aghy_thinned.csv?rlkey=4e1c3rd0nt21fn8oooydp49n3&dl=1"))
aghy_thinned<-aghy_occ_sampled_final_1990_2020
coordinates(aghy_thinned) <- ~ lon + lat
crs(aghy_thinned)
myCRS1 <- CRS("+init=epsg:4326") # WGS 84
crs(aghy_thinned) <- myCRS1 # Define the coordinate system that will be used

Aghy_bg <- sampleRandom(x=US_worldclim_recent_norm_stack_final,
                   size=10000,
                   na.rm=T, #removes the 'Not Applicable' points  
                   sp=T) # return spatial points 

aghy_selected <- sample(1:nrow(aghy_occ_sampled_final_1990_2020),nrow(aghy_occ_sampled_final_1990_2020)*0.75)
aghy_train <- aghy_occ_sampled_final_1990_2020[aghy_selected,] # this is the selection to be used for model training
aghy_test <- aghy_occ_sampled_final_1990_2020[-aghy_selected,] # this is the opposite of the selection which will be used for model testing


env_aghy_train <- raster::extract(US_worldclim_recent_norm_stack_final,aghy_train)
env_aghy_test <- raster::extract(US_worldclim_recent_norm_stack_final,aghy_test)
env_aghy_bg <- raster::extract(US_worldclim_recent_norm_stack_final,Aghy_bg)  
aghyPredictors <- rbind(env_aghy_train,env_aghy_bg)
aghyPredictors<-as.data.frame(aghyPredictors)
aghyResponse <- c(rep(1,nrow(env_aghy_train)),
                rep(0,nrow(env_aghy_bg))) 

mod_aghy <- dismo::maxent(aghyPredictors,aghyResponse)
test_aghy_p <- predict(mod_aghy, env_aghy_test)
test_aghy_a <- predict(mod_aghy, env_aghy_bg)
e_aghy <- evaluate(p=test_aghy_p, a=test_aghy_a)
threshold(e_aghy)

par(mfrow=c(1, 3))
# density(e_aghy)
boxplot(e_aghy, col=c('blue', 'red'))
plot(e_aghy, 'ROC')
# plot(e_aghy, 'TPR')
plot(e_aghy, 'kappa')
map_recent_aghy <- predict(US_worldclim_recent_norm_stack_final, mod_aghy,  progress='')

# writeRaster(map_recent_aghy, '/Users/jm200/Library/CloudStorage/Dropbox/Miller Lab/Herbaruim Project/Data/Model/raster_recent_aghy.tif')

par(mfrow=c(1,2))
plot(map_recent_aghy, main='2020')
plot(wrld_simpl, add=TRUE, border='dark grey')

tr_aghy <- threshold(e_aghy, 'spec_sens')
plot(map_recent_aghy > tr_aghy, main='2020')
plot(wrld_simpl, add=TRUE, border='dark grey')

dev.off()
# Agrostis perennans-----

## Download occurrence data from gbif for *Agrostis perennans*
# agpe_occ_raw <- gbif(genus="Agrostis",species="perennans",download=TRUE)
load("/Users/jm200/Library/CloudStorage/Dropbox/Miller Lab/Herbaruim Project/Data/Occurence/agpe_occ_raw.rdata")
head(agpe_occ_raw ) # to view the first few records the occurrence dataset use:

## Clean occurrence data from gbif
agpe_occ <- subset(agpe_occ_raw,(!is.na(lat))&(!is.na(lon))&!is.na(country)) # here we remove erroneous coordinates, where either the latitude or longitude is missing
cat(nrow(agpe_occ_raw)-nrow(agpe_occ), "records are removed") #Show the number of records that are removed from the dataset. 
agpe_occ %>% 
  dplyr::select(country,lon, lat, year,basisOfRecord)%>% 
  filter(country=="United States")->agpe_occ_clean

dups_agpe <- duplicated(agpe_occ_clean[c("lat","lon")]  )
agpe_occ_unique <- agpe_occ_clean[!dups_agpe,] # Remove duplicated data based on latitude and longitude
cat(nrow(agpe_occ_clean)-nrow(agpe_occ_unique), "records are removed")

table(agpe_occ_unique$basisOfRecord) # show the frequency table of “basisOfRecord”
hist(agpe_occ_unique$year,main="",xlim=c(1990,2020),ylim=c(0,500),col="orange2",xlab="years")

agpe_occ_unique %>% 
  dplyr::select(lon, lat, year,country)%>% 
  filter(year %in% (1990:2020) & as.numeric(lon >=-110.374160))->agpe_occ_1990_2020 #  filter from 1990 to 2020 to match climatic data.

summary(agpe_occ_1990_2020$year) # show a quick summary of years in the data

## Merge endophyte host data from  database to GBIF  data
endo_herb_georef %>% 
  dplyr::select(lon, lat, year,Spp_code,Country)%>% 
  dplyr::filter(Spp_code=="AGPE" & Country=="United States" & year %in% (1990:2020)) %>% 
  dplyr::rename(country=Country)->agpe_occgeoref_1990_2020 

agpe_occ_final_1990_2020<-rbind(agpe_occ_1990_2020[,-3],agpe_occgeoref_1990_2020[,-c(3,4)])
dupfinal1990_2020 <- duplicated(agpe_occ_final_1990_2020[,c(1,2)])
agpe_occ_unique_final_1990_2020 <- agpe_occ_final_1990_2020[!dupfinal1990_2020,]

## Cross-check coordinates by visual 
# library(maptools)
coordinates(agpe_occ_unique_final_1990_2020) <- ~lon+lat
crs(agpe_occ_unique_final_1990_2020) <- crs(wrld_simpl)
ovr_agpe_1990_2020 <- over(agpe_occ_unique_final_1990_2020, wrld_simpl)
head(ovr_agpe_1990_2020)
ctr_agpe_1990_2020 <- ovr_agpe_1990_2020$NAME
i_agpe_1990_2020 <- which(is.na(ctr_agpe_1990_2020)) # points (identified by their record numbers) do not match any country (that is, they are in an ocean)
# agpe_occ_unique_final_1990_2020[!i_agpe_1990_2020,]
j_agpe_1990_2020 <- which(ctr_agpe_1990_2020 != agpe_occ_unique_final_1990_2020$country)
cbind(ctr_agpe_1990_2020, agpe_occ_unique_final_1990_2020$country)[j_agpe_1990_2020,]
plot(agpe_occ_unique_final_1990_2020,col="grey")
plot(wrld_simpl, add=T, border='black', lwd=2)

## Spatial thinning of species occurence records
r_agpe_1990_2020 <- raster(agpe_occ_unique_final_1990_2020) # create a RasterLayer with the extent of aghy_occ_unique_final_1990_2020
res(r_agpe_1990_2020) <- 0.04166667 # set the resolution of the cells
r_agpe_1990_2020 <- extend(r_agpe_1990_2020, extent(r_agpe_1990_2020)+1) # expand (extend) the extent of the RasterLayer a little
agpe_occ_sampled_final_1990_2020 <- gridSample(agpe_occ_unique_final_1990_2020, r_agpe_1990_2020, n=1)
p_agpe_1990_2020 <- rasterToPolygons(r_agpe_1990_2020)
# points(agpe_occ_sampled_final_1990_2020, cex=1, col='red', pch='x')
agpe_occ_sampled_final_1990_2020<-as.data.frame(agpe_occ_sampled_final_1990_2020)
# write_csv(agpe_occ_sampled_final_1990_2020,"/Users/jm200/Library/CloudStorage/Dropbox/Miller Lab/Herbaruim Project/Data/Occurence/agpe_thinned.csv")
# agpe_thinned_recent<-read.csv("/Users/jm200/Library/CloudStorage/Dropbox/Miller Lab/Herbaruim Project/Data/Occurence/agpe_thinned_recent.csv", header=T)
# agpe_thinned<-agpe_thinned_recent[,-1]

## Projection on a map to double check if there is any point out of the US
plot(US_worldclim_recent_norm_stack_final,1)
plot(wrld_simpl, add=TRUE)
points(agpe_occ_sampled_final_1990_2020, col='#0072B2',pch=19)

## Maxent model for *Agrostis perennans*----
# aghy_thinned<-read_csv(url("https://www.dropbox.com/scl/fi/xviem0o588sns5c6gtg4s/aghy_thinned.csv?rlkey=4e1c3rd0nt21fn8oooydp49n3&dl=1"))
agpe_thinned<-agpe_occ_sampled_final_1990_2020
coordinates(agpe_thinned) <- ~ lon + lat
crs(agpe_thinned)
myCRS1 <- CRS("+init=epsg:4326") # WGS 84
crs(agpe_thinned) <- myCRS1 # Define the coordinate system that will be used

Agpe_bg <- sampleRandom(x=US_worldclim_recent_norm_stack_final,
                        size=10000,
                        na.rm=T, #removes the 'Not Applicable' points  
                        sp=T) # return spatial points 

agpe_selected <- sample(1:nrow(agpe_occ_sampled_final_1990_2020),nrow(agpe_occ_sampled_final_1990_2020)*0.75)
agpe_train <- agpe_occ_sampled_final_1990_2020[agpe_selected,] # this is the selection to be used for model training
agpe_test <- agpe_occ_sampled_final_1990_2020[-agpe_selected,] # this is the opposite of the selection which will be used for model testing


env_agpe_train <- raster::extract(US_worldclim_recent_norm_stack_final,agpe_train)
env_agpe_test <- raster::extract(US_worldclim_recent_norm_stack_final,agpe_test)
env_agpe_bg <- raster::extract(US_worldclim_recent_norm_stack_final,Agpe_bg)  
agpePredictors <- rbind(env_agpe_train,env_agpe_bg)
agpePredictors<-as.data.frame(agpePredictors)
agpeResponse <- c(rep(1,nrow(env_agpe_train)),
                  rep(0,nrow(env_agpe_bg))) 

mod_agpe <- dismo::maxent(agpePredictors,agpeResponse)
test_agpe_p <- predict(mod_agpe, env_agpe_test)
test_agpe_a <- predict(mod_agpe, env_agpe_bg)
e_agpe <- evaluate(p=test_agpe_p, a=test_agpe_a)
threshold(e_agpe)

par(mfrow=c(1, 3))
# density(e_agpe)
boxplot(e_agpe, col=c('blue', 'red'))
plot(e_agpe, 'ROC')
# plot(e_agpe, 'TPR')
plot(e_agpe, 'kappa')
dev.off()
map_recent_agpe <- predict(US_worldclim_recent_norm_stack_final, mod_agpe,  progress='')

# writeRaster(map_recent_agpe, '/Users/jm200/Library/CloudStorage/Dropbox/Miller Lab/Herbaruim Project/Data/Model/raster_recent_agpe.tif')

par(mfrow=c(1,2))
plot(map_recent_agpe, main='2020')
plot(wrld_simpl, add=TRUE, border='dark grey')

tr_agpe <- threshold(e_agpe, 'spec_sens')
plot(map_recent_agpe > tr_agpe, main='2020')
plot(wrld_simpl, add=TRUE, border='dark grey')
# points(pres_train_aghy, pch='+')
dev.off()

# Elymus virginicus-----

## Download occurrence data from gbif for *Elymus virginicus*
# elvi_occ_raw <- gbif(genus="Elymus",species="virginicus",download=TRUE) 
load("/Users/jm200/Library/CloudStorage/Dropbox/Miller Lab/Herbaruim Project/Data/Occurence/elvi_occ_raw.rdata")
head(elvi_occ_raw)

## Clean occurrence data from gbif

elvi_occ <- subset(elvi_occ_raw,(!is.na(lat))&(!is.na(lon))&!is.na(country)) # here we remove erroneous coordinates, where either the latitude or longitude is missing
cat(nrow(elvi_occ_raw)-nrow(elvi_occ), "records are removed") #Show the number of records that are removed from the dataset. 
elvi_occ %>% 
  dplyr::select(country,lon, lat, year,basisOfRecord)%>% 
  filter(country=="United States")->elvi_occ_clean

dups_elvi <- duplicated(elvi_occ_clean[c("lat","lon")]  )
elvi_occ_unique <- elvi_occ_clean[!dups_elvi,] # Remove duplicated data based on latitude and longitude
cat(nrow(elvi_occ_clean)-nrow(elvi_occ_unique), "records are removed")
table(elvi_occ_unique$basisOfRecord) # show the frequency table of “basisOfRecord”

hist(elvi_occ_unique$year,main="",xlim=c(1990,2020),ylim=c(0,2000),col="orange2",xlab="years")

elvi_occ_unique %>% 
  dplyr::select(lon, lat, year,country)%>% 
  filter(year %in% (1990:2020) & as.numeric(lon >=-110.374160))->elvi_occ_1990_2020 #  filter from 1990 to 2020 to match climatic data.

summary(elvi_occ_1990_2020$year) # show a quick summary of years in the data


## Import endophyte host data from Josh database and merge them to existing online data

endo_herb_georef %>% 
  dplyr::select(lon, lat, year,Spp_code,Country)%>% 
  dplyr::filter(Spp_code=="ELVI" & Country=="United States" & year %in% (1990:2020)) %>% 
  dplyr::rename(country=Country)->elvi_occgeoref_1990_2020 

elvi_occ_final_1990_2020<-rbind(elvi_occ_1990_2020[,-3],elvi_occgeoref_1990_2020[,-c(3,4)])
dupfinal1990_2020 <- duplicated(elvi_occ_final_1990_2020[,c(1,2)])
elvi_occ_unique_final_1990_2020 <- elvi_occ_final_1990_2020[!dupfinal1990_2020,]

## Cross-check coordinates by visual 
# library(maptools)
coordinates(elvi_occ_unique_final_1990_2020) <- ~lon+lat
crs(elvi_occ_unique_final_1990_2020) <- crs(wrld_simpl)
ovr_elvi_1990_2020 <- over(elvi_occ_unique_final_1990_2020, wrld_simpl)
head(ovr_elvi_1990_2020)
ctr_elvi_1990_2020 <- ovr_elvi_1990_2020$NAME
i_elvi_1990_2020 <- which(is.na(ctr_elvi_1990_2020)) # points (identified by their record numbers) do not match any country (that is, they are in an ocean)
elvi_occ_unique_final_1990_2020[!i_elvi_1990_2020,]
j_elvi_1990_2020 <- which(ctr_elvi_1990_2020 != elvi_occ_unique_final_1990_2020$country)
cbind(ctr_elvi_1990_2020, elvi_occ_unique_final_1990_2020$country)[j_elvi_1990_2020,]
plot(elvi_occ_unique_final_1990_2020,col="grey")
plot(wrld_simpl, add=T, border='black', lwd=2)

## Spatial thinning of species occurence records
r_elvi_1990_2020 <- raster(elvi_occ_unique_final_1990_2020) # create a RasterLayer with the extent of aghy_occ_unique_final_1990_2020
res(r_elvi_1990_2020) <- 0.04166667 # set the resolution of the cells
r_elvi_1990_2020 <- extend(r_elvi_1990_2020, extent(r_elvi_1990_2020)+1) # expand (extend) the extent of the RasterLayer a little
elvi_occ_sampled_final_1990_2020 <- gridSample(elvi_occ_unique_final_1990_2020, r_elvi_1990_2020, n=1)
p_elvi_1990_2020 <- rasterToPolygons(r_elvi_1990_2020)
# points(elvi_occ_sampled_final_1990_2020, cex=1, col='red', pch='x')
elvi_occ_sampled_final_1990_2020<-as.data.frame(elvi_occ_sampled_final_1990_2020)
# write_csv(elvi_occ_sampled_final_1990_2020,"/Users/jm200/Library/CloudStorage/Dropbox/Miller Lab/Herbaruim Project/Data/Occurence/elvi_thinned.csv")

# elvi_thinned_recent<-read.csv("/Users/jm200/Library/CloudStorage/Dropbox/Miller Lab/Herbaruim Project/Data/Occurence/elvi_thinned_recent.csv", header=T)
# elvi_thinned<-elvi_thinned_recent[,-1]
## Projection on a map to double check if there is any point out of the US
plot(US_worldclim_recent_norm_stack_final,1)
plot(wrld_simpl, add=TRUE)
points(elvi_occ_sampled_final_1990_2020, col='#0072B2',pch=19)

## Maxent model for *Elymus virginicus*----
# aghy_thinned<-read_csv(url("https://www.dropbox.com/scl/fi/xviem0o588sns5c6gtg4s/aghy_thinned.csv?rlkey=4e1c3rd0nt21fn8oooydp49n3&dl=1"))
elvi_thinned<-elvi_occ_sampled_final_1990_2020
coordinates(elvi_thinned) <- ~ lon + lat
crs(elvi_thinned)
myCRS1 <- CRS("+init=epsg:4326") # WGS 84
crs(elvi_thinned) <- myCRS1 # Define the coordinate system that will be used

elvi_bg <- sampleRandom(x=US_worldclim_recent_norm_stack_final,
                        size=10000,
                        na.rm=T, #removes the 'Not Applicable' points  
                        sp=T) # return spatial points 

elvi_selected <- sample(1:nrow(elvi_occ_sampled_final_1990_2020),nrow(elvi_occ_sampled_final_1990_2020)*0.75)
elvi_train <- elvi_occ_sampled_final_1990_2020[elvi_selected,] # this is the selection to be used for model training
elvi_test <- elvi_occ_sampled_final_1990_2020[-elvi_selected,] # this is the opposite of the selection which will be used for model testing


env_elvi_train <- raster::extract(US_worldclim_recent_norm_stack_final,elvi_train)
env_elvi_test <- raster::extract(US_worldclim_recent_norm_stack_final,elvi_test)
env_elvi_bg <- raster::extract(US_worldclim_recent_norm_stack_final,elvi_bg)  
elviPredictors <- rbind(env_elvi_train,env_elvi_bg)
elviPredictors<-as.data.frame(elviPredictors)
elviResponse <- c(rep(1,nrow(env_elvi_train)),
                  rep(0,nrow(env_elvi_bg))) 

mod_elvi <- dismo::maxent(elviPredictors,elviResponse)
test_elvi_p <- predict(mod_elvi, env_elvi_test)
test_elvi_a <- predict(mod_elvi, env_elvi_bg)
e_elvi <- evaluate(p=test_elvi_p, a=test_elvi_a)
threshold(e_elvi)


par(mfrow=c(1, 3))
# density(e_elvi)
boxplot(e_elvi, col=c('blue', 'red'))
plot(e_elvi, 'ROC')
# plot(e_elvi, 'TPR')
plot(e_elvi, 'kappa')
dev.off()
map_recent_elvi <- predict(US_worldclim_recent_norm_stack_final, mod_elvi,  progress='')

# writeRaster(map_recent_elvi, '/Users/jm200/Library/CloudStorage/Dropbox/Miller Lab/Herbaruim Project/Data/Model/raster_recent_elvi.tif')

par(mfrow=c(1,2))

plot(map_recent_elvi, main='2020')
plot(wrld_simpl, add=TRUE, border='dark grey')

tr_elvi <- threshold(e_elvi, 'spec_sens')

# points(pres_train_agpe, pch='+')
plot(map_recent_elvi > tr_elvi, main='2020')
plot(wrld_simpl, add=TRUE, border='dark grey')
# points(pres_train_aghy, pch='+')
dev.off()

# Final plot

# pdf("/Users/jm200/Library/CloudStorage/Dropbox/Miller Lab/Herbaruim Project/Figures/FigSDMMaxent.pdf",height = 5,width = 10)
# # layout(mat = layout.matrix)
# op<-par(mfrow = c(2,3), mar=c(3,3,2,1),oma = c(5, 5, 1, 3))
# plot(map_old_aghy,xaxt='n',main=substitute(paste(italic("A. hyemalis"))))
# text(x=-80,y=50,label="1925")
# # par(mar=c(2.2,1,3,1))
# plot(map_old_agpe,xaxt='n',yaxt='n', main=substitute(paste(italic("A. perennans"))))
# text(x=-80,y=50,label="1925")
# # par(mar=c(2.2,1,3,1))
# plot(map_old_elvi, xaxt='n',yaxt='n',main=substitute(paste(italic("E. virginicus"))))
# text(x=-80,y=50,label="1925")
# # par(mar=c(5,3,0,1))
# plot(map_recent_aghy, main="")
# text(x=-80,y=50,label="2020")
# # par(mar=c(5,1,0,1))
# plot(map_recent_agpe,yaxt='n', main="")
# text(x=-80,y=50,label="2020")
# # par(mar=c(5,1,0,1))
# plot(map_recent_elvi,yaxt='n', main="")
# text(x=-80,y=50,label="2020")
# title(xlab="Longitude", ylab="Latitude", outer=T, cex.lab=1.2,line = 1.5)
# 
# par(op)
# dev.off()
aghy<-raster::raster("/Users/jm200/Library/CloudStorage/Dropbox/Miller Lab/Herbaruim Project/Data/Model/1990-2020/aghy.tif")
agpe<-raster::raster("/Users/jm200/Library/CloudStorage/Dropbox/Miller Lab/Herbaruim Project/Data/Model/1990-2020/agpe.tif")
elvi<-raster::raster("/Users/jm200/Library/CloudStorage/Dropbox/Miller Lab/Herbaruim Project/Data/Model/1990-2020/elvi.tif")




plot(aghy)
# layout.matrix <- matrix(c(1, 2, 3, 4,5,6), nrow = 2, ncol = 3)
pdf("/Users/jm200/Library/CloudStorage/Dropbox/Miller Lab/Herbaruim Project/Figures/FigSDM_Presence_absence.pdf",height = 8,width = 5)
# layout(mat = layout.matrix)
op<-par(mfrow = c(3,1), mar=c(4,8,1.5,2))
# plot(map_old_aghy > tr_aghy,xaxt='n',main=substitute(paste(italic("A. hyemalis"))))
# text(x=-80,y=50,label="1925")
# # par(mar=c(2.2,1,3,1))
# plot(map_old_agpe > tr_agpe,xaxt='n',yaxt='n', main=substitute(paste(italic("A. perennans"))))
# text(x=-80,y=50,label="1925")
# # par(mar=c(2.2,1,3,1))
# plot(map_old_elvi > tr_elvi, xaxt='n',yaxt='n',main=substitute(paste(italic("E. virginicus"))))
# text(x=-80,y=50,label="1925")
# par(mar=c(5,3,0,1))
plot(aghy >  0.5102097 , main="",legend=FALSE,xlab="", ylab="Latitude",cex.lab=1.5)
mtext("A",side = 3, adj = 0,cex=1.25)
# par(mar=c(5,1,0,1))
plot(agpe > 0.52849, main="",xlab="", ylab="Latitude",legend=FALSE,cex.lab=1.5)
mtext("B",side = 3, adj = 0,cex=1.25)
# par(mar=c(5,1,0,1))
plot(elvi > 0.5422173, main="",legend=FALSE,xlab="Longitude", ylab="Latitude",cex.lab=1.5)
mtext("C",side = 3, adj = 0,cex=1.25)
par(op)
dev.off()


# plot(aghy)
# layout.matrix <- matrix(c(1, 2, 3, 4,5,6), nrow = 2, ncol = 3)
pdf("/Users/jm200/Library/CloudStorage/Dropbox/Miller Lab/Herbaruim Project/Figures/FigSDM.pdf",height = 8,width = 5)
# layout(mat = layout.matrix)
op<-par(mfrow = c(3,1), mar=c(4,8,1.5,2))
# plot(map_old_aghy > tr_aghy,xaxt='n',main=substitute(paste(italic("A. hyemalis"))))
# text(x=-80,y=50,label="1925")
# # par(mar=c(2.2,1,3,1))
# plot(map_old_agpe > tr_agpe,xaxt='n',yaxt='n', main=substitute(paste(italic("A. perennans"))))
# text(x=-80,y=50,label="1925")
# # par(mar=c(2.2,1,3,1))
# plot(map_old_elvi > tr_elvi, xaxt='n',yaxt='n',main=substitute(paste(italic("E. virginicus"))))
# text(x=-80,y=50,label="1925")
# par(mar=c(5,3,0,1))
plot(aghy , main="",xlab="", ylab="Latitude",cex.lab=1.5)
mtext("A",side = 3, adj = 0,cex=1.25)
# par(mar=c(5,1,0,1))
plot(agpe, main="",xlab="", ylab="Latitude",cex.lab=1.5)
mtext("B",side = 3, adj = 0,cex=1.25)
# par(mar=c(5,1,0,1))
plot(elvi, main="",xlab="Longitude", ylab="Latitude",cex.lab=1.5)
mtext("C",side = 3, adj = 0,cex=1.25)
par(op)
dev.off()

aghy_binary <- aghy >  0.5102097
agpe_binary <-agpe > 0.52849
elvi_binary <-elvi > 0.5422173

writeRaster(aghy_binary, '/Users/jm200/Library/CloudStorage/Dropbox/Miller Lab/github/EndoHerbarium/aghy_binary.tif')
writeRaster(agpe_binary, '/Users/jm200/Library/CloudStorage/Dropbox/Miller Lab/github/EndoHerbarium/agpe_binary.tif')
writeRaster(elvi_binary, '/Users/jm200/Library/CloudStorage/Dropbox/Miller Lab/github/EndoHerbarium/elvi_binary.tif')


graph<-list("A","B","C")
pdf("/Users/jm200/Library/CloudStorage/Dropbox/Miller Lab/Herbaruim Project/Figures/FigSDMMaxent_ROC.pdf",height = 5,width = 5)
par(mfrow=c(1, 3))
plot(e_aghy, 'ROC')
plot(e_agpe, 'ROC')
plot(e_elvi, 'ROC')
dev.off()
