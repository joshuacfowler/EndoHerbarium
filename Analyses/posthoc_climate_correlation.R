# Purpose: Perform post-hoc correlations between modelled spatially varying trends in endophyte prevalence and observed change in climate
# Authors: Joshua Fowler
# Updated: June 25, 2024

library(tidyverse) # for data manipulation and ggplot
library(INLA) # for fitting integrated nested Laplace approximation models
library(inlabru)
library(fmesher)

library(prism) # to import prism raster files

library(sf)
# library(rmapshaper)
library(terra)
library(tidyterra)


library(patchwork)
library(ggmap)

species_colors <- c("#1b9e77","#d95f02","#7570b3")
species_names <- c("A. hyemalis", "A. perennans", "E. virginicus")


################################################################################
############ Download and read in prism data rasters ###########################
################################################################################

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
tmean_annual_recent_cv <- tmean_annual_recent_sd/(tmean_annual_recent_norm+abs(terra::minmax(tmean_annual_recent_norm)[1,]))
tmean_spring_recent_cv <- tmean_spring_recent_sd/(tmean_spring_recent_norm+abs(terra::minmax(tmean_spring_recent_norm)[1,]))
tmean_summer_recent_cv <- tmean_summer_recent_sd/(tmean_summer_recent_norm+abs(terra::minmax(tmean_summer_recent_norm)[1,]))
tmean_autumn_recent_cv <- tmean_autumn_recent_sd/(tmean_autumn_recent_norm+abs(terra::minmax(tmean_autumn_recent_norm)[1,]))

tmean_annual_old_cv <- tmean_annual_old_sd/(tmean_annual_old_norm+abs(terra::minmax(tmean_annual_old_norm)[1,]))
tmean_spring_old_cv <- tmean_spring_old_sd/(tmean_spring_old_norm+abs(terra::minmax(tmean_spring_old_norm)[1,]))
tmean_summer_old_cv <- tmean_summer_old_sd/(tmean_summer_old_norm+abs(terra::minmax(tmean_summer_old_norm)[1,]))
tmean_autumn_old_cv <- tmean_autumn_old_sd/(tmean_autumn_old_norm+abs(terra::minmax(tmean_autumn_old_norm)[1,]))



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

epsg6703km <- paste(
  "+proj=aea +lat_0=23 +lon_0=-96 +lat_1=29.5",
  "+lat_2=45.5 +x_0=0 +y_0=0 +datum=NAD83",
  "+units=km +no_defs"
)
dd_crs <- "GEOGCRS[\"NAD83\",\n    DATUM[\"North American Datum 1983\",\n        ELLIPSOID[\"GRS 1980\",6378137,298.257222101004,\n            LENGTHUNIT[\"metre\",1]]],\n    PRIMEM[\"Greenwich\",0,\n        ANGLEUNIT[\"degree\",0.0174532925199433]],\n    CS[ellipsoidal,2],\n        AXIS[\"geodetic latitude (Lat)\",north,\n            ORDER[1],\n            ANGLEUNIT[\"degree\",0.0174532925199433]],\n        AXIS[\"geodetic longitude (Lon)\",east,\n            ORDER[2],\n            ANGLEUNIT[\"degree\",0.0174532925199433]],\n    ID[\"EPSG\",4269]]"

raster::crs(tmean_annual_difference) <- dd_crs
raster::crs(tmean_spring_difference) <- dd_crs
raster::crs(tmean_summer_difference) <- dd_crs
raster::crs(tmean_autumn_difference) <- dd_crs

raster::crs(tmean_annual_sd_difference) <- dd_crs
raster::crs(tmean_spring_sd_difference) <- dd_crs
raster::crs(tmean_summer_sd_difference) <- dd_crs
raster::crs(tmean_autumn_sd_difference) <- dd_crs

raster::crs(tmean_annual_cv_difference) <- dd_crs 
raster::crs(tmean_spring_cv_difference) <- dd_crs 
raster::crs(tmean_summer_cv_difference) <- dd_crs 
raster::crs(tmean_autumn_cv_difference) <- dd_crs 



raster::crs(ppt_annual_difference) <- dd_crs
raster::crs(ppt_spring_difference) <- dd_crs
raster::crs(ppt_summer_difference) <- dd_crs
raster::crs(ppt_autumn_difference) <- dd_crs

raster::crs(ppt_annual_sd_difference) <- dd_crs
raster::crs(ppt_spring_sd_difference) <- dd_crs
raster::crs(ppt_summer_sd_difference) <- dd_crs
raster::crs(ppt_autumn_sd_difference) <- dd_crs

raster::crs(ppt_annual_cv_difference) <- dd_crs 
raster::crs(ppt_spring_cv_difference) <- dd_crs 
raster::crs(ppt_summer_cv_difference) <- dd_crs 
raster::crs(ppt_autumn_cv_difference) <- dd_crs 





######### Next extracting these climate values at our pixels where we predict change in endophyte prevalence ####
# These are the objects we saved from the model output that include out spatially varying slopes

# getting the pixel values of slopes as a dataframe
svc.pred_aghy <- readRDS( file = "svc.pred_aghy.Rds")
svc.pred_agpe <- readRDS( file = "svc.pred_agpe.Rds")
svc.pred_elvi <- readRDS( file = "svc.pred_elvi.Rds")


vrt_aghy_dd <- svc.pred_aghy  %>% 
  spTransform(dd_crs)

vrt_agpe_dd <- svc.pred_agpe  %>% 
  spTransform(dd_crs)

vrt_elvi_dd <- svc.pred_elvi  %>% 
  spTransform(dd_crs)


aghy_prism_diff_pred_df <- tibble(lon = vrt_aghy_dd@coords[,1], lat = vrt_aghy_dd@coords[,2],
                                  coords.x1 = svc.pred_aghy@coords[,1], coords.x2 = svc.pred_aghy@coords[,2],
                                  tmean_annual_diff = as.numeric(unlist(terra::extract(tmean_annual_difference, vrt_aghy_dd@coords))),
                                  tmean_spring_diff = as.numeric(unlist(terra::extract(tmean_spring_difference, vrt_aghy_dd@coords))),
                                  tmean_summer_diff = as.numeric(unlist(terra::extract(tmean_summer_difference, vrt_aghy_dd@coords))),
                                  tmean_autumn_diff = as.numeric(unlist(terra::extract(tmean_autumn_difference, vrt_aghy_dd@coords))),
                                  ppt_annual_diff = as.numeric(unlist(terra::extract(ppt_annual_difference, vrt_aghy_dd@coords))),
                                  ppt_spring_diff = as.numeric(unlist(terra::extract(ppt_spring_difference, vrt_aghy_dd@coords))),
                                  ppt_summer_diff = as.numeric(unlist(terra::extract(ppt_summer_difference, vrt_aghy_dd@coords))),
                                  ppt_autumn_diff = as.numeric(unlist(terra::extract(ppt_autumn_difference, vrt_aghy_dd@coords))),
                                  
                                  tmean_annual_sd_diff = as.numeric(unlist(terra::extract(tmean_annual_sd_difference, vrt_aghy_dd@coords))),
                                  tmean_spring_sd_diff = as.numeric(unlist(terra::extract(tmean_spring_sd_difference, vrt_aghy_dd@coords))),
                                  tmean_summer_sd_diff = as.numeric(unlist(terra::extract(tmean_summer_sd_difference, vrt_aghy_dd@coords))),
                                  tmean_autumn_sd_diff = as.numeric(unlist(terra::extract(tmean_autumn_sd_difference, vrt_aghy_dd@coords))),
                                  ppt_annual_sd_diff = as.numeric(unlist(terra::extract(ppt_annual_sd_difference, vrt_aghy_dd@coords))),
                                  ppt_spring_sd_diff = as.numeric(unlist(terra::extract(ppt_spring_sd_difference, vrt_aghy_dd@coords))),
                                  ppt_summer_sd_diff = as.numeric(unlist(terra::extract(ppt_summer_sd_difference, vrt_aghy_dd@coords))),
                                  ppt_autumn_sd_diff = as.numeric(unlist(terra::extract(ppt_autumn_sd_difference, vrt_aghy_dd@coords))),
                                  
                                  tmean_annual_cv_diff = as.numeric(unlist(terra::extract(tmean_annual_cv_difference, vrt_aghy_dd@coords))),
                                  tmean_spring_cv_diff = as.numeric(unlist(terra::extract(tmean_spring_cv_difference, vrt_aghy_dd@coords))),
                                  tmean_summer_cv_diff = as.numeric(unlist(terra::extract(tmean_summer_cv_difference, vrt_aghy_dd@coords))),
                                  tmean_autumn_cv_diff = as.numeric(unlist(terra::extract(tmean_autumn_cv_difference, vrt_aghy_dd@coords))),
                                  ppt_annual_cv_diff = as.numeric(unlist(terra::extract(ppt_annual_cv_difference, vrt_aghy_dd@coords))),
                                  ppt_spring_cv_diff = as.numeric(unlist(terra::extract(ppt_spring_cv_difference, vrt_aghy_dd@coords))),
                                  ppt_summer_cv_diff = as.numeric(unlist(terra::extract(ppt_summer_cv_difference, vrt_aghy_dd@coords))),
                                  ppt_autumn_cv_diff = as.numeric(unlist(terra::extract(ppt_autumn_cv_difference, vrt_aghy_dd@coords))),
                                  svc.pred_aghy@data) 


write_csv(aghy_prism_diff_pred_df, file = "aghy_prism_diff_pred_df.csv")



agpe_prism_diff_pred_df <- tibble(lon = vrt_agpe_dd@coords[,1], lat = vrt_agpe_dd@coords[,2],
                                  coords.x1 = svc.pred_agpe@coords[,1], coords.x2 = svc.pred_agpe@coords[,2],
                                  tmean_annual_diff = as.numeric(unlist(terra::extract(tmean_annual_difference, vrt_agpe_dd@coords))),
                                  tmean_spring_diff = as.numeric(unlist(terra::extract(tmean_spring_difference, vrt_agpe_dd@coords))),
                                  tmean_summer_diff = as.numeric(unlist(terra::extract(tmean_summer_difference, vrt_agpe_dd@coords))),
                                  tmean_autumn_diff = as.numeric(unlist(terra::extract(tmean_autumn_difference, vrt_agpe_dd@coords))),
                                  ppt_annual_diff = as.numeric(unlist(terra::extract(ppt_annual_difference, vrt_agpe_dd@coords))),
                                  ppt_spring_diff = as.numeric(unlist(terra::extract(ppt_spring_difference, vrt_agpe_dd@coords))),
                                  ppt_summer_diff = as.numeric(unlist(terra::extract(ppt_summer_difference, vrt_agpe_dd@coords))),
                                  ppt_autumn_diff = as.numeric(unlist(terra::extract(ppt_autumn_difference, vrt_agpe_dd@coords))),
                                  
                                  tmean_annual_sd_diff = as.numeric(unlist(terra::extract(tmean_annual_sd_difference, vrt_agpe_dd@coords))),
                                  tmean_spring_sd_diff = as.numeric(unlist(terra::extract(tmean_spring_sd_difference, vrt_agpe_dd@coords))),
                                  tmean_summer_sd_diff = as.numeric(unlist(terra::extract(tmean_summer_sd_difference, vrt_agpe_dd@coords))),
                                  tmean_autumn_sd_diff = as.numeric(unlist(terra::extract(tmean_autumn_sd_difference, vrt_agpe_dd@coords))),
                                  ppt_annual_sd_diff = as.numeric(unlist(terra::extract(ppt_annual_sd_difference, vrt_agpe_dd@coords))),
                                  ppt_spring_sd_diff = as.numeric(unlist(terra::extract(ppt_spring_sd_difference, vrt_agpe_dd@coords))),
                                  ppt_summer_sd_diff = as.numeric(unlist(terra::extract(ppt_summer_sd_difference, vrt_agpe_dd@coords))),
                                  ppt_autumn_sd_diff = as.numeric(unlist(terra::extract(ppt_autumn_sd_difference, vrt_agpe_dd@coords))),
                                  
                                  tmean_annual_cv_diff = as.numeric(unlist(terra::extract(tmean_annual_cv_difference, vrt_agpe_dd@coords))),
                                  tmean_spring_cv_diff = as.numeric(unlist(terra::extract(tmean_spring_cv_difference, vrt_agpe_dd@coords))),
                                  tmean_summer_cv_diff = as.numeric(unlist(terra::extract(tmean_summer_cv_difference, vrt_agpe_dd@coords))),
                                  tmean_autumn_cv_diff = as.numeric(unlist(terra::extract(tmean_autumn_cv_difference, vrt_agpe_dd@coords))),
                                  ppt_annual_cv_diff = as.numeric(unlist(terra::extract(ppt_annual_cv_difference, vrt_agpe_dd@coords))),
                                  ppt_spring_cv_diff = as.numeric(unlist(terra::extract(ppt_spring_cv_difference, vrt_agpe_dd@coords))),
                                  ppt_summer_cv_diff = as.numeric(unlist(terra::extract(ppt_summer_cv_difference, vrt_agpe_dd@coords))),
                                  ppt_autumn_cv_diff = as.numeric(unlist(terra::extract(ppt_autumn_cv_difference, vrt_agpe_dd@coords))),
                                  svc.pred_agpe@data) 

write_csv(agpe_prism_diff_pred_df, file = "agpe_prism_diff_pred_df.csv")

elvi_prism_diff_pred_df <- tibble(lon = vrt_elvi_dd@coords[,1], lat = vrt_elvi_dd@coords[,2],
                                  coords.x1 = svc.pred_elvi@coords[,1], coords.x2 = svc.pred_elvi@coords[,2],
                                  tmean_annual_diff = as.numeric(unlist(terra::extract(tmean_annual_difference, vrt_elvi_dd@coords))),
                                  tmean_spring_diff = as.numeric(unlist(terra::extract(tmean_spring_difference, vrt_elvi_dd@coords))),
                                  tmean_summer_diff = as.numeric(unlist(terra::extract(tmean_summer_difference, vrt_elvi_dd@coords))),
                                  tmean_autumn_diff = as.numeric(unlist(terra::extract(tmean_autumn_difference, vrt_elvi_dd@coords))),
                                  ppt_annual_diff = as.numeric(unlist(terra::extract(ppt_annual_difference, vrt_elvi_dd@coords))),
                                  ppt_spring_diff = as.numeric(unlist(terra::extract(ppt_spring_difference, vrt_elvi_dd@coords))),
                                  ppt_summer_diff = as.numeric(unlist(terra::extract(ppt_summer_difference, vrt_elvi_dd@coords))),
                                  ppt_autumn_diff = as.numeric(unlist(terra::extract(ppt_autumn_difference, vrt_elvi_dd@coords))),
                                  
                                  tmean_annual_sd_diff = as.numeric(unlist(terra::extract(tmean_annual_sd_difference, vrt_elvi_dd@coords))),
                                  tmean_spring_sd_diff = as.numeric(unlist(terra::extract(tmean_spring_sd_difference, vrt_elvi_dd@coords))),
                                  tmean_summer_sd_diff = as.numeric(unlist(terra::extract(tmean_summer_sd_difference, vrt_elvi_dd@coords))),
                                  tmean_autumn_sd_diff = as.numeric(unlist(terra::extract(tmean_autumn_sd_difference, vrt_elvi_dd@coords))),
                                  ppt_annual_sd_diff = as.numeric(unlist(terra::extract(ppt_annual_sd_difference, vrt_elvi_dd@coords))),
                                  ppt_spring_sd_diff = as.numeric(unlist(terra::extract(ppt_spring_sd_difference, vrt_elvi_dd@coords))),
                                  ppt_summer_sd_diff = as.numeric(unlist(terra::extract(ppt_summer_sd_difference, vrt_elvi_dd@coords))),
                                  ppt_autumn_sd_diff = as.numeric(unlist(terra::extract(ppt_autumn_sd_difference, vrt_elvi_dd@coords))),
                                  
                                  tmean_annual_cv_diff = as.numeric(unlist(terra::extract(tmean_annual_cv_difference, vrt_elvi_dd@coords))),
                                  tmean_spring_cv_diff = as.numeric(unlist(terra::extract(tmean_spring_cv_difference, vrt_elvi_dd@coords))),
                                  tmean_summer_cv_diff = as.numeric(unlist(terra::extract(tmean_summer_cv_difference, vrt_elvi_dd@coords))),
                                  tmean_autumn_cv_diff = as.numeric(unlist(terra::extract(tmean_autumn_cv_difference, vrt_elvi_dd@coords))),
                                  ppt_annual_cv_diff = as.numeric(unlist(terra::extract(ppt_annual_cv_difference, vrt_elvi_dd@coords))),
                                  ppt_spring_cv_diff = as.numeric(unlist(terra::extract(ppt_spring_cv_difference, vrt_elvi_dd@coords))),
                                  ppt_summer_cv_diff = as.numeric(unlist(terra::extract(ppt_summer_cv_difference, vrt_elvi_dd@coords))),
                                  ppt_autumn_cv_diff = as.numeric(unlist(terra::extract(ppt_autumn_cv_difference, vrt_elvi_dd@coords))),
                                  svc.pred_elvi@data) 


write_csv(elvi_prism_diff_pred_df, file = "elvi_prism_diff_pred_df.csv")


######### reading in and looking at the prism data from above ####

aghy_prism_diff_pred_df <- read_csv("aghy_prism_diff_pred_df.csv")
agpe_prism_diff_pred_df <- read_csv("agpe_prism_diff_pred_df.csv")
elvi_prism_diff_pred_df <- read_csv("elvi_prism_diff_pred_df.csv")

# Plotting the change in climate at our pixel values

prism_diff_pred_df <- aghy_prism_diff_pred_df %>% 
  full_join(agpe_prism_diff_pred_df) %>% 
  full_join(elvi_prism_diff_pred_df) %>% 
  pivot_longer(contains("_diff")) %>% 
  separate(name, into = c("climate", "season"), sep = "_", extra = "drop", remove = FALSE) %>% 
  # filter(!grepl("_cv_", name) & ! grepl("annual", name)) %>% 
  mutate(season_f = factor(season, level = c("spring", "summer", "autumn")),
         moment = case_when(grepl("_sd_", name) ~ "sd",
                            TRUE ~ "mean"))
  
prism_summary <- prism_diff_pred_df %>% 
  filter(species == "A. hyemalis") %>% 
  ungroup() %>% 
  group_by(climate, season, moment) %>% 
  dplyr::summarise(mean_change = mean(value),
                   max_change = max(value),
                   min_change = min(value))


AGHY_tmean_change_plot <- ggplot(filter(prism_diff_pred_df, climate == "tmean" & species == "A. hyemalis"))+
  geom_tile(aes(x = coords.x1, y = coords.x2, fill = value)) +
  coord_sf()+
  facet_wrap(~moment + season_f)+
  scale_fill_viridis_c(option = "magma") + labs(x = "Lon.", y = "Lat.", fill = "Change in ºC")+
  theme(strip.text = element_text(size = rel(1)))
AGHY_tmean_change_plot

AGHY_ppt_change_plot <- ggplot(filter(prism_diff_pred_df, climate == "ppt" & species == "A. hyemalis"))+
  geom_tile(aes(x = coords.x1, y = coords.x2, fill = value)) +
  coord_sf()+
  facet_wrap(~ moment+season)+
  scale_fill_viridis_c(option = "magma") + labs(x = "Lon.", y = "Lat.", fill = "Change in mm.")+
  theme(strip.text = element_text(size = rel(1)))
AGHY_ppt_change_plot


AGHY_climate_change_plot <- AGHY_tmean_change_plot/AGHY_ppt_change_plot + plot_annotation(tag_levels = "A")

ggsave(AGHY_climate_change_plot, filename = "Plots/AGHY_climate_change_plot.png", width = 8, height = 10)



AGPE_tmean_change_plot <- ggplot(filter(prism_diff_pred_df, climate == "tmean" & species == "A. perennans"))+
  geom_tile(aes(x = coords.x1, y = coords.x2, fill = value)) +
  coord_sf()+
  facet_wrap(~moment + season_f)+
  scale_fill_viridis_c(option = "magma") + labs(x = "Lon.", y = "Lat.", fill = "Change in ºC")+
  theme(strip.text = element_text(size = rel(1)))
# AGPE_tmean_change_plot

AGPE_ppt_change_plot <- ggplot(filter(prism_diff_pred_df, climate == "ppt" & species == "A. perennans"))+
  geom_tile(aes(x = coords.x1, y = coords.x2, fill = value)) +
  coord_sf()+
  facet_wrap(~ moment+season)+
  scale_fill_viridis_c(option = "magma") + labs(x = "Lon.", y = "Lat.", fill = "Change in mm.")+
  theme(strip.text = element_text(size = rel(1)))
# AGPE_ppt_change_plot


AGPE_climate_change_plot <- AGPE_tmean_change_plot/AGPE_ppt_change_plot + plot_annotation( tag_levels = "A")

ggsave(AGPE_climate_change_plot, filename = "Plots/AGPE_climate_change_plot.png", width = 8, height = 10)



ELVI_tmean_change_plot <- ggplot(filter(prism_diff_pred_df, climate == "tmean" & species == "E. virginicus"))+
  geom_tile(aes(x = coords.x1, y = coords.x2, fill = value)) +
  coord_sf()+
  facet_wrap(~moment + season_f)+
  scale_fill_viridis_c(option = "magma") + labs(x = "Lon.", y = "Lat.", fill = "Change in ºC")+
  theme(strip.text = element_text(size = rel(1)))
# ELVI_tmean_change_plot

ELVI_ppt_change_plot <- ggplot(filter(prism_diff_pred_df, climate == "ppt" & species == "E. virginicus"))+
  geom_tile(aes(x = coords.x1, y = coords.x2, fill = value)) +
  coord_sf()+
  facet_wrap(~ moment+season)+
  scale_fill_viridis_c(option = "magma") + labs(x = "Lon.", y = "Lat.", fill = "Change in mm.")+
  theme(strip.text = element_text(size = rel(1)))
# ELVI_ppt_change_plot


ELVI_climate_change_plot <- ELVI_tmean_change_plot/ELVI_ppt_change_plot + plot_annotation( tag_levels = "A")

ggsave(ELVI_climate_change_plot, filename = "Plots/ELVI_climate_change_plot.png", width = 8, height = 10)






ggplot()+
  geom_tile(data = elvi_prism_diff_pred_df, aes(x = coords.x1, y = coords.x2, fill = tmean_autumn_diff))+
  lims(x = c(-1000,2200), y = c(500,3000))

ggplot()+
  geom_point(data = elvi_prism_diff_pred_df, aes(x = tmean_annual_sd_diff, y = tmean_summer_sd_diff))

svc.pred_climate <- aghy_prism_diff_pred_df %>% 
  full_join(agpe_prism_diff_pred_df) %>% 
  full_join(elvi_prism_diff_pred_df)

# Making a version of the dataset that subsamples locations
svc.pred_climate_subsample <- svc.pred_climate %>% 
group_by(species) %>%
sample_n(size =250) %>% 
  ungroup() %>% 
  mutate(Intercept = 1)

pixel_map <- ggplot(filter(svc.pred_climate_subsample, species == "A. hyemalis"))+
  coord_sf()+
  geom_tile(aes(x = coords.x1, y = coords.x2, fill = ppt_spring_sd_diff))+
  scale_fill_viridis_c(option = "magma") + labs(x = "Lon.", y = "Lat.", fill = "Change in mm.")+
  facet_wrap(~species)
ggsave(pixel_map, filename = "Plots/pixel_maps_AGHYprecip.png", width = 5, height = 5)

pixel_map <- ggplot(filter(svc.pred_climate_subsample, species == "A. hyemalis"))+
  coord_sf()+
  geom_tile(aes(x = coords.x1, y = coords.x2, fill = mean))+
  scale_fill_viridis_c(option = "turbo", na.value = "transparent")+
  labs(x = "Lon.", y = "Lat.", fill = "% change/year")+
  facet_wrap(~species)
ggsave(pixel_map, filename = "Plots/pixel_maps_AGHYtrends.png", width = 5, height = 5)


svc.pred_climate_subsample_long <- svc.pred_climate_subsample %>%
  pivot_longer(cols = contains("_diff")) %>%
  filter(!grepl("cv", name)) %>% 
  filter(!grepl("annual", name)) %>% 
  mutate(name = fct_relevel(name, c('tmean_spring_diff', 'tmean_summer_diff', "tmean_autumn_diff",
                                    'tmean_spring_sd_diff', 'tmean_summer_sd_diff', "tmean_autumn_sd_diff",
                                    'ppt_spring_diff', 'ppt_summer_diff', "ppt_autumn_diff",
                                    'ppt_spring_sd_diff', 'ppt_summer_sd_diff', "ppt_autumn_sd_diff"))) 
# 




# Making a version of the dataset that takes pixels only at locations with specimens
# list of unique locations
# endo_herb_AGHY$lat
# 
# svc.pred_climate_specimens <- svc.pred_climate %>% 
#   group_by(species) %>% 
#   
#   sample_n(size =200) %>% 
#   ungroup() %>% 
#   mutate(Intercept = 1)





# Looking at the data
ppt_regression_plot <- ggplot(filter(svc.pred_climate_subsample_long, grepl("ppt", name)&!grepl("sd", name)))+
  geom_point(aes(x = value, y = q0.5, color = species), alpha = .2)+
  geom_smooth(aes(x = value, y = q0.5, group = species, fill = species), color = "black", method = "lm")+
  scale_color_manual(values = species_colors) +
  facet_wrap(~species+ factor(name), scales = "free", nrow = 3) +
  theme(text = element_text(size = 8))+
  theme_light()+
  labs(y = "Trend in endophyte prevalence", x = expression("Change in Precip. (mm"^-1*")"))
ppt_regression_plot


ppt_sd_regression_plot <- ggplot(filter(svc.pred_climate_subsample_long, grepl("ppt", name)&grepl("sd", name) & !grepl("annual", name)))+
  geom_point(aes(x = value, y = mean,color = species), alpha = .2)+
  geom_smooth(aes(x = value, y = mean, group = species, fill = species), color = "black", method = "lm")+
  scale_color_manual(values = species_colors) +
  facet_wrap(~species+ factor(name), scales = "free", nrow = 3) +
  theme(strip.text = element_text(size = rel(400)))+
  theme_light()+
  labs(y = "Trend in endophyte prevalence", x = expression("Change in Precip. (mm"^-1*")"))
ppt_sd_regression_plot


tmean_regression_plot <- ggplot(filter(svc.pred_climate_subsample_long, grepl("tmean", name)&!grepl("sd", name)))+
  geom_point(aes(x = value, y = q0.5, color = species), alpha = .2)+
  geom_smooth(aes(x = value, y = q0.5, group = species, fill = species), color = "black", method = "lm")+
  scale_color_manual(values = species_colors) +
  facet_wrap(~species+ factor(name), scales = "free", nrow = 3) +
  theme(text = element_text(size = 8))+
  theme_light()+
  labs(y = "Trend in endophyte prevalence", x = expression("Change in Temp. (C)"))
tmean_regression_plot


tmean_sd_regression_plot <- ggplot(filter(svc.pred_climate_subsample_long, grepl("tmean", name)&grepl("sd", name) & !grepl("annual", name)))+
  geom_point(aes(x = value, y = mean,color = species), alpha = .2)+
  geom_smooth(aes(x = value, y = mean, group = species, fill = species), color = "black", method = "lm")+
  scale_color_manual(values = species_colors) +
  facet_wrap(~species+ factor(name), nrow = 3) +
  theme(strip.text = element_text(size = rel(400)))+
  theme_light()+
  labs(y = "Trend in endophyte prevalence", x = expression("Change in Temp. (C)"))
tmean_sd_regression_plot


######### perform a regression with INLA to test the relationship between change in prevalence and change in seasonal climate ####

data<- svc.pred_climate_subsample %>% 
  st_as_sf(coords = c("lon", "lat"), crs = 4326, remove = FALSE) %>% 
  st_transform(epsg6703km) %>% 
  mutate(
    easting = st_coordinates(.)[, 1],
    northing = st_coordinates(.)[, 2]
  ) 
fit_list <- list()

for(s in 1:3){
  data_temp <- data[data$species == species_names[s],]

  cmp <- ~  Intercept(1) + 
    spring_ppt(main = ppt_spring_diff, model = "linear", mean.linear=0, prec.linear=0.001) + summer_ppt(main = ppt_summer_diff, model = "linear", mean.linear=0, prec.linear=0.001) + autumn_ppt(main = ppt_autumn_diff, model = "linear", mean.linear=0, prec.linear=0.001) +
    spring_ppt_sd(main = ppt_spring_sd_diff, model = "linear", mean.linear=0, prec.linear=0.001) + summer_ppt_sd(main = ppt_summer_sd_diff, model = "linear", mean.linear=0, prec.linear=0.001) + autumn_ppt_sd(main = ppt_autumn_sd_diff, model = "linear", mean.linear=0, prec.linear=0.001) +
    spring_tmean(main = tmean_spring_diff, model = "linear", mean.linear=0, prec.linear=0.001) + summer_tmean(main = tmean_summer_diff, model = "linear", mean.linear=0, prec.linear=0.001) + autumn_tmean(main = tmean_autumn_diff, model = "linear", mean.linear=0, prec.linear=0.001) +
    spring_tmean_sd(main = tmean_spring_sd_diff, model = "linear", mean.linear=0, prec.linear=0.001) + summer_tmean_sd(main = tmean_summer_sd_diff, model = "linear", mean.linear=0, prec.linear=0.001) + autumn_tmean_sd(main = tmean_autumn_sd_diff, model = "linear", mean.linear=0, prec.linear=0.001)
    # spring_PxT(main = ppt_spring_diff*tmean_spring_diff, model = "linear", mean.linear=0, prec.linear=0.001) + 
    # summer_PxT(main = ppt_summer_diff*tmean_summer_diff, model = "linear", mean.linear=0, prec.linear=0.001) +
    # autumn_PxT(main = ppt_autumn_diff*tmean_autumn_diff, model = "linear", mean.linear=0, prec.linear=0.001)+
    # spring_Pxsd(main = ppt_spring_diff*ppt_spring_sd_diff, model = "linear", mean.linear=0, prec.linear=0.001) + 
    # summer_Pxsd(main = ppt_summer_diff*ppt_summer_sd_diff, model = "linear", mean.linear=0, prec.linear=0.001) +
    # autumn_Pxsd(main = ppt_autumn_diff*ppt_autumn_sd_diff, model = "linear", mean.linear=0, prec.linear=0.001)+
    # spring_Txsd(main = tmean_spring_diff*tmean_spring_sd_diff, model = "linear", mean.linear=0, prec.linear=0.001) + 
    # summer_Txsd(main = tmean_summer_diff*tmean_summer_sd_diff, model = "linear", mean.linear=0, prec.linear=0.001) +
    # autumn_Txsd(main = tmean_autumn_diff*tmean_autumn_sd_diff, model = "linear", mean.linear=0, prec.linear=0.001)
  
  
  

# setting up the model formula
fml.climate <- mean ~ 0 + Intercept + spring_ppt + summer_ppt + autumn_ppt + spring_ppt_sd + summer_ppt_sd + autumn_ppt_sd +
   spring_tmean + summer_tmean + autumn_tmean + spring_tmean_sd + summer_tmean_sd + autumn_tmean_sd  


lik_climate <- like(formula = fml.climate,
                    family = "gaussian",
                    data = data_temp)

fit <- bru(cmp,
           lik_climate,
           options = list(
             control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE),
             control.inla = list(int.strategy = "eb"),
             verbose = TRUE,
             bru_max_iter = 10))

fit_list[[species_names[s]]] <- fit 

}




# saveRDS(fit_list, file = "climate_fits.rds")
fit_list[[1]]$mode$mode.status
fit_list[[1]]$dic$dic
fit_list[[1]]$summary.fixed
fit_list[[1]]$summary.random



######### generate predictions for plotting ####


prediction_list <- list()
posterior_list <- list()
for(s in 1:3){
data_temp <- data[data$species == species_names[s],]

min_spring_ppt <- min(data_temp$ppt_spring_diff)
mean_spring_ppt <- mean(data_temp$ppt_spring_diff)
max_spring_ppt <- max(data_temp$ppt_spring_diff)
min_summer_ppt <- min(data_temp$ppt_summer_diff)
mean_summer_ppt <- mean(data_temp$ppt_summer_diff)
max_summer_ppt <- max(data_temp$ppt_summer_diff)
min_autumn_ppt <- min(data_temp$ppt_autumn_diff)
mean_autumn_ppt <- mean(data_temp$ppt_autumn_diff)
max_autumn_ppt <- max(data_temp$ppt_autumn_diff)


min_spring_ppt_sd <- min(data_temp$ppt_spring_sd_diff)
mean_spring_ppt_sd <- mean(data_temp$ppt_spring_sd_diff)
max_spring_ppt_sd <- max(data_temp$ppt_spring_sd_diff)
min_summer_ppt_sd <- min(data_temp$ppt_summer_sd_diff)
mean_summer_ppt_sd <- mean(data_temp$ppt_summer_sd_diff)
max_summer_ppt_sd <- max(data_temp$ppt_summer_sd_diff)
min_autumn_ppt_sd <- min(data_temp$ppt_autumn_sd_diff)
mean_autumn_ppt_sd <- mean(data_temp$ppt_autumn_sd_diff)
max_autumn_ppt_sd <- max(data_temp$ppt_autumn_sd_diff)


min_spring_tmean <- min(data_temp$tmean_spring_diff)
mean_spring_tmean <- mean(data_temp$tmean_spring_diff)
max_spring_tmean <- max(data_temp$tmean_spring_diff)
min_summer_tmean <- min(data_temp$tmean_summer_diff)
mean_summer_tmean <- mean(data_temp$tmean_summer_diff)
max_summer_tmean <- max(data_temp$tmean_summer_diff)
min_autumn_tmean <- min(data_temp$tmean_autumn_diff)
mean_autumn_tmean <- mean(data_temp$tmean_autumn_diff)
max_autumn_tmean <- max(data_temp$tmean_autumn_diff)


min_spring_tmean_sd <- min(data_temp$tmean_spring_sd_diff)
mean_spring_tmean_sd <- mean(data_temp$tmean_spring_sd_diff)
max_spring_tmean_sd <- max(data_temp$tmean_spring_sd_diff)
min_summer_tmean_sd <- min(data_temp$tmean_summer_sd_diff)
mean_summer_tmean_sd <- mean(data_temp$tmean_summer_sd_diff)
max_summer_tmean_sd <- max(data_temp$tmean_summer_sd_diff)
min_autumn_tmean_sd <- min(data_temp$tmean_autumn_sd_diff)
mean_autumn_tmean_sd <- mean(data_temp$tmean_autumn_sd_diff)
max_autumn_tmean_sd <- max(data_temp$tmean_autumn_sd_diff)

predict_vector <- c()

preddata <- tibble(species = species_names[s],
                   ppt_spring_diff = seq(min_spring_ppt, max_spring_ppt, length.out = 50),
                   ppt_summer_diff = seq(min_summer_ppt, max_summer_ppt, length.out = 50),
                   ppt_autumn_diff = seq(min_autumn_ppt, max_autumn_ppt, length.out = 50),
                   ppt_spring_sd_diff = seq(min_spring_ppt_sd, max_spring_ppt_sd, length.out = 50),
                   ppt_summer_sd_diff = seq(min_summer_ppt_sd, max_summer_ppt_sd, length.out = 50),
                   ppt_autumn_sd_diff = seq(min_autumn_ppt_sd, max_autumn_ppt_sd, length.out = 50),
                   tmean_spring_diff = seq(min_spring_tmean, max_spring_tmean, length.out = 50),
                   tmean_summer_diff = seq(min_summer_tmean, max_summer_tmean, length.out = 50),
                   tmean_autumn_diff = seq(min_autumn_tmean, max_autumn_tmean, length.out = 50),
                   tmean_spring_sd_diff = seq(min_spring_tmean_sd, max_spring_tmean_sd, length.out = 50),
                   tmean_summer_sd_diff = seq(min_summer_tmean_sd, max_summer_tmean_sd, length.out = 50),
                   tmean_autumn_sd_diff = seq(min_autumn_tmean_sd, max_autumn_tmean_sd, length.out = 50)
                   )

# gennerating predictions 

prediction <- predict(
  fit_list[[s]],
  newdata = preddata,
  formula =  ~ tibble("spring_ppt" =  Intercept + spring_ppt,
                      "summer_ppt" =  Intercept + summer_ppt,
                      "autumn_ppt" =  Intercept + autumn_ppt,
                      # 
                      "spring_ppt_sd" =  Intercept + spring_ppt_sd,
                      "summer_ppt_sd" =  Intercept + summer_ppt_sd,
                      "autumn_ppt_sd" =  Intercept + autumn_ppt_sd,
                      # 
                      "spring_tmean" =   Intercept + spring_tmean,
                      "summer_tmean" =   Intercept + summer_tmean,
                      "autumn_tmean" =   Intercept + autumn_tmean,
                      # 
                      "spring_tmean_sd" = Intercept + spring_tmean_sd,
                      "summer_tmean_sd" = Intercept + summer_tmean_sd,
                      "autumn_tmean_sd" = Intercept + autumn_tmean_sd),
  probs = c(0.025, 0.25, 0.5, 0.75, 0.975),
  n.samples = 100)
prediction_list[[species_names[s]]] <- prediction


posterior <- generate(
  fit_list[[s]],
  newdata = preddata,
  formula =  ~ c("spring_ppt" =  spring_ppt_latent,
                      "summer_ppt" =  summer_ppt_latent,
                      "autumn_ppt" =  autumn_ppt_latent,
                      # 
                      "spring_ppt_sd" =  spring_ppt_sd_latent,
                      "summer_ppt_sd" =  summer_ppt_sd_latent,
                      "autumn_ppt_sd" =  autumn_ppt_sd_latent,
                      # 
                      "spring_tmean" =   spring_tmean_latent,
                      "summer_tmean" =   summer_tmean_latent,
                      "autumn_tmean" =   autumn_tmean_latent,
                      # 
                      "spring_tmean_sd" = spring_tmean_sd_latent,
                      "summer_tmean_sd" = summer_tmean_sd_latent,
                      "autumn_tmean_sd" = autumn_tmean_sd_latent),
  n.samples = 500)
posterior_list[[species_names[s]]] <- posterior
}



prediction_df_aghy <-  prediction_list[[1]] %>% map_dfr(~ .x %>% as_tibble(), .id = "variable")
prediction_df_agpe <-  prediction_list[[2]] %>% map_dfr(~ .x %>% as_tibble(), .id = "variable")
prediction_df_elvi <-  prediction_list[[3]] %>% map_dfr(~ .x %>% as_tibble(), .id = "variable")

prediction_df_setup <- bind_rows(prediction_df_aghy, prediction_df_agpe, prediction_df_elvi) %>% 
  pivot_longer(cols =  c(ppt_spring_diff, ppt_summer_diff, ppt_autumn_diff, ppt_spring_sd_diff, ppt_summer_sd_diff, ppt_autumn_sd_diff, tmean_spring_diff, tmean_summer_diff, tmean_autumn_diff ,tmean_spring_sd_diff, tmean_summer_sd_diff ,tmean_autumn_sd_diff)) %>% 
  filter(variable == "spring_ppt" & name == "ppt_spring_diff"|
         variable == "summer_ppt" & name == "ppt_summer_diff"|
         variable == "autumn_ppt" & name == "ppt_autumn_diff"|
         variable == "spring_ppt_sd" & name == "ppt_spring_sd_diff"|
         variable == "summer_ppt_sd" & name == "ppt_summer_sd_diff"|
         variable == "autumn_ppt_sd" & name == "ppt_autumn_sd_diff"|

         variable == "spring_tmean" & name == "tmean_spring_diff"|
         variable == "summer_tmean" & name == "tmean_summer_diff"|
         variable == "autumn_tmean" & name == "tmean_autumn_diff"|
         variable == "spring_tmean_sd" & name == "tmean_spring_sd_diff"|
         variable == "summer_tmean_sd" & name == "tmean_summer_sd_diff"|
         variable == "autumn_tmean_sd" & name == "tmean_autumn_sd_diff") %>%
  mutate(name = fct_relevel(name, c('tmean_spring_diff', 'tmean_summer_diff', "tmean_autumn_diff",
                                    'tmean_spring_sd_diff', 'tmean_summer_sd_diff', "tmean_autumn_sd_diff",
                                    'ppt_spring_diff', 'ppt_summer_diff', "ppt_autumn_diff",
                                    'ppt_spring_sd_diff', 'ppt_summer_sd_diff', "ppt_autumn_sd_diff"))) %>% 
  mutate(moment = case_when(grepl("sd", name) ~ "sd",
                            !grepl("sd", name) ~ "mean"),
         variable_group = case_when(grepl("tmean", variable) ~ "Temperature",
                                    grepl("ppt", variable) ~ "Precipitation"))


climate_labels <- c(
  ppt_spring_diff = "Spring", 
  ppt_summer_diff = "Summer", 
  ppt_autumn_diff = "Autumn",
  ppt_spring_sd_diff = "Spring", 
  ppt_summer_sd_diff = "Summer", 
  ppt_autumn_sd_diff = "Autumn",
  tmean_spring_diff = "Spring", 
  tmean_summer_diff = "Summer", 
  tmean_autumn_diff = "Autumn",
  tmean_spring_sd_diff = "Spring", 
  tmean_summer_sd_diff = "Summer", 
  tmean_autumn_sd_diff = "Autumn",
  mean = "Mean",
  sd = "Standard Deviation"
)


# transforming the posterior samples to dataframe to caluclate summary statistics
colnames(posterior_list[[1]]) <- colnames(posterior_list[[2]]) <- colnames(posterior_list[[3]]) <- c( paste0("iter",1:500))


aghy_post_df <- as_tibble(posterior_list[[1]]) %>% 
  mutate(variable = rownames(posterior_list[[1]]),
         species = "A. hyemalis") 
agpe_post_df <- as_tibble(posterior_list[[2]]) %>% 
  mutate(variable = rownames(posterior_list[[2]]),
         species = "A. perennans")
elvi_post_df <- as_tibble(posterior_list[[3]]) %>% 
  mutate(variable = rownames(posterior_list[[3]]),
         species = "E. virginicus")


climate_post_df <- bind_rows(aghy_post_df, agpe_post_df, elvi_post_df) %>% 
  pivot_longer(cols = -c(species, variable), names_to = "iteration") %>% 
  mutate(positive = value>0) %>% 
  group_by(species, variable) %>% 
  dplyr::summarise(post_median = median(value),
                   prob_pos = (sum(positive)/500)*100) %>% 
  mutate(overlap = case_when(prob_pos >=95 ~ "significant",
                             prob_pos <95 & prob_pos >5 ~ "flat",
                             prob_pos <=5 ~ "significant"))

prediction_df <- prediction_df_setup %>% 
  left_join(climate_post_df)


tmean_trend <- ggplot(filter(prediction_df, grepl("tmean", name) & !grepl("sd", name))) +
  geom_vline(aes(xintercept = 0), color = "grey80")+geom_hline(aes(yintercept = 0), color = "grey80")+
  geom_point(data = filter(svc.pred_climate_subsample_long, grepl("tmean", name) & !grepl("annual", name) & !grepl("sd", name)), aes(x = value, y = mean, color = species), alpha = 0.2)+
  geom_line(aes(value, mean, linetype = overlap)) +
  geom_ribbon(aes(value, ymin = q0.025, ymax = q0.975, fill = species), alpha = 0.2) +
  geom_ribbon(aes(value, ymin = q0.25, ymax = q0.75, fill = species), alpha = 0.2) +
  facet_grid(species ~ name,
             labeller = labeller(name =climate_labels),
             scales = "free_x")+
  scale_color_manual(values = species_colors)+
  scale_fill_manual(values = species_colors)+
  scale_linetype_manual(values = c("dashed", "solid") )+
  guides(fill = "none", color = "none", linetype = "none")+
  labs(y = "Change in Prevalence (% per Year)", x = "Change in Temperature (ºC)")+
  theme_light()+
  theme(strip.background = element_blank(),
        strip.text = element_text(color = "black"),
        strip.text.x = element_text( size = rel(2)),
        strip.text.y = element_blank(),
        axis.text = element_text(size = rel(1)),
        axis.title = element_text(size = rel(1.5)),
        legend.text = element_text(face = "italic"),
        )
  # lims(y = c(0,1))

tmean_trend


tmean_sd_trend <- ggplot(filter(prediction_df, grepl("tmean", name) & grepl("sd", name))) +
  geom_vline(aes(xintercept = 0), color = "grey80")+geom_hline(aes(yintercept = 0), color = "grey80")+
  geom_point(data = filter(svc.pred_climate_subsample_long, grepl("tmean", name) & !grepl("annual", name) & grepl("sd", name)), aes(x = value, y = mean, color = species), alpha = 0.2)+
  geom_line(aes(value, mean, linetype = overlap)) +
  geom_ribbon(aes(value, ymin = q0.025, ymax = q0.975, fill = species), alpha = 0.2) +
  geom_ribbon(aes(value, ymin = q0.25, ymax = q0.75, fill = species), alpha = 0.2) +
  facet_grid(species ~ name,
             labeller = labeller(name =climate_labels),
             scales = "free_x")+
  scale_color_manual(values = species_colors)+
  scale_fill_manual(values = species_colors)+
  scale_linetype_manual(values = c("dashed", "solid") )+
  guides(fill = "none", color = "none", linetype = "none")+
  labs(y = "Change in Prevalence (% per Year)", x = "Change in SD(Temperature) (ºC)")+
  theme_light()+
  theme(strip.background = element_blank(),
        strip.text = element_text(color = "black"),
        strip.text.x = element_text( size = rel(2)),
        strip.text.y = element_blank(),
        axis.text = element_text(size = rel(1)),
        axis.title = element_text(size = rel(1.5)),
        legend.text = element_text(face = "italic"),
  )
# lims(y = c(0,1))



# tmean_sd_trend



ppt_trend <- ggplot(filter(prediction_df, grepl("ppt", name) & !grepl("sd", name))) +
  geom_vline(aes(xintercept = 0), color = "grey80")+geom_hline(aes(yintercept = 0), color = "grey80")+
  geom_point(data = filter(svc.pred_climate_subsample_long, grepl("ppt", name) & !grepl("annual", name) & !grepl("sd", name)), aes(x = value, y = mean, color = species), alpha = 0.2)+
  geom_line(aes(value, mean, linetype = overlap)) +
  geom_ribbon(aes(value, ymin = q0.025, ymax = q0.975, fill = species), alpha = 0.2) +
  geom_ribbon(aes(value, ymin = q0.25, ymax = q0.75, fill = species), alpha = 0.2) +
  facet_grid(species ~ name,
             labeller = labeller(name =climate_labels),
             scales = "free_x")+
  scale_color_manual(values = species_colors)+
  scale_fill_manual(values = species_colors)+
  scale_linetype_manual(values = c("dashed", "solid") )+
  guides(fill = "none", color = "none", linetype = "none")+
  labs(y = "Change in Prevalence (% per Year)", x = "Change in Precipitation (mm.)")+
  theme_light()+
  theme(strip.background = element_blank(),
        strip.text = element_text(color = "black"),
        strip.text.x = element_text( size = rel(2)),
        strip.text.y = element_text(face = "italic", size = rel(2), angle = 0),
        axis.text = element_text(size = rel(1)),
        axis.title = element_text(size = rel(1.5)),
        legend.text = element_text(face = "italic"),
  )
# lims(y = c(0,1))

# ppt_trend


ppt_sd_trend <- ggplot(filter(prediction_df, grepl("ppt", name) & grepl("sd", name))) +
  geom_vline(aes(xintercept = 0), color = "grey80")+geom_hline(aes(yintercept = 0), color = "grey80")+
  geom_point(data = filter(svc.pred_climate_subsample_long, grepl("ppt", name) & !grepl("annual", name) & grepl("sd", name)), aes(x = value, y = mean, color = species), alpha = 0.2)+
  geom_line(aes(value, mean, linetype = overlap)) +
  geom_ribbon(aes(value, ymin = q0.025, ymax = q0.975, fill = species), alpha = 0.2) +
  geom_ribbon(aes(value, ymin = q0.25, ymax = q0.75, fill = species), alpha = 0.2) +
  facet_grid(species ~ name,
             labeller = labeller(name =climate_labels),
             scales = "free_x")+
  scale_color_manual(values = species_colors)+
  scale_fill_manual(values = species_colors)+
  scale_linetype_manual(values = c("dashed", "solid") )+
  guides(fill = "none", color = "none", linetype = "none")+
  labs(y = "Change in Prevalence (% per Year)", x = "Change in SD(Precipitation) (mm.)")+
  theme_light()+
  theme(strip.background = element_blank(),
        strip.text = element_text(color = "black"),
        strip.text.x = element_text( size = rel(2)),
        strip.text.y = element_text(size = rel(2), face = "italic", angle = 0),
        axis.text = element_text(size = rel(1)),
        axis.title = element_text(size = rel(1.5)),
        legend.text = element_text(face = "italic"),
  )
# lims(y = c(0,1))

# ppt_sd_trend





# tmean_trend$layers[4:6] <- NULL ; tmean_trend$layers[[3]]$aes_params$alpha  <-0
# tmean_sd_trend$layers[4:6] <- NULL; tmean_sd_trend$layers[[3]]$aes_params$alpha  <-0
# ppt_trend$layers[4:6] <- NULL; ppt_trend$layers[[3]]$aes_params$alpha  <-0
# ppt_sd_trend$layers[4:6] <- NULL; ppt_sd_trend$layers[[3]]$aes_params$alpha  <-0

  
climate_trends_plot <- (tmean_trend + ppt_trend) /( tmean_sd_trend + ppt_sd_trend) + plot_annotation(tag_levels = "A") +  theme(plot.tag = element_text(face = 'bold'))

ggsave(climate_trends_plot, filename = "Plots/climate_trends_plot_intercept.png", width = 14, height = 10)

######### version of plot separating by species ####







aghy_tmean_trend <- ggplot(filter(prediction_df, grepl(species_names[1], species)  & grepl("tmean", name) & !grepl("sd", moment))) +
  geom_vline(aes(xintercept = 0), color = "grey80")+geom_hline(aes(yintercept = 0), color = "grey80")+
  geom_point(data = filter(svc.pred_climate_subsample_long, species == "A. hyemalis" & grepl("tmean", name) & !grepl("sd", name) & !grepl("annual", name)), aes(x = value, y = mean, color = species), alpha = 0.2)+
  geom_line(aes(value, mean, linetype = overlap)) +
  geom_ribbon(aes(value, ymin = q0.025, ymax = q0.975, fill = species), alpha = 0.2) +
  geom_ribbon(aes(value, ymin = q0.25, ymax = q0.75, fill = species), alpha = 0.2) +
  facet_grid(  moment~ name,
               labeller = labeller(name =climate_labels,
                                   moment = climate_labels),
               scales = "free_x")+
  scale_color_manual(values = species_colors)+
  scale_fill_manual(values = species_colors)+
  scale_linetype_manual(values = c("dashed", "solid") )+
  guides(fill = "none", color = "none", linetype = "none")+ 
  labs(y = "Change in Prevalence (% per Year)", x = "Change in Temperature (ºC)")+
  theme_light()+
  theme(strip.background = element_blank(),
        strip.text = element_text(color = "black", size = rel(.8)),
        strip.text.y = element_blank(),
        legend.text = element_text(face = "italic"),
  )
# lims(y = c(0,1))

# aghy_tmean_trend

aghy_tmean_sd_trend <- ggplot(filter(prediction_df, grepl(species_names[1], species)  & grepl("tmean", name) & grepl("sd", moment))) +
  geom_vline(aes(xintercept = 0), color = "grey80")+geom_hline(aes(yintercept = 0), color = "grey80")+
  geom_point(data = filter(svc.pred_climate_subsample_long,species == "A. hyemalis" & grepl("tmean", name) & grepl("sd", name) & !grepl("annual", name)), aes(x = value, y = mean, color = species), alpha = 0.2)+
  geom_line(aes(value, mean, linetype = overlap)) +
  geom_ribbon(aes(value, ymin = q0.025, ymax = q0.975, fill = species), alpha = 0.2) +
  geom_ribbon(aes(value, ymin = q0.25, ymax = q0.75, fill = species), alpha = 0.2) +
  facet_grid(  moment~ name,
               labeller = labeller(name =climate_labels,
                                   moment = climate_labels),
               scales = "free_x")+
  scale_color_manual(values = species_colors)+
  scale_fill_manual(values = species_colors)+
  scale_linetype_manual(values = c("dashed", "solid") )+
  guides(fill = "none", color = "none", linetype = "none")+ 
  labs(y = "Change in Prevalence (% per Year)", x = "Change in SD(Temperature) (mm.)")+
  theme_light()+
  theme(strip.background = element_blank(),
        strip.text = element_text(color = "black", size = rel(.8)),
        strip.text.y = element_blank(),
        legend.text = element_text(face = "italic"),
  )
# lims(y = c(0,1))

# aghy_tmean_sd_trend


aghy_ppt_trend <- ggplot(filter(prediction_df, grepl(species_names[1], species)  & grepl("ppt", name) & !grepl("summer", name) & !grepl("sd", moment))) +
  geom_vline(aes(xintercept = 0), color = "grey80")+geom_hline(aes(yintercept = 0), color = "grey80")+
  geom_point(data = filter(svc.pred_climate_subsample_long, grepl(species_names[1], species)  & grepl("ppt", name) & !grepl("summer", name) & !grepl("sd", name) & !grepl("annual", name)), aes(x = value, y = mean, color = species), alpha = 0.2)+
  geom_line(aes(value, mean, linetype = overlap)) +
  geom_ribbon(aes(value, ymin = q0.025, ymax = q0.975, fill = species), alpha = 0.2) +
  geom_ribbon(aes(value, ymin = q0.25, ymax = q0.75, fill = species), alpha = 0.2) +
  facet_grid(  moment~ name,
               labeller = labeller(name =climate_labels,
                                   moment = climate_labels),
               scales = "free_x")+
  scale_color_manual(values = species_colors)+
  scale_fill_manual(values = species_colors)+
  scale_linetype_manual(values = c("dashed", "solid") )+
  guides(fill = "none", color = "none", linetype = "none")+ 
  labs(y = "Change in Prevalence (% per Year)", x = "Change in Precipitation (mm.)")+
  theme_light()+
  theme(strip.background = element_blank(),
        strip.text = element_text(color = "black", size = rel(.8)),
        strip.text.y = element_blank(),
        legend.text = element_text(face = "italic"),
  )
# lims(y = c(0,1))
ggsave(aghy_ppt_trend, filename = "aghy_ppt_trend.png", width = 9, height = 4)
# aghy_ppt_trend

aghy_ppt_sd_trend <- ggplot(filter(prediction_df, grepl(species_names[1], species)  & grepl("ppt", name) & grepl("sd", moment))) +
  geom_vline(aes(xintercept = 0), color = "grey80")+geom_hline(aes(yintercept = 0), color = "grey80")+
  geom_point(data = filter(svc.pred_climate_subsample_long, species == "A. hyemalis" & grepl("ppt", name) & grepl("sd", name) & !grepl("annual", name)), aes(x = value, y = mean, color = species), alpha = 0.2)+
  geom_line(aes(value, mean, linetype = overlap)) +
  geom_ribbon(aes(value, ymin = q0.025, ymax = q0.975, fill = species), alpha = 0.2) +
  geom_ribbon(aes(value, ymin = q0.25, ymax = q0.75, fill = species), alpha = 0.2) +
  facet_grid(  moment~ name,
               labeller = labeller(name =climate_labels,
                                   moment = climate_labels),
               scales = "free_x")+
  scale_color_manual(values = species_colors)+
  scale_fill_manual(values = species_colors)+
  scale_linetype_manual(values = c("dashed", "solid") )+
  guides(fill = "none", color = "none", linetype = "none")+ 
  labs(y = "Change in Prevalence (% per Year)", x = "Change in SD(Precipitation) (mm.)")+
  theme_light()+
  theme(strip.background = element_blank(),
        strip.text = element_text(color = "black", size = rel(.8)),
        strip.text.y = element_blank(),
        legend.text = element_text(face = "italic"),
  )
# lims(y = c(0,1))

# aghy_ppt_sd_trend 

# aghy_tmean_trend$layers[4:6] <- NULL ; aghy_tmean_trend$layers[[3]]$aes_params$alpha  <-0
# aghy_tmean_sd_trend$layers[4:6] <- NULL; aghy_tmean_sd_trend$layers[[3]]$aes_params$alpha  <-0
# aghy_ppt_trend$layers[4:6] <- NULL; aghy_ppt_trend$layers[[3]]$aes_params$alpha  <-0
# aghy_ppt_sd_trend$layers[4:6] <- NULL; aghy_ppt_sd_trend$layers[[3]]$aes_params$alpha  <-0


AGHY_climate_trends_plot <- (aghy_tmean_trend + aghy_ppt_trend ) / (aghy_tmean_sd_trend+aghy_ppt_sd_trend)

ggsave(AGHY_climate_trends_plot, filename = "Plots/AGHY_climate_trends_plot.png", width = 12, height = 8)

# agpe plots
agpe_tmean_trend <- ggplot(filter(prediction_df, grepl(species_names[2], species)  & grepl("tmean", name) & grepl("autumn", name) & !grepl("sd", moment))) +
  geom_vline(aes(xintercept = 0), color = "grey80")+geom_hline(aes(yintercept = 0), color = "grey80")+
  geom_point(data = filter(svc.pred_climate_subsample_long, species == "A. perennans" & grepl("tmean", name) & grepl("autumn", name) & !grepl("sd", name) & !grepl("annual", name)), aes(x = value, y = mean, color = species), alpha = 0.2)+
  geom_line(aes(value, mean, linetype = overlap)) +
  geom_ribbon(aes(value, ymin = q0.025, ymax = q0.975, fill = species), alpha = 0.2) +
  geom_ribbon(aes(value, ymin = q0.25, ymax = q0.75, fill = species), alpha = 0.2) +
  facet_grid(  moment~ name,
               labeller = labeller(name =climate_labels,
                                   moment = climate_labels),
               scales = "free_x")+
  scale_color_manual(values = species_colors[2])+
  scale_fill_manual(values = species_colors[2])+
  scale_linetype_manual(values = c("flat"="dashed","significant"="solid") )+
  guides(fill = "none", color = "none", linetype = "none")+ 
  labs(y = "Change in Prevalence (% per Year)", x = "Change in Temperature (ºC)")+
  theme_light()+
  theme(strip.background = element_blank(),
        strip.text = element_text(color = "black", size = rel(.8)),
        strip.text.y = element_blank(),
        legend.text = element_text(face = "italic"),
  )
# lims(y = c(0,1))

# agpe_tmean_trend

agpe_tmean_sd_trend <- ggplot(filter(prediction_df, grepl(species_names[2], species)  & grepl("tmean", name) & grepl("autumn", name) & grepl("sd", moment))) +
  geom_vline(aes(xintercept = 0), color = "grey80")+geom_hline(aes(yintercept = 0), color = "grey80")+
  geom_point(data = filter(svc.pred_climate_subsample_long,species == "A. perennans" & grepl("tmean", name) & grepl("autumn", name)  & grepl("sd", name) & !grepl("annual", name)), aes(x = value, y = mean, color = species), alpha = 0.2)+
  geom_line(aes(value, mean, linetype = overlap)) +
  geom_ribbon(aes(value, ymin = q0.025, ymax = q0.975, fill = species), alpha = 0.2) +
  geom_ribbon(aes(value, ymin = q0.25, ymax = q0.75, fill = species), alpha = 0.2) +
  facet_grid(  moment~ name,
               labeller = labeller(name =climate_labels,
                                   moment = climate_labels),
               scales = "free_x")+
  scale_color_manual(values = species_colors[2])+
  scale_fill_manual(values = species_colors[2])+
  scale_linetype_manual(values = c("flat"="dashed","significant"="solid") )+
  guides(fill = "none", color = "none", linetype = "none")+ 
  labs(y = "Change in Prevalence (% per Year)", x = "Change in SD(Temperature) (mm.)")+
  theme_light()+
  theme(strip.background = element_blank(),
        strip.text = element_text(color = "black", size = rel(.8)),
        strip.text.y = element_blank(),
        legend.text = element_text(face = "italic"),
  )
# lims(y = c(0,1))

agpe_tmean_sd_trend


agpe_ppt_trend <- ggplot(filter(prediction_df, grepl(species_names[2], species)  & grepl("ppt", name) & grepl("autumn", name) & !grepl("sd", moment))) +
  geom_vline(aes(xintercept = 0), color = "grey80")+geom_hline(aes(yintercept = 0), color = "grey80")+
  geom_point(data = filter(svc.pred_climate_subsample_long, grepl(species_names[2], species)  & grepl("ppt", name) & grepl("autumn", name)  & !grepl("sd", name) & !grepl("annual", name)), aes(x = value, y = mean, color = species), alpha = 0.2)+
  geom_line(aes(value, mean, linetype = overlap)) +
  geom_ribbon(aes(value, ymin = q0.025, ymax = q0.975, fill = species), alpha = 0.2) +
  geom_ribbon(aes(value, ymin = q0.25, ymax = q0.75, fill = species), alpha = 0.2) +
  facet_grid(  moment~ name,
               labeller = labeller(name =climate_labels,
                                   moment = climate_labels),
               scales = "free_x")+
  scale_color_manual(values = species_colors[2])+
  scale_fill_manual(values = species_colors[2])+
  scale_linetype_manual(values = c("flat"="dashed","significant"="solid") )+
  guides(fill = "none", color = "none", linetype = "none")+ 
  labs(y = "Change in Prevalence (% per Year)", x = "Change in Precipitation (mm.)")+
  theme_light()+
  theme(strip.background = element_blank(),
        strip.text = element_text(color = "black", size = rel(.8)),
        strip.text.y = element_blank(),
        legend.text = element_text(face = "italic"),
  )
# lims(y = c(0,1))

# agpe_ppt_trend

agpe_ppt_sd_trend <- ggplot(filter(prediction_df, grepl(species_names[2], species)  & grepl("ppt", name) & grepl("autumn", name) & grepl("sd", moment))) +
  geom_vline(aes(xintercept = 0), color = "grey80")+geom_hline(aes(yintercept = 0), color = "grey80")+
  geom_point(data = filter(svc.pred_climate_subsample_long, species == "A. perennans" & grepl("ppt", name) & grepl("autumn", name) & grepl("sd", name) & !grepl("annual", name)), aes(x = value, y = mean, color = species), alpha = 0.2)+
  geom_line(aes(value, mean, linetype = overlap)) +
  geom_ribbon(aes(value, ymin = q0.025, ymax = q0.975, fill = species), alpha = 0.2) +
  geom_ribbon(aes(value, ymin = q0.25, ymax = q0.75, fill = species), alpha = 0.2) +
  facet_grid(  moment~ name,
               labeller = labeller(name =climate_labels,
                                   moment = climate_labels),
               scales = "free_x")+
  scale_color_manual(values = species_colors[2])+
  scale_fill_manual(values = species_colors[2])+
  scale_linetype_manual(values = c("flat"="dashed","significant"="solid") )+
  guides(fill = "none", color = "none", linetype = "none")+ 
  labs(y = "Change in Prevalence (% per Year)", x = "Change in SD(Precipitation) (mm.)")+
  theme_light()+
  theme(strip.background = element_blank(),
        strip.text = element_text(color = "black", size = rel(.8)),
        strip.text.y = element_blank(),
        legend.text = element_text(face = "italic"),
  )
# lims(y = c(0,1))

# agpe_ppt_sd_trend 

# agpe_tmean_trend$layers[4:6] <- NULL ; agpe_tmean_trend$layers[[3]]$aes_params$alpha  <-0
# agpe_tmean_sd_trend$layers[4:6] <- NULL; agpe_tmean_sd_trend$layers[[3]]$aes_params$alpha  <-0
# agpe_ppt_trend$layers[4:6] <- NULL; agpe_ppt_trend$layers[[3]]$aes_params$alpha  <-0
# agpe_ppt_sd_trend$layers[4:6] <- NULL; agpe_ppt_sd_trend$layers[[3]]$aes_params$alpha  <-0


AGPE_climate_trends_plot <- (agpe_tmean_trend + agpe_ppt_trend ) #/ (agpe_tmean_sd_trend+agpe_ppt_sd_trend)

ggsave(AGPE_climate_trends_plot, filename = "AGPE_climate_trends_plot.png", width = 9, height = 4)

# elvi plots
elvi_tmean_trend <- ggplot(filter(prediction_df, grepl(species_names[3], species)  & grepl("tmean", name) & grepl("summer", name) & !grepl("sd", moment))) +
  geom_vline(aes(xintercept = 0), color = "grey80")+geom_hline(aes(yintercept = 0), color = "grey80")+
  geom_point(data = filter(svc.pred_climate_subsample_long, species == "E. virginicus" & grepl("tmean", name) & grepl("summer", name) & !grepl("sd", name) & !grepl("annual", name)), aes(x = value, y = mean, color = species), alpha = 0.2)+
  geom_line(aes(value, mean, linetype = overlap)) +
  geom_ribbon(aes(value, ymin = q0.025, ymax = q0.975, fill = species), alpha = 0.2) +
  geom_ribbon(aes(value, ymin = q0.25, ymax = q0.75, fill = species), alpha = 0.2) +
  facet_grid(  moment~ name,
               labeller = labeller(name =climate_labels,
                                   moment = climate_labels),
               scales = "free_x")+
  scale_color_manual(values = species_colors[3])+
  scale_fill_manual(values = species_colors[3])+
  scale_linetype_manual(values = c("flat"="dashed","significant"="solid") )+
  guides(fill = "none", color = "none", linetype = "none")+ 
  labs(y = "Change in Prevalence (% per Year)", x = "Change in Temperature (ºC)")+
  theme_light()+
  theme(strip.background = element_blank(),
        strip.text = element_text(color = "black", size = rel(.8)),
        strip.text.y = element_blank(),
        legend.text = element_text(face = "italic"),
  )
# lims(y = c(0,1))

# elvi_tmean_trend

elvi_tmean_sd_trend <- ggplot(filter(prediction_df, grepl(species_names[3], species)  & grepl("tmean", name) & grepl("summer", name) & grepl("sd", moment))) +
  geom_vline(aes(xintercept = 0), color = "grey80")+geom_hline(aes(yintercept = 0), color = "grey80")+
  geom_point(data = filter(svc.pred_climate_subsample_long,species == "E. virginicus" & grepl("tmean", name) & grepl("summer", name) & grepl("sd", name) & !grepl("annual", name)), aes(x = value, y = mean, color = species), alpha = 0.2)+
  geom_line(aes(value, mean, linetype = overlap)) +
  geom_ribbon(aes(value, ymin = q0.025, ymax = q0.975, fill = species), alpha = 0.2) +
  geom_ribbon(aes(value, ymin = q0.25, ymax = q0.75, fill = species), alpha = 0.2) +
  facet_grid(  moment~ name,
               labeller = labeller(name =climate_labels,
                                   moment = climate_labels),
               scales = "free_x")+
  scale_color_manual(values = species_colors[3])+
  scale_fill_manual(values = species_colors[3])+
  scale_linetype_manual(values = c("flat"="dashed","significant"="solid") )+
  guides(fill = "none", color = "none", linetype = "none")+ 
  labs(y = "Change in Prevalence (% per Year)", x = "Change in SD(Temperature) (mm.)")+
  theme_light()+
  theme(strip.background = element_blank(),
        strip.text = element_text(color = "black", size = rel(.8)),
        strip.text.y = element_blank(),
        legend.text = element_text(face = "italic"),
  )
# lims(y = c(0,1))

# elvi_tmean_sd_trend


elvi_ppt_trend <- ggplot(filter(prediction_df, grepl(species_names[3], species)  & grepl("ppt", name) & grepl("summer", name) & !grepl("sd", moment))) +
  geom_vline(aes(xintercept = 0), color = "grey80")+geom_hline(aes(yintercept = 0), color = "grey80")+
  geom_point(data = filter(svc.pred_climate_subsample_long, grepl(species_names[3], species)  & grepl("ppt", name) & grepl("summer", name) & !grepl("sd", name) & !grepl("annual", name)), aes(x = value, y = mean, color = species), alpha = 0.2)+
  geom_line(aes(value, mean, linetype = overlap)) +
  geom_ribbon(aes(value, ymin = q0.025, ymax = q0.975, fill = species), alpha = 0.2) +
  geom_ribbon(aes(value, ymin = q0.25, ymax = q0.75, fill = species), alpha = 0.2) +
  facet_grid(  moment~ name,
               labeller = labeller(name =climate_labels,
                                   moment = climate_labels),
               scales = "free_x")+
  scale_color_manual(values = species_colors[3])+
  scale_fill_manual(values = species_colors[3])+
  scale_linetype_manual(values = c("flat"="dashed","significant"="solid") )+
  guides(fill = "none", color = "none", linetype = "none")+ 
  labs(y = "Change in Prevalence (% per Year)", x = "Change in Precipitation (mm.)")+
  theme_light()+
  theme(strip.background = element_blank(),
        strip.text = element_text(color = "black", size = rel(.8)),
        strip.text.y = element_blank(),
        legend.text = element_text(face = "italic"),
  )
# lims(y = c(0,1))

# elvi_ppt_trend

elvi_ppt_sd_trend <- ggplot(filter(prediction_df, grepl(species_names[3], species)  & grepl("ppt", name) & grepl("autumn", name) & grepl("sd", moment))) +
  geom_vline(aes(xintercept = 0), color = "grey80")+geom_hline(aes(yintercept = 0), color = "grey80")+
  geom_point(data = filter(svc.pred_climate_subsample_long, species == "E. virginicus" & grepl("ppt", name) & grepl("autumn", name) & grepl("sd", name) & !grepl("annual", name)), aes(x = value, y = mean, color = species), alpha = 0.2)+
  geom_line(aes(value, mean, linetype = overlap)) +
  geom_ribbon(aes(value, ymin = q0.025, ymax = q0.975, fill = species), alpha = 0.2) +
  geom_ribbon(aes(value, ymin = q0.25, ymax = q0.75, fill = species), alpha = 0.2) +
  facet_grid(  moment~ name,
               labeller = labeller(name =climate_labels,
                                   moment = climate_labels),
               scales = "free_x")+
  scale_color_manual(values = species_colors[3])+
  scale_fill_manual(values = species_colors[3])+
  scale_linetype_manual(values = c("flat"="dashed","significant"="solid") )+
  guides(fill = "none", color = "none", linetype = "none")+ 
  labs(y = "Change in Prevalence (% per Year)", x = "Change in SD(Precipitation) (mm.)")+
  theme_light()+
  theme(strip.background = element_blank(),
        strip.text = element_text(color = "black", size = rel(.8)),
        strip.text.y = element_blank(),
        legend.text = element_text(face = "italic"),
  )
# lims(y = c(0,1))

# elvi_ppt_sd_trend 

# elvi_tmean_trend$layers[4:6] <- NULL ; elvi_tmean_trend$layers[[3]]$aes_params$alpha  <-0
# elvi_tmean_sd_trend$layers[4:6] <- NULL; elvi_tmean_sd_trend$layers[[3]]$aes_params$alpha  <-0
# elvi_ppt_trend$layers[4:6] <- NULL; elvi_ppt_trend$layers[[3]]$aes_params$alpha  <-0
# elvi_ppt_sd_trend$layers[4:6] <- NULL; elvi_ppt_sd_trend$layers[[3]]$aes_params$alpha  <-0


ELVI_climate_trends_plot <- (elvi_tmean_trend + elvi_ppt_trend ) / (elvi_tmean_sd_trend+elvi_ppt_sd_trend)

ggsave(ELVI_climate_trends_plot, filename = "Plots/ELVI_climate_trends_plot.png", width = 9, height = 5)






###### generating some summaries for the paper ##############





###### Plotting the posterior distributions ###################################################

plot(fit, varname = "spring_ppt")

flist <- vector("list", NROW(fit$summary.fixed$ppt_spring))
for (i in seq_along(flist)) flist[[i]] <- plot(fit, "spring_ppt", index = i) +geom_vline(aes(xintercept = 0))
multiplot(plotlist = flist, cols = 3)








