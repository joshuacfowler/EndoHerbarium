# Purpose: Fits spatiotemporal model of change in endophyte prevalence in herbarium specimens
# Authors: Joshua Fowler and Tom Miller
# Updated: Feb 22, 2024

library(tidyverse) # for data manipulation and ggplot
library(INLA) # for fitting integrated nested Laplace approximation models
library(inlabru)
library(fmesher)

library(sf)
# library(rmapshaper)
library(terra)
library(tidyterra)
# library(ggpolypath)



library(patchwork)
library(ggmap)
library(pROC)

invlogit<-function(x){exp(x)/(1+exp(x))}
species_colors <- c("#1b9e77","#d95f02","#7570b3")
endophyte_colors <- c("#fdedd3","#f3c8a8", "#5a727b", "#4986c7", "#181914",  "#163381")


################################################################################
############ Read in the herbarium dataset ############################### 
################################################################################

endo_herb_georef <- read_csv(file = "~/Dropbox/Josh&Tom - shared/Endo_Herbarium/DigitizedHerbariumRecords/endo_herb_georef.csv") %>% 
  # filter(Country != "Canada") %>%
  mutate(sample_temp = Sample_id) %>% 
  separate(sample_temp, into = c("Herb_code", "spp_code", "specimen_code", "tissue_code")) %>% 
  mutate(species_index = as.factor(case_when(spp_code == "AGHY" ~ "1",
                                             spp_code == "AGPE" ~ "2",
                                             spp_code == "ELVI" ~ "3"))) %>% 
  mutate(species = case_when(spp_code == "AGHY" ~ "A. hyemalis",
                             spp_code == "AGPE" ~ "A. perennans",
                             spp_code == "ELVI" ~ "E. virginicus")) %>% 
  mutate(decade_bin = floor(year/10)*10) %>% 
  mutate(std_year = (year-mean(year, na.rm = T))) %>%  # I am mean centering but not scaling by standard deviation to preserve units for interpretation of the parameter values
  filter(scorer_factor != "Scorer26")
# mini_dataset <- endo_herb_georef %>% 
#   filter(year>1970 &year<1990 & species == "A. hyemalis") %>% 
#   mutate(presence = Endo_status_liberal) %>% 
#   filter(!is.na(Endo_status_liberal)) %>%
#   filter(!is.na(lon) & !is.na(year)) %>% 
#   filter(lon>-110 ) %>% 
#   filter(Country != "Canada" ) %>% 
#   dplyr::select(Sample_id, lat, lon, year, std_year, presence)
# 
# write_csv(mini_dataset, "2024-04-02_endophyte_data.csv")

endo_herb <- endo_herb_georef %>% 
  filter(!is.na(Endo_status_liberal)) %>%
  filter(!is.na(spp_code)) %>% 
  filter(!is.na(lon) & !is.na(year)) %>% 
  filter(lon>-110 ) %>% 
  filter(Country != "Canada" ) %>% 
  mutate(year_bin = case_when(year<1970 ~ "pre-1970",
                              year>=1970 ~ "post-1970")) %>% 
  mutate(scorer_id = case_when(scorer_id == "BellaGuttierez" ~ "BellaGutierrez",
                               TRUE ~ scorer_id)) %>% 
  mutate(endo_status_text = case_when(Endo_status_liberal == 0 ~ "E-",
                                      Endo_status_liberal == 1 ~ "E+"))  %>% 
  filter(spp_code %in% c("AGHY", "AGPE", "ELVI")) 

# converting the lat long to epsg 6703km in km
# define a crs
epsg6703km <- paste(
  "+proj=aea +lat_0=23 +lon_0=-96 +lat_1=29.5",
  "+lat_2=45.5 +x_0=0 +y_0=0 +datum=NAD83",
  "+units=km +no_defs"
)

endo_herb<- endo_herb %>% 
  st_as_sf(coords = c("lon", "lat"), crs = 4326, remove = FALSE) %>% 
  st_transform(epsg6703km) %>% 
  mutate(
    easting = st_coordinates(.)[, 1],
    northing = st_coordinates(.)[, 2]
  ) 



# Loading in contemporary survey data for AGHY and ELVI, which we will use for model validation
contemp_surveys <- read_csv(file = "contemp_surveys.csv") %>% 
  mutate(decade_bin = floor(Year/10)*10) %>% 
  mutate(std_year =   (Year-mean(endo_herb$year, na.rm = T))) %>% 
  # filter(SpeciesID == "AGHY") %>% 
  st_as_sf(coords = c("lon", "lat"), crs = 4326, remove = FALSE) %>% 
  st_transform(epsg6703km) %>% 
  mutate(
    easting = st_coordinates(.)[, 1],
    northing = st_coordinates(.)[, 2]
  ) 
  

contemp_random_sample <- read_csv(file = "contemp_random_sample.csv") %>% 
  st_as_sf(coords = c("lon", "lat"), crs = 4326, remove = FALSE) %>% 
  st_transform(epsg6703km) %>% 
  mutate(
    easting = st_coordinates(.)[, 1],
    northing = st_coordinates(.)[, 2]
  ) 


# register_google(key = "")
# map <- ggmap::get_map(zoom = 3, source = "google", maptype = c("satellite"))

outline_map <- map_data("world")
states_shape <- map_data("state")
counties_shape <- map_data("county")
# ggplot()+
#   geom_map(data = counties_shape, map = counties_shape, aes(long, lat, map_id = region))

collections_map <- ggplot()+
  geom_map(data = outline_map, map = outline_map, aes(map_id = region), color = "grey", linewidth = .1, fill = "#FAF9F6")+
  geom_map(data = states_shape, map = states_shape, aes(map_id = region), color = "grey", linewidth = .1, fill = NA)+
  geom_point(data = endo_herb, aes(x = lon, y = lat, color = species), alpha = .7, size = 1.2)+
  coord_sf(xlim = c(-109,-68), ylim = c(21,49))+
  scale_color_manual(values = species_colors)+
  theme_light()+
  theme(legend.text = element_text(face = "italic"))+
  labs(x = "Longitude", y = "Latitude", color = "Host Species")

# collections_map
# ggsave(collections_map, filename = "collections_map.png", width = 7, height = 4)

endo_status_map <- ggplot()+
  geom_map(data = outline_map, map = outline_map, aes( map_id = region), color = "grey", linewidth = .1, fill = "#FAF9F6")+
  geom_map(data = states_shape, map = states_shape, aes(map_id = region), color = "grey", linewidth = .1, fill = NA)+
  geom_point(data = endo_herb, aes(x = lon, y = lat, color = endo_status_text), alpha = .7, size = .8)+
  facet_wrap(~factor(year_bin, levels = c("pre-1970", "post-1970"))+species)+
  coord_sf(xlim = c(-109,-68), ylim = c(21,49))+
  scale_color_manual(values = c(endophyte_colors[2],endophyte_colors[6]))+
  theme_light()+
  theme(strip.background = element_blank(),
        strip.text  = element_text(face = "italic", color = "black"))+
  labs(x = "Longitude", y = "Latitude", color = "Endophyte Status")
# endo_status_map

# ggsave(endo_status_map, filename = "endo_status_map.png", width = 10, height = 5)



scorer_map <- ggplot()+
  geom_map(data = outline_map, map = outline_map, aes(map_id = region), color = "grey", linewidth = .1, fill = "#FAF9F6")+
  geom_map(data = states_shape, map = states_shape, aes(map_id = region), color = "grey", linewidth = .1, fill = NA)+
  geom_point(data = endo_herb, aes(x = lon, y = lat, color = scorer_factor, shape = species), alpha = .7, size = 1.2)+
  coord_sf(xlim = c(-109,-68), ylim = c(21,49))+
  theme_light()+
  theme(legend.text = element_text(face = "italic"))+
  labs(x = "Longitude", y = "Latitude", color = "Host Species")

# scorer_map
# ggsave(scorer_map, filename = "scorer_map.png", width = 7, height = 4)



# testing if the trends are still are present without oldest samples
# endo_herb <- endo_herb_georef %>% 
#   filter(!is.na(Endo_status_liberal_1)) %>%
#   filter(!is.na(Spp_code)) %>% 
#   filter(!is.na(lon) & !is.na(year)) %>% 
#   filter(lon>-110) %>% 
#   filter(year>1910)

# summary_endo_herb <- endo_herb %>% 
#   mutate(score_match = case_when(Endo_status_conservative == Endo_status_liberal ~ 1,
#                                  Endo_status_conservative != Endo_status_liberal ~ 0)) %>% 
#   group_by(Spp_code) %>%
#   summarize(n  = n(),
#             n_match = sum(score_match),
#             percent_match = n_match/n)
# 
mode <- function(codes){
  which.max(tabulate(codes))
}

summary_endo_herb <- endo_herb %>% 
  mutate(Sample_id_temp = Sample_id) %>% 
  filter(score_number == 1) %>%
  separate(Sample_id_temp, sep = "_", into = c("herbarium", "spp_code", "plant_no")) %>% select(-spp_code, -plant_no) %>% 
  # filter(seed_scored>0) %>% 
  filter(month<=12&month>0) %>% 
  group_by(species) %>% 
  summarize(n(),
            avg_seed = mean(seed_scored, na.rm = T),
            avg_month = mode(as.numeric(month)),
            min_year = min(year),
            max_year = max(year))



# making separate dataframes for each species

endo_herb_AGHY <- endo_herb %>% 
  filter(!is.na(Endo_status_liberal)) %>% 
  filter(spp_code == "AGHY") %>% 
  filter(!is.na(lon) & !is.na(year))

endo_herb_AGPE <- endo_herb %>% 
  filter(!is.na(Endo_status_liberal)) %>% 
  filter(spp_code == "AGPE") %>% 
  filter(!is.na(lon) & !is.na(year)) %>% 
  filter(lat>=20) # dropping one plant in mexico


endo_herb_ELVI <- endo_herb %>% 
  filter(!is.na(Endo_status_liberal)) %>% 
  filter(spp_code == "ELVI") %>% 
  filter(!is.na(lon) & !is.na(year)) 

endo_herb_list <- list()
endo_herb_list[[1]] <- endo_herb_AGHY
endo_herb_list[[2]] <- endo_herb_AGPE
endo_herb_list[[3]] <- endo_herb_ELVI


##########################################################################################
############ Setting up and running INLA model with inlabru ############################### 
##########################################################################################
# Running the model for each species separately
species_codes <- c("AGHY", "AGPE", "ELVI")
species_names <- c("A. hyemalis", "A. perennans", "E. virginicus")

fit_lists <- list()
mesh_list <- list()
bdry_polygon_list <- list()

for(s in 1:3){
##### Building a spatial mesh #####

data <- endo_herb_list[[s]]

# Build the spatial mesh from the coords for each species and a boundary around each species predicted distribution (eventually from Jacob's work ev)
coords <- cbind(data$easting, data$northing)

# non_convex_bdry <- inla.nonconvex.hull(coords, -0.03, -0.1, resolution = c(100, 100))


non_convex_bdry <- fm_extensions(
  data$geometry,
  convex = c(250, 500),
  concave = c(250, 500),
  crs = fm_crs(data)
)
# 
# plot(non_convex_bdry[[2]])


coastline <- st_make_valid(sf::st_as_sf(maps::map("world", regions = c("usa", "canada", "mexico"), plot = FALSE, fill = TRUE))) %>% st_transform(epsg6703km)
# plot(coastline)




bdry <- st_intersection(coastline$geom, non_convex_bdry[[1]])

# plot(bdry)

bdry_polygon <- st_cast(st_sf(bdry), "MULTIPOLYGON", group_or_split = TRUE) %>% st_union() %>% 
  as("Spatial")

non_convex_bdry[[1]] <- bdry_polygon



max.edge = diff(range(coords[,1]))/(50)


mesh <- fm_mesh_2d_inla(
  # loc = coords,
  boundary = non_convex_bdry, max.edge = c(max.edge, max.edge*2), # km inside and outside
  cutoff = 20,
  crs = fm_crs(data)
  # crs=CRS(proj4string(bdry_polygon))
) # cutoff is min edge
# plot it



# 
# # AGHY_non_convex_bdry <- inla.nonconvex.hull(AGHY_coords, -0.03, -0.05, resolution = c(100, 100))
# # AGPE_non_convex_bdry <- inla.nonconvex.hull(AGPE_coords, -0.04, -0.05, resolution = c(100, 100))
# # ELVI_non_convex_bdry <- inla.nonconvex.hull(ELVI_coords, -0.03, -0.05, resolution = c(100, 100))
# 
# sf::sf_use_s2(FALSE)
# bdry_st <- st_multipolygon(non_convex_bdry[[1]])
#                            mutate(easting = V1,  northing = V2) %>%
#                            st_as_sf(coords = c("easting", "northing"), crs = fm_crs(data)) %>%
#                            summarise(geometry = st_combine(geometry)) %>%
#                            st_cast("POLYGON"))
# # 
# # AGHY_bdry_st <- st_make_valid(as_tibble(AGHY_non_convex_bdry$loc)  %>% 
# #                                 mutate(lon = V1,  lat = V2) %>% 
# #                                 st_as_sf(coords = c("lon", "lat"), crs = 4326) %>% 
# #                                 summarise(geometry = st_combine(geometry)) %>% 
# #                                 st_cast("POLYGON"))
# # AGPE_bdry_st <- st_make_valid(as_tibble(AGPE_non_convex_bdry$loc)  %>% 
# #                                 mutate(lon = V1,  lat = V2) %>% 
# #                                 st_as_sf(coords = c("lon", "lat"), crs = 4326) %>% 
# #                                 summarise(geometry = st_combine(geometry)) %>% 
# #                                 st_cast("POLYGON"))
# # ELVI_bdry_st <- st_make_valid(as_tibble(ELVI_non_convex_bdry$loc)  %>% 
# #                                 mutate(lon = V1,  lat = V2) %>% 
# #                                 st_as_sf(coords = c("lon", "lat"), crs = 4326) %>% 
# #                                 summarise(geometry = st_combine(geometry)) %>% 
# #                                 st_cast("POLYGON"))
# 
# coastline <- st_make_valid(sf::st_as_sf(maps::map("world", regions = c("usa", "canada", "mexico"), plot = FALSE, fill = TRUE), crs = fm_crs(data)))
# # plot(coastline)
# 
# sf::st_transform(coastline, fm_crs(data))
# 
# 
# bdry <- st_intersection(coastline$geom, bdry_st)
# # AGHY_bdry <- st_intersection(coastline$geom, AGHY_bdry_st)
# # AGPE_bdry <- st_intersection(coastline$geom, AGPE_bdry_st)
# # ELVI_bdry <- st_intersection(coastline$geom, ELVI_bdry_st)
# 
# bdry_polygon <- st_cast(st_sf(bdry), "MULTIPOLYGON", group_or_split = TRUE) %>% st_union() %>% 
#   as("Spatial")
# # AGHY_bdry_polygon <- st_cast(st_sf(AGHY_bdry), "MULTIPOLYGON", group_or_split = TRUE) %>% st_union() %>% 
# #   as("Spatial")
# # AGPE_bdry_polygon <- st_cast(st_sf(AGPE_bdry), "MULTIPOLYGON", group_or_split = TRUE) %>% st_union() %>% 
# #   as("Spatial")
# # ELVI_bdry_polygon <- st_cast(st_sf(ELVI_bdry), "MULTIPOLYGON", group_or_split = TRUE) %>% st_union() %>% 
# #   as("Spatial")
# # plot(bdry_polygon)
# 
# max.edge = diff(range(coords[,2]))/(10)
# mesh10 <- inla.mesh.2d(loc = coords, max.edge = c(1,2)*max.edge,
#                        boundary = bdry_polygon,
#                        offset = c(1,4),
#                        cutoff = max.edge/(5))
# plot(mesh10)
# mesh5 <- inla.mesh.2d(loc = coords, max.edge = c(1,2)*(max.edge/2),
#                       boundary = bdry_polygon,
#                       offset = c(1,4),
#                       cutoff = max.edge/(5))
# plot(mesh5)
# # # 
# # # mesh2.5 <- inla.mesh.2d(loc = coords, max.edge = c(1,2)*(max.edge/2/2),
# # #                         boundary = bdry_polygon,
# # #                         offset = c(1,4),
# # #                         cutoff = max.edge/(10))
# # # # plot(mesh2.5)
# # 
# mesh1 <- inla.mesh.2d(loc = coords, max.edge = c(.5,2)*(max.edge/4),
#                       boundary = bdry_polygon,
#                       offset = c(1,4),
#                       cutoff = max.edge/(10))
# plot(mesh1)
# # 
# mesh <- inla.mesh.2d(loc = coords, max.edge = c(.5,2)*(max.edge/4),
#                      boundary = bdry_polygon,
#                      offset = c(1,4),
#                      cutoff = max.edge/(15),
#                      crs = 4326)

# ggplot() +
#   gg(data = mesh) +
#   geom_point(data = data, aes(x = easting, y = northing, col = species), size = 1) +
#   coord_sf()+
#   theme_bw() +
#   labs(x = "", y = "")


# make spde
# The priors from online tutorials are :   # P(practic.range < 0.05) = 0.01 # P(sigma > 1) = 0.01
# For ESA presentation, I used the following which at least "converged" but seem sensitive to choices
# for AGHY =  P(practic.range < 0.1) = 0.01 # P(sigma > 1) = 0.01
# for AGPE =  P(practic.range < 1) = 0.01 # P(sigma > 1) = 0.01
# for ELVI = P(practic.range < 1) = 0.01 # P(sigma > 1) = 0.01



# In general, the prior on the range of the spde should be bigger than the max edge of the mesh
prior_range <- max.edge*3
# the prior for the SPDE standard deviation is a bit trickier to explain, but since our data is binomial, I'm setting it to .5
prior_sigma <- 1

spde <- inla.spde2.pcmatern(
  mesh = mesh,
  prior.range = c(prior_range, 0.05),

  prior.sigma = c(prior_sigma, 0.05),
)

# spde_const <- inla.spde2.pcmatern(
#   mesh = mesh,
#   prior.range = c(prior_range, 0.05),
#   
#   prior.sigma = c(prior_sigma, 0.05),
#   # constr = TRUE
# )
# inlabru makes making weights spatial effects simpler because we don't have to make projector matrices for each effect. i.e we don't have to make an A-matrix for each spatially varying effect.
# this means we can go strat to making the components of the model


# This is the model formula with a spatial effect (spatially varying intercept) as well a spatially varying time slopes


# setting the random effects prior
pc_prec <- list(prior = "pcprec", param = c(1, 0.1))






# version of the model with overall year slope
  
year_components <- ~ Intercept(1)+ yearEffect(main = std_year, model = "linear") +
  space.int(geometry, model = spde) + # can add group argument to make specific to each species
#  time.slope(geometry, weights = std_year) + # can add constr = TRUE to make constrain to zero in the SPDE definition, which allows us to compare relative to mean condition
  collector(collector_factor, model = "iid", constr = TRUE, hyper = list(pc_prec))+
 scorer(scorer_factor, model = "iid", constr = TRUE, hyper = list(pc_prec))

# version of the model with spatially varying slopes
svc_components <- ~ 0+
  space.int(geometry, model = spde) + # can add group argument to make specific to each species
  time.slope(geometry, weights = std_year, model = spde) + # can add constr = TRUE to make constrain to zero in the SPDE definition, which allows us to compare relative to mean condition
  collector(collector_factor, model = "iid", constr = TRUE, hyper = list(pc_prec))+
  scorer(scorer_factor, model = "iid", constr = TRUE, hyper = list(pc_prec))




# formula, with "." meaning "add all the model components":
year_formula <- Endo_status_liberal ~ .
svc_formula <- Endo_status_liberal ~ .



# Now run the models


year.fit <- bru(year_components,
           like(
             formula = year_formula,
             family = "binomial",
             Ntrials = 1,
             data = endo_herb
           ),
           options = list(
             control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE),
             control.inla = list(int.strategy = "eb"),
             verbose = TRUE
           )
)

svc.fit <- bru(svc_components,
                     like(
                       formula = svc_formula,
                       family = "binomial",
                       Ntrials = 1,
                       data = endo_herb
                     ),
                     options = list(
                       control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE),
                       control.inla = list(int.strategy = "eb"),
                       verbose = TRUE
                     )
)

fit.list <- list(year_fit = year.fit, svc_fit = svc.fit)

fit_lists[[species_codes[s]]] <- list(fit.list)

mesh_list[[s]] <- mesh
bdry_polygon_list[[s]] <- bdry_polygon

}


saveRDS(fit_lists, file = "~/Dropbox/Josh&Tom - shared/Endo_Herbarium/Model_Output/fit_lists_withscorer.Rds")

saveRDS(mesh_list, file = "~/Dropbox/Josh&Tom - shared/Endo_Herbarium/Model_Output/mesh_list_withscorer.Rds")
saveRDS(bdry_polygon_list, file = "~/Dropbox/Josh&Tom - shared/Endo_Herbarium/Model_Output/boundary_polygon_list.Rds")

fit_lists_withoutscorer <- readRDS(fit_lists, file = "~/Dropbox/Josh&Tom - shared/Endo_Herbarium/Model_Output/fit_lists.Rds")


fit_lists <- readRDS(fit_lists, file = "~/Dropbox/Josh&Tom - shared/Endo_Herbarium/Model_Output/fit_lists_withscorer.Rds")
mesh_list <- readRDS( file = "~/Dropbox/Josh&Tom - shared/Endo_Herbarium/Model_Output/mesh_list_withscorer.Rds")
bdry_polygon_list <- readRDS( file = "~/Dropbox/Josh&Tom - shared/Endo_Herbarium/Model_Output/boundary_polygon_list.Rds")


fit_lists$AGHY[[1]]$year_fit$dic$dic
fit_lists$AGHY[[1]]$svc_fit$dic$dic

fit_lists$AGPE[[1]]$year_fit$dic$dic
fit_lists$AGPE[[1]]$svc_fit$dic$dic

fit_lists$ELVI[[1]]$year_fit$dic$dic
fit_lists$ELVI[[1]]$svc_fit$dic$dic


fit_lists_withoutscorer$AGHY[[1]]$year_fit$dic$dic
fit_lists_withoutscorer$AGHY[[1]]$svc_fit$dic$dic

fit_lists_withoutscorer$AGPE[[1]]$year_fit$dic$dic
fit_lists_withoutscorer$AGPE[[1]]$svc_fit$dic$dic

fit_lists_withoutscorer$ELVI[[1]]$year_fit$dic$dic
fit_lists_withoutscorer$ELVI[[1]]$svc_fit$dic$dic


fit_lists$AGHY[[1]]$year_fit$mode$mode.status
fit_lists$AGHY[[1]]$svc_fit$mode$mode.status

fit_lists$AGPE[[1]]$year_fit$mode$mode.status
fit_lists$AGPE[[1]]$svc_fit$mode$mode.status

fit_lists$ELVI[[1]]$year_fit$mode$mode.status
fit_lists$ELVI[[1]]$svc_fit$mode$mode.status

# 
# fit$summary.fixed
# fit$summary.random
# 
# fit$summary.hyperpar



################################################################################################################################
##########  Assessing model fit     ###############
################################################################################################################################

validation.pred <- list()
for(s in 1:3){
  for(m in 1:2){
data <- endo_herb_list[[s]]
# predicting the training data
validation.pred[[species_codes[s]]][[1]] <- predict(
  fit_lists[[species_codes[s]]][[1]]$year_fit,
  newdata = data,
  formula = ~ invlogit(Intercept + yearEffect + space.int + collector + scorer )) 
validation.pred[[species_codes[s]]][[2]] <- predict(
  fit_lists[[species_codes[s]]][[1]]$svc_fit,
  newdata = data,
  formula = ~ invlogit(time.slope + space.int + collector + scorer)) 
  }
}


rocobj <- pROC::roc(endo_herb_list[[1]]$Endo_status_liberal, validation.pred$AGHY[[1]]$mean)
rocobj <- pROC::roc(endo_herb_list[[1]]$Endo_status_liberal, validation.pred$AGHY[[2]]$mean)
rocobj <- pROC::roc(endo_herb_list[[2]]$Endo_status_liberal, validation.pred$AGPE[[1]]$mean)
rocobj <- pROC::roc(endo_herb_list[[2]]$Endo_status_liberal, validation.pred$AGPE[[2]]$mean)
rocobj <- pROC::roc(endo_herb_list[[3]]$Endo_status_liberal, validation.pred$ELVI[[1]]$mean)
rocobj <- pROC::roc(endo_herb_list[[3]]$Endo_status_liberal, validation.pred$ELVI[[2]]$mean)

ggroc(rocobj) 
  

# AUC values
rocobj$auc







################################################################################################################################
##########  Looking at the raw parameters from the SVC models   ###############
################################################################################################################################

# Mapping the coeffients

# get easting and northing limits


xlim <- range(mesh$loc[, 1])
ylim <- range(mesh$loc[, 2])
grd_dims <- round(c(x = diff(range(xlim))/50, y = diff(range(ylim))/50))

# make mesh projector to get model summaries from the mesh to the mapping grid
mesh_proj <- fm_evaluator(
  mesh,
  xlim = xlim, ylim = ylim, dims = grd_dims,
)



space_int <- data.frame(
  median = invlogit(fit$summary.random$space.int$"0.5quant"),
  range95 = invlogit(fit$summary.random$space.int$"0.975quant") -
    invlogit(fit$summary.random$space.int$"0.025quant")
)

time_slope <- data.frame(
  median = (exp(fit$summary.random$time.slope$"0.5quant")-1)*100,
  range95 = (exp(fit$summary.random$time.slope$"0.975quant") -
               exp(fit$summary.random$time.slope$"0.025quant"))*100
)

# # raw values
# time_slope <- data.frame(
#   median = (fit$summary.random$time_slope$"0.5quant"),
#   range95 = (fit$summary.random$time_slope$"0.975quant" -
#                fit$summary.random$time_slope$"0.025quant")
# )
# 
year_slope <- data.frame(
  median = ((fit$summary.fixed["yearEffect","0.5quant"] + fit$summary.random$time.slope$"0.5quant")),
  range95 = ((fit$summary.fixed["yearEffect","0.975quant"] + fit$summary.random$time.slope$"0.975quant") -
               (fit$summary.fixed["yearEffect","0.025quant"]+fit$summary.random$time.slope$"0.025quant"))
)

# 
# prev <- data.frame(
#   median = (fit$summary.fixed["year","0.5quant"]*2000 + fit$summary.random$space_int$"0.5quant"),
#   range95 = (fit$summary.fixed["year","0.975quant"]*2000 + fit$summary.random$time_slope$"0.975quant") -
#                (fit$summary.fixed["year","0.025quant"]*2000 + + fit$summary.random$time_slope$"0.025quant")
# )

# loop to get estimates on a mapping grid
pred_grids <- lapply(
  list(space_int = space_int, time_slope = time_slope, year_slope = year_slope),
  function(x) as.matrix(fm_evaluate(mesh_proj, x[,]))
)


# make a terra raster stack with the posterior median and range95
out_stk <- rast()
for (j in 1:length(pred_grids)) {
  mean_j <- cbind(expand.grid(x = mesh_proj$x, y = mesh_proj$y),
                  Z = c(matrix(pred_grids[[j]][, 1], grd_dims[1]))
  )
  mean_j <- rast(mean_j)
  range95_j <- cbind(expand.grid(X = mesh_proj$x, Y = mesh_proj$y),
                     Z = c(matrix(pred_grids[[j]][, 2], grd_dims[1]))
  )
  range95_j <- rast(range95_j)
  out_j <- c(mean_j, range95_j)
  terra::add(out_stk) <- out_j
}
names(out_stk) <- c(
  "space_median", "space_range95",  "time_median", "time_range95", "year_median", "year_range95"
)
#Masking the raster to our boundary


### out_stk <- raster::mask(out_stk, bdry_raster, touches = FALSE)

make_plot_field <- function(data_stk, scale_label) {
  ggplot() +
    # geom_map(data = outline_map, map = outline_map, aes(map_id = region), color = "grey", linewidth = .1, fill = "#FAF9F6")+
    # geom_map(data = states_shape, map = states_shape, aes( map_id = region), color = "grey", linewidth = .1, fill = NA)+
    geom_sf(fill = NA) +
    # coord_sf(xlim = c(-109,-68), ylim = c(21,49))+
    geom_spatraster(data = data_stk) +
    labs(x = "", y = "") +
    scale_fill_viridis_c(name = scale_label, option = "turbo", na.value = "transparent")+
    theme(text = element_text(size = 2))+
    theme_bw()
  }


pt <- make_plot_field(
  data_stk = out_stk[["time_median"]],
  scale_label = "Trend\nPosterior Median\n "
) 
   # coord_sf(xlim = c(-1000,2200), ylim = c(0,3000))
pt

# ggsave(pt, filename = "ELVI_SVC_trends_map.png", width = 10, height = 7)
  

pt_r <- make_plot_field(
  data_stk = out_stk[["time_range95"]],
  scale_label = "Trend\nPosterior 95 CI\n"
)
pt_r



py <- make_plot_field(
  data_stk = out_stk[["year_median"]],
  scale_label = "Trend\nPosterior Median\n "
) 
# coord_sf(xlim = c(-1000,2200), ylim = c(0,3000))
py

# ggsave(pt, filename = "ELVI_SVC_trends_map.png", width = 10, height = 7)


py_r <- make_plot_field(
  data_stk = out_stk[["year_range95"]],
  scale_label = "Trend\nPosterior 95 CI\n"
)
py_r




# ggsave(pt_r, filename = "AGPE_SVC_trends_CI_map.png", width = 10, height = 7)

ps <- make_plot_field(
  data_stk = out_stk[["space_median"]],
  scale_label = "Spatial\nPosterior Median"
)
ps

ps_r <- make_plot_field(
  data_stk = out_stk[["space_range95"]],
  scale_label = "Spatial\nPosterior 95 CI\n"
)
ps_r









# plotting the prevalence (without any temporal stuff)
pp <- make_plot_field(
  data_stk = out_stk[["prev_median"]],
  scale_label = "Prevalence\nPosterior Median"
)
pp

pp_r <- make_plot_field(
  data_stk = out_stk[["prev_range95"]],
  scale_label = "Spatial\nPosterior 95 CI\n"
)
pp_r

# Taking alot of material from this blog post: https://inlabru-org.github.io/inlabru/articles/2d_lgcp_covars.html

# plotting the spde range and sd posteriors
spde.range <- spde.posterior(fit, "space.int", what = "range")
spde.logvar <- spde.posterior(fit, "space.int", what = "log.variance")
spde.var <- spde.posterior(fit, "space.int", what = "variance")

spde.range <- spde.posterior(fit, "time.slope", what = "range")
spde.logvar <- spde.posterior(fit, "time.slope", what = "log.variance")
spde.var <- spde.posterior(fit, "time.slope", what = "variance")

range.plot <- plot(spde.range)
var.plot <- plot(spde.var)


range.plot
var.plot

# and plot the matern covariance (our spatial decay effect)

cov.plot <- plot(spde.posterior(fit, "space.int", what = "matern.covariance"))
cov.plot <- plot(spde.posterior(fit, "time.slope", what = "matern.covariance"))

cov.plot




################################################################################################################################
##########  Getting and plotting prediction from overall trend model  ###############
################################################################################################################################
year.pred <- list()

for (s in 1:3){
  data <- endo_herb_list[[s]]
min_year <- min(data$std_year)
max_year <- max(data$std_year)
preddata <- expand.grid(spp_code = species_codes[s], species = species_names[s], std_year = seq(min_year, max_year, length.out = 1000),
                        collector_factor = NA, scorer_factor = NA) 

# gennerating predictions and back-transforming the standardized year variable

mean_year <- mean(endo_herb$year)
sd_year <- sd(endo_herb$year)



year.pred[[s]] <- predict(
  fit_lists[[species_codes[s]]][[1]]$year_fit,
  newdata = preddata,
  formula = ~ invlogit(Intercept + yearEffect + collector +scorer)) %>% 
  mutate(year = std_year + mean_year)

}


year.pred_df <- bind_rows(year.pred, .id = "column_label")


# binning the data for plotting
endo_herb_binned <- endo_herb %>% 
  mutate(binned_year = cut(year, breaks = 10)) %>%
  group_by(Spp_code, species,binned_year) %>%   
  summarise(mean_year = mean(year),
            mean_endo = mean(Endo_status_liberal),
            mean_lon = mean(lon),
            mean_lat = mean(lat),
            sample = n(),
            se_endo = sd(Endo_status_liberal)/sqrt(sample)) %>% 
  mutate(lat_bin = case_when(mean_lat>=35 ~ paste("43"),
                             mean_lat<35 ~ paste("35")),
         lon_bin = case_when(mean_lon<=-94 ~ paste("-90"),
                             mean_lon>-94 ~ paste("-80") ))


year_trend <- ggplot(year.pred_df) +
  geom_point(data = endo_herb_binned, aes(x = mean_year, y = mean_endo, size = sample, color = species))+
  # geom_linerange(data = endo_herb_binned, aes(x = mean_year, ymin = mean_endo-se_endo, ymax = mean_endo+se_endo))+
  geom_point(data =endo_herb, aes(x = year, y = Endo_status_liberal), shape = "|")+
  geom_line(aes(year, mean)) +
  geom_ribbon(aes(year, ymin = q0.025, ymax = q0.975), alpha = 0.2) +
  geom_ribbon(aes(year, ymin = mean - 1 * sd, ymax = mean + 1 * sd), alpha = 0.2) +
  facet_wrap(~species)+
  scale_color_manual(values = species_colors)+
  labs(y = "Endophyte Prevalence", x = "Year", color = "Species", size = "Sample Size")+
  theme_light()+
  theme(strip.background = element_blank(),
        legend.text = element_text(face = "italic"))+
  lims(y = c(0,1))

year_trend


# Making a histogram of the data
year_hist <- ggplot()+
  geom_histogram(data = endo_herb_AGHY, aes(x = year), fill = species_colors[1], bins = 120, alpha = .6)+
  geom_histogram(data = endo_herb_ELVI, aes(x = year), fill = species_colors[3], bins = 120, alpha = .6)+
  geom_histogram(data = endo_herb_AGPE, aes(x = year),fill = species_colors[2],  bins = 120, alpha = .6)+
  facet_wrap(~species)+
  theme_classic()+
  theme(strip.background = element_blank(),
        strip.text = element_text(face = "italic"))+
  theme(axis.line = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank())+
  labs(y = "# of Specimens")+
  guides(fill = "none")
year_hist


year_plot <- year_hist + year_trend + 
  plot_layout(ncol = 1, heights = c(1,2)) + plot_annotation(tag_levels = "A")
  
year_plot
ggsave(year_plot, filename = "year_plot.png", width = 7, height = 5)
ggsave(year_plot, filename = "year_plot_withscorer.png", width = 7, height = 5)





################################################################################################################################
##########  Plotting the spatially varying trends ###############
################################################################################################################################
svc.pred <- list()



for (s in 1:3){
  vrt <- NA
  data <- endo_herb_list[[s]]
  mesh <- mesh_list[[s]]
  bdry_polygon <- bdry_polygon_list[[s]]

vrt <- inlabru::fm_pixels(mesh, mask = bdry_polygon, format = "sp", dims = c(40,40))# Note this is where we can mask the output according the whatever shape, such as the host distribution

vrt@data <- data.frame(std_year = rep(mean(data$std_year), length.out = length(vrt)),
                       spp_code = rep(species_codes[s], length.out = length(vrt)),
                       species = rep(species_names[s], length.out = length(vrt)))

# vrt@data <- data.frame(std_year = rep(NA, length.out = length(vrt)))


# ggplot()+
#   gg(mesh)+
#   gg(vrt, color = "red")



svc.pred[[s]] <- predict(fit_lists[[species_codes[s]]][[1]]$svc_fit, 
                      vrt, 
                      formula = ~ (exp(time.slope)-1)*100)

}

# make a base map
world_map <- st_as_sf(maps::map("world", plot = FALSE, fill = TRUE)) %>%
  st_transform(epsg6703km)
states_map <- st_as_sf(maps::map("state", plot = FALSE, fill = TRUE)) %>%
  st_transform(epsg6703km)

dim(vrt)
# dim(svc.pred$prev)
# 
# dim(svc.pred$space_pred)
dim(svc.pred[[1]])

# ggplot()+
#   gg(svc.pred$prev, aes(fill = mean))+
#   scale_fill_viridis_c(option = "plasma", na.value = "transparent", begin = 0, end = 1)
# 
# 
# ggplot()+
#   gg(svc.pred$space_pred, aes(fill = mean))+
#   scale_fill_viridis_c(option = "turbo", na.value = "transparent")
# SEtting a common color scale across
min_trend <- max(max(svc.pred[[1]]$mean),max(svc.pred[[2]]$mean), max(svc.pred[[3]]$mean))

max_trend <- min(min(svc.pred[[1]]$mean),min(svc.pred[[2]]$mean),min(svc.pred[[3]]$mean))

trendrange <- range(svc.pred[[1]]$mean, svc.pred[[2]]$mean, svc.pred[[3]]$mean)


space_x <- range(svc.pred[[1]]@coords[,1], svc.pred[[2]]@coords[,1], svc.pred[[3]]@coords[,1])
space_y <- range(svc.pred[[1]]@coords[,2], svc.pred[[2]]@coords[,2], svc.pred[[3]]@coords[,2])


svc_time_map_AGHY <- ggplot()+
  geom_sf(data = world_map, color = "grey", linewidth = .1, fill = "#FAF9F6") +
  geom_sf(data = states_map, color = "grey", linewidth = .1, fill = "#FAF9F6") +
  coord_sf(xlim = space_x, ylim = space_y)+
  gg(svc.pred[[1]], aes(fill = mean))+
  scale_fill_viridis_c(option = "turbo", na.value = "transparent", limits = trendrange)+
  labs(title = species_names[1], fill = "% change/year", y = "Latitude", x = "Longitude")+
  theme_light()+
  theme(plot.title = element_text(face = "italic"))
# svc_time_map_AGHY

svc_time_map_AGPE<- ggplot()+
  geom_sf(data = world_map, color = "grey", linewidth = .1, fill = "#FAF9F6") +
  geom_sf(data = states_map, color = "grey", linewidth = .1, fill = "#FAF9F6") +
  gg(svc.pred[[2]], aes(fill = mean))+
  coord_sf(xlim = space_x, ylim = space_y)+
  scale_fill_viridis_c(option = "turbo", na.value = "transparent", limits = trendrange)+
  labs(title = species_names[2], fill = "% change/year", y = "Latitude", x = "Longitude")+
  theme_light()+
  theme(plot.title = element_text(face = "italic"))
  

# svc_time_map_AGPE


svc_time_map_ELVI<- ggplot()+
  geom_sf(data = world_map, color = "grey", linewidth = .1, fill = "#FAF9F6") +
  geom_sf(data = states_map, color = "grey", linewidth = .1, fill = "#FAF9F6") +
  gg(svc.pred[[3]], aes(fill = mean))+
  coord_sf(xlim = space_x, ylim = space_y)+
  scale_fill_viridis_c(option = "turbo", na.value = "transparent", limits = trendrange)+
  labs(title = species_names[3], fill = "% change/year", y = "Latitude", x = "Longitude")+
  theme_light()+
  theme(plot.title = element_text(face = "italic"))
# svc_time_map_ELVI



svc_time_map <- svc_time_map_AGHY + svc_time_map_AGPE + svc_time_map_ELVI +
  plot_layout(ncol = 1, guides = 'collect') + plot_annotation(tag_levels = "A")
ggsave(svc_time_map, filename = "svc_time_map.png", width = 6, height = 12)





################################################################################################################################
##########  Plotting the spatial Intercepts ###############
################################################################################################################################
svc.pred_space <- list()



for (s in 1:3){
  vrt <- NA
  data <- endo_herb_list[[s]]
  mesh <- mesh_list[[s]]
  bdry_polygon <- bdry_polygon_list[[s]]
  
  vrt <- inlabru::fm_pixels(mesh, mask = bdry_polygon, format = "sp", dims = c(40,40))# Note this is where we can mask the output according the whatever shape, such as the host distribution
  
  vrt@data <- data.frame(std_year = rep(mean(data$std_year), length.out = length(vrt)),
                         spp_code = rep(species_codes[s], length.out = length(vrt)),
                         species = rep(species_names[s], length.out = length(vrt)))
  
  # vrt@data <- data.frame(std_year = rep(NA, length.out = length(vrt)))
  
  
  # ggplot()+
  #   gg(mesh)+
  #   gg(vrt, color = "red")
  # 
  
  
  svc.pred_space[[s]] <- predict(fit_lists[[species_codes[s]]][[1]]$svc_fit, 
                           vrt, 
                           formula = ~ invlogit(space.int))
  
}


dim(vrt)
# dim(svc.pred$prev)
# 
# dim(svc.pred$space_pred)
dim(svc.pred_space[[3]])


min_trend <- max(max(svc.pred_space[[1]]$mean),max(svc.pred_space[[2]]$mean), max(svc.pred_space[[3]]$mean))

max_trend <- min(min(svc.pred_space[[1]]$mean),min(svc.pred_space[[2]]$mean),min(svc.pred_space[[3]]$mean))

trendrange <- range(svc.pred_space[[1]]$mean, svc.pred_space[[2]]$mean, svc.pred_space[[3]]$mean)


space_x <- range(svc.pred_space[[1]]@coords[,1], svc.pred_space[[2]]@coords[,1], svc.pred_space[[3]]@coords[,1])
space_y <- range(svc.pred_space[[1]]@coords[,2], svc.pred_space[[2]]@coords[,2], svc.pred_space[[3]]@coords[,2])


svc_space_map_AGHY <- ggplot()+
  geom_sf(data = world_map, color = "grey", linewidth = .1, fill = "#FAF9F6") +
  geom_sf(data = states_map, color = "grey", linewidth = .1, fill = "#FAF9F6") +
  coord_sf(xlim = space_x, ylim = space_y)+
  gg(svc.pred_space[[1]], aes(fill = mean))+
  scale_fill_viridis_c(option = "turbo", na.value = "transparent", limits = trendrange)+
  labs(title = species_names[1], fill = "", y = "Latitude", x = "Longitude")+
  theme_light()+
  theme(plot.title = element_text(face = "italic"))
# svc_space_map_AGHY

svc_space_map_AGPE<- ggplot()+
  geom_sf(data = world_map, color = "grey", linewidth = .1, fill = "#FAF9F6") +
  geom_sf(data = states_map, color = "grey", linewidth = .1, fill = "#FAF9F6") +
  gg(svc.pred_space[[2]], aes(fill = mean))+
  coord_sf(xlim = space_x, ylim = space_y)+
  scale_fill_viridis_c(option = "turbo", na.value = "transparent", limits = trendrange)+
  labs(title = species_names[2], fill = "% change/year", y = "Latitude", x = "Longitude")+
  theme_light()+
  theme(plot.title = element_text(face = "italic"))


# svc_space_map_AGPE


svc_space_map_ELVI<- ggplot()+
  geom_sf(data = world_map, color = "grey", linewidth = .1, fill = "#FAF9F6") +
  geom_sf(data = states_map, color = "grey", linewidth = .1, fill = "#FAF9F6") +
  gg(svc.pred_space[[3]], aes(fill = mean))+
  coord_sf(xlim = space_x, ylim = space_y)+
  scale_fill_viridis_c(option = "turbo", na.value = "transparent", limits = trendrange)+
  labs(title = species_names[3], fill = "% change/year", y = "Latitude", x = "Longitude")+
  theme_light()+
  theme(plot.title = element_text(face = "italic"))
# svc_space_map_ELVI



svc_space_map <- svc_space_map_AGHY + svc_space_map_AGPE + svc_space_map_ELVI +
  plot_layout(ncol = 1, guides = 'collect') + plot_annotation(tag_levels = "A")
ggsave(svc_space_map, filename = "svc_space_map.png", width = 6, height = 12)




################################################################################################################################
##########  Making the contemp predictions ###############
################################################################################################################################

contemp.pred <- list()


for (s in 1:1){
  data <- contemp_surveys %>% 
    filter(SpeciesID == "AGHY") %>% 
    mutate(collector_factor = NA, scorer_factor = NA) %>% 
    st_transform(epsg6703km)
  
  # gennerating predictions and back-transforming the standardized year variable
  
  mean_year <- mean(endo_herb$year)
  sd_year <- sd(endo_herb$year)
  
  
  
  contemp.pred[[s]] <- predict(
    fit_lists[[species_codes[s]]][[1]]$svc_fit,
    newdata = data,
    formula = ~ invlogit(space.int+time.slope + collector + scorer)) %>% 
    mutate(year = std_year + mean_year)
  
}


contemp.pred_df <- bind_rows(contemp.pred, .id = "column_label")



contemp_lon <- ggplot(contemp.pred_df)+
  geom_point(aes(x = lon, y = endo_prev, size = sample_size), alpha = .2)+
  # geom_smooth(aes(x = lon, y = endo_prev))+#, method = "glm", method.args = list(family = "binomial"))+
  geom_point(aes(x = lon, y = mean), color = "darkgreen") +
  geom_linerange(aes(x = lon, ymax = mean+1.96*sd, ymin = mean-1.96*sd), color = "darkgreen")+
  ylim(0,1) + theme_classic()




contemp_lat <- ggplot(contemp.pred_df)+
  geom_point(aes(x = lat, y = endo_prev, size = sample_size), alpha = .2)+
  # geom_smooth(aes(x = lon, y = endo_prev))+#, method = "glm", method.args = list(family = "binomial"))+
  geom_point(aes(x = lat, y = mean), color = "darkgreen") +
  geom_linerange(aes(x = lat, ymax = mean+1.96*sd, ymin = mean-1.96*sd), color = "darkgreen")+
  ylim(0,1) + theme_classic()


contemp_obspred <- ggplot(contemp.pred_df)+
  geom_point(aes(x = mean, y = endo_prev),color = "darkgreen", alpha = .2)+
  geom_abline(intercept = 0, slope = 1)+
  # geom_smooth(aes(x = lon, y = endo_prev))+#, method = "glm", method.args = list(family = "binomial"))+
  geom_linerange(aes(y = endo_prev, xmax = mean+1.96*sd, xmin = mean-1.96*sd), color = "darkgreen")+
  lims(x = c(0,1), y = c(0,1)) + theme_classic()
contemp_obspred



contemp_test_plot <- contemp_obspred + contemp_lon + contemp_lat+
  plot_layout(nrow = 1, guides = "collect") + plot_annotation(tag_levels = "A")
contemp_test_plot
ggsave(contemp_test_plot, filename = "contemp_test_plot.png", width = 10, height = 4)
###


data(gorillas, package = "inlabru")




matern <- INLA::inla.spde2.pcmatern(gorillas$mesh,
                                    prior.sigma = c(0.1, 0.01),
                                    prior.range = c(0.01, 0.01)
)

# Define domain of the LGCP as well as the model components (spatial SPDE effect and Intercept)

cmp <- coordinates ~ mySmooth(main = coordinates, model = matern) + Intercept(1)

# Fit the model, with "eb" instead of full Bayes
fit_gorillas <- lgcp(cmp, gorillas$nests,
            samplers = gorillas$boundary,
            domain = list(coordinates = gorillas$mesh),
            options = list(control.inla = list(int.strategy = "eb"))
)




vrt_gorillas <- fm_pixels(gorillas$mesh, format = "sp")

# we obtain these vertices as a SpatialPointsDataFrame

ggplot() +
  gg(gorillas$mesh) +
  gg(vrt_gorillas, color = "red")



mySmooth <- predict(fit_gorillas, vrt_gorillas, ~mySmooth)


ggplot() +
  gg(gorillas$mesh) +
  gg(mySmooth, aes(color = mean), size = 3)




ggplot() +
  gg(gorillas$mesh, color = mySmooth$mean)




pxl <- fm_pixels(gorillas$mesh, format = "sp")
mySmooth2 <- predict(fit, pxl, ~mySmooth)

# This will give us a SpatialPixelDataFrame with the columns we are looking for

head(mySmooth2)
ggplot() +
  gg(mySmooth2)


# will be able to use "mask" argument to remove pixels outside of species distribution object
pred.df <- fm_pixels(mesh, mask = bdry_polygon, format = "sp", dims = c(10,10))
# pred.df$decade <- 1900
year <- rep(1900, n = dim(pred.df)^2)


# 
pred.df_m <- merge(as.data.frame(pred.df),year)



spdf <- SpatialPointsDataFrame(coords = pred.df, data = pred.df_m)
int1 <- predict(fit, 
                newdata = spdf, 
                formula = ~ list(space = space_int,
                                 space_time = Intercept + space_int + year))


ggplot(data = int1$space) +gg(aes(fill = mean))
  # geom_point(data =endo_herb, aes(x = year, y = Endo_status_liberal))+
  geom_spatraster(data = int1$space, aes(fill = mean))
  geom_line(aes(year, mean)) +
  geom_ribbon(aes(year, ymin = q0.025, ymax = q0.975), alpha = 0.2) +
  geom_ribbon(aes(year, ymin = mean - 1 * sd, ymax = mean + 1 * sd), alpha = 0.2) +
  lims(y = c(0,1))








###### Trying to get prediction from the year effect #####

decade.pred <- predict(
  fit,
  data.frame(decade = seq(1800, 2020, by = 10)),
  formula = ~ decade_eval(decade),
  include = character(0) # Not needed from version 2.8.0
)



decade.pred <- decade.pred$median


ggplot(decade.pred) +
  geom_line(aes(decade, mean))+
  geom_ribbon(aes(decade,ymin = q0.025,ymax = q0.975),alpha = 0.2) +
  geom_ribbon(aes(decade,ymin = mean - 1 * sd,ymax = mean + 1 * sd),alpha = 0.2)






pred.y_data <- pred.y_data[1:9] %>% 
  mutate(decade = year)

vrt <- fm_vertices(mesh, format = "sp")
pred <- predict(fit, vrt, ~ space_int)


class(pred)
head(vrt)
head(pred)


ggplot() +
  gg(mesh, color = mean)

##### Post-hoc correlations with climate drivers #####

# Read in the saved rasters which have the calculated change in climate.

# reading in the change in climate normals (1895-1925 and 1990-2020)
prism_diff_pred_df <- read_csv(file = "prism_diff_pred_df.csv") %>% 
  distinct() %>% 
  st_as_sf(coords = c("lon", "lat"), crs = 4326, remove = FALSE) %>% 
  st_transform(epsg6703km) %>% 
  mutate(
    coords.x1 = st_coordinates(.)[, 1],
    coords.x2 = st_coordinates(.)[, 2]
  ) 

# mutate(across(contains("ppt_"), ~.x*.1)) # changing the units of the change in precip. to be 10ths of millimeters
prism_for_plotting <- prism_diff_pred_df %>% 
  filter(tmean_annual_cv_diff <25 & tmean_annual_cv_diff>-2)  


ggplot()+
  gg(prism_diff_pred_df, aes(color = tmean_annual_diff ))

ggplot()+
  gg(svc.pred[[1]], aes(color = mean))
# Merging the climate date with our predicted temporal slopes

points <- prism_diff_pred_df %>%
  dplyr::select(lon, lat) 
  

# getting the pixel values of slopes as a dataframe


dim(svc.pred[[1]])
trend_df <- bind_rows(as_tibble(svc.pred[[1]]), as_tibble(svc.pred[[3]]), as_tibble(svc.pred[[2]]), .id = "column_label")

# 
# ggplot(data=trend_df)+
#   geom_tile(aes(x = coords.x1, y = coords.x2, fill = mean))


# pred_change_merge <- pred_ %>% 
pred_change_merge <- trend_df %>% 
  left_join(prism_diff_pred_df) %>% 
  filter(tmean_annual_cv_diff <25 & tmean_annual_cv_diff>-2) 
  distinct()  
  na.omit() %>% 
  mutate(lat = as.numeric(lat), lon = as.numeric(lon)) 

pred_change_merge_subsample <- pred_change_merge %>% 
  # group_by(species) %>% 
  sample_n(size =200)

ggplot(data=pred_change_merge_subsample)+
  geom_tile(aes(x = lon, y = lat, fill = time_median))

pred_change_merge_long <- pred_change_merge_subsample %>% 
  pivot_longer(cols = contains("_diff")) %>% 
  filter(!grepl("sd", name)) %>% 
  filter(!grepl("tmean_spring", name)) %>% 
  filter(!grepl("tmean_summer", name)) %>% 
  filter(!grepl("tmean_autumn", name)) 


# Looking at the data
ppt_regression_plot <- ggplot(filter(pred_change_merge_long, grepl("ppt", name)&!grepl("cv", name) & !grepl("annual", name)))+
  geom_point(aes(x = value, y = time_median), alpha = .2)+
  geom_smooth(aes(x = value, y = time_median), color = "black", method = "lm")+
  scale_color_manual(values = species_colors) +
  facet_wrap(~factor(name, levels = c("ppt_spring_diff", "ppt_summer_diff", "ppt_autumn_diff")), labeller = as_labeller(c("ppt_spring_diff" = "Change in Spring Ppt.", "ppt_summer_diff" = "Change in Summer Ppt.", "ppt_autumn_diff" = "Change in Autumn Ppt.")), scales = "free", nrow = 3) +
  theme(text = element_text(size = 8))+
  theme_light()+
  labs(y = "Trend in endophyte prevalence", x = expression("Change in Precip. (mm"^-1*")"))
ppt_regression_plot


ggsave(ppt_regression_plot, filename = "ppt_regression_plot_ESA.png", height = 6, width = 4)


ppt_CV_regression_plot <- ggplot(filter(pred_change_merge_long, grepl("ppt", name)&grepl("cv", name) & !grepl("annual", name)))+
  geom_point(aes(x = value, y = time_median), alpha = .2)+
  geom_smooth(aes(x = value, y = time_median), color = "black", method = "lm")+
  scale_color_manual(values = species_colors) +
  facet_wrap(~factor(name, levels = c("ppt_spring_cv_diff", "ppt_summer_cv_diff", "ppt_autumn_cv_diff")), labeller = as_labeller(c("ppt_spring_cv_diff" = "Change in Spring Ppt. CV", "ppt_summer_cv_diff" = "Change in Summer Ppt. CV", "ppt_autumn_cv_diff" = "Change in Autumn Ppt. CV ")), scales = "free", nrow = 3) +
  theme(strip.text = element_text(size = rel(400)))+
  theme_light()+
  labs(y = "Trend in endophyte prevalence", x = expression("Change in Precip. (mm"^-1*")"))
ppt_CV_regression_plot


ggsave(ppt_CV_regression_plot, filename = "ppt_CV_regression_plot_ESA.png", height = 6, width = 4)


tmean_regression_plot <- ggplot(filter(pred_change_merge_long, grepl("tmean", name)))+
  geom_point(aes(x = value, y = time_median), alpha = .2)+
  geom_smooth(aes(x = value, y = time_median), color = "black", method = "lm")+
  scale_color_manual(values = species_colors) +
  facet_wrap(~factor(name, levels = c("tmean_annual_diff", "tmean_annual_cv_diff")), labeller = as_labeller(c("tmean_annual_diff" = "Change in Annual Temp.", "tmean_annual_cv_diff" = "Change in Temp. CV ")),scales = "free", nrow = 3) +
  theme(text = element_text(size = 4))+
  theme_light()+
  labs(y = "Trend in endophyte prevalence", x = expression("Change in mean temp. ("*degree*"C)"))
tmean_regression_plot

ggsave(tmean_regression_plot, filename = "tmean_regression_plot_ESA.png", height = 5, width = 4)
