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



# Creating herbariumd levels
herbarium_levels <- levels(as.factor(endo_herb_georef$Herb_code))
herbarium_no <- paste0("Herbarium",1:nlevels(as.factor(endo_herb_georef$Herb_code)))

endo_herb_georef$herbarium_factor <- herbarium_no[match(as.factor(endo_herb_georef$Herb_code), herbarium_levels)]



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

# write_csv(endo_herb, file = "filtered_endo_herb_georef.csv")

# mini_dataset <- endo_herb %>%
#   filter(year>1970 &year<1990) %>%
#   mutate(presence = Endo_status_liberal) %>%
#   filter(!is.na(Endo_status_liberal)) %>%
#   filter(!is.na(easting) & !is.na(std_year)) %>%
#   filter(lon>-110) %>%
#   filter(Country != "Canada" ) %>%
#   dplyr::select(Sample_id, species_index, species, lat, lon, easting, northing, year, std_year, presence, scorer_factor)
# 
# write_csv(mini_dataset, "2024-05-16_endophyte_data.csv")


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
  geom_point(data = filter(endo_herb, scorer_factor == "Scorer7"), aes(x = lon, y = lat, color = scorer_factor, shape = species), alpha = .7, size = 1.2)+
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

endo_herb <- endo_herb %>% 
  filter(!is.na(Endo_status_liberal)) %>% 
  filter(!is.na(lon) & !is.na(year))# dropping one plant in mexico
  
  
##########################################################################################
############ Setting up and running INLA model with inlabru ############################### 
##########################################################################################
species_codes <- c("AGHY", "AGPE", "ELVI")
species_names <- c("A. hyemalis", "A. perennans", "E. virginicus")

##### Building a spatial mesh #####

data <- endo_herb

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



max.edge = diff(range(coords[,1]))/(30)

mesh1 <- fm_mesh_2d_inla(
  # loc = coords,
  boundary = non_convex_bdry, max.edge = c(max.edge, max.edge*2), # km inside and outside
  cutoff = max.edge/5,
  crs = fm_crs(data)
  # crs=CRS(proj4string(bdry_polygon))
) # cutoff is min edge


mesh2 <- fm_mesh_2d_inla(
  # loc = coords,
  boundary = non_convex_bdry, max.edge = c(max.edge/2, (max.edge/2)*2), # km inside and outside
  cutoff = max.edge/5,
  crs = fm_crs(data)
  # crs=CRS(proj4string(bdry_polygon))
) # cutoff is min edge



mesh3 <- fm_mesh_2d_inla(
  # loc = coords,
  boundary = non_convex_bdry, max.edge = c(max.edge/4, (max.edge/2)*2), # km inside and outside
  cutoff = max.edge/5,
  crs = fm_crs(data)
  # crs=CRS(proj4string(bdry_polygon))
) # cutoff is min edge

mesh_list <- list(mesh1, mesh2)


# plot it
# plot(mesh)

# ggplot() +
#   gg(data = mesh1) +
#   geom_point(data = data, aes(x = easting, y = northing, col = species), size = 1) +
#   coord_sf()+
#   theme_bw() +
#   labs(x = "", y = "")
# 
# 
# ggplot() +
#   gg(data = mesh2) +
#   geom_point(data = data, aes(x = easting, y = northing, col = species), size = 1) +
#   coord_sf()+
#   theme_bw() +
#   labs(x = "", y = "")
# 
# ggplot() +
#   gg(data = mesh3) +
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
fit_list <- list()
time_list <- list()

# for(m in 1:length(mesh_list)){
  for(m in 2){
  mesh <- mesh_list[[m]]


# In general, the prior on the range of the spde should be bigger than the max edge of the mesh
prior_range <- max.edge*3
# the prior for the SPDE standard deviation is a bit trickier to explain, but since our data is binomial, I'm setting it to .5
prior_sigma <- 1

spde1 <- inla.spde2.pcmatern(
  mesh = mesh,
  prior.range = c(prior_range, 0.5),

  prior.sigma = c(prior_sigma, 0.5),
)

spde2 <- inla.spde2.pcmatern(
  mesh = mesh,
  prior.range = c(prior_range/5, 0.5),
  
  prior.sigma = c(prior_sigma, 0.5),
)


spde3 <- inla.spde2.pcmatern(
  mesh = mesh,
  prior.range = c(prior_range*5, 0.5),
  
  prior.sigma = c(prior_sigma, 0.5),
)

spde_list <- list(spde1, spde2, spde3)
}
# inlabru makes making weights spatial effects simpler because we don't have to make projector matrices for each effect. i.e we don't have to make an A-matrix for each spatially varying effect.
# this means we can go strat to making the components of the model


# This is the model formula with a spatial effect (spatially varying intercept) as well a spatially varying time slopes


# setting the random effects prior
# set so that the probability of random effect sd exceeding 1 is 0.01
pc_prec1 <- list(prior = "pcprec", param = c(1, 0.1))
pc_prec2 <- list(prior = "pcprec", param = c(1, 0.01))



pc_prec_list <- list(pc_prec1, pc_prec2)



for(s in 1:length(spde_list)){
  spde <- spde_list[[s]]
  for(p in 1:length(pc_prec_list)){
    pc_prec <- pc_prec_list[[p]]

# version of the model with multiple likelihoods

cmp <- ~ space(geometry, model = spde) + space.species1(geometry, model = spde) + space.species2(geometry, model = spde) + + space.species3(geometry, model = spde) +
  time.species1(geometry, weights = std_year, model = spde) + time.species2(geometry, weights = std_year, model = spde) + time.species3(geometry, weights = std_year, model = spde) +
  + int.species1(1) + int.species2(1) + int.species3(1)+
  + year.species1(main = ~0 + std_year, model = "fixed") + year.species2(main = ~0 + std_year, model = "fixed") + year.species3(main = ~0 + std_year, model = "fixed")+
  herbarium(herbarium_factor, model = "iid", constr = TRUE, hyper = list(pc_prec))+
  scorer(scorer_factor, model = "iid", constr = TRUE, hyper = list(pc_prec)) +
  collector(collector_factor, model = "iid", constr = TRUE, hyper = list(pc_prec))

fml.aghy <- Endo_status_liberal ~ 0 + int.species1 + year.species1 + space + space.species1 + time.species1 + herbarium + scorer + collector
fml.agpe <- Endo_status_liberal ~ 0 + int.species2 + year.species2 + space + space.species2 + time.species2 + herbarium + scorer + collector
fml.elvi <- Endo_status_liberal ~ 0 + int.species3 + year.species3 + space + space.species3 + time.species3 + herbarium + scorer + collector

scorer_effect <- scorer.space(geometry, model = spde) + scorer(scorer_factor, model = "iid", constr = TRUE, hyper = list(pc_prec))

lik_aghy <- like(formula = fml.aghy,
                 family = "binomial",
                 Ntrials = 1,
                 data = data[data$species_index == 1,])

lik_agpe <- like(formula = fml.agpe,
                 family = "binomial",
                 Ntrials = 1,
                 data = data[data$species_index == 2,])

lik_elvi <- like(formula = fml.elvi,
                 family = "binomial",
                 Ntrials = 1,
                 data = data[data$species_index == 3,])
start=Sys.time()

          fit <- bru(cmp,
                     lik_aghy,
                     lik_agpe,
                     lik_elvi,
                     options = list(
                       control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE),
                       control.inla = list(int.strategy = "eb"),
                       verbose = FALSE,
                       bru_max_iter =3)
          )
          fit_list[[paste("m",m,"s",s,"p",p, sep = "")]] <- fit
          time_list[[paste("m",m,"s",s,"p",p, sep = "")]] <- Sys.time()-start
      }
    }
  }

saveRDS(fit_list, file = "~/Dropbox/Josh&Tom - shared/Endo_Herbarium/Model_Output/mesh_and_prior_sensitivity.Rds")
saveRDS(time_list, file = "~/Dropbox/Josh&Tom - shared/Endo_Herbarium/Model_Output/mesh_and_prior_timing.Rds")


fit_list <- readRDS(file = "~/Dropbox/Josh&Tom - shared/Endo_Herbarium/Model_Output/mesh_and_prior_sensitivity.Rds")
time_list <- readRDS(file = "~/Dropbox/Josh&Tom - shared/Endo_Herbarium/Model_Output/mesh_and_prior_timing.Rds")

###############################################################################
####### we can look at some model comparison and convergence diagnostics ######
###############################################################################

n <- length(fit_list)
for(i in 1:n){
print(fit_list[[i]]$dic$dic)
}

for(i in 1:n){
  print(fit_list[[i]]$mode$mode.status)
}
# This makes a plot of some of the model convergence info, but I only did a few iterations of the models, so it's doesn't look quite right, probably.
# pdf("convergence.pdf", width = 20, height = 20)
# inlabru:::make_track_plots(fit_list[[1]])[["default"]]
# dev.off()



###############################################################################
####### Examining the posteriors of the scorer effect across the different prior ######
###############################################################################

# only looking at if there is a difference for the p1 vs p2 for the first mesh and last mesh
random_scorers <- sample((1:NROW(fit_list[[1]]$summary.random$scorer)), 4)

slist <- vector("list", NROW(random_scorers))
for (i in seq_along(slist)){ slist[[i]] <- plot(fit_list[[1]], "scorer", index = i) + lims(x = c(-2.5,2.5))}
multiplot(plotlist = slist, cols = 3)

for (i in seq_along(slist)){ slist[[i]] <- plot(fit_list[[2]], "scorer", index = i) + lims(x = c(-2.5,2.5))}
  multiplot(plotlist = slist, cols = 3)

for (i in seq_along(slist)){ slist[[i]] <- plot(fit_list[[11]], "scorer", index = i) + lims(x = c(-2.5,2.5))}
  multiplot(plotlist = slist, cols = 3)
  
for (i in seq_along(slist)){ slist[[i]] <- plot(fit_list[[12]], "scorer", index = i) + lims(x = c(-2.5,2.5))}
  multiplot(plotlist = slist, cols = 3)
# These plots look basically identical


###############################################################################
####### Examining the spatial effects with different priors ######
###############################################################################
spde_plot_list <- list()
for(i in c(1,3,5,7,9,11)){
  # plotting the spde range and sd posteriors
  spde.range.overall <- plot(spde.posterior(fit_list[[i]], "space", what = "range"))+labs(title = "overall spatial")
  spde.var.overall <- plot(spde.posterior(fit_list[[i]], "space", what = "variance"))+labs(title = "overall spatial")
  spde.cov.overall <- plot(spde.posterior(fit_list[[i]], "space", what = "matern.covariance"))+labs(title = "overall spatial")
  
  spde.range.species1 <- plot(spde.posterior(fit_list[[i]], "space.species1", what = "range"))+labs(title = "species spatial")
  spde.var.species1 <- plot(spde.posterior(fit_list[[i]], "space.species1", what = "variance"))+labs(title = "species spatial")
  spde.cov.species1 <- plot(spde.posterior(fit_list[[i]], "space.species1", what = "matern.covariance"))+labs(title = "species spatial")
  
  spde.range.slope1 <- plot(spde.posterior(fit_list[[i]], "time.species1", what = "range"))+labs(title = "slope spatial")
  spde.var.slope1 <- plot(spde.posterior(fit_list[[i]], "time.species1", what = "variance"))+labs(title = "slope spatial")
  spde.cov.slope1 <- plot(spde.posterior(fit_list[[i]], "time.species1", what = "matern.covariance"))+labs(title = "slope spatial")

  spde_plot_list[[i]] <- (spde.range.overall|spde.var.overall|spde.cov.overall)/
                    (spde.range.species1|spde.var.species1|spde.cov.species1)/
                    (spde.range.slope1|spde.var.slope1|spde.cov.slope1) + plot_annotation(title = names(fit_list)[i])
}

  pdf("spde_covariances.pdf", width = 35, height = 15)
  wrap_elements(spde_plot_list[[1]])|wrap_elements(spde_plot_list[[3]])|wrap_elements(spde_plot_list[[5]])|wrap_elements(spde_plot_list[[7]])|wrap_elements(spde_plot_list[[9]])|wrap_elements(spde_plot_list[[11]])
  dev.off()
  
  
  
# The priors here do change the covariance function. I don't know that this is a huge problem, but I think it's part of why this has a difficult time.
# note that these are with the first mesh, so may be somewhat different with coarser meshs.

  
###############################################################################
####### looking at the overall trend prediction from the model ######
###############################################################################
  
  min_year <- min(data$std_year)
  max_year <- max(data$std_year)
  preddata <- expand.grid(species_index = c(1,2,3), std_year = seq(min_year, max_year, length.out = 1000),
                          collector_factor = NA, scorer_factor = NA, herbarium_factor = NA) %>% 
    mutate(species = case_when(species_index == 1 ~ species_names[1],
                               species_index == 2 ~ species_names[2],
                               species_index == 3 ~ species_names[3]))
  
  # gennerating predictions and back-transforming the standardized year variable
  
  mean_year <- mean(data$year)
  sd_year <- sd(data$year)
  
year.pred <- list()

# Using these functions to look at just the fits from one set of iid random effects priors
evens <- function(x) subset(x, x %% 2 = 0)
odds <- function(x) subset(x, x %% 2 != 0)

for(i in odds(1:n)){
  year.pred.aghy <- predict(
    fit_list[[i]],
    newdata = preddata[preddata$species_index == 1,],
    formula = ~ invlogit(int.species1 + year.species1)) %>% 
    mutate(year = std_year + mean_year)
  
  year.pred.agpe <- predict(
    fit_list[[i]],
    newdata = preddata[preddata$species_index == 2,],
    formula = ~ invlogit(int.species2 + year.species2)) %>% 
    mutate(year = std_year + mean_year)
  
  year.pred.elvi <- predict(
    fit_list[[i]],
    newdata = preddata[preddata$species_index == 3,],
    formula = ~ invlogit(int.species3 + year.species3 )) %>% 
    mutate(year = std_year + mean_year)
  
  year.pred[[i]] <- bind_rows(year.pred.aghy, year.pred.agpe, year.pred.elvi) %>% mutate(model = names(fit_list)[i])
} 

year.pred.df <- list_rbind(year.pred) %>% 
  mutate(mesh = case_when(grepl("m1", model) ~ "standard mesh",
                          grepl("m2", model) ~ "finer mesh"),
         range = case_when(grepl("s1", model) ~ "prior range = 342 km",
                           grepl("s2", model) ~ "prior range = 68 km",
                           grepl("s3", model) ~ "prior range = 1714 km"))
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
  

year_trend <- ggplot(year.pred.df) +
    geom_point(data = endo_herb_binned, aes(x = mean_year, y = mean_endo, size = sample, color = species))+
    # geom_linerange(data = endo_herb_binned, aes(x = mean_year, ymin = mean_endo-se_endo, ymax = mean_endo+se_endo))+
    geom_point(data =endo_herb, aes(x = year, y = Endo_status_liberal), shape = "|")+
    geom_line(aes(year, mean)) +
    geom_ribbon(aes(year, ymin = q0.025, ymax = q0.975), alpha = 0.2) +
    geom_ribbon(aes(year, ymin = mean - 1 * sd, ymax = mean + 1 * sd), alpha = 0.2) +
    facet_grid(range~mesh+species)+
    scale_color_manual(values = species_colors)+
    labs(y = "Endophyte Prevalence", x = "Year", color = "Species", size = "Sample Size")+
    theme_light()+
    theme(strip.text = element_text(size = rel(1.2)),
          legend.text = element_text(face = "italic"))+
    lims(y = c(0,1))

  # year_trend
  

  ggsave(year_trend, file = "Plots/prior_comparison_year_plot.png", width = 20, height = 20)

  # This looks very similar across models
  
  # Making a plot of the posteriors:
  
format <- theme(panel.grid.minor = element_line(linewidth = 0.1, linetype = 'dashed',colour = "grey"))
int_posterior <- list()
for(i in 1:n){
  int_aghy <- plot(fit_list[[i]], "int.species1")+ lims(x = c(-5,5)) + geom_vline(xintercept = fit_list[[i]]$summary.fixed$mean[1], color = species_colors[1], linewidth = 1.5) + geom_vline(xintercept = 0) + theme_classic() + format + labs(title = names(fit_list)[i], subtitle = "AGHY")
  int_agpe <- plot(fit_list[[i]], "int.species2")+ lims(x = c(-5,5)) + geom_vline(xintercept = fit_list[[i]]$summary.fixed$mean[2], color = species_colors[2], linewidth = 1.5) + geom_vline(xintercept = 0) + theme_classic() + format + labs(title = names(fit_list)[i], subtitle = "AGPE")
  int_elvi <- plot(fit_list[[i]], "int.species3")+ lims(x = c(-5,5)) + geom_vline(xintercept = fit_list[[i]]$summary.fixed$mean[3], color = species_colors[3], linewidth = 1.5) + geom_vline(xintercept = 0) + theme_classic() + format + labs(title = names(fit_list)[i], subtitle = "ELVI")
  int_posterior[[i]] <- int_aghy + int_agpe + int_elvi + plot_layout(ncol = 1)}
int_plot <- wrap_elements(int_posterior[[1]])|wrap_elements(int_posterior[[2]])|wrap_elements(int_posterior[[3]])|wrap_elements(int_posterior[[4]])|wrap_elements(int_posterior[[5]])|wrap_elements(int_posterior[[6]])|wrap_elements(int_posterior[[7]])|wrap_elements(int_posterior[[8]])|wrap_elements(int_posterior[[9]])|wrap_elements(int_posterior[[10]])|wrap_elements(int_posterior[[11]])|wrap_elements(int_posterior[[12]])
ggsave(int_plot, file = "comparison_intercept_plot.png", width = 35, height = 10)

year_posterior <- list()
for(i in 1:n){
  year_aghy <- plot(fit_list[[i]], "year.species1", index = 1) + lims(x = c(-.01,.025)) + geom_vline(xintercept = fit_list[[i]]$summary.random$year.species1$mean, color = species_colors[1], linewidth = 1.5) + geom_vline(xintercept = 0) + theme_classic() + format + labs(title = names(fit_list)[i], subtitle = "AGHY")
  year_agpe <- plot(fit_list[[i]], "year.species2", index = 1) + lims(x = c(-.01,.025)) + geom_vline(xintercept = fit_list[[i]]$summary.random$year.species2$mean, color = species_colors[2], linewidth = 1.5) + geom_vline(xintercept = 0) + theme_classic() + format + labs(title = names(fit_list)[i], subtitle = "AGPE")
  year_elvi <- plot(fit_list[[i]], "year.species3", index = 1) + lims(x = c(-.01,.025)) + geom_vline(xintercept = fit_list[[i]]$summary.random$year.species3$mean, color = species_colors[3], linewidth = 1.5) + geom_vline(xintercept = 0) + theme_classic() + format + labs(title = names(fit_list)[i], subtitle = "ELVI")
  year_posterior[[i]] <- year_aghy + year_agpe + year_elvi + plot_layout(ncol = 1)}
year_plot <- wrap_elements(year_posterior[[1]])|wrap_elements(year_posterior[[2]])|wrap_elements(year_posterior[[3]])|wrap_elements(year_posterior[[4]])|wrap_elements(year_posterior[[5]])|wrap_elements(year_posterior[[6]])|wrap_elements(year_posterior[[7]])|wrap_elements(year_posterior[[8]])|wrap_elements(year_posterior[[9]])|wrap_elements(year_posterior[[10]])|wrap_elements(year_posterior[[11]])|wrap_elements(year_posterior[[12]])

ggsave(year_plot, file = "comparison_yearslope_plot.png", width = 35, height = 10)

# The intercept posteriors look essentially the same across models. 
# The slopes are of the same direction and magnitude, but do have varying credible intervals

  
###############################################################################
####### now plotting the spatially varying slopes from each model ######
###############################################################################
# make a base map
world_map <- st_as_sf(maps::map("world", plot = FALSE, fill = TRUE)) %>%
  st_transform(epsg6703km)
states_map <- st_as_sf(maps::map("state", plot = FALSE, fill = TRUE)) %>%
  st_transform(epsg6703km)


vrt <- inlabru::fm_pixels(mesh_list[[1]], mask = bdry_polygon, format = "sp", dims = c(50,50))# Note this is where we can mask the output according the whatever shape, such as the host distribution

vrt@data <- expand.grid(std_year = rep(1, length.out = length(vrt)))

svc_combo <- list()
for(i in c(1,3,5,7,9,11)){


svc.pred1 <- predict(fit_list[[i]], 
                    vrt, 
                    formula = ~ ( exp(year.species1+ time.species1)-1)*100)
svc.pred1$species <- species_names[1]

svc.pred2 <- predict(fit_list[[i]], 
                     vrt, 
                     formula = ~ ( exp(year.species2+ time.species2)-1)*100)
svc.pred2$species <- species_names[2]

svc.pred3 <- predict(fit_list[[i]], 
                     vrt, 
                     formula = ~ ( exp(year.species3+ time.species3)-1)*100)
svc.pred3$species <- species_names[3]



min_trend <- min(svc.pred1$mean, svc.pred2$mean, svc.pred3$mean)

max_trend <- max(svc.pred1$mean, svc.pred2$mean, svc.pred3$mean)

trendrange <- range(svc.pred1$mean, svc.pred2$mean, svc.pred3$mean)

model_names <- data.frame(names(fit_list)) %>% 
  mutate(mesh = case_when(grepl("m1", names.fit_list.) ~ "standard mesh",
                          grepl("m2", names.fit_list.) ~ "finer mesh"),
         range = case_when(grepl("s1", names.fit_list.) ~ "prior range = 342 km",
                           grepl("s2", names.fit_list.) ~ "prior range = 68 km",
                           grepl("s3", names.fit_list.) ~ "prior range = 1714 km"))

space_x <- range(svc.pred1@coords[,1], svc.pred2@coords[,1], svc.pred3@coords[,1])
space_y <- range(svc.pred1@coords[,2], svc.pred2@coords[,2], svc.pred3@coords[,2])



svc_time_map_AGHY <- ggplot()+
  geom_sf(data = world_map, color = "grey", linewidth = .1, fill = "#FAF9F6") +
  geom_sf(data = states_map, color = "grey", linewidth = .1, fill = "#FAF9F6") +
  coord_sf(xlim = space_x, ylim = space_y)+
  gg(svc.pred1, aes(fill = mean))+
  scale_fill_viridis_c(option = "turbo", na.value = "transparent", limits = trendrange)+
  labs( title = paste0(model_names$mesh[i], " - ", model_names$range[i]), subtitle = species_names[1],fill = "% change/year", y = "Latitude", x = "Longitude")+
  theme_light()+
  theme(text = element_text(size = rel(4)),
        plot.title = element_text(size = rel(4)),
        plot.subtitle = element_text(face = "italic"), aspect.ratio = 1)
svc_time_map_AGHY

svc_time_map_AGPE <- ggplot()+
  geom_sf(data = world_map, color = "grey", linewidth = .1, fill = "#FAF9F6") +
  geom_sf(data = states_map, color = "grey", linewidth = .1, fill = "#FAF9F6") +
  coord_sf(xlim = space_x, ylim = space_y)+
  gg(svc.pred2, aes(fill = mean))+
  scale_fill_viridis_c(option = "turbo", na.value = "transparent", limits = trendrange)+
  labs( title = paste0(model_names$mesh[i], " - ", model_names$range[i]), subtitle = species_names[2],fill = "% change/year", y = "Latitude", x = "Longitude")+
  theme_light()+
  theme(text = element_text(size = rel(4)),
        plot.title = element_text(size = rel(4)),
        plot.subtitle = element_text(face = "italic"), aspect.ratio = 1)
svc_time_map_AGPE


svc_time_map_ELVI <- ggplot()+
  geom_sf(data = world_map, color = "grey", linewidth = .1, fill = "#FAF9F6") +
  geom_sf(data = states_map, color = "grey", linewidth = .1, fill = "#FAF9F6") +
  coord_sf(xlim = space_x, ylim = space_y)+
  gg(svc.pred3, aes(fill = mean))+
  scale_fill_viridis_c(option = "turbo", na.value = "transparent", limits = trendrange)+
  labs( title = paste0(model_names$mesh[i], " - ", model_names$range[i]), subtitle = species_names[3], fill = "% change/year", y = "Latitude", x = "Longitude")+
  theme_light()+
  theme(text = element_text(size = rel(4)),
        plot.title = element_text(size = rel(4)),
        plot.subtitle = element_text(face = "italic"), aspect.ratio = 1)
svc_time_map_ELVI


svc_combo[[i]] <- svc_time_map_AGHY|svc_time_map_AGPE|svc_time_map_ELVI + plot_layout(guides = "collect")
}

# svc_plot <- wrap_elements(svc_combo[[1]])/wrap_elements(svc_combo[[2]])/wrap_elements(svc_combo[[3]])/wrap_elements(svc_combo[[4]])/wrap_elements(svc_combo[[5]])/wrap_elements(svc_combo[[6]])/wrap_elements(svc_combo[[7]])/wrap_elements(svc_combo[[8]])/wrap_elements(svc_combo[[9]])/wrap_elements(svc_combo[[10]])/wrap_elements(svc_combo[[11]])/wrap_elements(svc_combo[[12]]) + plot_annotation(title = "SVC - year")
svc_plot_standard <- wrap_elements(svc_combo[[1]])/wrap_elements(svc_combo[[3]])/wrap_elements(svc_combo[[5]]) + plot_annotation(tag_levels = "A")

svc_plot_finer<- wrap_elements(svc_combo[[7]])/wrap_elements(svc_combo[[9]])/wrap_elements(svc_combo[[11]]) + plot_annotation(tag_levels = "A")

ggsave(svc_plot_standard, file = "Plots/standard_mesh_comparison_svc_plot.png", height = 15, width = 15, limitsize = FALSE)
ggsave(svc_plot_finer, file = "Plots/finer_mesh_comparison_svc_plot.png", height = 15, width = 15, limitsize = FALSE)


# Changing the mesh doesn't make a difference, changing the priors on scorer doesn't make a difference. Changing the spde prior does matter.
# Interestingly, the spde that is more fine has lower DIC values, but there isn't a big difference between the medium and the fine, but the coarser mesh does more poorly on DIC. 





###############################################################################
####### now plotting the species specific spatial effect from each model ######
###############################################################################
# use base map from last section
# use vertices from last section

space_combo <- list()
for(i in c(1,3,5,7,9,11)){
  
  
  space.pred1 <- predict(fit_list[[i]], 
                       vrt, 
                       formula = ~ space.species1)
  space.pred1$species <- species_names[1]
  
  space.pred2 <- predict(fit_list[[i]], 
                       vrt, 
                       formula = ~ space.species2)
  space.pred2$species <- species_names[2]
  
  space.pred3 <- predict(fit_list[[i]], 
                       vrt, 
                       formula = ~ space.species3)
  space.pred3$species <- species_names[3]
  
  
  
  min_trend <- min(space.pred1$mean, space.pred2$mean, space.pred3$mean)
  
  max_trend <- max(space.pred1$mean, space.pred2$mean, space.pred3$mean)
  
  trendrange <- range(space.pred1$mean, space.pred2$mean, space.pred3$mean)
  
  
  
  space_x <- range(space.pred1@coords[,1], space.pred2@coords[,1], space.pred3@coords[,1])
  space_y <- range(space.pred1@coords[,2], space.pred2@coords[,2], space.pred3@coords[,2])
  
  
  
  space_time_map_AGHY <- ggplot()+
    geom_sf(data = world_map, color = "grey", linewidth = .1, fill = "#FAF9F6") +
    geom_sf(data = states_map, color = "grey", linewidth = .1, fill = "#FAF9F6") +
    coord_sf(xlim = space_x, ylim = space_y)+
    gg(space.pred1, aes(fill = mean))+
    scale_fill_viridis_c(option = "turbo", na.value = "transparent", limits = trendrange)+
    labs( title = names(fit_list)[i], subtitle = species_names[3],fill = "mean", y = "Latitude", x = "Longitude")+
    theme_light()+
    theme(plot.subtitle = element_text(face = "italic"), aspect.ratio = 1)
  space_time_map_AGHY
  
  space_time_map_AGPE <- ggplot()+
    geom_sf(data = world_map, color = "grey", linewidth = .1, fill = "#FAF9F6") +
    geom_sf(data = states_map, color = "grey", linewidth = .1, fill = "#FAF9F6") +
    coord_sf(xlim = space_x, ylim = space_y)+
    gg(space.pred2, aes(fill = mean))+
    scale_fill_viridis_c(option = "turbo", na.value = "transparent", limits = trendrange)+
    labs( title = names(fit_list)[i], subtitle = species_names[3],fill = "mean", y = "Latitude", x = "Longitude")+
    theme_light()+
    theme(plot.subtitle = element_text(face = "italic"), aspect.ratio = 1)
  space_time_map_AGPE
  
  
  space_time_map_ELVI <- ggplot()+
    geom_sf(data = world_map, color = "grey", linewidth = .1, fill = "#FAF9F6") +
    geom_sf(data = states_map, color = "grey", linewidth = .1, fill = "#FAF9F6") +
    coord_sf(xlim = space_x, ylim = space_y)+
    gg(space.pred3, aes(fill = mean))+
    scale_fill_viridis_c(option = "turbo", na.value = "transparent", limits = trendrange)+
    labs( title = names(fit_list)[i], subtitle = species_names[3], fill = "mean", y = "Latitude", x = "Longitude")+
    theme_light()+
    theme(plot.subtitle = element_text(face = "italic"), aspect.ratio = 1)
  space_time_map_ELVI
  
  
  space_combo[[i]] <- space_time_map_AGHY|space_time_map_AGPE|space_time_map_ELVI + plot_layout(guides = "collect")
}

# space_plot <- wrap_elements(space_combo[[1]])/wrap_elements(space_combo[[2]])/wrap_elements(space_combo[[3]])/wrap_elements(space_combo[[4]])/wrap_elements(space_combo[[5]])/wrap_elements(space_combo[[6]])/wrap_elements(space_combo[[7]])/wrap_elements(space_combo[[8]])/wrap_elements(space_combo[[9]])/wrap_elements(space_combo[[10]])/wrap_elements(space_combo[[11]])/wrap_elements(space_combo[[12]]) + plot_annotation(title = "SVC - year")
space_plot <- wrap_elements(space_combo[[1]])/wrap_elements(space_combo[[3]])/wrap_elements(space_combo[[5]])/wrap_elements(space_combo[[7]])/wrap_elements(space_combo[[9]])/wrap_elements(space_combo[[11]]) + plot_annotation(title = "Spatial")

ggsave(space_plot, file = "comparison_space_plot.png", height = 60, width = 30, limitsize = FALSE)


# These look much more similar across species, which is good. The fine spatial effect with the first mesh (row 2) looks very different to the rest. It makes sense that the finer spatial effect maybe asking too much of the model so you get identifiability problems. But I didn't do the full number of iterations, so maybe it would have figured this out.


###############################################################################
####### plotting the shared spatial effect from each model ######
###############################################################################
# use base map from last section
# use vertices from last section

shared_combo <- list()
for(i in c(1,3,5,7,9,11)){
  
  
  shared.pred <- predict(fit_list[[i]], 
                         vrt, 
                         formula = ~ space)
  
  
  min_trend <- min(shared.pred$mean)
  
  max_trend <- max(shared.pred$mean)
  
  trendrange <- range(shared.pred$mean)
  
  
  
  space_x <- range(shared.pred@coords[,1])
  space_y <- range(shared.pred@coords[,2])
  
  
  
  shared_space_map<- ggplot()+
    geom_sf(data = world_map, color = "grey", linewidth = .1, fill = "#FAF9F6") +
    geom_sf(data = states_map, color = "grey", linewidth = .1, fill = "#FAF9F6") +
    coord_sf(xlim = space_x, ylim = space_y)+
    gg(shared.pred, aes(fill = mean))+
    scale_fill_viridis_c(option = "turbo", na.value = "transparent", limits = trendrange)+
    labs( title = names(fit_list)[i],fill = "mean", y = "Latitude", x = "Longitude")+
    theme_light()+
    theme(plot.subtitle = element_text(face = "italic"), aspect.ratio = 1)
  # shared_space_map
  

  shared_combo[[i]] <- shared_space_map + plot_layout(guides = "collect")
}

# shared_plot <- wrap_elements(shared_combo[[1]])/wrap_elements(shared_combo[[2]])/wrap_elements(shared_combo[[3]])/wrap_elements(shared_combo[[4]])/wrap_elements(shared_combo[[5]])/wrap_elements(shared_combo[[6]])/wrap_elements(shared_combo[[7]])/wrap_elements(shared_combo[[8]])/wrap_elements(shared_combo[[9]])/wrap_elements(shared_combo[[10]])/wrap_elements(shared_combo[[11]])/wrap_elements(shared_combo[[12]]) + plot_annotation(title = "SVC - year")
shared_plot <- wrap_elements(shared_combo[[1]])/wrap_elements(shared_combo[[3]])/wrap_elements(shared_combo[[5]])/wrap_elements(shared_combo[[7]])/wrap_elements(shared_combo[[9]])/wrap_elements(shared_combo[[11]]) + plot_annotation(title = "Spatial")

ggsave(shared_plot, file = "comparison_shared_plot.png", height = 40, width = 15, limitsize = FALSE)


# These look similar across meshes, but different with different spde


###############################################################################
####### Assessing model predictive ability on the training data ######
###############################################################################


data_train <- endo_herb %>% mutate(row_number = row_number())
data1 <- data_train %>% filter(species_index == 1)
data2 <- data_train %>% filter(species_index == 2)
data3 <- data_train %>% filter(species_index == 3)

# predicting the training data
validation.pred <- list()
for(i in c(1,3,5,7,9,11)){
validation.pred1 <- predict(
  fit_list[[i]],
  newdata = data1,
  formula = ~ invlogit(int.species1 + year.species1 + 
                         space +  space.species1 + 
                         time.species1 + 
                         herbarium + scorer + collector)) 

validation.pred2 <- predict(
  fit_list[[i]],
  newdata = data2,
  formula = ~ invlogit(int.species2 + year.species2 + 
                         space + space.species3 + 
                         time.species2 + 
                         herbarium + scorer + collector)) 

validation.pred3 <- predict(
  fit_list[[i]],
  newdata = data3,
  formula = ~ invlogit(int.species3 + year.species3 +
                         space + space.species3 + 
                         time.species3 +
                         herbarium + scorer + collector)) 

validation.pred[[i]] <- bind_rows(tibble(validation.pred1), tibble(validation.pred2), tibble(validation.pred3)) %>% 
  arrange(row_number) %>% 
  mutate(model = names(fit_list)[i])
}

rocobj <- list()
for(i in c(1,3,5,7,9,11)){
rocobj[[i]] <- pROC::roc(endo_herb$Endo_status_liberal, validation.pred[[i]]$mean)
}


# AUC values
for(i in c(1,3,5,7,9,11)){
print(rocobj[[i]]$auc)
}

# plotting ROC curves
multiplot(ggroc(rocobj) )

roclist <- vector("list", length(c(1,3,5,7,9,11)))
for (i in c(1,3,5,7,9,11)){ roclist[[i]] <- ggroc(rocobj[[i]])+labs(title = names(fit_list)[i]) + theme(aspect.ratio = .5)}
pdf("comparison_ROC_plots.pdf", width = 10, height = 10)
multiplot(plotlist = roclist, cols = 2)
dev.off()
# Technically one of these gives a higher predictive AUC value, but they are all basically the same.




################################################################################################################################
#######################################################################################################
####### Comparing model on the coarse mesh, but with with without different random effects ######
#######################################################################################################
################################################################################################################################



mod_list <- list()
time_list <- list()


mesh <- mesh_list[[1]]
  
  
  # In general, the prior on the range of the spde should be bigger than the max edge of the mesh
  prior_range <- max.edge*3
  # the prior for the SPDE standard deviation is a bit trickier to explain, but since our data is binomial, I'm setting it to .5
  prior_sigma <- 1
  
  spde1 <- inla.spde2.pcmatern(
    mesh = mesh,
    prior.range = c(prior_range, 0.5),
    
    prior.sigma = c(prior_sigma, 0.5),
  )
  
  # spde2 <- inla.spde2.pcmatern(
  #   mesh = mesh,
  #   prior.range = c(prior_range/5, 0.5),
  #   
  #   prior.sigma = c(prior_sigma, 0.5),
  # )
  # 
  # 
  # spde3 <- inla.spde2.pcmatern(
  #   mesh = mesh,
  #   prior.range = c(prior_range*5, 0.5),
  #   
  #   prior.sigma = c(prior_sigma, 0.5),
  # )
  # 
  spde_list <- list(spde1)
  
  # inlabru makes making weights spatial effects simpler because we don't have to make projector matrices for each effect. i.e we don't have to make an A-matrix for each spatially varying effect.
  # this means we can go strat to making the components of the model
  
  
  # This is the model formula with a spatial effect (spatially varying intercept) as well a spatially varying time slopes
  
  
  # setting the random effects prior
  # set so that the probability of random effect sd exceeding 1 is 0.01
  pc_prec1 <- list(prior = "pcprec", param = c(1, 0.1))
  # pc_prec2 <- list(prior = "pcprec", param = c(1, 0.01))
  
  
  
  pc_prec_list <- list(pc_prec1)
  

  for(s in 1:length(spde_list)){
    spde <- spde_list[[s]]
    for(p in 1:length(pc_prec_list)){
      pc_prec <- pc_prec_list[[p]]
      
      
      # version of the model with multiple likelihoods
      
      cmp <- ~ space(geometry, model = spde) + space.species1(geometry, model = spde) + space.species2(geometry, model = spde) + + space.species3(geometry, model = spde) +
        time.species1(geometry, weights = std_year, model = spde) + time.species2(geometry, weights = std_year, model = spde) + time.species3(geometry, weights = std_year, model = spde) +
        + int.species1(1) + int.species2(1) + int.species3(1)+
        + year.species1(main = ~0 + std_year, model = "fixed") + year.species2(main = ~0 + std_year, model = "fixed") + year.species3(main = ~0 + std_year, model = "fixed")+
        herbarium(herbarium_factor, model = "iid", constr = TRUE, hyper = list(pc_prec))+
        scorer(scorer_factor, model = "iid", constr = TRUE, hyper = list(pc_prec)) +
        collector(collector_factor, model = "iid", constr = TRUE, hyper = list(pc_prec))
      
      fml.aghy1 <- Endo_status_liberal ~ 0 + int.species1 + year.species1 + space + space.species1 + time.species1 + herbarium + scorer + collector
      fml.aghy2 <- Endo_status_liberal ~ 0 + int.species1 + year.species1 + space + space.species1 + time.species1 + scorer + collector
      fml.aghy3 <- Endo_status_liberal ~ 0 + int.species1 + year.species1 + space + space.species1 + time.species1 + collector

      
      fml.agpe1 <- Endo_status_liberal ~ 0 + int.species2 + year.species2 + space + space.species2 + time.species2 + herbarium + scorer + collector
      fml.agpe2 <- Endo_status_liberal ~ 0 + int.species2 + year.species2 + space + space.species2 + time.species2 + scorer + collector
      fml.agpe3 <- Endo_status_liberal ~ 0 + int.species2 + year.species2 + space + space.species2 + time.species2 + collector
      
      fml.elvi1 <- Endo_status_liberal ~ 0 + int.species3 + year.species3 + space + space.species3 + time.species3 + herbarium + scorer + collector
      fml.elvi2 <- Endo_status_liberal ~ 0 + int.species3 + year.species3 + space + space.species3 + time.species3 + scorer + collector
      fml.elvi3 <- Endo_status_liberal ~ 0 + int.species3 + year.species3 + space + space.species3 + time.species3 + collector
      
      fml_aghy_list <- list(fml.aghy1, fml.aghy1, fml.aghy3)
      fml_agpe_list <- list(fml.agpe1, fml.agpe1, fml.agpe3)
      fml_elvi_list <- list(fml.elvi1, fml.elvi1, fml.elvi3)
      for(f in 1:length(fml_aghy_list)){
      lik_aghy <- like(formula = fml_aghy_list[[f]],
                       family = "binomial",
                       Ntrials = 1,
                       data = data[data$species_index == 1,])
      
      lik_agpe <- like(formula = fml_agpe_list[[f]],
                       family = "binomial",
                       Ntrials = 1,
                       data = data[data$species_index == 2,])
      
      lik_elvi <- like(formula = fml_elvi_list[[f]],
                       family = "binomial",
                       Ntrials = 1,
                       data = data[data$species_index == 3,])
      start=Sys.time()
      
      fit <- bru(cmp,
                 lik_aghy,
                 lik_agpe,
                 lik_elvi,
                 options = list(
                   control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE),
                   control.inla = list(int.strategy = "eb"),
                   verbose = FALSE,
                   bru_max_iter =3)
      )
      mod_list[[paste("m",1,"s",s,"p",p,"f",f, sep = "")]] <- fit
      time_list[[paste("m",1,"s",s,"p",p,"f",f, sep = "")]] <- Sys.time()-start
    }
  }
}

saveRDS(mod_list, file = "~/Dropbox/Josh&Tom - shared/Endo_Herbarium/Model_Output/withwithoutscorer_sensitivity.Rds")
saveRDS(time_list, file = "~/Dropbox/Josh&Tom - shared/Endo_Herbarium/Model_Output/withwithoutscorer_timing.Rds")

mod_list <- readRDS( file = "~/Dropbox/Josh&Tom - shared/Endo_Herbarium/Model_Output/withwithoutscorer_sensitivity.Rds")
time_list <- readRDS( file = "~/Dropbox/Josh&Tom - shared/Endo_Herbarium/Model_Output/withwithoutscorer_timing.Rds")


##############################################################################
####### we can look at some model comparison and convergence diagnostics ######
###############################################################################

m <- length(mod_list)
for(i in 1:m){
  print(mod_list[[i]]$dic$dic)
}

for(i in 1:m){
  print(mod_list[[i]]$mode$mode.status)
}

###############################################################################
####### looking at the overall trend prediction from the model ######
###############################################################################

min_year <- min(data$std_year)
max_year <- max(data$std_year)
preddata <- expand.grid(species_index = c(1,2,3), std_year = seq(min_year, max_year, length.out = 1000),
                        collector_factor = NA, scorer_factor = "Scorer1000", herbarium_factor = NA) %>% 
  mutate(species = case_when(species_index == 1 ~ species_names[1],
                             species_index == 2 ~ species_names[2],
                             species_index == 3 ~ species_names[3]))

# gennerating predictions and back-transforming the standardized year variable

mean_year <- mean(data$year)
sd_year <- sd(data$year)
formula_list.aghy <- list("~ invlogit(int.species1 + year.species1 + herbarium + scorer_eval() + collector)", 
                     "~ invlogit(int.species1 + year.species1 + scorer_eval() + collector)", 
                     "~ invlogit(int.species1 + year.species1 + collector)")
formula_list.agpe <- list("~ invlogit(int.species2 + year.species2 + herbarium + scorer + collector)", 
                          "~ invlogit(int.species2 + year.species2 + scorer + collector)", 
                          "~ invlogit(int.species2 + year.species2 + collector)")
formula_list.elvi <- list("~ invlogit(int.species3 + year.species3 + herbarium + scorer + collector)", 
                          "~ invlogit(int.species3 + year.species3 + scorer + collector)", 
                          "~ invlogit(int.species3 + year.species3 + collector)")

year.pred <- list()
m <- length(mod_list)
for(i in 1:m){
  year.pred.aghy <- predict(
    mod_list[[i]],
    newdata = preddata[preddata$species_index == 1,],
    formula = formula(formula_list.aghy[[i]])) %>% 
    mutate(year = std_year + mean_year)
  
  year.pred.agpe <- predict(
    mod_list[[i]],
    newdata = preddata[preddata$species_index == 2,],
    formula = formula(formula_list.agpe[[i]])) %>% 
    mutate(year = std_year + mean_year)
  
  year.pred.elvi <- predict(
    mod_list[[i]],
    newdata = preddata[preddata$species_index == 3,],
    formula = formula(formula_list.elvi[[i]])) %>% 
    mutate(year = std_year + mean_year)
  
  year.pred[[i]] <- bind_rows(year.pred.aghy, year.pred.agpe, year.pred.elvi) %>% mutate(model = names(mod_list)[i])
} 

year.pred.df <- list_rbind(year.pred)
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


year_trend <- ggplot(year.pred.df) +
  geom_point(data = endo_herb_binned, aes(x = mean_year, y = mean_endo, size = sample, color = species))+
  # geom_linerange(data = endo_herb_binned, aes(x = mean_year, ymin = mean_endo-se_endo, ymax = mean_endo+se_endo))+
  geom_point(data =endo_herb, aes(x = year, y = Endo_status_liberal), shape = "|")+
  geom_line(aes(year, mean)) +
  geom_ribbon(aes(year, ymin = q0.025, ymax = q0.975), alpha = 0.2) +
  geom_ribbon(aes(year, ymin = mean - 1 * sd, ymax = mean + 1 * sd), alpha = 0.2) +
  facet_wrap(~model+species)+
  scale_color_manual(values = species_colors)+
  labs(y = "Endophyte Prevalence", x = "Year", color = "Species", size = "Sample Size")+
  theme_light()+
  theme(legend.text = element_text(face = "italic"))+
  lims(y = c(0,1))

# year_trend


ggsave(year_trend, file = "withandwithoutscorer_year_plot.png", width = 20, height = 20)

# This looks very similar across models

# Making a plot of the posteriors:

format <- theme(panel.grid.minor = element_line(linewidth = 0.1, linetype = 'dashed',colour = "grey"))
int_posterior <- list()
for(i in 1:n){
  int_aghy <- plot(mod_list[[i]], "int.species1")+ lims(x = c(-5,5)) + geom_vline(xintercept = fit_list[[i]]$summary.fixed$mean[1], color = species_colors[1], linewidth = 1.5) + geom_vline(xintercept = 0) + theme_classic() + format + labs(title = names(mod_list)[i], subtitle = "AGHY")
  int_agpe <- plot(mod_list[[i]], "int.species2")+ lims(x = c(-5,5)) + geom_vline(xintercept = fit_list[[i]]$summary.fixed$mean[2], color = species_colors[2], linewidth = 1.5) + geom_vline(xintercept = 0) + theme_classic() + format + labs(title = names(mod_list)[i], subtitle = "AGPE")
  int_elvi <- plot(mod_list[[i]], "int.species3")+ lims(x = c(-5,5)) + geom_vline(xintercept = fit_list[[i]]$summary.fixed$mean[3], color = species_colors[3], linewidth = 1.5) + geom_vline(xintercept = 0) + theme_classic() + format + labs(title = names(mod_list)[i], subtitle = "ELVI")
  int_posterior[[i]] <- int_aghy + int_agpe + int_elvi + plot_layout(ncol = 1)}
int_plot <- wrap_elements(int_posterior[[1]])|wrap_elements(int_posterior[[2]])|wrap_elements(int_posterior[[3]])|wrap_elements(int_posterior[[4]])|wrap_elements(int_posterior[[5]])|wrap_elements(int_posterior[[6]])|wrap_elements(int_posterior[[7]])|wrap_elements(int_posterior[[8]])|wrap_elements(int_posterior[[9]])|wrap_elements(int_posterior[[10]])|wrap_elements(int_posterior[[11]])|wrap_elements(int_posterior[[12]])
ggsave(int_plot, file = "withwithoutscorer_intercept_plot.png", width = 35, height = 10)

year_posterior <- list()
for(i in 1:n){
  year_aghy <- plot(mod_list[[i]], "year.species1", index = 1) + lims(x = c(-.01,.025)) + geom_vline(xintercept = mod_list[[i]]$summary.random$year.species1$mean, color = species_colors[1], linewidth = 1.5) + geom_vline(xintercept = 0) + theme_classic() + format + labs(title = names(mod_list)[i], subtitle = "AGHY")
  year_agpe <- plot(mod_list[[i]], "year.species2", index = 1) + lims(x = c(-.01,.025)) + geom_vline(xintercept = mod_list[[i]]$summary.random$year.species2$mean, color = species_colors[2], linewidth = 1.5) + geom_vline(xintercept = 0) + theme_classic() + format + labs(title = names(mod_list)[i], subtitle = "AGPE")
  year_elvi <- plot(mod_list[[i]], "year.species3", index = 1) + lims(x = c(-.01,.025)) + geom_vline(xintercept = mod_list[[i]]$summary.random$year.species3$mean, color = species_colors[3], linewidth = 1.5) + geom_vline(xintercept = 0) + theme_classic() + format + labs(title = names(mod_list)[i], subtitle = "ELVI")
  year_posterior[[i]] <- year_aghy + year_agpe + year_elvi + plot_layout(ncol = 1)}
year_plot <- wrap_elements(year_posterior[[1]])|wrap_elements(year_posterior[[2]])|wrap_elements(year_posterior[[3]])|wrap_elements(year_posterior[[4]])|wrap_elements(year_posterior[[5]])|wrap_elements(year_posterior[[6]])|wrap_elements(year_posterior[[7]])|wrap_elements(year_posterior[[8]])|wrap_elements(year_posterior[[9]])|wrap_elements(year_posterior[[10]])|wrap_elements(year_posterior[[11]])|wrap_elements(year_posterior[[12]])

ggsave(year_plot, file = "withwithoutscorer_yearslope_plot.png", width = 35, height = 10)



###############################################################################
####### now plotting the spatially varying slopes from each model ######
###############################################################################
# make a base map
world_map <- st_as_sf(maps::map("world", plot = FALSE, fill = TRUE)) %>%
  st_transform(epsg6703km)
states_map <- st_as_sf(maps::map("state", plot = FALSE, fill = TRUE)) %>%
  st_transform(epsg6703km)


vrt <- inlabru::fm_pixels(mesh_list[[1]], mask = bdry_polygon, format = "sp", dims = c(50,50))# Note this is where we can mask the output according the whatever shape, such as the host distribution

vrt@data <- expand.grid(std_year = rep(1, length.out = length(vrt)))

svc_combo <- list()
for(i in 1:m){
  
  
  svc.pred1 <- predict(mod_list[[i]], 
                       vrt, 
                       formula = ~ ( exp(year.species1+ time.species1)-1)*100)
  svc.pred1$species <- species_names[1]
  
  svc.pred2 <- predict(mod_list[[i]], 
                       vrt, 
                       formula = ~ ( exp(year.species2+ time.species2)-1)*100)
  svc.pred2$species <- species_names[2]
  
  svc.pred3 <- predict(mod_list[[i]], 
                       vrt, 
                       formula = ~ ( exp(year.species3+ time.species3)-1)*100)
  svc.pred3$species <- species_names[3]
  
  
  
  min_trend <- min(svc.pred1$mean, svc.pred2$mean, svc.pred3$mean)
  
  max_trend <- max(svc.pred1$mean, svc.pred2$mean, svc.pred3$mean)
  
  trendrange <- range(svc.pred1$mean, svc.pred2$mean, svc.pred3$mean)
  
  
  
  space_x <- range(svc.pred1@coords[,1], svc.pred2@coords[,1], svc.pred3@coords[,1])
  space_y <- range(svc.pred1@coords[,2], svc.pred2@coords[,2], svc.pred3@coords[,2])
  
  
  
  svc_time_map_AGHY <- ggplot()+
    geom_sf(data = world_map, color = "grey", linewidth = .1, fill = "#FAF9F6") +
    geom_sf(data = states_map, color = "grey", linewidth = .1, fill = "#FAF9F6") +
    coord_sf(xlim = space_x, ylim = space_y)+
    gg(svc.pred1, aes(fill = mean))+
    scale_fill_viridis_c(option = "turbo", na.value = "transparent", limits = trendrange)+
    labs( title = names(mod_list)[i], subtitle = species_names[1],fill = "% change/year", y = "Latitude", x = "Longitude")+
    theme_light()+
    theme(plot.subtitle = element_text(face = "italic"), aspect.ratio = 1)
  svc_time_map_AGHY
  
  svc_time_map_AGPE <- ggplot()+
    geom_sf(data = world_map, color = "grey", linewidth = .1, fill = "#FAF9F6") +
    geom_sf(data = states_map, color = "grey", linewidth = .1, fill = "#FAF9F6") +
    coord_sf(xlim = space_x, ylim = space_y)+
    gg(svc.pred2, aes(fill = mean))+
    scale_fill_viridis_c(option = "turbo", na.value = "transparent", limits = trendrange)+
    labs( title = names(mod_list)[i], subtitle = species_names[2],fill = "% change/year", y = "Latitude", x = "Longitude")+
    theme_light()+
    theme(plot.subtitle = element_text(face = "italic"), aspect.ratio = 1)
  svc_time_map_AGPE
  
  
  svc_time_map_ELVI <- ggplot()+
    geom_sf(data = world_map, color = "grey", linewidth = .1, fill = "#FAF9F6") +
    geom_sf(data = states_map, color = "grey", linewidth = .1, fill = "#FAF9F6") +
    coord_sf(xlim = space_x, ylim = space_y)+
    gg(svc.pred3, aes(fill = mean))+
    scale_fill_viridis_c(option = "turbo", na.value = "transparent", limits = trendrange)+
    labs( title = names(mod_list)[i], subtitle = species_names[3], fill = "% change/year", y = "Latitude", x = "Longitude")+
    theme_light()+
    theme(plot.subtitle = element_text(face = "italic"), aspect.ratio = 1)
  svc_time_map_ELVI
  
  
  svc_combo[[i]] <- svc_time_map_AGHY|svc_time_map_AGPE|svc_time_map_ELVI + plot_layout(guides = "collect")
}

# svc_plot <- wrap_elements(svc_combo[[1]])/wrap_elements(svc_combo[[2]])/wrap_elements(svc_combo[[3]])/wrap_elements(svc_combo[[4]])/wrap_elements(svc_combo[[5]])/wrap_elements(svc_combo[[6]])/wrap_elements(svc_combo[[7]])/wrap_elements(svc_combo[[8]])/wrap_elements(svc_combo[[9]])/wrap_elements(svc_combo[[10]])/wrap_elements(svc_combo[[11]])/wrap_elements(svc_combo[[12]]) + plot_annotation(title = "SVC - year")
svc_plot <- wrap_elements(svc_combo[[1]])/wrap_elements(svc_combo[[2]])/wrap_elements(svc_combo[[3]])+ plot_annotation(title = "SVC - year")

ggsave(svc_plot, file = "withwithoutscorer_svc_plot.png", height = 30, width = 30, limitsize = FALSE)





###############################################################################
####### Assessing model predictive ability on the training data ######
###############################################################################


data_train <- endo_herb %>% mutate(row_number = row_number())
data1 <- data_train %>% filter(species_index == 1)
data2 <- data_train %>% filter(species_index == 2)
data3 <- data_train %>% filter(species_index == 3)

formula_list <- list("~ invlogit(int.species1 + year.species1 + space +  space.species1 + time.species1 + herbarium + scorer + collector)", 
                     "~ invlogit(int.species1 + year.species1 + space +  space.species1 + time.species1 + scorer + collector)", 
                     "~ invlogit(int.species1 + year.species1 + space +  space.species1 + time.species1 + collector)")

# predicting the training data
validation.pred <- list()
for(i in 1:m){
  validation.pred1 <- predict(
    mod_list[[i]],
    newdata = data1,
    formula = formula(formula_list[[i]])) 
  
  validation.pred2 <- predict(
    mod_list[[i]],
    newdata = data2,
    formula = formula(formula_list[[i]])) 
  
  validation.pred3 <- predict(
    mod_list[[i]],
    newdata = data3,
    formula = formula(formula_list[[i]])) 
  
  validation.pred[[i]] <- bind_rows(tibble(validation.pred1), tibble(validation.pred2), tibble(validation.pred3)) %>% 
    arrange(row_number) %>% 
    mutate(model = names(mod_list)[i])
}

rocobj <- list()
for(i in c(1:m)){
  rocobj[[i]] <- pROC::roc(endo_herb$Endo_status_liberal, validation.pred[[i]]$mean)
}


# AUC values
for(i in c(1:m)){
  print(rocobj[[i]]$auc)
}

# plotting ROC curves
multiplot(ggroc(rocobj) )

pdf("withwithoutscorer_ROC_plots.pdf", width = 4, height = 4)
multiplot(ggroc(rocobj) )
dev.off()



###############################################################################
####### Plotting the spatial contribution of the iid random effects? ######
###############################################################################

scorer_pred_list <- list()
for(i in 1:2){
  
  
  scorer.pred <- predict(mod_list[[i]], 
                       data_train, 
                       formula = ~ scorer) 
    # mutate(model = names(mod_list)[i])

  scorer_pred_list[[i]] <- scorer.pred %>% 
    mutate(model = names(mod_list)[i])
}
  
  # min_trend <- min(scorer.pred[[1]]$mean, scorer.pred[[2]]$mean, scorer.pred[[3]]$mean)
  # 
  # max_trend <- max(scorer.pred[[1]]$mean, scorer.pred[[2]]$mean, scorer.pred[[3]]$mean)
  # 
  # trendrange <- range(scorer.pred[[1]]$mean, scorer.pred[[2]]$mean, scorer.pred[[3]]$mean)
  # 
  # 
  scorer.pred_df <- bind_rows(scorer_pred_list[[1]], scorer_pred_list[[2]]) %>% 
    mutate(year_bin = case_when(year>=1970 ~ "post-1970",
                                                                                                        year<1970 ~ "pre-1970"))

  scorer_contribution_map <- ggplot()+
    # geom_sf(data = world_map, color = "grey", linewidth = .1, fill = "#FAF9F6") +
    geom_sf(data = states_map, color = "grey", linewidth = .1, fill = "#FAF9F6") +
    # coord_sf(xlim = space_x, ylim = space_y)+
    gg(filter(scorer.pred_df, model == "m1s1p1f1"), aes(color = mean))+
    # scale_fill_viridis_c(option = "turbo", na.value = "transparent", limits = trendrange)+
    facet_grid(rows = vars(year_bin , species), cols = vars(scorer_factor))+
    labs(y = "Latitude", x = "Longitude")+
    theme_light()+
    theme(plot.subtitle = element_text(face = "italic"), aspect.ratio = 1)
  # scorer_contribution_map
  ggsave(scorer_contribution_map, filename = "scorer_contribution_map.png", width = 15, height = 25)

  
ggplot(filter(scorer.pred_df, model == "m1s1p1f1"))+
  geom_tile(aes(x = species, y =scorer_factor, fill = mean))



#
#
#
#
#
#
#
#
#
##
#
#
#
##
#
#
#
##
#
#
#
#
























################################################################################################################################
################################################################################################################################
################################################################################################################################
################################################################################################################################
################################################################################################################################
################################################################################################################################
################################################################################################################################
################################################################################################################################
################################################################################################################################
################################################################################################################################
################################################################################################################################
################################################################################################################################
################################################################################################################################
################################################################################################################################
################################################################################################################################



################################################################################################################################
##########  Making the contemp predictions ###############
################################################################################################################################

contemp.aghy <- contemp_surveys %>% 
    filter(SpeciesID == "AGHY") %>% 
    mutate(collector_factor = NA, scorer_factor = NA, species_index = 1, species = species_names[1]) %>% 
    st_transform(epsg6703km)

contemp.elvi <- contemp_surveys %>% 
  filter(SpeciesID == "ELVI") %>% 
  mutate(collector_factor = NA, scorer_factor = NA, species_index = 3, species = species_names[3]) %>% 
  st_transform(epsg6703km)


  # gennerating predictions and back-transforming the standardized year variable
  
  mean_year <- mean(endo_herb$year)
  sd_year <- sd(endo_herb$year)
  
  
  contemp.pred.aghy<- predict(
    fit,
    newdata = contemp.aghy,
    formula = ~ invlogit(int.species1 + year.species1 + space +  time.species1 + 
                           scorer + collector)) %>% 
    mutate(year = std_year + mean_year)

  
  contemp.pred.elvi<- predict(
    fit,
    newdata = contemp.elvi,
    formula = ~ invlogit(int.species3 + year.species3 +  time.species3 + 
                           scorer + collector)) %>% 
    mutate(year = std_year + mean_year)


contemp.pred <- bind_rows(contemp.pred.aghy, contemp.pred.elvi)




contemp_lon <- ggplot(contemp.pred)+
  geom_point(aes(x = lon, y = endo_prev, size = sample_size), alpha = .4)+
  geom_smooth(aes(x = lon, y = endo_prev, group = species), color = "black", method = "glm",  formula = "y ~ x", method.args = list(family = "binomial" ))+
  geom_linerange(aes(x = lon, ymin = `q0.025`, ymax = `q0.975`, color = species), alpha = .8)+
  geom_point(aes(x = lon, y = mean), shape = 4) + 
  scale_color_manual(values = c(species_colors[1], species_colors[3]))+
  ylim(0,1) + labs(y = "Endophyte Prevalance", x = "Longitude", color = "Species", size = "Sample Size")+
  facet_wrap(~species) + theme_classic()
# contemp_lon



contemp_lat <- ggplot(contemp.pred)+
  geom_point(aes(x = lat, y = endo_prev, size = sample_size), alpha = .4)+
  geom_smooth(aes(x = lat, y = endo_prev, group = species), color = "black", method = "glm",  formula = "y ~ x", method.args = list(family = "binomial" ))+
  geom_linerange(aes(x = lat, ymin = `q0.025`, ymax = `q0.975`, color = species))+
  geom_point(aes(x = lat, y = mean), shape = 4) +
  scale_color_manual(values = c(species_colors[1], species_colors[3]))+
  ylim(0,1) + labs(y = "Endophyte Prevalance", x = "Latitude",  color = "Species", size = "Sample Size")+
  facet_wrap(~species) + theme_classic()
# contemp_lat

contemp_obspred <- ggplot(contemp.pred)+
  geom_abline(intercept = 0, slope = 1)+
  geom_linerange(aes(y = endo_prev, xmin = `q0.025`, xmax = `q0.975`, color = species))+
  geom_point(aes(x = mean, y = endo_prev), shape = 4)+
  scale_color_manual(values = c(species_colors[1], species_colors[3]))+
  lims(x = c(0,1), y = c(0,1)) + labs(y = "Observed", x = "Predicted",  color = "Species")+
  facet_wrap(~species) + theme_classic()
# contemp_obspred



contemp_test_plot <- contemp_obspred + contemp_lon + contemp_lat+
  plot_layout(ncol = 1, guides = "collect") + plot_annotation(tag_levels = "A")
contemp_test_plot
ggsave(contemp_test_plot, filename = "contemp_test_plot.png", width = 10, height = 4)
###

# now looking at the ROC and AUC values for the contemporary dataset choosing only one plant from each population
# however we only have this information for AGHY
rocobj <- pROC::roc(contemp_random_sample$endo_status, contemp.pred.aghy$mean)

ggroc(rocobj) 


# AUC values
rocobj$auc


