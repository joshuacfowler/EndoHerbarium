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
library(tidybayes)

invlogit<-function(x){exp(x)/(1+exp(x))}
species_colors <- c("#1b9e77","#d95f02","#7570b3")
endophyte_colors <- c("#fdedd3","#f3c8a8", "#5a727b", "#4986c7", "#181914",  "#163381")


################################################################################
############ Read in the herbarium dataset ############################### 
################################################################################

endo_herb_georef <- read_csv(file = "~/Dropbox/Josh&Tom - shared/Endo_Herbarium/DigitizedHerbariumRecords/endo_herb_georef.csv") %>% 
  # filter(Country != "Canada") %>%
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
  mutate(year_bin = case_when(year<1969 ~ "pre-1969",
                              year>=1969 ~ "post-1969")) %>% 
  mutate(scorer_index = parse_number(scorer_factor),
         herbarium_index = parse_number(herbarium_factor),
         collector_index = parse_number(collector_factor)) %>% 
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

# calculating the min and max distance between points


# dist_vec <- cbind(endo_herb$easting, endo_herb$northing)
# dist_mat <- as.matrix(dist(dist_vec, method = "euclidean", diag = TRUE, upper = TRUE))
# min(dist_mat[dist_mat>0])
# max(dist_mat[dist_mat>0])
# median(dist_mat[dist_mat>0])
# 
# # calculating the width of georeferencing bounding boxes
# 
# endo_herb_boundbox_west <- endo_herb %>% 
#   select(north, south, east, west, lon, lat, easting, northing) %>% 
#   as_tibble() %>% 
#   st_as_sf(coords = c("west", "lat"), crs = 4326, remove = FALSE) %>% 
#   st_transform(epsg6703km) %>% 
#   mutate(
#     west_bb = st_coordinates(.)[, 1],
#     lat_bb = st_coordinates(.)[, 2]
#   ) %>% select(west_bb)
# endo_herb_boundbox_east <- endo_herb %>% 
#   select(north, south, east, west, lon, lat, easting, northing) %>% 
#   as_tibble() %>% 
#   st_as_sf(coords = c("east", "lat"), crs = 4326, remove = FALSE) %>% 
#   st_transform(epsg6703km) %>% 
#   mutate(
#     east_bb = st_coordinates(.)[, 1],
#     lat_bb = st_coordinates(.)[, 2]
#   ) %>% select(east_bb)
# 
# endo_herb_boundbox_width <- tibble(west_bb = endo_herb_boundbox_west$west_bb, east_bb = endo_herb_boundbox_east$east_bb) %>% 
#   mutate(bb_width = east_bb-west_bb)
# 
# 
# endo_herb_withbbwidth <- endo_herb %>% 
#   select(Sample_id,location_string, easting, northing, lat, lon) %>% 
#   cbind(bb_width = endo_herb_boundbox_width$bb_width) # I double checked a couple of these and the width I calculated seems to line up reasonable well with online measurements of the county East to West
# 
# bb_summary <- endo_herb_withbbwidth %>% 
#   summarize(mean_width = mean(bb_width),
#             max_width = max(bb_width),
#             min_width = min(bb_width),
#             median_width = median(bb_width)) %>% 
#   st_drop_geometry()

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
  coord_sf(xlim = c(-109,-68), ylim = c(19,49), crs = 4326)+
  # lims(x = c(-109,-68), y = c(21,49))+
  scale_color_manual(values = species_colors)+
  theme_light()+
  theme(legend.text = element_text(face = "italic"))+
  labs(x = "Longitude", y = "Latitude", color = "Host Species")

collections_map
ggsave(collections_map, filename = "Plots/collections_map.png", width = 7, height = 4)


# summarizing collections across states
# states_summary <- endo_herb %>% 
#   filter(!is.na(Endo_status_liberal), spp_code == "AGPE") %>% 
#   group_by(State) %>% 
#   summarize(n = n()) %>% arrange(-n)

endo_status_map <- ggplot()+
  geom_map(data = outline_map, map = outline_map, aes( map_id = region), color = "grey", linewidth = .1, fill = "#FAF9F6")+
  geom_map(data = states_shape, map = states_shape, aes(map_id = region), color = "grey", linewidth = .1, fill = NA)+
  geom_point(data = endo_herb, aes(x = lon, y = lat, color = endo_status_text), alpha = .7, size = .8)+
  facet_wrap(~factor(year_bin, levels = c("pre-1969", "post-1969"))+species)+
  coord_sf(xlim = c(-109,-68), ylim = c(19,49), crs = 4326)+
  scale_color_manual(values = c(endophyte_colors[2],endophyte_colors[6]))+
  theme_light()+
  theme(strip.background = element_blank(),
        strip.text  = element_text(face = "italic", color = "black"))+
  labs(x = "Longitude", y = "Latitude", color = "Endophyte Status")
endo_status_map
ggsave(endo_status_map, filename = "Plots/endo_status_map.png", width = 10, height = 5)



scorer_map <- ggplot()+
  geom_map(data = outline_map, map = outline_map, aes(map_id = region), color = "grey", linewidth = .1, fill = "#FAF9F6")+
  geom_map(data = states_shape, map = states_shape, aes(map_id = region), color = "grey", linewidth = .1, fill = NA)+
  geom_point(data = filter(endo_herb, scorer_factor == "Scorer7"), aes(x = lon, y = lat, color = scorer_factor, shape = species), alpha = .7, size = 1.2)+
  # coord_sf(xlim = c(-109,-68), ylim = c(21,49))+
  theme_light()+
  theme(legend.text = element_text(face = "italic"))+
  labs(x = "Longitude", y = "Latitude", color = "Host Species")

scorer_map
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


mesh <- fm_mesh_2d_inla(
  # loc = coords,
  boundary = non_convex_bdry, max.edge = c(max.edge, max.edge*2), # km inside and outside
  cutoff = max.edge/5,
  crs = fm_crs(data)
  # crs=CRS(proj4string(bdry_polygon))
) # cutoff is min edge
# plot it
# plot(mesh)


mesh_plot <- ggplot() +
  gg(data = mesh) +
  geom_point(data = data, aes(x = easting, y = northing, fill = species), color = "grey29", size  = 1, shape = 21) +
  coord_sf()+
  theme_bw() +
  scale_fill_manual(values = species_colors)+
  labs(x = "", y = "", color = "Species")
# mesh_plot
# ggsave(mesh_plot, filename = "Plots/mesh_plot.png", width = 6, height = 6)

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
  prior.range = c(prior_range, 0.5),
  
  prior.sigma = c(prior_sigma, 0.5),
)


# inlabru makes making weights spatial effects simpler because we don't have to make projector matrices for each effect. i.e we don't have to make an A-matrix for each spatially varying effect.
# this means we can go strat to making the components of the model


# This is the model formula with a spatial effect (spatially varying intercept) as well a spatially varying time slopes


# setting the random effects prior
pc_prec <- list(prior = "pcprec", param = c(1, 0.1))


# version of the model with multiple likelihoods

cmp <- ~  space(geometry, model = spde) + space.species1(geometry, model = spde) + space.species2(geometry, model = spde) + space.species3(geometry, model = spde) +
  time.species1(geometry, weights = std_year, model = spde) + time.species2(geometry, weights = std_year, model = spde) + time.species3(geometry, weights = std_year, model = spde) +
  + int.species1(1) + int.species2(1) + int.species3(1)+
  + year.species1(main = ~0 + std_year, model = "fixed") + year.species2(main = ~0 + std_year, model = "fixed") + year.species3(main = ~0 + std_year, model = "fixed")+
  # herbarium(herbarium_index, model = "iid", constr = TRUE, mapper = bru_mapper_index(max(data$herbarium_index)), hyper = list(pc_prec))+
  scorer(scorer_index, model = "iid", constr = TRUE, mapper = bru_mapper_index(max(data$scorer_index)), hyper = list(pc_prec)) +
  collector(collector_index, model = "iid", constr = TRUE, mapper = bru_mapper_index(max(data$collector_index, na.rm = T)), hyper = list(pc_prec))

fml.aghy <- Endo_status_liberal ~ 0 + int.species1 + year.species1 + space + space.species1 + time.species1 + scorer + collector 
fml.agpe <- Endo_status_liberal ~ 0 + int.species2 + year.species2 + space + space.species2 + time.species2 + scorer + collector 
fml.elvi <- Endo_status_liberal ~ 0 + int.species3 + year.species3 + space + space.species3 + time.species3 + scorer + collector



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




fit <- bru(cmp,
           lik_aghy,
           lik_agpe,
           lik_elvi,
           options = list(
             control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE),
             #control.inla = list(int.strategy = "eb"),
             verbose = TRUE,
             bru_max_iter = 10)
)


saveRDS(fit, file = "~/Dropbox/Josh&Tom - shared/Endo_Herbarium/Model_Output/multispecies_fit.Rds")
fit <- readRDS(file = "~/Dropbox/Josh&Tom - shared/Endo_Herbarium/Model_Output/multispecies_fit.Rds")

fit$dic$dic
fit$mode$mode.status

fit$summary.fixed
fit$summary.random


################################################################################################################################
##########  Assessing model fit     ###############
################################################################################################################################
# fml.aghy <- presence ~ 0 + int.species1 + year.species1 + space.species1 + time.species1 + scorer
# fml.agpe <- presence ~ 0 + int.species2 + year.species2 + space.species2 + time.species2 + scorer
# fml.elvi <- presence ~ 0 + int.species3 + year.species3 +  space.species3 + time.species3 + scorer


data <- endo_herb %>% mutate(row_number = row_number())
data1 <- data %>% filter(species_index == 1)
data2 <- data %>% filter(species_index == 2)
data3 <- data %>% filter(species_index == 3)

# predicting the training data

validation.pred1 <- predict(
  fit,
  newdata = data1,
  formula = ~ invlogit(int.species1 + year.species1 + 
                         space +  space.species1 + 
                         time.species1 + 
                         scorer_eval(scorer_index) + collector_eval(scorer_index)),
  n.samples = 100) 

validation.pred2 <- predict(
  fit,
  newdata = data2,
  formula = ~ invlogit(int.species2 + year.species2 + 
                         space + space.species3 + 
                         time.species2 + 
                         scorer_eval(scorer_index) + collector_eval(scorer_index)),
  n.samples = 100) 

validation.pred3 <- predict(
  fit,
  newdata = data3,
  formula = ~ invlogit(int.species3 + year.species3 +
                         space + space.species3 + 
                         time.species3 +
                         scorer_eval(scorer_index) + collector_eval(scorer_index)),
  n.samples = 100) 

validation.pred <- bind_rows(tibble(validation.pred1), tibble(validation.pred2), tibble(validation.pred3)) %>% 
  arrange(row_number)

rocobj <- pROC::roc(endo_herb$Endo_status_liberal, validation.pred$mean)

ROC_training_plot <- ggroc(rocobj) 
ggsave(ROC_training_plot, filename = "Plots/ROC_training_plot.png", width = 4, height = 4)

# AUC values
rocobj$auc
# 0.7892


# Taking alot of material from this blog post: https://inlabru-org.github.io/inlabru/articles/2d_lgcp_covars.html
post.pred1 <- generate(
  fit,
  newdata = data1,
  formula = ~ invlogit(int.species1 + year.species1 + 
                         space +  space.species1 + 
                         time.species1 + 
                         scorer+ collector),
  n.samples = 100) 
posterior_samples1 <- bind_cols(data1, post.pred1)


post.pred2 <- generate(
  fit,
  newdata = data2,
  formula = ~ invlogit(int.species2 + year.species2 + 
                         space + space.species3 + 
                         time.species2 + 
                         scorer+ collector),
  n.samples = 100) 
posterior_samples2 <- bind_cols(data2, post.pred2)

post.pred3 <- generate(
  fit,
  newdata = data3,
  formula = ~ invlogit(int.species3 + year.species3 +
                         space + space.species3 + 
                         time.species3 +
                         scorer+ collector),
  n.samples = 100) 

posterior_samples3 <- bind_cols(data3, post.pred3)

posterior_samples <- bind_rows(posterior_samples1, posterior_samples2, posterior_samples3) %>% 
  arrange(row_number) %>% 
  select(...55:...154) %>% st_drop_geometry %>% as.matrix()

# simulating datasets from the posterior samples
n_post_draws <- 100
y_sim <- matrix(NA,n_post_draws,length(data$Endo_status_liberal))

for(i in 1:n_post_draws){
y_sim[i,] <- rbinom(n = length(endo_herb$Endo_status_liberal), size = 1, prob = posterior_samples[,i])
}

saveRDS(y_sim, file = "y_sim.rds")
y_sim_df <- t(y_sim)
colnames(y_sim_df) <- paste("iter", 1:n_post_draws)
y_sim_df <- as_tibble(y_sim_df, .name_repair = "minimal") %>% 
  pivot_longer(cols = everything())

overlay_plot <- ggplot(data)+
  geom_line(data = y_sim_df, aes(x = value, group = name), stat="density", color = "red", alpha = .3) +
  geom_line(aes(x = Endo_status_liberal), stat="density") + 
  labs(x = "Endophyte Status", y = "Density")+
  theme_minimal()
overlay_plot
ggsave(overlay_plot, filename = "Plots/overlay_plot.png", width = 4, height = 4)

########## plotting the spde range and sd posteriors ################
spde.range <- spde.posterior(fit, "space", what = "range")
spde.var <- spde.posterior(fit, "space", what = "variance")

spde.range <- spde.posterior(fit, "time.species1", what = "range")
spde.var <- spde.posterior(fit, "time.species1", what = "variance")

range.plot <- plot(spde.range)
var.plot <- plot(spde.var)


range.plot
var.plot

# and plot the matern covariance (our spatial decay effect)

cov.plot <- plot(spde.posterior(fit, "space", what = "matern.covariance"))
cov.plot

cov.plot <- plot(spde.posterior(fit, "space.species2", what = "matern.covariance"))
cov.plot

cov.plot <- plot(spde.posterior(fit, "time.species2", what = "matern.covariance"))
cov.plot

cov.plot <- plot(spde.posterior(fit, "time.slope", what = "matern.covariance"))

cov.plot





##### Now I want to compare simulated data with observed data 


################################################################################################################################
##########  Getting and plotting prediction from overall trend model  ###############
################################################################################################################################

min_year <- min(data$std_year)
max_year <- max(data$std_year)
preddata <- expand.grid(species_index = c(1,2,3), std_year = seq(min_year, max_year, length.out = 1000),
                        collector_index = 9999, scorer_index = 9999, herbarium_index = 9999) %>% 
  mutate(species = case_when(species_index == 1 ~ species_names[1],
                             species_index == 2 ~ species_names[2],
                             species_index == 3 ~ species_names[3]))

# gennerating predictions and back-transforming the standardized year variable

mean_year <- mean(data$year)
sd_year <- sd(data$year)



year.pred.aghy <- predict(
  fit,
  newdata = preddata[preddata$species_index == 1,],
  formula = ~ invlogit(int.species1 + year.species1 + scorer_eval(scorer_index) + collector_eval(collector_index)),
  probs = c(0.025, 0.25, 0.5, 0.75, 0.975),
  n.samples = 100) %>% 
  mutate(year = std_year + mean_year)

year.pred.agpe <- predict(
  fit,
  newdata = preddata[preddata$species_index == 2,],
  formula = ~ invlogit(int.species2 + year.species2 + scorer_eval(scorer_index) + collector_eval(collector_index)),
  probs = c(0.025, 0.25, 0.5, 0.75, 0.975),
  n.samples = 100) %>% 
  mutate(year = std_year + mean_year)

year.pred.elvi <- predict(
  fit,
  newdata = preddata[preddata$species_index == 3,],
  formula = ~ invlogit(int.species3 + year.species3 + scorer_eval(scorer_index) + collector_eval(collector_index)),
  probs = c(0.025, 0.25, 0.5, 0.75, 0.975),
  n.samples = 100) %>% 
  mutate(year = std_year + mean_year)

year.pred <- bind_rows(year.pred.aghy, year.pred.agpe, year.pred.elvi)

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


year_trend <- ggplot(year.pred) +
  # geom_linerange(data = endo_herb_binned, aes(x = mean_year, ymin = mean_endo-se_endo, ymax = mean_endo+se_endo))+
  geom_point(data =endo_herb, aes(x = year, y = Endo_status_liberal), shape = "|")+
  geom_point(data = endo_herb_binned, aes(x = mean_year, y = mean_endo, size = sample, color = species))+
  geom_line(aes(year, mean)) +
  geom_ribbon(aes(year, ymin = q0.025, ymax = q0.975), alpha = 0.2) +
  geom_ribbon(aes(year, ymin = q0.25, ymax = q0.75), alpha = 0.2) +
  facet_wrap(~species)+
  scale_color_manual(values = species_colors)+
  labs(y = "Endophyte Prevalence", x = "Year", color = "Species", size = "Sample Size")+
  theme_light()+
  theme(strip.background = element_blank(),
        legend.text = element_text(face = "italic"),
        plot.margin = unit(c(0,.1,.1,.1), "line"))+
  lims(y = c(0,1))

# year_trend


# Making a histogram of the data
year_hist <- ggplot()+
  geom_histogram(data = endo_herb_AGHY, aes(x = year), fill = species_colors[1], bins = 120, alpha = .6)+
  geom_histogram(data = endo_herb_ELVI, aes(x = year), fill = species_colors[3], bins = 120, alpha = .6)+
  geom_histogram(data = endo_herb_AGPE, aes(x = year),fill = species_colors[2],  bins = 120, alpha = .6)+
  facet_wrap(~species)+
  theme_classic()+
  theme(strip.background = element_blank(),
        strip.text = element_text(face = "italic"),
        axis.line = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.ticks.length.x = unit(0,'lines'),
        axis.title.x = element_blank(),
        plot.margin = unit(c(0,1,0,1), "line"))+
  labs(y = "# of Specimens")+
  guides(fill = "none")
# year_hist





year_plot <- year_hist/year_trend +
 plot_layout(ncol = 1, heights = c(1, 2))+ plot_annotation(tag_levels = "A") &   theme(plot.tag.position = c(-0, .96))

year_plot
ggsave(year_plot, filename = "Plots/year_plot.png", width = 7, height = 5)



# Making a plot of the posteriors:

format <- theme(panel.grid.minor = element_line(linewidth = 0.1, linetype = 'dashed',colour = "grey"),
                axis.title = element_text(size = rel(2)),
                axis.text = element_text(size = rel(.5)),
                title = element_text(size = rel(.4)),
                plot.margin = unit(c(.1,.1,.1,.1), "cm"))

int_aghy <- plot(fit, "int.species1", size = 10)+ lims(x = c(-4.4,4)) + geom_vline(xintercept = fit$summary.fixed$mean[1], color = species_colors[1], linewidth = .5) + geom_vline(xintercept = 0) + theme_classic() + format+ labs(title = species_names[1], y = "Probability Density", x = "Intercept")
int_agpe <- plot(fit, "int.species2")+ lims(x = c(-4.4,4)) + geom_vline(xintercept = fit$summary.fixed$mean[2], color = species_colors[2], linewidth = .5) + geom_vline(xintercept = 0) + theme_classic() + format+ labs(title = species_names[2], y = "Probability Density", x = "Intercept")
int_elvi <- plot(fit, "int.species3")+ lims(x = c(-4.4,4)) + geom_vline(xintercept = fit$summary.fixed$mean[3], color = species_colors[3], linewidth = .5) + geom_vline(xintercept = 0) + theme_classic() + format+ labs(title = species_names[3], y = "Probability Density", x = "Intercept")
int_posterior <- int_aghy + int_agpe + int_elvi + plot_layout(ncol = 1, axis_titles = "collect")
# int_posterior



year_aghy <- plot(fit, "year.species1", index = 1) + lims(x = c(-.01,.025)) + geom_vline(xintercept = fit$summary.random$year.species1$mean, color = species_colors[1], linewidth = .5) + geom_vline(xintercept = 0) + theme_classic() + format + labs(title = " ", y = "Probability Density", x = "Slope", color  = "")
year_agpe <- plot(fit, "year.species2", index = 1) + lims(x = c(-.01,.025)) + geom_vline(xintercept = fit$summary.random$year.species2$mean, color = species_colors[2], linewidth = .5) + geom_vline(xintercept = 0) + theme_classic() + format + labs(title = " ", y = "Probability Density", x = "Slope", color  = "") 
year_elvi <- plot(fit, "year.species3", index = 1) + lims(x = c(-.01,.025)) + geom_vline(xintercept = fit$summary.random$year.species3$mean, color = species_colors[3], linewidth = .5) + geom_vline(xintercept = 0) + theme_classic() + format + labs(title = " ", y = "Probability Density", x = "Slope", color  = "") 
year_posterior <- year_aghy + year_agpe + year_elvi + plot_layout(ncol = 1, axis_titles = "collect")
# year_posterior




int_posterior_wrap <- int_posterior + plot_annotation(tag_levels = list(c("A", "", ""))) & theme(plot.tag = element_text(size = rel(2)))

year_posterior_wrap <- year_posterior + plot_annotation(tag_levels = list(c("B", "", ""))) & theme(plot.tag = element_text(size = rel(2)))

posterior_plot <- wrap_elements(int_posterior_wrap)|wrap_elements(year_posterior_wrap)
# posterior_plot
ggsave(posterior_plot, filename = "Plots/posterior_plot.png", width = 6, height = 5)


################################################################################################################################
##########  Plotting the posteriors of scorer and collector effects ###############
################################################################################################################################
# random_collectors <- sample((1:NROW(fit$summary.random$collector)), 12)

slist <- vector("list", NROW(fit$summary.random$scorer))
for (i in seq_along(slist)) slist[[i]] <- plot(fit, "scorer", index = i) + lims(x = c(-2.5,2.5))
multiplot(plotlist = slist, cols = 5)
pdf(file = "Plots/scorer_posterior_plot.pdf",  width = 8, height = 12)
multiplot(plotlist = slist, cols = 5)
dev.off()


random_collectors <- sample((1:NROW(fit$summary.random$collector)), 12)

clist <- vector("list", NROW(random_collectors))
for (i in 1:length(random_collectors)) clist[[i]] <- plot(fit, "collector", index = random_collectors[i]) + lims(x = c(-.1,.1))
multiplot(plotlist = clist, cols = 3)
pdf(file = "Plots/collector_posterior_plot.pdf",  width = 8, height = 12)
multiplot(plotlist = clist, cols = 3)
dev.off()

plot(fit, "scorer")
plot(fit, "collector")

plot(fit, "space")
plot(fit, "time.species1")



n_draws <- 500

posteriors.collector <- generate(
  fit,
  formula = ~ collector_latent,
  n.samples = n_draws) 

posteriors.scorer <- generate(
  fit,
  formula = ~ scorer_latent,
  n.samples = n_draws) 


colnames(posteriors.collector) <- c( paste0("iter",1:n_draws))
colnames(posteriors.scorer) <- c( paste0("iter",1:n_draws))

rownames(posteriors.collector) <- c(paste0("collector",1:nrow(posteriors.collector)))
rownames(posteriors.scorer) <- c(paste0("scorer",1:nrow(posteriors.scorer)))


posteriors <- rbind(posteriors.collector, posteriors.scorer)



posteriors_df <- as_tibble(posteriors, rownames = "param") %>% 
  mutate(param_type = case_when(grepl("collector", param) ~ "collector",
                                grepl("scorer", param) ~ "scorer"),
         number = parse_number(param)) %>% 
  pivot_longer( cols = -c(param, param_type, number), names_to = "iteration")

posteriors_summary <- posteriors_df %>% 
  group_by(param, number, param_type) %>% 
  dplyr::summarise(median = median(value),
                   Q97.5 = quantile(value, .975),
                   Q2.5 = quantile(value, .025))

collector_plot <- ggplot(data = filter(posteriors_summary, param_type == "collector"))+
  geom_vline(xintercept = 0)+
  geom_point(aes(x = median, y = number), alpha = .8)+
  geom_linerange(aes(xmin = Q2.5, xmax = Q97.5, y = number), alpha = .3)+
  labs(y = "Collector Number", y = "Posterior Estimate") + theme_minimal()
collector_plot
ggsave(collector_plot, filename = "Plots/collector_posterior.png", width = 5, height = 10)

scorer_plot <- ggplot(data = filter(posteriors_summary, param_type == "scorer"))+
  geom_vline(xintercept = 0)+
  geom_point(aes(x = median, y = number), alpha = .8)+
  geom_linerange(aes(xmin = Q2.5, xmax = Q97.5, y = number), alpha = .3)+
  labs(y = "Scorer Number", x = "Posterior Estimate") + theme_minimal()
scorer_plot
ggsave(scorer_plot, filename = "Plots/scorer_posterior.png", width = 5, height = 10)


################################################################################################################################
##########  Generating posterior samples for some summary statistics ###############
################################################################################################################################
n_draws <- 500

posteriors.aghy <- generate(
  fit,
  formula = ~ year.species1_latent,
  n.samples = n_draws) 
posteriors.agpe <- generate(
  fit,
  formula = ~ year.species2_latent,
  n.samples = n_draws) 
posteriors.elvi <- generate(
  fit,
  formula = ~ year.species3_latent,
  n.samples = n_draws) 

colnames(posteriors.aghy) <- c( paste0("iter",1:n_draws))
colnames(posteriors.agpe) <- c( paste0("iter",1:n_draws))
colnames(posteriors.elvi) <- c( paste0("iter",1:n_draws))

rownames(posteriors.aghy) <- c("year.aghy")
rownames(posteriors.agpe) <- c("year.agpe")
rownames(posteriors.elvi) <- c("year.elvi")


posteriors <- rbind(posteriors.aghy, posteriors.agpe, posteriors.elvi)



posteriors_df <- as_tibble(posteriors, rownames = "param") %>% 
  separate_wider_delim(param, delim = ".", names = c("param_type", "spp_label"), cols_remove = FALSE) %>% 
  pivot_longer( cols = -c(param, param_type, spp_label), names_to = "iteration")




posteriors_summary <- posteriors_df %>% 
  group_by(param, spp_label, param_type) %>% 
  summarize(mean = mean(value), 
            lwr = quantile(value, .025),
            upr = quantile(value, .975),
            prop_pos = sum(value>0)/500)

posterior_hist <- ggplot(posteriors_df)+
  stat_halfeye(aes(x = value, y = spp_label, fill = spp_label), breaks = 50, alpha = .6)+
  # stat_histinterval(aes(x = value, y = spp_label, fill = spp_label), breaks = 50, alpha = .6)+
  # geom_point(data = posteriors_summary, aes(x = mean, y = spp_label, color = spp_label))+
  # geom_linerange(data = posteriors_summary, aes(xmin = lwr, xmax = upr, y = spp_label, color = spp_label))+
  
  geom_vline(xintercept = 0)+
  facet_wrap(~param_type, scales = "free")+
  scale_color_manual(values = species_colors)+
  scale_fill_manual(values = species_colors)+
  scale_x_continuous(labels = scales::label_number(), guide = guide_axis(check.overlap = TRUE))+
  theme_bw()

posterior_hist
################################################################################################################################
##########  Plotting the spatially varying trends ###############
################################################################################################################################
# Making a mask for each species based on Jacob's host SDMs
aghy_raster <-terra::rast("aghy_binary.tif") 
aghy_raster_distribution <- mask(aghy_raster, mask = aghy_raster == 1, maskvalue = 0)

aghy_distribution_poly<-  terra::as.polygons(aghy_raster_distribution) 

aghy_distribution <- aghy_distribution_poly   %>% 
  st_as_sf(coords = c("lon", "lat"), crs = 4326, remove = FALSE) %>% 
  st_transform(epsg6703km)  




agpe_raster <- terra::rast("agpe_binary.tif")
agpe_raster_distribution <- mask(agpe_raster, mask = agpe_raster == 1, maskvalue = 0)

agpe_distribution_poly<-  terra::as.polygons(agpe_raster_distribution) 

agpe_distribution <- agpe_distribution_poly   %>% 
  st_as_sf(coords = c("lon", "lat"), crs = 4326, remove = FALSE) %>% 
  st_transform(epsg6703km)  


elvi_raster <- terra::rast("elvi_binary.tif")
elvi_raster_distribution <- mask(elvi_raster, mask = elvi_raster == 1, maskvalue = 0)

elvi_distribution_poly<-  terra::as.polygons(elvi_raster_distribution) 

elvi_distribution <- elvi_distribution_poly   %>% 
  st_as_sf(coords = c("lon", "lat"), crs = 4326, remove = FALSE) %>% 
  st_transform(epsg6703km)  


vrt_aghy <- inlabru::fm_pixels(mesh, mask = aghy_distribution, format = "sp", dims = c(60, 60))# Note this is where we can mask the output according the whatever shape, such as the host distribution
vrt_agpe <- inlabru::fm_pixels(mesh, mask = agpe_distribution, format = "sp", dims = c(60, 60))# Note this is where we can mask the output according the whatever shape, such as the host distribution
vrt_elvi <- inlabru::fm_pixels(mesh, mask = elvi_distribution, format = "sp", dims = c(60, 60))# Note this is where we can mask the output according the whatever shape, such as the host distribution

saveRDS(vrt_aghy, file = "aghy_distribution_df.rds")
saveRDS(vrt_agpe, file = "agpe_distribution_df.rds")
saveRDS(vrt_elvi, file = "elvi_distribution_df.rds")

vrt_aghy <- readRDS(file = "aghy_distribution_df.rds")
vrt_agpe <- readRDS(file = "agpe_distribution_df.rds")
vrt_elvi <- readRDS(file = "elvi_distribution_df.rds")


# ggplot()+
#   gg(mesh)+
#   gg(vrt_aghy, color = "red")


vrt_aghy@data <- expand.grid(std_year = rep(1, length.out = length(vrt_aghy)),
                             species_index = 1,
                             species = species_names[1])

vrt_agpe@data <- expand.grid(std_year = rep(1, length.out = length(vrt_agpe)),
                             species_index = 2,
                             species = species_names[2])

vrt_elvi@data <- expand.grid(std_year = rep(1, length.out = length(vrt_elvi)),
                             species_index = 3,
                             species = species_names[3])



svc.pred_aghy <- predict(fit, 
                         vrt_aghy, 
                         formula = ~ ( exp(year.species1+ time.species1)-1)*100)
svc.pred_agpe <- predict(fit, 
                         vrt_agpe, 
                         formula = ~ ( exp(year.species2+ time.species2)-1)*100)
svc.pred_elvi <- predict(fit, 
                         vrt_elvi, 
                         formula = ~ ( exp(year.species3+ time.species3)-1)*100)




saveRDS(svc.pred_aghy, file = "svc.pred_aghy.Rds")
saveRDS(svc.pred_agpe, file = "svc.pred_agpe.Rds")
saveRDS(svc.pred_elvi, file = "svc.pred_elvi.Rds")

svc.pred_aghy <- readRDS(file = "svc.pred_aghy.Rds")
svc.pred_agpe <- readRDS(file = "svc.pred_agpe.Rds")
svc.pred_elvi <- readRDS(file = "svc.pred_elvi.Rds")



# make a base map
world_map <- st_as_sf(maps::map("world", plot = FALSE, fill = TRUE)) %>%
  st_transform(epsg6703km)
states_map <- st_as_sf(maps::map("state", plot = FALSE, fill = TRUE)) %>%
  st_transform(epsg6703km)




min_trend <- max(svc.pred_aghy$mean, svc.pred_aghy$mean, svc.pred_aghy$mean)

max_trend <- min(svc.pred_aghy$mean, svc.pred_aghy$mean, svc.pred_aghy$mean)

trendrange <- range(svc.pred_aghy$mean, svc.pred_aghy$mean, svc.pred_aghy$mean)



space_x <- range(svc.pred_aghy@coords[,1],svc.pred_agpe@coords[,1],svc.pred_agpe@coords[,1])
space_y <- range(svc.pred_aghy@coords[,2],svc.pred_agpe@coords[,2],svc.pred_agpe@coords[,2])



svc_time_map_AGHY <- ggplot()+
  geom_sf(data = world_map, color = "grey", linewidth = .1, fill = "#FAF9F6") +
  geom_sf(data = states_map, color = "grey", linewidth = .1, fill = "#FAF9F6") +
  coord_sf(xlim = space_x, ylim = space_y)+
  gg(svc.pred_aghy, aes(fill = mean))+
  scale_fill_viridis_c(option = "turbo", na.value = "transparent", limits = trendrange)+
  labs(title = species_names[1], fill = "% change/year", y = "Latitude", x = "Longitude")+
  theme_light()+
  theme(plot.title = element_text(face = "italic"), aspect.ratio = 1)
# svc_time_map_AGHY


svc_time_map_AGPE <- ggplot()+
  geom_sf(data = world_map, color = "grey", linewidth = .1, fill = "#FAF9F6") +
  geom_sf(data = states_map, color = "grey", linewidth = .1, fill = "#FAF9F6") +
  coord_sf(xlim = space_x, ylim = space_y)+
  gg(svc.pred_agpe, aes(fill = mean))+
  scale_fill_viridis_c(option = "turbo", na.value = "transparent", limits = trendrange)+
  labs(title = species_names[2], fill = "% change/year", y = "Latitude", x = "Longitude")+
  theme_light()+
  theme(plot.title = element_text(face = "italic"), aspect.ratio = 1)
# svc_time_map_AGPE


svc_time_map_ELVI <- ggplot()+
  geom_sf(data = world_map, color = "grey", linewidth = .1, fill = "#FAF9F6") +
  geom_sf(data = states_map, color = "grey", linewidth = .1, fill = "#FAF9F6") +
  coord_sf(xlim = space_x, ylim = space_y)+
  gg(svc.pred_elvi, aes(fill = mean))+
  scale_fill_viridis_c(option = "turbo", na.value = "transparent", limits = trendrange)+
  labs(title = species_names[3], fill = "% change/year", y = "Latitude", x = "Longitude")+
  theme_light()+
  theme(plot.title = element_text(face = "italic"), aspect.ratio = 1)
# svc_time_map_ELVI



svc_time_map <- svc_time_map_AGHY + svc_time_map_AGPE + svc_time_map_ELVI + plot_layout(nrow = 1, guides = "collect") + plot_annotation(tag_levels = "A")
ggsave(svc_time_map, filename = "Plots/svc_time_map.png", width = 15, height = 5)



trendrange.CI <- range(svc.pred_aghy$`q0.975`-svc.pred_aghy$`q0.025`, svc.pred_agpe$`q0.975`-svc.pred_agpe$`q0.025`, svc.pred_elvi$`q0.975`-svc.pred_elvi$`q0.025`)




svc_time_map_AGHY.CI <- ggplot()+
  geom_sf(data = world_map, color = "grey", linewidth = .1, fill = "#FAF9F6") +
  geom_sf(data = states_map, color = "grey", linewidth = .1, fill = "#FAF9F6") +
  coord_sf(xlim = space_x, ylim = space_y)+
  gg(svc.pred_aghy, aes(fill = `q0.975`-`q0.025`))+
  scale_fill_gradient(low = "white", high = "deeppink4", na.value = "transparent", limits = trendrange.CI)+
  labs(title = species_names[1], y = "Latitude", x = "Longitude")+
  theme_light()+
  theme(plot.title = element_text(face = "italic"), aspect.ratio = 1)
# svc_time_map_AGHY.CI


svc_time_map_AGPE.CI <- ggplot()+
  geom_sf(data = world_map, color = "grey", linewidth = .1, fill = "#FAF9F6") +
  geom_sf(data = states_map, color = "grey", linewidth = .1, fill = "#FAF9F6") +
  coord_sf(xlim = space_x, ylim = space_y)+
  gg(svc.pred_agpe, aes(fill = `q0.975`-`q0.025`))+
  scale_fill_gradient(low = "white", high = "deeppink4", na.value = "transparent",  limits = trendrange.CI)+
  labs(title = species_names[2], y = "Latitude", x = "Longitude")+
  theme_light()+
  theme(plot.title = element_text(face = "italic"), aspect.ratio = 1)
# svc_time_map_AGPE.CI


svc_time_map_ELVI.CI <- ggplot()+
  geom_sf(data = world_map, color = "grey", linewidth = .1, fill = "#FAF9F6") +
  geom_sf(data = states_map, color = "grey", linewidth = .1, fill = "#FAF9F6") +
  coord_sf(xlim = space_x, ylim = space_y)+
  gg(svc.pred_elvi, aes(fill = `q0.975`-`q0.025`))+
  scale_fill_gradient(low = "white", high = "deeppink4", na.value = "transparent",  limits = trendrange.CI)+
  labs(title = species_names[3], y = "Latitude", x = "Longitude")+
  theme_light()+
  theme(plot.title = element_text(face = "italic"), aspect.ratio = 1)
# svc_time_map_ELVI.CI



svc_time_map.CI <- svc_time_map_AGHY.CI + svc_time_map_AGPE.CI + svc_time_map_ELVI.CI + plot_layout(nrow = 1, guides = "collect") + plot_annotation(tag_levels = "A")
ggsave(svc_time_map.CI, filename = "Plots/svc_time_map_CI.png", width = 15, height = 5)



################################################################################################################################
##########  Plotting the spatial Intercepts ###############
################################################################################################################################


# generating predicted prevalence in 1895
vrt_aghy@data <- expand.grid(std_year = rep(1895-mean_year, length.out = length(vrt_aghy)),
                             species_index = 1,
                             species = species_names[1])

vrt_agpe@data <- expand.grid(std_year = rep(1895-mean_year, length.out = length(vrt_agpe)),
                             species_index = 2,
                             species = species_names[2])

vrt_elvi@data <- expand.grid(std_year = rep(1895-mean_year, length.out = length(vrt_elvi)),
                             species_index = 3,
                             species = species_names[3])



svc.space.pred_aghy.1895 <- predict(fit, 
                         vrt_aghy, 
                         formula = ~ invlogit(int.species1 + space.species1 + year.species1+ time.species1))
svc.space.pred_agpe.1895 <- predict(fit, 
                         vrt_agpe, 
                         formula = ~ invlogit(int.species2 + space.species2 + year.species2+ time.species2))
svc.space.pred_elvi.1895 <- predict(fit, 
                         vrt_elvi, 
                         formula = ~ invlogit(int.species3 + space.species3 + year.species3+ time.species3))

# generating predicted prevalence in 2020

vrt_aghy@data <- expand.grid(std_year = rep(2020-mean_year, length.out = length(vrt_aghy)),
                             species_index = 1,
                             species = species_names[1])

vrt_agpe@data <- expand.grid(std_year = rep(2020-mean_year, length.out = length(vrt_agpe)),
                             species_index = 2,
                             species = species_names[2])

vrt_elvi@data <- expand.grid(std_year = rep(2020-mean_year, length.out = length(vrt_elvi)),
                             species_index = 3,
                             species = species_names[3])



svc.space.pred_aghy.2020 <- predict(fit, 
                                    vrt_aghy, 
                                    formula = ~ invlogit(int.species1 + space.species1 + year.species1+ time.species1))
svc.space.pred_agpe.2020 <- predict(fit, 
                                    vrt_agpe, 
                                    formula = ~ invlogit(int.species2 + space.species2 + year.species2+ time.species2))
svc.space.pred_elvi.2020 <- predict(fit, 
                                    vrt_elvi, 
                                    formula = ~ invlogit(int.species3 + space.species3 + year.species3+ time.species3))



# saving predictions
saveRDS(svc.space.pred_aghy.1895, file = "svc.space.pred_aghy.1895.Rds")
saveRDS(svc.space.pred_agpe.1895, file = "svc.space.pred_agpe.1895.Rds")
saveRDS(svc.space.pred_elvi.1895, file = "svc.space.pred_elvi.1895.Rds")

svc.space.pred_aghy.1895 <- readRDS(file = "svc.space.pred_aghy.1895.Rds")
svc.space.pred_agpe.1895 <- readRDS(file = "svc.space.pred_agpe.1895.Rds")
svc.space.pred_elvi.1895 <- readRDS(file = "svc.space.pred_elvi.1895.Rds")

saveRDS(svc.space.pred_aghy.2020, file = "svc.space.pred_aghy.2020.Rds")
saveRDS(svc.space.pred_agpe.2020, file = "svc.space.pred_agpe.2020.Rds")
saveRDS(svc.space.pred_elvi.2020, file = "svc.space.pred_elvi.2020.Rds")

svc.space.pred_aghy.2020 <- readRDS(file = "svc.space.pred_aghy.2020.Rds")
svc.space.pred_agpe.2020 <- readRDS(file = "svc.space.pred_agpe.2020.Rds")
svc.space.pred_elvi.2020 <- readRDS(file = "svc.space.pred_elvi.2020.Rds")




prevrange <- range(svc.space.pred_aghy.1895$mean, svc.space.pred_agpe.1895$mean, svc.space.pred_elvi.1895$mean,
                   svc.space.pred_aghy.2020$mean, svc.space.pred_agpe.2020$mean, svc.space.pred_elvi.2020$mean)


space_x <- range(svc.space.pred_aghy.1895@coords[,1],svc.space.pred_agpe.1895@coords[,1],svc.space.pred_agpe.1895@coords[,1])
space_y <- range(svc.space.pred_aghy.1895@coords[,2],svc.space.pred_agpe.1895@coords[,2],svc.space.pred_agpe.1895@coords[,2])

svc_space_map_AGHY.1895 <- ggplot()+
  geom_sf(data = world_map, color = "grey", linewidth = .1, fill = "#FAF9F6") +
  geom_sf(data = states_map, color = "grey", linewidth = .1, fill = "#FAF9F6") +
  coord_sf(xlim = space_x, ylim = space_y)+
  gg(svc.space.pred_aghy.1895, aes(fill = mean))+
  scale_fill_viridis_c(option = "turbo", na.value = "transparent", limits = prevrange)+
  labs(title = species_names[1], subtitle = "year: 1895", fill = "% prevalence", y = "Latitude", x = "Longitude")+
  theme_light()+
  theme(plot.title = element_text(face = "italic"))
svc_space_map_AGHY.1895

svc_space_map_AGHY.2020 <- ggplot()+
  geom_sf(data = world_map, color = "grey", linewidth = .1, fill = "#FAF9F6") +
  geom_sf(data = states_map, color = "grey", linewidth = .1, fill = "#FAF9F6") +
  coord_sf(xlim = space_x, ylim = space_y)+
  gg(svc.space.pred_aghy.2020, aes(fill = mean))+
  scale_fill_viridis_c(option = "turbo", na.value = "transparent", limits = prevrange)+
  labs(title = species_names[1], subtitle = "year: 2020", fill = "% prevalence", y = "Latitude", x = "Longitude")+
  theme_light()+
  theme(plot.title = element_text(face = "italic"))
svc_space_map_AGHY.2020

svc_space_map_AGPE.1895<- ggplot()+
  geom_sf(data = world_map, color = "grey", linewidth = .1, fill = "#FAF9F6") +
  geom_sf(data = states_map, color = "grey", linewidth = .1, fill = "#FAF9F6") +
  gg(svc.space.pred_agpe.1895, aes(fill = mean))+
  coord_sf(xlim = space_x, ylim = space_y)+
  scale_fill_viridis_c(option = "turbo", na.value = "transparent", limits = prevrange)+
  labs(title = species_names[2], subtitle = "year: 1895", fill = "% prevalence", y = "Latitude", x = "Longitude")+
  theme_light()+
  theme(plot.title = element_text(face = "italic"))

svc_space_map_AGPE.1895

svc_space_map_AGPE.2020<- ggplot()+
  geom_sf(data = world_map, color = "grey", linewidth = .1, fill = "#FAF9F6") +
  geom_sf(data = states_map, color = "grey", linewidth = .1, fill = "#FAF9F6") +
  gg(svc.space.pred_agpe.2020, aes(fill = mean))+
  coord_sf(xlim = space_x, ylim = space_y)+
  scale_fill_viridis_c(option = "turbo", na.value = "transparent", limits = prevrange)+
  labs(title = species_names[2], subtitle = "year: 2020", fill = "% prevalence", y = "Latitude", x = "Longitude")+
  theme_light()+
  theme(plot.title = element_text(face = "italic"))

svc_space_map_AGPE.2020


svc_space_map_ELVI.1895<- ggplot()+
  geom_sf(data = world_map, color = "grey", linewidth = .1, fill = "#FAF9F6") +
  geom_sf(data = states_map, color = "grey", linewidth = .1, fill = "#FAF9F6") +
  gg(svc.space.pred_elvi.1895, aes(fill = mean))+
  coord_sf(xlim = space_x, ylim = space_y)+
  scale_fill_viridis_c(option = "turbo", na.value = "transparent", limits = prevrange)+
  labs(title = species_names[3], subtitle = "year: 1895", fill = "% prevalence", y = "Latitude", x = "Longitude")+
  theme_light()+
  theme(plot.title = element_text(face = "italic"))
svc_space_map_ELVI.1895

svc_space_map_ELVI.2020<- ggplot()+
  geom_sf(data = world_map, color = "grey", linewidth = .1, fill = "#FAF9F6") +
  geom_sf(data = states_map, color = "grey", linewidth = .1, fill = "#FAF9F6") +
  gg(svc.space.pred_elvi.2020, aes(fill = mean))+
  coord_sf(xlim = space_x, ylim = space_y)+
  scale_fill_viridis_c(option = "turbo", na.value = "transparent", limits = prevrange)+
  labs(title = species_names[3], subtitle = "year: 2020", fill = "% prevalence", y = "Latitude", x = "Longitude")+
  theme_light()+
  theme(plot.title = element_text(face = "italic"))
svc_space_map_ELVI.2020


svc_space_map <- (svc_space_map_AGHY.1895 + svc_space_map_AGPE.1895 + svc_space_map_ELVI.1895) / (svc_space_map_AGHY.2020 + svc_space_map_AGPE.2020 + svc_space_map_ELVI.2020)+
  plot_layout(nrow = 2, guides = 'collect') + plot_annotation(tag_levels = "A")
ggsave(svc_space_map, filename = "Plots/svc_space_map_year.png", width = 12, height = 8)


# evaluating whether or not change in prevalence is associated with starting prevalence
vrt_aghy@data <- data.frame("trend.mean" = svc.pred_aghy$mean, "prev.mean" = svc.space.pred_aghy$mean)
vrt_agpe@data <- data.frame("trend.mean" = svc.pred_agpe$mean, "prev.mean" = svc.space.pred_agpe$mean)
vrt_elvi@data <- data.frame("trend.mean" = svc.pred_elvi$mean, "prev.mean" = svc.space.pred_elvi$mean)

initialprev.trend_aghy <- ggplot(data = vrt_aghy@data)+
  geom_smooth(aes(x = prev.mean, y = trend.mean), color = "black", method = 'lm')+
  geom_point(aes(x = prev.mean, y = trend.mean), alpha = .6)+
  labs(title = species_names[1], x = "% prevalence in 1895", y = "% change/year")+
  theme_light()+
  theme(plot.title = element_text(face = "italic"))
initialprev.trend_aghy

initialprev.trend_agpe <- ggplot(data = vrt_agpe@data)+
  geom_smooth(aes(x = prev.mean, y = trend.mean), color = "black", method = 'lm')+
  geom_point(aes(x = prev.mean, y = trend.mean), alpha = .6)+
  labs(title = species_names[2], x = "% prevalence in 1895", y = "% change/year")+
  theme_light()+
  theme(plot.title = element_text(face = "italic"))
initialprev.trend_agpe

initialprev.trend_elvi <- ggplot(data = vrt_elvi@data)+
  geom_smooth(aes(x = prev.mean, y = trend.mean),color = "black", method = 'lm')+
  geom_point(aes(x = prev.mean, y = trend.mean), alpha = .6)+
  labs(title = species_names[3], x = "% prevalence in 1895", y = "% change/year")+
  theme_light()+
  theme(plot.title = element_text(face = "italic"))
initialprev.trend_elvi


initialprev.trend_plot <- initialprev.trend_aghy+initialprev.trend_agpe+initialprev.trend_elvi+
  plot_layout(nrow = 1, axis_titles = "collect", guides = 'collect') + plot_annotation( tag_levels = "A")
ggsave(initialprev.trend_plot, filename = "Plots/initialprev_trend_plot.png", width = 8, height = 4)

################################################################################################################################
##########  Making the contemp predictions ###############
################################################################################################################################

contemp.aghy <- contemp_surveys %>% 
  filter(SpeciesID == "AGHY") %>% 
  mutate(herbarium_index = 9999, collector_index = 9999, scorer_index = NA, species_index = 1, species = species_names[1]) %>% 
  st_transform(epsg6703km)

contemp.elvi <- contemp_surveys %>% 
  filter(SpeciesID == "ELVI") %>% 
  mutate(herbarium_index = 9999, collector_index = 9999, scorer_index = NA, species_index = 3, species = species_names[3]) %>% 
  st_transform(epsg6703km)


# gennerating predictions and back-transforming the standardized year variable

mean_year <- mean(endo_herb$year)
sd_year <- sd(endo_herb$year)


contemp.pred.aghy<- predict(
  fit,
  newdata = contemp.aghy,
  formula = ~ invlogit(int.species1 + year.species1 + space + space.species1 + time.species1 + 
                         scorer_eval(scorer_index) + collector_eval(collector_index))) %>% 
  mutate(year = std_year + mean_year)


contemp.pred.elvi<- predict(
  fit,
  newdata = contemp.elvi,
  formula = ~ invlogit(int.species3 + year.species3 + space.species3 + time.species3 + 
                         scorer_eval(scorer_index) + collector_eval(collector_index))) %>% 
  mutate(year = std_year + mean_year)


contemp.pred <- bind_rows(contemp.pred.aghy, contemp.pred.elvi)




contemp_lon <- ggplot(contemp.pred)+
  geom_point(aes(x = lon, y = endo_prev, size = sample_size), alpha = .4)+
  geom_smooth(aes(x = lon, y = endo_prev, group = species), color = "black", method = "glm",  formula = "y ~ x", method.args = list(family = "binomial" ))+
  geom_linerange(aes(x = lon, ymin = `q0.025`, ymax = `q0.975`, color = species), alpha = .8)+
  geom_point(aes(x = lon, y = mean), shape = 4) + 
  scale_color_manual(values = c(species_colors[1], species_colors[3]))+
  guides(color = "none")+
  lims(x = c(min(contemp.aghy$lon, contemp.elvi$lon),max(contemp.aghy$lon, contemp.elvi$lon)), y = c(0,1)) +
  labs(y = "Endophyte Prevalance", x = "Longitude", color = "Species", size = "Sample Size")+
  facet_wrap(~species, ncol=1, scales = "free") + theme_classic()+ theme(strip.background = element_blank(),
                                                                         strip.text = element_blank())
# contemp_lon



contemp_lat <- ggplot(contemp.pred)+
  geom_point(aes(x = lat, y = endo_prev, size = sample_size), alpha = .4)+
  geom_smooth(aes(x = lat, y = endo_prev, group = species), color = "black", method = "glm",  formula = "y ~ x", method.args = list(family = "binomial" ))+
  geom_linerange(aes(x = lat, ymin = `q0.025`, ymax = `q0.975`, color = species))+
  geom_point(aes(x = lat, y = mean), shape = 4) +
  scale_color_manual(values = c(species_colors[1], species_colors[3]))+
  guides(color = "none")+
  lims(x = c(min(contemp.aghy$lat, contemp.elvi$lat),max(contemp.aghy$lat, contemp.elvi$lat)), y = c(0,1)) +
  labs(y = "Endophyte Prevalance", x = "Latitude",  color = "Species", size = "Sample Size")+
  facet_wrap(~species, ncol = 1, scales = "free") + theme_classic() + theme(strip.background = element_blank(),
                                                                            strip.text = element_blank())
# contemp_lat

contemp_obspred <- ggplot(contemp.pred)+
  geom_abline(intercept = 0, slope = 1)+
  geom_linerange(aes(y = endo_prev, xmin = `q0.025`, xmax = `q0.975`, color = species))+
  geom_point(aes(x = mean, y = endo_prev), shape = 4)+
  scale_color_manual(values = c(species_colors[1], species_colors[3]))+
  guides(color = "none")+
  lims(x = c(0,1), y = c(0,1)) + labs(y = "Observed", x = "Predicted",  color = "Species")+
  facet_wrap(~species, ncol = 1, scales = "free", strip.position = "right") + theme_classic()+ theme(strip.background = element_blank(),
                                                                           strip.text = element_text(size = rel(1)),
                                                                           strip.text.y = element_text(face = "italic", angle = 0))
# contemp_obspred



contemp_test_plot <- contemp_lon + contemp_lat+contemp_obspred +
  plot_layout(nrow = 1, guides = "collect") + plot_annotation(tag_levels = "A")
# contemp_test_plot
ggsave(contemp_test_plot, filename = "Plots/contemp_test_plot.png", width = 12, height = 5)
###

# now looking at the ROC and AUC values for the contemporary dataset choosing only one plant from each population
# however we only have this information for AGHY
rocobj <- pROC::roc(contemp_random_sample$endo_status, contemp.pred.aghy$mean)

ROC_test_plot <- ggroc(rocobj) 
ggsave(ROC_test_plot, filename = "Plots/ROC_test_plot.png", width = 4, height = 4)


# AUC values
rocobj$auc
#0.7463










###### Making a map to show where the contemp predictions ###########
contemp_surveys <- contemp_surveys %>% 
  mutate(species_name = case_when(SpeciesID == "AGHY" ~ "A. hyemalis",
                                  SpeciesID == "ELVI" ~ "E. virginicus"))
contemp_map <- ggplot()+
  geom_map(data = outline_map, map = outline_map, aes(map_id = region), color = "grey", linewidth = .1, fill = "#FAF9F6")+
  geom_map(data = states_shape, map = states_shape, aes(map_id = region), color = "grey", linewidth = .1, fill = NA)+
  geom_point(data = endo_herb, aes(x = lon, y = lat),shape = 4, alpha = .3, size = 1.2)+
  geom_point(data = contemp_surveys, aes(x = lon, y = lat, color = species_name))+
  coord_sf(xlim = c(-109,-68), ylim = c(19,49), crs = 4326)+
  scale_color_manual(values = c(species_colors[1], species_colors[3]))+
  theme_light()+
  theme(legend.text = element_text(face = "italic"))+
  labs(x = "Longitude", y = "Latitude", color = "Host Species")
contemp_map
ggsave(contemp_map, filename = "contemp_surveys_map.png", width = 7, height = 4)


