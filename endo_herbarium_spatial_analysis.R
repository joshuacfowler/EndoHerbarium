# Purpose: Fits spatial model of change in endophyte prevalence in herbarium specimens
# Authors: Joshua Fowler, Mallory Tucker, Ella Segal, Lani Dufresne, and Tom Miller
# Updated: Noov 11, 2022

library(tidyverse) # for data manipulation and ggplot
library(INLA) # for fitting integrated nested Laplace approximation models
library(inlabru)

library(sf)
library(rmapshaper)

library(patchwork)
library(ggmap)
library(ROCR)


invlogit<-function(x){exp(x)/(1+exp(x))}
species_colors <- c("#1b9e77","#d95f02","#7570b3")
endophyte_colors <- c("#fdedd3","#f3c8a8", "#5a727b", "#4986c7", "#181914",  "#163381")


################################################################################
############ Read in the herbarium dataset ############################### 
################################################################################

endo_herb_georef <- read_csv(file = "~/Dropbox/Josh&Tom - shared/Endo_Herbarium/DigitizedHerbariumRecords/endo_herb_georef.csv") %>%
  # filter(Country != "Canada") %>%
  mutate(Spp_code = case_when(grepl("AGHY", Sample_id) ~ "AGHY",
                              grepl("ELVI", Sample_id) ~ "ELVI",
                              grepl("AGPE", Sample_id) ~ "AGPE")) %>% 
  mutate(species = case_when(Spp_code == "AGHY" ~ "A. hyemalis",
                             Spp_code == "AGPE" ~ "A. perennans",
                             Spp_code == "ELVI" ~ "E. virginicus")) 

endo_herb_AGHY <- endo_herb_georef %>% 
  filter(!is.na(Endo_status_liberal)) %>% 
  filter(Spp_code == "AGHY") %>% 
  filter(!is.na(lon) & !is.na(year)) %>% 
  filter(lon>-110) %>% 
  filter(Country != "Canada" & !is.na(County)) 
  
endo_herb_AGPE <- endo_herb_georef %>% 
  filter(!is.na(Endo_status_liberal)) %>% 
  filter(Spp_code == "AGPE") %>% 
  filter(!is.na(lon) & !is.na(year)) %>% 
  filter(lon>-100) %>% 
  filter(Country != "Canada" & !is.na(County))  
  

endo_herb_ELVI <- endo_herb_georef %>% 
  filter(!is.na(Endo_status_liberal)) %>% 
  filter(Spp_code == "ELVI") %>% 
  filter(!is.na(lon) & !is.na(year)) %>% 
  filter(lon>-110) %>% 
  filter(Country != "Canada" & !is.na(County))
  

endo_herb <- endo_herb_georef %>% 
  filter(!is.na(Endo_status_liberal)) %>%
  filter(!is.na(Spp_code)) %>% 
  filter(!is.na(lon) & !is.na(year)) %>% 
  filter(lon>-110 ) %>% 
  filter(Country != "Canada" ) %>% 
  filter(!is.na(County) | is.na(Municipality)) %>% 
  mutate(year_bin = case_when(year<1970 ~ "pre-1970",
                              year>=1970 ~ "post-1970")) %>% 
  mutate(scorer_id = case_when(scorer_id == "BellaGuttierez" ~ "BellaGutierrez",
                               TRUE ~ scorer_id)) %>% 
  mutate(endo_status_text = case_when(Endo_status_liberal == 0 ~ "E-",
                                      Endo_status_liberal == 1 ~ "E+"))  
  # mutate(Endo_status_liberal = case_when(Spp_code == "AGPE" & Endo_status_liberal ==0 ~ 1,
  #                                         Spp_code == "AGPE" & Endo_status_liberal == 1~ 0, TRUE ~ Endo_status_liberal))




# updating the scorer labels
scorer_levels <- levels(as.factor(endo_herb$scorer_id))
scorer_no <- paste0("Scorer",1:nlevels(as.factor(endo_herb$scorer_id)))

endo_herb$scorer_id <- scorer_no[match(as.factor(endo_herb$scorer_id), scorer_levels)]

# updating the collector labels

endo_herb$collector_lastname <- str_replace_all(endo_herb$collector_lastname, "�", " ")
endo_herb$collector_lastname <- str_replace_all(endo_herb$collector_lastname, "xa0", " ")


collector_levels <- levels(as.factor(endo_herb$collector_lastname))
collector_no <- paste0("Collector",1:nlevels(as.factor(endo_herb$collector_lastname)))

endo_herb$collector_id <- collector_no[match(as.factor(endo_herb$collector_lastname), collector_levels)]


# endo_herb <- endo_herb_AGHY
# endo_herb <- endo_herb_AGPE
# endo_herb <- endo_herb_ELVI




# Loading in contemporary survey data for AGHY and ELVI, which we will use for model validation
contemp_surveys <- read_csv(file = "contemp_surveys.csv")
contemp_random_sample <- read_csv(file = "contemp_random_sample.csv")

contemp_random_sample <- contemp_random_sample %>% 
  left_join(contemp_surveys)


# register_google(key = "")
# map <- ggmap::get_map(zoom = 3, source = "google", maptype = c("satellite"))

outline_map <- map_data("world")
states_shape <- map_data("state")
counties_shape <- map_data("county")
# ggplot()+
#   geom_map(data = counties_shape, map = counties_shape, aes(long, lat, map_id = region))

collections_map <- ggplot()+
  geom_map(data = outline_map, map = outline_map, aes(long, lat, map_id = region), color = "grey", linewidth = .1, fill = "#FAF9F6")+
  geom_map(data = states_shape, map = states_shape, aes(long, lat, map_id = region), color = "grey", linewidth = .1, fill = NA)+
  geom_point(data = endo_herb, aes(x = lon, y = lat, color = species), alpha = .7, lwd = .5)+
  coord_sf(xlim = c(-109,-68), ylim = c(21,49))+
  scale_color_manual(values = species_colors)+
  theme_light()+
  theme(legend.text = element_text(face = "italic"))+
  labs(x = "Longitude", y = "Latitude", color = "Host Species")

collections_map
ggsave(collections_map, filename = "collections_map.png", width = 7, height = 4)

endo_status_map <- ggplot()+
  geom_map(data = outline_map, map = outline_map, aes(long, lat, map_id = region), color = "grey", linewidth = .1, fill = "#FAF9F6")+
  geom_map(data = states_shape, map = states_shape, aes(long, lat, map_id = region), color = "grey", linewidth = .1, fill = NA)+
  geom_point(data = endo_herb, aes(x = lon, y = lat, color = endo_status_text), alpha = .7, lwd = .8)+
  facet_wrap(~factor(year_bin, levels = c("pre-1970", "post-1970"))+species)+
  coord_sf(xlim = c(-109,-68), ylim = c(21,49))+
  scale_color_manual(values = c(endophyte_colors[2],endophyte_colors[6]))+
  theme_light()+
  theme(strip.background = element_blank(),
        strip.text  = element_text(face = "italic", color = "black"))+
  labs(x = "Longitude", y = "Latitude", color = "Endophyte Status")
endo_status_map

ggsave(endo_status_map, filename = "endo_status_map.png", width = 10, height = 5)



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
  separate(Sample_id_temp, sep = "_", into = c("herbarium", "spp_code", "plant_no")) %>% select(-spp_code, -plant_no) %>% 
  # filter(seed_scored>0) %>% 
  filter(month<=12&month>0) %>% 
  group_by(species) %>% 
  summarize(n(),
            avg_seed = mean(seed_scored, na.rm = T),
            avg_month = mode(as.numeric(month)))


# summary(lm(formula(Endo_status_liberal_1 ~ 0 + Spp_code + Spp_code:year + Spp_code:year:lat + Spp_code:year:lon + Spp_code:year:lat:lon), data = endo_herb))
################################################################################
############ Setting up and running INLA model ############################### 
################################################################################
# messing aruond with INLA model fitting
# generate mesh, which is used to fit the spatial effects
# using a coarser mesh to start out as this is easier computation and sufficient to capture spatial effects
coords <- cbind(endo_herb$lon, endo_herb$lat)
AGHY_coords <- cbind(endo_herb_AGHY$lon, endo_herb_AGHY$lat)
AGPE_coords <- cbind(endo_herb_AGPE$lon, endo_herb_AGPE$lat)
ELVI_coords <- cbind(endo_herb_ELVI$lon, endo_herb_ELVI$lat)

non_convex_bdry <- inla.nonconvex.hull(coords, -0.03, -0.05, resolution = c(100, 100))
AGHY_non_convex_bdry <- inla.nonconvex.hull(AGHY_coords, -0.03, -0.05, resolution = c(100, 100))
AGPE_non_convex_bdry <- inla.nonconvex.hull(AGPE_coords, -0.04, -0.05, resolution = c(100, 100))
ELVI_non_convex_bdry <- inla.nonconvex.hull(ELVI_coords, -0.03, -0.05, resolution = c(100, 100))

sf::sf_use_s2(FALSE)
bdry_st <- st_make_valid(as_tibble(non_convex_bdry$loc)  %>% 
  mutate(lon = V1,  lat = V2) %>% 
  st_as_sf(coords = c("lon", "lat"), crs = 4326) %>% 
  summarise(geometry = st_combine(geometry)) %>% 
  st_cast("POLYGON"))

AGHY_bdry_st <- st_make_valid(as_tibble(AGHY_non_convex_bdry$loc)  %>% 
                           mutate(lon = V1,  lat = V2) %>% 
                           st_as_sf(coords = c("lon", "lat"), crs = 4326) %>% 
                           summarise(geometry = st_combine(geometry)) %>% 
                           st_cast("POLYGON"))
AGPE_bdry_st <- st_make_valid(as_tibble(AGPE_non_convex_bdry$loc)  %>% 
                                mutate(lon = V1,  lat = V2) %>% 
                                st_as_sf(coords = c("lon", "lat"), crs = 4326) %>% 
                                summarise(geometry = st_combine(geometry)) %>% 
                                st_cast("POLYGON"))
ELVI_bdry_st <- st_make_valid(as_tibble(ELVI_non_convex_bdry$loc)  %>% 
                                mutate(lon = V1,  lat = V2) %>% 
                                st_as_sf(coords = c("lon", "lat"), crs = 4326) %>% 
                                summarise(geometry = st_combine(geometry)) %>% 
                                st_cast("POLYGON"))

coastline <- st_make_valid(sf::st_as_sf(maps::map("world", regions = c("usa", "canada", "mexico"), plot = FALSE, fill = TRUE)))
# plot(coastline)

bdry <- st_intersection(coastline$geom, bdry_st)
AGHY_bdry <- st_intersection(coastline$geom, AGHY_bdry_st)
AGPE_bdry <- st_intersection(coastline$geom, AGPE_bdry_st)
ELVI_bdry <- st_intersection(coastline$geom, ELVI_bdry_st)

bdry_polygon <- st_cast(st_sf(bdry), "MULTIPOLYGON", group_or_split = TRUE) %>% st_union() %>% 
  as("Spatial")
AGHY_bdry_polygon <- st_cast(st_sf(AGHY_bdry), "MULTIPOLYGON", group_or_split = TRUE) %>% st_union() %>% 
  as("Spatial")
AGPE_bdry_polygon <- st_cast(st_sf(AGPE_bdry), "MULTIPOLYGON", group_or_split = TRUE) %>% st_union() %>% 
  as("Spatial")
ELVI_bdry_polygon <- st_cast(st_sf(ELVI_bdry), "MULTIPOLYGON", group_or_split = TRUE) %>% st_union() %>% 
  as("Spatial")
# plot(bdry_polygon)

max.edge = diff(range(coords[,2]))/(10)
mesh10 <- inla.mesh.2d(loc = coords, max.edge = c(1,2)*max.edge,
                      boundary = bdry_polygon,
                      offset = c(1,4),
                      cutoff = max.edge/(5))
# plot(mesh10)
mesh5 <- inla.mesh.2d(loc = coords, max.edge = c(1,2)*(max.edge/2),
                       boundary = bdry_polygon,
                       offset = c(1,4),
                       cutoff = max.edge/(5))
# plot(mesh5)

mesh2.5 <- inla.mesh.2d(loc = coords, max.edge = c(1,2)*(max.edge/2/2),
                        boundary = bdry_polygon,
                        offset = c(1,4),
                        cutoff = max.edge/(10))
# plot(mesh2.5)

mesh1 <- inla.mesh.2d(loc = coords, max.edge = c(.5,2)*(max.edge/2/2),
                        boundary = bdry_polygon,
                        offset = c(1,4),
                        cutoff = max.edge/(10))
# plot(mesh1)

mesh <- inla.mesh.2d(loc = coords, max.edge = c(.5,2)*(max.edge/2/2),
                      boundary = bdry_polygon,
                      offset = c(1,4),
                      cutoff = max.edge/(20))

# mesh_plot <- ggplot()+
#   gg(mesh, exterior = FALSE, edge.color = "grey25", alpha = .1 )+
#   geom_point(data = endo_herb, aes(x = lon, y = lat), lwd = .5, color = "red")+
#   coord_sf()+
#   theme_minimal()
# ggsave(mesh_plot, filename = "mesh_plot.png", width = 7, height = 4)


mesh <- mesh
# points(coords, col = "red")
mesh$n # the number of mesh vertices


# defining the spatial random effect
A <- inla.spde.make.A(mesh, loc = coords)


spde <- inla.spde2.pcmatern(mesh = mesh,
                            prior.range = c(100, 0.01), # P(practic.range < 0.05) = 0.01
                            prior.sigma = c(100, 0.01)) # P(sigma > 1) = 0.01
spde.index <- inla.spde.make.index(name = "spatial.field", n.spde = spde$n.spde)

stack.fit <- inla.stack(data = list(Endo_status_liberal = endo_herb$Endo_status_liberal),
                    A = list(A,1),
                    effects = list(s = spde.index$spatial.field, 
                                   data.frame(Intercept = rep(1, nrow(endo_herb)),
                                              year = endo_herb$year,
                                              year2 = endo_herb$year^2,
                                              lat = endo_herb$lat,
                                              lat2 = endo_herb$lat^2,
                                              lon = endo_herb$lon,
                                              lon2 = endo_herb$lon^2,
                                              species = endo_herb$Spp_code,
                                              id = endo_herb$Sample_id,
                                              county = endo_herb$County,
                                              collector = endo_herb$collector_id,
                                              scorer = endo_herb$scorer_id)),
                    
                    tag = 'fit')
# making a data stack for the predicted value of endophyte prevalence
# need an even  grid of coordinates

# this is the set of data for which we want predictions across the whole space
pred_data <- data.frame(expand.grid(Intercept = 1, 
                        lon = seq(min(endo_herb$lon),max(endo_herb$lon), length.out = 100),
                        lat = seq(min(endo_herb$lat),max(endo_herb$lat), length.out = 100),
                        year = c(1870, 1895, 1920, 1925, 1970, 1990, 2020), 
                        # species = c("ELVI"),
                        species = c("AGHY", "AGPE", "ELVI"),
                        id = NA,
                        county = NA,
                        collector = NA,
                        scorer = NA)) %>% 
  mutate(year2 = year^2, lon2 = lon^2, lat2 = lat^2)


# keeping just the points that are part of the mesh's boundary
AGHY_simple_bdry <- st_as_sf(AGHY_bdry_polygon) %>% ms_simplify()  %>% st_make_valid() %>% as("Spatial") #%>% ms_filter_islands(min_area = 5000000000000)
AGPE_simple_bdry <- st_as_sf(AGPE_bdry_polygon) %>% ms_simplify()  %>% st_make_valid() %>% as("Spatial") #%>% ms_filter_islands(min_area = 5000000000000)
ELVI_simple_bdry <- st_as_sf(ELVI_bdry_polygon) %>% ms_simplify()  %>% st_make_valid() %>% as("Spatial") #%>% ms_filter_islands(min_area = 5000000000000)

# plot(simple_bdry)
aghy_ind <- point.in.polygon(
  pred_data$lon, pred_data$lat,
  raster::geom(AGHY_simple_bdry)[,5], raster::geom(AGHY_simple_bdry)[,6] 
)
agpe_ind <- point.in.polygon(
  pred_data$lon, pred_data$lat,
  raster::geom(AGPE_simple_bdry)[,5], raster::geom(AGPE_simple_bdry)[,6] 
)
elvi_ind <- point.in.polygon(
  pred_data$lon, pred_data$lat,
  raster::geom(ELVI_simple_bdry)[,5], raster::geom(ELVI_simple_bdry)[,6] 
)

aghy_pred_data <- pred_data %>% 
  filter(species == "AGHY")
aghy_coords_pred <- aghy_pred_data %>% 
  select(lon, lat) 
aghy_coords_pred <- as.matrix(aghy_coords_pred[which(aghy_ind == 1),])
# plot(coords_pred)
AGHY_pred_data <- aghy_pred_data[which(aghy_ind ==1),]

agpe_pred_data <- pred_data %>% 
  filter(species == "AGPE")
agpe_coords_pred <- agpe_pred_data %>% 
  select(lon, lat) 
agpe_coords_pred <- as.matrix(agpe_coords_pred[which(agpe_ind == 1),])
# plot(coords_pred)
AGPE_pred_data <- agpe_pred_data[which(agpe_ind ==1),]

elvi_pred_data <- pred_data %>% 
  filter(species == "ELVI")
elvi_coords_pred <- elvi_pred_data %>% 
  select(lon, lat) 
elvi_coords_pred <- as.matrix(elvi_coords_pred[which(elvi_ind == 1),])
# plot(coords_pred)
ELVI_pred_data <- elvi_pred_data[which(elvi_ind ==1),]
  
#Making a projection matrix for the predicted values
A_pred.aghy <- inla.spde.make.A(mesh = mesh, loc = aghy_coords_pred)
A_pred.agpe <- inla.spde.make.A(mesh = mesh, loc = agpe_coords_pred)
A_pred.elvi <- inla.spde.make.A(mesh = mesh, loc = elvi_coords_pred)


stack.pred.aghy <- inla.stack(data = list(Endo_status_liberal = NA),
                        A = list(A_pred.aghy,1),
                        effects = list(s = spde.index$spatial.field, 
                                       AGHY_pred_data),
                        tag = 'pred.aghy')
stack.pred.agpe <- inla.stack(data = list(Endo_status_liberal = NA),
                              A = list(A_pred.agpe,1),
                              effects = list(s = spde.index$spatial.field, 
                                             AGPE_pred_data),
                              tag = 'pred.agpe')
stack.pred.elvi <- inla.stack(data = list(Endo_status_liberal = NA),
                              A = list(A_pred.elvi,1),
                              effects = list(s = spde.index$spatial.field, 
                                             ELVI_pred_data),
                              tag = 'pred.elvi')


# this is the data to make a plot for the effect of year at the average lat/lon
# These are the mean locations for each species
mean_loc <- endo_herb %>% 
  group_by(Spp_code, species) %>% 
  summarize(#min_lon = min(lon) - .1*min(lon),
            mean_lon = mean(lon),
           # max_lon = max(lon) + .1*max(lon),
            min_lat = min(lat) + .1*min(lat),
            mean_lat = mean(lat),
            max_lat = max(lat) - .1*max(lat)) %>% 
  pivot_longer(cols = c(-Spp_code, -species, -mean_lon)) %>% 
  rename(lon = mean_lon,lat = value)
  # pivot_longer(cols = c(-Spp_code, -species), names_to = c("level", "coord"), names_sep = "_") %>% 
  # pivot_wider(names_from = coord, values_from = value)

pred.y_data <- data.frame(expand.grid(Intercept = 1, 
                                      # species = c("ELVI"),
                          species = c("AGHY", "AGPE", "ELVI"),
                          year = min(endo_herb$year):max(endo_herb$year), 
                          id = NA,
                          county = NA,
                          collector = NA,
                          scorer = NA)) %>% 
  full_join(mean_loc, by = c("species" = "Spp_code" ))

coords_pred.y <- pred.y_data %>% 
  select(lon, lat) %>% 
  as.matrix()

#Making a projection matrix for the predicted values
A_pred.y <- inla.spde.make.A(mesh = mesh, loc = coords_pred.y)


stack.pred.y <- inla.stack(data = list(Endo_status_liberal = NA),
                         A = list(A_pred.y,1),
                         effects = list(s = spde.index$spatial.field, 
                                        pred.y_data),
                         tag = 'pred.y')



# this is the data from contemporary surveys for model validation
# pred.test_data <- data.frame(Intercept = 1, 
#                                       lon = contemp_surveys$lon,
#                                       lat = contemp_surveys$lat,
#                                       year = contemp_surveys$Year, 
#                                       species = contemp_surveys$SpeciesID,
#                                       id = NA,
#                                       county = NA,
#                                       collector = NA)

pred.test_data <- data.frame(Intercept = 1, 
                             lon = contemp_surveys$lon,
                             lat = contemp_surveys$lat,
                             year = contemp_surveys$Year, 
                             species = contemp_surveys$SpeciesID,
                             id = "new",
                             county = NA,
                             collector = "new",
                             scorer = "new")

coords_pred.test <- pred.test_data %>% 
  select(lon, lat) %>% 
  as.matrix()

#Making a projection matrix for the predicted values
A_pred.test <- inla.spde.make.A(mesh = mesh, loc = coords_pred.test)


stack.pred.test <- inla.stack(data = list(Endo_status_liberal = NA),
                           A = list(A_pred.test,1),
                           effects = list(s = spde.index$spatial.field, 
                                          pred.test_data),
                           tag = 'test')



#
full_stack <- inla.stack(stack.fit, stack.pred.aghy, stack.pred.agpe, stack.pred.elvi, stack.pred.y, stack.pred.test)
# full_stack <- inla.stack(stack.fit, stack.pred.y)

# full_stack <- inla.stack(stack.fit)

# formula1 <- formula(Endo_status_liberal ~ 1 + year + year:lat + year:lon + year:lat:lon
#                     + f(s, model = spde) + f(collector, model = "iid") + f(scorer, model = "iid"))
# formula2 <- formula(Endo_status_liberal ~ 1 + year + year:lat + year:lon + year:lat:lon
#                     + f(s, model = spde) + f(county, model = "iid") + f(collector, model = "iid") + f(scorer, model = "iid"))
# formula3 <- formula(Endo_status_liberal ~ 1 + year + year:lat + year:lon + year:lat:lon
                    # + f(s, model = spde) + f(county,year, model = "iid") + f(collector, model = "iid") + f(scorer, model = "iid"))
# 
# formula <- formula(Endo_status_liberal_1 ~ 0 + year + year:lat + year:lon + year:lat:lon 
#                    + f(s, model = spde))
# 
# formula <- formula(Endo_status_liberal_1 ~ 1 + lon*lat*year+
#                    + f(s, model = spde))
formula1 <- formula(Endo_status_liberal ~ 0 + species + species:year + species:year:lat + species:year:lon + species:year:lat:lon
                    + f(s, model = spde) + f(collector, model = "iid") + f(scorer, model = "iid"))

formula2 <- formula(Endo_status_liberal ~ 0 + species + species:year + species:year:lat + species:year:lon + species:year:lat:lon
                    + f(s, model = spde) + f(county, year, model = "iid", const = TRUE) +  f(collector, model = "iid") + f(scorer, model = "iid"))

formula3 <- formula(Endo_status_liberal ~ 0 + species + species:year + species:year:lat + species:year:lon + species:year:lat:lon
                    + f(s, model = spde) + f(collector, model = "iid") + f(id, model = "iid") + f(scorer, model = "iid"))

formula4 <- formula(Endo_status_liberal ~ 0 + species + species:year +species:year:lat + species:year:lon + species:year:lat:lon + 
                      + f(s, model = spde) + f(collector, model = "iid"))

formula5 <- formula(Endo_status_liberal ~ 0 + species + species:year +species:year:lat + species:year:lon + species:year:lat:lon  + species:year:lat2 + species:year:lon2 + species:year:lat2:lon2 +
                      + f(s, model = spde) )

formula6 <- formula(Endo_status_liberal ~ 0 + species + f(year, model = "iid")
                      + f(s, model = spde) )

formula7 <- formula(Endo_status_liberal ~ 0 + species + species:year+
                    + f(s, model = spde) + f(county, year, model = "iid") +  f(collector, model = "iid") + f(scorer, model = "iid"))
formula8 <- formula(Endo_status_liberal ~ 0 + species*year*lat*lon)

formula9 <- formula(Endo_status_liberal ~ 0 + species*year*lat*lon + f(county, model = "iid"))

formula10 <- formula(Endo_status_liberal ~ 0 + species*year + f(county,year, model = "iid") +f(collector, model = "iid") + f(scorer, model = "iid"))



# formula <- formula(Endo_status_liberal_1 ~ 0 + species + year  + species:year + species:year:lat + species:year:lon + species:year:lat:lon
#                    + f(collector, model = "iid") + f(s, model = spde))

inla.setOption(inla.mode="experimental") 
I1 <- inla(formula = formula1, family = "binomial", Ntrials = 1,
          data = inla.stack.data(full_stack), 
          control.predictor = list(A = inla.stack.A(full_stack),
                                   link=1,compute=TRUE),
          control.family = list(link = "logit"),
          control.compute = list(config = TRUE, dic = TRUE),
          verbose = TRUE)

I2 <- inla(formula = formula2, family = "binomial", Ntrials = 1,
           data = inla.stack.data(full_stack), 
           control.predictor = list(A = inla.stack.A(full_stack),
                                    link=1,compute=TRUE),
           control.family = list(link = "logit"),
           control.compute = list(config = TRUE, dic = TRUE),
           verbose = TRUE)
I3 <- inla(formula = formula3, family = "binomial", Ntrials = 1,
           data = inla.stack.data(full_stack), 
           control.predictor = list(A = inla.stack.A(full_stack),
                                    link=1,compute=TRUE),
           control.family = list(link = "logit"),
           control.compute = list(config = TRUE, dic = TRUE),
           verbose = TRUE)
I4 <- inla(formula = formula4, family = "binomial", Ntrials = 1,
           data = inla.stack.data(full_stack), 
           control.predictor = list(A = inla.stack.A(full_stack),
                                    link=1,compute=TRUE),
           control.family = list(link = "logit"),
           control.compute = list(config = TRUE, dic = TRUE),
           verbose = TRUE)
I6 <- inla(formula = formula6, family = "binomial", Ntrials = 1,
           data = inla.stack.data(full_stack), 
           control.predictor = list(A = inla.stack.A(full_stack),
                                    link=1,compute=TRUE),
           control.family = list(link = "logit"),
           control.compute = list(config = TRUE, dic = TRUE),
           verbose = TRUE)
I10 <- inla(formula = formula10, family = "binomial", Ntrials = 1,
           data = inla.stack.data(full_stack), 
           control.predictor = list(A = inla.stack.A(full_stack),
                                    link=1,compute=TRUE),
           control.family = list(link = "logit"),
           control.compute = list(config = TRUE, dic = TRUE),
           verbose = TRUE)
I1$dic$dic
I2$dic$dic
I3$dic$dic
I4$dic$dic
I7$dic$dic
I8$dic$dic
I9$dic$dic
I10$dic$dic

I1$mode$mode.status # a 0 or low value indicates "convergence"

I <- I1
saveRDS(I, file = "~/Dropbox/Josh&Tom - shared/Endo_Herbarium/Model_Output/inla_spde.rds")
# improve estimates of hyper-parameters by running the following
# I2 <- inla.hyperpar(I2)

I$summary.fixed
I$summary.random


# pulling out the predicted values from the model and attaching them to our pred_data

AGHY_index <- inla.stack.index(full_stack, tag = "pred.aghy")$data
AGHY_pred_data$pred_mean <- I$summary.fitted.values[AGHY_index, "mean"]
AGHY_pred_data$pred_lwr <- I$summary.fitted.values[AGHY_index, "0.025quant"]
AGHY_pred_data$pred_upr <- I$summary.fitted.values[AGHY_index, "0.975quant"]

AGPE_index <- inla.stack.index(full_stack, tag = "pred.agpe")$data
AGPE_pred_data$pred_mean <- I$summary.fitted.values[AGPE_index, "mean"]
AGPE_pred_data$pred_lwr <- I$summary.fitted.values[AGPE_index, "0.025quant"]
AGPE_pred_data$pred_upr <- I$summary.fitted.values[AGPE_index, "0.975quant"]

ELVI_index <- inla.stack.index(full_stack, tag = "pred.elvi")$data
ELVI_pred_data$pred_mean <- I$summary.fitted.values[ELVI_index, "mean"]
ELVI_pred_data$pred_lwr <- I$summary.fitted.values[ELVI_index, "0.025quant"]
ELVI_pred_data$pred_upr <- I$summary.fitted.values[ELVI_index, "0.975quant"]

pred_data <- rbind(filter(AGHY_pred_data,!is.na(Intercept)), filter(AGPE_pred_data, !is.na(Intercept)), filter(ELVI_pred_data, !is.na(Intercept)))

# write_csv(pred_data, file = "pred_data.csv")

pred_data_long <- pred_data %>% 
  pivot_longer(cols = c(pred_mean, pred_lwr, pred_upr), names_to = "variable") %>% 
  mutate(species = case_when(species == "AGHY" ~ "A. hyemalis",
                             species == "AGPE" ~ "A. perennans",
                             species == "ELVI" ~ "E. virginicus"))


ggplot()+
  geom_map(data = outline_map, map = outline_map, aes(long, lat, map_id = region), color = "grey", linewidth = .1, fill = "#FAF9F6")+
  coord_sf(xlim = c(-109,-68), ylim = c(21,49)) +
  geom_tile(data = filter(pred_data_long, variable == "pred_mean" & year == 1920 | variable == "pred_mean" & year == 2020), aes(x = lon, y = lat, fill = value), alpha = .8)+
  facet_wrap(variable~year+species)+
  scale_fill_gradient(
    name = "Endophyte Prevalence",
    low = "blue", high = "orange",
    breaks = c(0,.25,.5,.75,1),
    limits = c(0,1)) +
  theme_light()+
  # theme(strip.background = element_blank())+
  labs(x = "Longitude", y = "Latitude")
  
# making maps for each species for manuscript
AGHY_prevalence <- ggplot()+
  geom_map(data = outline_map, map = outline_map, aes(long, lat, map_id = region), color = "grey", linewidth = .1, fill = "#FAF9F6")+
  coord_sf(xlim = c(-109,-68), ylim = c(21,49)) +
  geom_tile(data = filter(pred_data_long, species == "A. hyemalis" & variable == "pred_mean" & year == 1925 | species == "A. hyemalis" & variable == "pred_mean" & year == 2020), aes(x = lon, y = lat, fill = value), alpha = .8)+
  facet_wrap(~year, ncol = 1)+
  scale_fill_viridis_c(limits = c(0,1))+
  # scale_fill_gradient(
  #   name = "Endophyte Prevalence",
  #   low = endophyte_colors[1], high = endophyte_colors[6],
  #   breaks = c(0,.25,.5,.75,1),
  #   limits = c(0,1)) +
  theme_light()+
  theme(strip.background = element_blank(),
        strip.text = element_text(colour = "black"))+
  labs(title = expression(italic("A. hyemalis")), x = "Longitude", y = "Latitude", fill = "Endophyte Prevalence")
AGHY_prevalence

AGPE_prevalence <- ggplot()+
  geom_map(data = outline_map, map = outline_map, aes(long, lat, map_id = region), color = "grey", linewidth = .1, fill = "#FAF9F6")+
  coord_sf(xlim = c(-109,-68), ylim = c(21,49)) +
  geom_tile(data = filter(pred_data_long, species == "A. perennans" & variable == "pred_mean" & year == 1925 | species == "A. perennans" & variable == "pred_mean" & year == 2020), aes(x = lon, y = lat, fill = value), alpha = .8)+
  facet_wrap(~year, ncol = 1)+
  scale_fill_viridis_c(limits = c(0,1))+
  # scale_fill_gradient(
  #   name = "Endophyte Prevalence",
  #   low = "blue", high = "orange",
  #   breaks = c(0,.25,.5,.75,1),
  #   limits = c(0,1)) +
  theme_light()+
  theme(strip.background = element_blank(),
        strip.text = element_text(colour = "black"))+
  labs(title = expression(italic("A. perennans")), x = "Longitude", y = "", fill = "Endophyte Prevalence")
AGPE_prevalence

ELVI_prevalence <- ggplot()+
  geom_map(data = outline_map, map = outline_map, aes(long, lat, map_id = region), color = "grey", linewidth = .1, fill = "#FAF9F6")+
  coord_sf(xlim = c(-109,-68), ylim = c(21,49)) +
  geom_tile(data = filter(pred_data_long, species == "E. virginicus" & variable == "pred_mean" & year == 1925 | species == "E. virginicus" & variable == "pred_mean" & year == 2020), aes(x = lon, y = lat, fill = value), alpha = .8)+
  facet_wrap(~year, ncol = 1)+
  scale_fill_viridis_c(limits = c(0,1))+
  # scale_fill_gradient(
  #   name = "Endophyte Prevalence",
  #   low = "blue", high = "orange",
  #   breaks = c(0,.25,.5,.75,1),
  #   limits = c(0,1)) +
  theme_light()+
  theme(strip.background = element_blank(),
        strip.text = element_text(colour = "black"))+
  labs(title = expression(italic("E. virginicus")), x = "Longitude", y = "", fill = "Endophyte Prevalence")
ELVI_prevalence


prevalence_map <- AGHY_prevalence + AGPE_prevalence + ELVI_prevalence + plot_layout(nrow = 1, guides = "collect")
prevalence_map
ggsave(prevalence_map, filename = "prevalence_map.png", width = 11, height = 5)

pred_data_change <- pred_data %>% 
  dplyr::select(lon, lat, year, species,pred_mean, pred_lwr, pred_upr) %>% 
  ungroup() %>% 
  pivot_wider(names_from = year, values_from = c(pred_mean, pred_lwr, pred_upr)) %>% 
  mutate(mean_change = pred_mean_2020-pred_mean_1925,
         rel_change = mean_change/pred_mean_1925,
         log_1925 = log(pred_mean_1925),
         log_2020 = log(pred_mean_2020)) %>% 
  mutate(species_name = case_when(species == "AGHY" ~ "A. hyemalis",
                             species == "AGPE" ~ "A. perennans",
                             species == "ELVI" ~ "E. virginicus"))
write_csv(pred_data_change, file = "pred_data_change.csv")

rel_change_plot <- ggplot(data = pred_data_change)+
  geom_map(data = outline_map, map = outline_map, aes( map_id = region), color = "grey", linewidth = .1, fill = "#FAF9F6")+
  coord_sf(xlim = c(-109,-68), ylim = c(21,49)) +
  geom_tile( aes(x = lon, y = lat, fill = rel_change))+
  scale_fill_viridis_c(name = "Relative Change", option = "plasma")+
  facet_wrap(~species_name, ncol = 1)+
  theme_light()+
  theme(strip.background = element_blank(),
        strip.text = element_text(colour = "black", face = "italic"))
rel_change_plot

ggsave(rel_change_plot, filename = "rel_change_plot.png", width = 10, height = 15)

# ggplot(data = prism_diff_pred_df)+
#   geom_tile(aes(x = lon, y = lat, fill = tmean_spring_diff))

pred_change_summary <- pred_data_change %>% 
  group_by(species) %>% 
  summarize(avg_1920 = mean(pred_mean_1920),
            ave_2020 = mean(pred_mean_2020),
            avg_change = mean(mean_change),
            max_1920 = max(pred_mean_1920),
            max_2020 = max(pred_mean_2020),
            max_change = max(mean_change),
            min_1920 = min(pred_mean_1920),
            min_2020 = min(pred_mean_2020),
            min_change = min(mean_change))
######################################################################################################
############ Comparing the contemporary survey prevalence to predicted ############################### 
######################################################################################################
# pulling out the predicted values from the model and attaching them to our pred_data
index <- inla.stack.index(full_stack, tag = "test")$data
pred.test_data$pred_mean <- I$summary.fitted.values[index, "mean"]
pred.test_data$pred_lwr <- I$summary.fitted.values[index, "0.025quant"]
pred.test_data$pred_upr <- I$summary.fitted.values[index, "0.975quant"]

# pred.test_data$obs_endo <- contemp_random_sample$endo_status
pred.test_data$obs_prev <- contemp_surveys$endo_prev
pred.test_data$sample_size <- contemp_surveys$sample_size

# contemp_binned <- contemp_random_sample %>% 
#   mutate(lon_bin = cut(lon, breaks = 25)) %>% 
#   group_by(lon_bin) %>% 
#   summarize(mean_lon = mean(lon),
#             mean_prev = mean(endo_status))


obs_pred_plot <- ggplot(data = pred.test_data)+
  geom_abline(aes(intercept = 0, slope = 1))+
  geom_errorbar(aes(xmin = pred_lwr, xmax = pred_upr, y = obs_prev, col = species), width = .02, alpha = .9)+
  geom_point(aes(x = pred_mean, y = obs_prev, color = species), alpha = .9)+
  lims(x = c(0,1), y = c(0,1))+
  coord_sf()+
  scale_color_manual(values = species_colors)+
  labs(y = "Observed Prevalence", x = "Predicted Prevalence", color = "Species") +
  theme_light()+
  guides(color = "none")

obs_pred_plot
# ggsave(obs_pred_plot, filename = "obs_pred_plot_IDfx.png", width = 4, height = 4)



obs_pred_lon <- ggplot()+
  geom_point(data = pred.test_data, aes(x = lon, y = obs_prev, col = species, size = sample_size), alpha = .6)+
  geom_point(data = pred.test_data, aes(x = lon, y = pred_mean, col = species), size = 4, pch = 3) +
  geom_errorbar(data = pred.test_data, aes(x = lon, ymin = pred_lwr, ymax = pred_upr, col = species), alpha = .9)+
  scale_color_manual(values = species_colors)+
  facet_wrap(~species, ncol = 1)+
  theme_light()+
  theme(strip.background = element_blank(),
        strip.text = element_text(color = "black"))+
  labs(x = "Longitude", y = "Endophyte Prevalence", color = "Species") 
obs_pred_lon

obs_pred_lat <- ggplot()+
  geom_point(data = pred.test_data, aes(x = lat, y = obs_prev, col = species, size = sample_size), alpha = .6)+
  geom_point(data = pred.test_data, aes(x = lat, y = pred_mean, col = species),size = 4, pch = 3) +
  geom_errorbar(data = pred.test_data, aes(x = lat, ymin = pred_lwr, ymax = pred_upr, col = species), alpha = .9)+
  scale_color_manual(values = species_colors)+
  facet_wrap(~species, ncol= 1)+
  theme_light()+
  theme(strip.background = element_blank(),
        strip.text = element_text(color = "black"))+
  labs(x = "Latitude", y = "Endophyte Prevalence", color = "Species") 
obs_pred_lat



contemp_test_plot <- obs_pred_plot + obs_pred_lon + obs_pred_lat + 
  plot_layout(guides = "collect") + plot_annotation(tag_levels = "A")
ggsave(contemp_test_plot, filename = "contemp_test_plot.png", width = 12, height = 8)



ggplot(data = pred.test_data)+
  geom_point(aes(x = lat, y = pred_mean))+
  geom_point(aes(x = lat, y = obs_prev), col = "red")

ggplot(data = pred.test_data)+
  geom_point(aes(x = lon, y = pred_mean))+
  geom_point(aes(x = lon, y = obs_prev), col = "red")

ggplot(data = pred.test_data)+
  geom_point(aes(x = obs_prev, y = pred_mean))

# making a map of the contemp surveys vs historical
test_data_map <- ggplot()+
  geom_map(data = outline_map, map = outline_map, aes(map_id = region), color = "grey", linewidth = .1, fill = "#FAF9F6")+
  geom_map(data = states_shape, map = states_shape, aes(map_id = region), color = "grey", linewidth = .1, fill = NA)+
  geom_point(data = endo_herb, aes(x = lon, y = lat), color = "black",size = 1, pch = 3)+
  geom_point(data = contemp_surveys, aes(x = lon, y = lat), color = "red", size = 1)+
  coord_sf(xlim = c(-109,-68), ylim = c(21,49))+
  theme_light()+
  theme(strip.background = element_blank(),
        strip.text  = element_text(face = "italic", color = "black"))+
  labs(x = "Longitude", y = "Latitude", color = "Endophyte Status")
test_data_map

ggsave(test_data_map, filename = "test_data_map.png", width = 10, height = 8)




################################################################################
############ Plotting the overall trend over time ############################### 
################################################################################
# Binning the data by year for plotting

endo_herb_binned <- endo_herb %>% 
  mutate(binned_year = cut(year, breaks = 20)) %>%
  group_by(Spp_code, species,binned_year) %>%   
  summarise(mean_year = mean(year),
            mean_endo = mean(Endo_status_liberal),
            mean_lon = mean(lon),
            mean_lat = mean(lat),
            sample = n())
  
# mean_loc <- endo_herb_binned %>% 
#   group_by(Spp_code, species) %>% 
#   summarize(lon = mean(mean_lon),
#             lat = mean(mean_lat))

# pulling out the predicting effect across time
index <- inla.stack.index(full_stack, tag = "pred.y")$data
pred.y_data$pred_mean <- I$summary.fitted.values[index, "mean"]
pred.y_data$pred_lwr <- I$summary.fitted.values[index, "0.025quant"]
pred.y_data$pred_upr <- I$summary.fitted.values[index, "0.975quant"]

pred.y_data_neat <- pred.y_data %>% 
  mutate(species = species.y) %>% 
  filter(name != "min_lat") %>% 
  mutate(name  = case_when(name == "mean_lat" ~ paste("35"),
                           name == "max_lat" ~ paste("43"),
                           TRUE ~ name))

year_trend <- ggplot()+
  geom_ribbon(data = pred.y_data_neat, aes(x = year, ymin = pred_lwr, ymax = pred_upr, group = name, linetype = name), color = "grey40", fill = "grey", alpha = .1)+
  geom_line( data = pred.y_data_neat, aes(x = year, y = pred_mean, group = name, linetype = name), color = "black", linewidth = 1, alpha = .8)+
  geom_point(data = endo_herb, aes(x = year, y = Endo_status_liberal), pch = "|", alpha = .7)+
  geom_point(data = endo_herb_binned, aes(x = mean_year, y = mean_endo, col = species , size = sample), alpha = .9)+
  scale_color_manual(values = species_colors, guide = "none")+
  # scale_fill_manual(values = species_colors)+
  scale_linetype_manual(values = c("dotted", "solid"))+
  facet_wrap(~species, nrow = 1)+
  theme_classic()+
  theme(strip.background = element_blank(),
        strip.text = element_blank())+
  labs(y = "Endophyte Prevalence", x = "Year", color = "Species", fill = "Species", size = "Specimen Count", linetype = "Latitude")
year_trend
  
year_hist <- ggplot()+
  geom_histogram(data = filter(endo_herb, species == "A. hyemalis"), aes(x = year), fill = species_colors[1], bins = 120, alpha = .6)+
  geom_histogram(data = filter(endo_herb, species == "E. virginicus"), aes(x = year), fill = species_colors[3], bins = 120, alpha = .6)+
  geom_histogram(data = filter(endo_herb, species == "A. perennans"), aes(x = year),fill = species_colors[2],  bins = 120, alpha = .6)+
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
  plot_layout(ncol = 1, heights = c(1,2)) + plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text( size = 8))

year_plot
ggsave(year_plot, filename = "year_plot.png", width = 7, height = 4)

year_trend_summary <- pred.y_data %>% 
  group_by(species) %>% 
  summarize(minP = min(pred_mean),
            maxP = max(pred_mean)) %>% 
  summarize(min = mean(minP),
            max = mean(maxP))



####################################################################################################################
############ Overall trend across lat and across lon ############################### 
####################################################################################################################
endo_herb_lat_binned <- endo_herb %>% 
  mutate(binned_lat = cut(lat, breaks = 12)) %>%
  rename(species = Spp_code) %>% 
  group_by(species,binned_lat) %>%   
  summarise(mean_lat = mean(lat),
            mean_endo = mean(Endo_status_liberal),
            sample = n())
endo_herb_lon_binned <- endo_herb %>% 
  mutate(binned_lon = cut(year, breaks = 12)) %>%
  rename(species = Spp_code) %>% 
  group_by(species,binned_lon) %>%   
  summarise(mean_lon = mean(lon),
            mean_endo = mean(Endo_status_liberal),
            sample = n())


lat_trend <- ggplot()+
  geom_smooth(data = endo_herb_nice, aes(x = lat, y = Endo_status_liberal, color = species), method = "gam", method.args=list(family="binomial"))+
  geom_point(data = endo_herb_nice, aes(x = lat, y = Endo_status_liberal), pch = "|", alpha = .7)+
  geom_point(data = endo_herb_lat_binned, aes(x = mean_lat, y = mean_endo, col = species, lwd = sample))+
  theme_classic()+
  labs(y = "Endophyte Prevalence", x = "Latitude", lwd = "# of Specimens", color = "Species", fill = "Species")
lat_trend

lon_trend <- ggplot()+
  geom_smooth(data = endo_herb_nice, aes(x = lon, y = Endo_status_liberal, color = species), method = "gam", method.args=list(family="binomial"))+
  geom_point(data = endo_herb_nice, aes(x = lon, y = Endo_status_liberal), pch = "|", alpha = .7)+
  geom_point(data = endo_herb_lon_binned, aes(x = mean_lon, y = mean_endo, col = species, lwd = sample))+
  theme_classic()+
  labs(y = "Endophyte Prevalence", x = "Latitude", lwd = "# of Specimens", color = "Species", fill = "Species")
lon_trend
  
####################################################################################################################
############ Mapping the rate of change ############################### 
####################################################################################################################


pred_change <- pred_data %>% 
  distinct(lon,lat,lon2,lat2,species) 

pred_change$year_slope <- ifelse(pred_change$species == "AGPE", I$summary.fixed["speciesAGPE:year", "mean"], ifelse(pred_change$species == "ELVI", I$summary.fixed["speciesELVI:year", "mean"], ifelse(pred_change$species == "AGHY", I$summary.fixed["speciesAGHY:year", "mean"], NA)))
pred_change$year.lat_slope <- ifelse(pred_change$species == "AGPE", I$summary.fixed["speciesAGPE:year:lat", "mean"], ifelse(pred_change$species == "ELVI", I$summary.fixed["speciesELVI:year:lat", "mean"], ifelse(pred_change$species == "AGHY", I$summary.fixed["speciesAGHY:year:lat", "mean"], NA)))
pred_change$year.lon_slope <- ifelse(pred_change$species == "AGPE", I$summary.fixed["speciesAGPE:year:lon", "mean"], ifelse(pred_change$species == "ELVI", I$summary.fixed["speciesELVI:year:lon", "mean"], ifelse(pred_change$species == "AGHY", I$summary.fixed["speciesAGHY:year:lon", "mean"], NA)))
pred_change$year.lat.lon_slope <- ifelse(pred_change$species == "AGPE", I$summary.fixed["speciesAGPE:year:lat:lon", "mean"], ifelse(pred_change$species == "ELVI", I$summary.fixed["speciesELVI:year:lat:lon", "mean"], ifelse(pred_change$species == "AGHY", I$summary.fixed["speciesAGHY:year:lat:lon", "mean"], NA)))
pred_change$year.lat2_slope <- ifelse(pred_change$species == "AGPE", I$summary.fixed["speciesAGPE:year:lat2", "mean"], ifelse(pred_change$species == "ELVI", I$summary.fixed["speciesELVI:year:lat2", "mean"], ifelse(pred_change$species == "AGHY", I$summary.fixed["speciesAGHY:year:lat2", "mean"], NA)))
pred_change$year.lon2_slope <- ifelse(pred_change$species == "AGPE", I$summary.fixed["speciesAGPE:year:lon2", "mean"], ifelse(pred_change$species == "ELVI", I$summary.fixed["speciesELVI:year:lon2", "mean"], ifelse(pred_change$species == "AGHY", I$summary.fixed["speciesAGHY:year:lon2", "mean"], NA)))
pred_change$year.lat2.lon2_slope <- ifelse(pred_change$species == "AGPE", I$summary.fixed["speciesAGPE:year:lat2:lon2", "mean"], ifelse(pred_change$species == "ELVI", I$summary.fixed["speciesELVI:year:lat2:lon2", "mean"], ifelse(pred_change$species == "AGHY", I$summary.fixed["speciesAGHY:year:lat2:lon2", "mean"], NA)))

pred_change_df <- pred_change %>% 
  mutate(change = year_slope + year.lat_slope*lat + year.lon_slope*lon + year.lat.lon_slope*lat*lon)

  # mutate(change = year_slope + year.lat_slope*lat + year.lon_slope*lon + year.lat.lon_slope*lat*lon)

pred_change_means <- pred_change_df %>% 
  summarize(mean_change =  mean(change),
         sd_change = sd(change))

write_csv(pred_relative_change_df, "pred_change_df.csv")


register_google(key = "")
map <- ggmap::get_map(zoom = 4, maptype = c("satellite"))



ggmap(map)+
  geom_tile(data = pred_change_df, aes(x = lon, y = lat, fill = change), alpha = .3)+
  facet_wrap(~species)+
  scale_fill_gradient(
    name = "Temporal Trend",
    low = "blue", high = "orange")+
  theme_bw()

ggmap(map)+
  geom_tile(data = filter(pred_change_df, species == "AGHY"), aes(x = lon, y = lat, fill = change), alpha = .5)+
  # facet_wrap(~species)+
  scale_fill_gradient(
    name = "Temporal Trend",
    low = "blue", high = "orange")+
  theme_bw()

ggmap(map)+
  geom_tile(data = filter(pred_change_df, species == "ELVI"), aes(x = lon, y = lat, fill = change), alpha = .5)+
  # facet_wrap(~species)+
  scale_fill_gradient(
    name = "Temporal Trend",
    low = "blue", high = "orange")+
  theme_bw()

ggmap(map)+
  geom_tile(data = filter(pred_change_df, species == "AGPE"), aes(x = lon, y = lat, fill = change), alpha = .5)+
  # facet_wrap(~species)+
  scale_fill_gradient(
    name = "Temporal Trend",
    low = "blue", high = "orange")+
  theme_bw()

####################################################################################################################
############ correlating the temporal trend with the change in climate ############################### 
####################################################################################################################

# pred_change_df <- read_csv(file = "pred_change_df.csv") %>% 
#   mutate(lat = as.character(round(lat, digits = 3)), lon = as.character(round(lon, digits = 3)))

# reading in the change in climate normals (1895-1925 and 1990-2020)
prism_diff_pred_df <- read_csv(file = "prism_diff_pred_df.csv") %>% 
  distinct() 
  # mutate(across(contains("ppt_"), ~.x*.1)) # changing the units of the change in precip. to be 10ths of millimeters
prism_for_plotting <- prism_diff_pred_df %>% 
  filter(tmean_annual_cv_diff <25 & tmean_annual_cv_diff>-2)  

tmean_annual_cv_map <- ggplot(prism_for_plotting) +
  geom_map(data = outline_map, map = outline_map, aes(map_id = region), color = "grey", linewidth = .1, fill = "#FAF9F6")+
  geom_map(data = states_shape, map = states_shape, aes(map_id = region), color = "grey", linewidth = .1, fill = NA)+
  geom_tile(aes(x = lon, y = lat, fill = tmean_annual_cv_diff))+
  coord_sf(xlim = c(-109,-68), ylim = c(21,49))+
  theme_light()+
  theme(strip.background = element_blank(),
        strip.text  = element_text(face = "italic", color = "black"))+
  labs(title = "Change in mean annual temperature", x = "Longitude", y = "Latitude", fill = "Change in Coeff. of Var.")
  
tmean_annual_cv_map

tmean_annual_mean_map <- ggplot(prism_for_plotting) +
  geom_map(data = outline_map, map = outline_map, aes(map_id = region), color = "grey", linewidth = .1, fill = "#FAF9F6")+
  geom_map(data = states_shape, map = states_shape, aes(map_id = region), color = "grey", linewidth = .1, fill = NA)+
  geom_tile(aes(x = lon, y = lat, fill = tmean_annual_diff))+
  coord_sf(xlim = c(-109,-68), ylim = c(21,49))+
  theme_light()+
  theme(strip.background = element_blank(),
        strip.text  = element_text(face = "italic", color = "black"))+
  labs(title = "Change in CV of mean annual temperature", x = "Longitude", y = "Latitude", fill = expression("Change in " *degree*"C"))

tmean_annual_mean_map

ppt_annual_mean_map <- ggplot(prism_for_plotting) +
  geom_map(data = outline_map, map = outline_map, aes(map_id = region), color = "grey", linewidth = .1, fill = "#FAF9F6")+
  geom_map(data = states_shape, map = states_shape, aes(map_id = region), color = "grey", linewidth = .1, fill = NA)+
  geom_tile(aes(x = lon, y = lat, fill = ppt_annual_diff))+
  coord_sf(xlim = c(-109,-68), ylim = c(21,49))+
  theme_light()+
  theme(strip.background = element_blank(),
        strip.text  = element_text(face = "italic", color = "black"))+
  labs(title = "Change in total annual precipitation",x = "Longitude", y = "Latitude", fill = expression("Change in mm."))

ppt_annual_mean_map

ppt_annual_cv_map <- ggplot(prism_for_plotting) +
  geom_map(data = outline_map, map = outline_map, aes(map_id = region), color = "grey", linewidth = .1, fill = "#FAF9F6")+
  geom_map(data = states_shape, map = states_shape, aes(map_id = region), color = "grey", linewidth = .1, fill = NA)+
  geom_tile(aes(x = lon, y = lat, fill = ppt_annual_cv_diff))+
  coord_sf(xlim = c(-109,-68), ylim = c(21,49))+
  theme_light()+
  theme(strip.background = element_blank(),
        strip.text  = element_text(face = "italic", color = "black"))+
  labs(title = "Change in CV of total annual precipitation", x = "Longitude", y = "Latitude", fill = expression("Change in CV"))

ppt_annual_cv_map


ppt_spring_mean_map <- ggplot(prism_for_plotting) +
  geom_map(data = outline_map, map = outline_map, aes(map_id = region), color = "grey", linewidth = .1, fill = "#FAF9F6")+
  geom_map(data = states_shape, map = states_shape, aes(map_id = region), color = "grey", linewidth = .1, fill = NA)+
  geom_tile(aes(x = lon, y = lat, fill = ppt_spring_diff))+
  coord_sf(xlim = c(-109,-68), ylim = c(21,49))+
  theme_light()+
  theme(strip.background = element_blank(),
        strip.text  = element_text(face = "italic", color = "black"))+
  labs(title = "Change in total Spring precipitation", x = "Longitude", y = "Latitude", fill = expression("Change in mm."))

ppt_spring_mean_map

ppt_spring_cv_map <- ggplot(prism_for_plotting) +
  geom_map(data = outline_map, map = outline_map, aes(map_id = region), color = "grey", linewidth = .1, fill = "#FAF9F6")+
  geom_map(data = states_shape, map = states_shape, aes(map_id = region), color = "grey", linewidth = .1, fill = NA)+
  geom_tile(aes(x = lon, y = lat, fill = ppt_spring_cv_diff))+
  coord_sf(xlim = c(-109,-68), ylim = c(21,49))+
  theme_light()+
  theme(strip.background = element_blank(),
        strip.text  = element_text(face = "italic", color = "black"))+
  labs(title = "Change in CV of total Spring precipitation", x = "Longitude", y = "Latitude", fill = expression("Change in mm."))

ppt_spring_cv_map


ppt_summer_mean_map <- ggplot(prism_for_plotting) +
  geom_map(data = outline_map, map = outline_map, aes(map_id = region), color = "grey", linewidth = .1, fill = "#FAF9F6")+
  geom_map(data = states_shape, map = states_shape, aes(map_id = region), color = "grey", linewidth = .1, fill = NA)+
  geom_tile(aes(x = lon, y = lat, fill = ppt_summer_diff))+
  coord_sf(xlim = c(-109,-68), ylim = c(21,49))+
  theme_light()+
  theme(strip.background = element_blank(),
        strip.text  = element_text(face = "italic", color = "black"))+
  labs(title = "Change in total summer precipitation", x = "Longitude", y = "Latitude", fill = expression("Change in CV"))

ppt_summer_mean_map

ppt_summer_cv_map <- ggplot(prism_for_plotting) +
  geom_map(data = outline_map, map = outline_map, aes(map_id = region), color = "grey", linewidth = .1, fill = "#FAF9F6")+
  geom_map(data = states_shape, map = states_shape, aes(map_id = region), color = "grey", linewidth = .1, fill = NA)+
  geom_tile(aes(x = lon, y = lat, fill = ppt_summer_cv_diff))+
  coord_sf(xlim = c(-109,-68), ylim = c(21,49))+
  theme_light()+
  theme(strip.background = element_blank(),
        strip.text  = element_text(face = "italic", color = "black"))+
  labs(title = "Change in CV of total summer precipitation", x = "Longitude", y = "Latitude", fill = expression("Change in CV"))

ppt_summer_cv_map

ppt_autumn_mean_map <- ggplot(prism_for_plotting) +
  geom_map(data = outline_map, map = outline_map, aes(map_id = region), color = "grey", linewidth = .1, fill = "#FAF9F6")+
  geom_map(data = states_shape, map = states_shape, aes(map_id = region), color = "grey", linewidth = .1, fill = NA)+
  geom_tile(aes(x = lon, y = lat, fill = ppt_autumn_diff))+
  coord_sf(xlim = c(-109,-68), ylim = c(21,49))+
  theme_light()+
  theme(strip.background = element_blank(),
        strip.text  = element_text(face = "italic", color = "black"))+
  labs(title = "Change in total autumn precipitation",x = "Longitude", y = "Latitude", fill = expression("Change in mm."))

ppt_autumn_mean_map

ppt_autumn_cv_map <- ggplot(prism_for_plotting) +
  geom_map(data = outline_map, map = outline_map, aes(map_id = region), color = "grey", linewidth = .1, fill = "#FAF9F6")+
  geom_map(data = states_shape, map = states_shape, aes(map_id = region), color = "grey", linewidth = .1, fill = NA)+
  geom_tile(aes(x = lon, y = lat, fill = ppt_autumn_cv_diff))+
  coord_sf(xlim = c(-109,-68), ylim = c(21,49))+
  theme_light()+
  theme(strip.background = element_blank(),
        strip.text  = element_text(face = "italic", color = "black"))+
  labs(title = "Change in CV of total autumn precipitation",x = "Longitude", y = "Latitude", fill = expression("Change in CV"))

ppt_autumn_cv_map


ppt_change_maps <- ppt_annual_mean_map + ppt_annual_cv_map + ppt_spring_mean_map + ppt_spring_cv_map+ ppt_summer_mean_map + ppt_summer_cv_map + ppt_autumn_mean_map + ppt_autumn_cv_map + 
  plot_layout(ncol = 2) + plot_annotation(tag_levels = "A")

ggsave(ppt_change_maps, filename = "ppt_change_maps.png", width = 15, height = 20)


tmean_change_maps <- tmean_annual_mean_map + tmean_annual_cv_map +
  plot_layout(ncol = 2) + plot_annotation(tag_levels = "A")
ggsave(tmean_change_maps, filename = "tmean_change_maps.png", width = 15, height = 6)

#testing out pc's as a way to simplify the interpretation
# prism_pc <- prcomp(prism_diff_pred_df[,c(4,5,6,8,9,10)], scale = TRUE)
# summary(prism_pc)
# 
# prism_loadings <- as_tibble(cbind(rownames(prism_pc$rotation), prism_pc$rotation)) %>% 
#   rename(parameter = V1) %>% 
#   pivot_longer(-parameter) %>% 
#   mutate(value = as.numeric(value))
# 
# ggplot(prism_loadings)+
#   geom_bar(aes(x = name, y = value, fill = parameter), position = "dodge", stat = "identity") +
#   theme_light()
# 
# prism_pc_df <- as_tibble(cbind(prism_diff_pred_df$lon, prism_diff_pred_df$lat, prism_pc$x)) %>% 
#   rename(lon = V1, lat = V2)
# 
# pred_change_pc_merge <- pred_data_change %>% 
#   left_join(prism_pc_df) %>% 
#   filter(!is.na(PC1))
# ggplot(data = pred_change_pc_merge) +
#   geom_point(aes(x = lon, y = lat, color = PC3))
# 
# summary(lm(formula = formula(mean_change ~ 1 + PC1 + PC2 + PC3 + PC4), filter(pred_change_pc_merge,species == "AGHY")))

# Merging the climate date with our predicted change in mean prevalence across 1925 and 2020
# pred_change_merge <- pred_ %>% 
pred_change_merge <- pred_data_change %>% 
  left_join(prism_diff_pred_df) %>% 
  filter(tmean_annual_cv_diff <25 & tmean_annual_cv_diff>-2) %>% 
  distinct() %>% 
  na.omit() %>% 
  mutate(lat = as.numeric(lat), lon = as.numeric(lon)) %>% 
  select(-contains("_slope"))



pred_change_merge_long <- pred_change_merge %>% 
  pivot_longer(cols = contains("_diff")) %>% 
  filter(!grepl("sd", name)) %>% 
  filter(!grepl("tmean_spring", name)) %>% 
  filter(!grepl("tmean_summer", name)) %>% 
  filter(!grepl("tmean_autumn", name))

ppt_regression_plot <- ggplot(filter(pred_change_merge_long, grepl("ppt", name)))+
  geom_point(aes(x = value, y = rel_change, color = species), alpha = .2)+
  geom_smooth(aes(x = value, y = rel_change), color = "black", method = "lm")+
  scale_color_manual(values = species_colors) +
  facet_wrap(~name+species, scales = "free", nrow = 3) +
  theme_light()+
  labs(y = "Relative Change in Endophyte Prevalence (%)", x = expression("Change in Precip. (mm"^-1*")"))
ppt_regression_plot

tmean_regression_plot <- ggplot(filter(pred_change_merge_long, grepl("tmean", name)))+
  geom_point(aes(x = value, y = rel_change, color = species), alpha = .2)+
  geom_smooth(aes(x = value, y = rel_change), color = "black", method = "lm")+
  scale_color_manual(values = species_colors) +
  facet_wrap(~name+species, scales = "free", nrow = 3) +
  theme_light()+
  labs(y = "Change in Endophyte Prevalence (%)", x = expression("Change in mean temp. ("*degree*"C)"))
tmean_regression_plot

climate_regression_plot <- ppt_regression_plot + tmean_regression_plot +
  plot_layout(nrow = 1, guides = "collect") 
climate_regression_plot

climate_correlations <- pred_change_merge %>% 
  select(-contains("pred"), -lat, -lon) %>% 
  group_by(species) %>% 
  summarize(rel.ppt_annual_cor = cor(rel_change, ppt_annual_diff),
            rel.ppt_spring_cor = cor(rel_change, ppt_spring_diff),
            rel.ppt_summer_cor = cor(rel_change, ppt_summer_diff),
            rel.ppt_autumn_cor = cor(rel_change, ppt_autumn_diff),
            rel.tmean_annual_cor = cor(rel_change, tmean_annual_diff),
            rel.tmean_spring_cor = cor(rel_change, tmean_spring_diff),
            rel.tmean_summer_cor = cor(rel_change, tmean_summer_diff),
            rel.tmean_autumn_cor = cor(rel_change, tmean_autumn_diff),
            
            rel.ppt_annual_sd_cor = cor(rel_change, ppt_annual_sd_diff),
            rel.ppt_spring_sd_cor = cor(rel_change, ppt_spring_sd_diff),
            rel.ppt_summer_sd_cor = cor(rel_change, ppt_summer_sd_diff),
            rel.ppt_autumn_sd_cor = cor(rel_change, ppt_autumn_sd_diff),
            rel.tmean_annual_sd_cor = cor(rel_change, tmean_annual_sd_diff),
            rel.tmean_spring_sd_cor = cor(rel_change, tmean_spring_sd_diff),
            rel.tmean_summer_sd_cor = cor(rel_change, tmean_summer_sd_diff),
            rel.tmean_autumn_sd_cor = cor(rel_change, tmean_autumn_sd_diff),
            
            rel.ppt_annual_cv_cor = cor(rel_change, ppt_annual_cv_diff),
            rel.ppt_spring_cv_cor = cor(rel_change, ppt_spring_cv_diff),
            rel.ppt_summer_cv_cor = cor(rel_change, ppt_summer_cv_diff),
            rel.ppt_autumn_cv_cor = cor(rel_change, ppt_autumn_cv_diff),
            rel.tmean_annual_cv_cor = cor(rel_change, tmean_annual_cv_diff),
            rel.tmean_spring_cv_cor = cor(rel_change, tmean_spring_cv_diff),
            rel.tmean_summer_cv_cor = cor(rel_change, tmean_summer_cv_diff),
            rel.tmean_autumn_cv_cor = cor(rel_change, tmean_autumn_cv_diff),
            
            mean.ppt_annual_cor = cor(mean_change, ppt_annual_diff),
            mean.ppt_spring_cor = cor(mean_change, ppt_spring_diff),
            mean.ppt_summer_cor = cor(mean_change, ppt_summer_diff),
            mean.ppt_autumn_cor = cor(mean_change, ppt_autumn_diff),
            mean.tmean_annual_cor = cor(mean_change, tmean_annual_diff),
            mean.tmean_spring_cor = cor(mean_change, tmean_spring_diff),
            mean.tmean_summer_cor = cor(mean_change, tmean_summer_diff),
            mean.tmean_autumn_cor = cor(mean_change, tmean_autumn_diff),
            
            mean.ppt_annual_sd_cor = cor(mean_change, ppt_annual_sd_diff),
            mean.ppt_spring_sd_cor = cor(mean_change, ppt_spring_sd_diff),
            mean.ppt_summer_sd_cor = cor(mean_change, ppt_summer_sd_diff),
            mean.ppt_autumn_sd_cor = cor(mean_change, ppt_autumn_sd_diff),
            mean.tmean_annual_sd_cor = cor(mean_change, tmean_annual_sd_diff),
            mean.tmean_spring_sd_cor = cor(mean_change, tmean_spring_sd_diff),
            mean.tmean_summer_sd_cor = cor(mean_change, tmean_summer_sd_diff),
            mean.tmean_autumn_sd_cor = cor(mean_change, tmean_autumn_sd_diff),
            
            mean.ppt_annual_cv_cor = cor(mean_change, ppt_annual_cv_diff),
            mean.ppt_spring_cv_cor = cor(mean_change, ppt_spring_cv_diff),
            mean.ppt_summer_cv_cor = cor(mean_change, ppt_summer_cv_diff),
            mean.ppt_autumn_cv_cor = cor(mean_change, ppt_autumn_cv_diff),
            mean.tmean_annual_cv_cor = cor(mean_change, tmean_annual_cv_diff),
            mean.tmean_spring_cv_cor = cor(mean_change, tmean_spring_cv_diff),
            mean.tmean_summer_cv_cor = cor(mean_change, tmean_summer_cv_diff),
            mean.tmean_autumn_cv_cor = cor(mean_change, tmean_autumn_cv_diff),
            ) %>%  
  pivot_longer(cols = c( -species)) %>% 
  separate_wider_delim(name, ".", names = c("change", "name")) %>% 
  mutate(change = case_when(change == "mean" ~ "Change",
                            change == "rel" ~ "Relative Change")) %>% 
  mutate(parameter_name = factor(case_when(name == "ppt_annual_cor" ~ "Annual Precip.", name == "ppt_spring_cor" ~ "Spring Precip.",name == "ppt_summer_cor" ~ "Summer Precip.",name == "ppt_autumn_cor" ~ "Autumn Precip.",
                                           name == "ppt_annual_sd_cor" ~ "Annual Precip. SD", name == "ppt_spring_sd_cor" ~ "Spring Precip. SD",name == "ppt_summer_sd_cor" ~ "Summer Precip. SD",name == "ppt_autumn_sd_cor" ~ "Autumn Precip. SD",
                                           name == "ppt_annual_cv_cor" ~ "Annual Precip. CV", name == "ppt_spring_cv_cor" ~ "Spring Precip. CV",name == "ppt_summer_cv_cor" ~ "Summer Precip. CV",name == "ppt_autumn_cv_cor" ~ "Autumn Precip. CV",
                                           
                                           name == "tmean_annual_cor" ~ "Avg. Annual Temp.", name == "tmean_spring_cor" ~ "Avg. Spring Temp.",name == "tmean_summer_cor" ~ "Avg. Summer Temp.",name == "tmean_autumn_cor" ~ "Avg. Autumn Temp.",  
                                           name == "tmean_annual_sd_cor" ~ "Avg. Annual Temp. SD", name == "tmean_spring_sd_cor" ~ "Avg. Spring Temp. SD",name == "tmean_summer_sd_cor" ~ "Avg. Summer Temp. SD",name == "tmean_autumn_sd_cor" ~ "Avg. Autumn Temp. SD",  
                                           name == "tmean_annual_cv_cor" ~ "Avg. Annual Temp. CV", name == "tmean_spring_cv_cor" ~ "Avg. Spring Temp. CV",name == "tmean_summer_cv_cor" ~ "Avg. Summer Temp. CV",name == "tmean_autumn_cv_cor" ~ "Avg. Autumn Temp. CV",  
                                           TRUE ~ name),
                                 levels = c("Avg. Annual Temp.", "Avg. Autumn Temp.", "Avg. Summer Temp.", "Avg. Spring Temp.",  "Annual Precip.", "Autumn Precip.", "Summer Precip.", "Spring Precip.", 
                                            "Avg. Annual Temp. SD", "Avg. Autumn Temp. SD", "Avg. Summer Temp. SD", "Avg. Spring Temp. SD",  "Annual Precip. SD", "Autumn Precip. SD", "Summer Precip. SD", "Spring Precip. SD",
                                            "Avg. Annual Temp. CV", "Avg. Autumn Temp. CV", "Avg. Summer Temp. CV", "Avg. Spring Temp. CV",  "Annual Precip. CV", "Autumn Precip. CV", "Summer Precip. CV", "Spring Precip. CV"))) %>% 
  mutate(species = case_when(species == "AGHY" ~ "A. hyemalis",
                             species == "AGPE" ~ "A. perennans",
                             species == "ELVI" ~ "E. virginicus")) %>% 
  mutate(Season = case_when(grepl("annual",name) ~ "Annual",
                            grepl("spring",name) ~ "Spring",
                            grepl("summer",name) ~ "Summer",
                            grepl("autumn",name) ~ "Autumn")) %>% 
  mutate(mean_or_cv = case_when(grepl("cv", name) ~ "CV",
                                !grepl("cv", name) ~ "Mean")) %>% 
  mutate(driver = case_when(grepl("ppt", name) ~ "Precipitation",
                            grepl("tmean", name) ~ "Temperature")) %>% 
  mutate(weak = case_when(value > -.3 & value < .3 ~ "",
                          TRUE ~ ">.3")) %>% 
  filter(change == "Relative Change") %>% 
  filter(grepl("Precip", parameter_name) | grepl("Annual Temp.", parameter_name)) %>% 
  filter(!grepl("SD", parameter_name))


climate_corr_heatmap <- ggplot(data = climate_correlations)+
  geom_tile(aes(y = species, x = as.factor(parameter_name),  fill = value), color = "grey20")+
  geom_point(aes(y = species, x = as.factor(parameter_name), color = weak, alpha = weak), pch = "*", size = 7)+
  scale_color_manual(values = c("grey", "grey15"))+
  scale_alpha_manual(values = c(0,1))+
  scale_fill_distiller(palette = "RdBu", direction = "reverse", type = "div")+
  facet_wrap(factor(Season, levels = c("Annual", "Spring", "Summer", "Autumn")) ~ factor(driver, levels = c("Temperature", "Precipitation")), scales = "free_x", nrow = 1)+
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text.y = element_text(face = "italic")) +
  labs(x = "Observed change in...", y = "", fill = "Correlation Coefficient", color = "", alpha = "")
climate_corr_heatmap
ggsave(climate_corr_heatmap, filename = "climate_corr_heatmap.png", width = 10, height = 7)
# climate_corr_plot <- ggplot(data = climate_correlations)+
#   geom_vline(aes(xintercept = value, color = species)) +
#   facet_wrap(~parameter_name, ncol = 1)
# climate_corr_plot

aghy_pred_change_merge <- pred_change_merge %>% filter(species == "AGHY")
agpe_pred_change_merge <- pred_change_merge %>% filter(species == "AGPE")
elvi_pred_change_merge <- pred_change_merge %>% filter(species == "ELVI")

# bonferonni correction p value would be .005

cor.test(aghy_pred_change_merge$rel_change, aghy_pred_change_merge$ppt_annual_diff)$p.value #*
cor.test(aghy_pred_change_merge$rel_change, aghy_pred_change_merge$ppt_spring_diff)$p.value #*
cor.test(aghy_pred_change_merge$rel_change, aghy_pred_change_merge$ppt_summer_diff)$p.value #*
cor.test(aghy_pred_change_merge$rel_change, aghy_pred_change_merge$ppt_autumn_diff)$p.value
cor.test(aghy_pred_change_merge$rel_change, aghy_pred_change_merge$ppt_annual_cv_diff)$p.value #*
cor.test(aghy_pred_change_merge$rel_change, aghy_pred_change_merge$ppt_spring_cv_diff)$p.value 
cor.test(aghy_pred_change_merge$rel_change, aghy_pred_change_merge$ppt_summer_cv_diff)$p.value # marginal
cor.test(aghy_pred_change_merge$rel_change, aghy_pred_change_merge$ppt_autumn_cv_diff)$p.value #*
cor.test(aghy_pred_change_merge$rel_change, aghy_pred_change_merge$tmean_annual_diff)$p.value #*
cor.test(aghy_pred_change_merge$rel_change, aghy_pred_change_merge$tmean_annual_cv_diff)$p.value



cor.test(agpe_pred_change_merge$rel_change, agpe_pred_change_merge$ppt_annual_diff)$p.value #*
cor.test(agpe_pred_change_merge$rel_change, agpe_pred_change_merge$ppt_spring_diff)$p.value #*
cor.test(agpe_pred_change_merge$rel_change, agpe_pred_change_merge$ppt_summer_diff)$p.value #*
cor.test(agpe_pred_change_merge$rel_change, agpe_pred_change_merge$ppt_autumn_diff)$p.value#*
cor.test(agpe_pred_change_merge$rel_change, agpe_pred_change_merge$ppt_annual_cv_diff)$p.value #*
cor.test(agpe_pred_change_merge$rel_change, agpe_pred_change_merge$ppt_spring_cv_diff)$p.value 
cor.test(agpe_pred_change_merge$rel_change, agpe_pred_change_merge$ppt_summer_cv_diff)$p.value # 
cor.test(agpe_pred_change_merge$rel_change, agpe_pred_change_merge$ppt_autumn_cv_diff)$p.value
cor.test(agpe_pred_change_merge$rel_change, agpe_pred_change_merge$tmean_annual_diff)$p.value #*
cor.test(agpe_pred_change_merge$rel_change, agpe_pred_change_merge$tmean_annual_cv_diff)$p.value#*

cor.test(elvi_pred_change_merge$rel_change, elvi_pred_change_merge$ppt_annual_diff)$p.value #*
cor.test(elvi_pred_change_merge$rel_change, elvi_pred_change_merge$ppt_spring_diff)$p.value #*
cor.test(elvi_pred_change_merge$rel_change, elvi_pred_change_merge$ppt_summer_diff)$p.value #*
cor.test(elvi_pred_change_merge$rel_change, elvi_pred_change_merge$ppt_autumn_diff)$p.value#*
cor.test(elvi_pred_change_merge$rel_change, elvi_pred_change_merge$ppt_annual_cv_diff)$p.value #*
cor.test(elvi_pred_change_merge$rel_change, elvi_pred_change_merge$ppt_spring_cv_diff)$p.value 
cor.test(elvi_pred_change_merge$rel_change, elvi_pred_change_merge$ppt_summer_cv_diff)$p.value #*
cor.test(elvi_pred_change_merge$rel_change, elvi_pred_change_merge$ppt_autumn_cv_diff)$p.value #
cor.test(elvi_pred_change_merge$rel_change, elvi_pred_change_merge$tmean_annual_diff)$p.value #*
cor.test(elvi_pred_change_merge$rel_change, elvi_pred_change_merge$tmean_annual_cv_diff)$p.value#*

plot(filter(pred_change_merge, species == "AGHY")[,25:33])

summary(lm(formula = formula(mean_change ~ 0 + ppt_annual_diff + ppt_spring_diff + ppt_summer_diff + ppt_autumn_diff + ppt_winter_diff + tmean_annual_diff + tmean_spring_diff + tmean_summer_diff + tmean_autumn_diff +  tmean_winter_diff), filter(pred_change_merge,species == "AGPE")))
car::vif(lm(formula = formula(mean_change ~ 1 + ppt_spring_diff + ppt_summer_diff + ppt_autumn_diff + tmean_spring_diff + tmean_summer_diff + tmean_autumn_diff), filter(pred_change_merge,species == "AGPE")))
car::vif(lm(formula = formula(mean_change ~ 1 + ppt_spring_diff * ppt_summer_diff * ppt_autumn_diff *tmean_spring_diff * tmean_summer_diff * tmean_autumn_diff), filter(pred_change_merge,species == "AGPE")), type = "predictor")
cor(filter(pred_change_merge,species == "AGPE"))
summary(lm(formula = formula(mean_change ~  0 + ppt_spring_diff + ppt_summer_diff + ppt_autumn_diff + ppt_winter_diff + tmean_spring_diff + tmean_summer_diff + tmean_autumn_diff + tmean_winter_diff), filter(pred_change_merge,species == "AGPE")))
summary(lm(formula = formula(mean_change ~  0 + ppt_spring_diff*ppt_summer_diff *ppt_autumn_diff* ppt_winter_diff * tmean_spring_diff * tmean_summer_diff *tmean_autumn_diff * tmean_winter_diff), filter(pred_change_merge,species == "AGPE")))

summary(lm(formula = formula(mean_change ~  0 + ppt_annual_diff + tmean_annual_diff), filter(pred_change_merge,species == "AGPE")))
summary(lm(formula = formula(mean_change ~  0 + ppt_annual_diff*tmean_annual_diff), filter(pred_change_merge,species == "AGPE")))

hist(pred_change_merge$tmean_annual_diff)
cor.test(filter(pred_change_merge, species == "AGHY")$mean_change, filter(pred_change_merge, species == "AGHY")$tmean_summer_diff, method = c("pearson"))



# library(devtools)
# install_github('timcdlucas/INLAutils')
# library(INLAutils)
# stack <- inla.stack(data = list(mean_change = filter(pred_change_merge, species == c("ELVI"))$mean_change),
#                     A = list(1), effects = list(data.frame(Intercept = 1, filter(pred_change_merge, species == c("ELVI"))[,c(27,28,29,31,32,33)])))
# 
# I_step <- INLAstep(fam1 = "gaussian", data.frame(filter(pred_change_merge, species == c("ELVI"))),
#                    in_stack = stack,
#                    invariant = "0 + Intercept",
#                    direction = "backwards",
#                    include = c(27,28,29,31,32,33),
#                    inter = 3,
#                    y = "mean_change")
# 
# I_step$best_formula
# 
# 


# climate_formula <- formula(mean_change ~ 0 + species*ppt_spring_diff + species*ppt_summer_diff + species*ppt_autumn_diff + species*ppt_winter_diff + species*tmean_spring_diff + species*tmean_summer_diff + species*tmean_autumn_diff + species*tmean_winter_diff)
# climate_formula <- formula(mean_change ~ 0 + species*ppt_annual_diff + species*ppt_spring_diff + species*ppt_summer_diff + species*ppt_autumn_diff + species*ppt_winter_diff + species*tmean_annual_diff + species*tmean_spring_diff + species*tmean_summer_diff + species*tmean_autumn_diff + species*tmean_winter_diff)
# climate_formula <- formula(mean_change ~ 0 + ppt_spring_diff*ppt_summer_diff*ppt_autumn_diff*ppt_winter_diff*tmean_spring_diff*tmean_summer_diff*tmean_autumn_diff*tmean_winter_diff)
# climate_formula <- formula(mean_change ~ 0 + tmean_spring_diff+tmean_summer_diff+tmean_autumn_diff + ppt_spring_diff+ppt_summer_diff+ppt_autumn_diff)
climate_formula <- formula(mean_change ~ 1 + ppt_spring_diff+ppt_summer_diff+ppt_autumn_diff+tmean_spring_diff+tmean_summer_diff+tmean_autumn_diff)

climate_formula <- formula(mean_change ~ 1 + (tmean_spring_diff+tmean_summer_diff+tmean_autumn_diff + ppt_spring_diff+ppt_summer_diff+ppt_autumn_diff)^3)


climate_formula <- formula(mean_change ~ 1 + tmean_spring_diff*tmean_summer_diff*tmean_autumn_diff*ppt_spring_diff*ppt_summer_diff*ppt_autumn_diff)
climate_formula <- formula(mean_change ~ 1 + tmean_spring_diff)

# climate_formula <- formula(mean_change ~ 0 + ppt_spring_diff+ppt_summer_diff+ppt_autumn_diff+ppt_winter_diff)

# climate_formula <- formula(mean_change ~ 0 + species*ppt_annual_diff*tmean_annual_diff)
# climate_formula <- formula(mean_change ~ 0 + species*ppt_annual_diff + species*tmean_annual_diff)
# 
# climate_formula <- formula(mean_change ~ 0 + species + ppt_spring_diff + ppt_summer_diff + ppt_autumn_diff + ppt_winter_diff + tmean_spring_diff + tmean_summer_diff + tmean_autumn_diff + tmean_winter_diff)
climate_formula <- formula(rel_change ~ 1 + tmean_annual_diff + ppt_spring_diff +ppt_summer_diff + ppt_autumn_diff)

climate_formula <- formula(mean_change ~ 1 + pred_mean_1925 + tmean_annual_diff + ppt_spring_diff +ppt_summer_diff + ppt_autumn_diff)

I_climate <- list()

for(s in 1:3){
I_climate[[s]] <- inla(formula = climate_formula,
           data = filter(pred_change_merge, species == c("AGHY", "AGPE", "ELVI")[s]),
           control.compute = list(config = TRUE, dic = TRUE),
           verbose = FALSE)
}

I_climate[[2]]$summary.fixed
I_climate[[1]]$dic$dic

post_sampling <- list()
post_samps <- list()
for(s in 1:3){
post_sampling[[s]] <- inla.posterior.sample(n = 500, result = I_climate[[s]])
post_samps[[s]] <- inla.posterior.sample.eval(rownames(I_climate[[s]]$summary.fixed), post_sampling[[s]])
rownames(post_samps[[s]]) <- rownames(I_climate[[s]]$summary.fixed)
}


aghy_post_samps <- as_tibble(t(post_samps[[1]])) %>% 
  pivot_longer(cols = everything(), names_to = "parameter") %>% 
  mutate(species = "A. hyemalis")
agpe_post_samps <- as_tibble(t(post_samps[[2]])) %>% 
  pivot_longer(cols = everything(), names_to = "parameter") %>% 
  mutate(species = "A. perennans")
elvi_post_samps <- as_tibble(t(post_samps[[3]])) %>% 
  pivot_longer(cols = everything(), names_to = "parameter") %>% 
  mutate(species = "E. virginicus")

climate_post_samps <- aghy_post_samps %>% 
  bind_rows(agpe_post_samps) %>% 
  bind_rows(elvi_post_samps) %>% 
  mutate(parameter_name = factor(case_when(parameter == "ppt_spring_diff" ~ "Spring Precip.",parameter == "ppt_summer_diff" ~ "Summer Precip.",parameter == "ppt_autumn_diff" ~ "Autumn Precip.",
                                           parameter == "tmean_spring_diff" ~ "Avg. Spring Temp.",parameter == "tmean_summer_diff" ~ "Avg. Summer Temp.",parameter == "tmean_autumn_diff" ~ "Avg. Autumn Temp.",  TRUE ~ parameter),
                                 levels = c("(Intercept)", "Avg. Autumn Temp.", "Avg. Summer Temp.", "Avg. Spring Temp.",  "Autumn Precip.", "Summer Precip.", "Spring Precip.")))



library(ggridges)
tmean_post_plot <- ggplot()+
  geom_vline(aes(xintercept = 0))+
  geom_density_ridges(data = filter(climate_post_samps, grepl("tmean", parameter)), aes(y= parameter, x = value, fill = species), scale = 1, color = NA, alpha = .7)+
  # geom_density_ridges(data = filter(climate_post_samps, grepl("Temp.", parameter_name)), aes(y= parameter_name, x = value, fill = species), scale = 1, color = NA, alpha = .7)+
  scale_fill_manual(values = species_colors)+
  theme_classic()+
  theme(panel.grid.major.y = element_line(color = "darkgrey"),
        legend.text = element_text(face = "italic"))+
  xlab(expression("Rate of Prevalence Change "("%"/degree~C)))+
  ylab("")+ labs(fill = "Species")
tmean_post_plot

ppt_post_plot <- ggplot()+
  geom_vline(aes(xintercept = 0))+
  geom_density_ridges(data = filter(climate_post_samps, grepl("Precip.", parameter_name)), aes(y= parameter_name, x = value, fill = species), scale = 1, color = NA, alpha = .7)+
  scale_fill_manual(values = species_colors)+
  theme_classic()+
  theme(panel.grid.major.y = element_line(color = "darkgrey"),
        legend.text = element_text(face = "italic"))+
  xlab(expression("Rate of Prevalence Change " ("%"/"mm"^-1)))+
  ylab("") + labs(fill = "Species")
ppt_post_plot
  

climate_plot <- tmean_post_plot + ppt_post_plot + plot_layout(ncol = 1, guides = "collect") + plot_annotation(tag_levels = "A")
ggsave(climate_plot, filename = "climate_plot.png", height = 5, width = 6)



####################################################################################################################
############ Making a plot of the change in climate ############################### 
####################################################################################################################
prism_diff_long <- prism_diff_pred_df %>% 
  pivot_longer(c(-lon, -lat))
  
  
ppt_change_plot <- ggplot(filter(prism_diff_long, grepl("ppt", name) & grepl("annual", name)))+
  geom_tile(aes(x = lon, y = lat, fill =  value))+
  facet_wrap(~name)+
  coord_sf()+
  labs(fill = expression("mm"^-1)) +
  theme_light()
ppt_change_plot

tmean_change_plot <- ggplot(filter(prism_diff_long, grepl("tmean", name) & grepl("annual", name) & grepl("cv", name)))+
  geom_tile(aes(x = lon, y = lat, fill =  value))+
  facet_wrap(~name)+
  coord_sf()+
  labs(fill = expression(degree*"C"))+
  theme_light()
tmean_change_plot
  

climate_change_plot <- ppt_change_plot + tmean_change_plot + 
  plot_layout(ncol = 1) + plot_annotation(title = "Change in climate drivers between 1925 and 2020")
  
climate_change_plot
ggsave(climate_change_plot, filename = "climate_change_plot.png", width = 8, height = 4)


####################################################################################################################
############ Plotting the posterior estimates for the scorer and collector random effects ############################### 
####################################################################################################################
# we can sample from the posterior for the different variables
post_samp <- inla.posterior.sample(n = 500, result = I)

scorer_samps <- inla.posterior.sample.eval("scorer", post_samp)
rownames(scorer_samps) <- I$summary.random$scorer[,"ID"]
scorer_samps_df <- as_tibble(t(scorer_samps)) %>% 
  pivot_longer(cols = everything(), names_to = "Scorer")
  
collector_samps <- inla.posterior.sample.eval("collector", post_samp)
rownames(collector_samps) <- I$summary.random$collector[,"ID"]
collector_samps_df <- as_tibble(t(collector_samps)) %>% 
  pivot_longer(cols = everything(), names_to = "Collector")

scorer_summary <- scorer_samps_df %>% 
  group_by(Scorer) %>% 
  summarize(par_mean = mean(value), 
            par_lwr = quantile(value, probs = c(.025)),
            par_upr = quantile(value, probs = c(.975))) %>% 
  filter(Scorer != "new") %>% 
  mutate(Scorer_fct = as.factor(Scorer)) 

levels(scorer_summary$Scorer_fct) <- scorer_no

collector_summary <- collector_samps_df %>% 
  group_by(Collector) %>% 
  summarize(par_mean = mean(value), 
            par_lwr = quantile(value, probs = c(.025)),
            par_upr = quantile(value, probs = c(.975))) %>% 
  filter(Collector != "new") %>% 
  mutate(Collector_fct = as.factor(Collector)) 


levels(collector_summary$Collector_fct) <- collector_no


# hist(scorer_samps_df$value)
# 
# plot(inla.smarginal(I$marginals.random$scorer$index.1))


scorer_plot <- ggplot(scorer_summary)+
  geom_vline(xintercept = 0, lwd = .2)+
  geom_point(aes(x = par_mean,y = Scorer_fct))+
  geom_errorbarh(aes(xmin = par_lwr, xmax = par_upr, height = 0, y = Scorer_fct))+
  theme_light()+
  labs(y = "", x = "Posterior Estimate")
scorer_plot

collector_plot <- ggplot(collector_summary)+
  geom_vline(xintercept = 0, lwd = .2)+
  geom_point(aes(x = par_mean,y = Collector_fct), size = .3)+
  geom_errorbarh(aes(xmin = par_lwr, xmax = par_upr, y = Collector_fct), height = 0, alpha = .6)+
  theme_light()+
  theme(axis.text.y = element_blank())+
  labs(y = "Collector_ ID", x = "Posterior Estimate")

collector_plot

random_fx_plot <- collector_plot +  scorer_plot + 
  plot_layout(nrow = 1) + plot_annotation(tag_levels = "A")
ggsave(random_fx_plot, filename = "random_fx_plot.png", width = 8, height = 6)



# calculating summary overall posterior summary stuff for the manuscript
#
plot(I$marginals.fixed$`speciesAGHY:year`)
plot(I$marginals.fixed$`speciesAGPE:year`)
plot(I$marginals.fixed$`speciesELVI:year`)

# getting all the fixed effectts
fixed_samps <-inla.posterior.sample.eval(rownames(I$summary.fixed), post_samp)
rownames(fixed_samps) <- rownames(I$summary.fixed)
fixed_samps_df <- as_tibble(t(fixed_samps)) %>% 
  pivot_longer(cols = everything(), names_to = "parameter") %>% 
  filter()

# getting just the main year effects
year_samps <-inla.posterior.sample.eval(rownames(I$summary.fixed)[4:6], post_samp)
rownames(year_samps) <- rownames(I$summary.fixed)[4:6]
year_samps_df <- as_tibble(t(year_samps)) %>% 
  pivot_longer(cols = everything(), names_to = "Species") %>% 
  filter()

year_post_summary <- year_samps_df %>% 
  group_by(Species) %>% 
  summarize(mean_trend = mean(value),
            n_draws = n(),
            prop_pos = sum(value>0)/n_draws)




####################################################################################################################
############ Refitting the model with the Conservative sscores and plotting fixed effects ############################### 
####################################################################################################################
I
endo_herb <- endo_herb %>%
  mutate(Endo_status_conservative = case_when(Endo_status_conservative == 4 ~ 1, TRUE  ~ Endo_status_conservative))

stack.cons <- inla.stack(data = list(Endo_status_liberal = endo_herb$Endo_status_conservative),
                                  A = list(A,1),
                                  effects = list(s = spde.index$spatial.field, 
                                                 data.frame(Intercept = rep(1, nrow(endo_herb)),
                                                            year = endo_herb$year,
                                                            year2 = endo_herb$year^2,
                                                            lat = endo_herb$lat,
                                                            lat2 = endo_herb$lat^2,
                                                            lon = endo_herb$lon,
                                                            lon2 = endo_herb$lon^2,
                                                            species = endo_herb$Spp_code,
                                                            id = endo_herb$Sample_id,
                                                            county = endo_herb$County,
                                                            collector = endo_herb$collector_id,
                                                            scorer = endo_herb$scorer_id)),
                                  
                                  tag = 'fit')
I_cons <- inla(formula = formula1, family = "binomial", Ntrials = 1,
           data = inla.stack.data(stack.cons), 
           control.predictor = list(A = inla.stack.A(stack.cons),
                                    link=1,compute=TRUE),
           control.family = list(link = "logit"),
           control.compute = list(config = TRUE, dic = TRUE),
           verbose = TRUE)

I_cons$dic$dic
I_cons$mode$mode.status

post_samp.cons <- inla.posterior.sample(n = 500, result = I_cons)

fixed_samps.cons <-inla.posterior.sample.eval(rownames(I_cons$summary.fixed), post_samp.cons)
rownames(fixed_samps.cons) <- rownames(I_cons$summary.fixed)
fixed_samps.cons_df <- as_tibble(t(fixed_samps.cons)) %>% 
  pivot_longer(cols = everything(), names_to = "parameter") %>% 
  filter()

####################################################################################################################
############ Getting the posterior estimates for the temporal trend at each location ############################### 
####################################################################################################################

# we can sample from the posterior for the different variables
post_samp <- inla.posterior.sample(n = 1000, result = I)

AGHY_post_samps <- inla.posterior.sample.eval(c("speciesAGHY", "speciesAGHY:year", "speciesAGHY:year:lon", "speciesAGHY:year:lat", "speciesAGHY:year:lat:lon"), post_samp)


rownames(AGHY_post_samps) <- c("speciesAGHY", "speciesAGHY:year", "speciesAGHY:year:lon", "speciesAGHY:year:lat", "speciesAGHY:year:lat:lon") 

AGHY_post_df <- as_tibble(t(AGHY_post_samps)) %>% 
  pivot_longer(cols = everything(), names_to = "parameter")

ggplot(data = AGHY_post_df)+
  geom_histogram(aes(x = value), bins = 300)+
  geom_blank(aes(x = -value)) +
  facet_wrap(~parameter, scales = "free", ncol = 1, strip.position = "left")+
  geom_vline(xintercept = 0)+
  theme_classic()+
  theme(strip.background = element_blank(),
        strip.text.y.left = element_text(size = 8, angle = 0),
        axis.line.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank())+
  # axis.text.x = element_text(vjust = 1.5, hjust = 1.5, angle = 45))+
  labs(x = "Posterior Estimate",  y = "") + guides(fill = "none")


####################################################################################################################

# Plotting the overall trend across time
year_effect_df <- as_tibble(list(year=seq(from = min(endo_herb$year), to = max(endo_herb$year))))

year_effect_df$mean <- invlogit(I$summary.fixed["(Intercept)", "mean"] + 
                                  I$summary.fixed["year", "mean"]*year_effect_df$year )
                                  # I$summary.fixed["year:lat", "mean"]*year_effect_df$year*mean_lat + 
                                  # I$summary.fixed["year:lon", "mean"]*year_effect_df$year*mean_lon + 
                                  # I$summary.fixed["year:lat:lon", "mean"]*year_effect_df$year*mean_lat*mean_lon - 
                                  # mean(I$summary.random$s[,"mean"]))
year_effect_df$lwr <- invlogit(I$summary.fixed["(Intercept)", "0.025quant"] + 
                                 I$summary.fixed["year", "0.025quant"]*year_effect_df$year + 
                                 I$summary.fixed["year:lat", "0.025quant"]*year_effect_df$year*mean_lat + 
                                 I$summary.fixed["year:lon", "0.025quant"]*year_effect_df$year*mean_lon + 
                                 I$summary.fixed["year:lat:lon", "0.025quant"]*year_effect_df$year*mean_lat*mean_lon)
year_effect_df$upr <- invlogit(I$summary.fixed["(Intercept)", "0.975quant"] + 
                                 I$summary.fixed["year", "0.975quant"]*year_effect_df$year + 
                                 I$summary.fixed["year:lat", "0.975quant"]*year_effect_df$year*mean_lat + 
                                 I$summary.fixed["year:lon", "0.975quant"]*year_effect_df$year*mean_lon + 
                                 I$summary.fixed["year:lat:lon", "0.975quant"]*year_effect_df$year*mean_lat*mean_lon)

ggplot()+
  geom_point(data = endo_herb_binned, aes(x = mean_year, y = mean_endo, col = mean_lat, lwd = sample))+
  geom_point(data = endo_herb, aes(x = year, y = Endo_status_liberal_1))+
  geom_line( data = year_effect_df, aes(x = year, y = mean))+
  scale_color_viridis_b()+
  # geom_ribbon(aes(x = year, ymin = lwr, ymax = upr))+
  # lims(y = c(0,1))+
  theme_classic()

plot(I$summary.fixed["Intercept", "mean"] + I$summary.fixed["year", "mean"]*seq(min(endo_herb$year):max(endo_herb$year)))

#inla.make.lincomb can be used to make combo of parameters?





# Making the conditional predictions
year_effect <- inla.make.lincomb()

summary(I)
I$summary.fixed
I$summary.random
I$summary.fitted.values

# we can sample from the posterior for the different variables
post_samp <- inla.posterior.sample(n = 10, result = I)
post_samps <- inla.posterior.sample.eval(c("Intercept", "year", "year:lat:lon"), post_samp)
rownames(post_samps) <- c("Intercept", "year", "year:lat:lon")
post_samps_df <- as_tibble(t(post_samps))
hist(post_samps_df$Intercept)
hist(post_samps_df$year)

# pulling out the fitted values from the model
I <- I1
index <- inla.stack.index(full_stack, tag = "fit")$data
fit_mean <- I$summary.fitted.values[index, "mean"]
fit_lwr <- I$summary.fitted.values[index, "0.025quant"]
fit_upr <- I$summary.fitted.values[index, "0.975quant"]

yrep <- yrep_lwr <- yrep_upr <- c()
for(i in 1:length(index)){
yrep[i] <- rbinom(1,  n = 1, prob = fit_mean[i])
yrep_lwr[i] <- rbinom(1,  n = 1, prob = fit_lwr[i])
yrep_upr[i] <- rbinom(1,  n = 1, prob = fit_upr[i])
}

endo_herb$fit_mean <- fit_mean
endo_herb$yrep <- yrep
endo_herb$yrep_lwr <- yrep_lwr
endo_herb$yrep_upr <- yrep_upr

# plot to assess model fit
density_plot <- ggplot(endo_herb)+
  geom_density(aes((yrep)))+
  geom_density(aes(yrep_lwr), linetype = "dashed")+
  geom_density(aes(yrep_upr), linetype = "dashed")+
  geom_density(aes(Endo_status_liberal), color = "red")+
  theme_light()+
  labs(x = "Observed vs Predicted Density")
density_plot
ggsave(density_plot, filename = "density_plot.png",  width = 4, height = 3)

####################################################################################################################
############ Plotting ROC curves and AUC############################### 
####################################################################################################################
roc <- ROCR::prediction(endo_herb$fit_mean, endo_herb$Endo_status_liberal)
performance <- ROCR::performance(roc, "tpr", "fpr")
png("ROC_plot.png", width = 500, height = 500)
par(mfrow=c(1,1))
plot(performance, main = "ROC curve - Observed Data")

auc_roc <- performance(roc, measure = "auc")
auc_roc <- auc_roc@y.values[[1]]

dev.off()

roc_test <- ROCR::prediction(pred.test_data$pred_mean, pred.test_data$obs_endo)
performance_test <-  ROCR::performance(roc_test, "tpr", "fpr")
plot(performance_test, main = "ROC curve - Test Data")

auc_test <- performance(roc_test, measure = "auc")
auc_test <- auc_test@y.values[[1]]


auc_roc
auc_test

#####






fit_data <- rbind(
  data.frame(
    lon = coords[, 1], lat = coords[, 2],
    value = fit_mean, variable = "fit_mean"
  ),
  data.frame(
    lon = coords[, 1], lat = coords[, 2],
    value = fit_lwr, variable = "fit_lwr"
  ),
  data.frame(
    lon = coords[, 1], lat = coords[, 2], 
    value = fit_upr, variable = "fit_upr"
  )
)
fit_data$variable <- as.factor(fit_data$variable)

fit_data <- fit_data %>% 
  mutate(across(c(lat,lon),round, digits = 3))


ggplot(data = fit_data)+
  geom_tile(aes(lon, lat, fill = value), height = 1, width = 1)+
  # geom_point(data = filter(endo_herb, Spp_code == "AGHY"), aes(lon, lat, color = Endo_status_liberal_1), alpha = .6)+
  facet_wrap(~variable, nrow = 1) +
  coord_fixed(ratio = 1)+
  scale_fill_gradient(
    name = "Endophyte Prevalence",
    low = "blue", high = "orange",
    breaks = c(.25,.75),
    limits = c(0,1)) +
  theme_bw()





#####

fitted <- subset(endo_herb, select=c(Spp_code,year,Endo_status_liberal_1))
fitted$Endo_status_liberal_1 <- NA
year_seq <- expand.grid(species = c("AGHY"),year=seq(min(endo_herb$year, na.rm = T), max(endo_herb$year, na.rm = T)), Endo_status_liberal_1=NA)

pred <- invlogit(I$summary.fixed["speciesAGHY",1] + I$summary.fixed["year",1]*year_seq$year)
pred1 <- invlogit(I$summary.fixed["speciesAGPE",1] +  I$summary.fixed["year",1]*year_seq$year + I$summary.fixed["speciesAGPE:year",1]*year_seq$year)
pred2 <- invlogit(I$summary.fixed["speciesELVI",1] + I$summary.fixed["year",1]*year_seq$year + I$summary.fixed["speciesELVI:year",1]*year_seq$year)

plot(year_seq$year, pred,ylim = c(0,1) )
points( year_seq$year, pred)
points(year_seq$year, pred2)

# Trying too visualize the spatial effect
(stepsize <- .3)
nxy <- round(c(diff(range(coords[,1])), diff(range(coords[,2])))/stepsize)
projgrid <- inla.mesh.projector(spde$mesh, xlim=range(coords[,1]),
                                ylim=range(coords[,2]), dims=nxy)
xmean <- inla.mesh.project(projgrid, I$summary.random$s$mean)
xsd <- inla.mesh.project(projgrid, I$summary.random$s$sd)

# setting values outside of the mesh to NA
# require(splancs)
# This is cuttingtoo much out right now
# table(xy.in <- inout(projgrid$lattice$loc,
#                      cbind(coords[,1], coords[,2])))
# xmean[!xy.in] <- xsd[!xy.in] <- NA

library(gridExtra)
library(lattice)
do.call('grid.arrange',
        lapply(list(xmean, xsd),
               levelplot, col.regions=terrain.colors(16),
               xlab='', ylab='', scales=list(draw=FALSE)))



# version fitting time as temporally autocorrelated
rprior <- list(theta = list(prior = "pccor1", param = c(0, 0.9))) # setting a prior for the autocorrelation
n.years <- length(unique(endo_herb$year))
spde.index <- inla.spde.make.index(name = "spatial.field", n.spde = spde$n.spde, n.group =  n.years)

group <- endo_herb$year - min(endo_herb$year) + 1
A <- inla.spde.make.A(mesh = mesh5, loc = coords, group = group)
formula <-  formula(Endo_status_liberal_1 ~ 1 +
                      + f(s, model = spde,
                          group = year, control.group = list(model = "ar1",hyper = rprior)))
