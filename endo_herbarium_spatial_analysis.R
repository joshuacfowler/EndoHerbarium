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
species_colors <- c("#66c2a5", "#fc8d62", "#8da0cb")
endophyte_colors <- c("#fdedd3","#f3c8a8", "#5a727b", "#4986c7", "#181914",  "#163381")


################################################################################
############ Read in the herbarium dataset ############################### 
################################################################################

endo_herb_georef <- read_csv(file = "~/Dropbox/Josh&Tom - shared/Endo_Herbarium/DigitizedHerbariumRecords/endo_herb_georef.csv") %>%
  # filter(Country != "Canada") %>%
  mutate(Spp_code = case_when(grepl("AGHY", Sample_id) ~ "AGHY",
                              grepl("ELVI", Sample_id) ~ "ELVI",
                              grepl("AGPE", Sample_id) ~ "AGPE"))

endo_herb_AGHY <- endo_herb_georef %>% 
  filter(!is.na(Endo_status_liberal)) %>% 
  filter(Spp_code == "AGHY") %>% 
  filter(!is.na(lon) & !is.na(year)) %>% 
  filter(lon>-110)
endo_herb_AGPE <- endo_herb_georef %>% 
  filter(!is.na(Endo_status_liberal)) %>% 
  filter(Spp_code == "AGPE") %>% 
  filter(!is.na(lon) & !is.na(year)) %>% 
  filter(lon>-110)

endo_herb_ELVI <- endo_herb_georef %>% 
  filter(!is.na(Endo_status_liberal)) %>% 
  filter(Spp_code == "ELVI") %>% 
  filter(!is.na(lon) & !is.na(year)) 

endo_herb <- endo_herb_georef %>% 
  filter(!is.na(Endo_status_liberal)) %>%
  filter(!is.na(Spp_code)) %>% 
  filter(!is.na(lon) & !is.na(year)) %>% 
  filter(lon>-110 ) %>% 
  filter(Country != "Canada" & !is.na(County)) %>% 
  mutate(species = case_when(Spp_code == "AGHY" ~ "A. hyemalis",
                              Spp_code == "AGPE" ~ "A. perennans",
                              Spp_code == "ELVI" ~ "E. virginicus")) %>% 
  mutate(year_bin = case_when(year<1970 ~ "pre-1970",
                              year>=1970 ~ "post-1970")) %>% 
  mutate(endo_status_text = case_when(Endo_status_liberal == 0 ~ "E-",
                                      Endo_status_liberal == 1 ~ "E+"))

endo_herb$collector_lastname <- str_replace_all(endo_herb$collector_lastname, "�", " ")
endo_herb$collector_lastname <- str_replace_all(endo_herb$collector_lastname, "xa0", " ")
register_google(key = "")
map <- ggmap::get_map(zoom = 3, source = "google", maptype = c("satellite"))

outline_map <- map_data("world")
states_shape <- map_data("state")

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

# endo_herb <- endo_herb_AGHY
# endo_herb <- endo_herb_AGPE
# endo_herb <- endo_herb_ELVI

# Loading in contemporary survey data for AGHY and ELVI, which we will use for model validation
contemp_surveys <- read_csv(file = "contemp_surveys.csv")
contemp_random_sample <- read_csv(file = "contemp_random_sample.csv")



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
# summary(lm(formula(Endo_status_liberal_1 ~ 0 + Spp_code + Spp_code:year + Spp_code:year:lat + Spp_code:year:lon + Spp_code:year:lat:lon), data = endo_herb))
################################################################################
############ Setting up and running INLA model ############################### 
################################################################################
# messing aruond with INLA model fitting
# generate mesh, which is used to fit the spatial effects
# using a coarser mesh to start out as this is easier computation and sufficient to capture spatial effects
coords <- cbind(endo_herb$lon, endo_herb$lat)

non_convex_bdry <- inla.nonconvex.hull(coords, -0.03, -0.05, resolution = c(100, 100))
sf::sf_use_s2(FALSE)
bdry_st <- st_make_valid(as_tibble(non_convex_bdry$loc)  %>% 
  mutate(lon = V1,  lat = V2) %>% 
  st_as_sf(coords = c("lon", "lat"), crs = 4326) %>% 
  summarise(geometry = st_combine(geometry)) %>% 
  st_cast("POLYGON"))

coastline <- st_make_valid(sf::st_as_sf(maps::map("world", regions = c("usa", "canada", "mexico"), plot = FALSE, fill = TRUE)))
# plot(coastline)

bdry <- st_intersection(coastline$geom, bdry_st)

bdry_polygon <- st_cast(st_sf(bdry), "POLYGON") %>% st_union() %>% 
  as("Spatial")

plot(bdry_polygon)

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

mesh_plot <- ggplot()+
  gg(mesh, exterior = FALSE, edge.color = "grey25", alpha = .1 )+
  geom_point(data = endo_herb, aes(x = lon, y = lat), lwd = .5, color = "red")+
  coord_sf()+
  theme_minimal()
ggsave(mesh_plot, filename = "mesh_plot.png", width = 7, height = 4)


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
                                              collector = endo_herb$collector_lastname,
                                              scorer = endo_herb$scorer_id)),
                    
                    tag = 'fit')
# making a data stack for the predicted value of endophyte prevalence
# need an even  grid of coordinates

# this is the set of data for which we want predictions across the whole space
pred_data <- data.frame(expand.grid(Intercept = 1, 
                        lon = seq(min(endo_herb$lon),max(endo_herb$lon), length.out = 100),
                        lat = seq(min(endo_herb$lat),max(endo_herb$lat), length.out = 100),
                        year = c(1870, 1895, 1920, 1925, 1970, 1990, 2020), 
                        species = c("AGHY", "AGPE", "ELVI"),
                        id = NA,
                        county = NA,
                        collector = NA,
                        scorer = NA)) %>% 
  mutate(year2 = year^2, lon2 = lon^2, lat2 = lat^2)


# keeping just the points that are part of the mesh's boundary
simple_bdry <- st_as_sf(bdry_polygon) %>% ms_simplify()  %>% st_make_valid() %>% as("Spatial") #%>% ms_filter_islands(min_area = 5000000000000)
# plot(simple_bdry)
ind <- point.in.polygon(
  pred_data$lon, pred_data$lat,
  raster::geom(simple_bdry)[,5], raster::geom(simple_bdry)[,6] 
)

coords_pred <- pred_data %>% 
  select(lon, lat) 
coords_pred <- as.matrix(coords_pred[which(ind == 1),])
# plot(coords_pred)
pred_data <- pred_data[which(ind ==1),]
  
#Making a projection matrix for the predicted values
A_pred <- inla.spde.make.A(mesh = mesh, loc = coords_pred)



stack.pred <- inla.stack(data = list(Endo_status_liberal = NA),
                        A = list(A_pred,1),
                        effects = list(s = spde.index$spatial.field, 
                                       pred_data),
                        tag = 'pred')


# this is the data to make a plot for the effect of year at the average lat/lon
pred.y_data <- data.frame(expand.grid(Intercept = 1, 
                                    lon = mean(endo_herb$lon),
                                    lon2 = mean(endo_herb$lon)^2,
                                    lat = mean(endo_herb$lat),
                                    lat2 = mean(endo_herb$lat)^2,
                                    year = min(endo_herb$year):max(endo_herb$year), 
                                    year = (min(endo_herb$year):max(endo_herb$year))^2,
                                    species = c("AGHY", "AGPE", "ELVI"),
                                    id = NA,
                                    county = NA,
                                    collector = NA,
                                    scorer = NA))

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
                             lon = contemp_random_sample$lon,
                             lat = contemp_random_sample$lat,
                             year = contemp_random_sample$Year, 
                             species = contemp_random_sample$SpeciesID,
                             id = NA,
                             county = NA,
                             collector = NA,
                             scorer = NA)

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
full_stack <- inla.stack(stack.fit, stack.pred, stack.pred.y, stack.pred.test)
# full_stack <- inla.stack(stack.fit, stack.pred.y)

# full_stack <- inla.stack(stack.fit)

# formula <- formula(Endo_status_liberal ~ 1 + year + year:lat + year:lon + year:lat:lon
                    # + f(s, model = spde))
# 
# formula <- formula(Endo_status_liberal_1 ~ 0 + year + year:lat + year:lon + year:lat:lon 
#                    + f(s, model = spde))
# 
# formula <- formula(Endo_status_liberal_1 ~ 1 + lon*lat*year+
#                    + f(s, model = spde))
formula1 <- formula(Endo_status_liberal ~ 0 + species + species:year + species:year:lat + species:year:lon + species:year:lat:lon
                    + f(s, model = spde) + f(collector, model = "iid") + f(scorer, model = "iid"))

formula2 <- formula(Endo_status_liberal ~ 0 + species + species:year + species:year:lat + species:year:lon + species:year:lat:lon
                    + f(s, model = spde) + f(collector, model = "iid") + f(county, model = "iid") + f(id, model = "iid") + f(scorer, model = "iid"))

formula3 <- formula(Endo_status_liberal ~ 0 + species + species:year + species:year:lat + species:year:lon + species:year:lat:lon
                    + f(s, model = spde) + f(collector, model = "iid") + f(id, model = "iid") + f(scorer, model = "iid"))

formula4 <- formula(Endo_status_liberal ~ 0 + species + species:year +species:year:lat + species:year:lon + species:year:lat:lon + 
                      + f(s, model = spde) + f(collector, model = "iid"))

formula5 <- formula(Endo_status_liberal ~ 0 + species + species:year +species:year:lat + species:year:lon + species:year:lat:lon  + species:year:lat2 + species:year:lon2 + species:year:lat2:lon2 +
                      + f(s, model = spde) )

                    

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
I1$dic$dic
I2$dic$dic
I3$dic$dic
I4$dic$dic

I <- I1
saveRDS(I, file = "~/Dropbox/Josh&Tom - shared/Endo_Herbarium/Model_Output/inla_spde.rds")
# improve estimates of hyper-parameters by running the following
# I2 <- inla.hyperpar(I2)

I$summary.fixed
I$summary.random


# pulling out the predicted values from the model and attaching them to our pred_data
index <- inla.stack.index(full_stack, tag = "pred")$data
pred_data$pred_mean <- I$summary.fitted.values[index, "mean"]
pred_data$pred_lwr <- I$summary.fitted.values[index, "0.025quant"]
pred_data$pred_upr <- I$summary.fitted.values[index, "0.975quant"]

write_csv(pred_data, file = "pred_data.csv")

pred_data_long <- pred_data %>% 
  pivot_longer(cols = c(pred_mean, pred_lwr, pred_upr), names_to = "variable") %>% 
  mutate(species = case_when(species == "AGHY" ~ "A. hyemalis",
                             species == "AGPE" ~ "A. perennans",
                             species == "ELVI" ~ "E. virginicus"))

# endo_herb_neat <- endo_herb %>% 
#   mutate(binned_year= case_when(year<=1920))) %>% 
#   rename(species = Spp_code) %>% 
#   group_by(species, binned_year) %>% 
#     summarize(mean = mean(year))
#   mutate(year = case_when(binned_year == "(1.82e+03,1.92e+03]" ~ 1920,
#                           binned_year == "(1.97e+03,2.02e+03]" ~ 1970,
#                           binned_year == "(1.82e+03,1.92e+03]" ~ 2020))
  
ggplot()+
  geom_map(data = outline_map, map = outline_map, aes(long, lat, map_id = region), color = "grey", linewidth = .1, fill = "#FAF9F6")+
  coord_sf(xlim = c(-109,-68), ylim = c(21,49)) +
  geom_tile(data = filter(pred_data_long, variable == "pred_mean" & year == 1920 | variable == "pred_mean" & year == 2020), aes(x = lon, y = lat, fill = value), alpha = .8)+
  facet_wrap(variable~year+species, strip.position = "left")+
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
  geom_tile(data = filter(pred_data_long, species == "A. hyemalis" & variable == "pred_mean" & year == 1920 | species == "A. hyemalis" & variable == "pred_mean" & year == 2020), aes(x = lon, y = lat, fill = value), alpha = .8)+
  facet_wrap(~year, ncol = 1)+
  scale_fill_gradient(
    name = "Endophyte Prevalence",
    low = "blue", high = "orange",
    breaks = c(0,.25,.5,.75,1),
    limits = c(0,1)) +
  theme_light()+
  theme(strip.background = element_blank(),
        strip.text = element_text(colour = "black"))+
  labs(title = expression(italic("A. hyemalis")), x = "Longitude", y = "Latitude")
AGHY_prevalence

AGPE_prevalence <- ggplot()+
  geom_map(data = outline_map, map = outline_map, aes(long, lat, map_id = region), color = "grey", linewidth = .1, fill = "#FAF9F6")+
  coord_sf(xlim = c(-109,-68), ylim = c(21,49)) +
  geom_tile(data = filter(pred_data_long, species == "A. perennans" & variable == "pred_mean" & year == 1920 | species == "A. perennans" & variable == "pred_mean" & year == 2020), aes(x = lon, y = lat, fill = value), alpha = .8)+
  facet_wrap(~year, ncol = 1)+
  scale_fill_gradient(
    name = "Endophyte Prevalence",
    low = "blue", high = "orange",
    breaks = c(0,.25,.5,.75,1),
    limits = c(0,1)) +
  theme_light()+
  theme(strip.background = element_blank(),
        strip.text = element_text(colour = "black"))+
  labs(title = expression(italic("A. perennans")), x = "Longitude", y = "")
AGPE_prevalence

ELVI_prevalence <- ggplot()+
  geom_map(data = outline_map, map = outline_map, aes(long, lat, map_id = region), color = "grey", linewidth = .1, fill = "#FAF9F6")+
  coord_sf(xlim = c(-109,-68), ylim = c(21,49)) +
  geom_tile(data = filter(pred_data_long, species == "E. virginicus" & variable == "pred_mean" & year == 1920 | species == "E. virginicus" & variable == "pred_mean" & year == 2020), aes(x = lon, y = lat, fill = value), alpha = .8)+
  facet_wrap(~year, ncol = 1)+
  scale_fill_gradient(
    name = "Endophyte Prevalence",
    low = "blue", high = "orange",
    breaks = c(0,.25,.5,.75,1),
    limits = c(0,1)) +
  theme_light()+
  theme(strip.background = element_blank(),
        strip.text = element_text(colour = "black"))+
  labs(title = expression(italic("E. virginicus")), x = "Longitude", y = "")
ELVI_prevalence


prevalence_map <- AGHY_prevalence + AGPE_prevalence + ELVI_prevalence + plot_layout(nrow = 1, guides = "collect")
prevalence_map
ggsave(prevalence_map, filename = "prevalence_map.png", width = 12, height = 5)

pred_data_change <- pred_data %>% 
  dplyr::select(lon, lat, year, species,pred_mean, pred_lwr, pred_upr) %>% 
  ungroup() %>% 
  pivot_wider(names_from = year, values_from = c(pred_mean, pred_lwr, pred_upr)) %>% 
  mutate(mean_change = pred_mean_2020-pred_mean_1925)
write_csv(pred_data_change, file = "pred_data_change.csv")


ggplot(data = filter(pred_data_change, species == "ELVI"))+
  geom_tile(aes(x = lon, y = lat, fill = mean_change))+
  facet_wrap(~species)

ggplot(data = prism_diff_pred_df)+
  geom_tile(aes(x = lon, y = lat, fill = tmean_spring_diff))
######################################################################################################
############ Comparing the contemporary survey prevalence to predicted ############################### 
######################################################################################################
# pulling out the predicted values from the model and attaching them to our pred_data
index <- inla.stack.index(full_stack, tag = "test")$data
pred.test_data$pred_mean <- I$summary.fitted.values[index, "mean"]
pred.test_data$pred_lwr <- I$summary.fitted.values[index, "0.025quant"]
pred.test_data$pred_upr <- I$summary.fitted.values[index, "0.975quant"]


pred.test_data$obs_prev <- contemp_random_sample$endo_status

contemp_binned <- contemp_random_sample %>% 
  mutate(lon_bin = cut(lon, breaks = 25)) %>% 
  group_by(lon_bin) %>% 
  summarize(mean_lon = mean(lon),
            mean_prev = mean(endo_status))

ggplot(data = contemp_binned)+
  geom_point(data = contemp_random_sample, aes(x = lon, y = endo_status), col = "purple")+
  geom_point(aes(x = mean_lon, y = mean_prev), col = "red", lwd = 2)+
  geom_point(data = pred.test_data, aes(x = lon, y = pred_mean)) +
  geom_errorbar(data = pred.test_data, aes(x = lon, ymin = pred_lwr, ymax = pred_upr))+
  geom_point(data = filter(contemp_surveys, SpeciesID == "AGHY"), aes(x = lon, y = endo_prev), col = "blue")
  

ggplot(data = pred.test_data)+
  geom_point(aes(x = lat, y = lon, color = species))

ggplot(data = pred.test_data)+
  geom_point(aes(x = pred_mean, y = obs_prev, color = species, shape = as.factor(year)))+
  geom_errorbar(aes(xmin = pred_lwr, xmax = pred_upr, y = obs_prev, col = species))+
  # geom_abline(aes(intercept = 0, slope = 1))+
  lims(x = c(0,1), y = c(0,1))



ggplot(data = pred.test_data)+
  geom_point(aes(x = lat, y = pred_mean))+
  geom_point(aes(x = lat, y = obs_prev), col = "red")

ggplot(data = pred.test_data)+
  geom_point(aes(x = lon, y = pred_mean))+
  geom_point(aes(x = lon, y = obs_prev), col = "red")

################################################################################
############ Plotting the overall trend over time ############################### 
################################################################################
# Binning the data by year for plotting

endo_herb_binned <- endo_herb %>% 
  mutate(binned_year = cut(year, breaks = 20)) %>%
  group_by(species,binned_year) %>%   
  summarise(mean_year = mean(year),
            mean_endo = mean(Endo_status_liberal),
            sample = n())
  


# pulling out the predicting effect across time
index <- inla.stack.index(full_stack, tag = "pred.y")$data
pred.y_data$pred_mean <- I$summary.fitted.values[index, "mean"]
pred.y_data$pred_lwr <- I$summary.fitted.values[index, "0.025quant"]
pred.y_data$pred_upr <- I$summary.fitted.values[index, "0.975quant"]

pred.y_data_neat <- pred.y_data %>% 
  mutate(species = case_when(species == "AGHY" ~ "A. hyemalis",
                             species == "AGPE" ~ "A. perennans",
                             species == "ELVI" ~ "E. virginicus"))

year_trend <- ggplot()+
  geom_ribbon(data = pred.y_data_neat, aes(x = year, ymin = pred_lwr, ymax = pred_upr, group = species, fill = species), alpha = .2)+
  geom_line( data = pred.y_data_neat, aes(x = year, y = pred_mean, group = species, col = species), linewidth = 1)+
  geom_point(data = endo_herb, aes(x = year, y = Endo_status_liberal), pch = "|", alpha = .7)+
  geom_point(data = endo_herb_binned, aes(x = mean_year, y = mean_endo, col = species, lwd = sample))+
  scale_color_manual(values = species_colors)+
  scale_fill_manual(values = species_colors)+
  theme_classic()+
  labs(y = "Endophyte Prevalence", x = "Year", lwd = "# of Specimens", color = "Species", fill = "Species")

year_trend
  
year_hist <- ggplot()+
  geom_histogram(data = filter(endo_herb, species == "A. hyemalis"), aes(x = year), fill = species_colors[1], bins = 120, alpha = .6)+
  geom_histogram(data = filter(endo_herb, species == "E. virginicus"), aes(x = year), fill = species_colors[3], bins = 120, alpha = .6)+
  geom_histogram(data = filter(endo_herb, species == "A. perennans"), aes(x = year),fill = species_colors[2],  bins = 120, alpha = .6)+
  theme_classic()+
  theme(axis.line = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank())+
  labs(y = "# of Specimens")+
  guides(fill = "none")
year_hist


year_plot <- year_hist + year_trend + 
  plot_layout(ncol = 1, heights = c(1,2), guides = "collect") + plot_annotation(tag_levels = "A")

year_plot
ggsave(year_plot, filename = "year_plot.png", width = 7, height = 5)

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
  mutate(across(contains("ppt_"), ~.x*.1)) # changing the units of the change in precip. to be 10ths of millimeters


# Merging the climate date with our predicted change in mean prevalence across 1925 and 2020
# pred_change_merge <- pred_ %>% 
pred_change_merge <- pred_data_change %>% 
  left_join(prism_diff_pred_df) %>% 
  distinct() %>% 
  na.omit() %>% 
  mutate(lat = as.numeric(lat), lon = as.numeric(lon)) %>% 
  select(-contains("_slope"))

ggplot(pred_change_merge) +
  geom_point(aes(x = lon, y = lat))

pred_change_merge_long <- pred_change_merge %>% 
  pivot_longer(cols = contains("_diff"))

ggplot(pred_change_merge_long)+
  geom_point(aes(x = value, y = mean_change, color = species))+
  geom_smooth(aes(x = value, y = mean_change, fill = species), method = "glm")+
  facet_wrap(~species+name, scales = "free")

summary(lm(formula = formula(mean_change ~ 0 + ppt_annual_diff + ppt_spring_diff + ppt_summer_diff + ppt_autumn_diff + ppt_winter_diff + tmean_annual_diff + tmean_spring_diff + tmean_summer_diff + tmean_autumn_diff +  tmean_winter_diff), filter(pred_change_merge,species == "AGPE")))
summary(lm(formula = formula(mean_change ~  0 + ppt_spring_diff + ppt_summer_diff + ppt_autumn_diff + ppt_winter_diff + tmean_spring_diff + tmean_summer_diff + tmean_autumn_diff + tmean_winter_diff), filter(pred_change_merge,species == "AGPE")))
summary(lm(formula = formula(mean_change ~  0 + ppt_spring_diff*ppt_summer_diff *ppt_autumn_diff* ppt_winter_diff * tmean_spring_diff * tmean_summer_diff *tmean_autumn_diff * tmean_winter_diff), filter(pred_change_merge,species == "AGPE")))

summary(lm(formula = formula(mean_change ~  0 + ppt_annual_diff + tmean_annual_diff), filter(pred_change_merge,species == "AGPE")))
summary(lm(formula = formula(mean_change ~  0 + ppt_annual_diff*tmean_annual_diff), filter(pred_change_merge,species == "AGPE")))

hist(pred_change_merge$tmean_annual_diff)
cor.test(filter(pred_change_merge, species == "AGPE")$relative_change, filter(pred_change_merge, species == "AGPE")$tmean_annual_diff, method = c("pearson"))

# climate_formula <- formula(mean_change ~ 0 + species*ppt_spring_diff + species*ppt_summer_diff + species*ppt_autumn_diff + species*ppt_winter_diff + species*tmean_spring_diff + species*tmean_summer_diff + species*tmean_autumn_diff + species*tmean_winter_diff)
# climate_formula <- formula(mean_change ~ 0 + species*ppt_annual_diff + species*ppt_spring_diff + species*ppt_summer_diff + species*ppt_autumn_diff + species*ppt_winter_diff + species*tmean_annual_diff + species*tmean_spring_diff + species*tmean_summer_diff + species*tmean_autumn_diff + species*tmean_winter_diff)
# climate_formula <- formula(mean_change ~ 0 + ppt_spring_diff*ppt_summer_diff*ppt_autumn_diff*ppt_winter_diff*tmean_spring_diff*tmean_summer_diff*tmean_autumn_diff*tmean_winter_diff)
climate_formula <- formula(mean_change ~ 1 + ppt_spring_diff+ppt_summer_diff+ppt_autumn_diff+tmean_spring_diff+tmean_summer_diff+tmean_autumn_diff)
# climate_formula <- formula(mean_change ~ 0 + ppt_spring_diff+ppt_summer_diff+ppt_autumn_diff+ppt_winter_diff)

# climate_formula <- formula(mean_change ~ 0 + species*ppt_annual_diff*tmean_annual_diff)
# climate_formula <- formula(mean_change ~ 0 + species*ppt_annual_diff + species*tmean_annual_diff)
# 
# climate_formula <- formula(mean_change ~ 0 + species + ppt_spring_diff + ppt_summer_diff + ppt_autumn_diff + ppt_winter_diff + tmean_spring_diff + tmean_summer_diff + tmean_autumn_diff + tmean_winter_diff)

# making a data stack for the predicted value of endophyte 
I_climate <- list()
for(s in 1:3){
I_climate[[s]] <- inla(formula = climate_formula,
           data = filter(pred_change_merge, species == c("AGHY", "AGPE", "ELVI")[s]),
           control.compute = list(config = TRUE, dic = TRUE),
           verbose = FALSE)
}

I_climate[[1]]$summary.fixed
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



aghy_post_samps <- as_tibble(t(post_samps[[1]])) %>% 
  pivot_longer(cols = everything(), names_to = "parameter") %>% 
  mutate(parameter_name = factor(case_when(parameter == "ppt_spring_diff" ~ "Spring Precip.",parameter == "ppt_summer_diff" ~ "Summer Precip.",parameter == "ppt_autumn_diff" ~ "Autumn Precip.",
                               parameter == "tmean_spring_diff" ~ "Avg. Spring Temp.",parameter == "tmean_summer_diff" ~ "Avg. Summer Temp.",parameter == "tmean_autumn_diff" ~ "Avg. Autumn Temp.",  TRUE ~ parameter),
                               levels = c("(Intercept)", "Avg. Autumn Temp.", "Avg. Summer Temp.", "Avg. Spring Temp.",  "Autumn Precip.", "Summer Precip.", "Spring Precip.")))
agpe_post_samps <- as_tibble(t(post_samps[[2]])) %>% 
  pivot_longer(cols = everything(), names_to = "parameter") %>% 
  mutate(parameter_name = factor(case_when(parameter == "ppt_spring_diff" ~ "Spring Precip.",parameter == "ppt_summer_diff" ~ "Summer Precip.",parameter == "ppt_autumn_diff" ~ "Autumn Precip.",
                               parameter == "tmean_spring_diff" ~ "Avg. Spring Temp.",parameter == "tmean_summer_diff" ~ "Avg. Summer Temp.",parameter == "tmean_autumn_diff" ~ "Avg. Autumn Temp.",TRUE ~ parameter),
                               levels = c("(Intercept)", "Avg. Autumn Temp.", "Avg. Summer Temp.", "Avg. Spring Temp.",  "Autumn Precip.", "Summer Precip.", "Spring Precip.")))
elvi_post_samps <- as_tibble(t(post_samps[[3]])) %>% 
  pivot_longer(cols = everything(), names_to = "parameter") %>% 
  mutate(parameter_name = factor(case_when(parameter == "ppt_spring_diff" ~ "Spring Precip.",parameter == "ppt_summer_diff" ~ "Summer Precip.",parameter == "ppt_autumn_diff" ~ "Autumn Precip.",
                               parameter == "tmean_spring_diff" ~ "Avg. Spring Temp.",parameter == "tmean_summer_diff" ~ "Avg. Summer Temp.",parameter == "tmean_autumn_diff" ~ "Avg. Autumn Temp.", TRUE ~ parameter),
                               levels = c("(Intercept)", "Avg. Autumn Temp.", "Avg. Summer Temp.", "Avg. Spring Temp.",  "Autumn Precip.", "Summer Precip.", "Spring Precip.")))

library(ggridges)
tmean_post_plot <- ggplot()+
  geom_vline(aes(xintercept = 0))+
  geom_density_ridges(data = filter(climate_post_samps, grepl("Temp.", parameter_name)), aes(y= parameter_name, x = value, fill = species), scale = 1, color = NA, alpha = .7)+
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
png("ROC_plot.png", width = 500, height = 1000)
par(mfrow=c(2,1))
plot(performance, main = "ROC curve - Observed Data")

auc_roc <- performance(roc, measure = "auc")
auc_roc <- auc_roc@y.values[[1]]


roc_test <- ROCR::prediction(pred.test_data$pred_mean, pred.test_data$obs_prev)
performance_test <-  ROCR::performance(roc_test, "tpr", "fpr")
plot(performance_test, main = "ROC curve - Test Data")

auc_test <- performance(roc_test, measure = "auc")
auc_test <- auc_test@y.values[[1]]

dev.off()



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
