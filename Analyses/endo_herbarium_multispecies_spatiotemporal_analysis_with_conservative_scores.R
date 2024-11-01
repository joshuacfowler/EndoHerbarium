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
  filter(!is.na(Endo_status_conservative)) %>%
  filter(!is.na(spp_code)) %>% 
  filter(!is.na(lon) & !is.na(year)) %>% 
  filter(lon>-110 ) %>% 
  filter(Country != "Canada" ) %>% 
  mutate(year_bin = case_when(year<1970 ~ "pre-1970",
                              year>=1970 ~ "post-1970")) %>% 
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
  filter(!is.na(Endo_status_conservative)) %>% 
  filter(spp_code == "AGHY") %>% 
  filter(!is.na(lon) & !is.na(year))

endo_herb_AGPE <- endo_herb %>% 
  filter(!is.na(Endo_status_conservative)) %>% 
  filter(spp_code == "AGPE") %>% 
  filter(!is.na(lon) & !is.na(year)) %>% 
  filter(lat>=20) # dropping one plant in mexico


endo_herb_ELVI <- endo_herb %>% 
  filter(!is.na(Endo_status_conservative)) %>% 
  filter(spp_code == "ELVI") %>% 
  filter(!is.na(lon) & !is.na(year)) 

endo_herb <- endo_herb %>% 
  mutate(Endo_status_conservative = case_when(Endo_status_conservative == 4 ~ 1, 
                                              TRUE ~ Endo_status_conservative)) %>% 
  filter(!is.na(Endo_status_conservative)) %>% 
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

cmp <- ~ space(geometry, model = spde) + 
  space.species1(geometry, model = spde) + space.species2(geometry, model = spde) + space.species3(geometry, model = spde) +
  time.species1(geometry, weights = std_year, model = spde) + time.species2(geometry, weights = std_year, model = spde) + time.species3(geometry, weights = std_year, model = spde) +
  + int.species1(1) + int.species2(1) + int.species3(1)+
  + year.species1(main = ~0 + std_year, model = "fixed") + year.species2(main = ~0 + std_year, model = "fixed") + year.species3(main = ~0 + std_year, model = "fixed")+
  # herbarium(herbarium_index, model = "iid", constr = TRUE, mapper = bru_mapper_index(max(data$herbarium_index)), hyper = list(pc_prec))+
  scorer(scorer_index, model = "iid", constr = TRUE, mapper = bru_mapper_index(max(data$scorer_index)), hyper = list(pc_prec)) +
  collector(collector_index, model = "iid", constr = TRUE, mapper = bru_mapper_index(max(data$collector_index, na.rm = T)), hyper = list(pc_prec))

fml.aghy <- Endo_status_conservative ~ 0 + int.species1 + year.species1 + space + space.species1 + time.species1 + scorer + collector 
fml.agpe <- Endo_status_conservative ~ 0 + int.species2 + year.species2 + space + space.species2 + time.species2 + scorer + collector 
fml.elvi <- Endo_status_conservative ~ 0 + int.species3 + year.species3 + space + space.species3 + time.species3 + scorer + collector


data_aghy <- data[data$species_index == 1,]
data_agpe <- data[data$species_index == 2,]
data_elvi <- data[data$species_index == 3,]

lik_aghy <- like(formula = fml.aghy,
                 family = "binomial",
                 Ntrials = 1,
                 data = data_aghy)

lik_agpe <- like(formula = fml.agpe,
                 family = "binomial",
                 Ntrials = 1,
                 data = data_agpe)

lik_elvi <- like(formula = fml.elvi,
                 family = "binomial",
                 Ntrials = 1,
                 data = data_elvi)




fit <- bru(cmp,
           lik_aghy,
           lik_agpe,
           lik_elvi,
           options = list(
             control.compute = list(dic = FALSE, waic = FALSE, cpo = FALSE),
             control.inla = list(int.strategy = "eb"),
             verbose = TRUE,
             bru_max_iter = 10)
)


saveRDS(fit, file = "~/Dropbox/Josh&Tom - shared/Endo_Herbarium/Model_Output/multispecies_fit_conservative.Rds")
fit <- readRDS(file = "~/Dropbox/Josh&Tom - shared/Endo_Herbarium/Model_Output/multispecies_fit_conservative.Rds")

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
  formula = ~ invlogit(int.species1 + year.species1 + space+
                         space.species1 + 
                         time.species1 + 
                         scorer_eval(scorer_index) + collector_eval(scorer_index)),
  n.samples = 100) 

validation.pred2 <- predict(
  fit,
  newdata = data2,
  formula = ~ invlogit(int.species2 + year.species2 +  space+
                         space.species3 + 
                         time.species2 + 
                         scorer_eval(scorer_index) + collector_eval(scorer_index)),
  n.samples = 100) 

validation.pred3 <- predict(
  fit,
  newdata = data3,
  formula = ~ invlogit(int.species3 + year.species3 + space+
                         space.species3 + 
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
# 0.7892 for liberal
# 0.7467 for conservative


# Taking alot of material from this blog post: https://inlabru-org.github.io/inlabru/articles/2d_lgcp_covars.html

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
            mean_endo = mean(Endo_status_conservative),
            mean_lon = mean(lon),
            mean_lat = mean(lat),
            sample = n(),
            se_endo = sd(Endo_status_conservative)/sqrt(sample)) %>% 
  mutate(lat_bin = case_when(mean_lat>=35 ~ paste("43"),
                             mean_lat<35 ~ paste("35")),
         lon_bin = case_when(mean_lon<=-94 ~ paste("-90"),
                             mean_lon>-94 ~ paste("-80") ))


year_trend <- ggplot(year.pred) +
  # geom_linerange(data = endo_herb_binned, aes(x = mean_year, ymin = mean_endo-se_endo, ymax = mean_endo+se_endo))+
  geom_point(data =endo_herb, aes(x = year, y = Endo_status_conservative), shape = "|")+
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

year_trend


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



################################################################################################################################
##########  Generating posterior samples for some summary statistics ###############
################################################################################################################################
fit_lib <- readRDS(file = "~/Dropbox/Josh&Tom - shared/Endo_Herbarium/Model_Output/multispecies_fit.Rds")
fit_con <- readRDS(file = "~/Dropbox/Josh&Tom - shared/Endo_Herbarium/Model_Output/multispecies_fit_conservative.Rds")

n_draws <- 500

posteriors.aghy_lib <- generate(
  fit_lib,
  formula = ~ c(year.species1_latent,
                   int.species1_latent),
                   
  n.samples = n_draws) 
posteriors.agpe_lib <- generate(
  fit_lib,
  formula = ~ c(year.species2_latent,
                int.species2_latent),
  n.samples = n_draws) 
posteriors.elvi_lib <- generate(
  fit_lib,
  formula = ~ c(year.species3_latent,
                int.species3_latent),
  n.samples = n_draws) 

posteriors.aghy_con <- generate(
  fit_con,
  formula = ~ c(year.species2_latent,
                int.species2_latent),
  n.samples = n_draws) 
posteriors.agpe_con <- generate(
  fit_con,
  formula = ~ c(year.species2_latent,
                int.species2_latent),
  n.samples = n_draws) 
posteriors.elvi_con <- generate(
  fit_con,
  formula = ~ c(year.species2_latent,
                int.species2_latent),
  n.samples = n_draws) 

colnames(posteriors.aghy_lib) <- c( paste0("iter",1:n_draws))
colnames(posteriors.agpe_lib) <- c( paste0("iter",1:n_draws))
colnames(posteriors.elvi_lib) <- c( paste0("iter",1:n_draws))

rownames(posteriors.aghy_lib) <- c("year.AGHY", "int.AGHY")
rownames(posteriors.agpe_lib) <- c("year.AGPE", "int.AGPE")
rownames(posteriors.elvi_lib) <- c("year.ELVI", "int.ELVI")

colnames(posteriors.aghy_con) <- c( paste0("iter",1:n_draws))
colnames(posteriors.agpe_con) <- c( paste0("iter",1:n_draws))
colnames(posteriors.elvi_con) <- c( paste0("iter",1:n_draws))

rownames(posteriors.aghy_con) <- c("year.AGHY", "int.AGHY")
rownames(posteriors.agpe_con) <- c("year.AGPE", "int.AGPE")
rownames(posteriors.elvi_con) <- c("year.ELVI", "int.ELVI")


posteriors_lib <- rbind(posteriors.aghy_lib, posteriors.agpe_lib, posteriors.elvi_lib)
posteriors_con <- rbind(posteriors.aghy_con, posteriors.agpe_con, posteriors.elvi_con)




posteriors_lib_df <- as_tibble(posteriors_lib, rownames = "param") %>% 
  separate_wider_delim(param, delim = ".", names = c("param_type", "spp_label"), cols_remove = FALSE) %>% 
  pivot_longer( cols = -c(param, param_type, spp_label), names_to = "iteration") %>% 
  mutate(Data = "Liberal")

posteriors_con_df <- as_tibble(posteriors_con, rownames = "param") %>% 
  separate_wider_delim(param, delim = ".", names = c("param_type", "spp_label"), cols_remove = FALSE) %>% 
  pivot_longer( cols = -c(param, param_type, spp_label), names_to = "iteration") %>% 
  mutate(Data = "Conservative")

posteriors_df <- posteriors_lib_df %>% 
  bind_rows(posteriors_con_df)



posteriors_summary <- posteriors_df %>% 
  group_by(param, spp_label, param_type, Data) %>% 
  summarize(mean = mean(value), 
            lwr = quantile(value, .025),
            upr = quantile(value, .975),
            prop_pos = sum(value>0)/500)

posterior_hist <- ggplot(filter(posteriors_df, param_type == "year"))+
  stat_halfeye(aes(x = value, y = factor(spp_label, levels = c( "ELVI", "AGPE", "AGHY")), fill = Data), normalize = "groups", breaks = 50, alpha = .4)+
  # stat_histinterval(aes(x = value, y = spp_label, fill = spp_label), breaks = 50, alpha = .6)+
  # geom_point(data = posteriors_summary, aes(x = mean, y = spp_label, color = spp_label))+
  # geom_linerange(data = posteriors_summary, aes(xmin = lwr, xmax = upr, y = spp_label, color = spp_label))+
  geom_vline(xintercept = 0)+
  facet_wrap(~param_type, scales = "free", labeller = as_labeller(c(year = "Global Year Trend")))+
  scale_fill_manual(values = c("royalblue", "grey20"))+
  scale_color_manual(values = c("royalblue", "grey20"))+
  scale_x_continuous(labels = scales::label_number(), guide = guide_axis(check.overlap = TRUE))+
  labs(x = "Posterior Estimate", y = "", fill = "Data")+
  theme_bw()

posterior_hist
# ggsave(posterior_hist, "Plots/posterior_comparison.png", width = 5, height = 5)
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

# saveRDS(vrt_aghy, file = "aghy_distribution_df.rds")
# saveRDS(vrt_agpe, file = "agpe_distribution_df.rds")
# saveRDS(vrt_elvi, file = "elvi_distribution_df.rds")
# 
# vrt_aghy <- readRDS(file = "aghy_distribution_df.rds")
# vrt_agpe <- readRDS(file = "agpe_distribution_df.rds")
# vrt_elvi <- readRDS(file = "elvi_distribution_df.rds")


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



svc.pred_aghy <- predict(fit_con, 
                         vrt_aghy, 
                         formula = ~ ( exp(year.species1+ time.species1)-1)*100)
svc.pred_agpe <- predict(fit_con, 
                         vrt_agpe, 
                         formula = ~ ( exp(year.species2+ time.species2)-1)*100)
svc.pred_elvi <- predict(fit_con, 
                         vrt_elvi, 
                         formula = ~ ( exp(year.species3+ time.species3)-1)*100)




# saveRDS(svc.pred_aghy, file = "svc.pred_aghy.Rds")
# saveRDS(svc.pred_agpe, file = "svc.pred_agpe.Rds")
# saveRDS(svc.pred_elvi, file = "svc.pred_elvi.Rds")

# svc.pred_aghy <- readRDS(file = "svc.pred_aghy.Rds")
# svc.pred_agpe <- readRDS(file = "svc.pred_agpe.Rds")
# svc.pred_elvi <- readRDS(file = "svc.pred_elvi.Rds")



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



svc_time_map <- svc_time_map_AGHY /svc_time_map_AGPE / svc_time_map_ELVI + plot_layout(ncol = 1, guides = "collect") 
# ggsave(svc_time_map, filename = "Plots/svc_time_map.png", width = 15, height = 5)


conservative_comparison_plot <- posterior_hist + svc_time_map + plot_annotation(tag_levels = "A")

ggsave(conservative_comparison_plot, filename = "Plots/conservative_comparison_plot.png", width = 10, height = 7)




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
ggsave(svc_time_map.CI, filename = "Plots/svc_time_map.CI.png", width = 15, height = 5)



################################################################################################################################
##########  Plotting the spatial Intercepts ###############
################################################################################################################################
svc.pred_space <- list()
svc.pred_space_year <- list()



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
  
  
  svc.pred_space_year[[s]] <- predict(fit_lists[[species_codes[s]]][[1]]$year_fit, 
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
  gg(svc.pred_space_year[[1]], aes(fill = mean))+
  scale_fill_viridis_c(option = "turbo", na.value = "transparent", limits = trendrange)+
  labs(title = species_names[1], fill = "", y = "Latitude", x = "Longitude")+
  theme_light()+
  theme(plot.title = element_text(face = "italic"))
# svc_space_map_AGHY

svc_space_map_AGPE<- ggplot()+
  geom_sf(data = world_map, color = "grey", linewidth = .1, fill = "#FAF9F6") +
  geom_sf(data = states_map, color = "grey", linewidth = .1, fill = "#FAF9F6") +
  gg(svc.pred_space_year[[2]], aes(fill = mean))+
  coord_sf(xlim = space_x, ylim = space_y)+
  scale_fill_viridis_c(option = "turbo", na.value = "transparent", limits = trendrange)+
  labs(title = species_names[2], fill = "% change/year", y = "Latitude", x = "Longitude")+
  theme_light()+
  theme(plot.title = element_text(face = "italic"))


# svc_space_map_AGPE


svc_space_map_ELVI<- ggplot()+
  geom_sf(data = world_map, color = "grey", linewidth = .1, fill = "#FAF9F6") +
  geom_sf(data = states_map, color = "grey", linewidth = .1, fill = "#FAF9F6") +
  gg(svc.pred_space_year[[3]], aes(fill = mean))+
  coord_sf(xlim = space_x, ylim = space_y)+
  scale_fill_viridis_c(option = "turbo", na.value = "transparent", limits = trendrange)+
  labs(title = species_names[3], fill = "% change/year", y = "Latitude", x = "Longitude")+
  theme_light()+
  theme(plot.title = element_text(face = "italic"))
# svc_space_map_ELVI



svc_space_map <- svc_space_map_AGHY + svc_space_map_AGPE + svc_space_map_ELVI +
  plot_layout(ncol = 1, guides = 'collect') + plot_annotation(tag_levels = "A")
ggsave(svc_space_map, filename = "svc_space_map_year.png", width = 6, height = 12)




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
  formula = ~ invlogit(int.species3 + year.species3 + space + space.species3 + time.species3 + 
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
  ylim(0,1) + labs(y = "Endophyte Prevalance", x = "Longitude", color = "Species", size = "Sample Size")+
  facet_wrap(~species) + theme_classic()+ theme(strip.background = element_blank(),
                                                strip.text = element_text(size = rel(1.2)))
# contemp_lon



contemp_lat <- ggplot(contemp.pred)+
  geom_point(aes(x = lat, y = endo_prev, size = sample_size), alpha = .4)+
  geom_smooth(aes(x = lat, y = endo_prev, group = species), color = "black", method = "glm",  formula = "y ~ x", method.args = list(family = "binomial" ))+
  geom_linerange(aes(x = lat, ymin = `q0.025`, ymax = `q0.975`, color = species))+
  geom_point(aes(x = lat, y = mean), shape = 4) +
  scale_color_manual(values = c(species_colors[1], species_colors[3]))+
  guides(color = "none")+
  ylim(0,1) + labs(y = "Endophyte Prevalance", x = "Latitude",  color = "Species", size = "Sample Size")+
  facet_wrap(~species) + theme_classic() + theme(strip.background = element_blank(),
                                                 strip.text = element_blank())
# contemp_lat

contemp_obspred <- ggplot(contemp.pred)+
  geom_abline(intercept = 0, slope = 1)+
  geom_linerange(aes(y = endo_prev, xmin = `q0.025`, xmax = `q0.975`, color = species))+
  geom_point(aes(x = mean, y = endo_prev), shape = 4)+
  scale_color_manual(values = c(species_colors[1], species_colors[3]))+
  lims(x = c(0,1), y = c(0,1)) + labs(y = "Observed", x = "Predicted",  color = "Species")+
  facet_wrap(~species) + theme_classic()+ theme(strip.background = element_blank(),
                                                strip.text = element_blank())
# contemp_obspred



contemp_test_plot <- contemp_lon + contemp_lat+contemp_obspred +
  plot_layout(ncol = 1, guides = "collect") + plot_annotation(tag_levels = "A")
contemp_test_plot
# ggsave(contemp_test_plot, filename = "Plots/contemp_test_plot.png", width = 10, height = 10)
###

# now looking at the ROC and AUC values for the contemporary dataset choosing only one plant from each population
# however we only have this information for AGHY
rocobj <- pROC::roc(contemp_random_sample$endo_status, contemp.pred.aghy$mean)

ROC_test_plot <- ggroc(rocobj) 
# ggsave(ROC_test_plot, filename = "Plots/ROC_test_plot.png", width = 4, height = 4)


# AUC values
rocobj$auc
# Liberal = 0.7463
# Conservative = 0.734









###### Making a map to show where the contemp predictions ###########
contemp_map <- ggplot()+
  geom_map(data = outline_map, map = outline_map, aes(map_id = region), color = "grey", linewidth = .1, fill = "#FAF9F6")+
  geom_map(data = states_shape, map = states_shape, aes(map_id = region), color = "grey", linewidth = .1, fill = NA)+
  geom_point(data = endo_herb, aes(x = lon, y = lat),shape = 4, alpha = .3, size = 1.2)+
  geom_point(data = contemp_surveys, aes(x = lon, y = lat, color = SpeciesID))+
  coord_sf(xlim = c(-109,-68), ylim = c(21,49))+
  scale_color_manual(values = c(species_colors[1], species_colors[3]))+
  theme_light()+
  theme(legend.text = element_text(face = "italic"))+
  labs(x = "Longitude", y = "Latitude", color = "Host Species")
contemp_map
ggsave(contemp_map, filename = "contemp_surveys_map.png", width = 7, height = 4)


