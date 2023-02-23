# Purpose: Fits spatial model of change in endophyte prevalence in herbarium specimens
# Authors: Joshua Fowler, Mallory Tucker, Ella Segal, Lani Dufresne, and Tom Miller
# Updated: Noov 11, 2022

library(tidyverse) # for data manipulation and ggplot
library(INLA) # for fitting integrated nested Laplace approximation models
library(inlabru)
library(sf)
library(patchwork)
library(ggmap)

invlogit<-function(x){exp(x)/(1+exp(x))}

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
  filter(lon>-110) 

endo_herb$collector_lastname <- str_replace_all(endo_herb$collector_lastname, "�", " ")
endo_herb$collector_lastname <- str_replace_all(endo_herb$collector_lastname, "xa0", " ")


endo_herb <- endo_herb_AGHY
endo_herb <- endo_herb_AGPE
endo_herb <- endo_herb_ELVI



# testing if the trends are still are present without oldest samples
# endo_herb <- endo_herb_georef %>% 
#   filter(!is.na(Endo_status_liberal_1)) %>%
#   filter(!is.na(Spp_code)) %>% 
#   filter(!is.na(lon) & !is.na(year)) %>% 
#   filter(lon>-110) %>% 
#   filter(year>1910)


summary(lm(formula(Endo_status_liberal_1 ~ 0 + Spp_code + Spp_code:year + Spp_code:year:lat + Spp_code:year:lon + Spp_code:year:lat:lon), data = endo_herb))
################################################################################
############ Setting up and running INLA model ############################### 
################################################################################
# messing aruond with INLA model fitting
# generate mesh, which is used to fit the spatial effects
# using a coarser mesh to start out as this is easier computation and sufficient to capture spatial effects
coords <- cbind(endo_herb$lon, endo_herb$lat)

non_convex_bdry <- inla.nonconvex.hull(coords, -0.03, -0.05, resolution = c(100, 100))

max.edge = diff(range(coords[,2]))/(10)
mesh10 <- inla.mesh.2d(loc = coords, max.edge = c(1,2)*max.edge,
                      boundary = non_convex_bdry,
                      offset = c(1,4),
                      cutoff = max.edge/(5))
# plot(mesh10)

mesh5 <- inla.mesh.2d(loc = coords, max.edge = c(1,2)*(max.edge/2),
                       boundary = non_convex_bdry,
                       offset = c(1,4),
                       cutoff = max.edge/(5))
# plot(mesh5)

mesh2.5 <- inla.mesh.2d(loc = coords, max.edge = c(1,2)*(max.edge/2/2),
                        boundary = non_convex_bdry,
                        offset = c(1,4),
                        cutoff = max.edge/(10))
# plot(mesh2.5)

mesh1 <- inla.mesh.2d(loc = coords, max.edge = c(.5,2)*(max.edge/2/2),
                        boundary = non_convex_bdry,
                        offset = c(1,4),
                        cutoff = max.edge/(10))
# plot(mesh1)
# points(coords, col = "red")
mesh5$n # the number of mesh vertices


# defining the spatial random effect
A <- inla.spde.make.A(mesh5, loc = coords)
spde <- inla.spde2.pcmatern(mesh = mesh5,
                            prior.range = c(.05, 0.01), # P(practic.range < 0.05) = 0.01
                            prior.sigma = c(1, 0.01)) # P(sigma > 1) = 0.01
spde.index <- inla.spde.make.index(name = "spatial.field", n.spde = spde$n.spde)

  
stack.fit <- inla.stack(data = list(Endo_status_liberal = endo_herb$Endo_status_liberal),
                    A = list(A,1),
                    effects = list(s = spde.index$spatial.field, 
                                   data.frame(Intercept = rep(1, nrow(endo_herb)),
                                              year = endo_herb$year,
                                              year2 = endo_herb$year^2,
                                              lat = endo_herb$lat,
                                              lon = endo_herb$lon,
                                              species = endo_herb$Spp_code,
                                              county = endo_herb$County,
                                              collector = endo_herb$collector_lastname)),
                    
                    tag = 'fit')
# making a data stack for the predicted value of endophyte prevalence
# need an even  grid of coordinates

# this is the set of data for which we want predictions
pred_data <- data.frame(expand.grid(Intercept = 1, 
                        lon = seq(min(coords[,1]),max(coords[,1]), length.out = 75),
                        lat = seq(min(coords[,2]),max(coords[,2]), length.out = 75),
                        year = c(1920, 1970, 2020), 
                        year2 = c(1920^2, 1970^2, 2020^2),
                        species = c("AGHY", "ELVI", "AGPE"),
                        county = NA,
                        collector = NA))

# keeping just the points that are part of the mesh's boundary
ind <- point.in.polygon(
  pred_data$lon, pred_data$lat,
  non_convex_bdry$loc[, 1], non_convex_bdry$loc[, 2]
)
coords_pred <- pred_data %>% 
  select(lon, lat)
coords_pred <- as.matrix(coords_pred[which(ind == 1),])

pred_data <- pred_data[which(ind ==1),]
  
#Making a projection matrix for the predicted values
A_pred <- inla.spde.make.A(mesh = mesh5, loc = coords_pred)



stack.pred <- inla.stack(data = list(Endo_status_liberal = NA),
                        A = list(A_pred,1),
                        effects = list(s = spde.index$spatial.field, 
                                       pred_data),
                        tag = 'pred')


# this is the data to make a plot for the effect of year at the average lat/lon
pred.y_data <- data.frame(expand.grid(Intercept = 1, 
                                    lon = mean(endo_herb$lon),
                                    lat = mean(endo_herb$lat),
                                    year = min(endo_herb$year):max(endo_herb$year), 
                                    year = (min(endo_herb$year):max(endo_herb$year))^2,
                                    species = c("AGHY", "AGPE", "ELVI"),
                                    county = NA,
                                    collector = NA))

# keeping just the points that are part of the mesh's boundary
ind.y <- point.in.polygon(
  pred.y_data$lon, pred.y_data$lat,
  non_convex_bdry$loc[, 1], non_convex_bdry$loc[, 2]
)
coords_pred.y <- pred.y_data %>% 
  select(lon, lat)
coords_pred.y <- as.matrix(coords_pred.y[which(ind.y == 1),])

pred.y_data <- pred.y_data[which(ind.y ==1),]

#Making a projection matrix for the predicted values
A_pred.y <- inla.spde.make.A(mesh = mesh5, loc = coords_pred.y)


stack.pred.y <- inla.stack(data = list(Endo_status_liberal = NA),
                         A = list(A_pred.y,1),
                         effects = list(s = spde.index$spatial.field, 
                                        pred.y_data),
                         tag = 'pred.y')



#
full_stack <- inla.stack(stack.fit, stack.pred, stack.pred.y)
full_stack <- inla.stack(stack.fit, stack.pred.y)

# full_stack <- inla.stack(stack.fit)

# formula <- formula(Endo_status_liberal_1 ~ 1 + year + year:lat + year:lon + year:lat:lon 
#                     + f(s, model = spde))
# 
# formula <- formula(Endo_status_liberal_1 ~ 0 + year + year:lat + year:lon + year:lat:lon 
#                    + f(s, model = spde))
# 
# formula <- formula(Endo_status_liberal_1 ~ 1 + lon*lat*year+
#                    + f(s, model = spde))


formula1 <- formula(Endo_status_liberal ~ 0 + species + species:year + species:year:lat + species:year:lon + species:year:lat:lon
                  + f(s, model = spde))# + f(collector, model = "iid"))
formula2 <- formula(Endo_status_liberal ~ 0 + species + species:year + species:year:lat + species:year:lon + species:year:lat:lon
                    + f(s, model = spde) + f(collector, model = "iid"))

formula3 <- formula(Endo_status_liberal ~ 0 + species + species:year +species:year:lat + species:year:lon + species:year:lat:lon + species:year2 + species:year2:lat + species:year2:lon + species:year2:lat:lon +
                    + f(s, model = spde) + f(collector, model = "iid"))
                    
                    

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
I <- I2
saveRDS(I, file = "~/Dropbox/Josh&Tom - shared/Endo_Herbarium/Model_Output/inla_spde.rds")
# improve estimates of hyper-parameters by running the following
I <- inla.hyperpar(I)

I$summary.fixed
I$summary.random


# pulling out the predicted values from the model and attaching them to our pred_data
index <- inla.stack.index(full_stack, tag = "pred")$data
pred_data$pred_mean <- I2$summary.fitted.values[index, "mean"]
pred_data$pred_lwr <- I2$summary.fitted.values[index, "0.025quant"]
pred_data$pred_upr <- I2$summary.fitted.values[index, "0.975quant"]

write_csv(pred_data, file = "pred_data.csv")

pred_data_long <- pred_data %>% 
  pivot_longer(cols = c(pred_mean, pred_lwr, pred_upr), names_to = "variable")

endo_herb_neat <- endo_herb %>% 
  mutate(binned_year= case_when(year<=1920))) %>% 
  rename(species = Spp_code) %>% 
  group_by(species, binned_year) %>% 
    summarize(mean = mean(year))
  mutate(year = case_when(binned_year == "(1.82e+03,1.92e+03]" ~ 1920,
                          binned_year == "(1.97e+03,2.02e+03]" ~ 1970,
                          binned_year == "(1.82e+03,1.92e+03]" ~ 2020))
  
write_csv(pred_data)
ggplot(data = filter(pred_data_long, variable == "pred_mean"))+
  geom_tile(aes(x = lon, y = lat, fill = value))+
  facet_wrap(variable~year+species)+
  scale_fill_gradient(
    name = "Endophyte Prevalence",
    low = "blue", high = "orange",
    breaks = c(.25,.5,.75),
    limits = c(0,1)) +
  theme_bw()

################################################################################
############ Plotting the overall trend over time ############################### 
################################################################################
# Binning the data by year for plotting

endo_herb_binned <- endo_herb %>% 
  mutate(binned_year = cut(year, breaks = 12)) %>%
  rename(species = Spp_code) %>% 
  group_by(species,binned_year) %>%   
  summarise(mean_year = mean(year),
            mean_endo = mean(Endo_status_liberal),
            sample = n())
  
  # mutate(binned_lat = cut(lat, breaks = 1), binned_year = cut(year, breaks = 12)) %>%  
  # rename(species = Spp_code) %>% 
  # group_by(species, binned_lat, binned_year) %>%   
  # summarise(mean_lat = mean(lat),
  #           mean_year = mean(year),
  #           mean_endo = mean(Endo_status_liberal_1),
  #           sample = n())

endo_herb_nice <- endo_herb %>% 
  rename(species = Spp_code) 
  
mean_lon <- min(endo_herb$lon)
mean_lat <- min(endo_herb$lat)


# pulling out the predicting effect across time
index <- inla.stack.index(full_stack, tag = "pred.y")$data
pred.y_data$pred_mean <- I1$summary.fitted.values[index, "mean"]
pred.y_data$pred_lwr <- I1$summary.fitted.values[index, "0.025quant"]
pred.y_data$pred_upr <- I1$summary.fitted.values[index, "0.975quant"]

year_trend <- ggplot()+
  geom_ribbon(data = pred.y_data, aes(x = year, ymin = pred_lwr, ymax = pred_upr, group = species, fill = species), alpha = .2)+
  geom_line( data = pred.y_data, aes(x = year, y = pred_mean, group = species, col = species), linewidth = 1)+
  geom_point(data = endo_herb_nice, aes(x = year, y = Endo_status_liberal), pch = "|", alpha = .7)+
  geom_point(data = endo_herb_binned, aes(x = mean_year, y = mean_endo, col = species, lwd = sample))+
  # facet_wrap(~species)+
  # scale_color_viridis_b()+
  # geom_ribbon(aes(x = year, ymin = lwr, ymax = upr))+
  # lims(y = c(0,1))+
  theme_classic()+
  labs(y = "Endophyte Prevalence", x = "Year", lwd = "# of Specimens", color = "Species", fill = "Species")

year_trend

year_hist <- ggplot()+
  geom_histogram(data = endo_herb_nice, aes(x = year, fill = species, group = species), bins = 150, alpha = .4)+
  # facet_wrap(~species)+
  # scale_color_viridis_b()+
  # geom_ribbon(aes(x = year, ymin = lwr, ymax = upr))+
  # lims(y = c(0,1))+
  theme_classic()+
  theme(axis.line = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank())+
  guides(fill = "none")+
year_hist



year_plot <- year_hist + year_trend + plot_layout(ncol = 1, heights = c(1,2), guides = "collect")

year_plot

####################################################################################################################
############ Mapping the rate of change ############################### 
####################################################################################################################


pred_change <- pred_data %>% 
  distinct(lon,lat,species) 

pred_change$year_slope <- ifelse(pred_change$species == "AGPE", I$summary.fixed["speciesAGPE:year", "mean"], ifelse(pred_change$species == "ELVI", I$summary.fixed["speciesELVI:year", "mean"], ifelse(pred_change$species == "AGHY", I$summary.fixed["speciesAGHY:year", "mean"], NA)))
pred_change$year.lat_slope <- ifelse(pred_change$species == "AGPE", I$summary.fixed["speciesAGPE:year:lat", "mean"], ifelse(pred_change$species == "ELVI", I$summary.fixed["speciesELVI:year:lat", "mean"], ifelse(pred_change$species == "AGHY", I$summary.fixed["speciesAGHY:year:lat", "mean"], NA)))
pred_change$year.lon_slope <- ifelse(pred_change$species == "AGPE", I$summary.fixed["speciesAGPE:year:lon", "mean"], ifelse(pred_change$species == "ELVI", I$summary.fixed["speciesELVI:year:lon", "mean"], ifelse(pred_change$species == "AGHY", I$summary.fixed["speciesAGHY:year:lon", "mean"], NA)))
pred_change$year.lat.lon_slope <- ifelse(pred_change$species == "AGPE", I$summary.fixed["speciesAGPE:year:lat:lon", "mean"], ifelse(pred_change$species == "ELVI", I$summary.fixed["speciesELVI:year:lat:lon", "mean"], ifelse(pred_change$species == "AGHY", I$summary.fixed["speciesAGHY:year:lat:lon", "mean"], NA)))

pred_change_df <- pred_change %>% 
  mutate(change = year_slope + year.lat_slope*lat + year.lon_slope*lon + year.lat.lon_slope*lat*lon) 

pred_change_means <- pred_change_df %>% 
  summarize(mean_change =  mean(change),
         sd_change = sd(change))
pred_relative_change_df <- pred_change_df %>% 
  mutate(relative_change = (change-pred_change_means$mean_change)/pred_change_means$sd_change)

write_csv(pred_relative_change_df, "pred_change_df.csv")


register_google(key = "AIzaSyDuAdpozRqmb8Sms-XfivxLi3tzlifJdMw")
map <- ggmap::get_map(zoom = 4, maptype = c("satellite"))

ggmap(map)+
  geom_tile(data = pred_change_df, aes(x = lon, y = lat, fill = change), alpha = .3)+
  facet_wrap(~species)+
  scale_fill_gradient(
    name = "Temporal Trend",
    low = "blue", high = "orange")+
  theme_bw()



####################################################################################################################
############ correlating the temporal trend with the change in climate ############################### 
####################################################################################################################

pred_change_df <- read_csv(file = "pred_change_df.csv") %>% 
  mutate(lat = as.character(round(lat, digits = 3)), lon = as.character(round(lon, digits = 3)))
prism_diff_pred_df <- read_csv(file = "prism_diff_pred_df.csv") %>% 
  mutate(lat = as.character(round(lat, digits = 3)), lon = as.character(round(lon, digits = 3)))



pred_change_merge <- pred_change_df %>% 
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
  geom_point(aes(x = value, y = relative_change, color = species))+
  geom_smooth(aes(x = value, y = relative_change, fill = species), method = "glm")+
  facet_wrap(~name, scales = "free")

summary(lm(formula = formula(relative_change ~ 0 + ppt_annual_diff + ppt_spring_diff + ppt_summer_diff + ppt_autumn_diff + ppt_winter_diff + tmean_annual_diff + tmean_spring_diff + tmean_summer_diff + tmean_autumn_diff +  tmean_winter_diff), filter(pred_change_merge,species == "AGHY")))
summary(lm(formula = formula(relative_change ~  0 + ppt_spring_diff + ppt_summer_diff + ppt_autumn_diff + ppt_winter_diff + tmean_spring_diff + tmean_summer_diff + tmean_autumn_diff + tmean_winter_diff), filter(pred_change_merge,species == "ELVI")))

hist(pred_change_merge$tmean_annual_diff)
cor.test(filter(pred_change_merge, species == "AGPE")$relative_change, filter(pred_change_merge, species == "AGPE")$tmean_annual_diff, method = c("pearson"))
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

endo_herb$yrep <- yrep
endo_herb$yrep_lwr <- yrep_lwr
endo_herb$yrep_upr <- yrep_upr

# plot to assess model fit
ggplot(endo_herb)+
  geom_density(aes((yrep)))+
  geom_density(aes(yrep_lwr), linetype = "dashed")+
  geom_density(aes(yrep_upr), linetype = "dashed")+
  geom_density(aes(Endo_status_liberal), color = "red")+
  theme_bw()+
  labs(x = "Observed vs Predicted Density")


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
