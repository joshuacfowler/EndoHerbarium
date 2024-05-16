## Plotting spde model predictions
## Apr 2, 2024


library(tidyverse) # for data manipulation and ggplot
library(INLA) # for fitting integrated nested Laplace approximation models
library(inlabru)
library(fmesher)

library(sf)
# library(rmapshaper)
library(terra)
library(tidyterra)
# library(ggpolypath)

invlogit<-function(x){exp(x)/(1+exp(x))}

data <- read_csv(file = "2024-04-02_endophyte_data.csv")

# View(data)

# converting the lat long to epsg 6703km in km
# define a crs
epsg6703km <- paste(
  "+proj=aea +lat_0=23 +lon_0=-96 +lat_1=29.5",
  "+lat_2=45.5 +x_0=0 +y_0=0 +datum=NAD83",
  "+units=km +no_defs"
)

data<- data %>% 
  st_as_sf(coords = c("lon", "lat"), crs = 4326, remove = FALSE) %>% 
  st_transform(epsg6703km) %>% 
  mutate(
    easting = st_coordinates(.)[, 1],
    northing = st_coordinates(.)[, 2]
  ) 




##########################################################################################
############ Setting up and running INLA model with inlabru ############################### 
##########################################################################################

##### Building a spatial mesh #####

# Build the spatial mesh from the coords for each species and a boundary around each species predicted distribution
coords <- cbind(data$easting, data$northing)



non_convex_bdry <- fm_extensions(
  data$geometry,
  convex = c(250, 500),
  concave = c(250, 500)
)

max.edge = diff(range(coords[,1]))/(50)

# For now I'm just trying to create a standard mesh, but in the full script, I add a inner boundary based on the US coastline
mesh <- fm_mesh_2d_inla(
  # loc = coords,
  boundary = non_convex_bdry, max.edge = c(max.edge, max.edge*2), # km inside and outside
  cutoff = 20,
  crs = fm_crs(data)
  # crs=CRS(proj4string(bdry_polygon))
) # cutoff is min edge


ggplot() +
  gg(data = mesh) +
  geom_point(data = data, aes(x = easting, y = northing), size = 1) +
  coord_sf()+
  theme_bw() +
  labs(x = "", y = "")





# In general, the prior on the range of the spde should be bigger than the max edge of the mesh
prior_range <- max.edge*3
# the prior for the SPDE standard deviation is a bit trickier to explain, but since our data is binomial, I'm setting it to .5
prior_sigma <- 1

spde <- inla.spde2.pcmatern(
  mesh = mesh,
  prior.range = c(prior_range, 0.05),
  
  prior.sigma = c(prior_sigma, 0.05),
  # constr = TRUE
)



# setting the random effects prior
pc_prec <- list(prior = "pcprec", param = c(1, 0.1))


# version of the model without spatially varying time slope
svc_components <- ~ Intercept(1)+ yearEffect(main = std_year, model = "linear") + # can add Intercept(1) or not because we have the spatial intercept
  space.int(geometry, model = spde)   


# formula, with "." meaning "add all the model components":
svc_formula <- presence ~ .



# Now run the model

fit <- bru(svc_components,
           like(
             formula = svc_formula,
             family = "binomial",
             Ntrials = 1,
             data = data
           ),
           options = list(
             control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE),
             control.inla = list(int.strategy = "eb"),
             verbose = TRUE
           )
)



###### Getting and plotting prediction from model #####
# Here I am showing the trend over time

min_year <- min(data$std_year)
max_year <- max(data$std_year)
preddata <- data.frame(std_year = seq(min_year,max_year, length.out = 10))

# gennerating predictions and back-transforming the standardized year variable


year.pred <- predict(
  fit,
  newdata = preddata,
  formula = ~ invlogit(Intercept + yearEffect)) 


ggplot(year.pred) +
  geom_point(data =data, aes(x = std_year, y = presence), shape = "|")+
  geom_line(aes(std_year, mean)) +
  geom_ribbon(aes(std_year, ymin = q0.025, ymax = q0.975), alpha = 0.2) +
  geom_ribbon(aes(std_year, ymin = mean - 1 * sd, ymax = mean + 1 * sd), alpha = 0.2) +
  lims(y = c(0,1))



# Up to this point, the model predictions seems to work


# Making spatial predictions

pxl <- inlabru::fm_pixels(mesh, format = "sp")# Note this is where we can mask the output according the whatever shape, such as the host distribution

pxl$std_year <- 3

ggplot()+
  gg(mesh)+
  gg(pxl, color = "red")


space_prediction <- predict(fit, 
                                pxl, formula = ~invlogit(Intercept + space.int + std_year))

dim(pxl)
dim(data)
dim(space_prediction) # The prediction is the same number of rows as fitted data, not the new pixels object


class(pxl) # The pxl object is class "SpatialPixelsDataFrame"


# and plotting it gives a mean 1:1 line and squiggly uncertainty ribbon 
ggplot()+
  gg(space_prediction, aes(fill = mean))




####### Instead just predicting on the original dataset

space_prediction <- predict(fit, 
                                data, formula = ~invlogit(Intercept + space.int))



# which gives predictions at the original data points, but I haven't been super successful at recreating this.
ggplot()+
  gg(space_prediction, aes(color = mean))




##########################################################################################
############ Running the gorillas vignette ############################### 
##########################################################################################

# I am able to recreate the gorillas example from the website

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




pxl_gorillas <- fm_pixels(gorillas$mesh, format = "sp")

# we obtain these vertices as a SpatialPointsDataFrame

ggplot() +
  gg(gorillas$mesh) +
  gg(pxl_gorillas, color = "red")



mySmooth <- predict(fit_gorillas, pxl_gorillas, ~mySmooth)


ggplot() +
  gg(gorillas$mesh) +
  gg(mySmooth, aes(color = mean), size = 3)


