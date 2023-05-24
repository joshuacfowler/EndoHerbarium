library(tidyverse)
library(MASS)
library(INLA)
library(fields)

invlogit<-function(x){exp(x)/(1+exp(x))}

##### Simulating data to see if the INLA SPDE model can pick up both increasing and decreasing trends in spatial data #####





## time and locations
t <- 12
time <- 1:t
loc <- as.matrix(expand.grid(seq(0, 3, 0.2), seq(0, 3, 0.2)))
locdist <- as.matrix(dist(loc))

## Space covariance -- sigma = 1, range = 0.3
## Matern covariance # to be compared with inla: smoothness(nu) = alpha - d/2, d = 2 here and default in inla alpha = 2.
Vs <- Matern(d =locdist, range = .1, smoothness = .1, phi = 10) # low values here are essentially more local spatial autocorrelation, whereas bigger values are larger autocorrelation


## Time covariance -- r = 0 
# setting this to 0 is essentially no temporal autocorrelation
Vt <- diag(t) 
Vt <- 1* .0^abs(row(Vt)-col(Vt))

## Cross covariance
Vc <- kronecker(Vs, Vt)

## simulate the data
set.seed(10)
xx <- crossprod(chol(Vc), rep(rnorm(n = nrow(loc) * t)))
lon <- rep(loc[,1], each = t)
lat <-  rep(loc[,2], each = t)
time <- rep(1:t,  nrow(loc))
prob_x <-  xx + 0*lon +0*lat + 0*lat*lon
prob_xt <- xx + 0*lon +0*lat + 0*lat*lon + .1*time + -.1*time*lat + .1*time*lon + -.1*time*lat*lon

## Create the time spatial data frame
simdf <- data.frame(lon = lon, lat = lat, 
                    prob= xx, 
                    prob_x = prob_x,
                    prob_xt = prob_xt,
                    datn = rbinom(n = length(xx), size = 1,prob = invlogit(prob_x)),
                    datnt = rbinom(n = length(xx), size = 1,prob = invlogit(prob_xt)),
                    time = time )

simdf_samp <- simdf %>% 
  slice_sample(n = 600)

## plot

ggplot(data = simdf)+
  geom_tile(aes(x = lon, y = lat, fill = prob)) +coord_sf()+ facet_wrap(~time) 

ggplot(data = simdf)+
  geom_tile(aes(x = lon, y = lat, fill = prob_x)) +coord_sf()+ facet_wrap(~time) 

ggplot(data = simdf)+
  geom_tile(aes(x = lon, y = lat, fill = invlogit(prob_xt))) +coord_sf()+ facet_wrap(~time) 


ggplot(data = simdf_samp)+
  geom_tile(aes(x = lon, y = lat, fill = datn)) +coord_sf()+ facet_wrap(~time) 

ggplot(data = simdf_samp)+
  geom_tile(aes(x = lon, y = lat, fill = datnt)) +coord_sf()+ facet_wrap(~time) 




ggplot(data = filter(simdf, lon == 0.0 & lat == 0.0))+
  geom_point(aes(x = time, y = datn ), method = "glm")
ggplot(data = filter(simdf, lon == 0.0 & lat == 0.0))+
  geom_smooth(aes(x = time, y = datn ), method = "glm")
#










# setting up the data to be sampled from the generate points for each species
set.seed(123)
n_samp <- 100 #number of samples per species
n_slope <- 5 # number of temporal trend scenarios
n_space <-  5 # number of spatial scenarios
n_regions <- c(2,3,4,5,6) # max number of spatial regions (squared)

year <- 0:200
slopes <- LETTERS[1:n_slope]

# longitude <- seq(-75,-100, length.out = n_samp*n_slope) + rnorm(n= n_samp*n_slope, mean = 0, sd = 1)
# latitude <- seq(30,40, length.out = n_samp*n_slope) + rnorm(n= n_samp*n_slope, mean = 0, sd = 1)


longitude <- seq(0,100, length.out = n_samp*n_slope) + rnorm(n= n_samp*n_slope, mean = 0, sd = 1)
latitude <- seq(0,100, length.out = n_samp*n_slope) + rnorm(n= n_samp*n_slope, mean = 0, sd = 1)

# sampling from those to make a set of coordinates for our fake data
coords <- data.frame(long=sample(longitude, replace = TRUE, size = n_samp*n_slope), lat=sample(latitude, replace = TRUE, size = n_samp*n_slope), year = sample(year, replace = TRUE, size = n_samp*n_slope))

coords_binned <- list()

for(i in 1:n_space){
coords_binned[[i]]<- coords %>% 
  mutate(one.region = "A",
         bin.long = cut(long, breaks = n_regions[i], labels = LETTERS[1:n_regions[i]]),
         bin.lat = as.character(as.numeric(cut(lat, breaks = n_regions[i], labels =1:n_regions[i]))))  %>% 
  mutate(region = paste0(bin.long, bin.lat),
         space_treatment = paste0("space",i))
}

coords_df <- as_tibble(rbind(coords_binned[[1]],
                            coords_binned[[2]],
                            coords_binned[[3]],
                            coords_binned[[4]],
                            coords_binned[[5]]))


ggplot(data = coords_df)+
  geom_point(aes(x = long, y = lat, color = year))+
  facet_wrap(~space_treatment)
ggplot(data = coords_df)+
  geom_point(aes(x = long, y = lat, color = region))+
  facet_wrap(~space_treatment)

ggplot(data = coords)+
  geom_point(aes(x = long, y = lat, color = slopes))

ggplot(data = coords)+
  geom_histogram(aes(x = year))+
  facet_wrap(~slopes)

coords_group <- coords %>% 
  group_by(slopes) %>% 
  summarize(n = n())


X <- model.matrix(~ 1+year, coords)
Z <- model.matrix(~ region -1, coords)
ZZ <- model.matrix(~spp:year - 1, coords)


alpha <- c(0, 0) # intercept and overall slope
gamma <- seq(from = 1, to = 0, length.out = n_spp) # species intercepts

# gamma <- rep(.5, n_spp) # species intercepts

beta <- seq(from = -.02, to= .02, length.out = n_spp) # species slopes

linpred <-  X %*% alpha + Z%*% gamma + ZZ%*%beta # + rnorm(n_samp*n_spp, 0, 100)

y <- rbinom(n = n_samp*n_spp, size = 1, prob = invlogit(linpred))

coords$y <- y

ggplot(data = coords)+
  geom_point(aes(x = year, y = y))+
  facet_wrap(~spp)

ggplot(data = filter(coords))+
  geom_point(aes(x = long, y = lat, color = y))+
  facet_wrap(~spp)
ggplot(data = filter(coords, year > 50))+
  geom_point(aes(x = long, y = lat, color = y))+
  facet_wrap(~spp)





library(lme4)
model <- glm(y ~ 0 + spp*year, data = coords, family = "binomial")
summary(model)







####### TEsting out gstat package
library(gstat)
v.sim <- vgm(200, "Mat", 30, 20, kappa = 0.1)
g.sim <- gstat(NULL, "z", z~long+lat+year,
               locations=~long+lat, beta = c(0.10, 0.20, 0.40),
               dummy=TRUE, model = v.sim, nmax=50)
yy <- predict(g.sim, newdata = coords, nsim = 1) # mydata is a data.table and has the longitude, latitude, and covariate_x
coordinates(yy) <- c("long", "lat")



model_matrix <- model.matrix(~ 1 + as.factor(spp) + spp:year, coords)

# now we will use a multivariate normal to generate spatial dependency between the points in the outcome (0/1)
# These parameters define the matern covariance matrix
kappa <- 7
nu <- 1
coords_m <- cbind(coords$long, coords$lat)
dmat <- dist(coords_m)
cor_m <- as.matrix((2^(1 - nu)/gamma(nu)) * ((kappa * dmat)^nu) * (besselK(dmat *kappa, nu)))
cor_m[1:5, 1:5]

# adding a bit of noise to the covariance matrix as sort of "measurement error"


sigma2e <- 0.3
sigma2x <- 5

diag(cor_m) <- 1
cov_m <- sigma2e * diag(n) + sigma2x * cor_m

L <- chol(cov_m)
set.seed(234)

# This is the linear predictor
betaspp <- 1
betayear <- seq(from = -100, to = 100, length.out = 9)
mu <- model_matrix[,1:9]*betaspp + model_matrix[,10:18]%*%betayear
y <-rnorm(n = n, mean = mu, sd = 50000)

coords$y <- y
head(coords)

ggplot(data = filter(coords, year>1960))+
  geom_point(aes(x = long, y = lat, color = y))+
  facet_wrap(~spp)

ggplot(data = filter(coords))+
  geom_point(aes(x = year, y = y))+
  facet_wrap(~spp)






groups <- 1:5
N <- 250
g <- factor(sample(groups, replace = TRUE, size = N), levels = groups)
x <- rnorm(N)
X <- model.matrix(~ x)
Z <- model.matrix(~ g - 1)

beta <- c(10, 2)
gamma <- rnorm(length(groups), 0, 10)

y = X %*% beta + Z%*% gamma + rnorm(N, 0, 0.3)


library(lme4)
d <- data.frame(y = y, x = x, g = g)
model = lmer(y ~ x + (1|g), data = d)

summary(model)

Grid <- data.frame(x.easting = sample(longitude, size = 100), x.northing = sample(latitude, size= 100))
n <- nrow(Grid)

# Explanatory variables and coefficients
x1 <- rnorm(n) %>% round(2)
x2 <- rnorm(n) %>% round(2)




# Spatial field
distance <- as.matrix(dist(Grid))
omega <- MASS::mvrnorm(n     = 1, 
                       mu    = rep(0,n), 
                       Sigma = 0.4 * exp(-0.1 * distance))

eta <- x1 + x2 + omega

d <- Grid %>% 
  mutate(Y_normal = rnorm(n, eta, sd = 0.1) %>% round(2),
         Y_pois   = rpois(n, exp(eta)),
         trials   = 10,
         Y_binom  = rbinom(n = n, size = trials, prob = 1/(1+exp(-eta))),
         Y_bern = rbinom(n = n, size = 1, prob = 1/(1+exp(-eta))),
         x1       = x1,
         x2       = x2) %>% 
  sample_n(size = 1000)
  

ggplot(d, aes(x = x.easting, y = x.northing, color = as.factor(Y_bern))) + 
  geom_point() +
  theme_bw()

##### Simulating data for space-time interaction models #####

pred_data <- data.frame(expand.grid(Intercept = 1, 
                                    lon = seq(0,3, .2),
                                    lat = seq(0,3,.2),
                                    year = seq(1900,2010, by = 1), 
                                    species = c("AGHY", "AGPE", "ELVI")))

time <- c(1:20)
loc <- as.matrix(expand.grid(seq(0,3, .2),seq(0,3,.2)))
locdist <- as.matrix(dist(loc))
# Simulate matern covariance
library(geoR)
Vs <- matern(locdist, phi = 0.3, kappa = 1)
# image(Vs)

# Covariance across time
Vt <- diag(20)
Vt <- 1 + .1^abs(row(Vt)-col(Vt))

# cross co-variance (Kronecker)

Vc <- kronecker(Vs, Vt) 

## simulate the data
set.seed(123)
xx <- crossprod(chol(Vc), rep(rnorm(nrow(loc) * 20)))


## Create the time spatial data frame
simdf <- data.frame(lat = rep(loc[,1], each = 20), lon = rep(loc[,2], each = 20), 
                    dat= xx, datn = as.numeric(rbernoulli(length(xx), p = invlogit(xx))), year = rep(1:20,  nrow(loc)) )

lattice::levelplot(datn ~ lat + lon|year, data = filter(simdf, year == 1:10))

subset_simdf <- simdf[sample(nrow(simdf), 100),]

ggplot(subset_simdf)+
  geom_point(aes(x = lon, y = lat, color = as.factor(datn)))


##### Building the SPDE model #####

## Select time knots and generate time mesh
knots = seq(1, 20, length = 20)
mesh1 = inla.mesh.1d(loc = knots, degree = 2, boundary = "free")

## generate space mesh 
coords <- cbind(subset_simdf$lon, subset_simdf$lat)
mesh2 <- inla.mesh.2d(loc = coords, offset = c(0.1, 0.4), max.edge = 0.3 )
plot(mesh2)

# make the matern SPDE 
A <- inla.spde.make.A(mesh2, loc = coords)

spde <- inla.spde2.pcmatern(mesh = mesh2,
                            prior.range = c(0.05, 0.01), # P(practic.range < 0.05) = 0.01
                            prior.sigma = c(1, 0.01)) # P(sigma > 1) = 0.01
spde.index <- inla.spde.make.index(name = "spatial.field", n.spde = spde$n.spde)


stack.fit <- inla.stack(data = list(Endo_status_liberal_1 = endo_herb$Endo_status_liberal_1),
                        A = list(A,1),
                        effects = list(s = spde.index$spatial.field, 
                                       data.frame(Intercept = rep(1, nrow(endo_herb)),
                                                  year = endo_herb$year,
                                                  lat = endo_herb$lat,
                                                  lon = endo_herb$lon,
                                                  species = endo_herb$Spp_code,
                                                  collector = endo_herb$collector_lastname)),
                        
                        tag = 'fit')

