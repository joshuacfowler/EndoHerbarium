# Purpose: Fits spatial model of change in endophyte prevalence in herbarium specimens
# Authors: Joshua Fowler, Mallory Tucker, Ella Segal, Lani Dufresne, and Tom Miller
# Updated: Noov 11, 2022

library(tidyverse) # for data manipulation and ggplot
library(INLA) # for fitting integrated nested Laplace approximation models
library(sf)


################################################################################
############ Read in the herbarium dataset ############################### 
################################################################################

endo_herb_georef <- read_csv(file = "~/Dropbox/Josh&Tom - shared/Endo_Herbarium/DigitizedHerbariumRecords/endo_herb_georef.csv") %>%
  filter(Country != "Canada") %>%
  mutate(Spp_code = case_when(grepl("AGHY", Sample_id) ~ "AGHY",
                              grepl("ELVI", Sample_id) ~ "ELVI",
                              grepl("AGPE", Sample_id) ~ "AGPE"))

endo_herb_AGHY <- endo_herb_georef %>% 
  filter(!is.na(Endo_status_liberal_1)) %>% 
  filter(Spp_code == "AGHY") %>% 
  filter(!is.na(lon) & !is.na(year)) %>% 
  filter(lon>-110)

endo_herb_ELVI <- endo_herb_georef %>% 
  filter(!is.na(Endo_status_liberal_1)) %>% 
  filter(Spp_code == "ELVI") %>% 
  filter(!is.na(lon) & !is.na(year)) 

endo_herb <- endo_herb_georef %>% 
  filter(!is.na(Endo_status_liberal_1)) %>%
  filter(!is.na(lon) & !is.na(year)) %>% 
  filter(lon>-110)


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
                            prior.range = c(0.05, 0.01), # P(practic.range < 0.05) = 0.01
                            prior.sigma = c(1, 0.01)) # P(sigma > 1) = 0.01


stack <- inla.stack(data = list(Endo_status_liberal_1 = endo_herb$Endo_status_liberal_1),
                    A = list(A,1),
                    effects = list(s = 1:spde$n.spde, 
                                   data.frame(Intercept = rep(1, nrow(endo_herb)),
                                              year = endo_herb$year,
                                              species = endo_herb$Spp_code)),
                    tag = 'dat')

formula <- formula(Endo_status_liberal_1 ~ 0+species*year + f(s, model = spde))
I <- inla(formula = formula, family = "binomial", Ntrials = 1,
          data = inla.stack.data(stack), 
          control.predictor = list(A = inla.stack.A(stack),
                                   link=1,compute=TRUE),
          control.family = list(link = "logit"),
          verbose = FALSE,
          )

I$summary.fixed
I$summary.random
I$summary.fitted.values

# plotting the year effect
fitted <- subset(endo_herb, select=c(Spp_code,year,Endo_status_liberal_1))
fitted$Endo_status_liberal_1 <- NA
year_seq <- expand.grid(species = c("AGHY","ELVI","POAU"),year=seq(min(endo_herb_georef$year, na.rm = T), max(endo_herb_georef$year, na.rm = T)), Endo_status_liberal_1=NA)



# Trying too visualize the spatial effect
(stepsize <- 4*1/111)
nxy <- round(c(diff(range(coords[,1])), diff(range(coords[,2])))/stepsize)
projgrid <- inla.mesh.projector(spde$mesh, xlim=range(coords[,1]),
                                ylim=range(coords[,2]), dims=nxy)
xmean <- inla.mesh.project(projgrid, I$summary.random$s$mean)
xsd <- inla.mesh.project(projgrid, I$summary.random$s$sd)

# setting values outside of the mesh to NA
require(splancs)
table(xy.in <- inout(projgrid$lattice$loc,
                     cbind(coords[,1], coords[,2])))
xmean[!xy.in] <- xsd[!xy.in] <- NA

library(gridExtra)
library(lattice)
do.call('grid.arrange',
        lapply(list(xmean, xsd),
               levelplot, col.regions=terrain.colors(16),
               xlab='', ylab='', scales=list(draw=FALSE)))

