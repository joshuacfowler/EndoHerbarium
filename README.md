# EndoHerbarium
This repository houses analysis and manuscript for an investigation of the response of grass-fungal symbionts to global climate change drivers across two centuries using surveys of historic herbarium specimens
### Repository Authors: 
Josh Fowler, Jacob Moutouama, and Tom Miller

## README file last updated: 
Oct. 7, 2025

### Project Overview:
Fungal endophytes are widespread symbionts of grasses that have been shown to provide a variety of mutualistic benefits under environmental stress. These symbionts are transmitted vertically (from maternal plant to offspring) through seeds. We are quantifying how the prevalence of this interaction has shifted across the distributions of three grass host species (\emph{A. hyemalis}, \emph{A. perennans}, and \emph{Elymus virginicus}), and how the rate of those shifts is associated with rates of change in seasonal climate drivers.

The data include 2,346 herbarium specimens collected between 1824 and 2019. We sampled seed tissue from these specimens between 2017 and 2019, and then quantified endophyte presence/absence. Combined with a spatially-varying coefficient model implemented in an approximate Bayesian framework, we flexibly estimate trends across the geographic distribution of each host.

Cleaned survey data can be found at the following Dryad data repository: [https://doi.org/10.5061/dryad.rn8pk0pn0]{https://doi.org/10.5061/dryad.rn8pk0pn0} 


## Repository Folder Description:
This repository is set up with three folders:

### Analyses 
holds scripts to cobble data and run analyses
File Name  | Description
------------- | -------------
ecolab_validation_data_prep.R|compiles contemporary endophyte survey data used as test data
endo_herbarium_mesh_and_prior_sensitivity_analysis.R|performs sensitivity analyses for model choices included SPDE mesh and priors
endo_herbarium_multispecies_spatiotemporal_analysis_with_conservative_scores.R|replicates central analysis of spatially-varying temporal trends with 'Conservative' endophyte scores
endo_herbarium_multispecies_spatiotemporal_analysis|Central analysis of spatially-varying temporal trends in endophyte prevalence; performs model fitting in INLA and generates figures
endo_herbarium_records_merge.R|script to combine herbarium endophyte survey data with various specimen locality databases, along with use of ggmap package to generate locality coordinates
posthoc_climate_correlation.R|downloads climate data from PRISM, calculates observed change in climate drivers, and tests for relationships between modeled trends in endophyte prevalence and change in climate drivers
sample_size_sensitivity_analysis.R|analysis to test ability of model to capture trends in data with spatially-biased missing data. performs this task by generating simulated data, and then fitting models with know proportion of data missing in one quadrant

### Manuscript 
holds drafts and figures of manuscript. Final manuscript files are appended with "_GCB"

File Name  | Description
------------- | -------------
endo_herbarium.bib | bib file containing references
endo_herbarium_GCB.tex | tex file for compiling final submission to journal.
Other misc. files | Various aux and log files for latex, based around formatting template from \emph{The American Naturalist}.



### Plots 
holds figures for manuscript.

File Name  | Description
------------- | -------------
AGHY_climate_change_plot.png|observed change in mean and sd of seasonal temp and preciptation across distribution of \emph{A. hyemalis}
AGPE_climate_change_plot.png|observed change in mean and sd of seasonal temp and preciptation across distribution of \emph{A. perennans}
ELVI_climate_change_plot.png|observed change in mean and sd of seasonal temp and preciptation across distribution of \emph{E. virginicus}
ROC_test_plot.png|AUC performance of model on contemporary test data
ROC_training_plot.png|AUC performance of model on historic training data
climate_trends_plot_intercept.png|relationships between change in prevalence and change in seasonal climate drivers (Fig. 5)
collections_map.png|maps of historic herbarium specimen locations (Fig. 1)
collector_posterior.png|posterior estimates of collector random effects
comparison_svc_plot.png|result of sensitivity analysis of different spde priors
conservative_comparison_plot.png|result of sensitivity analysis of use of "Conservative" vs. "Liberal" scores
contemp_surveys_map.png|map of contemporary test data
contemporary_test_plot.png|map of predictions across contemporary test data (Fig. 6)
endo_status_map.png|map of herbarium specimens colored by endophyte status
finer_mesh_comparison_svc_plot.png|result of sensitivity analysis of finer mesh
initialprev_trend_plot.png|relationship between starting endophyte prevalence and rate of change in endophyte prevalence
mesh_plot.png|mesh used during model fitting
overlay_plot.png|graphical posterior predictive check
posterior_plot.png|graph of intercept and slope parameter posterior distributions
prior_comparison_year_plot.png|prior sensitivity analysis
sample_size_fulldata_plot.png|plot of simulated data as part of sensitivity analysis of spatially biased sampling
sample_size_species2data_plot.png|plot of simulated data for "Species 2" across time periods as part of sensitivity analysis of spatially biased sampling
scorer_posterior.png|posterior_estimates of scorer random effects
sim_svc_CI_plot.png|uncertainty width in posterior estimates as part of sensitivity analysis of spatially biased sampling
sim_svc_plot.png|posterior estimates as part of sensitivity analysis of spatially biased sampling
standard_mesh_comparison_svc_plot.png|result of sensitivity analysis of finer mesh
svc_space_map_year.png|map of predicted endophyte prevalence for each host species across two time points (Fig. 4)
svc_time_map.png|map of predicted trend in endophyte prevalence across each host distribution (Fig. 3)
svc_time_map_CI.png|map of posterior credible interval width of trends in endophyte prevalence
year_plot.png|global trend in endophyte prevalence for each hosts species




### Misc. 
Other files including some used for Jacob's SDM analysis to generate host range maps.














