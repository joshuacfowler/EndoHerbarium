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
2022-07-13_README_for_year_random_effects_in_MPM.txt | notes about meeting between Josh and Tom clarifying how we are lining up the year random effects of the vital rate models within the population model. summary is that surv/growth occur in year t->t1, and reproduction occurs in year t1, but this reproduction informs new recruits in the following year.
LTREB_endodemog_2021_Plant_locations_and_maps | script to compile maps of surviving plants in the plots during census based on XY coordinates

### Manuscript 
holds drafts and figures of manuscript. Final manuscript files are appended with "_GCB"

File Name  | Description
------------- | -------------
.bib | bib file containing references
endo_herbarium_GCB.tex | tex file for compiling final submission to journal.
Other misc. files | Various aux and log files for latex, based around formatting template from \emph{The American Naturalist}.



### Plots 
holds figures for manuscript.

File Name  | Description
------------- | -------------



