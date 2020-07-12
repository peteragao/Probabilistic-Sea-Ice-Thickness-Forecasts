# Probabilistic-Sea-Ice-Thickness-Forecasts/R/main.R
#
# Author: Peter Gao (petergao@uw.edu) 
#
# License: (c) Peter Gao 2020, MIT
#
# This script
#   - loads and cleans data
#   - fits and generates forecasts from contour model
#   - fits and generates forecasts from conditional thickness model
#   - generates plots for paper
# ------------------------------------------------------------------------------
# FILE MANAGEMENT --------------------------------------------------------------
# PATH TO REPO
project_home <- "~/Dropbox/Probabilistic-Sea-Ice-Thickness-Forecasts/" 

# CLEAN DATA -------------------------------------------------------------------
source(paste0(project_home, "R/00_cleaning/clean_data.R"))

# CONTOUR MODEL ----------------------------------------------------------------
source(paste0(project_home, "R/01_contour/get_contour_forecasts.R"))

# THICKNESS MODEL --------------------------------------------------------------
source(paste0(project_home, "R/02_cdn_thickness/get_thickness_forecasts.R"))

# SUMMARIZE RESULTS ------------------------------------------------------------
source(paste0(project_home, "R/03_evaluation/summarize_results.R"))

# GENERATE FIGURES ------------------------------------------------------------
source(paste0(project_home, "R/03_evaluation/make_paper_figures.R"))