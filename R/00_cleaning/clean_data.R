# Probabilistic-Sea-Ice-Thickness-Forecasts/R/00_cleaning/clean_data.R
#
# Author: Peter Gao (petergao@uw.edu) 
#
# License: (c) Peter Gao 2020, MIT
#
# This script loads, cleans, and formats data from the following sources:
#  - NSIDC Polar Stereographic Grid and Region Maps
#  - AWI-CM thickness reanalyses
#  - ECMWF SEAS5 forecasts
#
# in: data/raw_NSIDC/region_n.msk
# in: data/raw_AWI_CM/Exp_CTD_T_[YEAR]_SIC_SIT_SIV.nc for [YEAR] in 2007:2018
# in: data/raw_AWI_CM/nod2d.out
# in: data/raw_AWI_CM/fesom.initial.mesh.diag.nc
# in: data/raw_ECMWF/seas5_ens_YEAR/seas5_ens_[ENS]_sit_[YEAR][MTH]01.nc")
#     for [YEAR] in 2015:2018, [MTH] in 01:12, [ENS] in 0:24
#
# out: data/clean_NSIDC/regions_raster.rds
# out: data/clean_NSIDC/regions_key.rds
# out: data/clean_AWI_CM/AWI_CM_locs_sf.rds
# out: data/clean_AWI_CM/AWI_CM_locs_PS_sf.rds
# out: data/clean_AWI_CM/AWI_CM_area_sf.rds
# out: data/clean_AWI_CM/AWI_CM_tbl_[YEAR].rds for [YEAR] in 2007:2018
# out: data/clean_ECMWF/sit_val[YEAR][MTH].rds
#      for [YEAR] in 2015:2018, [MTH] in 01:12
# out: data/clean_ECMWF/ECMWF_AWI_CM_tbl_lead_[LEAD].rds
#      for [LEAD] in 1:7
# ------------------------------------------------------------------------------

rm(list = ls())
# LIBRARIES -------------------------------------------------------------------
library(ncdf4)
library(sp)
library(raster)
library(proj4)
library(sf)
library(lubridate)
library(tidyr)
library(dplyr)
library(spheRlab)
library(akima)
library(stringr)

# FILE MANAGEMENT --------------------------------------------------------------
# PATH TO REPO
project_home <- "~/Dropbox/Probabilistic-Sea-Ice-Thickness-Forecasts/" 

cleaning_script_path <- paste0(project_home, "R/00_cleaning/") # SCRIPTS HERE

data_path <- paste0(project_home, "data/") # STORE DATA HERE

in_NSIDC_path <- paste0(data_path, "raw_NSIDC/") # Store raw NSIDC files here
in_ECMWF_path <- paste0(data_path, "raw_ECMWF/") # Store raw ECMWF files here
in_AWI_CM_path <- paste0(data_path, "raw_AWI_CM/") # Store raw AWI-CM files here

out_NSIDC_path <- paste0(data_path, "clean_NSIDC/") # cleaned NSIDC files 
out_ECMWF_path <- paste0(data_path, "clean_ECMWF/") # cleaned ECMWF files 
out_AWI_CM_path <- paste0(data_path, "clean_AWI_CM/") # cleaned AWI-CM files 

# CONSTANTS --------------------------------------------------------------------

# Polar Stereographic Projection used by NSIDC
PS_crs <- paste("+init=epsg:3413 +units=km +proj=stere +lat_0=90 +lat_ts=70",
                "+lon_0=-45 +k=1 +x_0=0 +y_0=0 +datum=WGS84 +no_defs",
                "+ellps=WGS84 +towgs84=0,0,0")

# Polar Stereographic grid used by NSIDC
PS_ra <- raster(nrows = 448, ncols = 304, 
                xmn = -3850, xmx = 3750,
                ymn = -5350, ymx = 5850,
                crs = PS_crs)

# Northern hemisphere lat/lon
LL_crop_extent <- extent(c(0, 360, 0, 90))

# CLEAN NSIDC INFO -------------------------------------------------------------
# in: data/raw_NSIDC/region_n.msk
# out: data/clean_NSIDC/regions_raster.rds
# out: data/clean_NSIDC/regions_key.rds
source(paste0(cleaning_script_path, "clean_NSIDC.R"))

# CLEAN AWI-CM DATA -------------------------------------------------------------
# in: data/raw_AWI_CM/Exp_CTD_T_[YEAR]_SIC_SIT_SIV.nc for [YEAR] in 2007:2018
# in: data/raw_AWI_CM/nod2d.out
# in: data/raw_AWI_CM/fesom.initial.mesh.diag.nc

# out: data/clean_AWI_CM/AWI_CM_locs_sf.rds
# out: data/clean_AWI_CM/AWI_CM_locs_PS_sf.rds
# out: data/clean_AWI_CM/AWI_CM_area_sf.rds
# out: data/clean_AWI_CM/AWI_CM_tbl_[YEAR].rds for [YEAR] in 2007:2018
source(paste0(cleaning_script_path, "clean_AWI_CM.R"))

# CLEAN ECMWF DATA -------------------------------------------------------------
# in: data/clean_AWI_CM/AWI_CM_locs_PS_sf.rds
# in: data/clean_AWI_CM/AWI_CM_tbl_[YEAR].rds for [YEAR] in 2007:2018
# in: data/raw_ECMWF/seas5_ens_YEAR/seas5_ens_[ENS]_sit_[YEAR][MTH]01.nc")
#     for [YEAR] in 2015:2018, [MTH] in 01:12, [ENS] in 0:24

# out: data/clean_ECMWF/sit_val[YEAR][MTH].rds
#      for [YEAR] in 2015:2018, [MTH] in 01:12
# out: data/clean_ECMWF/ECMWF_AWI_CM_tbl_lead_[LEAD].rds
#      for [LEAD] in 1:7
source(paste0(cleaning_script_path, "clean_ECMWF.R"))