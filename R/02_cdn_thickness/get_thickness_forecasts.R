# Probabilistic-Sea-Ice-Thickness-Forecasts/R/02_cdn_thickness/get_thickness_forecasts.R
#
# Author: Peter Gao (petergao@uw.edu) 
#
# License: (c) Peter Gao 2020, MIT
#
# For each initialization month, this model fits thickness modeland generates 
# thickness field forecasts conditional on contour forecasts
#
# out: output/01_contour/observed_level_15_2011_2018.rds
# out: output/01_contour/observed_level_15_2011_2018_sf.rds
#
# INLA OUTPUT
# out: output/02_cdn_thickness/[MODEL]/INLA/[INIT_TM].rds
#      for [MODEL] in c("cdn_on_edge", "uncdn"),
#          [INIT_TM] in 97:141
#
# FULL PREDICTIONS
# out: output/02_cdn_thickness/[MODEL]/pred/[INIT_TM].rds
#      for [MODEL] in c("cdn_on_edge", "uncdn"),
#          [INIT_TM] in 97:141
#
# SUMMARY OF PREDICTIONS
# out: output/02_cdn_thickness/[MODEL]/summary/[INIT_TM].rds
#      for [MODEL] in c("cdn_on_edge", "uncdn"),
#          [INIT_TM] in 97:141
#
# ------------------------------------------------------------------------------
rm(list = ls())
# CONSTANTS --------------------------------------------------------------------
THRESHOLD <- .15 # below this, set thickness = 0
INIT_TMS <- 97:141 # internal time index; time 97 is Jan 2015
MAX_EDGE <- 60 # mesh parameter
N_CONTOURS <- 25 # number of contours to condition upon
CLUSTER = T # run on cluster?
TEST_DATA = F # true for testing with small dataset
CDN_ON_EDGE <- T # toggle for Spatio-Temporal+ model

N_CLIM_TRN_YR <- 4 # climatology training period
N_TREND_TRN_MTH <- 3 # trend training period
N_TEST_MTH <- 3 # testing period

# FILE MANAGEMENT --------------------------------------------------------------
# PATH TO REPO
project_home <- "~/Dropbox/Probabilistic-Sea-Ice-Thickness-Forecasts/" 

data_path <- paste0(project_home, "data/") # STORE DATA HERE
NSIDC_path <- paste0(data_path, "clean_NSIDC/") # cleaned NSIDC files 
ECMWF_path <- paste0(data_path, "clean_ECMWF/") # cleaned ECMWF files 
AWI_CM_path <- paste0(data_path, "clean_AWI_CM/") # cleaned AWI-CM files 

if (CDN_ON_EDGE) {
  out_path <- paste0(project_home, "results/02_cdn_thickness/cdn_on_edge/")
} else {
  out_path <- paste0(project_home, "results/02_cdn_thickness/uncdn/")
}

contour_path <- paste0(project_home, "output/01_contour")
# LIBRARIES --------------------------------------------------------------------
print("Loading libraries...")
library(dplyr)
library(tidyr)
library(sf)
library(gstat)
library(raster)
library(sp)
library(INLA)
library(stringr)
library(ggplot2)
INLA:::inla.dynload.workaround()

source(paste0(project_home, "R/02_cdn_thickness/thickness_model_functions.R"))

main <- function(init_tms, n_clim_trn_yr, n_trend_trn_mth, n_test_mth, out_dir,
                 sample_from_posterior = F) {
  
  for (init_tm in init_tms) {
    last_time_pt <- Sys.time()
    start_time <- Sys.time()
    
    init_mth <- ifelse(init_tm %% 12 == 0, 12, init_tm %% 12)
    init_yr <- floor((init_tm - 1) / 12) + 2007
    load_yrs <- (init_yr - n_clim_trn_yr):init_yr 
    if ((init_mth - n_trend_trn_mth) < 0) {
      load_yrs <- c((init_yr - n_clim_trn_yr - 1), load_yrs)
    } else if (init_mth + n_test_mth > 12) {
      load_yrs <- c(load_yrs, init_yr + 1)
    }
    print("loading and cleaning data...")
    
    dt <- load_AWI_CM_data(yrs = load_yrs, test = TEST_DATA)
    
    load_time <- Sys.time() - last_time_pt
    last_time_pt <- Sys.time()
    
    AWI_CM_locs_sf <- dt$AWI_CM_locs_sf
    AWI_CM_tbl <- dt$AWI_CM_tbl
    print("fitting models...")
    INLA_out_file <- paste0(out_dir, "INLA/", init_tm, ".rds") 
    if (file.exists(INLA_out_file)) {
      print(paste0("loading ",  INLA_out_file))
      trn_res <- readRDS(INLA_out_file)
    } else {
      trn_res <- fit_models(AWI_CM_tbl, AWI_CM_locs_sf, init_tm, 
                            n_clim_trn_yr, n_trend_trn_mth,
                            n_test_mth, long_res = T)
      saveRDS(trn_res, INLA_out_file)
      print(paste0("saved ",  INLA_out_file))
    }
    
    train_time_tbl<- trn_res$time_tbl
    last_time_pt <- Sys.time()
    
    print("making predictions...")
    pred_out_file <- paste0(out_dir, "pred/", init_tm, ".rds")
    if (file.exists(pred_out_file)) {
      print(paste0("loading ",  pred_out_file))
      forecast_obj <- readRDS(pred_out_file)
    } else {
      out <- list()
      out$full_marginal_mse_tbl <- list()
      out$full_marginal_mae_tbl <- list()
      out$full_marginal_cov_tbl <- list()
      out$in_contour_marginal_mse_tbl <- list()
      out$in_contour_marginal_mae_tbl <- list()
      out$in_contour_marginal_cov_tbl <- list()
      out$vol_cov_tbl <- list()
      out$NS_marginal_cov_tbl <- list()
      out$S_marginal_cov_tbl <- list()
      for (i in 1:n_test_mth) {
        forecast_obj <- get_forecasts(trn_res, init_tm + i,
                                      n_contours = N_CONTOURS,
                                      qs = c(.7, .8, .9, .95, .99))
        
        out$full_marginal_mse_tbl[[i]] <- 
          forecast_obj$full_marginal_mse_tbl
        out$full_marginal_mae_tbl[[i]] <- 
          forecast_obj$full_marginal_mae_tbl
        out$full_marginal_cov_tbl[[i]] <- 
          forecast_obj$full_marginal_cov_tbl
        out$in_contour_marginal_mse_tbl[[i]] <- 
          forecast_obj$in_contour_marginal_mse_tbl
        out$in_contour_marginal_mae_tbl[[i]] <- 
          forecast_obj$in_contour_marginal_mae_tbl
        out$in_contour_marginal_cov_tbl[[i]] <- 
          forecast_obj$in_contour_marginal_cov_tbl
        out$vol_cov_tbl[[i]] <- forecast_obj$vol_cov_tbl
        out$NS_marginal_cov_tbl[[i]] <- forecast_obj$NS_marginal_cov_tbl
        out$S_marginal_cov_tbl[[i]] <- forecast_obj$S_marginal_cov_tbl
      }
      forecast_obj <- list()
      forecast_obj$full_marginal_mse_tbl <- 
        do.call(bind_rows, out$full_marginal_mse_tbl)
      forecast_obj$full_marginal_mae_tbl <- 
        do.call(bind_rows, out$full_marginal_mae_tbl)
      forecast_obj$full_marginal_cov_tbl <-
        do.call(bind_rows, out$full_marginal_cov_tbl)
      forecast_obj$in_contour_marginal_mse_tbl <- 
        do.call(bind_rows, out$in_contour_marginal_mse_tbl)
      forecast_obj$in_contour_marginal_mae_tbl <-
        do.call(bind_rows, out$in_contour_marginal_mae_tbl)
      forecast_obj$in_contour_marginal_cov_tbl <- 
        do.call(bind_rows, out$in_contour_marginal_cov_tbl)
      forecast_obj$vol_cov_tbl <- do.call(bind_rows, out$vol_cov_tbl)
      forecast_obj$NS_marginal_cov_tbl <-  do.call(bind_rows, out$NS_marginal_cov_tbl)
      forecast_obj$S_marginal_cov_tbl <- do.call(bind_rows, out$S_marginal_cov_tbl)
      saveRDS(forecast_obj, pred_out_file)
      print(paste0("saved ", pred_out_file))
    }
    pred_time_tbl <- forecast_obj$time_tbl
    last_time_pt <- Sys.time()
    
    print("summarizing results...")
    summary_out_file <- paste0(out_dir, "summary/", init_tm, ".rds")
    out_summary <- get_summary(trn_res, forecast_obj)
    
    gather_time <- Sys.time() - last_time_pt
    last_time_pt <- Sys.time()
    
    time_tbl <- tibble(init_tm = init_tm,
                       load = load_time,
                       gather = gather_time,
                       total = Sys.time() - start_time)
    #time_tbl <- cbind(time_tbl, train_time_tbl, pred_time_tbl)
  }
  out_summary$time_tbl <- time_tbl
  saveRDS(out_summary, summary_out_file)
  print(paste0("saved ", summary_out_file))
  return(out_summary)
}
if (CLUSTER) {
  args = commandArgs(TRUE)
  # supplied at the commandlines
  INIT_TMS = as.numeric(args[1])
}

print(INIT_TMS)
test <- main(INIT_TMS, N_CLIM_TRN_YR, N_TREND_TRN_MTH,
             N_TEST_MTH, paste0(out_path))
# if (init_tms == 69) {
#   source(paste0(path, "cluster_summarize.R"))
# }
