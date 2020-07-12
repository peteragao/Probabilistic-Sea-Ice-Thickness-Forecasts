# Probabilistic-Sea-Ice-Thickness-Forecasts/R/03_evaluation/summarize_results.R
#
# Author: Peter Gao (petergao@uw.edu) 
#
# License: (c) Peter Gao 2020, MIT
#
# This script is designed to gather results for evaluating thickness 
# forecasts from the spatiotemporal two-stage thickness model
#
# out: output/02_cdn_thickness/[MODEL]/summary.rds
#      for [MODEL] in c("cdn_on_edge", "uncdn")
# out: output/ECMWF_vol_cov_tbl.rds
#
# FILE MANAGEMENT --------------------------------------------------------------
rm(list = ls())
project_home <- "~/Dropbox/Probabilistic-Sea-Ice-Thickness-Forecasts/" 
data_path <- paste0(project_home, "data/") # STORE DATA HERE
ECMWF_path <- paste0(data_path, "clean_ECMWF/") # cleaned ECMWF files 
AWI_CM_path <- paste0(data_path, "clean_AWI_CM/") # cleaned AWI-CM files 
output_path <- paste0(project_home, "output/")
obs_contour_path <- paste0(project_home, "output/01_contour/")
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
library(stars)

# CONSTANTS --------------------------------------------------------------------
THRESHOLD <- .15 # below this, set thickness = 0
# HELPFUL FUNCTIONS ------------------------------------------------------------

main <- function(init_tms) {
  summary <- list()
  ic_mse_tbl_list <- list()
  ic_mae_tbl_list <- list()
  ic_cov_tbl_list <- list()
  mse_tbl_list <- list()
  mae_tbl_list <- list()
  cov_tbl_list <- list()

  vol_mse_tbl_list <- list()
  vol_mae_tbl_list <- list()
  vol_cov_tbl_list <- list()
  NS_param_list <- list()
  S_param_list <- list()

  for (init_tm in init_tms) {
    print(init_tm)
    summary_path <- paste0(res_path, "summary/", init_tm, ".rds")
    if (file.exists(summary_path)) {
      in_file <- readRDS(summary_path)
      
      test_tms <- unique(in_file$NS_pred_tbl$AWI_CM_tm)
      test_regs <- 6:8
      
      model_mse_tbl_list <- list()
      model_mae_tbl_list <- list()
      ic_model_mse_tbl_list <- list()
      ic_model_mae_tbl_list <- list()
      
      # compare with reference forecasts
      for (ld in 1:3) {
        # add ECMWF data for comparison
        ECMWF <- 
          readRDS(paste0(AWI_CM_path,
                         "ECMWF_AWI_CM_tbl_lead_", ld, ".rds")) %>%
          dplyr::select(AWI_CM_CELL, AWI_CM_tm, ECMWF_heff) 
        test_tbl <- in_file$S_pred_tbl %>%
          ungroup() %>%
          filter(in_obs_contour | pred_y_mean > 0) %>%
          filter(lead == ld) %>%
          left_join(ECMWF, by = c("AWI_CM_tm" , "AWI_CM_CELL")) %>%
          dplyr::select(AWI_CM_tm, ECMWF_heff, AWI_CM_heff, lead) %>%
          drop_na() %>%
          pivot_longer(cols = -c(AWI_CM_tm, lead, AWI_CM_heff),
                       names_to = "model", values_to = "heff")
        model_mse_tbl_list[[ld]] <- test_tbl %>%
          group_by(AWI_CM_tm, lead, model) %>%
          summarize(mse = mean((AWI_CM_heff - heff)^2))
        model_mae_tbl_list[[ld]] <- test_tbl %>%
          group_by(AWI_CM_tm, lead, model) %>%
          summarize(mae = mean(abs(AWI_CM_heff - heff)))
        
        ic_test_tbl <- in_file$S_pred_tbl %>%
          ungroup() %>%
          filter(in_obs_contour) %>%
          filter(lead == ld) %>%
          left_join(ECMWF, by = c("AWI_CM_tm" , "AWI_CM_CELL")) %>%
          dplyr::select(AWI_CM_tm, ECMWF_heff, AWI_CM_heff, lead) %>%
          drop_na() %>%
          pivot_longer(cols = -c(AWI_CM_tm, lead, AWI_CM_heff),
                       names_to = "model", values_to = "heff")
        ic_model_mse_tbl_list[[ld]] <- test_tbl %>%
          group_by(AWI_CM_tm, lead, model) %>%
          summarize(mse = mean((AWI_CM_heff - heff)^2))
        ic_model_mae_tbl_list[[ld]] <- test_tbl %>%
          group_by(AWI_CM_tm, lead, model) %>%
          summarize(mae = mean(abs(AWI_CM_heff - heff)))
      }
      
      # evaluation for all locations 
      mse_tbl_list <-
        c(mse_tbl_list,
          list(in_file$full_marginal_mse_tbl),
          model_mse_tbl_list)
        
      mae_tbl_list <-
        c(mae_tbl_list,
          list(in_file$full_marginal_mae_tbl),
          model_mae_tbl_list)
      cov_tbl_list <- c(cov_tbl_list, list(in_file$full_marginal_cov_tbl))
      
      # evaluation only for locations in observation contour
      ic_mse_tbl_list <- c(ic_mse_tbl_list,
                           c(list(in_file$in_contour_marginal_mse_tbl),
                             ic_model_mse_tbl_list))
      ic_mae_tbl_list <- c(ic_mae_tbl_list,
                           c(list(in_file$in_contour_marginal_mae_tbl),
                             ic_model_mae_tbl_list))
      ic_cov_tbl_list <- c(ic_cov_tbl_list, 
                           list(in_file$in_contour_marginal_cov_tbl))
      
      # get params
      NS_param_tbl <- tibble(AWI_CM_tm = init_tm,
                             model = "Non-spatial",
                             s_e_med = 1/sqrt(in_file$summary_hyperpar_NS[1,4]),
                             s_w_med = 1/sqrt(in_file$summary_hyperpar_NS[2,4]),
                             rho_w_med =in_file$summary_hyperpar_NS[3,4])
      NS_param_list <- c(NS_param_list, list(NS_param_tbl))
      S_param_tbl <- tibble(AWI_CM_tm = init_tm,
                              model = "Spatial",
                              s_e_med = 1/sqrt(in_file$summary_hyperpar_S[1,4]),
                              r_med = in_file$summary_hyperpar_S[2,4],
                              s_z_med = in_file$summary_hyperpar_S[3,4],
                              rho_ST_med = in_file$summary_hyperpar_S[4,4])
      S_param_list <- c(S_param_list, list(S_param_tbl))
      
      # volume comparisons
      clim_vol_mse_tbl <- in_file$vol_cov_tbl %>%
        filter(q == min(in_file$vol_cov_tbl$q)) %>%
        filter(model == "Spatial") %>%
        mutate(model = "Climatology") %>%
        group_by(region, lead, AWI_CM_tm, model) %>%
        summarize(mse = mean((CLIM_vol - vol_mean)^2))
      clim_vol_mae_tbl <- in_file$vol_cov_tbl %>%
        filter(q == min(in_file$vol_cov_tbl$q)) %>%
        filter(model == "Spatial") %>%
        mutate(model = "Climatology") %>%
        group_by(region, lead, AWI_CM_tm, model) %>%
        summarize(mae = mean(abs(CLIM_vol - vol_mean)))
      vol_mse_tbl <- in_file$vol_cov_tbl %>%
        filter(q == min(in_file$vol_cov_tbl$q)) %>%
        group_by(region, lead, AWI_CM_tm, model) %>%
        summarize(mse = mean((AWI_CM_vol - vol_mean)^2))
      vol_mae_tbl <- in_file$vol_cov_tbl %>%
        filter(q == min(in_file$vol_cov_tbl$q)) %>%
        group_by(region, lead, AWI_CM_tm, model) %>%
        summarize(mae = mean(abs(AWI_CM_vol - vol_mean)))
      vol_cov_tbl <- in_file$vol_cov_tbl %>%
        ungroup() %>%
        group_by(region, lead, AWI_CM_tm, model, q) %>%
        summarize(cov_rate = mean(cov))
      
      vol_mse_tbl_list <- c(vol_mse_tbl_list, 
                            list(vol_mse_tbl, clim_vol_mse_tbl))
      vol_mae_tbl_list <- c(vol_mae_tbl_list, 
                            list(vol_mae_tbl, clim_vol_mae_tbl))
      vol_cov_tbl_list <- c(vol_cov_tbl_list, list(vol_cov_tbl))
    }
  } 

  summary$NS_param_tbl <- bind_rows(NS_param_list)
  summary$S_param_tbl <- bind_rows(S_param_list)
  summary$mse_tbl <- bind_rows(mse_tbl_list)
  summary$mae_tbl <- bind_rows(mae_tbl_list)
  summary$cov_tbl <- bind_rows(cov_tbl_list)
  summary$ic_mse_tbl <- bind_rows(ic_mse_tbl_list)
  summary$ic_mae_tbl <- bind_rows(ic_mae_tbl_list)
  summary$ic_cov_tbl <- bind_rows(ic_cov_tbl_list)
  summary$vol_mse_tbl <- bind_rows(vol_mse_tbl_list)
  summary$vol_mae_tbl <- bind_rows(vol_mae_tbl_list)
  summary$vol_cov_tbl <- bind_rows(vol_cov_tbl_list)
  saveRDS(summary, paste0(res_path, "summary.rds"))
  return(summary)
}


add_vol_comp_tbl <- function(init_tms) {
  summary <- readRDS(paste0(res_path, "summary.rds"))
  area_sf <- readRDS(paste0(AWI_CM_path,
                            "AWI_CM_area_sf.rds")) %>%
    mutate(area = area / 1e9) 
  AWI_CM_locs_sf <- readRDS(paste0(AWI_CM_path, 
                                   "AWI_CM_locs_PS_sf.rds")) %>%
    left_join(area_sf %>% st_set_geometry(NULL), by = "AWI_CM_CELL")
  vol_comp_tbl_list <- list()
  # compare with reference forecasts
  for (ld in 1:3) {
    ECMWF <- 
      readRDS(paste0(AWI_CM_path, "ECMWF_AWI_CM_tbl_lead_", ld, ".rds"))  %>%
      dplyr::select(AWI_CM_CELL, AWI_CM_tm, ECMWF_heff) 
    for (init_tm in init_tms) {
      print(init_tm)
      summary_path <- paste0(res_path, "summary/", init_tm, ".rds")
      if (file.exists(summary_path)) {
        in_file <- readRDS(summary_path)
        if (ld == 1) {
          model_vol_tbl <- in_file$vol_cov_tbl %>%
            filter(q == min(in_file$S_pred_tbl$q)) %>%
            group_by(region, AWI_CM_tm, lead, model) %>%
            summarize(vol_mean = mean(vol_mean),
                      AWI_CM_vol = mean(AWI_CM_vol),
                      CLIM_vol = mean(CLIM_vol)) %>%
            ungroup() %>%
            dplyr::select(region, AWI_CM_tm, lead, model, vol_mean) %>%
            rename(vol =  vol_mean) %>%
            rename(source = model)
          vol_comp_tbl_list <- c(vol_comp_tbl_list, list(model_vol_tbl))
        }

        vol_comp_tbl <- in_file$S_pred_tbl %>% 
          filter(lead == ld) %>%
          filter(q == min(in_file$S_pred_tbl$q)) %>%
          filter(in_obs_contour) %>%
          left_join(AWI_CM_locs_sf %>% st_set_geometry(NULL),
                    by = "AWI_CM_CELL") %>%
          left_join(ECMWF, by = c("AWI_CM_tm" , "AWI_CM_CELL")) %>%
          filter(region %in% 6:8) %>%
          dplyr::select(AWI_CM_CELL, area, region, AWI_CM_tm,  lead,
                        CLIM_heff, AWI_CM_heff, ECMWF_heff) 
        full_vol_comp_tbl <- vol_comp_tbl %>%
          pivot_longer(values_to = "heff", names_to = "source",
                       -c(AWI_CM_tm, area, AWI_CM_CELL, region, lead)) %>%
          group_by(AWI_CM_tm, source, lead) %>%
          summarize(vol = sum(area * heff, na.rm = T)) %>%
          ungroup() %>%
          mutate(region = 0) %>%
          mutate(source = str_sub(source, 1, -6))
        vol_comp_tbl <- vol_comp_tbl %>%
          pivot_longer(values_to = "heff", names_to = "source",
                       -c(AWI_CM_tm, area, AWI_CM_CELL, region, lead)) %>%
          group_by(region, AWI_CM_tm, source, lead) %>%
          summarize(vol = sum(area * heff, na.rm = T)) %>%
          ungroup() %>%
          mutate(source = str_sub(source, 1, -6))
        
        vol_comp_tbl_list <- c(vol_comp_tbl_list, list(full_vol_comp_tbl,
                                                       vol_comp_tbl))
        
      }
    }
  }
  
  vol_comp_tbl <- bind_rows(vol_comp_tbl_list)
  summary$vol_comp_tbl <- vol_comp_tbl
  saveRDS(summary, paste0(res_path, "summary.rds"))
  return(summary)
}

# compare all fitted results
for (model in c("uncdn", "cdn_on_edge")) {
  res_path <- paste0(project_home, "output/02_cdn_thickness/", model, "/")
  test <- main(97:141)
  test2 <- add_vol_comp_tbl(97:141)
}


# EVALUATE ECMWF COVERAGE -------------------------------------------------------

obs_contour_sf <- readRDS(paste0(obs_contour_path,
                                  "observed_level_15_2011_2018_sf.rds"))
locs_sf <- readRDS(paste0(AWI_CM_path, "AWI_CM_locs_PS_sf.rds")) %>%
  select(AWI_CM_CELL, region)
area_sf <-  readRDS(paste0(AWI_CM_path, "AWI_CM_area_sf.rds")) %>%
  st_set_geometry(NULL) %>%
  dplyr::select(AWI_CM_CELL, area)

qs <- c(.80,.9,.95)
cov_tbl_list <- list()
for (val_yr in 2015:2018) {
  for (val_mth in 1:12) {
    val_mth_str <- stringr::str_pad(val_mth, width = 2, pad = "0")
    in_file <- readRDS(paste0(ECMWF_path, "sit_val",
                              val_yr, val_mth_str, ".rds"))
    in_file_list <- list()
    for (i in unique(in_file$AWI_CM_tm)) {
      temp <- locs_sf %>% 
        left_join(in_file %>% filter(AWI_CM_tm == i),
                  by = "AWI_CM_CELL")
      obs_contour <- obs_contour_sf %>%
        filter(AWI_CM_tm == i)
      in_obs_contour = st_intersects(temp, obs_contour, sparse = F)[,1]
      temp <- temp %>%
        mutate(in_obs_contour = in_obs_contour) %>%
        st_set_geometry(NULL)
      in_file_list <- c(in_file_list, 
                        list(temp))
    }
    in_file <- bind_rows(in_file_list) %>%
      filter(!is.na(ECMWF_heff)) %>%
      filter(in_obs_contour) %>%
      filter(region %in% 6:8) %>%
      left_join(area_sf, by = "AWI_CM_CELL")
    print(paste0("loaded ", val_yr, val_mth_str, ".rds"))
    for (q in qs) {
      print(q)
      cov_tbl <- in_file %>%
        group_by(lead, region, mth, yr, AWI_CM_tm, ens_i) %>%
        summarize(ECMWF_vol = sum(ECMWF_heff * area),
                  AWI_CM_vol = sum(AWI_CM_heff * area)) %>%
        ungroup() %>%
        group_by(lead, region, mth, yr, AWI_CM_tm, AWI_CM_vol) %>%
        summarize(lb = quantile(ECMWF_vol, (1 - q) / 2),
                  ub = quantile(ECMWF_vol, 1 - (1 - q) / 2)) %>%
        ungroup() %>% mutate(cov = lb < AWI_CM_vol & ub > AWI_CM_vol,
                             q = q)
      all_cov_tbl <- in_file %>%
        group_by(lead, mth, yr, AWI_CM_tm, ens_i) %>%
        summarize(ECMWF_vol = sum(ECMWF_heff * area),
                  AWI_CM_vol = sum(AWI_CM_heff * area)) %>%
        ungroup() %>%
        group_by(lead, mth, yr, AWI_CM_tm, AWI_CM_vol) %>%
        summarize(lb = quantile(ECMWF_vol, (1 - q) / 2),
                  ub = quantile(ECMWF_vol, 1 - (1 - q) / 2)) %>%
        ungroup() %>% mutate(cov = lb < AWI_CM_vol & ub > AWI_CM_vol,
                             q = q) %>%
        mutate(region = 0)
      cov_tbl_list <- c(cov_tbl_list,
                        list(cov_tbl),
                        list(all_cov_tbl))
    }
  }
}
vol_cov_tbl <- bind_rows(cov_tbl_list)
saveRDS(vol_cov_tbl, paste0(output_path, "ECMWF_vol_cov_tbl.rds"))

cov_tbl_list <- list()
for (val_yr in 2015:2018) {
  for (val_mth in 1:12) {
    val_mth_str <- stringr::str_pad(val_mth, width = 2, pad = "0")
    in_file <- readRDS(paste0(ECMWF_path, "sit_val",
                              val_yr, val_mth_str, ".rds"))
    
    in_file_list <- list()
    for (i in unique(in_file$AWI_CM_tm)) {
      temp <- locs_sf %>%
        left_join(in_file %>% filter(AWI_CM_tm == i),
                  by = "AWI_CM_CELL")
      obs_contour <- obs_contour_sf %>%
        filter(AWI_CM_tm == i)
      in_obs_contour = st_intersects(temp, obs_contour, sparse = F)[,1]
      temp <- temp %>%
        mutate(in_obs_contour = in_obs_contour) %>%
        st_set_geometry(NULL)
      in_file_list <- c(in_file_list,
                        list(temp))
    }
    in_file <- bind_rows(in_file_list) %>%
      filter(!is.na(ECMWF_heff)) %>%
      filter(in_obs_contour) %>%
      filter(region %in% 6:8)
    
    for (q in qs) {
      cov_tbl <- in_file %>%
        group_by(lead, AWI_CM_CELL, mth, yr, AWI_CM_tm, AWI_CM_heff) %>%
        summarize(lb = quantile(ECMWF_heff, (1 - q) / 2),
                  ub = quantile(ECMWF_heff, 1 - (1 - q) / 2)) %>%
        ungroup() %>% mutate(cov = lb < AWI_CM_heff & ub > AWI_CM_heff,
                             q = q)
      cov_tbl_list <- c(cov_tbl_list,
                        list(cov_tbl))
    }
  }
}
full_cov_tbl <- bind_rows(cov_tbl_list)
summary_cov_tbl <- full_cov_tbl %>%
  filter(AWI_CM_tm >= 100) %>%
  group_by(lead, q) %>%
  summarize(cov = mean(cov)) %>%
  ungroup() %>%
  mutate(model = "ECMWF")

summary_cov_tbl_by_mth <- full_cov_tbl %>%
  filter(AWI_CM_tm >= 100) %>%
  group_by(lead, q, mth) %>%
  summarize(cov = mean(cov)) %>%
  ungroup() %>%
  mutate(model = "ECMWF")

saveRDS(summary_cov_tbl,
        paste0(output_path, "ECMWF_marginal_cov_tbl.rds"))
saveRDS(summary_cov_tbl_by_mth,
        paste0(output_path, "ECMWF_marginal_cov_tbl_by_mth.rds"))
saveRDS(full_cov_tbl,
        paste0(output_path, "ECMWF_marginal_cov_tbl_by_cell.rds"))

