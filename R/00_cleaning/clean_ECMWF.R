# Probabilistic-Sea-Ice-Thickness-Forecasts/R/00_cleaning/clean_ECMWF.R
#
# Author: Peter Gao (petergao@uw.edu) 
#
# License: (c) Peter Gao 2020, MIT
#
# This script loads ECMWF data from .nc files obtained via 
# correspondence with Steffen Tietsche from ECMWF.
# Generates raster objects for ECMWF monthly SIT for multiple validation dates
# ------------------------------------------------------------------------------

val_yrs <- 2015:2018
# USE FOR RUNNING FROM COMMAND LINE:
# args = commandArgs(TRUE)
# # supplied at the command line
# val_yrs = (2015:2018)[as.numeric(args[1])]


file_str <- paste0(in_ECMWF_path, "seas5_ens_")

AWI_CM_locs <- readRDS(paste0(out_AWI_CM_path, "AWI_CM_locs_PS_sf.rds")) #sf
AWI_CM_locs_sp <- as_Spatial(st_geometry(AWI_CM_locs))
for (val_yr in val_yrs) {
  init_yr <- val_yr
  AWI_CM_tbl <- readRDS(paste0(out_AWI_CM_path, "AWI_CM_tbl_", val_yr, ".rds"))
  for (val_mth in 1:12) {
    val_mth_str <- stringr::str_pad(val_mth, width = 2, pad = "0")
    out_list <- list()
    for (lead in 1:4) {
      init_mth <- val_mth - lead + 1
      init_mth <- ifelse(init_mth <= 0, init_mth + 12, init_mth)
      init_mth_str <- stringr::str_pad(init_mth, width = 2, pad = "0")
      if (init_mth > val_mth) {
        init_yr <- val_yr - 1
      }
      tbl_list <- list()
      for (ens_i in 0:24) {
        in_brick <- brick(paste0(file_str, init_yr, "/seas5_ens_", ens_i,
                                 "_sit_", init_yr, init_mth_str, "01.nc"))
        val_mths <- lubridate::month(as.Date(names(in_brick), format = "X%Y.%m.%d"))
        mthly_preds <- stackApply(in_brick, val_mths, fun = mean) # calc monthly mean
        mthly_preds <- crop(mthly_preds, LL_crop_extent)
        mthly_preds_PS <- projectRaster(rotate(mthly_preds), PS_ra)
        extracted <- raster::extract(mthly_preds[[lead]],
                                     AWI_CM_locs_sp, method = 'bilinear')
        ECMWF_sub_tbl <- 
          tibble(ECMWF_heff = extracted,
                 ens_i = ens_i, 
                 AWI_CM_CELL = AWI_CM_locs$AWI_CM_CELL,
                 mth = val_mth, yr = val_yr, lead = lead)
        
        tbl_list[[ens_i + 1]] <- ECMWF_sub_tbl
      }
      lead_tbl <- bind_rows(tbl_list) %>%
        left_join(AWI_CM_tbl, by = c("AWI_CM_CELL" = "AWI_CM_CELL",
                                     "mth" = "mth", "yr" = "yr")) %>%
        filter(!is.na(AWI_CM_heff))
      out_list[[lead]] <- lead_tbl
    }
    out_tbl <- bind_rows(out_list) %>%
      filter(!is.na(ECMWF_heff))
    saveRDS(out_tbl,
            paste0(out_ECMWF_path, "sit_val",  val_yr, val_mth_str, ".rds"))
    print(paste0("saved ", out_ECMWF_path, "sit_val",  val_yr, val_mth_str, ".rds"))
  }
}


load_path <- paste0(in_ECMWF_path, "seas5_ens_")
for (yr in 2010:2018) {
  for (init_mth in 1:12) {
    init_mth_str <- stringr::str_pad(init_mth, width = 2, pad = "0")
    mthly_preds_list <- list()
    for (ens_i in 0:24) {
      print(ens_i)
      in_brick <- brick(paste0(load_dir, yr, "/seas5_ens_", ens_i,
                               "_sit_", yr, init_mth_str, "01.nc"))
      val_mths <- lubridate::month(as.Date(names(in_brick), format = "X%Y.%m.%d"))
      mthly_preds <- stackApply(in_brick, val_mths, fun = mean)
      mthly_preds <- dropLayer(mthly_preds, 8)
      mthly_preds_list[[ens_i + 1]] <- mthly_preds
    }
    full_mthly_preds <- stack(mthly_preds_list)
    ens_mean_mthly_preds <- stackApply(full_mthly_preds, indices = 1:7, fun = mean)

    ens_mean_mthly_preds <- crop(ens_mean_mthly_preds, LL_crop_extent)
    ens_mean_mthly_preds_PS <- projectRaster(rotate(ens_mean_mthly_preds), PS_ra)
    saveRDS(ens_mean_mthly_preds_PS,
            paste0(out_ECMWF_path, "ens_mean_sit_", yr, init_mth_str, ".rds"))
    print(paste0("saved ", out_ECMWF_path, "ens_mean_sit_", yr, init_mth_str, ".rds"))
  }
}

for (yr in 2011:2018) {
  for (val_mth in 1:12) {
    print(paste0(yr, val_mth))
    val_mth_str <- stringr::str_pad(val_mth, width = 2, pad = "0")
    init_mths <- ((val_mth-7+1):val_mth)%%12
    init_mths[init_mths == 0] <- 12
    layer_list <- list()
    for (i in 1:7) {
      init_mth <- init_mths[i]
      layer_i <- 8 - i
      init_mth_str <- stringr::str_pad(init_mth, width = 2, pad = "0")
      if (init_mth > val_mth) {
        in_pth <- paste0(out_ECMWF_path, "ens_mean_sit_", yr-1,
                         init_mth_str, ".rds")
      } else {
        in_pth <- paste0(out_ECMWF_path, "ens_mean_sit_", yr,
                         init_mth_str, ".rds")
      }
      layer_list[[i]] <- readRDS(in_pth)[[layer_i]]
    }
    restacked <- stack(layer_list)
    names(restacked) <- paste0("lead", 7:1)
    saveRDS(restacked,
            paste0(out_ECMWF_path, "ens_mean_sit_val",
                   yr, val_mth_str, ".rds"))
    print(paste0("saved ", out_ECMWF_path, "ens_mean_sit_val",
                 yr, val_mth_str, ".rds"))
  }
}



## COMBINE AWI_CM and ECMWF DATASETS into tbls
# Load AWI_CM tibble
AWI_CM_locs <- readRDS(paste0(AWI_CM_path, "AWI_CM_locs_sf.rds")) #sf
AWI_CM_locs_sp <- as_Spatial(st_geometry(AWI_CM_locs))
for (lead_tm in 1:3) {
  print(lead_tm)
  ECMWF_tbl_list <- list()
  for (yr in 2015:2018) {
    AWI_CM_tbl <- readRDS(paste0(AWI_CM_path, "cleaned/AWI_CM_tbl_", yr, ".rds"))
    for (mth in 1:12) {
      print(paste0(yr, "/", mth))
      mth_str <- stringr::str_pad(mth, width = 2, pad = "0")
      in_ECMWF <- readRDS(paste0(ECMWF_path, "by_val_mth/ens_mean_sit_val", 
                                 yr, mth_str, ".rds"))
      extracted <- raster::extract(in_ECMWF[[8 - lead_tm]],
                                   AWI_CM_locs_sp, method = 'bilinear')
      ECMWF_sub_tbl <- 
        tibble(ECMWF_heff = extracted,
               AWI_CM_CELL = AWI_CM_locs$AWI_CM_CELL,
               mth = mth, yr = yr, lead = lead_tm) %>%
        left_join(AWI_CM_tbl, by = c("AWI_CM_CELL" = "AWI_CM_CELL",
                                     "mth" = "mth", "yr" = "yr"))
      ECMWF_tbl_list <- c(ECMWF_tbl_list,
                          list(ECMWF_sub_tbl %>% filter(!is.na(AWI_CM_heff))))
    }
  }
  ECMWF_tbl <- do.call(bind_rows, ECMWF_tbl_list)
  saveRDS(ECMWF_tbl, paste0(out_ECMWF_path, "ECMWF_AWI_CM_tbl_lead_", lead_tm, ".rds"))
  print(paste0("saved ", out_ECMWF_path, "ECMWF_AWI_CM_tbl_lead_", lead_tm, ".rds"))
}











