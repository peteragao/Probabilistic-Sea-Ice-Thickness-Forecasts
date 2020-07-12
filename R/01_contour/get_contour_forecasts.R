# Probabilistic-Sea-Ice-Thickness-Forecasts/R/01_contour/get_contour_forecasts.R
#
# Author: Peter Gao (petergao@uw.edu) 
#
# License: (c) Peter Gao 2020, MIT
#
# This script generates sea ice contour forecasts for 2015-2018
#
# in: data/bootstrapV3_1/bt_[YEAR][MTH]_f11_v3.1_n.bin
#     for [YEAR] in 2015:2018, [MTH] in 01:12
#
# out: output/01_contour/observed_level_15_2011_2018.rds
# out: output/01_contour/observed_level_15_2011_2018_sf.rds
#
# MAPS OF PROBABILITY OF ICE
# out: output/01_contour/cont_prob_init_[YEAR]_[MTH].rds
#      for [YEAR] in 2015:2018, [MTH] in 01:12
#
# SAMPLED CONTOURS
# out: output/01_contour/sample_contours_init_[YEAR]_[MTH].rds
#      for [YEAR] in 2015:2018, [MTH] in 01:12
#
# CONDITIONING POINTS
# out: output/01_contour/ice_edges_init_[YEAR]_[MTH].rds
#      for [YEAR] in 2015:2018, [MTH] in 01:12
#
# ------------------------------------------------------------------------------
rm(list = ls())
# LIBRARIES --------------------------------------------------------------------
library(dplyr)
library(tidyr)
library(sf)
library(gstat)
library(raster)
library(sp)
library(INLA)
library(stringr)
library(ggplot2)
library(IceCast) 
INLA:::inla.dynload.workaround()
# FILE MANAGEMENT --------------------------------------------------------------
project_home <- "~/Dropbox/Probabilistic-Sea-Ice-Thickness-Forecasts/" 
data_path <- paste0(project_home, "data/") # STORE DATA HERE
in_BS_path <- paste0(data_path, "bootstrapV3_1/")
output_path <- paste0(project_home, "output/01_contour")
# CONSTANTS --------------------------------------------------------------------
obs_start_year <- 2011
obs_end_year <- 2018
obs_yrs <- 2011:2018
level = 15
bs_version <- 3.1
dat_type_obs <- "bootstrap"
PS_crs <- paste("+init=epsg:3413 +units=km +proj=stere +lat_0=90 +lat_ts=70",
                "+lon_0=-45 +k=1 +x_0=0 +y_0=0 +datum=WGS84 +no_defs",
                "+ellps=WGS84 +towgs84=0,0,0")

# GET OBSERVED CONTOURS --------------------------------------------------------
# pull all Bootstrap SIC observations 
obs <- read_monthly_BS(start_year = obs_start_year, end_year = obs_end_year,
                       version = bs_version, file_folder = in_BS_path)
obs_reg_list <- list()

# loop through observation array and convert to sf objects
for (i in 1:length(obs_yrs)) {
  for (mth in 1:12) {
    print(paste0(mth, "/", i))
    # IceCast get_region function
    obs_reg <- get_region(dat = obs[i, mth,,], 
                          dat_type = "bootstrap", level = level)
    obs_reg <- st_sf(st_combine(st_as_sf(obs_reg)))
    obs_reg_list  <- c(obs_reg_list,  list(obs_reg))
  }
}
saveRDS(obs_reg_list, paste0(output_path, "observed_level_", 
                             level, "_", obs_start_year, 
                             "_", obs_end_year, ".rds"))

# convert to sf object for convenience
obs_reg_sf <- do.call(rbind, obs_reg_list) %>%
  rename(contour = 'st_combine.st_as_sf.obs_reg..') %>%
  mutate(mth = rep(1:12,  length(obs_yrs)),
         yr = rep(obs_yrs, each = 12),
         AWI_CM_tm = (yr - 2007) * 12 + mth)
st_crs(obs_reg_sf) <- PS_crs
saveRDS(obs_reg_sf, paste0(output_path, "observed_level_", 
                           level, "_", obs_start_year, 
                           "_", obs_end_year, "_sf.rds"))

# FIT CONTOUR MODEL AND GENERATE FORECASTS -------------------------------------
# CONSTANTS --------------------------------------------------------------------
set.seed(563)

n_train_mth <- 1
n_test_mth <-  3
n_train_yrs <- 10
test_yrs <- 2015:2018
init_months <- 1:12

# USE FOR RUNNING FROM COMMAND LINE:
# args = commandArgs(TRUE)
# # supplied at the command lines
# init_months <- as.numeric(args[1])

eps = .01 # tolerance for logit transform
buff = 0
level <- 15 # grid cells with SIC < level are ice free
n_gen <- 50 # NUMBER OF FORECAST

bs_version <- 3.1
dat_type_obs <- "bootstrap"

# USE INLA SAMPLE TO GENERATE CONTOURS -----------------------------------------
# samples based on fitted inla_res object and sample object
# obtained from inla.posterior.sample
gen_cont_INLA <- function(r, sample, inla_res, test_mth, test_yr,
                          reg_info, n_gen, rand, zeroes, ones) {
  
  stopifnot(n_gen >= 2)
  
  # initial info
  n_lines <- length(rand)
  angs_r <- reg_info$angs[[r]]
  y_max <- sapply(reg_info$sec_lengths[[r]], sum) #max lengths (with land)
  
  # identify which latent variables to pull out
  sample_ind <- which(inla_res$.args$data$inla_yr == test_yr & 
                        inla_res$.args$data$inla_mth == test_mth &
                        inla_res$.args$data$icecast_line_id %in% which(rand))
  
  # pull out transformed latent variables (z) 
  z_gen <- do.call(rbind, lapply(sample, function(x) x$latent[sample_ind]))
  
  # add measurement error to obtain sample predictions w = z + delta
  w_gen <- z_gen + rnorm(length(sample_ind), 
                         sd = sqrt(1/inla_res$summary.hyperpar[1,1]))
  
  # transform back to proportions
  prop_gen <- ilogit(w_gen)
  
  # convert to appropriate line segment lengths using IceCast
  y_gen <- apply(prop_gen, 1, function(x){prop_to_y(prop = x, r, reg_info,
                                                    inds = which(rand))})
  
  # drop lines that are very close to zero length
  y_gen[y_gen <= 12.5] <- 0
  if (sum(rand) == 1) {
    y_gen <- matrix(y_gen, ncol = n_gen)
  }
  fixed_pts <- matrix(NA, nrow = n_lines, ncol = 2)
  fixed_pts[zeroes, ] <- reg_info$start_coords[[r]][zeroes,]
  fixed_pts[ones, ] <- reg_info$start_coords[[r]][ones,] +
    cbind(y_max[ones]*cos(angs_r[ones]),
          y_max[ones]*sin(angs_r[ones]))
  # make contours out of x's and y's
  # - conts_r list holds generated contours
  # - edge_coords list holds points on the ice edge 
  #               where ice should be close to zero
  conts_r <- list()
  edge_coords <- list()
  
  for (k in 1:n_gen) {
    
    #random points
    angs_r <- reg_info$angs[[r]]
    
    new_pts <- reg_info$start_coords[[r]][rand,] +
      cbind(y_gen[,k]*cos(angs_r[rand]),
            y_gen[,k]*sin(angs_r[rand]))
    # transform new_pts to sf object
    prop_k <- y_gen[, k] / y_max[rand]
    edge_coords[[k]] <- new_pts %>%
      as.data.frame %>%
      filter(prop_k > .1 & prop_k < .9) %>%
      sf::st_as_sf(coords = c(1,2)) %>%
      mutate(traj_id = k)
    
    out_pts <- fixed_pts
    out_pts[rand, ] <- new_pts
    
    # make new polygons
    conts_r[[k]] <- make_polygons(r, my_end = out_pts, poly_name = "new")
  }
  
  return(list(conts_r = conts_r,
              edge_coords = do.call(rbind, edge_coords)))
}

lines_to_fit_INLA <- function(prop_tilde, eps, buff = 0) {
  #find indices on boundaries
  n_lines <- nrow(prop_tilde)
  at_lb <- apply(prop_tilde, 1, function(x){all(x == logit(eps))})
  at_ub <- apply(prop_tilde, 1, function(x){all(x == logit(1 - eps))})
  
  rand <- !at_lb & !at_ub
  
  return(list("rand" = rand, "zeroes" = at_lb, "ones" = at_ub))
}
# LOOP THROUGH TEST YEARS FOR SELECTED MONTH -----------------------------------
for (init_month in init_months)
  mth_str <- str_pad(init_month, 2, "left", "0")
  for (test_yr in test_yrs) {
    
    print(test_yr)
    
    # which months to condition on when making forecasts
    train_mths <- (init_month - n_train_mth + 1):(init_month)
    
    # which months to generate forecasts for
    test_mths <- (init_month + 1):(init_month + n_test_mth)
    
    # years of the training months
    trend_train_yrs <- rep(test_yr, n_train_mth)
    trend_train_yrs <- 
      ifelse(train_mths <= 0, trend_train_yrs - 1, trend_train_yrs)
    
    # years of the testing months
    trend_test_yrs <- rep(test_yr, n_test_mth)
    trend_test_yrs <- 
      ifelse(test_mths > 12, trend_test_yrs + 1, trend_test_yrs)
    
    # loop over the starts/ends of years
    test_mths <- ifelse(test_mths > 12, test_mths - 12, test_mths)
    train_mths <- ifelse(train_mths <= 0, train_mths + 12, train_mths)
    
    # combine all months/years needed for fitting/generating contour forecasts
    model_mths <- c(train_mths, test_mths)
    model_yrs <- c(trend_train_yrs, trend_test_yrs)
    n_model_mth <- length(model_mths)
    
    # w_clim           list containing mean of observed transformed proportions 
    #                  over all years
    # w_train          list of observed transformed proportions over
    #                  all years
    # regs_to_fit_list list of information about which regions to fit for each
    #                  month
    w_clim <- list()
    w_train <- list()
    regs_to_fit_list <- list()
    for (i in 1:n_model_mth) {
      w_clim[[i]] <- list()
      w_train[[i]] <- list()
      
      train_start_year <- model_yrs[i] - n_train_yrs 
      train_end_year <- model_yrs[i] - 1
      forecast_year <- model_yrs[i]
      
      # get observations for desired model month (IceCast)
      obs <- read_monthly_BS(start_year = train_start_year, 
                             end_year = forecast_year,
                             version = bs_version, file_folder = in_BS_path)
      
      # convert observed sea ice concentration to y lengths (IceCast)
      y_train <- find_y(start_year = train_start_year, end_year = forecast_year,
                        obs_start_year = train_start_year,
                        pred_start_year = NULL, 
                        observed = obs[, model_mths[i], ,],
                        predicted = NULL, reg_info, month, level,
                        dat_type_obs, dat_type_pred, obs_only = TRUE)
      
      # Determine which regions to fit, and the 'full' polygon
      # convert observations to proportions and logit of proportions
      y_train_only <- lapply(y_train$obs, function(x) x[1:n_train_yrs])
      temp <- to_fit(y_obs = y_train_only, reg_info)
      regs_to_fit_list[[i]] <- list()
      regs_to_fit_list[[i]]$to_fit <- temp$regs_to_fit
      regs_to_fit_list[[i]]$full <- temp$full
      full <- temp$full
      
      #convert lengths to proportions and transformed proportions
      prop_train <- y_to_prop(y = y_train$obs, 1:5, reg_info)
      prop_train <- lapply(prop_train, function(y){sapply(y, function(x){x})})
      
      # loop over all regions
      for (r in 1:5) {
        ub_ind <- which(prop_train[[r]] >= 1 - eps)
        prop_train[[r]][ub_ind] <- 1 - eps
        lb_ind <- which(prop_train[[r]] <= eps)
        prop_train[[r]][lb_ind] <- eps
        w_train[[i]][[r]] <- logit(prop_train[[r]])
        w_clim[[i]][[r]] <- rowMeans(logit(prop_train[[r]])[, 1:(n_train_yrs)])
        # include NA values for regions not fit
        if (!(r %in% regs_to_fit_list[[i]]$to_fit)) {
          w_train[[i]][[r]][,] <- NA
        }
      }
      
    }
    
    # prepare data for use with INLA
    # for each region, obs_arrays_list contains array representing observed w val
    # for each region, clim_arrays_list contains array representing avg w val
    obs_arrays_list <- list()
    clim_arrays_list <- list()
    for (r in 1:5) {
      n_lines <- nrow(w_train[[1]][[r]])
      obs_arrays_list[[r]] <- 
        array(NA, dim = c(n_lines,
                          n_model_mth, 
                          n_train_yrs + 1))
      clim_arrays_list[[r]] <- 
        array(NA, dim = c(n_lines,
                          n_model_mth, 
                          n_train_yrs + 1))
      for (i in 1:n_model_mth) {
        for (j in 1:(n_train_yrs + 1)) {
          obs_arrays_list[[r]][, i, j] <- w_train[[i]][[r]][, j] 
          clim_arrays_list[[r]][, i, j] <- w_clim[[i]][[r]]
        }
      }
    }
    
    # store info about which regions to fit and which lines are random each month
    fit_info_list <- list()
    for (r in 1:5)  {
      fit_info_list[[r]] <- list()
      for (i in 1:n_model_mth) {
        if (r %in% regs_to_fit_list[[i]]$to_fit) {
          fit_info_list[[r]][[i]] <- 
            lines_to_fit_INLA(prop_tilde = w_train[[i]][[r]],
                              eps = eps, buff = buff)
          # don't fit region if only one random line
          if (sum(fit_info_list[[r]][[i]]$rand) == 1) {
            regs_to_fit_list[[i]]$to_fit <- 
              setdiff(regs_to_fit_list[[i]]$to_fit, r)
            fit_info_list[[r]][[i]] <- list(rand = rep(F, nrow(w_train[[i]][[r]])))
          }
        } else {
          fit_info_list[[r]][[i]] <- list(rand = rep(F, nrow(w_train[[i]][[r]])))
        }
      }
      # identify lines to fit by combining all random line indices
      fit_info_list[[r]]$inla_fit_ind <- 
        sort(unique(unlist(lapply(fit_info_list[[r]], function(x) which(x$rand)))))
    }
    
    print("starting INLA fit")
    res <- list()
    # loop over all regions that we might want to fit
    # model fit separately for each region r
    for (r in unique(unlist(lapply(regs_to_fit_list, function(x) x$to_fit))))  {
      start_time <- proc.time()
      # which lines to fit 
      fit_ind <- fit_info_list[[r]]$inla_fit_ind 
      n_lines <- length(fit_ind)
      print('formatting data')
      
      # pull data for this region
      clim_fit <- clim_arrays_list[[r]][fit_ind, ,]
      obs_fit <- obs_arrays_list[[r]][fit_ind, ,]
      
      # generate tibble for use with INLA
      inla_df <- tibble(clim = as.vector(clim_fit),
                        obs = as.vector(obs_fit),
                        inla_yr = rep(1:(n_train_yrs + 1), 
                                      each = n_lines * n_model_mth),
                        inla_line_id = rep(1:n_lines,
                                           (n_train_yrs + 1) * n_model_mth),
                        icecast_line_id = rep(fit_ind, 
                                              (n_train_yrs + 1) * n_model_mth),
                        inla_mth = rep(rep(1:n_model_mth, each = n_lines), 
                                       (n_train_yrs + 1)))
      
      # remove test data (to be forecasted) from view
      inla_df <- inla_df %>%
        mutate(obs = ifelse((inla_yr == (n_train_yrs + 1) & 
                               inla_mth > n_train_mth), NA, obs)) 
      is_variable <- sapply(1:nrow(inla_df), function(i)
        inla_df$icecast_line_id[i] %in% which(fit_info_list[[r]][[inla_df$inla_mth[i]]]$rand))
      # remove invariant lines from view
      inla_df <- inla_df %>%
        mutate(obs = ifelse(is_variable, obs, NA))
      
      print('running INLA')
      
      # pc priors chosen so P(sigma > 1) = 0.01
      # and                 P(rho > 0) = 0.9
      prec_pc_prior= list(prec = list(prior = "pc.prec", param = c(1, 0.01)))
      ar1_pc_prior <- list(theta = list(prior = 'pccor1', param = c(0, 0.9)))
      
      # RUN INLA
      inla_res <-
        inla(obs ~ 0 + offset(clim) + 
               f(inla_line_id, model = "ar1", 
                 group = inla_mth, # ar1 across months
                 hyper = ar1_pc_prior, 
                 control.group = list(model = 'ar1', hyper = ar1_pc_prior),
                 replicate = inla_yr), # iid across years
             data = inla_df,
             control.family = list(hyper = prec_pc_prior), # prior for obs model
             control.predictor = list(compute=TRUE), # compute pred with offset
             control.compute = list(config=TRUE)) # for inla.posterior.sample
      
      print('saving results')
      out <- list()
      out$res <- inla_res
      res[[r]] <- out
      end_time <- proc.time()
      elapse_time <- end_time - start_time
      print(sprintf("INLA for region %i finished, elapsed time %f", 
                    r, elapse_time[3]))
    }
    # not stored locally to save space
    # uncomment if needed
    # saveRDS(res, paste0(output_path, "INLA_res", test_yr, "_", mth_str, ".rds"))
    # print(paste0("saved ", output_path, "INLA_res", test_yr, "_", mth_str, ".rds"))
    
    # sample n_gen contour trajectories
    sample_list <- list()
    for (r in unique(unlist(lapply(regs_to_fit_list, function(x) x$to_fit)))) {
      sample_list[[r]] <- inla.posterior.sample(n = n_gen, res[[r]]$res, 
                                                selection = list(Predictor = 0))
    }
    
    # generate ice edge points sf object
    #          contour sf object
    #          ice probability raster object
    edge_sf_list <- list()
    cont_sf_list <- list()
    cont_prob_list <- list()
    
    for (i in 1:n_model_mth)  {
      test_mth <- model_mths[i]
      
      #generate contours
      indiv_conts <- list()
      edge_coords <- list()
      for (r in regs_to_fit_list[[i]]$to_fit) {
        fit_info <- fit_info_list[[r]][[i]]
        
        # get contours based on sample_list
        out <- gen_cont_INLA(r, sample_list[[r]], res[[r]]$res,
                             test_mth = i, test_yr = n_train_yrs + 1,
                             reg_info, n_gen, rand = fit_info$rand, 
                             zeroes = fit_info$zeroes, ones = fit_info$ones)
        
        indiv_conts[[r]] <- out$conts_r
        edge_coords[[r]] <- out$edge_coords
        print(sprintf("contours for region %i generated", r))
      }
      # merge contours combining results from each region r
      conts <- merge_conts(conts = indiv_conts, full = regs_to_fit_list[[i]]$full)
      # compute probability map
      cont_prob_list[[i]] <- prob_map(merged = conts)
      
      # create sf objects
      cont_sf <- lapply(conts, function(x) st_sf(st_combine(st_as_sf(x))))
      edge_sf <- do.call(rbind, edge_coords) %>%
        mutate(yr = test_yr,
               mth =  test_mth)
      edge_sf_list[[i]] <- edge_sf
      cont_sf <- do.call(rbind, cont_sf) %>%
        rename(contour = st_combine.st_as_sf.x..) %>%
        mutate(yr = test_yr,
               mth =  test_mth,
               traj_id = 1:n_gen)
      cont_sf_list[[i]] <- cont_sf
    }
    edge_sf <- do.call(rbind, edge_sf_list)
    cont_sf <- do.call(rbind, cont_sf_list)
    
    # save desired objects
    saveRDS(cont_prob_list, paste0(output_path, "cont_prob_init_",
                                   test_yr, "_", mth_str, ".rds"))
    print(paste0("saved ", output_path, "cont_prob_init_",
                 test_yr, "_", mth_str, ".rds"))
    saveRDS(edge_sf, paste0(output_path, "ice_edges_init_", 
                            test_yr, "_", mth_str, ".rds"))
    print(paste0("saved ", output_path, "ice_edges_init_", 
                 test_yr, "_", mth_str, ".rds"))
    saveRDS(cont_sf, paste0(output_path, "sample_contours_init_",
                            test_yr, "_", mth_str, ".rds"))
    print(paste0("saved ", output_path, "sample_contours_init_",
                 test_yr, "_", mth_str, ".rds"))
  }


