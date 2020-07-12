
# HELPFUL FUNCTIONS ------------------------------------------------------------
#
# COMPUTE QUANTILES OF MIXTURE OF NORMALS --------------------------------------
#' Mixture of normal distributions
#'
#' @param q vector of quantiles
#' @param wts vector of weights of mixture (sum to 1)
#' @param means vector of mean parameters
#' @param sds vector of sd parameters
#'
#' @return distribution function
#' @export
#'
#' @examples 
#' pnormmix(.5, c(.2, .8), c(0, 1), c(0, 1))
pnormmix <- function(q, wts, means, sds) {
  sum(wts * pnorm(q, mean = means, sd = sds))
}

#' Mixture of normal distributions
#'
#' @param p vector of probabilities
#' @param wts vector of weights of mixture (sum to 1)
#' @param means vector of mean parameters
#' @param sds vector of sd parameters
#'
#' @return quantile function
#' @export
#'
#' @examples
#' qnormmix(.5, c(.2, .8), c(0, 1), c(0, 1))
qnormmix <- function(p, wts, means, sds) {
  df <- function(q) pnormmix(q, wts, means, sds) - p
  return( uniroot(df, c(-5000, 20000))$root ) 
}

#
# LOAD AWI-CM DATA -------------------------------------------------------------
# test <- load_AWI_CM_data(test = T)
# test <- load_AWI_CM_data(test = F)
load_AWI_CM_data <- function(yrs = 2010:2016,  test = F) {
  
  # area sf
  area_sf <- readRDS(paste0(AWI_CM_path, "cleaned/AWI_CM_area_sf.rds")) %>%
    mutate(area = area / 1e9) 
  
  # tibble of CMST estimate locations
  AWI_CM_locs_sf <- readRDS(paste0(AWI_CM_path, "cleaned/AWI_CM_locs_PS_sf.rds")) %>%
    left_join(area_sf %>% st_set_geometry(NULL), by = "AWI_CM_CELL")
  
  # tibble of CMST estimates 2010-2016
  tbl_list <- list()
  for (i in yrs) {
    tbl_list <- 
      c(tbl_list, list(readRDS(paste0(AWI_CM_path, "cleaned/AWI_CM_tbl_", i, ".rds"))))
  }
  AWI_CM_tbl <- bind_rows(tbl_list)
  
  AWI_CM_locs_sf <- AWI_CM_locs_sf %>%
    filter(AWI_CM_CELL %in% AWI_CM_tbl$AWI_CM_CELL)
  
  # test set for trying out code
  if (test) {
    AWI_CM_locs_sf <- AWI_CM_locs_sf %>%
      st_crop(y = c(xmin = -100, ymin = -100, xmax = 100, ymax = 100))
  }
  
  # identify locations where ice appears in the record
  all_z_tbl <- AWI_CM_tbl %>% 
    group_by(AWI_CM_CELL) %>%
    summarize(all_z = all(AWI_CM_heff <= THRESHOLD))
  
  AWI_CM_locs_sf <- AWI_CM_locs_sf %>%
    left_join(all_z_tbl, by = "AWI_CM_CELL") 
  
  AWI_CM_tbl <- AWI_CM_locs_sf %>% 
    st_set_geometry(NULL) %>%
    left_join(AWI_CM_tbl, by = "AWI_CM_CELL") %>%
    mutate(AWI_CM_heff = ifelse(AWI_CM_heff > THRESHOLD, AWI_CM_heff, 0)) %>%
    filter(!is.na(region)) %>%
    filter(!all_z) %>%
    filter(!is.na(AWI_CM_heff)) %>%
    filter(region %in% 6:8)
  
  # return data and locations
  out <- list()
  out$AWI_CM_tbl <- AWI_CM_tbl 
  out$AWI_CM_locs_sf <- AWI_CM_locs_sf %>%
    dplyr::select(-area) %>%
    filter(AWI_CM_CELL %in% AWI_CM_tbl$AWI_CM_CELL)
  out
}

# LOAD TRAINING DATA AND COMPUTE CLIMATOLOGY------------------------------------
load_training_data <- 
  function(AWI_CM_tbl, AWI_CM_locs_sf, obs_contour_sf,
           init_tm, n_clim_trn_yr, n_trend_trn_mth, n_test_mth) {
    
    # test timepoints
    test_tms <- (init_tm + 1):(init_tm + n_test_mth)
    
    # trend training timepoints
    trend_trn_tms <- (init_tm - n_trend_trn_mth + 1):init_tm
    
    spde_model_tms <- c(trend_trn_tms, test_tms)
    
    # climatology training timepoints
    clim_trn_tms <- as.vector(sapply(spde_model_tms, 
                                     function(i) i - (1:n_clim_trn_yr) * 12))
    
    # estimate climatological average
    clim_trn_tbl <- AWI_CM_tbl %>%
      filter(AWI_CM_tm %in% clim_trn_tms) %>%
      group_by(AWI_CM_CELL, mth) %>%
      summarise(CLIM_heff = mean(AWI_CM_heff))
    
    trend_trn_tbl_list <- list()
    contour_poly_sf_list <- list()
    for (i in spde_model_tms) {
      print(i)
      current_mth <- i -  floor((i - 1) / 12) * 12
      current_yr <- 2007 + floor((i - 1)/ 12)
      
      AWI_CM_trend_tbl <-  AWI_CM_tbl %>%
        filter(AWI_CM_tm == i) %>%
        filter(!all_z) %>%
        left_join(clim_trn_tbl, by = c("AWI_CM_CELL", "mth")) %>% 
        mutate(anom = AWI_CM_heff - CLIM_heff,
               group_tm = i + 1 - min(spde_model_tms))
      AWI_CM_model_sf <- AWI_CM_trend_tbl %>% 
        dplyr::select(-region, -all_z) %>%
        left_join(AWI_CM_locs_sf, by = "AWI_CM_CELL") %>%
        st_as_sf() 
      
      contour_poly_sf <- obs_contour_sf %>%
        filter(AWI_CM_tm == i)
      
      contour_poly_sf_list <- c(contour_poly_sf_list, list(contour_poly_sf))
      
      in_obs_contour = st_intersects(AWI_CM_model_sf, contour_poly_sf, sparse = F)[,1]
      AWI_CM_model_sf <- AWI_CM_model_sf %>% 
        mutate(in_obs_contour = in_obs_contour)
      trend_trn_tbl_list <- c(trend_trn_tbl_list, list(AWI_CM_model_sf))
    }
    trend_trn_tbl <- do.call(rbind, trend_trn_tbl_list)
    out <- list()
    out$clim_trn_tbl <- clim_trn_tbl
    out$test_tbl <- trend_trn_tbl %>% filter(AWI_CM_tm > max(trend_trn_tms))
    out$trn_tbl <- trend_trn_tbl %>% filter(AWI_CM_tm <= max(trend_trn_tms))
    out$contour_poly_sf_list <- contour_poly_sf_list
    out
  }


# SETUP SPDE APPROXIMATION -----------------------------------------------------
setup_models <- function(trn_tbl, test_tbl = NULL) {
  spde_trn_tbl <- trn_tbl
  LC_ID_tbl <- NULL
  lcs <- NULL
  
  if (!is.null(test_tbl)) {
    spde_trn_tbl <- rbind(trn_tbl %>% 
                            mutate(in_pred_contour = T,
                                   on_pred_contour = F,
                                   traj_id = NA),
                          test_tbl %>%
                            mutate(anom = ifelse(in_pred_contour, NA, anom))) %>%
      mutate(group_tm = AWI_CM_tm - max(trn_tbl$AWI_CM_tm) + 1)
  }
  
  # mesh parameters
  ## WARNING: MAY WANT TO TWEAK THESE PARAMETERS
  obs_bnd = inla.nonconvex.hull(st_coordinates(spde_trn_tbl), 
                                convex = 80, resolution = c(150,150))
  mesh = inla.mesh.2d(boundary = obs_bnd, max.edge = c(MAX_EDGE, 200),
                      offset = c(200, 200))
  
  # stationary model
  s_spde <- inla.spde2.pcmatern(
    # Mesh and smoothness parameter
    mesh = mesh, alpha = 2,
    prior.range = c(2000, .5),
    prior.sigma = c(2, .01))
  
  n_group <- max(spde_trn_tbl$group_tm)
  
  #test_regs <- na.omit(unique(test_tbl$region))
  # space-time index
  i_x <- inla.spde.make.index(name = 'i_x', 
                              n.spde = s_spde$n.spde, 
                              n.group = n_group)
  
  # projection matrix
  A <- inla.spde.make.A(mesh, 
                        loc = st_coordinates(spde_trn_tbl),
                        group = spde_trn_tbl$group_tm)
  
  
  # data "stack" used by INLA
  stk <- inla.stack(
    data = list(y = spde_trn_tbl$anom),
    A = list(A), 
    effects = list(list(i_x = i_x$i_x, 
                        i_x.group = i_x$i_x.group)),
    tag = 'trn'
  )
  if (!is.null(test_tbl)) {
    test_tms <- unique(test_tbl$AWI_CM_tm)
    test_regs <- 6:8
    lc_row_list <- list()
    lc_id_list <- list()
    current <- 1
    for (i in test_tms) {
      for (reg in test_regs) {
        prop_lc <- as.vector(((spde_trn_tbl$region == reg & 
                                 spde_trn_tbl$AWI_CM_tm == i & 
                                 spde_trn_tbl$in_pred_contour) * 
                                spde_trn_tbl$area) %*% A)
        if (!all(prop_lc == 0)) {
          lc_row_list[[current]] <- prop_lc
          lc_id_list[[current]] <- 
            tibble(ID = current, AWI_CM_tm = i, region = reg)
          current <- current + 1
        }
      }
      ## ADD FULL VOLUME LC
      lc_row_list <- c(lc_row_list, 
                       list(as.vector(((spde_trn_tbl$AWI_CM_tm == i &
                                          spde_trn_tbl$in_pred_contour) * 
                                         spde_trn_tbl$area) %*% A)))
      lc_id_list <- c(lc_id_list, 
                      list(tibble(ID = current, AWI_CM_tm = i, region = 0)))
      current <- current + 1
    }
    
    lc_mat <- do.call(rbind, lc_row_list)
    LC_ID_tbl  <- bind_rows(lc_id_list)
    lcs <- inla.make.lincombs(i_x = lc_mat)
  }
  out <- list()
  out$LC_ID_tbl <- LC_ID_tbl
  out$lcs <- lcs
  out$s_spde <- s_spde
  out$i_x <- i_x # index for GMRF
  out$stk <- stk
  return(out)
  
}
# TRAIN MODELS BASED ON TREND TRAINING DATA ------------------------------------
fit_models <- function(AWI_CM_tbl, AWI_CM_locs_sf, init_tm, 
                       n_clim_trn_yr, n_trend_trn_mth,
                       n_test_mth, long_res = F) {
  
  last_time_pt <- Sys.time()
  start_time <- last_time_pt
  obs_contour_sf <- readRDS( paste0(contour_path,
                                    "observed_level_15_2011_2018_sf.rds"))
  split <- load_training_data(AWI_CM_tbl, AWI_CM_locs_sf, 
                              obs_contour_sf, init_tm,
                              n_clim_trn_yr, n_trend_trn_mth, n_test_mth)
  trn_tbl <- split$trn_tbl
  test_tbl <- split$test_tbl
  clim_trn_tbl <- split$clim_trn_tbl
  
  sub_trn_tbl <- trn_tbl %>% filter(in_obs_contour)
  trn_models <- setup_models(sub_trn_tbl)
  
  load_trn_time <- Sys.time() - last_time_pt
  last_time_pt <- Sys.time()
  
  # PRIORS ETC.
  # P(rho_T>.1)=.9
  rho_prior <- list(theta = list(prior = 'pccor1', param = c(0.1, 0.9)))
  sigma_w_prior <- list(prec = list(prior = "pc.prec", param = c(1, 0.01)))
  sigma_e_prior <- list(prec = list(prior = "pc.prec", param = c(1, 0.01)))
  # INLA SETTINGS
  ctr_group_w <- list(model = 'ar1', hyper = rho_prior)
  ctr_group_x <- list(model = 'ar1', hyper = rho_prior)
  
  # RUN INLA
  # stationary model
  ptm <- proc.time()
  S_trn_res <- 
    inla(y ~ 0 + 
           f(i_x, model = trn_models$s_spde,
             group = i_x.group, control.group = ctr_group_x),
         data = inla.stack.data(trn_models$stk), 
         control.family = list(hyper = sigma_e_prior),
         control.predictor = list(A = inla.stack.A(trn_models$stk), compute = T))
  
  S_trn_theta <- S_trn_res$internal.summary.hyperpar$mode
  
  S_trn_time <- Sys.time() - last_time_pt
  last_time_pt <- Sys.time()
  
  # non-spatial model
  NS_trn_res <- inla(anom ~ 0 + f(AWI_CM_CELL, model = 'iid', 
                                  hyper = sigma_w_prior,
                                  group = group_tm, control.group = ctr_group_w),
                     data = sub_trn_tbl, 
                     control.family = list(hyper = sigma_e_prior),
                     control.predictor = list(compute = T),
                     control.compute = list(config = TRUE))
  
  NS_trn_theta <- NS_trn_res$internal.summary.hyperpar$mode
  
  NS_trn_time <- Sys.time() - last_time_pt
  last_time_pt <- Sys.time()
  
  time_tbl <- tibble(load_train = load_trn_time,
                     S_train = S_trn_time,
                     NS_train = NS_trn_time)
  
  
  out <- list()
  out$time_tbl <- time_tbl
  out$NS_trn_res <- NS_trn_res
  out$NS_trn_theta <- NS_trn_theta
  out$S_trn_res <- S_trn_res
  out$S_trn_theta <- S_trn_theta
  out$trn_models <- trn_models
  out$trn_tbl <- trn_tbl
  out$clim_trn_tbl <- clim_trn_tbl
  out$test_tbl <- test_tbl
  out$AWI_CM_locs_sf <- AWI_CM_locs_sf
  out
}
# ADD FORECAST CONTOURS TO TEST_TBL OBJECT FOR FORECASTING ---------------------
augment_test_tbl <- function(test_tbl, pred_contour_sf, 
                             AWI_CM_tbl, AWI_CM_locs_sf, clim_trn_tbl,
                             ice_edge_coords = NULL) {
  test_tms <- unique(test_tbl$AWI_CM_tm)
  
  group_tms <- unique(test_tbl$group_tm)
  out_tbl_list <- list()
  for (i in 1:length(test_tms)) {
    
    current_tm <- test_tms[i]
    print(current_tm)
    current_mth <- current_tm -  floor((current_tm - 1) / 12) * 12
    current_yr <- 2007 + floor((current_tm - 1)/ 12)
    out_test_tbl <- test_tbl %>%
      filter(AWI_CM_tm == test_tms[i]) %>%
      mutate(on_pred_contour = FALSE,
             traj_id = NA)
    combined_test_sf <- out_test_tbl
    if (!is.null(ice_edge_coords)) {
      ice_edge_pts_sf <- ice_edge_coords %>%
        filter(AWI_CM_tm == test_tms[i]) %>%
        mutate(on_pred_contour = TRUE,
               in_obs_contour = FALSE,
               all_z = F,
               region = NA,
               area = 0,
               mth = current_mth,
               yr = current_yr,
               AWI_CM_CELL = NA,
               group_tm = group_tms[i],
               CMST_tm = AWI_CM_tm - 45)
      
      # interpolation of climatological average at contour pts
      # first, get climatology at all test locations (sp object)
      clim_all_points_sp <- AWI_CM_locs_sf %>%
        left_join(clim_trn_tbl %>%
                    filter(mth == current_mth), by = "AWI_CM_CELL") %>%
        as("Spatial")
      # set up formula for interpolation via gstat
      clim_fm <- gstat(formula = CLIM_heff ~ 1, data = clim_all_points_sp)
      ice_edge_pts_sp <- SpatialPoints(st_coordinates(ice_edge_pts_sf),
                                       proj4string = CRS(proj4string(clim_all_points_sp)))
      # add climatology column to contour pts sf
      ice_edge_pts_sf <- ice_edge_pts_sf %>%
        mutate(CLIM_heff = predict(clim_fm, ice_edge_pts_sp)$var1.pred,
               AWI_CM_heff = .15,
               anom = AWI_CM_heff - CLIM_heff) 
      combined_test_sf <- rbind(combined_test_sf, ice_edge_pts_sf)
    }
    current_pred_contour_sf <- pred_contour_sf %>% 
      filter(AWI_CM_tm == test_tms[i])
    
    in_pred_contour = st_intersects(combined_test_sf,
                                    current_pred_contour_sf, sparse = F)[,1]
    combined_test_sf <- combined_test_sf %>% 
      mutate(in_pred_contour = in_pred_contour & !on_pred_contour)
    out_tbl_list <- c(out_tbl_list, list(combined_test_sf))
  }
  out_test_tbl <- do.call(rbind, out_tbl_list)
}


# GENERATE FORECASTS CONDITIONAL ON CONTOURS -----------------------------------
get_forecasts <- function(trn_res, test_tm,
                          n_contours = NULL,
                          qs = c(.7, .8, .9, .95, .99)) {
  last_time_pt <- Sys.time()
  start_time <- last_time_pt
  trn_tbl <- trn_res$trn_tbl
  sub_trn_tbl <- trn_res$trn_tbl %>% 
    filter(in_obs_contour) %>%
    filter(AWI_CM_tm == max(trn_tbl$AWI_CM_tm)) 
  init_tm <- max(trn_tbl$AWI_CM_tm)
  init_mth <- init_tm -  floor((init_tm - 1) / 12) * 12
  init_yr <- 2007 + floor((init_tm - 1) / 12)
  mth_str <- str_pad(init_mth, 2, "left", "0")
  print(paste0(contour_path, "sample_contours_init_", 
               init_yr, "_", mth_str, ".rds"))
  pred_contour_sf <- readRDS(paste0(contour_path, "sample_contours_init_", 
                                    init_yr, "_", mth_str, ".rds"))
  if (CDN_ON_EDGE) {
    ice_edge_coords_sf <- readRDS(paste0(contour_path, "ice_edges_init_",
                                         init_yr, "_", mth_str, ".rds"))
  }
  test_tbl <- trn_res$test_tbl
  test_regs <- 6:8
  AWI_CM_locs_sf <- trn_res$AWI_CM_locs_sf
  clim_trn_tbl <- trn_res$clim_trn_tbl
  
  NS_marginal_pred_list <- list()
  S_marginal_pred_list <- list()
  NS_volume_pred_list <- list()
  S_volume_pred_list <- list()
  
  for (k in 1:n_contours) {
    
    current_tm <- test_tm
    current_mth <- current_tm -  floor((current_tm - 1) / 12) * 12
    current_yr <- 2007 + floor((current_tm - 1) / 12)
    current_contour_sf <- 
      pred_contour_sf %>%
      dplyr::filter(mth == current_mth) %>%
      dplyr::filter(traj_id == k) %>%
      mutate(AWI_CM_tm = current_tm)
    st_crs(current_contour_sf) <- st_crs(test_tbl)
    
    if (CDN_ON_EDGE) {
      test_edge_coords_sf <- ice_edge_coords_sf %>%
        dplyr::filter(mth == current_mth) %>%
        dplyr::filter(traj_id == k) %>% 
        mutate(AWI_CM_tm = current_tm) 
      st_crs(test_edge_coords_sf) <- st_crs(test_tbl)
      aug_test_tbl <- augment_test_tbl(test_tbl %>% filter(AWI_CM_tm == current_tm),
                                       current_contour_sf,
                                       AWI_CM_tbl, AWI_CM_locs_sf, clim_trn_tbl,
                                       ice_edge_coords = test_edge_coords_sf) 
    } else {
      aug_test_tbl <- augment_test_tbl(test_tbl %>% filter(AWI_CM_tm == current_tm),
                                       current_contour_sf,
                                       AWI_CM_tbl, AWI_CM_locs_sf, clim_trn_tbl,
                                       ice_edge_coords = NULL) 
    }
    test_models <- setup_models(sub_trn_tbl, 
                                aug_test_tbl %>% 
                                  filter(in_pred_contour | on_pred_contour))
    
    # PRIORS ETC.
    # P(rho_T>.1)=.9
    rho_prior <- list(theta = list(prior = 'pccor1', param = c(0.1, 0.9)))
    sigma_w_prior <- list(prec = list(prior = "pc.prec", param = c(1, 0.01)))
    sigma_e_prior <- list(prec = list(prior = "pc.prec", param = c(1, 0.01)))
    # INLA SETTINGS
    ctr_group_w <- list(model = 'ar1', hyper = rho_prior)
    ctr_group_x <- list(model = 'ar1', hyper = rho_prior)
    
    # GET POSTERIORS FOR SPATIAL MODEL, GIVEN SAMPLED CONTOUR k
    # PULL THETA FROM TRN_RES
    S_res <- 
      inla(y ~ 0 + 
             f(i_x, model = test_models$s_spde,
               group = i_x.group, control.group = ctr_group_x),
           data = inla.stack.data(test_models$stk), 
           lincomb = test_models$lcs,
           control.family = list(hyper = sigma_e_prior),
           control.predictor = list(A = inla.stack.A(test_models$stk), compute = T),
           control.compute=list(config = TRUE),
           control.mode = list(theta = trn_res$S_trn_theta, restart = F))
    
    NS_trn_tbl <- rbind(sub_trn_tbl%>%
                          mutate(in_pred_contour = T,
                                 on_pred_contour = F,
                                 traj_id = NA),
                        aug_test_tbl %>% 
                          filter(in_pred_contour) %>%
                          mutate(anom = NA)) %>%
      mutate(group_tm = AWI_CM_tm - max(sub_trn_tbl$AWI_CM_tm) + 1)
    lc_row_list <- list()
    lc_id_list <- list()
    current <- 1
    for (reg in test_regs) {
      prop_lc <-
        as.vector((NS_trn_tbl$region == reg & 
                     NS_trn_tbl$AWI_CM_tm == test_tm & 
                     NS_trn_tbl$in_pred_contour) * NS_trn_tbl$area)
      if (!all(prop_lc == 0)) {
        lc_row_list[[current]] <- prop_lc
        lc_id_list[[current]] <- 
          tibble(ID = current, AWI_CM_tm = test_tm, region = reg)
        current <- current + 1
      }
    }
    ## ADD FULL VOLUME LC
    lc_row_list <- c(lc_row_list, 
                     list(as.vector(((NS_trn_tbl$AWI_CM_tm == test_tm &
                                        NS_trn_tbl$in_pred_contour) * NS_trn_tbl$area))))
    lc_id_list <- c(lc_id_list, 
                    list(tibble(ID = current, AWI_CM_tm = test_tm, region = 0)))
    current <- current + 1
    
    
    lc_mat <- do.call(rbind, lc_row_list)
    LC_ID_tbl  <- bind_rows(lc_id_list)
    lcs <- inla.make.lincombs(AWI_CM_CELL = lc_mat)
    
    # GET POSTERIORS FOR NON-SPATIAL MODEL, GIVEN SAMPLED CONTOUR k
    # PULL THETA FROM TRN_RES
    NS_res <- inla(anom ~ 0 + f(AWI_CM_CELL, model = 'iid', 
                                hyper = sigma_w_prior,
                                group = group_tm, control.group = ctr_group_w),
                   data = NS_trn_tbl, 
                   lincomb = lcs,
                   control.family = list(hyper = sigma_e_prior),
                   control.predictor = list(compute = T),
                   control.compute = list(config = TRUE),
                   control.mode = list(theta = trn_res$NS_trn_theta, restart = F))
    NS_fc <- get_NS_forecasts(NS_res, trn_tbl, aug_test_tbl, LC_ID_tbl)
    S_fc <- get_S_forecasts(S_res, trn_tbl, aug_test_tbl,  LC_ID_tbl)
    
    NS_marginal_pred_list[[k]] <- NS_fc$marginal_pred_tbl 
    S_marginal_pred_list[[k]] <-  S_fc$marginal_pred_tbl 
    NS_volume_pred_list[[k]] <- NS_fc$volume_pred_tbl
    S_volume_pred_list[[k]]<- S_fc$volume_pred_tbl 
  }
  pred_time <- Sys.time() - last_time_pt
  last_time_pt <- Sys.time()
  print("finished with posteriors, now summarizing results")
  NS_pred_tbl_list <- list()
  S_pred_tbl_list <- list()
  NS_vol_pred_tbl_list <- list()
  S_vol_pred_tbl_list <- list()
  for (k in 1:n_contours) {
    NS_pred_tbl_list[[k]] <- test_tbl %>% 
      filter(AWI_CM_tm == current_tm) %>%
      dplyr::select(AWI_CM_CELL, AWI_CM_tm, CLIM_heff,
                    AWI_CM_heff, in_obs_contour) %>%
      st_set_geometry(NULL) %>%
      left_join(NS_marginal_pred_list[[k]], 
                by = c("AWI_CM_tm", "AWI_CM_CELL")) %>%
      mutate(pred_y_mean = ifelse(is.na(pred_y_mean), 0, pred_y_mean),
             pred_y_sd = ifelse(is.na(pred_y_sd), 0, pred_y_sd),
             contour = k)
    S_pred_tbl_list[[k]] <- test_tbl %>% 
      filter(AWI_CM_tm == current_tm) %>%
      dplyr::select(AWI_CM_CELL, AWI_CM_tm, CLIM_heff,
                    AWI_CM_heff, in_obs_contour) %>%
      st_set_geometry(NULL) %>%
      left_join(S_marginal_pred_list[[k]], 
                by = c("AWI_CM_tm", "AWI_CM_CELL")) %>%
      mutate(pred_y_mean = ifelse(is.na(pred_y_mean), 0, pred_y_mean),
             pred_y_sd = ifelse(is.na(pred_y_sd), 0, pred_y_sd),
             contour = k)
    NS_vol_pred_tbl_list[[k]] <- NS_volume_pred_list[[k]] %>%
      mutate(contour =  k)
    S_vol_pred_tbl_list[[k]] <- S_volume_pred_list[[k]] %>%
      mutate(contour =  k)
  }
  print("calculating marginal coverages")
  NS_marginal_cov_tbl_list <- list()
  S_marginal_cov_tbl_list <- list()
  NS_vol_cov_tbl_list <- list()
  S_vol_cov_tbl_list <- list()
  true_vol_tbl <- test_tbl %>% 
    st_set_geometry(NULL) %>%
    filter(region %in% 6:8) %>%
    group_by(region, AWI_CM_tm) %>%
    summarize(
      CLIM_vol = sum(area * CLIM_heff),
      AWI_CM_vol = sum(area * AWI_CM_heff)) %>%
    ungroup()
  true_full_vol_tbl <- test_tbl %>% 
    st_set_geometry(NULL) %>%
    filter(in_obs_contour) %>%
    mutate(region = 0) %>%
    group_by(region, AWI_CM_tm) %>%
    summarize(CLIM_vol = sum(area * CLIM_heff),
              AWI_CM_vol = sum(area * AWI_CM_heff)) %>%
    ungroup()
  true_vol_tbl <- bind_rows(true_vol_tbl,
                            true_full_vol_tbl)
  for (which_q in qs) {
    NS_cov_tbl <- bind_rows(NS_pred_tbl_list) %>%
      group_by(AWI_CM_CELL, AWI_CM_tm, AWI_CM_heff) %>%
      summarize(pred_y_mean = mean(pred_y_mean),
                ub = qnormmix(p = (2 * which_q - 1), wts = rep(1 / n()),
                              means = pred_y_mean, sds = pred_y_sd),
                lb = qnormmix(p = ((1 - which_q) / 2), wts = rep(1 / n()),
                              means = pred_y_mean, sds = pred_y_sd)) %>%
      ungroup() %>%
      left_join(test_tbl %>%
                  dplyr::select(AWI_CM_CELL, AWI_CM_tm, AWI_CM_heff,
                                CLIM_heff, in_obs_contour),
                by = c("AWI_CM_CELL", "AWI_CM_tm")) %>%
      mutate(q = which_q) %>%
      ungroup()
    NS_marginal_cov_tbl_list <- c(NS_marginal_cov_tbl_list, 
                                  list(NS_cov_tbl))
    NS_vol_cov_tbl <- bind_rows(NS_vol_pred_tbl_list) %>%
      group_by(region, AWI_CM_tm) %>%
      summarize(vol_mean = mean(vol_mean),
                ub = qnormmix(p = (2 * which_q - 1), wts = rep(1 / n()),
                              means = vol_mean, sds = vol_sd),
                lb = qnormmix(p = ((1 - which_q) / 2), wts = rep(1 / n()),
                              means = vol_mean, sds = vol_sd))%>%
      ungroup() %>%
      left_join(true_vol_tbl %>% # get true volume
                  dplyr::select(region, AWI_CM_tm, 
                                AWI_CM_vol, CLIM_vol),
                by = c("region", "AWI_CM_tm")) %>%
      mutate(q = which_q) %>%
      ungroup()
    NS_vol_cov_tbl_list <- c(NS_vol_cov_tbl_list,
                             list(NS_vol_cov_tbl))
    S_cov_tbl <- bind_rows(S_pred_tbl_list) %>%
      group_by(AWI_CM_CELL, AWI_CM_tm) %>%
      summarize(pred_y_mean = mean(pred_y_mean),
                ub = qnormmix(p = (2 * which_q - 1), wts = rep(1 / n()),
                              means = pred_y_mean, sds = pred_y_sd),
                lb = qnormmix(p = ((1 - which_q) / 2), wts = rep(1 / n()),
                              means = pred_y_mean, sds = pred_y_sd))%>%
      ungroup() %>%
      left_join(test_tbl %>%
                  dplyr::select(AWI_CM_CELL, AWI_CM_tm, AWI_CM_heff,
                                CLIM_heff, in_obs_contour),
                by = c("AWI_CM_CELL", "AWI_CM_tm")) %>%
      mutate(q = which_q) %>%
      ungroup()
    S_marginal_cov_tbl_list <- c(S_marginal_cov_tbl_list, 
                                 list(S_cov_tbl))
    S_vol_cov_tbl <- bind_rows(S_vol_pred_tbl_list) %>%
      group_by(region, AWI_CM_tm) %>%
      summarize(vol_mean = mean(vol_mean),
                ub = qnormmix(p = (2 * which_q - 1), wts = rep(1 / n()),
                              means = vol_mean, sds = vol_sd),
                lb = qnormmix(p = ((1 - which_q) / 2), wts = rep(1 / n()),
                              means = vol_mean, sds = vol_sd)) %>%
      ungroup() %>%
      left_join(true_vol_tbl %>% # get true volume
                  dplyr::select(region, AWI_CM_tm, 
                                AWI_CM_vol, CLIM_vol),
                by = c("region", "AWI_CM_tm")) %>%
      mutate(q = which_q) %>%
      ungroup()
    S_vol_cov_tbl_list <- c(S_vol_cov_tbl_list,
                            list(S_vol_cov_tbl))
  }
  NS_marginal_cov_tbl <- bind_rows(NS_marginal_cov_tbl_list) %>%
    mutate(ub = ifelse(ub < 0, 0, ub),
           lb = ifelse(lb < 0, 0, lb),
           pred_y_mean = ifelse(pred_y_mean < 0, 0, pred_y_mean),
           cov = (AWI_CM_heff <= ub) & (AWI_CM_heff >= lb)) %>%
    mutate(lead = AWI_CM_tm - max(trn_tbl$AWI_CM_tm)) 
  NS_vol_cov_tbl <- bind_rows(NS_vol_cov_tbl_list) %>%
    mutate(cov = (AWI_CM_vol <= ub) & (AWI_CM_vol >= lb)) %>%
    mutate(lead = AWI_CM_tm - max(trn_tbl$AWI_CM_tm))
  S_marginal_cov_tbl <- bind_rows(S_marginal_cov_tbl_list) %>%
    mutate(ub = ifelse(ub < 0, 0, ub),
           lb = ifelse(lb < 0, 0, lb),
           pred_y_mean = ifelse(pred_y_mean < 0, 0, pred_y_mean),
           cov = (AWI_CM_heff <= ub) & (AWI_CM_heff >= lb)) %>%
    mutate(lead = AWI_CM_tm - max(trn_tbl$AWI_CM_tm)) 
  S_vol_cov_tbl <- bind_rows(S_vol_cov_tbl_list) %>%
    mutate(cov = (AWI_CM_vol <= ub) & (AWI_CM_vol >= lb)) %>%
    mutate(lead = AWI_CM_tm - max(trn_tbl$AWI_CM_tm)) 
  print("generating summaries")
  CLIM_mse_tbl <- test_tbl %>% 
    filter(AWI_CM_tm == current_tm) %>%
    st_set_geometry(NULL) %>%
    dplyr::select(AWI_CM_CELL, AWI_CM_tm, CLIM_heff, 
                  AWI_CM_heff, in_obs_contour) %>%
    mutate(lead = AWI_CM_tm - max(trn_tbl$AWI_CM_tm)) %>%
    group_by(AWI_CM_tm, lead) %>%
    summarize(mse = mean((AWI_CM_heff - CLIM_heff)^2)) %>%
    mutate(model = "Climatology") %>%
    ungroup()
  CLIM_mae_tbl <- test_tbl %>% 
    filter(AWI_CM_tm == current_tm) %>%
    st_set_geometry(NULL) %>%
    dplyr::select(AWI_CM_CELL, AWI_CM_tm, CLIM_heff, 
                  AWI_CM_heff, in_obs_contour) %>%
    mutate(lead = AWI_CM_tm - max(trn_tbl$AWI_CM_tm)) %>%
    group_by(AWI_CM_tm, lead) %>%
    summarize(mae = mean(abs(AWI_CM_heff - CLIM_heff))) %>%
    mutate(model = "Climatology") %>%
    ungroup()
  CLIM_in_contour_mse_tbl <- test_tbl %>% 
    filter(AWI_CM_tm == current_tm) %>%
    st_set_geometry(NULL) %>%
    dplyr::select(AWI_CM_CELL, AWI_CM_tm, CLIM_heff,
                  AWI_CM_heff, in_obs_contour) %>%
    filter(in_obs_contour) %>%
    mutate(lead = AWI_CM_tm - max(trn_tbl$AWI_CM_tm)) %>%
    group_by(AWI_CM_tm, lead) %>%
    summarize(mse = mean((AWI_CM_heff - CLIM_heff)^2)) %>%
    mutate(model = "Climatology") %>%
    ungroup()
  CLIM_in_contour_mae_tbl <- test_tbl %>% 
    filter(AWI_CM_tm == current_tm) %>%
    st_set_geometry(NULL) %>%
    dplyr::select(AWI_CM_CELL, AWI_CM_tm, CLIM_heff, 
                  AWI_CM_heff, in_obs_contour) %>%
    filter(in_obs_contour) %>%
    mutate(lead = AWI_CM_tm - max(trn_tbl$AWI_CM_tm)) %>%
    group_by(AWI_CM_tm, lead) %>%
    summarize(mae = mean(abs(AWI_CM_heff - CLIM_heff))) %>%
    mutate(model = "Climatology") %>%
    ungroup()
  
  NS_summary <- summarize_preds(NS_marginal_cov_tbl %>%
                                  filter(in_obs_contour | pred_y_mean > 0),
                                "Non-spatial")
  S_summary <- summarize_preds(S_marginal_cov_tbl %>%
                                 filter(in_obs_contour | pred_y_mean > 0),
                               "Spatial")
  NS_in_contour_summary <- 
    summarize_preds(NS_marginal_cov_tbl %>% filter(in_obs_contour), "Non-spatial")
  S_in_contour_summary <- 
    summarize_preds(S_marginal_cov_tbl %>% filter(in_obs_contour), "Spatial")
  
  full_marginal_mse_tbl <- bind_rows(CLIM_mse_tbl,
                                     NS_summary$mse_tbl,
                                     S_summary$mse_tbl)
  full_marginal_mae_tbl <- bind_rows(CLIM_mae_tbl,
                                     NS_summary$mae_tbl,
                                     S_summary$mae_tbl)
  full_marginal_cov_tbl <- bind_rows(NS_summary$cov_tbl,
                                     S_summary$cov_tbl)
  in_contour_marginal_mse_tbl <- bind_rows(CLIM_in_contour_mse_tbl,
                                           NS_in_contour_summary$mse_tbl,
                                           S_in_contour_summary$mse_tbl)
  in_contour_marginal_mae_tbl <- bind_rows(CLIM_in_contour_mae_tbl,
                                           NS_in_contour_summary$mae_tbl,
                                           S_in_contour_summary$mae_tbl)
  in_contour_marginal_cov_tbl <- bind_rows(NS_in_contour_summary$cov_tbl,
                                           S_in_contour_summary$cov_tbl)
  vol_cov_tbl <- bind_rows(NS_vol_cov_tbl %>% mutate(model = "Non-spatial"),
                           S_vol_cov_tbl %>% mutate(model = "Spatial"))
  
  summary_time <- Sys.time() - last_time_pt
  last_time_pt <- Sys.time()
  time_tbl <- tibble(pred = pred_time,
                     summary = summary_time)
  out <- list()
  out$time_tbl <- time_tbl
  out$full_marginal_mse_tbl <- full_marginal_mse_tbl
  out$full_marginal_mae_tbl <- full_marginal_mae_tbl
  out$full_marginal_cov_tbl <- full_marginal_cov_tbl
  out$in_contour_marginal_mse_tbl <- in_contour_marginal_mse_tbl
  out$in_contour_marginal_mae_tbl <- in_contour_marginal_mae_tbl
  out$in_contour_marginal_cov_tbl <- in_contour_marginal_cov_tbl
  out$vol_cov_tbl <- vol_cov_tbl
  out$NS_marginal_cov_tbl <- NS_marginal_cov_tbl
  out$S_marginal_cov_tbl <- S_marginal_cov_tbl
  return(out)
}

# NONSPATIAL MODEL FORECASTS ---------------------------------------------------
# get forecasts from non-spatial model (based on INLA marginal posteriors)
get_NS_forecasts <- function(NS_res, trn_tbl, aug_test_tbl, LC_ID_tbl) {
  sub_trn_tbl <- trn_tbl %>%
    filter(in_obs_contour) %>%
    filter(AWI_CM_tm == max(trn_tbl$AWI_CM_tm))
  sub_test_tbl <- aug_test_tbl %>% 
    filter(in_pred_contour) %>%
    mutate(anom = NA)
  test_ind <- (nrow(sub_trn_tbl) + 1):(nrow(sub_trn_tbl) + nrow(sub_test_tbl))
  
  # posterior means of hyperparameters
  est_sigma_e <- sqrt(1 / NS_res$summary.hyperpar$mean[1])
  
  pred_w_sd <- NS_res$summary.linear.predictor$sd[test_ind]
  pred_y_sd <- sqrt(pred_w_sd^2 + est_sigma_e^2)
  
  marginal_pred_tbl <- sub_test_tbl %>% 
    st_set_geometry(NULL) %>%
    dplyr::select(AWI_CM_CELL, AWI_CM_tm, in_pred_contour,
                  on_pred_contour, CLIM_heff) %>%
    mutate(pred_y_mean = CLIM_heff + 
             NS_res$summary.linear.predictor$mean[test_ind],
           pred_y_sd = pred_y_sd)
  marginal_pred_tbl <- marginal_pred_tbl %>%
    filter(!on_pred_contour) %>%
    dplyr::select(-c(in_pred_contour, on_pred_contour, CLIM_heff))
  
  volume_pred_tbl <- sub_test_tbl %>% 
    st_set_geometry(NULL) %>%
    mutate(var_err = est_sigma_e ^ 2) %>%
    filter(region %in% 6:8) %>%
    group_by(region, AWI_CM_tm) %>%
    summarize(var_err = sum(var_err)) %>%
    ungroup()
  
  full_volume_pred_tbl <- sub_test_tbl %>% 
    st_set_geometry(NULL) %>%
    mutate(var_err =  est_sigma_e ^ 2) %>%
    mutate(region = 0) %>%
    group_by(region, AWI_CM_tm) %>%
    summarize(var_err = sum(var_err)) %>%
    ungroup()
  volume_pred_tbl <- bind_rows(volume_pred_tbl,
                               full_volume_pred_tbl)
  
  volume_pred_tbl <- volume_pred_tbl %>%
    left_join(LC_ID_tbl, by = c("region", "AWI_CM_tm")) %>%
    left_join(NS_res$summary.lincomb.derived %>%
                dplyr::select(ID, mean, sd), by = "ID") %>%
    mutate(vol_mean = CLIM_vol + mean,
           vol_sd = sqrt(var_err + sd ^ 2)) %>%
    dplyr::select(-c(var_err, mean, sd, ID))
  
  out <- list()
  out$marginal_pred_tbl <- marginal_pred_tbl
  out$volume_pred_tbl <- volume_pred_tbl
  out
  
}

# SPATIAL STATIONARY MODEL FORECASTS -------------------------------------------
get_S_forecasts <- function(S_res, trn_tbl, aug_test_tbl, LC_ID_tbl) {
  sub_trn_tbl <- trn_tbl %>%
    filter(in_obs_contour) %>%
    filter(AWI_CM_tm == max(trn_tbl$AWI_CM_tm))
  sub_test_tbl <- aug_test_tbl %>% 
    filter(in_pred_contour | on_pred_contour) %>%
    mutate(anom = NA)
  
  est_sigma_e <- sqrt(1 / S_res$summary.hyperpar$mean[1])
  est_rho_x <- S_res$summary.hyperpar$mean[4]
  est_sigma_z <- S_res$summary.hyperpar$mean[3]
  
  test_ind <- (nrow(sub_trn_tbl) + 1):(nrow(sub_trn_tbl) + nrow(sub_test_tbl))
  pred_wz_sd <- S_res$summary.linear.predictor$sd[test_ind]
  pred_y_sd <- sqrt(pred_wz_sd^2 + est_sigma_e^2)
  marginal_pred_tbl <- sub_test_tbl %>% 
    st_set_geometry(NULL) %>%
    dplyr::select(AWI_CM_CELL, AWI_CM_tm, in_pred_contour, 
                  on_pred_contour, CLIM_heff) %>%
    mutate(pred_y_mean = CLIM_heff + 
             S_res$summary.linear.predictor$mean[test_ind],
           pred_y_sd = pred_y_sd)
  
  marginal_pred_tbl <- marginal_pred_tbl %>%
    filter(!on_pred_contour) %>%
    dplyr::select(-c(in_pred_contour, on_pred_contour, CLIM_heff))
  
  volume_pred_tbl <- sub_test_tbl %>% 
    st_set_geometry(NULL) %>%
    mutate(var_err = est_sigma_e ^ 2) %>%
    filter(region %in% 6:8) %>%
    group_by(region, AWI_CM_tm) %>%
    summarize(var_err = sum(var_err)) %>%
    ungroup()
  
  full_volume_pred_tbl <- sub_test_tbl %>% 
    st_set_geometry(NULL) %>%
    mutate(var_err =  est_sigma_e ^ 2) %>%
    mutate(region = 0)%>%
    group_by(region, AWI_CM_tm) %>%
    summarize(var_err = sum(var_err)) %>%
    ungroup()
  
  volume_pred_tbl <- bind_rows(volume_pred_tbl,
                               full_volume_pred_tbl)
  
  volume_pred_tbl <- volume_pred_tbl %>%
    left_join(LC_ID_tbl, by = c("region", "AWI_CM_tm")) %>%
    left_join(S_res$summary.lincomb.derived %>%
                dplyr::select(ID, mean, sd), by = "ID") %>%
    mutate(vol_mean = CLIM_vol + mean,
           vol_sd = sqrt(var_err + sd ^ 2)) %>%
    dplyr::select(-c(var_err, mean, sd, ID))
  
  out <- list()
  out$marginal_pred_tbl <- marginal_pred_tbl
  out$volume_pred_tbl <- volume_pred_tbl
  out
  
}

# summarize coverage rates across locations
summarize_preds <- function(cov_tbl, model_name) {
  cov_summary_tbl <- cov_tbl %>% 
    group_by(q, AWI_CM_tm, lead) %>%
    summarize(cov_rate = mean(cov))  %>%
    mutate(model = model_name) %>%
    ungroup()
  mse_tbl <- cov_tbl %>%
    group_by(AWI_CM_tm, lead) %>%
    summarize(mse = mean((AWI_CM_heff - pred_y_mean)^2)) %>%
    mutate(model = model_name) %>%
    ungroup()
  mae_tbl <- cov_tbl %>%
    group_by(AWI_CM_tm, lead) %>%
    summarize(mae = mean(abs(AWI_CM_heff - pred_y_mean))) %>%
    mutate(model = model_name) %>%
    ungroup()
  summary <- list()
  summary$cov_tbl <- cov_summary_tbl
  summary$mse_tbl <- mse_tbl
  summary$mae_tbl <- mae_tbl
  summary
}


get_summary <- function(trn_res, forecast_obj) {
  out_summary <- forecast_obj$summary
  
  out_summary$S_pred_tbl <- forecast_obj$S_marginal_cov_tbl %>%
    filter(q == .95)
  out_summary$NS_pred_tbl <- forecast_obj$NS_marginal_cov_tbl %>%
    filter(q == .95)
  
  out_summary$summary_hyperpar_NS <- trn_res$NS_trn_res$summary.hyperpar
  out_summary$summary_hyperpar_S <- trn_res$S_trn_res$summary.hyperpar
  
  out_summary$vol_cov_tbl <- forecast_obj$vol_cov_tbl
  out_summary$full_marginal_mse_tbl <- 
    forecast_obj$full_marginal_mse_tbl
  out_summary$full_marginal_mae_tbl <- 
    forecast_obj$full_marginal_mae_tbl  
  out_summary$full_marginal_cov_tbl <- 
    forecast_obj$full_marginal_cov_tbl
  out_summary$in_contour_marginal_mse_tbl <-
    forecast_obj$in_contour_marginal_mse_tbl
  out_summary$in_contour_marginal_mae_tbl <-
    forecast_obj$in_contour_marginal_mae_tbl
  out_summary$in_contour_marginal_cov_tbl <- 
    forecast_obj$in_contour_marginal_cov_tbl
  return(out_summary)
}