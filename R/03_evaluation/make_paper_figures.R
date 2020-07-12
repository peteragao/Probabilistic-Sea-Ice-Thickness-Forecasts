# Probabilistic-Sea-Ice-Thickness-Forecasts/R/03_evaluation/make_paper_figures.R
#
# Author: Peter Gao (petergao@uw.edu) 
#
# License: (c) Peter Gao 2020, MIT
#
# This model takes generated forecasts and summaries and creates plots for
# the paper
#

project_home <- "~/Dropbox/Probabilistic-Sea-Ice-Thickness-Forecasts/" 
data_path <- paste0(project_home, "data/") # STORE DATA HERE
ECMWF_path <- paste0(data_path, "clean_ECMWF/") # cleaned ECMWF files 
AWI_CM_path <- paste0(data_path, "clean_AWI_CM/") # cleaned AWI-CM files 

model <- "uncdn" # default choice for plots
res_path <- paste0(project_home, "results/02_cdn_thickness/", model, "/")
contour_path <- paste0(project_home, "output/01_contour")
tab_path <- paste0(project_home, "output/paper/tables/")
fig_path <- paste0(project_home, "output/paper/figures/")
output_path <- paste0(project_home, "output/")
# LIBRARIES --------------------------------------------------------------------
print("Loading libraries...")
library(dplyr)
library(tidyr)
library(sf)
library(gstat)
library(raster)
library(sp)
library(ggplot2)
library(rasterVis)
library(INLA)
library(RColorBrewer)
library(gstat)
library(latex2exp)
library(xtable)
library(stars)
library(viridis)
library(stringr)
#inla.setOption(pardiso.license = "~/sea-ice/analysis/pardiso.lic")
INLA:::inla.dynload.workaround()
#inla.pardiso.check()

# ANALYSIS CONSTANTS -----------------------------------------------------------
THRESHOLD <- .15 # below this, set thickness = 0
PS_crs <- paste("+init=epsg:3413 +units=km +proj=stere +lat_0=90 +lat_ts=70",
                "+lon_0=-45 +k=1 +x_0=0 +y_0=0 +datum=WGS84 +no_defs",
                "+ellps=WGS84 +towgs84=0,0,0 )")
PS_ra <- raster(nrows = 448, ncols = 304, 
                xmn = -3850, xmx = 3750,
                ymn = -5350, ymx = 5850,
                crs = PS_crs)
load_polar_map <- function(crs = CRS("+init=epsg:3413 +units=km")) {
  data("wrld_simpl", package = "maptools")      
  wm <- crop(wrld_simpl, extent(-180, 180, 0, 90)) 
  wm.proj <- spTransform(wm, crs)
  return(wm.proj)
}
polar_map <- st_as_sf(load_polar_map(crs = PS_crs))

# PLOT CONSTANTS ---------------------------------------------------------------
plt_tms <- 97:144
plt_yr <- 2007 + floor((plt_tms - 1)/ 12)
plt_mth <- plt_tms - floor((plt_tms - 1) / 12) * 12
plt_tm_labels <- paste0(plt_mth, "/", plt_yr)
max_tm <- max(plt_tms)
n_tm <- length(plt_tms)

quartzFonts(helvetica = c("Helvetica Neue", "Helvetica Neue Bold", 
                          "Helvetica Neue Light Italic", 
                          "Helvetica Neue Bold Italic"))
par(family = "helvetica")
col_tbl <- tibble(source = c("Adj. non-spatial", "Adj. spatial", "AWI-CM", "CLIM", 
                             "CMST", "ECMWF", "Nominal", "Non-spatial",
                             "Spatio-temporal+", "Spatio-temporal"),
                  cols = c( "#D55E00", "#CC79A7","black", "gray70", "#0072B2",
                            "#009E73","black", "#E69F00", "#F0E442", "#56B4E9"))

my_theme <- theme_classic() + 
  theme(strip.background = element_rect(colour=NA, fill="#e6f1ff"),
        title =element_text(size=8, face='bold'),
        legend.text = element_text(size = 8),
        axis.title.x =element_text(size= 8, face='bold'), 
        axis.text.x= element_text(size= 8),
        axis.title.y =element_text(size= 8, face='bold'), 
        axis.text.y= element_text(size= 8))
# LOAD FUNCTIONS ---------------------------------------------------------------
source(paste0(project_home, "R/02_cdn_thickness/thickness_model_functions.R"))

# LOAD SUMMARY -----------------------------------------------------------------
summary <- readRDS(paste0(res_path, "summary.rds"))

# FIGURE 1: example SIT maps ---------------------------------------------------
test <- load_AWI_CM_data()
AWI_CM_locs_sf <- test$AWI_CM_locs_sf
AWI_CM_tbl <- test$AWI_CM_tbl

# which months to plot (internal time index)
map_tms <- c(85,88, 91, 94)+12
mth <- map_tms - floor((map_tms - 1) / 12) * 12
mth_str <- str_pad(mth, 2, "left", "0")
yr <- 2007 + floor((map_tms - 1)/ 12) 
out_ra_list <- list()

# load relevant data
for (tm in map_tms) {
  in_file <- 
    readRDS(paste0(project_home, "results/02_cdn_thickness/uncdn/summary/", 
                   tm, ".rds"))$S_pred_tbl %>%
    mutate(AWI_CM_heff = ifelse(in_obs_contour, AWI_CM_heff, 0)) 
  ex_map_sf <- AWI_CM_locs_sf %>%
    left_join(in_file %>% dplyr::filter(AWI_CM_tm == tm + 1)) 
  ex_map_ra <- rasterize(ex_map_sf  %>% dplyr::select(AWI_CM_heff), PS_ra, fun = mean)
  ex_map_ra <- raster::aggregate(ex_map_ra[[2]], 2)
  out_ra_list <- c(out_ra_list, list(ex_map_ra))
}
ex_map_ra <- stack(out_ra_list)
ex_map_stars <- st_as_stars(ex_map_ra) %>%
  st_set_dimensions(3, values = paste0(yr, "/", mth_str))
# load contours
obs_contour_sf <- readRDS(paste0(contour_path, 
                                 "observed_level_15_2011_2018_sf.rds"))
# generate plot
pdf(paste0(fig_path, "example_SIT_maps.pdf"), width = 7.5, height = 2.5)
ggplot()+geom_stars(data = ex_map_stars)+
  facet_grid(~band) + 
  geom_sf(data = polar_map, color = "gray40", 
          fill = "gray50", size = .25,na.rm = T) +
  coord_sf(xlim = c(-3200, 2800), ylim = c(-4100,2500), expand = F) + 
  scale_fill_gradientn(name = "Thickness (m.)", guide = "colourbar", 
                       limits = c(0,6),
                       colors =  c("white", viridis(6, direction = -1)),
                       values = c(0, seq(.01, 1, length.out = 6)), 
                       na.value = "gray90") +
  my_theme + xlab("") + ylab("")
dev.off()

# FIGURE 2: region map ---------------------------------------------------------

# load locations and region info
AWI_CM_locs_PS_sf <- readRDS(paste0(AWI_CM_path, "AWI_CM_locs_locs_PS_sf.rds"))
reg_key <- readRDS(paste0(data_path, "clean_NSIDC/regions_key.rds"))

plt_sf <- AWI_CM_locs_PS_sf  %>% filter(region %in% 6:8)
reg_labels <- reg_key$reg_str[match(unique(plt_sf$region), reg_key$region)]
plt_sf <- plt_sf %>%
  mutate(region = factor(reg_key$reg_str[match(region, reg_key$region)])) 
plt_ra <- rasterize(plt_sf  %>% dplyr::select(region), PS_ra, fun='last')
plt_ra  <- raster::aggregate(plt_ra[[2]], 2, fun = min)
plt_ra  <-ratify(plt_ra)

# convert to stars object
plt_stars <- st_as_stars(plt_ra) 
plt_stars$region <- ifelse(is.na(plt_stars$region), NA, plt_stars$region)
plt_stars$region <- as.factor(plt_stars$region)

pdf(paste0(fig_path, "region_map.pdf"), width = 3.5, height = 2.25)
ggplot(plt_sf) +geom_stars(data = plt_stars) + 
  geom_sf(data = polar_map, color = "gray40", 
          fill = "gray50", size = .15, na.rm = T) +
  scale_color_discrete(name = "region")+
  guides(color = guide_legend(override.aes = list(size=5))) + 
  scale_fill_manual(name = "Region", na.value = "gray90",
                    values = c("#1b9e77", "#d95f02", "#7570b3", "gray90"),
                    labels = c(reg_key$reg_str[7:9], "Other ocean"),
                    breaks = c(6:8, NA)) +
  coord_sf(xlim = c(-3200, 2800), ylim = c(-4100,2500), expand = F) + 
  my_theme + xlab("") + ylab("") 

# FIGURE 2: contour schematic --------------------------------------------------
r <- 1
start_pts <- reg_info$start_coords[[r]]

# generate fake data
y_max <- sapply(reg_info$sec_lengths[[r]], sum)
y_gen <- prop_to_y(prop = c(rep(1, 15), 
                            ilogit(rnorm(length(y_max) - 35,
                                         mean = 2, sd = .3)),
                            rep(1, 20)),  
                   r, reg_info, inds = 1:length(y_max))
angs_r <- reg_info$angs[[r]]
new_pts <- reg_info$start_coords[[r]] + 
  cbind(y_gen*cos(angs_r), y_gen*sin(angs_r))

# maximum end points
max_pts <- reg_info$start_coords[[r]] + 
  cbind(y_max*cos(angs_r), y_max*sin(angs_r))
all_pts <- rbind(new_pts)

# fixed points
B <- start_pts %>% 
  as.data.frame() %>%
  st_as_sf(coords = c(1,2))
st_crs(B) <- st_crs(PS_wm)

M <- new_pts %>% 
  as.data.frame() %>%
  st_as_sf(coords = c(1,2))
st_crs(M) <- st_crs(PS_wm)

# loop last point to first point to construct polygon
Lobs <- st_sf(geo = st_sfc(st_linestring(rbind(start_pts[35,],
                                               new_pts[35,]))))
st_crs(Lobs) <- st_crs(PS_wm)
Lmax <- st_sf(geo = st_sfc(st_linestring(rbind(start_pts[37,],
                                               max_pts[37,]))))
st_crs(Lmax) <- st_crs(PS_wm)

# create line segments (max length)
Lm <- list()
for (i in 1:nrow(new_pts)) {
  Lm[[i]] <- st_linestring(rbind(start_pts[i,], max_pts[i,]))
}

Lm <- st_sf(geo = st_sfc(st_multilinestring(Lm)))
st_crs(Lm) <- st_crs(PS_wm)

# create line segments (observed length)
L <- list()
for (i in 1:nrow(new_pts)) {
  L[[i]] <- st_linestring(rbind(start_pts[i,] , new_pts[i,]))
}
L <- st_sf(geo = st_sfc(st_multilinestring(L)))
st_crs(L) <- st_crs(PS_wm)

# create polygon (observed)
S <- rbind(all_pts, all_pts[1,])  %>%
  st_linestring() %>%
  st_cast("MULTILINESTRING")
S <- st_sf(geo = st_sfc(S))
st_crs(S) <- st_crs(PS_wm)
S_poly <- st_polygonize(S)

# create polygon (maximum)
Sm <- rbind(max_pts, max_pts[1,])  %>%
  st_linestring() %>%
  st_cast("MULTILINESTRING")
Sm <- st_sf(geo = st_sfc(Sm))
Sm_poly <- st_polygonize(Sm)
st_crs(Sm_poly) <- st_crs(PS_wm)

# plot for central Arctic
p0 <- ggplot() + 
  geom_sf(data = Sm_poly, aes(fill = 'outside region'), col = NA)+ 
  geom_sf(data = Sm_poly, aes(fill = 'not ice'), col = NA)+ 
  geom_sf(data = S_poly, aes(fill = 'ice')) + 
  geom_sf(data = PS_wm, col = NA) + 
  geom_sf(data = Lm, col = "gray80", linetype = 'dashed', size = .3) + 
  geom_sf(data = L, col = "gray60", size = .5) + 
  geom_sf(data = S, aes(col = 'S')) + 
  geom_sf(data = Lobs, aes(col = 'Y'), size = .9) + 
  geom_sf(data = Lmax, aes(col = 'M'), size = .9)+
  geom_sf(data = PS_wm, col = NA, fill = 'gray50') +
  geom_sf(data = B, aes(col = 'B')) + 
  theme_classic() + theme(panel.background = element_rect(fill = "gray90")) + 
  xlim(-2800, 2100) + ylim(-1300,2400) + 
  scale_colour_discrete(name = "Legend",
                        breaks=c("S", "Y", "M", "B"),
                        labels=c(unname(TeX("\\mathbf{c}")), 
                                 unname(TeX("\\mathbf{y_{35}}")), 
                                 unname(TeX("\\mathbf{M_{37}}")),
                                 unname(TeX("\\mathbf{b}"))))+
  scale_fill_manual(name = "",
                    breaks=c("ice", "not ice", "outside region"),
                    values=c("#cffffd","white", "gray90"))+
  ggtitle("Central Arctic")

# plot for Greenland
r <- 5
start_pts <- reg_info$start_coords[[r]]
y_max <- sapply(reg_info$sec_lengths[[r]], sum)
y_gen <- prop_to_y(prop = ilogit(sin(1:length(y_max)/length(y_max))),  
                   r, reg_info, inds = 1:length(y_max))
angs_r <- reg_info$angs[[r]]
new_pts <- reg_info$start_coords[[r]] + 
  cbind(y_gen*cos(angs_r), y_gen*sin(angs_r))
max_pts <- reg_info$start_coords[[r]] + 
  cbind(y_max*cos(angs_r), y_max*sin(angs_r))
all_pts <- rbind(start_pts, new_pts[nrow(new_pts):1,])
B <- start_pts %>% 
  as.data.frame() %>%
  st_as_sf(coords = c(1,2))
st_crs(B) <- st_crs(PS_wm)
M <- new_pts %>% 
  as.data.frame() %>%
  st_as_sf(coords = c(1,2))
st_crs(M) <- st_crs(PS_wm)

Lobs <-
  st_sf(geo = st_sfc(st_linestring(rbind(start_pts[21,], new_pts[21,]))))
st_crs(Lobs) <- st_crs(PS_wm)
Lmax <- 
  st_sf(geo = st_sfc(st_linestring(rbind(start_pts[20,], max_pts[20,]))))
st_crs(Lmax) <- st_crs(PS_wm)

# create line segments (max length)
Lm <- list()
for (i in 1:nrow(new_pts)) {
  Lm[[i]] <- st_linestring(rbind(start_pts[i,] , max_pts[i,]))
}

Lm <- st_sf(geo = st_sfc(st_multilinestring(Lm)))
st_crs(Lm) <- st_crs(PS_wm)

# create line segments (observed length)
L <- list()
for (i in 1:nrow(new_pts)) {
  L[[i]] <- st_linestring(rbind(start_pts[i,] , new_pts[i,]))
}
L <- st_sf(geo = st_sfc(st_multilinestring(L)))
st_crs(L) <- st_crs(PS_wm)

# create polygon (observed)
S <- rbind(all_pts, all_pts[1, ])  %>%
  st_linestring() %>%
  st_cast("MULTILINESTRING")
S <- st_sf(geo = st_sfc(S))
st_crs(S) <- st_crs(PS_wm)
S_poly <- st_polygonize(S)

# create polygon (maximum)
Sm <- rbind(start_pts, max_pts[nrow(new_pts):1,], start_pts[1,])  %>%
  st_linestring() %>%
  st_cast("MULTILINESTRING")
Sm <- st_sf(geo = st_sfc(Sm))
Sm_poly <- st_polygonize(Sm)
st_crs(Sm_poly) <- st_crs(PS_wm)
p1 <- ggplot() +
  geom_sf(data = Sm_poly, aes(fill = 'outside region'), col = NA)+ 
  geom_sf(data = Sm_poly, aes(fill = 'not ice'), col = NA)+ 
  geom_sf(data = S_poly, aes(fill = 'ice')) + 
  geom_sf(data = PS_wm, col = NA) + 
  geom_sf(data = Lm, col = "gray80", linetype = 'dashed', size = .3) + 
  geom_sf(data = L, col = "gray60", size = .5) + 
  geom_sf(data = S, aes(col = 'S')) + 
  geom_sf(data = Lobs, aes(col = 'Y'), size = .9) + 
  geom_sf(data = Lmax, aes(col = 'M'), size = .9)+
  geom_sf(data = PS_wm, col = NA, fill = 'gray50') +
  geom_sf(data = B, aes(col = 'B')) + 
  theme_classic() + theme(panel.background = element_rect(fill = "gray90")) + 
  scale_colour_discrete(name = "Legend",
                        breaks=c("S", "Y", "M", "B"),
                        labels=c(unname(TeX("\\mathbf{c}")), 
                                 unname(TeX("\\mathbf{y_{21}}")), 
                                 unname(TeX("\\mathbf{M_{20}}")),
                                 unname(TeX("\\mathbf{b}")))) +
  scale_fill_manual(name = "",
                    breaks=c("ice", "not ice", "outside region"),
                    values=c("#cffffd","white", "gray90"))+
  ggtitle("Greenland Sea") + coord_sf(xlim = c(0, 2000), ylim = c(-3550,-100)) +
  scale_x_continuous(breaks = 0:50)
ggsave(paste0(fig_path, "paper/contour_schematic.pdf"), p0+p1)

# FIGURE 3: example contour preds ----------------------------------------------
# Example contour predictions and errors
init_tm <- 105
test_tm <- init_tm + 1
mth <- init_tm- floor((init_tm - 1) / 12) * 12
mth_str <- str_pad(mth, 2, "left", "0")
yr <- 2007 + floor((init_tm - 1)/ 12) 
test_mth <- test_tm - floor((test_tm - 1) / 12) * 12
test_mth_str <- str_pad(test_mth, 2, "left", "0")
test_yr <- 2007 + floor((test_tm - 1)/ 12) 
# observed contours
obs_contour_sf <- readRDS(paste0(contour_path, 
                                 "observed_level_15_2011_2018_sf.rds"))
# predicted contours
pred_contour_sf <- readRDS(paste0(contour_path, "sample_contours_init_",
                                  yr, "_", mth_str, ".rds"))

cont_prob_list <- readRDS(paste0(contour_path, "cont_prob_init_",
                                 test_yr, "_", mth_str, ".rds"))
cont_prob <- cont_prob_list[[2]]
cont_prob_ra <- flip(raster(t(cont_prob), template = PS_ra),2)
cont_prob_stars <- st_as_stars(cont_prob_ra) 
st_crs(pred_contour_sf) <- st_crs(obs_contour_sf)
ice_edge_sf <- readRDS(paste0(contour_path, "ice_edges_init_",
                              yr, "_", mth_str, ".rds"))
st_crs(ice_edge_sf) <- st_crs(obs_contour_sf) 

# sample contours
g_sample <- ggplot() + 
  geom_sf(data = polar_map, color = "gray40",
          fill = "gray50", size = .25,na.rm = T) + 
  geom_sf(data = pred_contour_sf %>%
            filter(traj_id %in% c(2,3,5)) %>% 
            filter(mth == test_mth), fill = NA, col = "tomato") + 
  geom_sf(data = ice_edge_sf %>%
            filter(traj_id %in% c(2,3,5)) %>% 
            filter(mth == test_mth), col= '#7570b3', size = .8) + 
  facet_wrap(~traj_id, nrow = 1) + 
  ggtitle(paste0("Sample contour forecasts (",
                 month.name[test_mth], " ", yr, ")")) +
  coord_sf(xlim = c(-3200, 2800), ylim = c(-4100,2500), expand = F) + 
  my_theme + 
  theme(
    strip.background = element_blank(),
    strip.text.x = element_blank()
  )

# Probability of ice map
g_poi <- ggplot() + 
   geom_sf(data = polar_map, color = "gray40", 
           fill = "gray50", size = .25,na.rm = T) + 
  geom_stars(data = cont_prob_stars) + 
  geom_sf(data = obs_contour_sf %>%
            filter(AWI_CM_tm == test_tm), fill = NA, col = "dodgerblue") + 
  scale_fill_gradientn(name = "Probability of ice", 
                       colors = c("white",
                                  viridis(6, option = "C", direction = -1)), 
                       guide = "legend") + 
  ylab("") + xlab("") +
  ggtitle(paste0("Probability of Ice (", month.name[test_mth], " ", yr, ")")) +
  coord_sf(xlim = c(-3200, 2800), ylim = c(-4100,2500), expand = F) + my_theme

pdf(paste0(fig_path, "example_contour_preds.pdf"), width = 6, height = 2.5)
g_sample
dev.off()
pdf(paste0(fig_path, "example_contour_poi.pdf"), width = 2.5, height = 2.5)
g_poi
dev.off()
# FIGURE 4: example thickness predictions --------------------------------------
# Example predictions and errors
init_tm <- 105
ld <- 1
test_tm <- init_tm + ld
mth <- init_tm- floor((init_tm - 1) / 12) * 12
mth_str <- str_pad(mth, 2, "left", "0")
yr <- 2007 + floor((init_tm - 1)/ 12) 
test_mth <- test_tm - floor((test_tm - 1) / 12) * 12
test_mth_str <- str_pad(test_mth, 2, "left", "0")
# load results
# area sf
area_sf <- readRDS(paste0(AWI_CM_path, "AWI_CM_area_sf.rds")) %>%
  mutate(area = area / 1e9) 

# tibble of AWI_CM estimate locations
in_file <- readRDS(paste0(project_home, 
                          "results/02_cdn_thickness/uncdn/summary/", 
                          tm, ".rds"))
cdn_file <- readRDS(paste0(project_home, 
                           "results/02_cdn_thickness/cdn_on_edge/summary/",
                           tm, ".rds"))
ECMWF_tbl <- readRDS(paste0(AWI_CM_path, "ECMWF_AWI_CM_tbl_lead_1.rds")) %>%
  filter(AWI_CM_tm == test_tm) %>%
  rename("ECMWF_mean" = ECMWF_heff) %>%
  mutate("ECMWF_err" = AWI_CM_heff - ECMWF_mean) %>%
  select(AWI_CM_CELL, "ECMWF_mean", "ECMWF_err")
AWI_CM_locs_sf <-
  readRDS(paste0(AWI_CM_path, "AWI_CM_locs_PS_sf.rds")) %>%
  left_join(area_sf %>% st_set_geometry(NULL), by = "AWI_CM_CELL") %>%
  filter(AWI_CM_CELL %in% in_file$S_pred_tbl$AWI_CM_CELL)
NS_pred_tbl <- in_file$NS_pred_tbl %>% filter(lead == ld) %>%
  mutate(pred_y_mean = ifelse(in_obs_contour, pred_y_mean, 0)) %>%
  rename("NS_uncdn_mean" = pred_y_mean) %>%
  mutate("NS_uncdn_len" = ub - lb,
         "NS_uncdn_err" = AWI_CM_heff - NS_uncdn_mean) %>%
  select(AWI_CM_CELL, "NS_uncdn_mean", "NS_uncdn_len", "NS_uncdn_err")
S_cdn_tbl <- in_file$S_pred_tbl %>% filter(lead == ld) %>%
  mutate(pred_y_mean = ifelse(in_obs_contour, pred_y_mean, 0)) %>%
  mutate(CLIM_heff = ifelse(in_obs_contour, CLIM_heff, 0)) %>%
  rename("S_cdn_mean" = pred_y_mean) %>%
  rename("CLIM_mean" = CLIM_heff) %>%
  mutate("CLIM_err" = AWI_CM_heff - CLIM_mean) %>%
  mutate("S_cdn_len" = ub - lb,
         "S_cdn_err" = AWI_CM_heff - S_cdn_mean) %>%
  select(AWI_CM_CELL,"CLIM_mean" , "CLIM_err" , 
         "S_cdn_mean", "S_cdn_len", "S_cdn_err")
S_pred_sf <- AWI_CM_locs_sf %>% 
  left_join(in_file$S_pred_tbl %>% filter(lead == 1)) %>%
  mutate(pred_y_mean = ifelse(in_obs_contour, pred_y_mean, 0)) %>%
  mutate(AWI_CM_heff = ifelse(in_obs_contour, AWI_CM_heff, 0)) %>%
  rename("S_uncdn_mean" = pred_y_mean) %>%
  rename("AWI_CM" = AWI_CM_heff) %>%
  mutate("S_uncdn_len" = ub - lb,
         "S_uncdn_err" = AWI_CM - S_uncdn_mean) %>%
  select(AWI_CM_CELL, "S_uncdn_mean", "S_uncdn_len",
         "S_uncdn_err", "AWI_CM", in_obs_contour) %>%
  left_join(NS_pred_tbl) %>%
  left_join(S_cdn_tbl) %>% 
  left_join(ECMWF_tbl) %>%
  mutate(ECMWF_mean = ifelse(in_obs_contour, ECMWF_mean, 0),
         ECMWF_err = ifelse(in_obs_contour, ECMWF_err, 0),
         S_cdn_err = ifelse(in_obs_contour, S_cdn_err, 0),
         CLIM_err = ifelse(in_obs_contour, CLIM_err, 0),
         NS_uncdn_err = ifelse(in_obs_contour,  NS_uncdn_err, 0),
         S_uncdn_err = ifelse(in_obs_contour,  S_uncdn_err, 0))
pred_ra <- rasterize(S_pred_sf %>%
                       dplyr::select( "NS_uncdn_mean", "S_uncdn_mean",
                                      "S_cdn_mean","ECMWF_mean", "CLIM_mean" ,"AWI_CM"),
                     PS_ra, fun = mean)
pred_ra <- raster::aggregate(pred_ra[[2:7]], 2)
pred_stars <- st_as_stars(pred_ra) %>%
  st_set_dimensions(3, values = c("Non-spatial", "Spatio-temporal", 
                                  "Spatio-temporal+", "ECMWF", "Climatology","AWI-CM"))
g_pred <- ggplot()+geom_stars(data = pred_stars)+
  facet_wrap(~factor(band)) + 
  geom_sf(data = polar_map, color = "gray40", 
          fill = "gray50", size = .25, na.rm = T) +
  coord_sf(xlim = c(-3200, 2800), ylim = c(-4100,2500), expand = F) + 
  scale_fill_gradientn(name = "Thickness (m.)", guide = "colourbar", 
                       limits = c(0,6),
                       colors =  c("white", viridis(6, direction = -1)),
                       values = c(0, seq(.01, 1, length.out = 6)), na.value = "gray90") +
  my_theme + xlab("") + ylab("")+ 
  theme(strip.text = element_text(size = 7),
        legend.key.size = unit(.4, "cm"))
error_ra <- rasterize(S_pred_sf %>% dplyr::select("NS_uncdn_err",
                                                  "S_uncdn_err", "S_cdn_err",
                                                  "ECMWF_err", "CLIM_err"),
                      PS_ra, fun = mean)
error_ra <- raster::aggregate(error_ra[[2:6]], 2)
error_stars <- st_as_stars(error_ra) %>%
  st_set_dimensions(3, values = c("Non-spatial", "Spatio-temporal", 
                                  "Spatio-temporal+", "ECMWF", "Climatology"))

g_error <- ggplot()+geom_stars(data = error_stars)+
  facet_wrap(~factor(band)) + 
  geom_sf(data = polar_map, color = "gray40", fill = "gray50", size = .25,na.rm = T) +
  coord_sf(xlim = c(-3200, 2800), ylim = c(-4100,2500), expand = F) + 
  scale_fill_gradient2(name = "Error (m.)    ", guide = "colorbar",
                       low = "#a50026",  mid = "white",  high = "#313695",
                       midpoint = 0, limits = c(-2,2),
                       na.value = "gray90") +
  my_theme + xlab("") + ylab("")+ 
  theme(strip.text = element_text(size = 7),
        legend.key.size = unit(.4, "cm"))


len_ra <- rasterize(S_pred_sf %>% dplyr::select("NS_uncdn_len",
                                                "S_uncdn_len", "S_cdn_len"),
                    PS_ra, fun = mean)
len_ra <- raster::aggregate(len_ra[[2:4]], 2)
len_stars <- st_as_stars(len_ra) %>%
  st_set_dimensions(3, values = c("Non-spatial CI length", 
                                  "Spatial CI length", "Spatial cdn CI length"))
g_len <- ggplot()+geom_stars(data = len_stars)+
  facet_grid(~band) + 
  geom_sf(data = polar_map, color = "gray40", fill = "gray50", size = .25,na.rm = T) +
  coord_sf(xlim = c(-3200, 2800), ylim = c(-4100,2500), expand = F) + 
  scale_fill_gradientn(name = "CI length (m.)",
                       colors =  c("white", "#f7fcf0", "#084081"),
                       values = c(0, .05,  1), na.value = "gray90") +  
  my_theme + xlab("") + ylab("") + 
  theme(strip.text = element_text(size = 7),
        legend.key.size = unit(.4, "cm"))
pdf(paste0(fig_path, "example_SIT_preds.pdf"), width = 7, height = 5)
g_pred
dev.off()
pdf(paste0(fig_path, "example_SIT_errs.pdf"), width = 7, height = 5)
g_error
dev.off()
pdf(paste0(fig_path, "example_SIT_CI_lengths.pdf"), width = 7, height = 2.5)
g_len
dev.off()
# TABLE 1: marginal coverage rates -----------------------------------------------
summary_uncdn <- readRDS(paste0(res_path, "summary.rds"))

summary_cdn <- readRDS(paste0(project_home, 
                              "results/02_cdn_thickness/cdn_on_edge/summary.rds"))

st_plus_mse <- summary_cdn$mse_tbl %>%
  filter(AWI_CM_tm > min(plt_tms) + 3 & AWI_CM_tm < max(plt_tms)) %>%
  filter(model == "Spatial") %>%
  mutate(model = "Spatio-temporal+")

mse_tbl <- summary_uncdn$mse_tbl %>% 
  filter(AWI_CM_tm > min(plt_tms) + 3 & AWI_CM_tm < max(plt_tms)) %>%
  mutate(model = recode(as.factor(model),"ECMWF_heff" = "ECMWF",
                        "Spatial" = "Spatio-temporal")) %>%
  bind_rows(st_plus_mse)

ECMWF_cov <- readRDS(paste0(output_path,
                            "ECMWF_marginal_cov_tbl.rds")) %>%
  filter(q == .9) %>%
  dplyr::select(-q) %>%
  rename(cov_rate = cov)  %>%
  filter(lead < 4)

st_plus_cov <- summary_cdn$cov_tbl %>%
  filter(AWI_CM_tm > min(plt_tms) + 3 & AWI_CM_tm < max(plt_tms)) %>%
  filter(q == 0.9) %>%
  filter(model == "Spatial") %>%
  mutate(model = "Spatio-temporal+")
cov_tbl <- summary$cov_tbl %>% 
  filter(AWI_CM_tm > min(plt_tms) + 3 & AWI_CM_tm < max(plt_tms)) %>%
  mutate(model = recode(as.factor(model),"ECMWF_heff" = "ECMWF",
                        "Spatial" = "Spatio-temporal")) %>%
  filter(q == 0.9) %>%
  bind_rows(st_plus_cov) %>%
  group_by(lead, model) %>%
  summarize(cov_rate = sqrt(mean(cov_rate))) %>%
  bind_rows(ECMWF_cov)

out_tbl <- mse_tbl %>%
  group_by(lead, model) %>%
  summarize(rmse = sqrt(mean(mse))) %>%
  left_join(cov_tbl) %>% 
  arrange(model) %>%
  select(model, everything())
kable(out_tbl, format = 'latex', digits = c(1, 2, 3, 3)) %>%
  collapse_rows(columns = 1)  %>%
  kable_styling(font_size = 7) 
write(out_kable, file = paste0(tab_path, "marginal_cov_tbl.tex"))

# TABLE 2: volume coverage rates -----------------------------------------------
# Volume coverage
reg_key <- readRDS(paste0(data_path, "clean_NSIDC/regions_key.rds"))
reg_key$reg_code[1] <- "ALL"
reg_key$reg_str[1] <- "All selected regions"
reg_labels <-
  reg_key$reg_str[match(unique(summary$vol_cov_tbl$region), reg_key$region)]
ECMWF_vol_cov <- readRDS(paste0(output_path,
                                "ECMWF_vol_cov_tbl.rds")) %>%
  group_by(region, lead) %>%
  summarize(cov_rate = mean(cov)) %>%
  ungroup() %>%
  mutate(model = "ECMWF") 

reg_labels <-
  reg_key$reg_str[match(unique(ECMWF_vol_cov$region), reg_key$region)]
ECMWF_vol_cov <- ECMWF_vol_cov %>%
  mutate(region = factor(region, labels = reg_labels)) %>%
  filter(lead < 4)
reg_labels <- 
  reg_key$reg_str[match(unique(summary$vol_cov_tbl$region), reg_key$region)]

out_tbl_plus <- summary_cdn$vol_cov_tbl %>%
  ungroup() %>%
  mutate(region = factor(region, labels = reg_labels)) %>%
  filter(q == .9) %>%
  filter(model == "Spatial") %>%
  mutate(model = "Spatio-temporal+") %>%
  group_by(region, lead, model) %>%
  summarize(cov_rate = mean(cov_rate)) %>%
  ungroup() %>%
  pivot_wider(names_from = region, values_from = cov_rate) 


out_tbl <- summary$vol_cov_tbl %>%
  ungroup() %>%
  mutate(region = factor(region, labels = reg_labels)) %>%
  filter(q == .9) %>%
  group_by(region, lead, model) %>%
  summarize(cov_rate = mean(cov_rate)) %>%
  ungroup() %>%
  bind_rows(ECMWF_vol_cov) %>%
  pivot_wider(names_from = region, values_from = cov_rate) %>%
  bind_rows(out_tbl_plus) %>% 
  mutate(model = recode(as.factor(model),
                        "Spatial" = "Spatio-temporal")) %>%
  arrange(model) %>%
  select(model, everything())

out_kable <- kable(out_tbl, format = 'latex', digits = c(1, 2, 2, 2, 2,2)) %>%
  collapse_rows(columns = 1)  %>%
  kable_styling(font_size = 7) 
write(out_kable, file = paste0(tab_path, "vol_comp_tbl.tex"))
# TABLE 2: volume coverage rates -----------------------------------------------

