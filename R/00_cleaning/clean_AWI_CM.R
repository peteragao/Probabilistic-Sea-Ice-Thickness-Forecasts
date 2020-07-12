# Probabilistic-Sea-Ice-Thickness-Forecasts/R/00_cleaning/clean_AWI_CM.R
#
# Author: Peter Gao (petergao@uw.edu) 
#
# License: (c) Peter Gao 2020, MIT
#
# This script loads monthly AWI_CM reanalysis data published online by
# Mu et al. 2020 (10.5281/zenodo.3597428) in the file
# Exp_CTD_T_SIC_SIT_SIV.tar.bz
# 
# Information on the AWI-CM grid is loaded from fesom.initial.mesh.diag.nc
# and nod2d.out
#
# Output:
#   - AWI_CM_locs_sf: sf object of the AWI-CM observation locations
#   - AWI_CM_locs_PS_sf: sf object of the projected AWI-CM observation locations
#   - area_sf:  sf object of the areas of AWI-CM cells
#   - AWI_CM_tbl: tibble of AWI-CM observations (for each year)
#
# AWI_CM_locs_sf ----------------------------------------------------------------

load_AWI_CM_locs_sf <- function(path = paste0(in_AWI_CM_path, "nod2d.out")) {
  coords <- read.table(path, skip = 1)[, 2:3] %>%
    rename(lon = V2, lat = V3)
  # parameters taken from pyfesom docs at 
  # https://buildmedia.readthedocs.org/media/pdf/pyfesom/latest/pyfesom.pdf
  coords <- sl.rot(lon = coords$lon,
                   lat = coords$lat,
                   alpha=50,beta=15,gamma=-90, invert = T)
  locs_tbl <- tibble(lon = coords$lon, 
                     lat = coords$lat,
                     AWI_CM_CELL = 1:length(coords$lat))
  return(st_as_sf(locs_tbl, coords = c("lon", "lat"),
                  crs = CRS("+init=EPSG:4326")))
}

AWI_CM_locs_sf <- load_AWI_CM_locs_sf()
saveRDS(AWI_CM_locs_sf, paste0(out_AWI_CM_path, "AWI_CM_locs_sf.rds"))
print(paste0("saved ", out_AWI_CM_path, "AWI_CM_locs_sf.rds"))

# Add regions and project to PS ------------------------------------------------
regions_raster <- readRDS(paste0(out_NSIDC_path, "regions_raster.rds"))
AWI_CM_locs_PS_sf <- AWI_CM_locs_sf %>%
  st_transform(crs = PS_crs) 
AWI_CM_locs_PS_sf <- AWI_CM_locs_PS_sf %>%
  mutate(region = raster::extract(regions_raster,
                                  as_Spatial(AWI_CM_locs_PS_sf)))
saveRDS(AWI_CM_locs_PS_sf, paste0(out_AWI_CM_path, "AWI_CM_locs_PS_sf.rds"))
print(paste0("saved ", out_AWI_CM_path, "AWI_CM_locs_PS_sf.rds"))

# AWI_CM_area_sf ---------------------------------------------------------------
ncin <- nc_open(paste0(in_AWI_CM_path, "fesom.initial.mesh.diag.nc"))
area = ncvar_get(ncin, "cluster_area")
area_tbl <- tibble(area = area,
                   AWI_CM_CELL = 1:length(area))
area_sf <- AWI_CM_locs_sf %>%
  left_join(area_tbl, by = "AWI_CM_CELL")
saveRDS(area_sf, paste0(out_AWI_CM_path, "AWI_CM_area_sf.rds"))
print(paste0("saved ", out_AWI_CM_path, "AWI_CM_area_sf.rds"))

# annual versions of AWI_CM_tbl  -----------------------------------------------

load_AWI_CM_yr <- function(yr) {
  ncin <- nc_open(paste0(in_AWI_CM_path,
                         "Exp_CTD_T_", yr, "_SIC_SIT_SIV.nc"))
  heff <- ncvar_get(ncin, "hice_a")
  AWI_CM_tbl <- tibble(AWI_CM_heff = as.vector(heff),
                       yr = yr,
                       mth = rep(1:12, each = nrow(heff)),
                       AWI_CM_CELL = rep(1:nrow(heff), 12))
  return(AWI_CM_tbl)
}

AWI_CM_locs_sf <- load_AWI_CM_locs_sf()
lats <- st_coordinates(AWI_CM_locs_sf)[,2]
N_HEMI_CELL <- AWI_CM_locs_sf$AWI_CM_CELL[which(lats > 0)]
N_AWI_CM_locs_sf <- AWI_CM_locs_sf[N_HEMI_CELL,]
for (i in 2007:2018) {
  print(i)
  AWI_CM_tbl <- load_AWI_CM_yr(i) %>%
    mutate(CMST_tm = (yr - 2010) * 12 + mth - 9) %>%
    mutate(AWI_CM_tm = (yr - 2007) * 12 + mth) %>%
    filter(AWI_CM_CELL %in% N_HEMI_CELL)
  saveRDS(AWI_CM_tbl, paste0(out_AWI_CM_path, "AWI_CM_tbl_", i, ".rds"))
  print(paste0("saved ", out_AWI_CM_path, "AWI_CM_tbl_", i, ".rds"))
}


