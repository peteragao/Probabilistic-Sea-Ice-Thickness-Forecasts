# Probabilistic-Sea-Ice-Thickness-Forecasts/R/00_cleaning/clean_NSIDC.R
#
# Author: Peter Gao (petergao@uw.edu) 
#
# License: (c) Peter Gao 2020, MIT
# This script loads NSIDC region information and stores it in raster format
#
# REGION IDS -------------------------------------------------------------------
# generate region labels on NSIDC's 25km Polar Stereographic grid
# as of 2020/04/21, file info available at
# https://nsidc.org/data/polar-stereo/tools_masks.html#region_masks

rasterize_region_masks <- function(verbose = F) {
  # binary file obtained from NSIDC
  to_read <- file(paste0(in_NSIDC_path, "region_n.msk"), 'rb')
  in_bin <- readBin(to_read, integer(), n = 448 * 304 + 300, size = 1)
  in_bin <- in_bin[301:length(in_bin)]
  regions <- matrix(in_bin, nrow = 448, ncol = 304, byrow = T)
  
  
  # NSIDC 25km polar stereographic grid
  regions_raster <- raster(regions, crs = PS_crs,
                           xmn=-3850, xmx=3750,
                           ymn=-5350, ymx=5850)
  
  out <- regions_raster
  # optional: print out regions key
  if (verbose) {
    cat("1: Non-regional ocean\n",
        "2: Sea of Okhotsk and Japan\n",	
        "3: Bering Sea\n",
        "4: Hudson Bay\n",
        "5: Baffin Bay/Davis Strait/Labrador Sea\n",
        "6: Greenland Sea\n",
        "7: Kara and Barents Seas\n",
        "8: Arctic Ocean\n",
        "9: Canadian Archipelago\n",
        "10: Gulf of St. Lawrence\n",
        "11: Land\n",
        "12: Coast\n",
        "0: Lakes, extended coast\n", sep = '')
  }
  close(to_read)
  return(out)
}
regions_raster <- rasterize_region_masks()
saveRDS(regions_key, paste0(out_NSIDC_path, "regions_key.rds"))
print(paste0("saved ", out_NSIDC_path, "regions_raster.rds"))
      
# key for region codes
regions_key <- tibble(region = 0:12,
                      reg_code = c("LAK", "NRO", "OKH", "BER", 
                                   "HUD", "BAF", "GRE", "BAR", 
                                   "CEN", "CAA", "GSL", "LAN", "COA"),
                      reg_str = c("Lakes", "Non-regional ocean",
                                  "Okhotsk/Japan", "Bering Sea", 
                                  "Hudson Bay","Baffin Bay",
                                  "Greenland Sea", "Kara/Barents", 
                                  "Central Arctic", "Canadian Arch.",
                                  "St. Lawrence","Land", "Coast"))

saveRDS(regions_key, paste0(out_NSIDC_path, "regions_key.rds"))
print(paste0("saved ", out_NSIDC_path, "regions_key.rds"))
