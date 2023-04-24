################################################################################
###       Visualizing Planet Mosaics (before and after pre-processing)       ###
###                        Code by HGR 1/2023                                ###
################################################################################

### Load Libraries #############################################################
library(terra)
library(sf)
library(magick)
library(RStoolbox)
library(terrainr)
library(ggthemes)
library(tidyverse)
################################################################################

### Load Data ##################################################################
# Sentinel Data
s2_0 <- rast('~/Documents/permafrost_pathways/rts_mapping/planet_processing_test/data/sentinel_data/yamal_gydan_polygons/Sentinel-2_YG_Mosaic_bgrnir_UTM42N_0.tif')
# prior to pre-processing
pre_cal_files_0 <- list.files('~/Documents/permafrost_pathways/rts_mapping/planet_processing_test/data/planet_data/yamal_gydan_polygons/polygon_id_0',
                        pattern = 'SR_clip\\.tif$',
                        full.names = TRUE,
                        recursive = TRUE)
pre_cal_images_0 <- map(pre_cal_files_0,
                  ~ rast(.x))
names(pre_cal_images_0) <- map(pre_cal_images_0,
                             ~ .x@ptr$get_sourcenames())
# post pre-processing
cal_files_0 <- list.files('~/Documents/permafrost_pathways/rts_mapping/planet_processing_test/data/automated_download/mad_output/calibrated_files',
                            pattern = '^0_.+SR_clip_arosics_mad\\.tif$',
                            full.names = TRUE,
                            recursive = TRUE)
cal_images_0 <- map(cal_files_0,
                      ~ rast(.x))
names(cal_images_0) <- map(cal_images_0,
                           ~ .x@ptr$get_sourcenames())
# CRS
local_crs <- crs(cal_images_0[[1]])
################################################################################

### Functions ##################################################################
# To convert raster to scaled dataframe
rast_to_df <- function(img, scale_values = NULL, band_names) {
  filepath <- img@ptr$filenames()
  img_id <- str_extract(img@ptr$get_sourcenames(), '[:digit:]{8}_[:digit:]{6}_([:digit:]{2}_)?[:alnum:]{4}')
  print(img_id)
  n_bands <- nlyr(img)
  names(img) <- band_names[1:n_bands]
  print(filepath)
  
  df <- img %>%
    as.data.frame(xy = TRUE)
  
  if (is.null(scale_values)) {
    scale_values <- c()
    for (band in 1:nlyr(img)) {
      scale_values <- c(scale_values, 
                        min(select(df, band_names[band]), na.rm = TRUE),
                        max(select(df, band_names[band]), na.rm = TRUE))
    }
    scale_values = matrix(scale_values, byrow = TRUE, ncol = 2)
  }
  print(scale_values)
  
  if (dim(scale_values)[1] == 1) {
    scale_values <- matrix(rep(scale_values, n_bands), ncol = 2, byrow = TRUE)
    warning('The set of scale values you provided are being used across all bands.')
  } else if (dim(scale_values)[1] > n_bands) {
    warning('More sets of scale values given than number of bands. Ignoring extra values.')
  } else if (dim(scale_values)[1] < n_bands) {
    stop('The number of scale value sets is greater than 1 and less than the number of bands.')
  }
  
  cal_status <- str_extract(filepath, 'calibrated|planet_data')
  cal_status <- case_when(cal_status == 'calibrated' ~ 'Calibrated',
                          cal_status == 'planet_data' ~ 'Uncalibrated')
  
  df <- df %>%
    mutate(across(c(3:ncol(.)), 
                  ~ case_when(
                    .x < scale_values[match(cur_column(), band_names), ][1] ~ 
                      0,
                    .x > scale_values[match(cur_column(), band_names), ][2] ~ 
                      1,
                    TRUE ~ (.x - scale_values[match(cur_column(), band_names), ][1])/
                      (scale_values[match(cur_column(), band_names), ][2] - scale_values[match(cur_column(), band_names), ][1])
                    )
                  )) %>%
    mutate(img = img_id,
           cal = factor(cal_status, levels = c('Uncalibrated', 'Calibrated')))
  
  return(df)
}
################################################################################

### Format Data for plotting ###################################################
bands = c('blue', 'green', 'red', 'nir')
scale_values <- matrix(c(280, 415, 470, 0, 750, 930, 1180, 2000), 
                       ncol = 2)

pid_0_df <- map_dfr(append(pre_cal_images_0[c('20170804_063123_1006_3B_AnalyticMS_SR_clip',
                                              '20180719_044304_104a_3B_AnalyticMS_SR_clip',
                                              '20180721_064457_103b_3B_AnalyticMS_SR_clip',
                                              '20210720_043317_32_106e_3B_AnalyticMS_SR_clip')],
                           cal_images_0[c('0_20180719_044304_104a_3B_AnalyticMS_SR_clip_arosics_mad',
                                          '0_20180721_064457_103b_3B_AnalyticMS_SR_clip_arosics_mad',
                                          '0_20210720_043317_32_106e_3B_AnalyticMS_SR_clip_arosics_mad')]),
                    ~ rast_to_df(.x, scale_values, bands))
################################################################################

### Faceted Figure #############################################################
# ggplot() +
#   geom_spatial_rgb(data = filter(yg_mosaics_pre_df, mosaic == '2017'),
#                    aes(x = x, y = y, r = red, g = green, b = blue),
#                    scale = 1) +
#   theme_bw() +
#   theme(axis.title = element_blank(),
#         axis.text.x = element_text(angle = 90, vjust = 0.5)) +
#   coord_sf(crs = crs(yg_mosaic_pre_2017),
#            datum = crs(yg_mosaic_pre_2017),
#            expand = FALSE)

pid_0_fig <- ggplot() +
  geom_spatial_rgb(data = pid_0_df,
                   aes(x = x, y = y, r = red, g = green, b = blue),
                   scale = 1) +
  facet_grid(cal ~ img) +
  theme_bw() +
  theme(axis.title = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5)) +
  coord_sf(crs = local_crs,
           datum = local_crs,
           expand = FALSE)
pid_0_fig
# ggsave('figures/pid_0_examples.jpg',
#        pid_0_fig,
#        height = 5,
#        width = 9,
#        bg = 'white')
################################################################################
