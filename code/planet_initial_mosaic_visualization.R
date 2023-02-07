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
s2 <- rast('~/Documents/permafrost_pathways/rts_mapping/planet_processing_test/data/sentinel_data/Sentinel-2_2021_Yamal_GS_Mosaic_bgrnir_UTM42N.tif')
# prior to pre-processing
yg_mosaic_pre_2017 <- rast('/home/hrodenhizer/Documents/permafrost_pathways/rts_mapping/planet_processing_test/data/unprocessed_composite/output/Planet_Median_Composite_Uncalibrated_2017.tif')
yg_mosaic_pre_2018 <- rast('/home/hrodenhizer/Documents/permafrost_pathways/rts_mapping/planet_processing_test/data/unprocessed_composite/output/Planet_Median_Composite_Uncalibrated_2018.tif')
yg_mosaic_pre_2019 <- rast('/home/hrodenhizer/Documents/permafrost_pathways/rts_mapping/planet_processing_test/data/unprocessed_composite/output/Planet_Median_Composite_Uncalibrated_2019.tif')
yg_mosaic_pre_2020 <- rast('/home/hrodenhizer/Documents/permafrost_pathways/rts_mapping/planet_processing_test/data/unprocessed_composite/output/Planet_Median_Composite_Uncalibrated_2020.tif')
yg_mosaic_pre_2021 <- rast('/home/hrodenhizer/Documents/permafrost_pathways/rts_mapping/planet_processing_test/data/unprocessed_composite/output/Planet_Median_Composite_Uncalibrated_2021.tif')
yg_mosaic_pre_all <- rast('/home/hrodenhizer/Documents/permafrost_pathways/rts_mapping/planet_processing_test/data/unprocessed_composite/output/Planet_Median_Composite_Uncalibrated_All.tif')
# post pre-processing
yg_mosaic_2017 <- rast('/home/hrodenhizer/Documents/permafrost_pathways/rts_mapping/planet_processing_test/data/calibrated_composite/output/Planet_Median_Composite_Calibrated_2017.tif')
yg_mosaic_2018 <- rast('/home/hrodenhizer/Documents/permafrost_pathways/rts_mapping/planet_processing_test/data/calibrated_composite/output/Planet_Median_Composite_Calibrated_2018.tif')
yg_mosaic_2019 <- rast('/home/hrodenhizer/Documents/permafrost_pathways/rts_mapping/planet_processing_test/data/calibrated_composite/output/Planet_Median_Composite_Calibrated_2019.tif')
yg_mosaic_2020 <- rast('/home/hrodenhizer/Documents/permafrost_pathways/rts_mapping/planet_processing_test/data/calibrated_composite/output/Planet_Median_Composite_Calibrated_2020.tif')
yg_mosaic_2021 <- rast('/home/hrodenhizer/Documents/permafrost_pathways/rts_mapping/planet_processing_test/data/calibrated_composite/output/Planet_Median_Composite_Calibrated_2021.tif')
yg_mosaic_all <- rast('/home/hrodenhizer/Documents/permafrost_pathways/rts_mapping/planet_processing_test/data/calibrated_composite/output/Planet_Median_Composite_Calibrated_All.tif')
# subset area for zooming in
zoom_area <-  st_sf(
  st_sfc(
    st_polygon(list(cbind(
      c(477000, 477000, 479000, 479000, 477000),
      c(7724000, 7725000, 7725000, 7724000, 7724000)
    ))),
    crs = crs(yg_mosaic_pre_2017)))
################################################################################

### Functions ##################################################################
# To convert raster to scaled dataframe
rast_to_df <- function(img, scale_values = NULL, band_names) {
  filepath <- img@ptr$filenames()
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
  
  if (dim(scale_values)[1] == n_bands) {
  } else if (dim(scale_values)[1] == 1) {
    scale_values <- matrix(rep(scale_values, n_bands), ncol = 2, byrow = TRUE)
  } else if (dim(scale_values)[1] > n_bands) {
    warning('More sets of scale values given than number of bands. Ignoring extra values.')
  } else if (dim(scale_values)[1] < n_bands) {
    stop('Incorrect number of scale values for the number of bands')
  }
  
  year <- str_extract(filepath, '[:digit:]{4}|All')
  cal_status <- str_extract(filepath, 'Calibrated|Uncalibrated')
  
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
    mutate(mosaic = year,
           cal = cal_status)
  
  return(df)
}
################################################################################

### Format Data for plotting ###################################################
yg_mosaics_pre <- list(yg_mosaic_pre_2017,
                       yg_mosaic_pre_2018,
                       yg_mosaic_pre_2019,
                       yg_mosaic_pre_2020,
                       yg_mosaic_pre_2021,
                       yg_mosaic_pre_all)
bands = c('blue', 'green', 'red')
scale_values <- matrix(c(200, 100, 100, 580, 1100, 1000), 
                       ncol = 2)

yg_mosaics_pre_df <- map_dfr(yg_mosaics_pre,
                         ~ rast_to_df(.x, scale_values, bands))

yg_mosaics_cal <- list(yg_mosaic_2017,
                       yg_mosaic_2018,
                       yg_mosaic_2019,
                       yg_mosaic_2020,
                       yg_mosaic_2021,
                       yg_mosaic_all)
bands = c('blue', 'green', 'red', 'nir')
scale_values <- matrix(c(200, 100, 100, 0, 580, 1100, 1000, 2000), 
                       ncol = 2)

yg_mosaics_cal_df <- map_dfr(yg_mosaics_cal,
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

mosaic_pre_cal_fig <- ggplot() +
  geom_spatial_rgb(data = yg_mosaics_pre_df,
                   aes(x = x, y = y, r = red, g = green, b = blue),
                   scale = 1) +
  facet_grid(cal ~ mosaic) +
  theme_bw() +
  theme(axis.title = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5)) +
  coord_sf(crs = crs(yg_mosaic_pre_2017),
           datum = crs(yg_mosaic_pre_2017),
           expand = FALSE)
mosaic_pre_cal_fig
# ggsave('figures/mosaics_uncalibrated.jpg',
#        mosaic_pre_cal_fig,
#        height = 5,
#        width = 6.5,
#        bg = 'white')

mosaic_cal_fig <- ggplot() +
  geom_spatial_rgb(data = yg_mosaics_cal_df,
                   aes(x = x, y = y, r = red, g = green, b = blue),
                   scale = 1) +
  facet_grid(cal ~ mosaic) +
  theme_bw() +
  theme(axis.title = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5)) +
  coord_sf(crs = crs(yg_mosaic_pre_2017),
           datum = crs(yg_mosaic_pre_2017),
           expand = FALSE)
mosaic_cal_fig
# ggsave('figures/mosaics_calibrated.jpg',
#        mosaic_cal_fig,
#        height = 5,
#        width = 6.5,
#        bg = 'white')
################################################################################

### Gif ########################################################################
img <- image_graph(400, 1000, res = 300)
map2(yg_mosaics_pre,
     year,
    ~ print(
      ggRGB(.x,
            r = 3,
            g = 2,
            b = 1,
            limits = matrix(c(100, 100, 300, 900, 900, 850), 
                            ncol = 2)) +
        ggtitle(.y) +
        theme_map() +
        theme(plot.title = element_text(hjust = 0.5))
    )
)
dev.off()
image_write_video(img,
                  delay = 2,
                  '~/Documents/permafrost_pathways/rts_mapping/planet_processing_test/data/unprocessed_composite/output/planet_composite.gif')
################################################################################