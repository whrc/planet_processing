################################################################################
###     Create Polygons to Improve Coverage of Nitze Regions                 ###
###                      Code by HGR 3/2023                                  ###
################################################################################

### Load Libraries #############################################################
library(sf)
library(tidyverse)
################################################################################

### Load Data ##################################################################
nitze_bbox <- read_sf('data/nitze_regions/bboxes/nitze_bbox_water_removed.shp')
image_footprints <- read_sf('data/nitze_regions/planet_images_footprints.shp')
################################################################################

### Functions ##################################################################
long2UTM <- function(long) {
  zone = (floor((long + 180)/6) %% 60) + 1
  return(paste0('EPSG:', 32600 + zone)) # would be 32700 for Southern Hemisphere
}

extract_utm <- function(data_sf, subset) {
  long <- data_sf %>%
    filter(region == subset) %>% 
    st_coordinates() %>% 
    as.data.frame() %>% 
    summarise(long = mean(X))
  
  utm_zone <- long2UTM(long)
  
  return(utm_zone)
}

plot_missing_area <- function(bbox, footprint, missing_area) {
  return(
    ggplot(bbox) +
      geom_sf(color = 'black', 
              fill = 'transparent',
              linewidth = 1) +
      geom_sf(data = footprints %>%
                filter(region == 2) %>%
                st_transform(crs = st_crs(bbox)),
              color = 'black',
              fill = 'transparent',
              linewidth = 0.5) +
      geom_sf(data = missing_area, 
              color = 'red', 
              fill = 'transparent',
              linewidth = 1)
  )
}
################################################################################

### Merge footprints by region #################################################
# add region information to nitze_bbox
nitze_bbox <- nitze_bbox %>%
  st_cast('POLYGON') %>%
  mutate(region = seq(1, 7))

# join the nitze_bbox region information to the footprints
image_footprints <- image_footprints %>%
  st_join(nitze_bbox)
# merge and dissolve the footprints by region
image_footprints <- image_footprints %>%
  aggregate(by = list(image_footprints$region), FUN = first) %>%
  select(region)

# Intersect the footprints with the original bboxes
counter <- 1
for (nregion in seq(1, 7)) {
  bbox <- nitze_bbox %>%
    filter(region == nregion) %>%
    st_transform(crs = extract_utm(., nregion))
  print(bbox)
  footprints <- image_footprints %>%
    filter(region == nregion) %>%
    st_transform(crs = st_crs(bbox))
  name <- paste0('nitze_bbox_r', nregion, '_improve_coverage')
  missing_area <- st_difference(bbox, footprints)
  
  if (nrow(missing_area) >= 1) {
    missing_area <- missing_area %>%
      mutate(geometry = st_as_sfc(st_bbox(geometry))) %>%
      st_intersection(bbox) %>%
      select(region)
    assign(name, missing_area)
    st_write(missing_area,
             paste0('data/nitze_regions/bboxes/', name, '.shp'),
             delete_dsn = TRUE)
    if (counter == 1) {
      nitze_bbox_improve_coverage <- missing_area %>%
        st_transform(crs = 4326)
    } else {
      nitze_bbox_improve_coverage = rbind(nitze_bbox_improve_coverage,
                                          missing_area %>%
                                            st_transform(crs = 4326))
    }
    counter <- counter + 1
  }
  plot <- plot_missing_area(bbox, footprints, missing_area)
  print(plot)
}
ggplot() +
  geom_sf(data = nitze_bbox_improve_coverage)
st_write(nitze_bbox_improve_coverage,
         'data/nitze_regions/bboxes/nitze_bbox_improve_coverage.shp',
         delete_dsn = TRUE)
################################################################################
