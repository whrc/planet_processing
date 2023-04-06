################################################################################
###              Reproject Nitze RTS Polygons to UtM                         ###
###                      Code by HGR 2/2023                                  ###
################################################################################

### Load Libraries #############################################################
library(sf)
library(tidyverse)
################################################################################

### Load Data ##################################################################
nitze_polygons <- read_sf('data/nitze_regions/nitze_polygons.shp')
nitze_bbox <- read_sf('data/nitze_regions/nitze_bbox.shp')
nitze_github <- read_sf('data/nitze_regions/ImageFootprints_RTS_PlanetScope_Nitze_iteration001.shp')
nitze_github <- nitze_github %>%
  st_make_valid() %>%
  filter(str_detect(image_id, '^[:digit:]{8}'))
################################################################################

### Reproject to Appropriate UTM Zone ##########################################
# name the regions
nitze_bbox <- nitze_bbox %>%
  mutate(region = seq(1, 6))

nitze_polygons <- nitze_polygons %>%
st_join(nitze_bbox, join = st_within) %>%
  arrange(region)

# take a look at the data on a map
map(nitze_bbox$region,
    ~ print(ggplot(filter(nitze_polygons,
                    region == .x)) +
              geom_sf() +
              ggtitle(paste0('Region ', .x))))

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

# reproject the polygons
nitze_polygons_by_region <- map(
  nitze_bbox$region,
  ~ nitze_polygons %>%
    filter(region == .x) %>%
    st_transform(crs = extract_utm(nitze_polygons, .x))
)

names(nitze_polygons_by_region) <- paste0('r', nitze_bbox$region)

# save the reprojected polygons
map2(nitze_polygons_by_region,
    names(nitze_polygons_by_region),
    ~ write_sf(.x,
               paste0('data/nitze_regions/nitze_polygons_', .y, '.shp')))

# reproject the bboxes
nitze_bbox_by_region <- map(
  nitze_bbox$region,
  ~ nitze_bbox %>%
    slice(.x) %>%
    st_transform(crs = st_crs(nitze_polygons_by_region[[.x]]))
)

names(nitze_bbox_by_region) <- paste0('r', nitze_bbox$region)

# save the reprojected bboxes
map2(nitze_bbox_by_region,
     names(nitze_bbox_by_region),
     ~ write_sf(.x,
                paste0('data/nitze_regions/nitze_bbox_', .y, '.shp')))
################################################################################

### Get Image Years ############################################################
bbox <- map_dfr(nitze_github %>%
                  st_intersects(nitze_bbox),
                ~ as.data.frame(list(.x), col.names = c('bbox')))
nitze_github <- nitze_github %>%
  mutate(bbox,
         year = as.numeric(str_sub(image_id, 1, 4)))

years <- nitze_github %>%
  group_by(bbox, year) %>%
  count()

# 2019 is most common for all regions
best_year <- years %>%
  group_by(bbox) %>%
  filter(n == max(n))
################################################################################