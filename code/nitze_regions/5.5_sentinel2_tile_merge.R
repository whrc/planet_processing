################################################################################
###                    Merge Tiled Sentinel Images                           ###
###                         Code by HGR 2/2023                               ###
################################################################################

# Load libraries
library(stringr)
library(purrr)
library(terra)

# Read in data, merge, write file, and delete inputs
files <- list.files('/explore/nobackup/people/hrodenhi/rts_mapping/planet_processing_test/data/nitze_regions/sentinel_data',
                    full.names = TRUE)
for (region in seq(0, 6)) {
  
  print('#####################################################################')
  print(paste0('Starting on region ', region))
  
  # get file names within region
  region_files <- files[which(str_detect(files, paste0('Region_', region)))]
  
  utm_zones <- unique(c(str_extract(region_files, 'UTM[:digit:]{1,2}N')))
  
  for (utm_zone in utm_zones) {
    
    print(paste0('Starting on  ', utm_zone))
    
    # get file names within region in correct UTM zone
    current_files <- region_files[which(str_detect(region_files, utm_zone))]
    
    # import raster files
    raster_list <- map(current_files,
                       ~ rast(.x))
    print(raster_list[[1]])

    # merge raster files into one raster
    if (length(raster_list) > 1) {
      print(c('Merging files: ', current_files))
      merged <- raster_list[[1]]
      for (image in raster_list[2:length(raster_list)]) {
        merged <- merge(merged, image)
      }
      print(merged)
      # save merged raster
      output_name <- paste0(
        str_extract(current_files[1], 
                    '^/.+UTM[:digit:]{1,2}N'),
        '_merged.tif')
      print(paste0('Saving merged raster to: ', output_name))
      writeRaster(merged, output_name)
      
      # remove files that were just merged
      print('Removing files.')
      file.remove(current_files)
      
    }
    
  }
  
}
