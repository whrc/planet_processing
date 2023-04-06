#!/usr/bin/env python
# coding: utf-8

# # Fix COGs exported from GEE with empty overviews
# https://issuetracker.google.com/issues/186412979?pli=1
# ## Import Libraries
import xarray as xr
import rioxarray as rxr
import os

composite_dir = '/explore/nobackup/people/hrodenhi/rts_mapping/planet_processing_test/data/nitze_regions/calibrated_composites/'
composite_files = []
for root, subdirs, files in os.walk(composite_dir):
    for file in files:
        composite_files.append(os.path.join(root, file))

print(len(composite_files))
for idx, file in enumerate(composite_files):
  print(str(idx) + ': ' + file)
  data = rxr.open_rasterio(file)
  data.rio.to_raster(file, driver = 'COG')
