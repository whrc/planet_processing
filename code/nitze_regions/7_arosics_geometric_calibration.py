#!/usr/bin/env python
# coding: utf-8

# # Run AROSICS Geometric Calibration on Planet Images
# ## TODO:
# - delete masked file after running AROSICS

# ## Import Libraries
import xarray as xr
import rioxarray as rxr
from geoarray import GeoArray
from arosics import COREG_LOCAL
import os
import re
from pprint import pprint


# ## Define Functions
def makemydir(dir_path):
    try:
        os.makedirs(dir_path)
    except OSError:
        pass
    os.chdir(dir_path) # this changes the working directory - necessary?


# ## Define Directories and Load Data
# Read in data
ref_dir = '/explore/nobackup/people/hrodenhi/rts_mapping/planet_processing_test/data/nitze_regions/sentinel_data'
base_dir = '/explore/nobackup/people/hrodenhi/rts_mapping/planet_processing_test/data/nitze_regions/planet_data'
ref_files = []
for root, subdirs, files in os.walk(ref_dir):
    for file in files:
        if re.match('.*tif$', file):
            ref_files.append(os.path.join(root, file))
print(ref_files)
        
tgt_files = []
for root, subdirs, files in os.walk(base_dir):
    for file in files:
        if re.match('.*SR.*tif$', file):
            tgt_files.append(os.path.join(root, file))
            
tgt_files = sorted(tgt_files)

udm_files = []
for root, subdirs, files in os.walk(base_dir):
    for file in files:
        if re.match('.*udm2.*tif$', file):
            udm_files.append(os.path.join(root, file))
            
udm_files = sorted(udm_files)

pprint('Length of tgt_files: ' + str(len(tgt_files)))
pprint('Length of udm_files: ' + str(len(udm_files)))

mask_dir = '/explore/nobackup/people/hrodenhi/rts_mapping/planet_processing_test/data/nitze_regions/arosics_output/masked_files/'
makemydir(mask_dir)
out_dir = '/explore/nobackup/people/hrodenhi/rts_mapping/planet_processing_test/data/nitze_regions/arosics_output/calibrated_files/'
makemydir(out_dir)
out_fig_dir = '/explore/nobackup/people/hrodenhi/rts_mapping/planet_processing_test/data/nitze_regions/arosics_output/figures/'
makemydir(out_fig_dir)


# use this to find the start index if the script stops running midway:
# tgt_files.index('name_of_image_that_broke_it.tif')
start_idx = 0

# ## Run AROSICS
for tgt_file, udm_file in zip(tgt_files[start_idx:], udm_files[start_idx:]):
    
    print('###########################################')
    print('Beginning with target file: ' + tgt_file)
    
    masked = os.listdir('/explore/nobackup/people/hrodenhi/rts_mapping/planet_processing_test/data/nitze_regions/arosics_output/masked_files/')
    masked = [string[2:-11] for string in masked]
    region = tgt_file.split('/')[10].split('_')[1]

    if tgt_file.split('/')[-1].split('.')[0] not in masked:
        print('Cloud masking the target image.')
        # Mask imagery using clear mask
        tgt = rxr.open_rasterio(tgt_file)
        udm = rxr.open_rasterio(udm_file).sel(band = 1)
        tgt = tgt.where(tgt != 0)
        tgt = tgt.where(udm != 0)
        tgt = tgt.astype('uint16')
        
        # save masked file
        tgt.rio.to_raster(os.path.join(mask_dir, region + '_' + tgt_file.split('/')[-1][0:-4] + '_masked.tif'))
        print('Done masking the target image.')
        
    calibrated = os.listdir('/explore/nobackup/people/hrodenhi/rts_mapping/planet_processing_test/data/nitze_regions/arosics_output/calibrated_files/')
    calibrated = [string[2:].split('_arosics')[0] for string in calibrated]
    
    if tgt_file.split('/')[-1].split('.')[0] not in calibrated:
        print('Starting AROSICS.')
        
        # get crs of tgt file to determine correct ref file
        crs = 'UTM' + str(rxr.open_rasterio(tgt_file).rio.crs).split(':')[1][3:5].lstrip('0')
        
        # get correct ref file
        # /explore/nobackup/people/hrodenhi/rts_mapping/planet_processing_test/data/nitze_regions/planet_data/region_x/unique_id/PSScene/file.tif
        ref = [file for file in ref_files if re.match('.*Region_' + region + '_' + crs + 'N.tif$', file)][0]
        print('Reference file: ' + ref)
        
        # Apply AROSICS correction from red band to all bands
        out_path = os.path.join(out_dir, region + '_' + tgt_file.split('/')[-1][0:-4] + '_arosics.tif')
        out_fig_path = os.path.join(out_fig_dir, region + '_' + tgt_file.split('/')[-1][0:-4] + '_coreg_points.tif')
        kwargs = {
            'grid_res'     : 200, # can make this larger to improve run time
            'window_size'  : (200, 200), # A larger window size could slow things down, but too small might lead to less accurate alignment
            'q'            : False,
            'r_b4match'    : 3,
            's_b4match'    : 3,
            'path_out'     : out_path,
            'fmt_out'      : 'GTiff',
            'CPUs'         : 4
        }
        
        try:
            tgt_cal = COREG_LOCAL(ref, tgt_file, **kwargs)
            tgt_cal.correct_shifts(min_points_local_corr = 4)
            tgt_cal.view_CoRegPoints(savefigPath = out_fig_path)
        except IndexError:
            print("There's not enough data to match the target and reference.")
            out_path = os.path.join(out_dir, region + '_' + tgt_file.split('/')[-1][0:-4] + '_arosics_failed.tif')
            # save masked file in AROSICS folder
            tgt.rio.to_raster(out_path)
        except Exception:
            print("Something went wrong in AROSICS.")
            out_path = os.path.join(out_dir, region + '_' + tgt_file.split('/')[-1][0:-4] + '_arosics_failed.tif')
            # save masked file in AROSICS folder
            tgt.rio.to_raster(out_path)
        
    print('\n')

# If AROSICS removed all tie points in flagging and returned a warning, the 
# try-except statment did not catch it and it did not save an output file.
# This saves the masked file in the AROSICS output folder.
masked_files = []
for root, subdirs, files in os.walk('/explore/nobackup/people/hrodenhi/rts_mapping/planet_processing_test/data/nitze_regions/arosics_output/masked_files/'):
    for file in files:
        masked_files.append(os.path.join(root, file))
masked_files = sorted(masked_files)
masked_pattern = [string.split('/')[-1].split('_masked')[0] for string in masked_files]

calibrated = os.listdir('/explore/nobackup/people/hrodenhi/rts_mapping/planet_processing_test/data/nitze_regions/arosics_output/calibrated_files/')
calibrated = [string.split('_arosics')[0] for string in calibrated]
missing_files = [file for file in masked_pattern if file not in calibrated]

for missing_file in missing_files:
    print('Copying ' + missing_file + ' to:')
    masked_file = [file for file in masked_files 
                   if re.match('^.*' + missing_file + '_masked\.tif$', 
                               file)][0]
    masked = rxr.open_rasterio(masked_file)
    out_path = masked_file.split('masked_files/')[0] + 'calibrated_files/' + masked_file.split('masked_files/')[1].split('_masked')[0] + '_arosics_failed.tif'
    print(out_path)
    masked.rio.to_raster(out_path)

