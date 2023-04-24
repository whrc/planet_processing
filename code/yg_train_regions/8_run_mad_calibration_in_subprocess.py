#!/usr/bin/env python
# coding: utf-8

### TODO:
# change data dir back

# 7_mad_radiometric_calibration.py has a memory leak caused by reading in data
# using (rio)xarray. It seemed easiest to workaround this issue by splitting
# up the image processing into smaller subprocesses, thereby freeing the memory
# after each subprocess. Since it seems unlikely that we will need to do this
# kind of processing long-term, it did not make sense to fix the code by 
# adapting it to use a different package. If we do end up using this code in the
# future, a real fix is needed.

# Import Libraries
import gc
import os
import re
import subprocess
import math

# function to create a directory
def makemydir(dir_path):
    try:
        os.makedirs(dir_path)
    except OSError:
        pass

# Find target files
tgt_dir = '/explore/nobackup/people/hrodenhi/rts_mapping/planet_processing_test/data/yg_train_regions/arosics_output/calibrated_files/'
# tgt_dir = 'data/yg_train_regions/arosics_output/calibrated_files/'
tgt_files = []
for root, subdirs, files in os.walk(tgt_dir):
    for file in files:
        if re.match('.*SR.*tif$', file):
            tgt_files.append(os.path.join(root, file))
            
tgt_files = sorted(tgt_files)

# check for files that have already been calibrated
out_dir = '/explore/nobackup/people/hrodenhi/rts_mapping/planet_processing_test/data/yg_train_regions/mad_output/calibrated_files/'
makemydir(out_dir)
cal_files = os.listdir(out_dir)
cal_files = [file.split('_mad')[0] for file in cal_files]
            
tgt_files = [file for file in tgt_files if file.split('/')[11].split('.')[0] not in cal_files]

# chunk up tgt_files for subprocesses of 5 images
chunk_size = 5
iterations = math.ceil(len(tgt_files)/chunk_size)
# iterate over chunks, running 7_mad_radiometric_calibration.py each time
for chunk in range(0, iterations):
    print('Starting subprocess ' + str(chunk), flush = True)
    result = subprocess.run(['python', '7_mad_radiometric_calibration.py'], 
                            capture_output=True, 
                            text=True,
                            timeout = 3600)
    print(result.stdout)
    print(result.stderr)
    gc.collect()
