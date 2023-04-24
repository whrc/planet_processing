#!/usr/bin/env python
# coding: utf-8


# ## Import Libraries
import numpy as np
import pandas as pd
import xarray as xr
import rioxarray as rxr
from osgeo import gdal, gdalconst
from numpy.linalg import inv, eig
from scipy.stats import chi2
import time
from sklearn.cluster import KMeans
import imageio
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.linear_model import LinearRegression
import os
import re
import psutil
import gc

gc.set_debug(gc.DEBUG_SAVEALL)

# IRMAD functions from https://github.com/ChenHongruixuan/ChangeDetectionRepository
'''
Python implementation of IRMAD
A. A. Nielsen, “The regularized iteratively reweighted MAD method for change detection in multi- and hyperspectral data,” IEEE Trans. Image Process., vol. 16, no. 2, pp. 463–478, 2007.
'''
def covw(center_X, center_Y, w):
    n = w.shape[1]
    sqrt_w = np.sqrt(w)
    sum_w = w.sum()
    V = np.concatenate((center_X, center_Y), axis=0)
    V = sqrt_w * V
    dis = np.dot(V, V.T) / sum_w * (n / (n - 1))

    return dis

# max_iter = 1 is the same as MAD
def IRMAD(img_X, img_Y, max_iter=50, epsilon=1e-3):
    bands_count_X, num = img_X.shape

    weight = np.ones((1, num))  # (1, height * width)
    can_corr = 100 * np.ones((bands_count_X, 1))
    for _iter in range(max_iter):
        print(_iter)
        mean_X = np.sum(weight * img_X, axis=1, keepdims=True) / np.sum(weight)
        mean_Y = np.sum(weight * img_Y, axis=1, keepdims=True) / np.sum(weight)

        # centralization
        center_X = img_X - mean_X
        center_Y = img_Y - mean_Y
        
        # also can use np.cov, but the result would be sightly different with author' result acquired by MATLAB code
        cov_XY = covw(center_X, center_Y, weight)
        size = cov_XY.shape[0]
        sigma_11 = cov_XY[0:bands_count_X, 0:bands_count_X]  # + 1e-4 * np.identity(3)
        sigma_22 = cov_XY[bands_count_X:size, bands_count_X:size]  # + 1e-4 * np.identity(3)
        sigma_12 = cov_XY[0:bands_count_X, bands_count_X:size]  # + 1e-4 * np.identity(3)
        sigma_21 = sigma_12.T

        tgt_mat = np.dot(np.dot(np.dot(inv(sigma_11), sigma_12), inv(sigma_22)), sigma_21)
        eigenvalue, eigenvector_X = eig(tgt_mat)  # the eigenvalue and eigenvector of image X
        # sort eigenvector based on the size of eigenvalue
        eigenvalue = np.sqrt(eigenvalue)

        idx = eigenvalue.argsort()
        eigenvalue = eigenvalue[idx]

        if (_iter + 1) == 1:
            print('Canonical correlations')
        print(eigenvalue)

        eigenvector_X = eigenvector_X[:, idx]

        eigenvector_Y = np.dot(np.dot(inv(sigma_22), sigma_21), eigenvector_X)  # the eigenvector of image Y

        # tune the size of X and Y, so the constraint condition can be satisfied
        norm_X = np.sqrt(1 / np.diag(np.dot(eigenvector_X.T, np.dot(sigma_11, eigenvector_X))))
        norm_Y = np.sqrt(1 / np.diag(np.dot(eigenvector_Y.T, np.dot(sigma_22, eigenvector_Y))))
        eigenvector_X = norm_X * eigenvector_X
        eigenvector_Y = norm_Y * eigenvector_Y

        mad_variates = np.dot(eigenvector_X.T, center_X) - np.dot(eigenvector_Y.T, center_Y)  # (6, width * height)

        if np.max(np.abs(can_corr - eigenvalue)) < epsilon:
            break
        can_corr = eigenvalue
        # calculate chi-square distance and probility of unchanged
        mad_var = np.reshape(2 * (1 - can_corr), (bands_count_X, 1))
        chi_square_dis = np.sum(mad_variates * mad_variates / mad_var, axis=0, keepdims=True)
        weight = 1 - chi2.cdf(chi_square_dis, bands_count_X)

    if (_iter + 1) == max_iter:
        print('the canonical correlation may not be converged')
    else:
        print('the canonical correlation is converged, the iteration is %d' % (_iter + 1))

    return mad_variates, can_corr, mad_var, eigenvector_X, eigenvector_Y, \
           sigma_11, sigma_22, sigma_12, chi_square_dis, weight


def get_binary_change_map(data):
    """
    get binary change map
    :param data:
    :param method: cluster method
    :return: binary change map
    """

    cluster_center = KMeans(n_clusters=2, max_iter=1500).fit(data.T).cluster_centers_.T  # shape: (1, 2)
    # cluster_center = k_means_cluster(weight, cluster_num=2)
    print('k-means cluster is done, the cluster center is ', cluster_center)
    dis_1 = np.linalg.norm(data - cluster_center[0, 0], axis=0, keepdims=True)
    dis_2 = np.linalg.norm(data - cluster_center[0, 1], axis=0, keepdims=True)

    bcm = np.copy(data)  # binary change map
    if cluster_center[0, 0] > cluster_center[0, 1]:
        bcm[dis_1 > dis_2] = 1
        bcm[dis_1 <= dis_2] = 0
    else:
        bcm[dis_1 > dis_2] = 0
        bcm[dis_1 <= dis_2] = 1

    return bcm # 1 = invariant pixel


# reprojection function
def reproj_to_sentinel(tgt, ref):
    return(tgt.rio.reproject_match(ref))
    
# linear model function
def lm_coefs(df):
    x = np.array(df[ ~np.isnan(df['tgt'])]['tgt']).reshape(-1, 1)
    y = np.array(df[ ~np.isnan(df['ref'])]['ref']).reshape(-1, 1)
    return pd.DataFrame({'tgt_file': pd.Series(df['tgt_file'][0]),
                         'tgt_band': pd.Series(df['tgt_band'][0]),
                         'intercept': pd.Series(float(LinearRegression().fit(x, y).intercept_)),
                         'slope': pd.Series(float(LinearRegression().fit(x, y).coef_))})


# function to create a directory
def makemydir(dir_path):
    try:
        os.makedirs(dir_path)
    except OSError:
        pass


# Read in data
tgt_dir = '/explore/nobackup/people/hrodenhi/rts_mapping/planet_processing_test/data/yg_val_regions/arosics_output/calibrated_files/'
tgt_files = []
for root, subdirs, files in os.walk(tgt_dir):
    for file in files:
        if re.match('.*SR.*tif$', file):
            tgt_files.append(os.path.join(root, file))
            
tgt_files = sorted(tgt_files)
print('There are ' + str(len(tgt_files)) + ' target files.')

ref_dir = '/explore/nobackup/people/hrodenhi/rts_mapping/planet_processing_test/data/yg_val_regions/sentinel_data/'
ref_files = []
for root, subdirs, files in os.walk(ref_dir):
    for file in files:
        if re.match('.*tif$', file):
            ref_files.append(os.path.join(root, file))
        
mad_out_dir = '/explore/nobackup/people/hrodenhi/rts_mapping/planet_processing_test/data/yg_val_regions/mad_output/mad_files/'
makemydir(mad_out_dir)
out_dir = '/explore/nobackup/people/hrodenhi/rts_mapping/planet_processing_test/data/yg_val_regions/mad_output/calibrated_files/'
makemydir(out_dir)
lm_data_out = '/explore/nobackup/people/hrodenhi/rts_mapping/planet_processing_test/data/yg_val_regions/mad_output/linear_models.csv'

# check for files that have already been calibrated
cal_files = os.listdir(out_dir)
cal_files = [file.split('_mad')[0] for file in cal_files if re.match('^.*tif$', file)]
tgt_files = [file for file in tgt_files if file.split('/')[11].split('.')[0] not in cal_files]

try:
  failed_files = list(pd.read_csv(failed_files_out).failed_files)
  tgt_files = [file for file in tgt_files if file.split('/')[11].split('.')[0] not in failed_files]
except:
  pass

print('There are ' + str(len(tgt_files)) + ' target files that have not been calibrated.')

# Run MAD calibration
for tgt_file in tgt_files[0:5]:
    print('Reading in: ' + tgt_file)
     
    # read in target file
    bands = [1, 2, 3, 4]
    tgt = rxr.open_rasterio(tgt_file).sel(band = bands)
    
    # find and read in reference file (sentinel-2)
    crs = str(tgt.rio.crs).split(':')[1][3:5].lstrip('0')
    region = tgt_file.split('/')[-1].split('_')[0]
    ref_path = [file for file in ref_files if re.match('.*Region_' + region + '_UTM' + crs + 'N\\.tif$', file)][0]
    print('Reading in: ' + ref_path)
    ref = rxr.open_rasterio(ref_path)
    
    print('Reprojecting and formatting data for MAD')
    # reproject target file to match Sentinel-2 composite (important for grid not projection)
    tgt_reproj = reproj_to_sentinel(tgt, ref)
    
    # mask target file to ensure nodata is represented by np.nan
    tgt_reproj = tgt_reproj.where(tgt_reproj != 0)
    
    # extract image shape
    img_bands, img_height, img_width = tgt_reproj.shape 
    
    # reshape data to be a 1D array for each band
    tgt_values = tgt_reproj.values.reshape(img_bands, -1)
    ref_values = ref.values.reshape(img_bands, -1)
    
    # remove missing data for MAD algorithm
    non_missing_values = np.where(~np.isnan(tgt_values) & 
                                  ~np.isnan(ref_values))
    tgt_values = tgt_values[non_missing_values].reshape(img_bands, -1)
    ref_values = ref_values[non_missing_values].reshape(img_bands, -1)
    
    
    # run MAD
    try:
      print('Running MAD')
      mad, can_coo, mad_var, ev_1, ev_2, sigma_11, sigma_22, sigma_12, chi_2, noc_weight = IRMAD(tgt_values, 
                                                                                                 ref_values,
                                                                                                 max_iter=20,
                                                                                                 epsilon=1e-3)
      k_means_bcm = get_binary_change_map(np.sqrt(chi_2)).reshape(-1)
      
      print('Reshaping and saving MAD output')
      # reshape MAD output to reflect shape of input including missing data locations
      output = np.empty((img_height*img_width))
      output[:] = np.nan
      output[non_missing_values[1][np.where(non_missing_values[0] == 0)]] = k_means_bcm
      output = np.reshape(output, (img_height, img_width))
      
      # create a geotiff from MAD results
      mad_output = tgt_reproj.sel(band = 1)
      mad_output.values = output
      mad_output.attrs['long_name'] = 'mad_mask'
      mad_output.astype('uint16')
      
      # Save MAD output
      mad_output.rio.to_raster(mad_out_dir + tgt_file.split('/')[-1][0:-4] + '_mad.tif')
      
      print('Running linear models to calibrate radiometry')
      # Run linear models on invariant pixels identified with MAD and apply linear correction to each band
      tgt_masked = tgt_reproj.where(mad_output == 1)
      ref_masked = ref.where(mad_output == 1)
      
      tgt_cal = tgt.where(tgt != 0)
      calibrated_values = list()
      lm_output = pd.DataFrame()
      for band in bands:
          # reshape by band
          band_name = tgt_reproj.attrs['long_name'][band-1]
          tgt_temp = tgt_masked.sel(band = band).values.reshape(-1)
          ref_temp = ref_masked.sel(band = band).values.reshape(-1)
          tgt_temp_clean = tgt_temp[~np.isnan(tgt_temp) & ~np.isnan(ref_temp)]
          ref_temp_clean = ref_temp[~np.isnan(tgt_temp) & ~np.isnan(ref_temp)]
          invariant_pixels = pd.DataFrame({
              'tgt_file': tgt_file,
              'tgt_band': band_name,
              'tgt': tgt_temp_clean,
              'ref': ref_temp_clean
          })
          
          # run linear model
          lm_current = lm_coefs(invariant_pixels)
          print(lm_current)
          
          # add linear model output to dataframe for export
          lm_output = pd.concat([lm_output, lm_current], axis = 0)
          
          # calculate calibrated raster values
          intercept = lm_current['intercept'][0]
          slope = lm_current['slope'][0]
          calibrated_values.append([
              intercept + slope * x 
              for x in tgt_cal.sel(band = band).values.reshape(-1)
          ])
      
      # reshape calibrated values to reflect output raster shape
      calibrated_values = np.array(calibrated_values)
      calibrated_values = calibrated_values.reshape(tgt_cal.shape)
      
      # overwrite raster values with calibrated values
      tgt_cal.values = calibrated_values
      tgt_cal.astype('uint16')
      
      # save calibrated file
      print('Saving calibrated file')
      tgt_cal.rio.to_raster(
        out_dir + tgt_file.split('/')[-1][0:-4] + '_mad.tif', 
        driver = 'COG'
        )
      
      # update linear model metadata file
      lm_output.to_csv(lm_data_out,
                       index = False,
                       mode = 'a',
                       header = not os.path.exists(lm_data_out))
  
      print('Memory usage: ' 
            + str(psutil.Process().memory_info().rss/(1024*1024)/1000) 
            + 'GB')
    except Exception:
      print('MAD Failed')
      failed_files = pd.DataFrame({'failed_files': [tgt_file.split('/')[-1][0:-4]]})
      failed_files.to_csv(failed_files_out,
                          index = False,
                          mode = 'a',
                          header = not os.path.exists(failed_files_out))
      
    print('\n')
