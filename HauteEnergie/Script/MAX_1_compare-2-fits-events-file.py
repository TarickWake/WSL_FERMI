import os, glob, sys

# Set up matplotlib
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

from astropy.io import fits
from astropy.table import Table
import numpy as np				### numarray or numpy depends on the python install

fitsfile1='HauteEnergie/Data/source_observation/lat_photon_weekly_w015_p305_v001.fits'  # path to your file number one
f1 = fits.open(fitsfile1)
fits.info(fitsfile1)
events1 = Table.read(fitsfile1,hdu=1)

fitsfile2='HauteEnergie/Data/source_observation/lat_photon_weekly_w016_p305_v001.fits'  # path to your file number one
f2 = fits.open(fitsfile2)
fits.info(fitsfile2)
events2 = Table.read(fitsfile2,hdu=1)

# print various info (examples)
#print(events1.columns)
#print(events1['ENERGY'].unit)
#print(events1['ENERGY'][0])

fig, axs = plt.subplots(3, 1, figsize=(8, 8))

energy1 = events1['ENERGY']
emin_mask1 = energy1 > 50000.00
zenith_angle1 = events1['ZENITH_ANGLE']
zenith_angle_mask1 = zenith_angle1 < 100
time1 = events1['TIME']
time1_min = time1 > 3.108e8 
time1_max = time1 < 3.120e8
time1_cut = time1_min & time1_max
#time1_cut = time1_max
l1 = events1['L']
b1 = events1['B']
#mask1_tot = zenith_angle_mask1 & time1_cut
mask1_tot = zenith_angle_mask1
l1_cut = l1[mask1_tot] 
b1_cut = b1[mask1_tot] 
print(len(l1_cut))


energy2 = events2['ENERGY']
emin_mask2 = energy2 > 50000.00
zenith_angle2 = events2['ZENITH_ANGLE']
zenith_angle_mask2 = zenith_angle2 < 100
time2 = events2['TIME']
time2_min = time2 > 3.108e8 
time2_max = time2 < 3.120e8
time2_cut = time2_min & time2_max
#time2_cut = time2_max
mask2_tot = zenith_angle_mask2
l2 = events2['L']
b2 = events2['B']
l2_cut = l2[mask2_tot] 
b2_cut = b2[mask2_tot] 
print(len(l2_cut))

k=1

axs[0].set_xlabel('L (deg.)', fontsize=10)
axs[0].set_ylabel('B (deg.)', fontsize=10)
counts0, xedges0, yedges0, im0 = axs[0].hist2d(l1_cut,b1_cut,bins=k*[36,18],norm=LogNorm())
fig.colorbar(im0,ax=axs[0])

axs[1].set_xlabel('L (deg.)', fontsize=10)
axs[1].set_ylabel('B (deg.)', fontsize=10)
counts1, xedges1, yedges1, im1 = axs[1].hist2d(l2_cut,b2_cut,bins=k*[36,18],norm=LogNorm())
fig.colorbar(im1,ax=axs[1])

diff = counts0 - counts1
ratio = counts0 / counts1
#from PIL import Image
#pil_image=Image.fromarray(diff)
#pil_image.show()

#fig = plt.figure()
#plt.imshow(ratio)
#plt.title("Plot 2D array")
#plt.show()

axs[2].set_xlabel('L (deg.)', fontsize=10)
axs[2].set_ylabel('B (deg.)', fontsize=10)
axs[2] = diff, xedges1, yedges1, im0
# fig.colorbar(im2,ax=axs[2])

plt.show()