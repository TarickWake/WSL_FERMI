import os, glob, sys

# Set up matplotlib
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

from astropy.io import fits
from astropy.table import Table
import numpy as np				### numarray or numpy depends on the python install


#fitsfile='/mnt/e/data2/fermi_analysis/all-sky-10years-5gev.fits'    # put here the path to your favorite file
fitsfile=r'/Users/tiffon/OneDrive - etu.u-bordeaux.fr/S10/HauteEnergie/w749-750/sourcefiles/lat_photon_weekly_w749_p305_v001.fits' 
f1 = fits.open(fitsfile)

fits.info(fitsfile)
events = Table.read(fitsfile,hdu=1)
print(events.columns)
print(events['ENERGY'].unit)
print(events['ENERGY'][0])

fig, axs = plt.subplots(2, 1, figsize=(8, 8))

energy = events['ENERGY']
emin_mask = energy > 50000.00
l = events['L']
b = events['B']

axs[0].set_xlabel('L (deg.)', fontsize=10)
axs[0].set_ylabel('B (deg.)', fontsize=10)
counts0, xedges0, yedges0, im0 = axs[0].hist2d(l,b,bins=[360,180],norm=LogNorm())
fig.colorbar(im0,ax=axs[0])

zenith_angle = events['ZENITH_ANGLE']
zenith_angle_mask = zenith_angle < 80
l_cut = l[zenith_angle_mask] 
b_cut = b[zenith_angle_mask] 
axs[1].set_xlabel('L (deg.)', fontsize=10)
axs[1].set_ylabel('B (deg.)', fontsize=10)
counts1, xedges1, yedges1, im1 = axs[1].hist2d(l_cut,b_cut,bins=[360,180],norm=LogNorm())
fig.colorbar(im1,ax=axs[1])
plt.show()

exit(0)
