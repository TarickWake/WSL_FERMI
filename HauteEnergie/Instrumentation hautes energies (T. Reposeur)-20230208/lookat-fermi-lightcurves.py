import os, glob, sys

# Set up matplotlib
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

from astropy.io import fits
from astropy.table import Table
import numpy as np				### numarray or numpy depends on the python install

fitsfile1='anonymous-exemple.fits'   # put here your favorite file
f1 = fits.open(fitsfile1)
fits.info(fitsfile1)
events1 = Table.read(fitsfile1,hdu=1)


# print various info (examples)
print(events1.columns)
#print(events1['ENERGY'].unit)

fig, ax = plt.subplots()


time1 = events1['TIME']
n1 = events1['COUNTS']
dn1 = events1['ERROR']
expos1 = events1['EXPOSURE']

tmin = time1 > 3.115e8 
tmax = time1 < 3.122e8
tmin = time1 > 3.115e8 
tmax = time1 < 9.122e8
tcut = tmin & tmax

n1_cut = n1[tcut]
dn1_cut = dn1[tcut]
expos1_cut = expos1[tcut]

flux1 = n1_cut/expos1_cut
dflux1 = dn1_cut/expos1_cut

ax.set_xlabel('Time (Mission Elapsed Time)', fontsize=10)
ax.set_ylabel('Flux (photon/cm2/s)', fontsize=10)
#ax.plot(time1,flux1)
ax.errorbar(time1[tcut],flux1,yerr=dflux1,fmt='o')
#plt.yscale('log')
plt.show()

exit(0)
