import os
import sys

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

from astropy.io import fits
from astropy.table import Table
import numpy as np

class GammaFits:
    def __init__(self, fitsfile):
        self.fitsfile = fitsfile
        self.f1 = fits.open(self.fitsfile)
        self.events = Table.read(self.fitsfile, hdu=1)
        self.energy = self.events['ENERGY']
        self.L = self.events['L']
        self.B = self.events['B']
        self.zenith_angle = self.events['ZENITH_ANGLE']
    def show_info(self):
        fits.info(self.fitsfile)
        print(self.events.columns)
        print(self.events['ENERGY'].unit)
        print(self.events['ENERGY'][0])
    
    def plot_hist2d(self):
        fig, axs = plt.subplots(2, 1, figsize=(8, 8))

        emin_mask = self.energy > 50000.00

        l = self.events['L']
        b = self.events['B']
        axs[0].set_xlabel('L (deg.)', fontsize=10)
        axs[0].set_ylabel('B (deg.)', fontsize=10)
        counts0, xedges0, yedges0, im0 = axs[0].hist2d(self.L, self.B, bins=[360,180], norm=LogNorm())
        fig.colorbar(im0, ax=axs[0])


        zenith_angle_mask = self.zenith_angle < 80
        l_cut = self.L[zenith_angle_mask] 
        b_cut = self.B[zenith_angle_mask] 

        axs[1].set_xlabel('L (deg.)', fontsize=10)
        axs[1].set_ylabel('B (deg.)', fontsize=10)
        counts1, xedges1, yedges1, im1 = axs[1].hist2d(l_cut, b_cut, bins=[360,180], norm=LogNorm())
        fig.colorbar(im1, ax=axs[1])
        plt.show()

fitsfile = "/Users/tiffon/OneDrive - etu.u-bordeaux.fr/S10/HauteEnergie/w749-750/ sourcefiles/lat_photon_weekly_w749_p305_v001.fits"

gf = GammaFits(fitsfile)
gf.show_info()
gf.plot_hist2d()