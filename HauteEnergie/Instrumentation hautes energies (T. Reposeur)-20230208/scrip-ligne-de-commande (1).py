import os, sys
import numpy as np

emin=100
emax=300000
roi=15
tstart=311130500
tstop=312338900
dtime=10800

# infile = '../lat_photon_weekly_w128-129-from-data-server.fits'
infile = "/home/clem/WSL_FERMI/HauteEnergie/lat_photon_weekly_w148_p305_v001.fits"
scfile = '../lat_spacecraft_weekly_w128-129-from-data-server.fits'
outprefix = 'lat_photon_3c454.3'

cmd = 'gtselect infile=%s outfile=%s_e%s-%s_roi%s.fits tmin=INDEF tmax=INDEF ra=342.894699 dec=18.81113535 rad=%s emin=%s emax=%s zmax=105' % (infile, outprefix,emin, emax, roi, roi, emin, emax)
os.system(cmd)

cmd = 'gtbin evfile=%s_e%s-%s_roi%s.fits scfile=%s outfile=%s_e%s-%s_roi%s_LC.fits algorithm=LC tbinalg=LIN tstart=%s tstop=%s dtime=%s ' %(outprefix, emin, emax, roi, scfile, outprefix, emin, emax, roi, tstart, tstop, dtime)
os.system(cmd)

cmd = 'gtexposure infile=%s_e%s-%s_roi%s_LC.fits scfile=%s irfs=CALDB srcmdl=none specin=-2.1' %(outprefix, emin, emax, roi, scfile)
os.system(cmd)
