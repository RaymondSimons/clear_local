import numpy as np
from numpy import *
import astropy
from astropy.io import fits
import glob
from glob import glob
gals = ['VELA%.2i'%i for i in arange(1, 36)]


runs_dir = '/astro/snyder_lab2/New_HydroART_images/VELA_v2'
for g, gal in enumerate(gals):
    if gal != 'VELA18':
        print g, gal
        fls = glob(runs_dir + '/%s/%s_a*_sunrise/images/broadbandz.fits'%(gal, gal))
        for f, fl in enumerate(fls):
            a = fl.split('/')[-3].split('_')[-2].lstrip('a')
            print a, fl