import numpy as np
from numpy import *
import astropy
from astropy.io import fits
import glob
from glob import glob
gals = ['VELA%.2i'%i for i in arange(1, 36)]



runs_dir = '/astro/snyder_lab2/New_HydroART_images/VELA_v2'
for g, gal in enumerate(gals):
    print gal
    if gal != 'VELA18':
        cat = open('../../sunrise_centers/%s_sunrise_centers.cat'%gal, 'w+')
        cat.write('#%s, galaxy centers used for SUNRISE calculations\n'%gal)        
        cat.write('simname\tscale\tx\ty\tz\n')
        fls = glob(runs_dir + '/%s/%s_a*_sunrise/images/broadbandz.fits'%(gal, gal))
        for f, fl in enumerate(fls):
            a = fl.split('/')[-3].split('_')[-2].lstrip('a')
            print a
            data = fits.open(fl)
            x = data['SFRHIST'].header['translate_originX']
            y = data['SFRHIST'].header['translate_originY']
            z = data['SFRHIST'].header['translate_originZ']
            cat.write('%s\t%s\t%.2f\t%.2f\t%.2f\n'%(gal, a, x, y, z))

        cat.close()