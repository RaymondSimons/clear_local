import astropy
from astropy.io import fits
import argparse
import glob
from glob import glob
import numpy as np
from numpy import *


def make_metal_profile(fl):
    zfits = fits.open(fl)
    zmap = zfits['Z_ALL'].data

    x1 = 40.
    y1 = 40.
    x = (np.arange(0, shape(zmap)[0]) - x1 + 0.5) * pix_scale
    y = (np.arange(0, shape(zmap)[1]) - y1 + 0.5) * pix_scale

    xv, yv = np.meshgrid(x, y)
    r = sqrt(xv**2. + yv**2.)

if __name__ == '__main__':

    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--field', default = 'GS1')
    args = vars(parser.parse_args())

    field = args['field']

    fls = glob('/user/rsimons/metal_maps/%s_*_metals.fits'%field)
    
    for fl in fls:
        make_metal_profile(fl = fl)
