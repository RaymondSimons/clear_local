import astropy
from astropy.io import fits
import argparse
import glob
from glob import glob
import numpy as np
from numpy import *
import photutils
from photutils import detect_sources


def make_metal_profile(fl):
    zfits = fits.open(fl)
    try: zmap = zfits['Z_R3'].data
    except: return
    imR = zfits['R3'].data
    eimR = zfits['eR3'].data
    imR[xmn_pix:xmx_pix, xmn_pix:xmx_pix] = nan
    eimR[xmn_pix:xmx_pix, xmn_pix:xmx_pix] = nan
    imR[(imR<0) | (eimR/imR/log(10) > 0.4)] = nan


    segm = detect_sources(imR, -99, npixels=5)

    x1 = 40.
    y1 = 40.
    pix_scale = 0.1
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
