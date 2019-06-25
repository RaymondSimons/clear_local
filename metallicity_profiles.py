import astropy
from astropy.io import fits
import argparse
import glob
from glob import glob
import numpy as np
from numpy import *
import photutils
from photutils import detect_sources
import matplotlib
import matplotlib.pyplot as plt

def make_metal_profile(fl):
    zfits = fits.open(fl)
    try: zmap = zfits['Z_R3'].data
    except: return
    zmap = zfits['Z_R3'].data[:,:,0]
    ezmap = zfits['Z_R3'].data[:,:,1] - zmap
    imR = zfits['R3'].data
    eimR = zfits['eR3'].data
    xmn_pix = 26
    xmx_pix = 54

    imR[xmn_pix:xmx_pix, xmn_pix:xmx_pix] = nan
    eimR[xmn_pix:xmx_pix, xmn_pix:xmx_pix] = nan
    imR[(imR<0) | (eimR/imR/log(10) > 0.4)] = nan

    segm = detect_sources(imR, -99, npixels=5)

    x1 = 40.
    y1 = 40.
    pix_scale = 0.1
    x = (np.arange(0, shape(zmap)[0]) - x1) * pix_scale
    y = (np.arange(0, shape(zmap)[1]) - y1) * pix_scale

    xv, yv = np.meshgrid(x, y)
    r = sqrt(xv**2. + yv**2.)

    segm = detect_sources(imR, -99, npixels=5)

    lbl_interest = array(segm.data)[40, 40]
    if lbl_interest == 0:
        small_box = array(segm.data)[36:46, 36:46].ravel() 
        if len(small_box[small_box > 0]) > 0:  lbl_interest = min(small_box[small_box > 0])
        else: zmap[:,:] = nan

    zmap[segm.data != lbl_interest] = nan

    fig = plt.figure(figsize = (8, 8))
    ax = fig.add_subplot(1,1,1)
    print (zmap.shape)
    ax.imshow(zmap, vmin = 8, vmax = 9,  interpolation = 'nearest', cmap = 'viridis')

    fig_name = fl.split('/')[-1].replace('.fits', '.png')
    fig.savefig('/home/rsimons/figures/metal_maps/metal_maps_%s'%fig_name, dpi = 300)



if __name__ == '__main__':

    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--field', default = 'GS1')
    args = vars(parser.parse_args())

    field = args['field']

    fls = glob('/user/rsimons/metal_maps/%s_*_metals.fits'%field)
    
    for fl in fls:
        make_metal_profile(fl = fl)
