import astropy
from astropy.io import fits
from astropy.cosmology import Planck15 as cosmo
import argparse
import glob
from glob import glob
import numpy as np
from numpy import *
import photutils
from photutils import detect_sources
import matplotlib
import matplotlib.pyplot as plt
#from joblib import Parallel, delayed

plt.ioff()
plt.close('all')
def make_metal_profile(fl, grizli_cat):
    di = fl.split('/')[-1].split('_')[1]
    fld = fl.split('/')[-1].split('_')[0]
    z = grizli_cat[1].data['z_50'][grizli_cat[1].data['ID'] == int(di)][0]
    fl_check = glob('/Users/rsimons/Desktop/clear/bad/*{}*{}*'.format(fld, di))
    if len(fl_check) > 0: return

    zfits = fits.open(fl)
    try: zmap = zfits['Z_R3'].data
    except: return
    xmn_pix = 30
    xmx_pix = 50
    midx = (xmx_pix - xmn_pix)/2

    zmap = zfits['Z_R3'].data[xmn_pix:xmx_pix, xmn_pix:xmx_pix,0]
    ezmap = zfits['Z_R3'].data[xmn_pix:xmx_pix, xmn_pix:xmx_pix,2]
    imR = zfits['R3'].data[xmn_pix:xmx_pix, xmn_pix:xmx_pix]
    eimR = zfits['eR3'].data[xmn_pix:xmx_pix, xmn_pix:xmx_pix]



    O3 = zfits['OIII'].data[xmn_pix:xmx_pix, xmn_pix:xmx_pix]
    Hb = zfits['HB'].data[xmn_pix:xmx_pix, xmn_pix:xmx_pix]

    fit_hdu = fits.open('/Users/rsimons/Dropbox/metal_maps_full/{}_{}_reduced_full.fits'.format(fld, di))
    direct_im = fit_hdu['DSCI'].data[xmn_pix:xmx_pix, xmn_pix:xmx_pix]
    seg_im = fit_hdu['SEG'].data[xmn_pix:xmx_pix, xmn_pix:xmx_pix]
    direct_im_temp = direct_im.copy()
    direct_im_temp[seg_im != seg_im[int(shape(seg_im)[0]/2.), int(shape(seg_im)[0]/2.)]] = nan
    midx, midy = photutils.centroid_2dg(direct_im_temp)

    #midx, midy = photutils.centroid_2dg(O3)
    if (midx < 0) | (midy < 0): return
    if (midx > xmx_pix - xmn_pix) | (midy > xmx_pix - xmn_pix): return


    imR[(imR<0) | (eimR/imR/log(10) > 0.30)] = nan

    segm = detect_sources(imR, -99, npixels=5)

    #x1 = 14.
    #y1 = 14.
    pix_scale = 0.1 / cosmo.arcsec_per_kpc_proper(z)
    x = (np.arange(0, shape(zmap)[0]) - midx) * pix_scale.value
    y = (np.arange(0, shape(zmap)[1]) - midy) * pix_scale.value

    xv, yv = np.meshgrid(x, y)
    r = sqrt(xv**2. + yv**2.)

    segm = detect_sources(imR, -99, npixels=5)

    lbl_interest = array(segm.data)[int(midx), int(midy)]
    if lbl_interest == 0:
        small_box = array(segm.data)[int(midx - 1):int(midx +1), int(midy - 1):int(midy +1)].ravel() 
        if len(small_box[small_box > 0]) > 0:  lbl_interest = min(small_box[small_box > 0])
        else: lbl_interest = 99

    zmap[segm.data != lbl_interest] = nan
    ezmap[segm.data != lbl_interest] = nan
    r[segm.data != lbl_interest] = nan

    '''
    bad = where(seg_im != seg_im[int(shape(seg_im)[0]/2.), int(shape(seg_im)[0]/2.)])
    zmap[bad] = nan
    ezmap[bad] = nan
    r[bad] = nan
    '''

    rad_z = {}
    gd = where(~isnan(zmap.ravel()))[0]
    if (len(gd) > 15.) & (len(gd) < 50.):    
        r_gd = r.ravel()[gd]
        z_gd = zmap.ravel()[gd]
        ez_gd = ezmap.ravel()[gd]
        rad_z['r'] = r_gd
        rad_z['z'] = z_gd
        ez_gd = array([max([ez, 0.3]) for ez in ez_gd])
        rad_z['ez'] = ez_gd
        if nanmax(r_gd) > 3.5:
            outname = fl.split('/')[-1].replace('.fits', '.npy')
            p, V = np.polyfit(r_gd, z_gd, deg = 1., w = 1./(ez_gd**2.), cov = True)
            rad_z['p'] = p
            rad_z['V'] = V

            fig2 = plt.figure(figsize = (30, 5))
            #ax2 = fig2.add_subplot(1,1,1 )
            axim = plt.subplot2grid((1, 6), (0, 0), colspan=1)   
            axo3 = plt.subplot2grid((1, 6), (0, 1), colspan=1)   
            axhb = plt.subplot2grid((1, 6), (0, 2), colspan=1)                
            ax = plt.subplot2grid((1, 6), (0, 3), colspan=1)        
            ax2 = plt.subplot2grid((1, 6), (0, 4), colspan=2)        
            #ax = fig.add_subplot(1,1,1)
            
            cmap = plt.cm.viridis
            cmap.set_bad('k')
            axim.imshow(direct_im,  interpolation = 'nearest', cmap = 'Greys_r')
            axo3.imshow(O3,  interpolation = 'nearest', cmap = cmap)
            axhb.imshow(Hb,  interpolation = 'nearest', cmap = cmap)

            ax.imshow(zmap, vmin = 8, vmax = 9,  interpolation = 'nearest', cmap = cmap)
            for axx in [axim, ax, axo3, axhb]:
                axx.plot(midx, midy, 'x', color = 'grey')



            draws = np.random.multivariate_normal(p, V, size = 100)
            x = linspace(0, 10, 100)
            ax2.errorbar(r_gd, z_gd, yerr = ez_gd, color = 'blue', fmt = 'o', alpha = 1.0, markersize = 8, linewidth = 1., zorder = 2)

            for d in draws:
                ax2.plot(x, x*d[0] + d[1], color = 'blue', alpha = 0.1)

            ax2.set_xlim(0,10)
            ax2.set_xlabel('distance from center, kpc')
            ax2.set_ylabel('log (O/H) + 12')


            np.save('/Users/rsimons/Desktop/clear/z_radius_new/{}'.format(outname), rad_z)
            fig_name = fl.split('/')[-1].replace('.fits', '.png')

            fig2.savefig('/Users/rsimons/Desktop/clear/mm_figs/metal_radius_%s'%fig_name, dpi = 300)


if __name__ == '__main__':

    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--field', default = 'GS1')
    args = vars(parser.parse_args())

    field = args['field']

    #fields = ['GS1', 'GS2', 'GS3', 'GS4', 'GS5', 'ERSPRIME']
    fields = ['GN1', 'GN2', 'GN3', 'GN4', 'GN5', 'GN7']
    #fields = ['GS1']#, 'GN2', 'GN3', 'GN4', 'GN5', 'GN7']


    for field in fields:
        print (field)
        fls = glob('/Users/rsimons/Dropbox/metal_maps_reduced/%s_*_metals.fits'%field)
        grizli_cat = fits.open('/Users/rsimons/Desktop/clear/Catalogs/grizli_v2.1_cats/{}_lines_grizli.fits'.format(field))
        for fl in fls:#3:4]: 
            make_metal_profile(fl = fl, grizli_cat = grizli_cat)













