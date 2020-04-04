import astropy
from astropy.io import fits, ascii
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
import os
from astropy.table import Table
from joblib import Parallel, delayed


plt.ioff()
plt.close('all')





import astropy
from astropy.io import fits, ascii
import glob
from glob import glob
import numpy as np
from numpy import *
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy.ma as ma
import photutils
from photutils import detect_sources
plt.ioff()
plt.close('all')


def clean_metal_maps(fl):
    mm_fits = fits.open(fl)

    np.random.seed(1)

    zmax = np.inf
    figdir = '/Users/rsimons/Dropbox/clear/figures/metals/metal_maps'
    fitdir = '/Users/rsimons/Dropbox/clear/products/metals/metal_maps_cleaned'

    lines = [('OII', 'oii3726;oii3729', 3727.),
             ('OIII', 'oiii4959;oiii5007', 5007.),
             ('Hb', 'hbeta', 4863.),
             ('Ha', 'nii6548;halpha;nii6584', 6563.),
             ('SII', 'sii6717;sii6731', 6725.)
            ]

    line_cnt = 0
    for l, (line, izi_line, wave) in enumerate(lines):
        try: im = mm_fits['%s_S_EC'%line].data
        except: continue
        line_cnt+=1

    nzrows = 8
    nrows = line_cnt + nzrows
    ncols = len(fit_types)

    fig = plt.figure(figsize = (8., 2.*nrows))
    res = {}
    cmap = plt.cm.viridis
    cmap.set_bad('k')
    cmap_peaks = plt.cm.terrain
    cmap_peaks.set_bad('k')

    vmin = 8.5
    vmax = 9.5
    flux_clip = 1.

    master_hdulist = []
    prihdr = fits.Header()
    prihdr['COMMENT'] = "Storing the cleaned metallicity maps in this FITS file."
    prihdu = fits.PrimaryHDU(header=prihdr)    
    master_hdulist.append(prihdu)

    colhdr = fits.Header()

    master_hdulist[0].header['nlines'] = 0



    for ft, fit_type in enumerate(fit_types):
        ct = 0
        Zax1 = plt.subplot2grid((nrows, ncols), (0, ft))
        Zax2 = plt.subplot2grid((nrows, ncols), (1, ft))
        Zax3 = plt.subplot2grid((nrows, ncols), (2, ft))
        Zax4 = plt.subplot2grid((nrows, ncols), (3, ft))
        Zax5 = plt.subplot2grid((nrows, ncols), (4, ft))
        Zax6 = plt.subplot2grid((nrows, ncols), (5, ft))
        Dax1 = plt.subplot2grid((nrows, ncols), (6, ft))
        Sax1 = plt.subplot2grid((nrows, ncols), (7, ft))

        dir_im = mm_fits['DSCI'].data[20:60, 20:60]
        dir_im_header = mm_fits['DSCI'].header

        dir_seg = mm_fits['SEG'].data[20:60, 20:60]
        dir_seg_header = mm_fits['SEG'].header

        dir_im[dir_seg != dir_seg[20,20]] = 0.


        #create new seg
        
        segm = detect_sources(data      = dir_im, 
                              threshold = dir_im[20,20]/20., 
                              npixels   = 5,           
                              connectivity = 8)
        dir_seg = segm.data

        di_seg = dir_seg[20,20]
        dir_im_seg = dir_im.copy()
        dir_im_seg[dir_seg != di_seg] = 0.
        midx, midy = photutils.centroid_2dg(dir_im_seg)

        dir_im_header['midx'] = midx
        dir_im_header['midy'] = midy

        dir_seg_clip = dir_seg.copy()
        dir_seg_clip[dir_seg_clip != di_seg ] = 0
        cmap_dir = plt.cm.Greys_r
        cmap.set_bad('k')
        Dax1.imshow(dir_im,  cmap = cmap_dir, vmin = dir_im[20,20]/20., vmax = dir_im[20,20])
        Sax1.imshow(dir_seg_clip, cmap = cmap_dir)

        master_hdulist.append(fits.ImageHDU(data = dir_im, header = dir_im_header, name = 'DSCI'))
        master_hdulist.append(fits.ImageHDU(data = dir_seg, header = dir_seg_header, name = 'SEG'))
        master_hdulist.append(fits.ImageHDU(data = dir_seg_clip, header = dir_seg_header, name = 'SEG_CLIP'))


        Zhdr =  mm_fits['Z%s'%(fit_type)].header
        Zim = mm_fits['Z%s'%(fit_type)].data
        Z_mode = Zim[:,:,0].copy()
        Z_l68 = Zim[:,:,1]
        Z_u68 = Zim[:,:,2]
        Z_np = Zim[:,:,3]
        Z_lerr = abs(Z_l68 - Z_mode)
        Z_uerr = abs(Z_u68 - Z_mode)


        x = (np.arange(0, shape(Z_mode)[0]) - midx)
        y = (np.arange(0, shape(Z_mode)[1]) - midy)

        xv, yv = np.meshgrid(x, y)
        r = sqrt(xv**2. + yv**2.)

        res['r{}'.format(fit_type)] = r
        res['midx{}'.format(fit_type)] = midx
        res['midy{}'.format(fit_type)] = midy

        master_hdulist.append(fits.ImageHDU(data = r, header = dir_im_header, name = 'r'))


        Zax1.imshow(Z_mode, cmap = cmap, vmin = vmin,  vmax = vmax)
        Zax2.imshow(Z_lerr, cmap = cmap, vmin = 0.,    vmax = 0.5)

        Z_seg = ma.masked_where((dir_seg != di_seg ) | (Z_lerr > zmax) | (Z_uerr > zmax) | (np.isnan(Z_mode)), Z_mode)
        Zax4.imshow(Z_seg,  cmap = cmap, vmin = vmin,  vmax = vmax)


        for ax in [Zax1, Zax2, Zax3, Zax4, Zax5,Zax6, Dax1, Sax1]:
            ax.axis('off')
            ax.set_xticklabels([''])
            ax.set_yticklabels([''])
            ax.plot(midx, midy, 'x', color = 'black')




        fs = 8
        if ft == 1: Zax1.annotate('smoothed',(0.1, 0.95),va = 'top', xycoords = 'axes fraction', color = 'white', fontweight = 'bold', fontsize = fs)
        if ft == 2: Zax1.annotate('exctinction-corrected',(0.1, 0.95),va = 'top', xycoords = 'axes fraction', color = 'white', fontweight = 'bold', fontsize = fs)
        if ft == 3: Zax1.annotate('smoothed\nexctinction-corrected',(0.1, 0.95), va = 'top', xycoords = 'axes fraction', color = 'white', fontweight = 'bold', fontsize = fs)

        fs = 15
        if ft == 0:
            Zax1.annotate('Z',(0.1, 0.95), va = 'top', xycoords = 'axes fraction', color = 'white', fontweight = 'bold', fontsize = fs)
            Zax2.annotate('Z$_{err}$',(0.1, 0.95), va = 'top', xycoords = 'axes fraction', color = 'white', fontweight = 'bold', fontsize = fs)
            Zax3.annotate('Z$_{Z, seg}$',(0.1, 0.95), va = 'top', xycoords = 'axes fraction', color = 'white', fontweight = 'bold', fontsize = fs)
            Zax4.annotate('Z$_{dir, seg}$',(0.1, 0.95), va = 'top', xycoords = 'axes fraction', color = 'white', fontweight = 'bold', fontsize = fs)
            Zax5.annotate('Z$_{both, seg}$',(0.1, 0.95), va = 'top', xycoords = 'axes fraction', color = 'white', fontweight = 'bold', fontsize = fs)
            Zax6.annotate('Z$_{both, peaks}$',(0.1, 0.95), va = 'top', xycoords = 'axes fraction', color = 'white', fontweight = 'bold', fontsize = fs)
            Dax1.annotate('direct$_{F105W}$',(0.1, 0.95), va = 'top', xycoords = 'axes fraction', color = 'white', fontweight = 'bold', fontsize = fs)
            Sax1.annotate('seg$_{F105W}$',(0.1, 0.95), va = 'top', xycoords = 'axes fraction', color = 'white', fontweight = 'bold', fontsize = fs)



        try: 
            flux_im =   mm_fits['%s%s'%('OII', fit_type)].data
            flux_eim =  mm_fits['E%s%s'%('OII', fit_type)].data
        except: 
            flux_im =   mm_fits['%s%s'%('HB', fit_type)].data
            flux_eim =  mm_fits['E%s%s'%('HB', fit_type)].data

        Z_modeforseg = Z_mode.copy()
        Z_modeforseg[(Z_lerr > zmax) |\
                     (Z_uerr > zmax) |\
                     (flux_im/flux_eim < flux_clip) ] = np.nan


        #Make segmentation map off of the Zmap
        segm = detect_sources(data      = Z_modeforseg, 
                              threshold = -99, 
                              npixels   = 5,           
                              connectivity = 4)
        if not segm: 
            print ('SEGMENTATION FAILED %s %s'%(fit_type, fl.split('/')[-1]))
            Z_temp = Z_mode * np.nan
            Zax3.imshow(Z_temp, cmap = cmap, vmin = vmin,  vmax = vmax)
            Zax5.imshow(Z_temp,  cmap = cmap, vmin = vmin, vmax = vmax)
            Zax6.imshow(Z_temp,  cmap = cmap_peaks, vmin = 1, vmax = 2)

            master_hdulist.append(fits.ImageHDU(data = Z_temp, header = Zhdr, name = 'Z_clean'))
            master_hdulist.append(fits.ImageHDU(data = Z_temp, header = Zhdr, name = 'npeaks_clean'))
            master_hdulist.append(fits.ImageHDU(data = Z_temp, header = Zhdr, name = 'leZ_clean'))
            master_hdulist.append(fits.ImageHDU(data = Z_temp, header = Zhdr, name = 'ueZ_clean'))

            res['Z{}'.format(fit_type)] = Z_temp
            res['Npeaks{}'.format(fit_type)] = Z_temp
            res['elZ{}'.format(fit_type)] = Z_temp
            res['euZ{}'.format(fit_type)] = Z_temp
            continue

        else:
            lbl_interest = array(segm.data)[int(midx), int(midy)]   
            
            if lbl_interest == 0:
                small_box = array(segm.data)[int(midx - 1):int(midx +1), int(midy - 1):int(midy +1)].ravel() 
                if len(small_box[small_box > 0]) > 0:  lbl_interest = min(small_box[small_box > 0])
                else: lbl_interest = 99

        Z_show = ma.masked_where((segm.data != lbl_interest) | (Z_lerr > zmax) | (Z_uerr > zmax) | (np.isnan(Z_mode)), Z_mode)
        mask_both = (dir_seg != di_seg ) | (segm.data != lbl_interest) | (Z_lerr > zmax) | (Z_uerr > zmax) | (flux_im/flux_eim < flux_clip) |  (np.isnan(Z_mode))
        Z_both = ma.masked_where(mask_both,   Z_mode)
        eZl_both = ma.masked_where(mask_both, Z_lerr)
        eZu_both = ma.masked_where(mask_both, Z_uerr)
        Z_np_both = ma.masked_where(mask_both, Z_np)






        Zax3.imshow(Z_show, cmap = cmap, vmin = vmin,  vmax = vmax)
        Zax5.imshow(Z_both,  cmap = cmap, vmin = vmin, vmax = vmax)
        Zax6.imshow(Z_np_both,  cmap = cmap_peaks, vmin = 1, vmax = 2)
        norm = matplotlib.colors.Normalize(vmin=1, vmax=2)
        Zax6.annotate('1 peak', (0.9, 0.15), xycoords = 'axes fraction',
                      ha = 'right', va = 'bottom', color = cmap_peaks(norm(1)))
        Zax6.annotate('2 peaks',(0.9, 0.05), xycoords = 'axes fraction',
                      ha = 'right', va = 'bottom',  color = cmap_peaks(norm(2)))




        res['Z{}'.format(fit_type)] = Z_both
        res['Npeaks{}'.format(fit_type)] = Z_np_both

        res['elZ{}'.format(fit_type)] = eZl_both
        res['euZ{}'.format(fit_type)] = eZu_both


        Z_mode_save = Z_mode.copy()
        Z_np_save   = Z_np.copy()
        Z_lerr_save = Z_lerr.copy()
        Z_uerr_save = Z_uerr.copy()

        Z_mode_save[mask_both] = np.nan
        Z_np_save  [mask_both] = np.nan
        Z_lerr_save[mask_both] = np.nan
        Z_uerr_save[mask_both] = np.nan



        master_hdulist.append(fits.ImageHDU(data = Z_mode_save,  header = Zhdr, name = 'Z_clean'))
        master_hdulist.append(fits.ImageHDU(data = Z_np_save  ,  header = Zhdr, name = 'npeaks_clean'))
        master_hdulist.append(fits.ImageHDU(data = Z_lerr_save,  header = Zhdr, name = 'leZ_clean'))
        master_hdulist.append(fits.ImageHDU(data = Z_uerr_save,  header = Zhdr, name = 'ueZ_clean'))


        ct = 0
        for l, (line, izi_line, wave) in enumerate(lines):
            try: im = mm_fits['%s%s'%(line, fit_type)].data
            except: continue

            ax = plt.subplot2grid((nrows, ncols), (ct+nzrows, ft))
            ct+=1
            ax.set_xticklabels([''])
            ax.set_yticklabels([''])
            ax.axis('off')
            #if ft == 0: 
            pl = ax.imshow(im, cmap = 'viridis')
            #else: ax.imshow(im, cmap = 'viridis', vmin = pl.get_clim()[0], vmax = pl.get_clim()[1])

            fs = 15
            if ft == 0:
                ax.annotate(  '%s'%(line),(0.1, 0.95), va = 'top', xycoords = 'axes fraction', color = 'white', fontweight = 'bold', fontsize = fs)

            ax.plot(midx, midy, 'x', color = 'grey')


            line_data = mm_fits['%s%s'%(line, fit_type)].data
            line_header = mm_fits['%s%s'%(line, fit_type)].header
            master_hdulist.append(fits.ImageHDU(data = line_data,  header = line_header, name = '%s%s'%(line, fit_type)))

            line_data_new = line_data.copy()
            line_data_new[mask_both] = np.nan
            master_hdulist.append(fits.ImageHDU(data = line_data_new,  header = line_header, name = '%s%s_seg'%(line, fit_type)))

        master_hdulist[0].header['nlines'] = ct


    figname = figdir + '/' + fl.split('/')[-1].replace('fits', 'png')
    fitsname = fitdir + '/' + fl.split('/')[-1].replace('.fits', '_cleaned.fits')
    print ('saving %s'%figname)
    fig.tight_layout()
    fig.savefig(figname)


    print ('\tSaving to ' + fitsname)
    thdulist = fits.HDUList(master_hdulist)
    thdulist.writeto(fitsname, overwrite = True)


    plt.close('all')
    return res


def make_clean_metal_map_figure(fld, di, fit_types):
    indir = '/Users/rsimons/Dropbox/clear/products/metals/metal_maps_cleaned'
    fl = indir + '/{}_{}_metals_highZbranch_cleaned.fits'.format(fld, di)
    mm = fits.open(fl)

    midx = mm[1].header['midx']
    midy = mm[1].header['midy']
    ncols =  5#mm[0].header['nlines']
    if ncols == 0.: return
    nrows =  2
    fig = plt.figure(figsize = (ncols*2., nrows * 2.))

    names = []
    for b in mm:
        if 'S_EC_SEG' in b.name:
            names.append(b.name.rstrip('_S_EC_SEG'))

    ln_names = ['OII', 'HB', 'OIII', 'HA', 'SII']
    for l, ln in enumerate(ln_names):
        if ln in names: 
            ax = plt.subplot2grid((nrows, ncols), (0, int(l)))
            im = mm[ln + '_S_EC'].data
            ax.set_xticklabels([''])
            ax.set_yticklabels([''])
            ax.axis('off')
            ax.imshow(im, cmap = 'viridis')
            fs = 15

            ax.annotate(  '%s'%(ln),(0.1, 0.95), va = 'top', xycoords = 'axes fraction', color = 'white', fontweight = 'bold', fontsize = fs)
            ax.plot(midx, midy, 'x', color = 'grey')

    cmap = plt.cm.viridis
    cmap.set_bad('k')
    cmap_dir = plt.cm.Greys_r
    cmap_dir.set_bad('k')

    vmin = 8.5
    vmax = 9.5
    Dax = plt.subplot2grid((nrows, ncols), (1, 0))
    Zax = plt.subplot2grid((nrows, ncols), (1, 1))

    Dax.imshow(mm['DSCI'].data, cmap = cmap_dir)
    Zax.imshow(mm['Z_clean'].data, cmap = cmap, vmin = vmin,  vmax = vmax)
    for ax in [Dax, Zax]:
        ax.set_xticklabels([''])
        ax.set_yticklabels([''])
        ax.axis('off')
        ax.plot(midx, midy, 'x', color = 'grey')

    figdir = '/Users/rsimons/Dropbox/clear/figures/metals/metal_maps_cleaned'
    fig.tight_layout()
    fig.savefig(figdir + '/' + fl.split('/')[-1].replace('.fits', '.png'), dpi = 400)

def write_metal_profile(fld, di, fit_types):
    indir = '/Users/rsimons/Dropbox/clear/products/metals/metal_maps'
    fl = indir + '/{}_{}_metals_highZbranch.fits'.format(fld, di)
    res = clean_metal_maps(fl)
    np.save('/Users/rsimons/Dropbox/clear/products/metals/metal_profiles/{}_{}_highZbranch.npy'.format(fld, di), res)
def fit_metal_profile(fld, di, fit_types):
    fl = '/Users/rsimons/Dropbox/clear/products/metals/metal_profiles/{}_{}_highZbranch.npy'.format(fld, di)
    if not os.path.isfile(fl): return
    a = np.load(fl, allow_pickle = True)[()]
    try: x = a['Z_S_EC'].mask
    except: return
    nrows = len(fit_types)
    ncols = 4
    fig = plt.figure(figsize = (ncols*2., nrows*2.))
    res = {}


    for ft, fit_type in enumerate(fit_types):
        res['p{}'.format(fit_type)]  = np.array([np.nan, np.nan] )
        res['V{}'.format(fit_type)]  = np.array([[np.nan, np.nan],[np.nan, np.nan]])
        res['r{}'.format(fit_type)]  = np.array([np.nan, np.nan])
        res['Z{}'.format(fit_type)]  = np.array([np.nan, np.nan])
        res['eZ{}'.format(fit_type)] = np.array([np.nan, np.nan])

        ax_Z =  plt.subplot2grid((nrows, ncols), (ft, 0))
        ax_nP =  plt.subplot2grid((nrows, ncols), (ft, 1))
        ax_p =  plt.subplot2grid((nrows, ncols), (ft, 2))
        ax_pf =  plt.subplot2grid((nrows, ncols), (ft, 3))

        if ft == 3:
            ax_p.set_xlabel('distance from center (pix)')
            ax_pf.set_xlabel('distance from center (pix)')
        ax_p.set_ylabel(r'$\log$ (O/H) + 12')


        cmap = plt.cm.viridis
        cmap.set_bad('k')
        cmap_peaks = plt.cm.terrain
        cmap_peaks.set_bad('k')

        vmin = 8.5
        vmax = 9.5

        Zmap = a['Z{}'.format(fit_type)]
        Npeaks_im = a['Npeaks{}'.format(fit_type)]
        Npeaks = Npeaks_im.ravel()

        ax_Z.imshow(Zmap, cmap = cmap, vmin = vmin,  vmax = vmax)

        ax_nP.imshow(Npeaks_im,  cmap = cmap_peaks, vmin = 1, vmax = 2)

        Zmap_mask = Zmap.mask

        elZmap = a['elZ{}'.format(fit_type)]
        euZmap = a['euZ{}'.format(fit_type)]

        crit1 = (~Zmap_mask.ravel()) & (~isnan(Zmap.data.ravel()))

        r = a['r{}'.format(fit_type)].ravel()[crit1]
        Z = Zmap.data.ravel()[crit1]
        elZ = elZmap.data.ravel()[crit1]
        euZ = euZmap.data.ravel()[crit1]

        elZ = array([max(eZ, 0.3) for eZ in elZ])
        seuZ = array([max(eZ, 0.3) for eZ in euZ])
        ax_p.errorbar(r, Z, yerr = [elZ, euZ], linestyle = 'None', fmt = 'o', color = 'grey', alpha = 0.2, zorder = 0,  markersize = 2)


        ax_p.set_ylim(6.9, 9.6)


        for ax in [ax_Z, ax_nP]:
            ax.plot(a['midx{}'.format(fit_type)], a['midy{}'.format(fit_type)], 'x', color = 'red')
            ax.axis('off')

        if len(r) == 0: 
            ax_p.axis('off')
            ax_pf.axis('off')

        if len(Z) > 5:
            outl = nanstd(Z)

            def reject_outliers(data_for_mdev, data, m = 2.):
                d = np.abs(data_for_mdev - np.median(data_for_mdev))
                mdev = np.median(d)
                d2 = np.abs(data - np.median(data_for_mdev))
                s = d2/mdev if mdev else 0.
                return s<m

            def reject_by_number(Z_use, Z):
                above_8 = len(where(Z_use > 8)[0])
                below_8 = len(where(Z_use < 8)[0])

 
                if above_8/(len(Z_use)) > 0.7: 
                    gd = where(Z > 8)
                elif below_8/(len(Z_use)) > 0.7:
                    gd = where(Z < 8)
                else:
                    gd = []
                return gd

            #gd_tofit = reject_outliers(Z[r < 4], Z)
            gd_tofit = reject_by_number(Z[r < 10], Z)

            r   = r  [gd_tofit]
            Z   = Z  [gd_tofit]
            elZ = elZ[gd_tofit]
            euZ = euZ[gd_tofit]


            if len(Z) > 5.:
                ax_p.errorbar(r, Z, yerr = [elZ, euZ], linestyle = 'None', fmt = 'o', color = 'blue', alpha = 0.5, zorder = 1 , markersize = 2)
                ax_pf.errorbar(r, Z, yerr = [elZ, euZ], linestyle = 'None', fmt = 'o', color = 'blue', alpha = 0.5, zorder = 1, markersize = 2)

                ax_pf.set_xlim(ax_p.get_xlim())
                if (max(r) - min(r) > 2) & (len(where(r < 4.)[0]) > 3.):
                    eZ = np.mean((elZ, euZ), axis= 0 )

                    p, V = np.polyfit(r, Z, deg = 1., w = 1./(eZ**2.), cov = True)
                    draws = np.random.multivariate_normal(p, V, size = 100)
                    x = linspace(0, max(r) + 1., 100)
                    res['p{}'.format(fit_type)]  = p
                    res['V{}'.format(fit_type)]  = V
                    res['r{}'.format(fit_type)]  = r
                    res['Z{}'.format(fit_type)]  = Z
                    res['eZ{}'.format(fit_type)] = eZ
                    for d in draws:
                        ax_pf.plot(x, x*d[0] + d[1], color = 'black', alpha = 0.05, zorder = 2)



                    ylm_min = max(min(Z-elZ) - 0.2, 7.0)
                    ylm_max = min(max(Z+euZ) + 0.2, 9.5)

                    ax_pf.set_ylim(ylm_min, ylm_max)
                    ax_p.axhline(y = ylm_min, xmin = 0.8, xmax = 1.0, linestyle = '--', color = 'grey')
                    ax_p.axhline(y = ylm_max, xmin = 0.8, xmax = 1.0, linestyle = '--', color = 'grey')

    figdir = '/Users/rsimons/Dropbox/clear/figures/metals/metal_profiles'



    npsavename = fl.split('/')[-1].replace('.npy', '_profile_fit.npy')
    np.save('/Users/rsimons/Dropbox/clear/products/metals/metal_profiles/fits/{}'.format(npsavename), res)
    figname = figdir + '/' + fl.split('/')[-1].replace('npy', 'png')


    print ('saving %s'%figname)
    fig.tight_layout()
    fig.savefig(figname, dpi = 400)

    plt.close('all')




if __name__ == '__main__':

    izi_cat = ascii.read('/Users/rsimons/Dropbox/clear/catalogs/good_izi.cat', header_start = 0)
    fit_types = array(['', '_S', '_EC', '_S_EC'])
    fit_types = array(['_S_EC'])


    izi_cat = izi_cat
    if True:
        if False:        
            for f, (fld, di) in enumerate(zip(izi_cat['field'], izi_cat['id'])):
                if True:#(fld == 'GS1') & (di == 46685):
                    write_metal_profile(fld, di, fit_types)
                    make_clean_metal_map_figure(fld, di, fit_types)
                    #fit_metal_profile(fld, di, fit_types)
        else:
            Parallel(n_jobs = -1)(delayed(write_metal_profile)(fld, di, fit_types) for f, (fld, di) in enumerate(zip(izi_cat['field'], izi_cat['id'])))        
            #Parallel(n_jobs = -1)(delayed(make_clean_metal_map_figure)(fld, di, fit_types) for f, (fld, di) in enumerate(zip(izi_cat['field'], izi_cat['id'])))        
            #Parallel(n_jobs = -1)(delayed(fit_metal_profile)(fld, di, fit_types) for f, (fld, di) in enumerate(zip(izi_cat['field'], izi_cat['id'])))

    if False:
        #Write metal profile catalog
        catalog_dic = {}
        catalog_dic['field'] = []
        catalog_dic['id'] = []

        for fit_type in fit_types:
            catalog_dic['m{}'.format(fit_type)] = []
            catalog_dic['m{}_err'.format(fit_type)] = []
            catalog_dic['b{}'.format(fit_type)] = []
            catalog_dic['b{}_err'.format(fit_type)] = []
        izi_cat = ascii.read('/Users/rsimons/Dropbox/clear/catalogs/good_izi_new.cat', header_start = 0)

        for f, (fld, di) in enumerate(zip(izi_cat['field'], izi_cat['id'])):
            fl = '/Users/rsimons/Dropbox/clear/products/metals/metal_profiles/fits/{}_{}_highZbranch_profile_fit.npy'.format(fld, di)
            if not os.path.isfile(fl): continue
            res = np.load(fl, allow_pickle = True)[()]
            catalog_dic['field'].append(fld)
            catalog_dic['id'].append(di)
            for fit_type in fit_types:
                catalog_dic['m{}'.format(fit_type)].append(res['p{}'.format(fit_type)][0])
                catalog_dic['m{}_err'.format(fit_type)].append(sqrt(res['V{}'.format(fit_type)][0,0]))
                catalog_dic['b{}'.format(fit_type)].append(res['p{}'.format(fit_type)][1])
                catalog_dic['b{}_err'.format(fit_type)].append(sqrt(res['V{}'.format(fit_type)][1,1]))

        data= Table(catalog_dic)
        ascii.write(catalog_dic, '/Users/rsimons/Dropbox/clear/catalogs/metal_highZbranch_profile_fits.cat', format = 'commented_header')

        


















